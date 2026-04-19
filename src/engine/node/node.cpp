#include "engine/node/node.h"
#include "engine/plant.h"
#include "engine/sugar.h"
#include "engine/world_params.h"
#include "engine/chemical/chemical_registry.h"
#include "engine/vascular.h"
#include "engine/perf_log.h"
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/meristems/meristem_types.h"
#include <algorithm>
#include <cmath>
#include <glm/geometric.hpp>

namespace botany {

// --- Node (base) ---

Node::Node(uint32_t id, NodeType type, glm::vec3 position, float radius)
    : id(id)
    , parent(nullptr)
    , offset(position)
    , position(position)
    , radius(radius)
    , type(type)
    , age(0)
{
    // Initialize all chemical map entries to zero
    chemicals[ChemicalID::Auxin] = 0.0f;
    chemicals[ChemicalID::Cytokinin] = 0.0f;
    chemicals[ChemicalID::Gibberellin] = 0.0f;
    chemicals[ChemicalID::Sugar] = 0.0f;
    chemicals[ChemicalID::Ethylene] = 0.0f;
    chemicals[ChemicalID::Stress] = 0.0f;
    chemicals[ChemicalID::Water] = 0.0f;
}

void Node::add_child(Node* child) {
    children.push_back(child);
    child->parent = this;
}

void Node::replace_child(Node* old_child, Node* new_child) {
    for (auto& c : children) {
        if (c == old_child) {
            c = new_child;
            new_child->parent = this;
            old_child->parent = nullptr;

            // Transfer canalization biases to new child entry
            auto it_flow = auxin_flow_bias.find(old_child);
            if (it_flow != auxin_flow_bias.end()) {
                auxin_flow_bias[new_child] = it_flow->second;
                auxin_flow_bias.erase(it_flow);
            }
            return;
        }
    }
}

void Node::tick(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    age++;
    sync_world_position();

    float s0 = chemical(ChemicalID::Sugar);
    deduct_maintenance_sugar(world);
    tick_sugar_maintenance = s0 - chemical(ChemicalID::Sugar);

    float s1 = chemical(ChemicalID::Sugar);
    update_tissue(plant, world);
    tick_sugar_activity = chemical(ChemicalID::Sugar) - s1;

    sync_world_position();  // re-sync after tissue may have changed offset
    if (update_physics(plant, g, world)) return;

    float s2 = chemical(ChemicalID::Sugar);
    transport_chemicals(g);
    tick_sugar_transport = chemical(ChemicalID::Sugar) - s2;

    // Flush chemicals received from parent this tick (anti-teleportation:
    // received chemicals couldn't cascade during this node's own transport
    // because they were held in the buffer, not in chemical()).
    for (auto& [id, amount] : transport_received) {
        chemical(id) += amount;
    }
    transport_received.clear();

    // Starvation check runs after the flush so that vascular sugar delivered
    // this tick (held in transport_received until now) is credited before we
    // decide whether the node is dying. A node that receives enough sugar
    // each tick to cover maintenance must not accumulate starvation ticks.
    if (check_starvation(plant, world)) return;

    // Clear growth reserve — update_tissue() has already used the sugar it needed.
    sugar_reserved_for_growth = 0.0f;
}

// --- Tick helpers ---

void Node::sync_world_position() {
    position = parent ? parent->position + offset : offset;
    if (!std::isfinite(position.x) || !std::isfinite(position.y) || !std::isfinite(position.z)) {
        position = parent ? parent->position : glm::vec3(0.0f);
    }
    // Shoot-side tissue can't penetrate the ground
    bool is_underground_type = (type == NodeType::ROOT || type == NodeType::ROOT_APICAL);
    if (!is_underground_type && position.y < 0.0f) {
        offset.y -= position.y;  // adjust offset so position lands at y=0
        position.y = 0.0f;
    }
}

bool Node::pay_maintenance(Plant& plant, const WorldParams& world) {
    deduct_maintenance_sugar(world);
    return check_starvation(plant, world);
}

void Node::transport_chemicals(const Genome& g) {
    transport_with_children(g);
    update_canalization(g);
    decay_chemicals(g);
}

void Node::deduct_maintenance_sugar(const WorldParams& world) {
    float cost = maintenance_cost(world);
    chemical(ChemicalID::Sugar) = std::max(0.0f, chemical(ChemicalID::Sugar) - cost);
}


bool Node::check_starvation(Plant& plant, const WorldParams& world) {
    if (chemical(ChemicalID::Sugar) <= 0.0f) starvation_ticks++;
    else starvation_ticks = 0;

    // Type-specific death threshold.  Meristems (APICAL, ROOT_APICAL) use the
    // base world threshold and are further protected by quiescence (see
    // update_tissue).  Non-meristem woody tissue (STEM, ROOT) survives much
    // longer on stored parenchyma reserves — modeled here as a flat longer
    // threshold until starch reserves are properly simulated.  LEAF uses the
    // base threshold since leaves don't store starch reserves meaningfully.
    const Genome& g = plant.genome();
    uint32_t death_threshold = world.starvation_ticks_max;
    if (type == NodeType::STEM)      death_threshold = g.starvation_ticks_max_stem;
    else if (type == NodeType::ROOT) death_threshold = g.starvation_ticks_max_root;

    if (starvation_ticks >= death_threshold && parent) {
        die(plant);
        return true;
    }
    return false;
}

bool Node::update_physics(Plant& plant, const Genome& g, const WorldParams& world) {
    compute_mass(g, world);
    compute_stress(g, world);
    return apply_droop_and_break(plant, g, world);
}

void Node::compute_mass(const Genome& g, const WorldParams& world) {
    float self_mass = 0.0f;
    bool is_underground = (type == NodeType::ROOT || type == NodeType::ROOT_APICAL);

    if (!is_underground) {
        if (type == NodeType::LEAF) {
            auto* leaf = as_leaf();
            self_mass = leaf ? (leaf->leaf_size * leaf->leaf_size * world.leaf_mass_density) : 0.0f;
        } else if (type == NodeType::STEM) {
            float length = std::max(glm::length(offset), 0.01f);
            self_mass = 3.14159f * radius * radius * length * g.wood_density;
        } else if (is_meristem()) {
            self_mass = world.meristem_mass;
        }
    }

    total_mass = self_mass;
    mass_moment = self_mass * position;
    for (const Node* child : children) {
        total_mass += child->total_mass;
        mass_moment += child->mass_moment;
    }
}

void Node::compute_stress(const Genome& g, const WorldParams& world) {
    bool is_underground = (type == NodeType::ROOT || type == NodeType::ROOT_APICAL);
    stress = 0.0f;

    if (is_underground || position.y <= world.ground_support_height) return;

    // Compute bending stress from off-center child mass
    float self_mass = total_mass;
    for (const Node* child : children) self_mass -= child->total_mass;
    // ^ recovers self_mass without recomputing — total_mass = self + children

    float child_mass = total_mass - self_mass;
    if (child_mass > 1e-6f && radius > 1e-6f) {
        glm::vec3 child_com = (mass_moment - self_mass * position) / child_mass;
        float dx = child_com.x - position.x;
        float dz = child_com.z - position.z;
        float lever_arm = std::sqrt(dx * dx + dz * dz);
        float torque = child_mass * world.gravity * lever_arm;
        float cross_section = 3.14159f * radius * radius;
        stress = torque / cross_section;
    }

    // Stress hormone production
    if (stress > 0.0f) {
        float break_stress = g.wood_density * world.break_strength_factor;
        float stress_ratio = stress / break_stress;
        if (stress_ratio > g.stress_hormone_threshold) {
            float excess = (stress_ratio - g.stress_hormone_threshold)
                         / (1.0f - g.stress_hormone_threshold);
            chemical(ChemicalID::Stress) += excess * g.stress_hormone_production_rate;
        }
    }
}

bool Node::apply_droop_and_break(Plant& plant, const Genome& g, const WorldParams& world) {
    if (type != NodeType::STEM) return false;
    bool is_underground = (type == NodeType::ROOT || type == NodeType::ROOT_APICAL);
    if (is_underground) return false;

    float break_stress = g.wood_density * world.break_strength_factor;
    float droop_threshold = break_stress * g.wood_flexibility;

    // Ground support — stems at or below ground level rest on the surface,
    // transferring their load to the ground. No droop, no break.
    bool ground_supported = position.y <= world.ground_support_height;

    // Break check — need at least a grandparent (don't snap the first internode)
    if (!ground_supported && stress >= break_stress && parent && parent->parent) {
        die(plant);
        return true;
    }

    // Droop — rotate offset toward gravity (ground-supported stems don't droop)
    if (!ground_supported && stress > droop_threshold) {
        float excess = (stress - droop_threshold) / (break_stress - droop_threshold);
        float droop_angle = std::min(excess * world.droop_rate, world.droop_rate);
        float len = glm::length(offset);
        if (len > 1e-4f) {
            glm::vec3 dir = offset / len;
            glm::vec3 down(0.0f, -1.0f, 0.0f);
            glm::vec3 axis = glm::cross(dir, down);
            float axis_len = glm::length(axis);
            if (axis_len > 1e-6f) {
                axis /= axis_len;
                float c = std::cos(droop_angle);
                float s = std::sin(droop_angle);
                glm::vec3 new_dir = dir * c
                    + glm::cross(axis, dir) * s
                    + axis * glm::dot(axis, dir) * (1.0f - c);
                offset = glm::normalize(new_dir) * len;
            }
        }
    }

    // Elastic recovery — spring back toward rest direction when understressed
    float rest_len = glm::length(rest_offset);
    if (rest_len > 1e-4f && stress < droop_threshold) {
        glm::vec3 rest_dir = rest_offset / rest_len;
        float len = glm::length(offset);
        if (len > 1e-4f) {
            glm::vec3 cur_dir = offset / len;
            float cos_angle = glm::dot(cur_dir, rest_dir);
            if (cos_angle < 0.9999f) {
                float angle_gap = std::acos(std::min(cos_angle, 1.0f));
                float stiffness = std::min(radius * 10.0f, 2.0f);
                float recovery = std::min(g.elastic_recovery_rate * stiffness, angle_gap);

                glm::vec3 axis = glm::cross(cur_dir, rest_dir);
                float axis_len = glm::length(axis);
                if (axis_len > 1e-6f) {
                    axis /= axis_len;
                    float c = std::cos(recovery);
                    float s = std::sin(recovery);
                    glm::vec3 new_dir = cur_dir * c
                        + glm::cross(axis, cur_dir) * s
                        + axis * glm::dot(axis, cur_dir) * (1.0f - c);
                    offset = glm::normalize(new_dir) * len;
                }
            }
        }
    }

    return false;
}

void Node::update_tissue(Plant& /*plant*/, const WorldParams& /*world*/) {}

void Node::compute_growth_reserve(const Genome& /*g*/, const WorldParams& /*world*/) {
    sugar_reserved_for_growth = 0.0f;
}

float Node::get_bias_multiplier(Node* child, const Genome& g) const {
    float flow = 0.0f;
    auto it_f = auxin_flow_bias.find(child);
    if (it_f != auxin_flow_bias.end()) flow = it_f->second;
    return 1.0f + g.canalization_weight * flow;
}

float Node::get_parent_auxin_flow_bias() const {
    if (!parent) {
        float max_bias = 0.0f;
        for (const auto& [child, bias] : auxin_flow_bias)
            max_bias = std::max(max_bias, bias);
        return max_bias;
    }
    auto it = parent->auxin_flow_bias.find(const_cast<Node*>(this));
    if (it == parent->auxin_flow_bias.end()) return 0.0f;
    return it->second;
}

void Node::die(Plant& plant) {
    // Detach from parent
    if (parent) {
        parent->auxin_flow_bias.erase(this);
        auto& siblings = parent->children;
        siblings.erase(std::remove(siblings.begin(), siblings.end(), this), siblings.end());
        parent = nullptr;
    }

    // Collect self + all descendants, queue for deferred removal
    std::vector<Node*> to_die = {this};
    size_t i = 0;
    while (i < to_die.size()) {
        for (Node* child : to_die[i]->children) {
            to_die.push_back(child);
        }
        i++;
    }

    // Clear children so recursive tick snapshot is empty
    children.clear();

    for (Node* n : to_die) {
        plant.queue_removal(n);
    }
}

float Node::maintenance_cost(const WorldParams& /*world*/) const {
    return 0.0f;
}

void Node::transport_with_children(const Genome& g) {
    if (children.empty()) return;

    bool is_seed = (parent == nullptr);

    // NOTE: last_auxin_flux.clear() was here. Moved to end of update_canalization()
    // so that PIN-recorded flux (written before the DFS walk) is preserved and
    // accumulated together with diffusion flux before canalization reads it.

    float ref_radius = g.initial_radius;
    float parent_cap_sugar = sugar_cap(*this, g);
    float parent_cap_water = water_cap(*this, g);

    for (const auto& dp : diffusion_params(g)) {
        bool has_cap = (dp.id == ChemicalID::Sugar || dp.id == ChemicalID::Water);

        auto cap_for = [&](const Node& n) -> float {
            if (dp.id == ChemicalID::Water) return water_cap(n, g);
            if (dp.id == ChemicalID::Sugar) return sugar_cap(n, g);
            return 0.0f;
        };
        float parent_cap = has_cap ? (dp.id == ChemicalID::Water ? parent_cap_water : parent_cap_sugar) : 0.0f;

        // --- Per-child info needed for both phases ---
        struct ChildInfo {
            Node* child;
            float child_radius;
            float child_cap;
            float bias_mult;
            ChemicalDiffusionParams child_dp;
            float desired;   // positive = child wants to receive, negative = child wants to give
        };
        std::vector<ChildInfo> infos;
        infos.reserve(children.size());

        // Track whether this parent has vasculature (for skipping vascular chemicals)
        bool parent_vascular = has_vasculature(*this, g);

        for (Node* child : children) {
            // Skip vascular chemicals on mature-to-mature edges — the global
            // vascular pass already handled bulk flow. Local diffusion still
            // handles last-mile delivery to leaves, meristems, and young nodes.
            if (is_vascular_chemical(dp.id) && parent_vascular && has_vasculature(*child, g)) {
                continue;
            }

            float child_radius = child->radius;
            if (child->type == NodeType::LEAF) {
                auto* leaf = child->as_leaf();
                if (leaf) {
                    child_radius = std::max(child->radius, leaf->leaf_size * ref_radius);
                }
            }

            float child_cap = has_cap ? cap_for(*child) : 0.0f;

            ChemicalDiffusionParams child_dp = dp;
            if (dp.bias != 0.0f && (child->type == NodeType::ROOT ||
                                     child->type == NodeType::ROOT_APICAL)) {
                child_dp.bias = -dp.bias;
            }

            float desired = compute_transport_flow(
                child->chemical(dp.id), chemical(dp.id),
                child_cap, parent_cap,
                child_radius, radius, ref_radius, child_dp);

            float bm = get_bias_multiplier(child, g);
            infos.push_back({child, child_radius, child_cap, bm, child_dp, desired});
        }

        // --- Phase 1: children giving to parent (desired < 0) ---
        float parent_val = chemical(dp.id);
        float parent_headroom = has_cap ? std::max(0.0f, parent_cap - parent_val) : 1e30f;
        float total_inflow = 0.0f;

        for (auto& ci : infos) {
            if (ci.desired >= 0.0f) continue;
            float give = std::min(-ci.desired, ci.child->chemical(dp.id));
            ci.desired = -give;
            total_inflow += give;
        }

        // Bias-weighted redistribution when total exceeds parent's capacity
        if (total_inflow > parent_headroom && total_inflow > 1e-8f) {
            float total_weighted = 0.0f;
            for (auto& ci : infos) {
                if (ci.desired >= 0.0f) continue;
                total_weighted += (-ci.desired) * ci.bias_mult;
            }
            if (total_weighted > 1e-8f) {
                for (auto& ci : infos) {
                    if (ci.desired >= 0.0f) continue;
                    float raw_give = -ci.desired;
                    float share = parent_headroom * (raw_give * ci.bias_mult / total_weighted);
                    ci.desired = -std::min(share, raw_give);
                }
            }
            total_inflow = parent_headroom;
        }

        // Apply inflows
        for (auto& ci : infos) {
            if (ci.desired >= 0.0f) continue;
            float give = -ci.desired;
            ci.child->chemical(dp.id) -= give;
            parent_val += give;
            if (dp.id == ChemicalID::Auxin) {
                last_auxin_flux[ci.child] += give;
            }
        }
        chemical(dp.id) = parent_val;

        // --- Phase 2: parent giving to children (desired > 0) ---
        // Seed pass-through: recompute desired flows using updated parent level
        // after Phase 1 inflows, so receivers see the full transit amount.
        // Non-seed: use pre-computed desired flows.
        if (is_seed) {
            for (auto& ci : infos) {
                if (ci.desired <= 0.0f) continue; // only recompute receivers
                ci.desired = compute_transport_flow(
                    ci.child->chemical(dp.id), chemical(dp.id),
                    ci.child_cap, parent_cap,
                    ci.child_radius, radius, ref_radius, ci.child_dp);
                if (ci.desired < 0.0f) ci.desired = 0.0f; // was a receiver, stays a receiver
            }
        }

        float available = chemical(dp.id);

        struct Receiver {
            Node* child;
            float raw_weight;
            float weight;
            float headroom;
        };
        std::vector<Receiver> receivers;
        float sum_receiver_caps = 0.0f;
        float sum_receiver_vals = 0.0f;
        for (auto& ci : infos) {
            if (ci.desired <= 0.0f) continue;
            float headroom = has_cap
                ? std::max(0.0f, ci.child_cap - ci.child->chemical(dp.id))
                : 1e30f;
            if (headroom > 1e-8f) {
                receivers.push_back({ci.child, ci.desired, ci.desired * ci.bias_mult, headroom});
                sum_receiver_caps += ci.child_cap;
                sum_receiver_vals += ci.child->chemical(dp.id);
            }
        }

        // Multi-way equalization cap: prevents oscillation when per-child
        // equalization clamps sum to more than the parent has. Only kicks in
        // when the total demand actually exceeds the equilibrium budget —
        // in chain topologies (single child), the per-pair clamp is sufficient
        // and more permissive, so sugar flows through without piling up.
        float max_give = available;
        if (!receivers.empty()) {
            // First compute what per-pair flows would demand in total
            float raw_total_demand = 0.0f;
            for (const auto& r : receivers) raw_total_demand += r.raw_weight;

            // Compute the multi-way equilibrium budget
            float equil_budget = available;
            float total_system = available + sum_receiver_vals;
            if (has_cap) {
                float total_cap = parent_cap + sum_receiver_caps;
                float parent_equil = (total_cap > 1e-8f)
                    ? total_system * parent_cap / total_cap : 0.0f;
                equil_budget = std::max(0.0f, available - parent_equil);
            } else {
                float parent_equil = total_system / (1.0f + static_cast<float>(receivers.size()));
                equil_budget = std::max(0.0f, available - parent_equil);
            }

            // Only apply the cap when per-pair totals would exceed it
            if (raw_total_demand > equil_budget) {
                max_give = equil_budget;
            }
        }

        while (!receivers.empty() && available > 1e-8f && max_give > 1e-8f) {
            float raw_total = 0.0f;
            float bias_total = 0.0f;
            for (const auto& r : receivers) {
                raw_total += r.raw_weight;
                bias_total += r.weight;
            }
            if (bias_total < 1e-8f) break;

            float to_distribute = std::min({available, raw_total, max_give});
            bool any_filled = false;

            for (auto& r : receivers) {
                float share = to_distribute * (r.weight / bias_total);
                float actual = std::min(share, r.headroom);
                // All nodes (including seed) defer to received buffer.
                // Anti-teleportation: child's update_tissue() must not see
                // transit chemicals — the buffer is flushed after transport.
                r.child->transport_received[dp.id] += actual;
                available -= actual;
                max_give -= actual;
                r.headroom -= actual;
                if (r.headroom < 1e-8f) any_filled = true;
                if (dp.id == ChemicalID::Auxin) {
                    last_auxin_flux[r.child] += actual;
                }
            }

            if (!any_filled) break;

            receivers.erase(
                std::remove_if(receivers.begin(), receivers.end(),
                    [](const Receiver& r) { return r.headroom < 1e-8f; }),
                receivers.end());
        }

        chemical(dp.id) = available;
    }
}

void Node::update_canalization(const Genome& g) {
    for (Node* child : children) {
        // Get auxin flux for this child this tick (PIN + diffusion combined).
        float flux = 0.0f;
        auto it = last_auxin_flux.find(child);
        if (it != last_auxin_flux.end()) flux = it->second;

        // Saturation = fraction of this connection's PIN capacity currently used.
        // r² is the child's cross-sectional area; pin_capacity_per_area is the max
        // flux density. Clamped to [0, 1] — physical PIN can't exceed full saturation.
        float r2 = child->radius * child->radius;
        float capacity = r2 * g.pin_capacity_per_area;
        float saturation = (capacity > 1e-8f)
            ? std::min(flux / capacity, 1.0f) : 0.0f;

        // Exponential smoothing: bias chases saturation at smoothing_rate.
        // Natural decay: when flux stops, saturation→0 and bias lerps toward 0
        // over ~1/smoothing_rate ticks. No separate decay param needed.
        float& flow_bias = auxin_flow_bias[child];
        flow_bias += (saturation - flow_bias) * g.smoothing_rate;
    }

    // NOTE: last_auxin_flux is cleared at the start of Plant::tick_tree() rather
    // than here, so flux values survive the full tick and can be read by the
    // post-tick canalization logger in engine.cpp.
}

void Node::decay_chemicals(const Genome& g) {
    for (const auto& dp : diffusion_params(g)) {
        if (dp.decay_rate > 0.0f) {
            chemical(dp.id) *= (1.0f - dp.decay_rate);
        }
    }
}

bool Node::is_meristem() const {
    return type == NodeType::APICAL || type == NodeType::ROOT_APICAL;
}

// Downcasting helpers
StemNode*       Node::as_stem()       { return type == NodeType::STEM ? static_cast<StemNode*>(this) : nullptr; }
const StemNode* Node::as_stem() const { return type == NodeType::STEM ? static_cast<const StemNode*>(this) : nullptr; }
RootNode*       Node::as_root()       { return type == NodeType::ROOT ? static_cast<RootNode*>(this) : nullptr; }
const RootNode* Node::as_root() const { return type == NodeType::ROOT ? static_cast<const RootNode*>(this) : nullptr; }
LeafNode*       Node::as_leaf()       { return type == NodeType::LEAF ? static_cast<LeafNode*>(this) : nullptr; }
const LeafNode* Node::as_leaf() const { return type == NodeType::LEAF ? static_cast<const LeafNode*>(this) : nullptr; }
ApicalNode*       Node::as_apical()       { return type == NodeType::APICAL ? static_cast<ApicalNode*>(this) : nullptr; }
const ApicalNode* Node::as_apical() const { return type == NodeType::APICAL ? static_cast<const ApicalNode*>(this) : nullptr; }
RootApicalNode*       Node::as_root_apical()       { return type == NodeType::ROOT_APICAL ? static_cast<RootApicalNode*>(this) : nullptr; }
const RootApicalNode* Node::as_root_apical() const { return type == NodeType::ROOT_APICAL ? static_cast<const RootApicalNode*>(this) : nullptr; }

} // namespace botany
