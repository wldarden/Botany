#include "engine/node/node.h"
#include "engine/plant.h"
#include "engine/sugar.h"
#include "engine/world_params.h"
#include "engine/chemical/chemical_registry.h"
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
            auto it_struct = structural_flow_bias.find(old_child);
            if (it_struct != structural_flow_bias.end()) {
                structural_flow_bias[new_child] = it_struct->second;
                structural_flow_bias.erase(it_struct);
            }

            return;
        }
    }
}

void Node::tick(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    age++;
    sync_world_position();
    if (handle_energy_cost(plant, world)) return;
    tissue_tick(plant, world);
    if (update_physics(plant, g, world)) return;
    update_chemicals(g);
}

// --- Tick helpers ---

void Node::sync_world_position() {
    position = parent ? parent->position + offset : offset;
    if (!std::isfinite(position.x) || !std::isfinite(position.y) || !std::isfinite(position.z)) {
        position = parent ? parent->position : glm::vec3(0.0f);
    }
}

bool Node::handle_energy_cost(Plant& plant, const WorldParams& world) {
    pay_maintenance(world);
    return check_starvation(plant, world);
}

void Node::update_chemicals(const Genome& g) {
    transport_with_children(g);
    decay_chemicals(g);
}

void Node::pay_maintenance(const WorldParams& world) {
    float cost = maintenance_cost(world);
    chemical(ChemicalID::Sugar) = std::max(0.0f, chemical(ChemicalID::Sugar) - cost);
}


bool Node::check_starvation(Plant& plant, const WorldParams& world) {
    if (chemical(ChemicalID::Sugar) <= 0.0f) starvation_ticks++;
    else starvation_ticks = 0;

    if (starvation_ticks >= world.starvation_ticks_max && parent) {
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

    // Break check — ground-anchored stems can't snap
    bool can_break = position.y > world.ground_support_height && parent && parent->parent;
    if (can_break && stress >= break_stress) {
        die(plant);
        return true;
    }

    // Droop — rotate offset toward gravity
    if (stress > droop_threshold) {
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

void Node::tissue_tick(Plant& /*plant*/, const WorldParams& /*world*/) {}

float Node::get_bias_multiplier(Node* child, const Genome& g) const {
    float flow = 0.0f, structural = 0.0f;
    auto it_f = auxin_flow_bias.find(child);
    if (it_f != auxin_flow_bias.end()) flow = it_f->second;
    auto it_s = structural_flow_bias.find(child);
    if (it_s != structural_flow_bias.end()) structural = it_s->second;
    return 1.0f + g.canalization_weight * (flow + structural);
}

void Node::die(Plant& plant) {
    // Detach from parent
    if (parent) {
        parent->auxin_flow_bias.erase(this);
        parent->structural_flow_bias.erase(this);
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

    float ref_radius = g.initial_radius;
    float parent_cap_sugar = sugar_cap(*this, g);

    for (const auto& dp : diffusion_params(g)) {
        bool has_cap = (dp.id == ChemicalID::Sugar);

        // --- Compute desired flow for each child ---
        struct ChildFlow {
            Node* child;
            float desired;   // positive = child wants to receive, negative = child wants to give
            float child_cap;
            float bias_mult; // canalization weight multiplier for this connection
        };
        std::vector<ChildFlow> flows;
        flows.reserve(children.size());

        for (Node* child : children) {
            // Leaf petiole radius: bigger leaves have thicker stalks
            float child_radius = child->radius;
            if (child->type == NodeType::LEAF) {
                auto* leaf = child->as_leaf();
                if (leaf) {
                    child_radius = std::max(child->radius, leaf->leaf_size * ref_radius);
                }
            }

            float child_cap = has_cap ? sugar_cap(*child, g) : 0.0f;
            float parent_cap = has_cap ? parent_cap_sugar : 0.0f;

            float desired = compute_transport_flow(
                child->chemical(dp.id), chemical(dp.id),
                child_cap, parent_cap,
                child_radius, radius, ref_radius, dp);

            if (std::abs(desired) > 1e-8f) {
                float bm = get_bias_multiplier(child, g);
                flows.push_back({child, desired, child_cap, bm});
            }
        }

        if (flows.empty()) continue;

        // --- Phase 1: children giving to parent (desired < 0) ---
        // Each child gives proportional to |desired|, limited by its supply.
        // Parent receives, limited by its remaining capacity.
        // When headroom is constrained, bias-weighted proportional redistribution
        // determines who gets priority (redistribution, not amplification).
        float parent_val = chemical(dp.id);
        float parent_headroom = has_cap ? std::max(0.0f, parent_cap_sugar - parent_val) : 1e30f;
        float total_inflow = 0.0f;

        for (auto& cf : flows) {
            if (cf.desired >= 0.0f) continue; // skip receivers
            float give = std::min(-cf.desired, cf.child->chemical(dp.id));
            cf.desired = -give; // store actual (negative = giving)
            total_inflow += give;
        }

        // Bias-weighted redistribution when total exceeds parent's capacity
        if (total_inflow > parent_headroom && total_inflow > 1e-8f) {
            float total_weighted = 0.0f;
            for (auto& cf : flows) {
                if (cf.desired >= 0.0f) continue;
                total_weighted += (-cf.desired) * cf.bias_mult;
            }
            if (total_weighted > 1e-8f) {
                for (auto& cf : flows) {
                    if (cf.desired >= 0.0f) continue;
                    float raw_give = -cf.desired;
                    float share = parent_headroom * (raw_give * cf.bias_mult / total_weighted);
                    cf.desired = -std::min(share, raw_give);
                }
            }
            total_inflow = parent_headroom;
        }

        // Apply inflows
        for (auto& cf : flows) {
            if (cf.desired >= 0.0f) continue;
            float give = -cf.desired;
            cf.child->chemical(dp.id) -= give;
            parent_val += give;
        }
        chemical(dp.id) = parent_val;

        // --- Phase 2: parent giving to children (desired > 0) ---
        // Distribute proportionally by bias-adjusted weight. Budget uses raw
        // total (no amplification). Bias only changes who gets what share.
        // If a child fills up, redistribute its remainder among unfilled children.
        float available = chemical(dp.id);

        // Collect receivers with bias-adjusted weights
        struct Receiver {
            Node* child;
            float raw_weight;   // unbiased desired flow
            float weight;       // bias-adjusted for proportional split
            float headroom;     // how much this child can still accept
        };
        std::vector<Receiver> receivers;
        for (auto& cf : flows) {
            if (cf.desired <= 0.0f) continue;
            float headroom = has_cap
                ? std::max(0.0f, cf.child_cap - cf.child->chemical(dp.id))
                : 1e30f;
            if (headroom > 1e-8f) {
                receivers.push_back({cf.child, cf.desired, cf.desired * cf.bias_mult, headroom});
            }
        }

        // Iteratively distribute — redistribute remainder when a child fills
        while (!receivers.empty() && available > 1e-8f) {
            float raw_total = 0.0f;
            float bias_total = 0.0f;
            for (const auto& r : receivers) {
                raw_total += r.raw_weight;
                bias_total += r.weight;
            }
            if (bias_total < 1e-8f) break;

            float to_distribute = std::min(available, raw_total);
            bool any_filled = false;

            for (auto& r : receivers) {
                float share = to_distribute * (r.weight / bias_total);
                float actual = std::min(share, r.headroom);
                r.child->chemical(dp.id) += actual;
                available -= actual;
                r.headroom -= actual;
                if (r.headroom < 1e-8f) any_filled = true;
            }

            if (!any_filled) break; // everyone got their share, done

            // Remove filled children and retry with remainder
            receivers.erase(
                std::remove_if(receivers.begin(), receivers.end(),
                    [](const Receiver& r) { return r.headroom < 1e-8f; }),
                receivers.end());
        }

        chemical(dp.id) = available;
    }
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
