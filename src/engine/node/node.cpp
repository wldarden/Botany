#include "engine/node/node.h"
#include "engine/plant.h"
#include "engine/sugar.h"
#include "engine/world_params.h"
#include "engine/chemical/chemical_registry.h"
#include "engine/perf_log.h"
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/node/leaf_node.h"
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
            return;
        }
    }
}

void Node::tick(Plant& plant, const WorldParams& world) {
    PerfStats* perf = plant.perf();

    {   // --- Position ---
        ScopedTimer t(perf ? perf->node_position_ms : ScopedTimer::dummy);
        position = parent ? parent->position + offset : offset;
        if (!std::isfinite(position.x) || !std::isfinite(position.y) || !std::isfinite(position.z)) {
            position = parent ? parent->position : glm::vec3(0.0f);
        }
    }

    age++;
    const Genome& g = plant.genome();

    {   // --- Maintenance ---
        ScopedTimer t(perf ? perf->node_maintenance_ms : ScopedTimer::dummy);
        float cost = maintenance_cost(world);
        chemical(ChemicalID::Sugar) = std::max(0.0f, chemical(ChemicalID::Sugar) - cost);
    }

    // Cap clamp
    float cap = sugar_cap(*this, g);
    chemical(ChemicalID::Sugar) = std::min(chemical(ChemicalID::Sugar), cap);

    // Starvation tracking + death
    if (chemical(ChemicalID::Sugar) <= 0.0f) starvation_ticks++;
    else starvation_ticks = 0;

    if (starvation_ticks >= world.starvation_ticks_max && parent) {
        die(plant);
        return;
    }

    {   // --- Produce ---
        produce(plant, world);
    }

    {   // --- Grow ---
        ScopedTimer t(perf ? perf->node_grow_ms : ScopedTimer::dummy);
        grow(plant, world);
    }

    {   // --- Mass & stress ---
        ScopedTimer t(perf ? perf->node_mass_stress_ms : ScopedTimer::dummy);
        float self_mass = 0.0f;
        bool is_underground = (type == NodeType::ROOT || type == NodeType::ROOT_APICAL || type == NodeType::ROOT_AXILLARY);

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

        stress = 0.0f;
        if (!is_underground && position.y > world.ground_support_height) {
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
        }

        if (!is_underground && stress > 0.0f) {
            float break_stress = g.wood_density * world.break_strength_factor;
            float stress_ratio = stress / break_stress;  // 0 = no load, 1 = breaking
            if (stress_ratio > g.stress_hormone_threshold) {
                float excess = (stress_ratio - g.stress_hormone_threshold)
                             / (1.0f - g.stress_hormone_threshold);  // normalize 0-1
                chemical(ChemicalID::Stress) += excess * g.stress_hormone_production_rate;
            }
        }
    }

    {   // --- Droop and break ---
        ScopedTimer t(perf ? perf->node_droop_break_ms : ScopedTimer::dummy);
        bool is_underground = (type == NodeType::ROOT || type == NodeType::ROOT_APICAL || type == NodeType::ROOT_AXILLARY);
        if (type == NodeType::STEM && !is_underground && stress > 0.0f) {
            float break_stress = g.wood_density * world.break_strength_factor;
            float droop_threshold = break_stress * g.wood_flexibility;

            // Ground-anchored stems can't snap — the root system holds them.
            // Only break stems above ground support height with a parent that also has a parent (not trunk base).
            bool can_break = position.y > world.ground_support_height && parent && parent->parent;
            if (can_break && stress >= break_stress) {
                die(plant);
                return;
            }

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
        }
    }

    {   // --- Transport ---
        ScopedTimer t(perf ? perf->node_transport_ms : ScopedTimer::dummy);
        transport_chemicals(g);
    }

    // Re-clamp sugar after transport — diffusion can overfill small nodes
    chemical(ChemicalID::Sugar) = std::min(chemical(ChemicalID::Sugar), cap);
}

void Node::produce(Plant& /*plant*/, const WorldParams& /*world*/) {}

void Node::grow(Plant& /*plant*/, const WorldParams& /*world*/) {}

void Node::die(Plant& plant) {
    // Detach from parent
    if (parent) {
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

void Node::transport_chemicals(const Genome& g) {
    if (parent) {
        float ref_radius = g.initial_radius;

        // Leaves use an effective petiole radius for transport.
        // This guarantees they can export their full sugar production each tick.
        // Petiole scales with leaf size — bigger leaves have thicker stalks.
        float effective_radius = radius;
        if (type == NodeType::LEAF) {
            auto* leaf = as_leaf();
            if (leaf) {
                effective_radius = std::max(radius, leaf->leaf_size * ref_radius);
            }
        }

        for (const auto& dp : diffusion_params(g)) {
            // Sugar uses volume-based capacity; hormones have no cap (0 = raw values).
            float my_cap = 0.0f, parent_cap = 0.0f;
            if (dp.id == ChemicalID::Sugar) {
                my_cap = sugar_cap(*this, g);
                parent_cap = sugar_cap(*parent, g);
            }

            transport_chemical(
                chemical(dp.id), parent->chemical(dp.id),
                my_cap, parent_cap,
                effective_radius, parent->radius, ref_radius,
                dp);

            // Decay (sugar doesn't decay)
            if (dp.decay_rate > 0.0f) {
                chemical(dp.id) *= (1.0f - dp.decay_rate);
            }
        }
    } else {
        // Seed: just decay
        for (const auto& dp : diffusion_params(g)) {
            if (dp.decay_rate > 0.0f) {
                chemical(dp.id) *= (1.0f - dp.decay_rate);
            }
        }
    }
}

bool Node::is_meristem() const {
    return type == NodeType::SHOOT_APICAL || type == NodeType::SHOOT_AXILLARY
        || type == NodeType::ROOT_APICAL  || type == NodeType::ROOT_AXILLARY;
}

// Downcasting helpers
StemNode*       Node::as_stem()       { return type == NodeType::STEM ? static_cast<StemNode*>(this) : nullptr; }
const StemNode* Node::as_stem() const { return type == NodeType::STEM ? static_cast<const StemNode*>(this) : nullptr; }
RootNode*       Node::as_root()       { return type == NodeType::ROOT ? static_cast<RootNode*>(this) : nullptr; }
const RootNode* Node::as_root() const { return type == NodeType::ROOT ? static_cast<const RootNode*>(this) : nullptr; }
LeafNode*       Node::as_leaf()       { return type == NodeType::LEAF ? static_cast<LeafNode*>(this) : nullptr; }
const LeafNode* Node::as_leaf() const { return type == NodeType::LEAF ? static_cast<const LeafNode*>(this) : nullptr; }

ShootApicalNode*       Node::as_shoot_apical()       { return type == NodeType::SHOOT_APICAL ? static_cast<ShootApicalNode*>(this) : nullptr; }
const ShootApicalNode* Node::as_shoot_apical() const { return type == NodeType::SHOOT_APICAL ? static_cast<const ShootApicalNode*>(this) : nullptr; }
ShootAxillaryNode*       Node::as_shoot_axillary()       { return type == NodeType::SHOOT_AXILLARY ? static_cast<ShootAxillaryNode*>(this) : nullptr; }
const ShootAxillaryNode* Node::as_shoot_axillary() const { return type == NodeType::SHOOT_AXILLARY ? static_cast<const ShootAxillaryNode*>(this) : nullptr; }
RootApicalNode*       Node::as_root_apical()       { return type == NodeType::ROOT_APICAL ? static_cast<RootApicalNode*>(this) : nullptr; }
const RootApicalNode* Node::as_root_apical() const { return type == NodeType::ROOT_APICAL ? static_cast<const RootApicalNode*>(this) : nullptr; }
RootAxillaryNode*       Node::as_root_axillary()       { return type == NodeType::ROOT_AXILLARY ? static_cast<RootAxillaryNode*>(this) : nullptr; }
const RootAxillaryNode* Node::as_root_axillary() const { return type == NodeType::ROOT_AXILLARY ? static_cast<const RootAxillaryNode*>(this) : nullptr; }

} // namespace botany
