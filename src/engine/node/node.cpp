#include "engine/node/node.h"
#include "engine/plant.h"
#include "engine/sugar.h"
#include "engine/world_params.h"
#include "engine/chemical/chemical_registry.h"
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
    age++;
    const Genome& g = plant.genome();

    // Maintenance sugar consumption
    float cost = maintenance_cost(g);
    chemical(ChemicalID::Sugar) = std::max(0.0f, chemical(ChemicalID::Sugar) - cost);

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

    // Grow before transport: leaves use sugar for expansion first, then
    // export the surplus. Otherwise transport drains leaves dry and they
    // never mature — the "tragedy of the commons" problem.
    grow(plant, world);

    // --- Mass & stress computation ---
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

    // Accumulate subtree mass from direct children (their values are one tick stale)
    total_mass = self_mass;
    mass_moment = self_mass * position;
    for (const Node* child : children) {
        total_mass += child->total_mass;
        mass_moment += child->mass_moment;
    }

    // Stress computation (above-ground only, skip if near ground)
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

    // Stress hormone production (above-ground only)
    if (!is_underground && stress > 0.0f) {
        chemical(ChemicalID::Stress) += stress * g.stress_hormone_production_rate;
    }

    // --- Droop and break (above-ground stems only) ---
    if (type == NodeType::STEM && !is_underground && stress > 0.0f) {
        float break_stress = g.wood_density * world.break_strength_factor;
        float droop_threshold = break_stress * g.wood_flexibility;

        if (stress >= break_stress) {
            // Branch snaps — remove this node and entire subtree
            die(plant);
            return;  // node is dead, skip transport
        }

        if (stress > droop_threshold) {
            // Branch droops toward gravity
            float excess = (stress - droop_threshold) / (break_stress - droop_threshold);
            float droop_angle = std::min(excess * world.droop_rate, world.droop_rate);
            float len = glm::length(offset);
            if (len > 1e-4f) {
                glm::vec3 dir = offset / len;
                glm::vec3 down(0.0f, -1.0f, 0.0f);
                // Rotate dir toward down by droop_angle using Rodrigues' rotation
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

    transport_chemicals(g);

    // Re-clamp sugar after transport — diffusion can overfill small nodes
    chemical(ChemicalID::Sugar) = std::min(chemical(ChemicalID::Sugar), cap);
}

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

float Node::maintenance_cost(const Genome& /*g*/) const {
    return 0.0f;
}

void Node::transport_chemicals(const Genome& g) {
    if (parent) {
        for (const auto& dp : diffusion_params(g)) {
            if (dp.id == ChemicalID::Sugar) {
                // Sugar transport: cap-aware to prevent annihilation.
                // Limit flow by receiver's headroom so sugar isn't destroyed by cap clamp.
                float my_cap = sugar_cap(*this, g);
                float parent_cap = sugar_cap(*parent, g);
                float flow = (chemical(dp.id) - parent->chemical(dp.id)) * dp.diffusion_rate;
                if (flow > 0.0f) {
                    float headroom = std::max(0.0f, parent_cap - parent->chemical(dp.id));
                    flow = std::min({flow, chemical(dp.id), headroom});
                } else {
                    float headroom = std::max(0.0f, my_cap - chemical(dp.id));
                    flow = std::max({flow, -parent->chemical(dp.id), -headroom});
                }
                chemical(dp.id) -= flow;
                parent->chemical(dp.id) += flow;
            } else {
                diffuse(chemical(dp.id), parent->chemical(dp.id), dp.diffusion_rate);
                chemical(dp.id) *= (1.0f - dp.decay_rate);
            }
        }
    } else {
        // Seed: just decay
        for (const auto& dp : diffusion_params(g)) {
            chemical(dp.id) *= (1.0f - dp.decay_rate);
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
