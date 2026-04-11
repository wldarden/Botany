#include "engine/node.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include "engine/meristems/helpers.h"
#include <algorithm>
#include <cmath>
#include <glm/geometric.hpp>

namespace botany {

using meristem_helpers::sugar_growth_fraction;

// --- Node (base) ---

Node::Node(uint32_t id, NodeType type, glm::vec3 position, float radius)
    : id(id)
    , parent(nullptr)
    , offset(position)
    , position(position)
    , radius(radius)
    , type(type)
    , age(0)
    , auxin(0.0f)
    , cytokinin(0.0f)
    , meristem(nullptr)
{}

void Node::add_child(Node* child) {
    children.push_back(child);
    child->parent = this;
}

void Node::tick(Plant& /*plant*/, const WorldParams& /*world*/) {
    age++;
}

// Downcasting helpers
StemNode*       Node::as_stem()       { return type == NodeType::STEM ? static_cast<StemNode*>(this) : nullptr; }
const StemNode* Node::as_stem() const { return type == NodeType::STEM ? static_cast<const StemNode*>(this) : nullptr; }
RootNode*       Node::as_root()       { return type == NodeType::ROOT ? static_cast<RootNode*>(this) : nullptr; }
const RootNode* Node::as_root() const { return type == NodeType::ROOT ? static_cast<const RootNode*>(this) : nullptr; }
LeafNode*       Node::as_leaf()       { return type == NodeType::LEAF ? static_cast<LeafNode*>(this) : nullptr; }
const LeafNode* Node::as_leaf() const { return type == NodeType::LEAF ? static_cast<const LeafNode*>(this) : nullptr; }

// --- StemNode ---

StemNode::StemNode(uint32_t id, glm::vec3 position, float radius)
    : Node(id, NodeType::STEM, position, radius)
{}

void StemNode::tick(Plant& plant, const WorldParams& world) {
    Node::tick(plant, world);
    const Genome& g = plant.genome();

    // Secondary growth: interior nodes (no active tip meristem) thicken.
    bool is_active_tip = meristem && meristem->is_tip() && meristem->active;
    if (!is_active_tip) {
        float max_cost = g.thickening_rate * world.sugar_cost_thickening;
        float gf = sugar_growth_fraction(sugar, g.sugar_save_stem, max_cost);
        if (gf > 1e-6f) {
            float actual_rate = g.thickening_rate * gf;
            sugar -= actual_rate * world.sugar_cost_thickening;
            radius += actual_rate;
        }
    }

    // Intercalary growth: young interior nodes elongate their internode.
    if (!meristem && parent) {
        if (age < g.internode_maturation_ticks && g.internode_elongation_rate > 1e-8f) {
            // GA boosts elongation rate
            float ga_boost = 1.0f + gibberellin * g.ga_elongation_sensitivity;
            // Ethylene inhibits elongation
            float eth_inhibit = std::max(0.0f, 1.0f - ethylene * g.ethylene_elongation_inhibition);
            float effective_rate = g.internode_elongation_rate * ga_boost * eth_inhibit;

            // GA-modulated max internode length
            float max_len = g.max_internode_length * (1.0f + gibberellin * g.ga_length_sensitivity);
            float current_len = glm::length(offset);
            if (current_len < max_len) {
                float max_cost = effective_rate * world.sugar_cost_elongation;
                float gf = sugar_growth_fraction(sugar, g.sugar_save_stem, max_cost);
                if (gf > 1e-6f) {
                    float actual_rate = effective_rate * gf;
                    sugar -= actual_rate * world.sugar_cost_elongation;
                    float len = glm::length(offset);
                    if (len > 1e-4f) {
                        offset += (offset / len) * actual_rate;
                    }
                }
            }
        }
    }
}

// --- RootNode ---

RootNode::RootNode(uint32_t id, glm::vec3 position, float radius)
    : Node(id, NodeType::ROOT, position, radius)
{}

void RootNode::tick(Plant& plant, const WorldParams& world) {
    Node::tick(plant, world);
    const Genome& g = plant.genome();

    // Secondary growth: interior nodes thicken.
    bool is_active_tip = meristem && meristem->is_tip() && meristem->active;
    if (!is_active_tip) {
        float max_cost = g.thickening_rate * world.sugar_cost_thickening;
        float gf = sugar_growth_fraction(sugar, g.sugar_save_root, max_cost);
        if (gf > 1e-6f) {
            float actual_rate = g.thickening_rate * gf;
            sugar -= actual_rate * world.sugar_cost_thickening;
            radius += actual_rate;
        }
    }

    // Intercalary growth: young interior nodes elongate.
    if (!meristem && parent) {
        if (age < g.root_internode_maturation_ticks && g.root_internode_elongation_rate > 1e-8f) {
            float ga_boost = 1.0f + gibberellin * g.ga_elongation_sensitivity;
            float eth_inhibit = std::max(0.0f, 1.0f - ethylene * g.ethylene_elongation_inhibition);
            float effective_rate = g.root_internode_elongation_rate * ga_boost * eth_inhibit;

            float max_len = g.root_max_internode_length * (1.0f + gibberellin * g.ga_length_sensitivity);
            float current_len = glm::length(offset);
            if (current_len < max_len) {
                float max_cost = effective_rate * world.sugar_cost_elongation;
                float gf = sugar_growth_fraction(sugar, g.sugar_save_root, max_cost);
                if (gf > 1e-6f) {
                    float actual_rate = effective_rate * gf;
                    sugar -= actual_rate * world.sugar_cost_elongation;
                    float len = glm::length(offset);
                    if (len > 1e-4f) {
                        offset += (offset / len) * actual_rate;
                    }
                }
            }
        }
    }
}

// --- LeafNode ---

LeafNode::LeafNode(uint32_t id, glm::vec3 position, float radius)
    : Node(id, NodeType::LEAF, position, radius)
{}

void LeafNode::tick(Plant& plant, const WorldParams& world) {
    Node::tick(plant, world);
    const Genome& g = plant.genome();
    glm::vec3 up(0.0f, 1.0f, 0.0f);

    // Phototropism: rotate leaf offset toward sky
    if (g.leaf_phototropism_rate > 1e-8f) {
        float offset_len = glm::length(offset);
        if (offset_len > 1e-4f) {
            glm::vec3 dir = offset / offset_len;
            float cos_angle = glm::dot(dir, up);
            if (cos_angle < 0.999f) {
                glm::vec3 axis = glm::cross(dir, up);
                float axis_len = glm::length(axis);
                if (axis_len > 1e-6f) {
                    axis /= axis_len;
                    float angle_to_up = std::acos(std::min(cos_angle, 1.0f));
                    float turn = std::min(g.leaf_phototropism_rate, angle_to_up);

                    float cost = turn * world.sugar_cost_phototropism;
                    if (sugar >= cost) {
                        sugar -= cost;
                        float c = std::cos(turn);
                        float s = std::sin(turn);
                        glm::vec3 new_dir = dir * c
                            + glm::cross(axis, dir) * s
                            + axis * glm::dot(axis, dir) * (1.0f - c);
                        offset = glm::normalize(new_dir) * offset_len;
                    }
                }
            }
        }
    }

    // Size growth toward max_leaf_size
    if (leaf_size >= g.max_leaf_size) return;

    float max_growth = g.leaf_growth_rate;
    float remaining = g.max_leaf_size - leaf_size;
    float growth = std::min(max_growth, remaining);

    float cost = growth * world.sugar_cost_leaf_growth;
    if (sugar < cost) {
        growth = sugar / world.sugar_cost_leaf_growth;
        cost = sugar;
    }
    if (growth < 1e-7f) return;

    leaf_size += growth;
    sugar -= cost;
}

} // namespace botany
