#include "engine/sugar.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include "engine/node/node.h"
#include "engine/node/leaf_node.h"
#include <algorithm>
#include <cmath>
#include <glm/geometric.hpp>

namespace botany {

float sugar_cap(const Node& node, const Genome& g) {
    float length = std::max(glm::length(node.offset), 0.01f);

    switch (node.type) {
        case NodeType::LEAF: {
            float ls = node.as_leaf()->leaf_size;
            float area = ls * ls;
            return std::max(g.sugar_cap_minimum, area * g.sugar_storage_density_leaf);
        }
        case NodeType::STEM:
        case NodeType::ROOT: {
            float volume = 3.14159f * node.radius * node.radius * length;
            float cap = volume * g.sugar_storage_density_wood;
            if (!node.parent) {
                cap = std::max(cap, g.seed_sugar);
            }
            return std::max(g.sugar_cap_minimum, cap);
        }
        case NodeType::SHOOT_APICAL:
        case NodeType::SHOOT_AXILLARY:
        case NodeType::ROOT_APICAL:
        case NodeType::ROOT_AXILLARY:
            return g.sugar_cap_meristem;
    }
    return g.sugar_cap_minimum;
}

// Light exposure: engine/light.cpp (world-level, called from Engine::tick)
// Sugar consumption + starvation + death: Node::tick()

void grow_leaves(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();
    glm::vec3 up(0.0f, 1.0f, 0.0f);

    plant.for_each_node_mut([&](Node& node) {
        auto* leaf = node.as_leaf();
        if (!leaf) return;

        // Phototropism: rotate leaf offset toward sky
        if (g.leaf_phototropism_rate > 1e-8f) {
            float offset_len = glm::length(node.offset);
            if (offset_len > 1e-4f) {
                glm::vec3 dir = node.offset / offset_len;
                float cos_angle = glm::dot(dir, up);
                if (cos_angle < 0.999f) { // not already pointing up
                    glm::vec3 axis = glm::cross(dir, up);
                    float axis_len = glm::length(axis);
                    if (axis_len > 1e-6f) {
                        axis /= axis_len;
                        float angle_to_up = std::acos(std::min(cos_angle, 1.0f));
                        float turn = std::min(g.leaf_phototropism_rate, angle_to_up);

                        float cost = turn * world.sugar_cost_phototropism;
                        if (node.chemical(ChemicalID::Sugar) >= cost) {
                            node.chemical(ChemicalID::Sugar) -= cost;
                            float c = std::cos(turn);
                            float s = std::sin(turn);
                            glm::vec3 new_dir = dir * c
                                + glm::cross(axis, dir) * s
                                + axis * glm::dot(axis, dir) * (1.0f - c);
                            node.offset = glm::normalize(new_dir) * offset_len;
                        }
                    }
                }
            }
        }

        // Size growth
        if (leaf->leaf_size >= g.max_leaf_size) return;

        float max_growth = g.leaf_growth_rate;
        float remaining = g.max_leaf_size - leaf->leaf_size;
        float growth = std::min(max_growth, remaining);

        float cost = growth * world.sugar_cost_leaf_growth;
        if (node.chemical(ChemicalID::Sugar) < cost) {
            growth = node.chemical(ChemicalID::Sugar) / world.sugar_cost_leaf_growth;
            cost = node.chemical(ChemicalID::Sugar);
        }
        if (growth < 1e-7f) return;

        leaf->leaf_size += growth;
        node.chemical(ChemicalID::Sugar) -= cost;
    });
}

void transport_sugar(Plant& plant, const WorldParams& world) {
    // Light exposure: computed by Engine before plant ticks (engine/light.cpp)
    // Sugar production: LeafNode::tick()
    // Sugar diffusion: Node::transport_chemicals()
    // Sugar consumption + starvation death: Node::tick()
}

} // namespace botany
