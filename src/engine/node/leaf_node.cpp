#include "engine/node/leaf_node.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include <algorithm>
#include <cmath>
#include <glm/geometric.hpp>

namespace botany {

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
