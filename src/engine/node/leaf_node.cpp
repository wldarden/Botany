#include "engine/node/leaf_node.h"
#include "engine/plant.h"
#include "engine/sugar.h"
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

    photosynthesize(g, world);
    phototropism(g, world);
    grow(g, world);
}

void LeafNode::photosynthesize(const Genome& g, const WorldParams& world) {
    if (leaf_size <= 1e-6f || senescence_ticks != 0) return;

    float cap = sugar_cap(*this, g);
    if (sugar >= cap) return;

    float angle_efficiency = 1.0f;
    float offset_len = glm::length(offset);
    if (offset_len > 1e-4f) {
        glm::vec3 leaf_normal = offset / offset_len;
        angle_efficiency = std::max(0.0f, leaf_normal.y);
    }

    sugar += light_exposure * angle_efficiency
           * world.light_level * leaf_size
           * g.sugar_production_rate;
    sugar = std::min(sugar, cap);
}

void LeafNode::phototropism(const Genome& g, const WorldParams& world) {
    if (g.leaf_phototropism_rate <= 1e-8f) return;

    float offset_len = glm::length(offset);
    if (offset_len <= 1e-4f) return;

    glm::vec3 dir = offset / offset_len;
    glm::vec3 up(0.0f, 1.0f, 0.0f);
    float cos_angle = glm::dot(dir, up);
    if (cos_angle >= 0.999f) return;

    glm::vec3 axis = glm::cross(dir, up);
    float axis_len = glm::length(axis);
    if (axis_len <= 1e-6f) return;
    axis /= axis_len;

    float angle_to_up = std::acos(std::min(cos_angle, 1.0f));
    float turn = std::min(g.leaf_phototropism_rate, angle_to_up);

    float cost = turn * world.sugar_cost_phototropism;
    if (sugar < cost) return;

    sugar -= cost;
    float c = std::cos(turn);
    float s = std::sin(turn);
    glm::vec3 new_dir = dir * c
        + glm::cross(axis, dir) * s
        + axis * glm::dot(axis, dir) * (1.0f - c);
    offset = glm::normalize(new_dir) * offset_len;
}

void LeafNode::grow(const Genome& g, const WorldParams& world) {
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

float LeafNode::maintenance_cost(const Genome& g) const {
    return g.sugar_maintenance_leaf * leaf_size * leaf_size;
}

} // namespace botany
