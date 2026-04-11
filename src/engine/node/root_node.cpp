#include "engine/node/root_node.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include "engine/node/meristems/helpers.h"
#include <algorithm>
#include <glm/geometric.hpp>

namespace botany {

using meristem_helpers::sugar_growth_fraction;

RootNode::RootNode(uint32_t id, glm::vec3 position, float radius)
    : Node(id, NodeType::ROOT, position, radius)
{}

void RootNode::tick(Plant& plant, const WorldParams& world) {
    Node::tick(plant, world);
    const Genome& g = plant.genome();

    thicken(g, world);
    elongate(g, world);
}

void RootNode::thicken(const Genome& g, const WorldParams& world) {
    float max_cost = g.thickening_rate * world.sugar_cost_thickening;
    float gf = sugar_growth_fraction(sugar, g.sugar_save_root, max_cost);
    if (gf <= 1e-6f) return;

    float actual_rate = g.thickening_rate * gf;
    sugar -= actual_rate * world.sugar_cost_thickening;
    radius += actual_rate;
}

void RootNode::elongate(const Genome& g, const WorldParams& world) {
    if (!parent) return;
    if (age >= g.root_internode_maturation_ticks) return;
    if (g.root_internode_elongation_rate <= 1e-8f) return;

    float ga_boost = 1.0f + gibberellin * g.ga_elongation_sensitivity;
    float eth_inhibit = std::max(0.0f, 1.0f - ethylene * g.ethylene_elongation_inhibition);
    float effective_rate = g.root_internode_elongation_rate * ga_boost * eth_inhibit;

    float max_len = g.root_max_internode_length * (1.0f + gibberellin * g.ga_length_sensitivity);
    float current_len = glm::length(offset);
    if (current_len >= max_len) return;

    float max_cost = effective_rate * world.sugar_cost_elongation;
    float gf = sugar_growth_fraction(sugar, g.sugar_save_root, max_cost);
    if (gf <= 1e-6f) return;

    float actual_rate = effective_rate * gf;
    sugar -= actual_rate * world.sugar_cost_elongation;
    if (current_len > 1e-4f) {
        offset += (offset / current_len) * actual_rate;
    }
}

float RootNode::maintenance_cost(const Genome& g) const {
    float length = std::max(glm::length(offset), 0.01f);
    float volume = 3.14159f * radius * radius * length;
    return g.sugar_maintenance_root * volume;
}

} // namespace botany
