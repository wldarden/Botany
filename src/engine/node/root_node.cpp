#include "engine/node/root_node.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include <algorithm>
#include <glm/geometric.hpp>

namespace botany {

RootNode::RootNode(uint32_t id, glm::vec3 position, float radius)
    : Node(id, NodeType::ROOT, position, radius)
{}

void RootNode::grow(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();
    thicken(g, world);
    elongate(g, world);
}

void RootNode::thicken(const Genome& g, const WorldParams& world) {
    float effective_rate = g.thickening_rate;

    // Auxin-driven cambial growth: auxin flows basipetally from the canopy
    // down into the root system — main root gets the most, thickens the most.
    float auxin_gf = std::min(chemical(ChemicalID::Auxin) / std::max(g.auxin_threshold, 1e-6f), 1.0f);
    effective_rate *= auxin_gf;

    // Sugar funds the growth
    float max_cost = effective_rate * world.sugar_cost_stem_growth;
    float sugar_gf = (max_cost > 1e-6f) ? std::min(chemical(ChemicalID::Sugar) / max_cost, 1.0f) : 1.0f;
    if (sugar_gf <= 1e-6f) return;

    float actual_rate = effective_rate * sugar_gf;
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_stem_growth;
    radius += actual_rate;
}

void RootNode::elongate(const Genome& g, const WorldParams& world) {
    if (!parent) return;
    if (age >= g.root_internode_maturation_ticks) return;
    if (g.root_internode_elongation_rate <= 1e-8f) return;

    float ga_boost = 1.0f + chemical(ChemicalID::Gibberellin) * g.ga_elongation_sensitivity;
    float eth_inhibit = std::max(0.0f, 1.0f - chemical(ChemicalID::Ethylene) * g.ethylene_elongation_inhibition);
    float effective_rate = g.root_internode_elongation_rate * ga_boost * eth_inhibit;

    float max_len = g.max_internode_length * (1.0f + chemical(ChemicalID::Gibberellin) * g.ga_length_sensitivity);
    float current_len = glm::length(offset);
    if (current_len >= max_len) return;

    // Elongation is sugar-gated only (not cytokinin) — same reasoning as stems
    float max_cost = effective_rate * world.sugar_cost_stem_growth;
    float sugar_gf = (max_cost > 1e-6f) ? std::min(chemical(ChemicalID::Sugar) / max_cost, 1.0f) : 1.0f;
    if (sugar_gf <= 1e-6f) return;

    float actual_rate = effective_rate * sugar_gf;
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_stem_growth;
    if (current_len > 1e-4f) {
        offset += (offset / current_len) * actual_rate;
    }
}

float RootNode::maintenance_cost(const WorldParams& world) const {
    float length = std::max(glm::length(offset), 0.01f);
    float volume = 3.14159f * radius * radius * length;
    return world.sugar_maintenance_root * volume;
}

} // namespace botany
