#include "engine/node/root_node.h"
#include "engine/node/meristems/helpers.h"
#include "engine/plant.h"
#include "engine/sugar.h"
#include "engine/world_params.h"
#include <algorithm>
#include <glm/geometric.hpp>

namespace botany {

RootNode::RootNode(uint32_t id, glm::vec3 position, float radius)
    : Node(id, NodeType::ROOT, position, radius)
{}

void RootNode::update_tissue(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();
    absorb_water(g, world);
    thicken(g, world);
    elongate(g, world);
}

void RootNode::absorb_water(const Genome& g, const WorldParams& world) {
    float length = std::max(glm::length(offset), 0.01f);
    float surface_area = 2.0f * 3.14159f * radius * length;
    float cap = water_cap(*this, g);
    float fill_fraction = (cap > 1e-6f) ? chemical(ChemicalID::Water) / cap : 1.0f;
    float gradient = std::max(0.0f, world.soil_moisture - fill_fraction);
    float absorbed = g.water_absorption_rate * surface_area * gradient;
    chemical(ChemicalID::Water) = std::min(chemical(ChemicalID::Water) + absorbed, cap);
}

void RootNode::thicken(const Genome& g, const WorldParams& world) {
    // Same vascular-driven model as StemNode. Root connections accumulate
    // structural_flow_bias from sugar and cytokinin transport (not auxin —
    // real root polar auxin transport governs patterning/gravitropism, not
    // cambial signaling). Well-used root connections thicken; lateral roots
    // with low flux stay thin.
    float bias = get_parent_structural_bias();
    if (bias < 1e-6f) return;

    float max_cost = g.cambium_responsiveness * bias * world.sugar_cost_stem_growth;
    float sugar_gf = (max_cost > 1e-6f) ? std::min(chemical(ChemicalID::Sugar) / max_cost, 1.0f) : 1.0f;
    if (sugar_gf <= 1e-6f) return;

    float actual_rate = g.cambium_responsiveness * bias * sugar_gf;
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_stem_growth;
    radius += actual_rate;
}

void RootNode::elongate(const Genome& g, const WorldParams& world) {
    if (!parent) return;
    if (age >= g.root_internode_maturation_ticks) return;
    if (g.root_internode_elongation_rate <= 1e-8f) return;

    float ga_boost = 1.0f + chemical(ChemicalID::Gibberellin) * g.ga_elongation_sensitivity;
    float eth_inhibit = std::max(0.0f, 1.0f - chemical(ChemicalID::Ethylene) * g.ethylene_elongation_inhibition);
    float auxin_boost = meristem_helpers::auxin_growth_factor(
        chemical(ChemicalID::Auxin), g.root_auxin_max_boost, g.root_auxin_half_saturation);
    float effective_rate = g.root_internode_elongation_rate * ga_boost * eth_inhibit * auxin_boost;

    float max_len = g.max_internode_length * (1.0f + chemical(ChemicalID::Gibberellin) * g.ga_length_sensitivity);
    float current_len = glm::length(offset);
    if (current_len >= max_len) return;

    // Elongation requires sugar (construction) and water (turgor pressure).
    float max_cost = effective_rate * world.sugar_cost_stem_growth;
    float sugar_gf = (max_cost > 1e-6f) ? std::min(chemical(ChemicalID::Sugar) / max_cost, 1.0f) : 1.0f;
    if (sugar_gf <= 1e-6f) return;
    float water_gf = meristem_helpers::turgor_fraction(chemical(ChemicalID::Water), water_cap(*this, g));
    if (water_gf <= 1e-6f) return;

    float actual_rate = effective_rate * sugar_gf * water_gf;
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
