#include "engine/node/stem_node.h"
#include "engine/node/meristems/helpers.h"
#include "engine/plant.h"
#include "engine/sugar.h"
#include "engine/world_params.h"
#include <algorithm>
#include <glm/geometric.hpp>

namespace botany {

StemNode::StemNode(uint32_t id, glm::vec3 position, float radius)
    : Node(id, NodeType::STEM, position, radius)
{}

void StemNode::update_tissue(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();
    thicken(g, world);
    elongate(g, world);
}

void StemNode::thicken(const Genome& g, const WorldParams& world) {
    // Cambium activity is driven by PIN canalization history (auxin_flow_bias), not age.
    // Connections where auxin has flowed repeatedly have higher PIN saturation and thus
    // stronger bias; those that haven't stay thin. This creates the self-reinforcing
    // loop: main axis carries most auxin flux → highest bias → fastest thickening →
    // widest pipe → more vascular flow → more flux. Lateral branches stay thin.
    float bias;
    if (!parent) {
        // Seed node: use the max of its children's auxin_flow_bias entries.
        bias = 0.0f;
        for (auto& [child, b] : auxin_flow_bias) {
            bias = std::max(bias, b);
        }
    } else {
        // Get the bias this parent has recorded for this child connection.
        float b = 0.0f;
        auto it = parent->auxin_flow_bias.find(this);
        if (it != parent->auxin_flow_bias.end()) b = it->second;
        bias = b;
    }
    if (bias < 1e-6f) return;

    // Dense wood costs proportionally more sugar per unit of radial growth.
    float density_scale = g.wood_density / world.reference_wood_density;

    // Sugar affordability: can we pay for the full bias-driven rate this tick?
    float max_cost = g.cambium_responsiveness * bias * world.sugar_cost_stem_growth * density_scale;
    float sugar_gf = (max_cost > 1e-6f) ? std::min(chemical(ChemicalID::Sugar) / max_cost, 1.0f) : 1.0f;
    if (sugar_gf <= 1e-6f) return;

    // Thigmomorphogenesis: mechanical stress boosts cambial activity (KEEP — real biology).
    float stress_boost = 1.0f + chemical(ChemicalID::Stress) * g.stress_thickening_boost;
    float actual_rate = g.cambium_responsiveness * bias * sugar_gf * stress_boost;
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_stem_growth * density_scale;
    radius += actual_rate;
}

void StemNode::elongate(const Genome& g, const WorldParams& world) {
    if (!parent) return;
    if (age >= g.internode_maturation_ticks) return;
    if (g.internode_elongation_rate <= 1e-8f) return;

    float ga_boost = 1.0f + chemical(ChemicalID::Gibberellin) * g.ga_elongation_sensitivity;
    float eth_inhibit = std::max(0.0f, 1.0f - chemical(ChemicalID::Ethylene) * g.ethylene_elongation_inhibition);
    float stress_inhibit = std::max(0.0f, 1.0f - chemical(ChemicalID::Stress) * g.stress_elongation_inhibition);
    float auxin_boost = meristem_helpers::auxin_growth_factor(
        chemical(ChemicalID::Auxin), g.stem_auxin_max_boost, g.stem_auxin_half_saturation);
    float effective_rate = g.internode_elongation_rate * ga_boost * eth_inhibit * stress_inhibit * auxin_boost;

    float max_len = g.max_internode_length * (1.0f + chemical(ChemicalID::Gibberellin) * g.ga_length_sensitivity);
    float current_len = glm::length(offset);
    if (current_len >= max_len) return;

    // Elongation requires sugar (construction) and water (turgor pressure).
    float density_scale = g.wood_density / world.reference_wood_density;
    float max_cost = effective_rate * world.sugar_cost_stem_growth * density_scale;
    float sugar_gf = (max_cost > 1e-6f) ? std::min(chemical(ChemicalID::Sugar) / max_cost, 1.0f) : 1.0f;
    if (sugar_gf <= 1e-6f) return;
    float water_gf = meristem_helpers::turgor_fraction(chemical(ChemicalID::Water), water_cap(*this, g));
    if (water_gf <= 1e-6f) return;

    float actual_rate = effective_rate * sugar_gf * water_gf;
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_stem_growth * density_scale;
    if (current_len > 1e-4f) {
        offset += (offset / current_len) * actual_rate;
    }
}

float StemNode::maintenance_cost(const WorldParams& world) const {
    float length = std::max(glm::length(offset), 0.01f);
    float volume = 3.14159f * radius * radius * length;
    return world.sugar_maintenance_stem * volume;
}

} // namespace botany
