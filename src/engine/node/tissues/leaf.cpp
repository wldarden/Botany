#include "engine/node/tissues/leaf.h"
#include "engine/node/meristems/helpers.h"
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

void LeafNode::update_tissue(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    produce_gibberellin(g);
    float net_sugar = photosynthesize(plant, g, world);
    transpire(g, world);
    // Ethylene-triggered senescence: high ethylene overrides carbon balance
    if (senescence_ticks == 0 && chemical(ChemicalID::Ethylene) >= g.ethylene_abscission_threshold) {
        senescence_ticks = 1;
    }
    check_carbon_balance(g, world, net_sugar);
    if (advance_senescence(plant, g)) return;
    phototropism(g, world);
    expand(g, world);
}

void LeafNode::produce_gibberellin(const Genome& g) {
    if (age < g.ga_leaf_age_max && leaf_size > 1e-6f && senescence_ticks == 0) {
        chemical(ChemicalID::Gibberellin) += leaf_size * g.ga_production_rate;
    }
}

void LeafNode::transpire(const Genome& g, const WorldParams& world) {
    float leaf_area = leaf_size * leaf_size;
    float wcap = water_cap(*this, g);
    float stomatal = wcap > 1e-6f
        ? std::clamp(chemical(ChemicalID::Water) / wcap, 0.2f, 1.0f)
        : 1.0f;
    float transpired = g.transpiration_rate * leaf_area * light_exposure * stomatal;
    chemical(ChemicalID::Water) = std::max(0.0f, chemical(ChemicalID::Water) - transpired);
}

void LeafNode::check_carbon_balance(const Genome& g, const WorldParams& world, float net_sugar) {
    if (senescence_ticks != 0 || age < g.min_leaf_age_before_abscission) return;

    float cost = maintenance_cost(world);
    if (net_sugar < cost) deficit_ticks++;
    else                  deficit_ticks = 0;

    if (deficit_ticks >= g.leaf_abscission_ticks) {
        senescence_ticks = 1;
    }
}

bool LeafNode::advance_senescence(Plant& plant, const Genome& g) {
    if (senescence_ticks == 0) return false;
    senescence_ticks++;
    if (senescence_ticks >= g.senescence_duration) {
        die(plant);
        return true;
    }
    return false;
}

float LeafNode::photosynthesize(Plant& plant, const Genome& g, const WorldParams& world) {
    if (leaf_size <= 1e-6f || senescence_ticks != 0) return 0.0f;

    float cap = sugar_cap(*this, g);
    float sugar_before = chemical(ChemicalID::Sugar);
    if (sugar_before >= cap) return 0.0f;

    float angle_efficiency = 1.0f;
    float facing_len = glm::length(facing);
    if (facing_len > 1e-4f) {
        glm::vec3 leaf_normal = facing / facing_len;
        angle_efficiency = std::max(0.0f, leaf_normal.y);
    }

    float leaf_area = leaf_size * leaf_size;

    // Stomatal conductance: water deficit partially closes stomata
    float wcap = water_cap(*this, g);
    float stomatal = wcap > 1e-6f
        ? std::clamp(chemical(ChemicalID::Water) / wcap, 0.2f, 1.0f)
        : 1.0f;

    float sugar_produced = light_exposure * angle_efficiency
           * world.light_level * leaf_area
           * world.sugar_production_rate
           * stomatal;
    chemical(ChemicalID::Sugar) += sugar_produced;
    chemical(ChemicalID::Sugar) = std::min(chemical(ChemicalID::Sugar), cap);
    float delta = chemical(ChemicalID::Sugar) - sugar_before;

    // Photosynthesis water cost: small deduction proportional to sugar produced
    float water_cost = sugar_produced * g.photosynthesis_water_ratio;
    chemical(ChemicalID::Water) = std::max(0.0f, chemical(ChemicalID::Water) - water_cost);

    if (delta > 0.0f) plant.add_sugar_produced(delta);
    return delta;
}

void LeafNode::phototropism(const Genome& g, const WorldParams& world) {
    if (g.leaf_phototropism_rate <= 1e-8f) return;

    float facing_len = glm::length(facing);
    if (facing_len <= 1e-4f) return;

    glm::vec3 dir = facing / facing_len;
    // sun_direction points toward ground; negate to get the toward-light vector leaves track.
    glm::vec3 light_dir = -world.sun_direction;
    float cos_angle = glm::dot(dir, light_dir);
    if (cos_angle >= 0.999f) return;

    glm::vec3 axis = glm::cross(dir, light_dir);
    float axis_len = glm::length(axis);
    if (axis_len <= 1e-6f) return;
    axis /= axis_len;

    float angle_to_up = std::acos(std::min(cos_angle, 1.0f));
    float turn = std::min(g.leaf_phototropism_rate, angle_to_up);

    float cost = turn * world.sugar_cost_phototropism;
    if (chemical(ChemicalID::Sugar) < cost) return;

    chemical(ChemicalID::Sugar) -= cost;
    float c = std::cos(turn);
    float s = std::sin(turn);
    glm::vec3 new_dir = dir * c
        + glm::cross(axis, dir) * s
        + axis * glm::dot(axis, dir) * (1.0f - c);
    facing = glm::normalize(new_dir);
}

void LeafNode::expand(const Genome& g, const WorldParams& world) {
    if (leaf_size >= g.max_leaf_size) return;

    float auxin_boost = meristem_helpers::auxin_growth_factor(
        chemical(ChemicalID::Auxin), g.leaf_auxin_max_boost, g.leaf_auxin_half_saturation);
    float max_growth = g.leaf_growth_rate * auxin_boost;
    float remaining = g.max_leaf_size - leaf_size;
    float growth = std::min(max_growth, remaining);

    // Turgor pressure gates leaf expansion — cells can't inflate without water
    float water_gf = meristem_helpers::turgor_fraction(chemical(ChemicalID::Water), water_cap(*this, g));
    if (water_gf < 1e-6f) return;
    growth *= water_gf;

    float cost = growth * world.sugar_cost_leaf_growth;
    if (chemical(ChemicalID::Sugar) < cost) {
        growth = chemical(ChemicalID::Sugar) / world.sugar_cost_leaf_growth;
        cost = chemical(ChemicalID::Sugar);
    }
    if (growth < 1e-7f) return;

    leaf_size += growth;
    chemical(ChemicalID::Sugar) -= cost;

    // Auxin production: growing leaves produce auxin proportional to growth rate.
    // No growth (full size, stressed, starved) → zero auxin.
    float growth_fraction = growth / g.leaf_growth_rate;
    float produced = growth_fraction * g.leaf_growth_auxin_multiplier * g.leaf_auxin_baseline;
    chemical(ChemicalID::Auxin) += produced;
    tick_auxin_produced += produced;

    // Extend petiole proportionally as leaf grows
    // Target offset length = base_offset + petiole_length * (leaf_size / max_leaf_size)
    float olen = glm::length(offset);
    if (olen > 1e-4f) {
        float petiole_growth = g.leaf_petiole_length * (growth / g.max_leaf_size);
        offset += (offset / olen) * petiole_growth;
    }
}

void LeafNode::compute_growth_reserve(const Genome& g, const WorldParams& world) {
    sugar_reserved_for_growth = 0.0f;
    // Senescing leaves don't expand; full-grown leaves don't need a reserve.
    if (senescence_ticks != 0 || leaf_size >= g.max_leaf_size) return;

    // Mirror expand(): compute expected sugar cost this tick before vascular runs.
    float auxin_boost = meristem_helpers::auxin_growth_factor(
        chemical(ChemicalID::Auxin), g.leaf_auxin_max_boost, g.leaf_auxin_half_saturation);
    float max_growth = g.leaf_growth_rate * auxin_boost;
    float remaining  = g.max_leaf_size - leaf_size;
    float growth     = std::min(max_growth, remaining);

    float water_gf = meristem_helpers::turgor_fraction(
        chemical(ChemicalID::Water), water_cap(*this, g));
    growth *= water_gf;

    float cost = growth * world.sugar_cost_leaf_growth;
    // Reserve only what we already have — can't protect sugar we don't possess.
    sugar_reserved_for_growth = std::min(cost, chemical(ChemicalID::Sugar));
}

float LeafNode::maintenance_cost(const WorldParams& world) const {
    return world.sugar_maintenance_leaf * leaf_size * leaf_size;
}

} // namespace botany
