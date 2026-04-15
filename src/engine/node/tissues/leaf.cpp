#include "engine/node/tissues/leaf.h"
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

void LeafNode::tissue_tick(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    // Gibberellin production: young leaves produce GA on themselves
    if (age < g.ga_leaf_age_max && leaf_size > 1e-6f && senescence_ticks == 0) {
        chemical(ChemicalID::Gibberellin) += leaf_size * g.ga_production_rate;
    }

    // Photosynthesis: sugar + cytokinin production
    float sugar_before = chemical(ChemicalID::Sugar);
    photosynthesize(g, world);
    float delta = chemical(ChemicalID::Sugar) - sugar_before;
    if (delta > 0.0f) plant.add_sugar_produced(delta);

    // Carbon-balance abscission: if production < maintenance for too long, senesce.
    // Young leaves are immune — they're still growing and not yet net producers.
    if (senescence_ticks == 0 && age >= g.min_leaf_age_before_abscission) {
        float cost = maintenance_cost(world);
        if (delta < cost) deficit_ticks++;
        else              deficit_ticks = 0;

        if (deficit_ticks >= g.leaf_abscission_ticks) {
            senescence_ticks = 1;
        }
    }

    // Senescence countdown → leaf drop
    if (senescence_ticks > 0) {
        senescence_ticks++;
        if (senescence_ticks >= g.senescence_duration) {
            die(plant);
            return;
        }
    }

    // Growth
    phototropism(g, world);
    grow_size(g, world);
}

void LeafNode::photosynthesize(const Genome& g, const WorldParams& world) {
    if (leaf_size <= 1e-6f || senescence_ticks != 0) return;

    float cap = sugar_cap(*this, g);
    if (chemical(ChemicalID::Sugar) >= cap) return;

    float angle_efficiency = 1.0f;
    float facing_len = glm::length(facing);
    if (facing_len > 1e-4f) {
        glm::vec3 leaf_normal = facing / facing_len;
        angle_efficiency = std::max(0.0f, leaf_normal.y);
    }

    float leaf_area = leaf_size * leaf_size;
    float sugar_produced = light_exposure * angle_efficiency
           * world.light_level * leaf_area
           * world.sugar_production_rate;
    chemical(ChemicalID::Sugar) += sugar_produced;
    chemical(ChemicalID::Sugar) = std::min(chemical(ChemicalID::Sugar), cap);

    // Cytokinin production: proportional to actual photosynthetic output.
    // This is the "I have producing leaves" signal that gates all growth.
    chemical(ChemicalID::Cytokinin) += sugar_produced * g.cytokinin_production_rate;
}

void LeafNode::phototropism(const Genome& g, const WorldParams& world) {
    if (g.leaf_phototropism_rate <= 1e-8f) return;

    float facing_len = glm::length(facing);
    if (facing_len <= 1e-4f) return;

    glm::vec3 dir = facing / facing_len;
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
    if (chemical(ChemicalID::Sugar) < cost) return;

    chemical(ChemicalID::Sugar) -= cost;
    float c = std::cos(turn);
    float s = std::sin(turn);
    glm::vec3 new_dir = dir * c
        + glm::cross(axis, dir) * s
        + axis * glm::dot(axis, dir) * (1.0f - c);
    facing = glm::normalize(new_dir);
}

void LeafNode::grow_size(const Genome& g, const WorldParams& world) {
    if (leaf_size >= g.max_leaf_size) return;

    float max_growth = g.leaf_growth_rate;
    float remaining = g.max_leaf_size - leaf_size;
    float growth = std::min(max_growth, remaining);

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
    chemical(ChemicalID::Auxin) += growth_fraction * g.leaf_growth_auxin_multiplier * g.leaf_auxin_baseline;

    // Extend petiole proportionally as leaf grows
    // Target offset length = base_offset + petiole_length * (leaf_size / max_leaf_size)
    float olen = glm::length(offset);
    if (olen > 1e-4f) {
        float petiole_growth = g.leaf_petiole_length * (growth / g.max_leaf_size);
        offset += (offset / olen) * petiole_growth;
    }
}

float LeafNode::maintenance_cost(const WorldParams& world) const {
    return world.sugar_maintenance_leaf * leaf_size * leaf_size;
}

} // namespace botany
