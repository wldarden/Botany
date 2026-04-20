#include "engine/vascular_sub_stepped.h"
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/node/node.h"
#include <glm/geometric.hpp>
#include <cmath>

namespace botany {

void vascular_sub_stepped(Plant& /*plant*/, const Genome& /*g*/, const WorldParams& /*world*/) {
    // Empty stub — per-sub-step logic is added in Tasks 10-18.
}

float radial_permeability_sugar(float radius, const Genome& g) {
    const float base   = g.base_radial_permeability_sugar;
    const float floor  = g.radial_floor_fraction_sugar;
    const float r_half = g.radial_half_radius_sugar;
    const float ratio  = radius / r_half;
    return base * (floor + (1.0f - floor) / (1.0f + ratio * ratio));
}

float radial_permeability_water(float radius, const Genome& g) {
    const float base   = g.base_radial_permeability_water;
    const float floor  = g.radial_floor_fraction_water;
    const float r_half = g.radial_half_radius_water;
    const float ratio  = radius / r_half;
    return base * (floor + (1.0f - floor) / (1.0f + ratio * ratio));
}

// Returns the length of the internode (distance from parent to this node).
// For the seed (no parent), returns 1.0 so the seed has a meaningful capacity.
static float node_length(const Node& n) {
    if (!n.parent) return 1.0f;
    return glm::length(n.offset);
}

float phloem_capacity(const Node& n, const Genome& g) {
    if (!n.phloem()) return 0.0f;
    const float L = node_length(n);
    return 3.14159265358979f * n.radius * n.radius * L * g.phloem_fraction;
}

float xylem_capacity(const Node& n, const Genome& g) {
    if (!n.xylem()) return 0.0f;
    const float L = node_length(n);
    return 3.14159265358979f * n.radius * n.radius * L * g.xylem_fraction;
}

} // namespace botany
