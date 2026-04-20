#include "engine/vascular_sub_stepped.h"
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/world_params.h"

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

} // namespace botany
