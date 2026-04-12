// src/engine/world_params.h
#pragma once

#include <cstdint>
#include <string>
#include <glm/vec3.hpp>

namespace botany {

// Unit system:
//   Distance: 1 unit = 1 dm (10 cm)
//   Time:     1 tick = 1 hour
//   Sugar:    grams of glucose (g)
//   Light:    fraction of full sun (1.0 = standard clear-sky PAR ~2000 µmol/m²/s)

struct WorldParams {
    float light_level = 1.0f;                   // fraction of full sun
    uint32_t starvation_ticks_max = 1200;       // hours (~50 days) until starvation death

    // Construction costs — energy required per unit of biological growth.
    // Derived from ~1.2-1.4 g glucose / g dry mass (includes 25% growth respiration
    // overhead), scaled by typical cross-section area at initial radius.
    // These are per-distance, not per-time — unchanged by tick rate.
    float sugar_cost_shoot_growth = 2.0f;       // g glucose / dm of shoot extension
    float sugar_cost_root_growth  = 1.5f;       // g glucose / dm of root extension (roots ~15% cheaper)
    float sugar_cost_thickening   = 0.5f;       // g glucose / dm of radial thickening
    float sugar_cost_elongation   = 1.0f;       // g glucose / dm of internode elongation (cell expansion, cheaper than new growth)
    float sugar_cost_activation   = 0.3f;       // g glucose per meristem activation (bud break)
    float sugar_cost_leaf_growth  = 1.5f;       // g glucose / dm of leaf expansion (leaf tissue is expensive — 1.5 g/g dry mass)
    float sugar_cost_phototropism = 0.001f;      // g glucose / radian of leaf turning
    float light_cell_size         = 0.075f;     // dm — shadow map cell size (smaller = higher resolution, more cells)
    glm::vec3 light_direction     = glm::vec3(0.0f, 1.0f, 0.0f); // unit vector pointing TOWARD light source
};

inline WorldParams default_world_params() {
    return WorldParams{
        .light_level = 1.0f,
        .starvation_ticks_max = 1200,
        .sugar_cost_shoot_growth = 2.0f,
        .sugar_cost_root_growth  = 1.5f,
        .sugar_cost_thickening   = 0.5f,
        .sugar_cost_elongation   = 1.0f,
        .sugar_cost_activation   = 0.3f,
        .sugar_cost_leaf_growth  = 1.5f,
        .sugar_cost_phototropism = 0.001f,
        .light_cell_size         = 0.075f,
        .light_direction         = glm::vec3(0.0f, 1.0f, 0.0f),
    };
}

WorldParams load_world_params(const std::string& path);

} // namespace botany
