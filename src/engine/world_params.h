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
    float soil_moisture = 1.0f;                  // fraction of field capacity (1.0 = saturated)
    uint32_t starvation_ticks_max = 2200;       // hours (~50 days) until starvation death

    // Construction costs — energy required per unit of biological growth.
    // Derived from ~1.2-1.4 g glucose / g dry mass (includes 25% growth respiration
    // overhead), scaled by typical cross-section area at initial radius.
    // These are per-distance, not per-time — unchanged by tick rate.
    float sugar_production_rate    = 0.02f;      // g glucose / (dm² leaf area · hr) at full sun — photosynthetic constant
    float sugar_cost_meristem_growth = 1.5f;    // g glucose / dm of shoot meristem tip extension
    float sugar_cost_root_growth  = 1.5f;       // g glucose / dm of root tip extension
    float sugar_cost_stem_growth  = 1.0f;       // g glucose / dm of stem/root internode growth (thickening + elongation)
    float sugar_cost_activation   = 0.3f;       // g glucose per meristem activation (bud break)
    float sugar_cost_leaf_growth  = 1.5f;       // g glucose / dm of leaf expansion (leaf tissue is expensive — 1.5 g/g dry mass)
    float sugar_cost_phototropism = 0.001f;      // g glucose / radian of leaf turning

    // Maintenance costs — respiration energy per unit tissue per hour.
    // Leaf: ~10% of gross photosynthesis. Stem: mid-range sapwood (mostly dead cells).
    // Root: slightly below stem (fine roots are active but less dense).
    // Meristem: metabolically active growing tissue.
    float sugar_maintenance_leaf     = 0.002f;   // g glucose / (dm² leaf area · hr)
    float sugar_maintenance_stem     = 0.01f;    // g glucose / (dm³ stem volume · hr)
    float sugar_maintenance_root     = 0.004f;   // g glucose / (dm³ root volume · hr)
    float sugar_maintenance_meristem = 0.0005f;  // g glucose / hr per active meristem tip
    float sugar_meristem_photosynthesis = 1.0f;  // shoot meristem sugar production at full light, as multiple of maintenance cost

    float light_cell_size         = 0.075f;     // dm — shadow map cell size (smaller = higher resolution, more cells)
    glm::vec3 light_direction     = glm::vec3(0.50f, 1.0f, 0.0f); // unit vector pointing TOWARD light source
    uint32_t light_update_interval = 10;  // ticks between shadow map recomputation
    bool cpu_light_enabled        = true;  // set false when GPU LightSystem is active

    // Stress physics
    float gravity = 9.81f;                  // m/s² — gravitational acceleration
    float break_strength_factor = 5000.0f;   // stress units per (g/dm³) of wood density — wood is strong
    float reference_wood_density = 50.0f;   // g/dm³ — density at which sugar costs are calibrated
    float leaf_mass_density = 1.0f;         // g/dm² of leaf area (~2g for a 15cm leaf)
    float meristem_mass = 0.1f;             // g — fixed mass for meristem tips
    float ground_support_height = 0.5f;     // dm — below this Y, stress is zeroed
    float droop_rate = 0.01f;               // radians/tick — max angular droop when overstressed
};

inline WorldParams default_world_params() {
    return WorldParams{
        .light_level = 1.0f,
        .starvation_ticks_max = 1200,
        .sugar_production_rate    = 0.02f,
        .sugar_cost_meristem_growth = 2.0f,
        .sugar_cost_root_growth  = 1.5f,
        .sugar_cost_stem_growth  = 1.0f,
        .sugar_cost_activation   = 0.3f,
        .sugar_cost_leaf_growth  = 1.5f,
        .sugar_cost_phototropism = 0.001f,
        .sugar_maintenance_leaf     = 0.002f,
        .sugar_maintenance_stem     = 0.01f,
        .sugar_maintenance_root     = 0.004f,
        .sugar_maintenance_meristem = 0.0005f,
        .sugar_meristem_photosynthesis = 1.0f,
        .light_cell_size         = 0.075f,
        .light_direction         = glm::vec3(0.0f, 1.0f, 0.0f),
        .gravity                 = 9.81f,
        .break_strength_factor   = 5.0f,
        .reference_wood_density  = 50.0f,
        .leaf_mass_density       = 5.0f,
        .meristem_mass           = 0.1f,
        .ground_support_height   = 0.5f,
        .droop_rate              = 0.01f,
    };
}

WorldParams load_world_params(const std::string& path);

} // namespace botany
