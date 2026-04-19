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
    float sugar_cost_activation   = 0.01f;      // g glucose per meristem activation (bud break) — low; dormant buds accumulate little sugar
    float sugar_cost_leaf_growth  = 1.5f;       // g glucose / dm of leaf expansion (leaf tissue is expensive — 1.5 g/g dry mass)
    float sugar_cost_phototropism = 0.001f;      // g glucose / radian of leaf turning

    // Maintenance costs — respiration energy per unit tissue per hour.
    // Leaf: ~10% of gross photosynthesis. Stem/root: living ring (cambium, inner bark,
    // ray parenchyma) around a dead wood core — scales with surface area, not volume.
    // Meristem: metabolically active growing tissue.
    float sugar_maintenance_leaf     = 0.002f;   // g glucose / (dm² leaf area · hr)
    float sugar_maintenance_stem     = 0.002f;   // g glucose / (dm² living bark surface · hr)
                                                 //   formula: π × radius × length × rate
                                                 //   young stem (r=0.015, L=0.5): 4.7×10⁻⁵ g/tick
                                                 //   mature trunk (r=0.1, L=1.0): 6.3×10⁻⁴ g/tick
                                                 //   real range: 0.001–0.005 g/(dm²·hr); 0.002 is midpoint
    float sugar_maintenance_root     = 0.002f;   // g glucose / (dm² living root surface · hr)
                                                 //   same formula as stem — living ring around stele
    float sugar_maintenance_meristem = 0.0005f;  // g glucose / hr per active meristem tip
    float sugar_meristem_photosynthesis = 0.0f;  // shoot meristem self-photosynthesis as a multiple of maintenance cost (0 = heterotrophic, real apex meristems depend entirely on leaf sugar)

    float light_cell_size         = 0.075f;     // dm — shadow map cell size (smaller = higher resolution, more cells)
    glm::vec3 light_direction     = glm::vec3(0.50f, 1.0f, 0.0f); // unit vector pointing TOWARD light source
    // Sun direction: unit vector pointing TOWARD ground (same convention as LightSystem::sun_direction).
    // (0,-1,0) = overhead sun. Negated to get the toward-light vector used in phototropism.
    glm::vec3 sun_direction       = glm::vec3(0.0f, -1.0f, 0.0f);
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

    // Phloem (Münch pressure flow) parameters
    float max_phloem_velocity  = 10.0f;    // dm/tick (= 1 m/hr) — physical upper bound on phloem sap velocity.
                                           //   Caps: velocity = min(dp × conductance_per_pressure, max_phloem_velocity).
                                           //   Real phloem: 0.3–1.5 m/hr; 10 dm/tick is the high end.
    uint32_t phloem_iterations = 3;        // inner Jacobi iterations per phloem_resolve call.
                                           //   Each iteration propagates sugar one pipe section.
                                           //   3 = conservative; 5 = aggressive (travels further per tick).
    float phloem_reference_radius = 0.015f;// dm — reference radius (superseded by velocity-capped model,
                                           //   retained for reference; not used in phloem_resolve).
    float phloem_ring_thickness   = 0.002f;// dm (= 0.2 mm) — thickness of the living phloem+cambium ring.
                                           //   Used in phloem_ring_area(r, t) = π × (2r×t − t²).
                                           //   Used ONLY for pipe capacity (flow_vol = velocity × ring_area).
                                           //   Real dicots: ~0.1–0.5 mm; constant regardless of stem diameter.
    float max_sugar_concentration = 300.0f;// g/dm³ — upper bulk-concentration cap (~30% sucrose solution).
                                           //   Sugar cap per node = node_volume(n) × max_sugar_concentration.
    float leaf_thickness          = 0.003f;// dm — leaf mesophyll depth for node_volume of LEAF nodes.
                                           //   node_volume(leaf) = leaf_size² × leaf_thickness.

    // Debug / diagnostics
    bool vascular_debug_log = false;        // write per-junction vascular flow to debug/vascular_log.csv
    bool phloem_debug_log   = false;        // write per-node Münch phloem data to debug/phloem_log.csv
                                            //   columns: tick, node_id, node_type, parent_id,
                                            //     sugar (pre-phloem), volume, concentration, pressure,
                                            //     water_fraction, sugar_loaded_from_leaf,
                                            //     sugar_unloaded_to_meristem, sugar_flow_in,
                                            //     sugar_flow_out, net_flow
                                            //   + a SUMMARY row per tick with conservation_error check
    uint32_t current_tick = 0;              // set by caller each tick; used to label vascular log rows
};

inline WorldParams default_world_params() {
    return WorldParams{
        .light_level = 1.0f,
        .starvation_ticks_max = 1200,
        .sugar_production_rate    = 0.02f,
        .sugar_cost_meristem_growth = 2.0f,
        .sugar_cost_root_growth  = 1.5f,
        .sugar_cost_stem_growth  = 1.0f,
        .sugar_cost_activation   = 0.01f,
        .sugar_cost_leaf_growth  = 1.5f,
        .sugar_cost_phototropism = 0.001f,
        .sugar_maintenance_leaf     = 0.002f,
        .sugar_maintenance_stem     = 0.002f,
        .sugar_maintenance_root     = 0.002f,
        .sugar_maintenance_meristem = 0.0005f,
        .sugar_meristem_photosynthesis = 0.0f,
        .light_cell_size         = 0.075f,
        .light_direction         = glm::vec3(0.0f, 1.0f, 0.0f),
        .sun_direction           = glm::vec3(0.0f, -1.0f, 0.0f),
        .gravity                 = 9.81f,
        .break_strength_factor   = 5.0f,
        .reference_wood_density  = 50.0f,
        .leaf_mass_density       = 5.0f,
        .meristem_mass           = 0.1f,
        .ground_support_height   = 0.5f,
        .droop_rate              = 0.01f,
        .max_phloem_velocity     = 10.0f,
        .phloem_iterations       = 3,
        .phloem_reference_radius = 0.015f,
        .phloem_ring_thickness   = 0.002f,
        .max_sugar_concentration = 300.0f,
        .leaf_thickness          = 0.003f,
    };
}

WorldParams load_world_params(const std::string& path);

} // namespace botany
