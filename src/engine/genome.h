// src/engine/genome.h
#pragma once

#include <cstdint>

namespace botany {

// Unit system: distance = dm (10cm), time = hours, sugar = grams glucose.
// See world_params.h for full unit documentation.

struct Genome {
    // Hormone production & sensitivity (dimensionless signaling units)
    float auxin_production_rate;
    float auxin_diffusion_rate;       // fraction diffused per tick
    float auxin_decay_rate;
    float auxin_threshold;

    float cytokinin_production_rate;    // cytokinin per g sugar produced by leaves (leaf productivity signal)
    float cytokinin_diffusion_rate;    // fraction diffused per tick
    float cytokinin_decay_rate;
    float cytokinin_threshold;
    float cytokinin_growth_threshold;  // cytokinin level for full-speed growth (gates all growth processes)

    // Shoot growth
    float growth_rate;                // dm/hr — shoot tip extension speed
    uint32_t shoot_plastochron;       // ticks between node creation (time-based, like real meristems)
    float branch_angle;               // radians
    float thickening_rate;            // dm/hr — radial growth of interior nodes
    float internode_elongation_rate;  // dm/hr — intercalary stretch of young internodes
    float max_internode_length;       // dm — max internode length (elongation target)
    uint32_t internode_maturation_ticks; // hours until internode stops elongating

    // Root growth
    float root_growth_rate;           // dm/hr — root tip extension speed
    uint32_t root_plastochron;        // ticks between root node creation
    float root_branch_angle;          // radians
    float root_internode_elongation_rate;  // dm/hr
    uint32_t root_internode_maturation_ticks; // hours
    float root_gravitropism_strength; // how strongly roots turn downward near surface
    float root_gravitropism_depth;    // dm — depth at which gravitropism correction begins

    // Geometry
    float max_leaf_size;              // dm — maximum leaf side-length
    float leaf_growth_rate;           // dm/hr — how fast leaves grow from bud to max
    float leaf_bud_size;              // dm — initial size of a leaf bud
    float leaf_petiole_length;        // dm — stalk length holding leaf away from stem
    float initial_radius;             // dm — stem radius at creation (~5mm)
    float root_initial_radius;        // dm — root radius at creation (~2.5mm)
    float tip_offset;                 // dm — small forward offset when chaining nodes
    float growth_noise;               // radians — max angular perturbation per segment

    // Sugar / photosynthesis (g glucose)
    float sugar_production_rate;      // g glucose / (dm² leaf area · hr) at full sun
    float sugar_diffusion_rate;       // fraction diffused per tick
    float sugar_maintenance_leaf;     // g glucose / (dm² leaf area · hr)
    float sugar_maintenance_stem;     // g glucose / (dm³ volume · hr)
    float sugar_maintenance_root;     // g glucose / (dm³ volume · hr)
    float sugar_maintenance_meristem; // g glucose / hr per active meristem tip

    float seed_sugar;                 // g glucose — initial reserves in the seed

    // Sugar storage caps — maximum sugar a node can hold, proportional to tissue volume.
    // Prevents unbounded accumulation; represents finite starch storage capacity.
    float sugar_storage_density_wood; // g glucose max / dm³ of stem/root tissue
    float sugar_storage_density_leaf; // g glucose max / dm² of leaf area
    float sugar_cap_minimum;          // g glucose — floor for tiny/new stem/root/leaf nodes
    float sugar_cap_meristem;         // g glucose — cap for meristem nodes (must hold growth sugar)

    float leaf_phototropism_rate;     // radians/hr — how fast leaves turn toward light

    // Activation thresholds — parent node sugar needed for bud break (g glucose)
    // Gibberellin — promotes internode elongation, produced by young leaves
    float ga_production_rate;         // GA produced per dm of leaf_size per tick
    uint32_t ga_leaf_age_max;         // ticks — only leaves younger than this produce GA
    float ga_elongation_sensitivity;  // how strongly GA boosts elongation rate
    float ga_length_sensitivity;      // how strongly GA boosts target internode length
    float ga_diffusion_rate;          // fraction diffused per tick
    float ga_decay_rate;              // exponential decay per tick

    // Leaf abscission — carbon-balance driven
    uint32_t leaf_abscission_ticks;        // ticks of net sugar deficit before senescence
    uint32_t min_leaf_age_before_abscission; // ticks — young leaves are immune (still growing)

    // Ethylene — stress/crowding gas signal, triggers leaf abscission
    float ethylene_starvation_rate;       // production when sugar = 0
    float ethylene_shade_rate;            // production from low light
    float ethylene_shade_threshold;       // light_exposure below which shade-ethylene kicks in
    float ethylene_age_rate;              // production ramp from old age
    uint32_t ethylene_age_onset;          // tick age when age-ethylene starts
    float ethylene_crowding_rate;         // production per nearby node
    float ethylene_crowding_radius;       // dm — radius for crowding density count
    float ethylene_diffusion_radius;      // dm — gas cloud spread distance
    float ethylene_abscission_threshold;  // ethylene level triggering leaf senescence
    float ethylene_elongation_inhibition; // strength of elongation suppression
    uint32_t senescence_duration;         // ticks from senescence start to leaf drop

    // Stress — mechanical load response
    float wood_density;                   // g/dm³ — mass per volume, also determines strength
    float wood_flexibility;               // 0-1 — droop threshold as fraction of break threshold
    float stress_hormone_production_rate; // hormone per unit stress
    float stress_hormone_diffusion_rate;  // fraction diffused per tick
    float stress_hormone_decay_rate;      // fraction decayed per tick
    float stress_thickening_boost;        // thickening multiplier per unit stress hormone
    float stress_elongation_inhibition;   // elongation suppression per unit stress hormone
    float stress_gravitropism_boost;      // gravitropism pull per unit stress hormone
};

inline Genome default_genome() {
    return Genome{
        .auxin_production_rate = 0.15f,
        .auxin_diffusion_rate = 0.3f,
        .auxin_decay_rate = 0.15f,
        .auxin_threshold = 0.15f,

        .cytokinin_production_rate = 5.0f,   // cytokinin per g sugar produced by leaves
        .cytokinin_diffusion_rate = 0.3f,
        .cytokinin_decay_rate = 0.05f,
        .cytokinin_threshold = 0.15f,
        .cytokinin_growth_threshold = 0.1f,

        .growth_rate = 0.008f,              // ~2 cm/day = 0.8 mm/hr
        .shoot_plastochron = 24,            // 1 day between node creation (like real meristems)
        .branch_angle = 0.785f,             // ~45 degrees
        .thickening_rate = 0.00004f,        // ~3.5 mm radius/year
        .internode_elongation_rate = 0.004f, // dm/hr — intercalary stretch after creation
        .max_internode_length = 1.0f,       // 10 cm — elongation target
        .internode_maturation_ticks = 72,    // 3 days until lockout

        .root_growth_rate = 0.004f,         // ~1 cm/day = 0.4 mm/hr
        .root_plastochron = 24,             // 1 day between root node creation
        .root_branch_angle = 0.35f,         // ~20 degrees
        .root_internode_elongation_rate = 0.002f, // dm/hr
        .root_internode_maturation_ticks = 48,    // 2 days
        .root_gravitropism_strength = .20f,
        .root_gravitropism_depth = 0.5f,

        .max_leaf_size = 1.5f,              // 15 cm side-length at maturity (realistic broad leaf)
        .leaf_growth_rate = 0.005f,         // ~0.5 mm/hr — full size (1.5dm) in ~300 hrs (~12 days)
        .leaf_bud_size = 0.02f,             // 2 mm bud
        .leaf_petiole_length = 0.5f,        // 5 cm petiole (realistic for broad leaves)
        .initial_radius = 0.05f,            // 5 mm
        .root_initial_radius = 0.025f,      // 2.5 mm
        .tip_offset = 0.01f,
        .growth_noise = 0.26f,              // ~15 degrees

        .sugar_production_rate = 0.02f,      // g glucose / (dm² leaf area · hr) — maintenance is ~5% of full-sun production
        .sugar_diffusion_rate = 0.8f,        // high base rate — radius scaling is the real throttle
        .sugar_maintenance_leaf = 0.001f,    // g / (dm² · hr) — ~25% of gross production (realistic leaf respiration)
        .sugar_maintenance_stem = 0.028f,   // g / (dm³ · hr) — wood is cheap per volume
        .sugar_maintenance_root = 0.005f,   // g / (dm³ · hr) — fine roots expensive (high turnover)
        .sugar_maintenance_meristem = 0.0001f, // tiny organ — costs ~1/3 of one mature leaf

        .seed_sugar = 48.0f,                 // ~15 days heterotrophic growth

        .sugar_storage_density_wood = 500.0f, // g glucose max / dm³ — high cap so stems can pass sugar through
        .sugar_storage_density_leaf = 2.0f,   // g glucose max / dm² — enough buffer for export
        .sugar_cap_minimum = 0.1f,            // floor for tiny/new nodes
        .sugar_cap_meristem = 2.0f,           // meristem tips — must hold enough sugar for active growth

        .leaf_phototropism_rate = 0.02f,    // ~1.1 deg/hr — full correction in ~3 days

        // Gibberellin
        .ga_production_rate = 0.5f,
        .ga_leaf_age_max = 168,               // 7 days
        .ga_elongation_sensitivity = 2.0f,
        .ga_length_sensitivity = 1.5f,
        .ga_diffusion_rate = 0.2f,
        .ga_decay_rate = 0.15f,

        // Leaf abscission
        .leaf_abscission_ticks = 500,         // ~21 days of net deficit
        .min_leaf_age_before_abscission = 240, // ~10 days — let new leaves grow to size

        // Ethylene
        .ethylene_starvation_rate = 0.3f,
        .ethylene_shade_rate = 0.2f,
        .ethylene_shade_threshold = 0.3f,
        .ethylene_age_rate = 0.05f,
        .ethylene_age_onset = 720,            // 30 days
        .ethylene_crowding_rate = 0.1f,
        .ethylene_crowding_radius = 0.5f,     // dm
        .ethylene_diffusion_radius = 1.0f,    // dm
        .ethylene_abscission_threshold = 0.5f,
        .ethylene_elongation_inhibition = 1.0f,
        .senescence_duration = 48,            // 2 days

        // Stress
        .wood_density = 50.0f,                    // g/dm³ — light deciduous wood
        .wood_flexibility = 0.5f,                 // droop starts at 50% of break stress
        .stress_hormone_production_rate = 0.001f,  // low signaling — stress values are large numbers
        .stress_hormone_diffusion_rate = 0.15f,   // moderate local diffusion
        .stress_hormone_decay_rate = 0.2f,        // fades quickly — local signal
        .stress_thickening_boost = 1.0f,          // 1:1 hormone-to-thickening boost
        .stress_elongation_inhibition = 1.0f,     // 1:1 hormone-to-elongation suppression
        .stress_gravitropism_boost = 0.5f,        // moderate upward correction
    };
}

} // namespace botany
