// src/engine/genome.h
#pragma once

#include <cstdint>

namespace botany {

// Unit system: distance = dm (10cm), time = hours, sugar = grams glucose.
// See world_params.h for full unit documentation.

struct Genome {
    // Hormone production & sensitivity (dimensionless signaling units)
    float auxin_production_rate;
    float auxin_transport_rate;
    float auxin_directional_bias;     // -1=basipetal only, 0=gradient, +1=acropetal only
    float auxin_decay_rate;
    float auxin_threshold;

    float cytokinin_production_rate;
    float cytokinin_transport_rate;
    float cytokinin_directional_bias; // -1=basipetal only, 0=gradient, +1=acropetal only
    float cytokinin_decay_rate;
    float cytokinin_threshold;

    // Shoot growth
    float growth_rate;                // dm/hr — shoot tip extension speed
    float max_internode_length;       // dm — max internode length (well-fed, lucky)
    float min_internode_length;       // dm — min internode length (starved/unlucky)
    float branch_angle;               // radians
    float thickening_rate;            // dm/hr — radial growth of interior nodes
    float internode_elongation_rate;  // dm/hr — intercalary stretch of young internodes
    uint32_t internode_maturation_ticks; // hours until internode stops elongating

    // Root growth
    float root_growth_rate;           // dm/hr — root tip extension speed
    float root_max_internode_length;  // dm — max root internode length
    float root_min_internode_length;  // dm — min root internode length
    float root_branch_angle;          // radians
    float root_internode_elongation_rate;  // dm/hr
    uint32_t root_internode_maturation_ticks; // hours
    float root_gravitropism_strength; // how strongly roots turn downward near surface
    float root_gravitropism_depth;    // dm — depth at which gravitropism correction begins

    // Geometry
    float max_leaf_size;              // dm — maximum leaf side-length
    float leaf_growth_rate;           // dm/hr — how fast leaves grow from bud to max
    float leaf_bud_size;              // dm — initial size of a leaf bud
    float initial_radius;             // dm — stem radius at creation (~5mm)
    float root_initial_radius;        // dm — root radius at creation (~2.5mm)
    float tip_offset;                 // dm — small forward offset when chaining nodes
    float growth_noise;               // radians — max angular perturbation per segment

    // Sugar / photosynthesis (g glucose)
    float sugar_production_rate;      // g glucose / (dm leaf_size · hr) at full sun
    float sugar_transport_conductance;// conductance scaling for diffusion between nodes
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

    // Sugar save thresholds — minimum reserve before growth occurs (g glucose)
    float sugar_save_shoot;           // reserve for shoot apical meristems
    float sugar_save_root;            // reserve for root apical meristems
    float sugar_save_stem;            // reserve for stem thickening

    float leaf_phototropism_rate;     // radians/hr — how fast leaves turn toward light

    // Activation thresholds — parent node sugar needed for bud break (g glucose)
    float sugar_activation_shoot;     // parent sugar needed for shoot axillary activation
    float sugar_activation_root;      // parent sugar needed for root axillary activation

    // Gibberellin — promotes internode elongation, produced by young leaves
    float ga_production_rate;         // GA produced per dm of leaf_size per tick
    uint32_t ga_leaf_age_max;         // ticks — only leaves younger than this produce GA
    float ga_elongation_sensitivity;  // how strongly GA boosts elongation rate
    float ga_length_sensitivity;      // how strongly GA boosts target internode length
    float ga_transport_rate;          // fraction transported per tick (biased transport)
    float ga_directional_bias;        // -1=basipetal, 0=gradient, +1=acropetal
    float ga_decay_rate;              // exponential decay per tick

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
};

inline Genome default_genome() {
    return Genome{
        .auxin_production_rate = 1.0f,
        .auxin_transport_rate = 0.3f,
        .auxin_directional_bias = -0.9f,
        .auxin_decay_rate = 0.05f,
        .auxin_threshold = 0.15f,

        .cytokinin_production_rate = 1.0f,
        .cytokinin_transport_rate = 0.3f,
        .cytokinin_directional_bias = 0.9f,
        .cytokinin_decay_rate = 0.05f,
        .cytokinin_threshold = 0.15f,

        .growth_rate = 0.008f,              // ~2 cm/day = 0.8 mm/hr
        .max_internode_length = 1.0f,       // 10 cm
        .min_internode_length = 0.3f,       // 3 cm — starved/unlucky minimum
        .branch_angle = 0.785f,             // ~45 degrees
        .thickening_rate = 0.00004f,        // ~3.5 mm radius/year
        .internode_elongation_rate = 0.004f, // dm/hr — half of tip growth rate
        .internode_maturation_ticks = 72,    // 3 days until lockout

        .root_growth_rate = 0.004f,         // ~1 cm/day = 0.4 mm/hr
        .root_max_internode_length = 0.8f,  // 8 cm
        .root_min_internode_length = 0.2f,  // 2 cm
        .root_branch_angle = 0.35f,         // ~20 degrees
        .root_internode_elongation_rate = 0.002f, // dm/hr
        .root_internode_maturation_ticks = 48,    // 2 days
        .root_gravitropism_strength = 2.0f,
        .root_gravitropism_depth = 0.5f,

        .max_leaf_size = 0.3f,              // 3 cm side-length at maturity
        .leaf_growth_rate = 0.001f,         // ~0.24 mm/hr — full size in ~250 hrs (~10 days)
        .leaf_bud_size = 0.02f,             // 2 mm bud
        .initial_radius = 0.05f,            // 5 mm
        .root_initial_radius = 0.025f,      // 2.5 mm
        .tip_offset = 0.01f,
        .growth_noise = 0.26f,              // ~15 degrees

        .sugar_production_rate = 0.012f,    // g glucose / (dm leaf · hr) — ~7 g/m²/day
        .sugar_transport_conductance = 1.0f,  // ~5% gradient transfer per iter at 5mm radius
        .sugar_maintenance_leaf = 0.013f,   // g / (dm² · hr) — leaves dominate maintenance budget
        .sugar_maintenance_stem = 0.028f,   // g / (dm³ · hr) — wood is cheap per volume
        .sugar_maintenance_root = 0.135f,   // g / (dm³ · hr) — fine roots expensive (high turnover)
        .sugar_maintenance_meristem = 0.001f, // high per mass, small organ

        .seed_sugar = 28.0f,                 // ~15 days heterotrophic growth

        .sugar_storage_density_wood = 50.0f,  // g glucose max / dm³ — ~5% of dry mass as starch
        .sugar_storage_density_leaf = 0.5f,   // g glucose max / dm² — thin tissue, less storage
        .sugar_cap_minimum = 0.05f,           // floor for tiny/new nodes
        .sugar_cap_meristem = 2.0f,           // meristem tips — must hold enough sugar for active growth

        .sugar_save_shoot = 0.01f,          // buffer before growth
        .sugar_save_root = 0.005f,
        .sugar_save_stem = 0.02f,

        .leaf_phototropism_rate = 0.02f,    // ~1.1 deg/hr — full correction in ~3 days

        .sugar_activation_shoot = 0.5f,     // well-fed branch needed for shoot bud break
        .sugar_activation_root = 0.3f,      // roots branch with less sugar

        // Gibberellin
        .ga_production_rate = 0.5f,
        .ga_leaf_age_max = 168,               // 7 days
        .ga_elongation_sensitivity = 2.0f,
        .ga_length_sensitivity = 1.5f,
        .ga_transport_rate = 0.2f,            // moderate transport
        .ga_directional_bias = -0.7f,         // mostly basipetal (leaf -> parent -> trunk)
        .ga_decay_rate = 0.15f,               // decays faster than auxin — short-range signal

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
    };
}

} // namespace botany
