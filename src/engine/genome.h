// src/engine/genome.h
#pragma once

#include <cstdint>

namespace botany {

// Unit system: distance = dm (10cm), time = hours, sugar = grams glucose.
// See world_params.h for full unit documentation.

struct Genome {
    // Hormone production & sensitivity (dimensionless signaling units)
    // was: auxin_production_rate
    float apical_auxin_baseline;
    float apical_growth_auxin_multiplier; // growth-scaled bonus: total = baseline * (1 + multiplier * gf)
    float auxin_diffusion_rate;       // fraction diffused per tick
    float auxin_decay_rate;
    float auxin_threshold;
    float auxin_shade_boost;           // shade-avoidance: production multiplier boost in low light (0 = none)
    float auxin_sugar_half_saturation; // g glucose — sugar level for half-max production (Michaelis-Menten)
    float auxin_water_half_saturation;      // ml — water level at which auxin production hits half-max. Floor 0.1 is encoded in the formula, not evolvable.
    float cytokinin_sugar_half_saturation;  // g glucose — sugar level at which CK production hits half-max.
    float cytokinin_water_half_saturation;  // ml — water level at which CK production hits half-max.
    float quiescence_threshold;             // ticks — active meristem reverts to dormant after this many consecutive low-sugar ticks (before starvation_ticks_max death).
    uint32_t starvation_ticks_max_stem;     // ticks — STEM nodes (non-meristem woody tissue) survive this many consecutive zero-sugar ticks before death. Stand-in until starch reserves are modeled.
    uint32_t starvation_ticks_max_root;     // ticks — ROOT nodes (non-meristem root tissue) survive this many consecutive zero-sugar ticks before death. Fine roots have less reserve than trunk.
    float auxin_bias;                  // equilibrium shift for basipetal flow (negative = toward root)
    float leaf_auxin_baseline;            // scaling constant for leaf auxin production (decoupled from apical)
    float leaf_growth_auxin_multiplier;   // fraction of leaf_auxin_baseline produced at max leaf growth

    // Auxin growth sensitivity — saturating Michaelis-Menten per tissue type.
    // max_boost: signed ceiling (positive = promotes, negative = inhibits growth).
    // half_saturation: auxin level for half-max effect.
    float stem_auxin_max_boost;
    float stem_auxin_half_saturation;
    float root_auxin_max_boost;
    float root_auxin_half_saturation;
    float leaf_auxin_max_boost;
    float leaf_auxin_half_saturation;
    float apical_auxin_max_boost;
    float apical_auxin_half_saturation;
    float root_apical_auxin_max_boost;
    float root_apical_auxin_half_saturation;

    float cytokinin_production_rate;    // cytokinin per g sugar produced by leaves (leaf productivity signal)
    float cytokinin_diffusion_rate;    // fraction diffused per tick
    float cytokinin_decay_rate;
    float cytokinin_threshold;
    float cytokinin_growth_threshold;  // cytokinin level for full-speed growth (gates all growth processes)
    float cytokinin_bias;              // equilibrium shift for acropetal flow (positive = toward tips)

    // Transport capacity — limits chemical throughput per connection per tick
    float hormone_base_transport;      // throughput floor for hormone signals
    float hormone_transport_scale;     // radius amplification for hormone throughput
    float sugar_base_transport;        // throughput floor for sugar
    float sugar_transport_scale;       // radius amplification for sugar throughput

    // Corticular photosynthesis (green stem / bark development)
    float stem_photosynthesis_rate;       // g glucose / (dm² stem surface · hr) at full light
    float stem_green_radius_threshold;    // dm — stems thicker than this have bark, no photosynthesis

    // Shoot growth
    float growth_rate;                // dm/hr — shoot tip extension speed
    uint32_t shoot_plastochron;       // ticks between node creation (time-based, like real meristems)
    float branch_angle;               // radians
    float cambium_responsiveness;     // dm/hr·bias — thickening rate per unit auxin_flow_bias (PIN saturation 0→1).
                                      // delta_radius = cambium_responsiveness × auxin_flow_bias × sugar_gf × stress_boost.
                                      // Zero = monocot (no secondary thickening).
    float internode_elongation_rate;  // dm/hr — intercalary stretch of young internodes
    float max_internode_length;       // dm — max internode length (elongation target)
    uint32_t internode_maturation_ticks; // hours until internode stops elongating (visual: stiff cylinders)

    // Root growth
    float root_growth_rate;           // dm/hr — root tip extension speed
    uint32_t root_plastochron;        // ticks between root node creation
    float root_branch_angle;          // radians
    float root_internode_elongation_rate;  // dm/hr
    uint32_t root_internode_maturation_ticks; // hours until root internode stops elongating
    float root_gravitropism_strength; // how strongly roots turn downward near surface
    float root_gravitropism_depth;    // dm — depth at which gravitropism correction begins
    float root_cytokinin_production_rate; // baseline cytokinin produced per tick (mirrors apical_auxin_baseline)
    float root_tip_auxin_production_rate; // auxin produced by root tips (PIN recycling local maximum)
    float root_auxin_growth_threshold;    // Km for auxin-gated root growth fraction
    float root_ck_growth_floor;            // Km for CK-gated root elongation — low so trace CK permits growth; zero CK stops it
    float root_auxin_activation_threshold; // min auxin to activate dormant root meristem
    float root_cytokinin_inhibition_threshold; // cytokinin above this inhibits root activation
    uint32_t primary_root_lateral_delay_internodes; // primary RA skips spawning dormant lateral buds for its first N internodes (head start for main tap root)

    // Geometry
    float max_leaf_size;              // dm — maximum leaf side-length
    float leaf_growth_rate;           // dm/hr — how fast leaves grow from bud to max
    float leaf_bud_size;              // dm — initial size of a leaf bud
    float leaf_petiole_length;        // dm — stalk length holding leaf away from stem
    float leaf_opacity;               // 0–1 — fraction of light blocked per leaf (0.85 blocks 85%, transmits 15%)
    float initial_radius;             // dm — stem radius at creation (~5mm)
    float root_initial_radius;        // dm — root radius at creation (~2.5mm)
    float tip_offset;                 // dm — small forward offset when chaining nodes
    float growth_noise;               // radians — max angular perturbation per segment

    // Sugar economy (g glucose)
    float sugar_diffusion_rate;       // fraction diffused per tick
    float seed_sugar;                 // g glucose — initial reserves in the seed

    // Sugar storage caps — maximum sugar a node can hold, proportional to tissue volume.
    // Prevents unbounded accumulation; represents finite starch storage capacity.
    float sugar_storage_density_wood; // g glucose max / dm³ of stem/root tissue
    float sugar_storage_density_leaf; // g glucose max / dm² of leaf area
    float sugar_cap_minimum;          // g glucose — floor for tiny/new stem/root/leaf nodes
    float sugar_cap_meristem;         // g glucose — cap for meristem nodes. Kept small (0.1g) to match
                                          // microscopic tissue volume; prevents outsized phloem demand.

    // Water economy (ml)
    float water_absorption_rate;          // ml / (dm² root surface · hr) per unit soil_moisture
    float transpiration_rate;             // ml / (dm² leaf area · hr) at full light
    float photosynthesis_water_ratio;     // ml water consumed per g sugar produced
    float water_storage_density_stem;     // ml / dm³ of stem/root tissue
    float water_storage_density_leaf;     // ml / dm² of leaf area
    float water_cap_meristem;             // ml — fixed cap for meristem nodes
    float water_diffusion_rate;           // fraction diffused per tick
    float water_bias;                     // upward equilibrium shift (positive = toward tips)
    float water_base_transport;           // throughput floor
    float water_transport_scale;          // radius scaling on throughput

    float leaf_phototropism_rate;     // radians/hr — how fast leaves turn toward light
    float meristem_gravitropism_rate; // base pull toward set-point angle (always-on plagiotropism)
    float meristem_phototropism_rate; // shade-scaled pull toward light direction

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
    float stress_hormone_threshold;       // stress ratio (0-1) below which no hormone produced
    float stress_hormone_production_rate; // hormone per unit excess stress ratio
    float stress_hormone_diffusion_rate;  // fraction diffused per tick
    float stress_hormone_decay_rate;      // fraction decayed per tick
    float stress_thickening_boost;        // thickening multiplier per unit stress hormone
    float stress_elongation_inhibition;   // elongation suppression per unit stress hormone
    float stress_gravitropism_boost;      // gravitropism pull per unit stress hormone
    float elastic_recovery_rate;          // radians/tick — spring-back toward rest direction

    // Canalization — PIN transport history biases auxin flow and cambium
    float smoothing_rate;                 // lerp rate for auxin_flow_bias toward current PIN saturation (~20 tick response)
    float canalization_weight;            // global scaling on auxin_flow_bias effect in local diffusion (0 = disabled)

    // PIN transport
    float pin_capacity_per_area;          // AU/(dm²·tick) — max auxin transport per unit cross-section at full efficiency.
                                          // Also the denominator in saturation = flux / (r² × pin_capacity_per_area).
    float pin_base_efficiency;            // [0–1] — cold-start PIN efficiency when auxin_flow_bias = 0 (constitutively active PINs)

    // Vascular transport
    float meristem_sink_fraction;         // fraction of sugar_cap_meristem a meristem can demand per tick.
                                          // Caps demand at cap×fraction instead of cap−current so meristems
                                          // can't out-compete leaves for limited phloem sugar. At 0.05 and
                                          // sugar_cap_meristem=0.1, max demand = 0.005g/tick — well below
                                          // typical leaf production (~0.025g/tick per leaf).
    float vascular_radius_threshold;      // dm — minimum radius for bulk vascular admission.
                                          // Set below initial_radius (0.015 dm) so newly spawned
                                          // internodes qualify from birth. Only pre-initial-radius
                                          // manually-created nodes stay excluded.

    // --- Compartmented vascular model (2026-04-19) ---
    // Radial flow permeability curve (per chemical class): perm(r) = base ×
    // (floor + (1 - floor) / (1 + (r / half_radius)²)).  Young thin stems
    // leak freely; mature thick trunks asymptote to `floor × base`.  See
    // spec section 6 for full rationale.

    // Phloem radial (sugar between stem's own local_env and own phloem):
    float base_radial_permeability_sugar;   // max permeability at r = 0
    float radial_floor_fraction_sugar;      // asymptote as r → ∞, fraction of base
    float radial_half_radius_sugar;         // dm — curve inflection

    // Xylem radial (water and cytokinin between stem's own local_env and own xylem):
    float base_radial_permeability_water;
    float radial_floor_fraction_water;
    float radial_half_radius_water;

    // Pipe cross-section fractions (what portion of stem cross-section is
    // sieve tubes vs vessel elements):
    float phloem_fraction;                  // 0–1, portion of π·r² that's phloem
    float xylem_fraction;                   // 0–1, portion of π·r² that's xylem

    // Source/sink targets:
    float leaf_reserve_fraction_sugar;      // leaf keeps this fraction of its
                                            // local sugar before loading to phloem
    float meristem_sink_target_fraction;    // meristem refills to this fraction
                                            // of sugar_cap per tick
    float leaf_turgor_target_fraction;      // leaves refill to this fraction of
                                            // water_cap via xylem pull
    float root_water_reserve_fraction;      // roots keep this fraction of their
                                            // local water before actively pumping
                                            // surplus into xylem.  The active pump
                                            // models root pressure — the osmotic
                                            // mechanism that pushes xylem sap
                                            // upward in young/small-leaf plants
                                            // when transpiration pull is weak.
                                            // 0.3 matches the leaf sugar reserve.
};

inline Genome default_genome() {
    return Genome{
        .apical_auxin_baseline = 0.15f,
        .apical_growth_auxin_multiplier = 2.0f,  // total = baseline * 3 at max growth
        .auxin_diffusion_rate = 0.05f,         // very slow polar transport — cell-to-cell only, ~1 cm/hr in real plants
        .auxin_decay_rate = 0.12f,             // moderate decay — half-life ~5.4 hours
        .auxin_threshold = 0.15f,
        .auxin_shade_boost = 0.5f,           // shade can increase production by 50%
        .auxin_sugar_half_saturation = 0.3f, // modest sugar needed for decent production
        .auxin_water_half_saturation = 0.1f,     // matches SA stomatal threshold magnitude — half-max at moderate turgor
        .cytokinin_sugar_half_saturation = 0.05f, // CK more sensitive to sugar than auxin; sharper response
        .cytokinin_water_half_saturation = 0.1f,  // matches auxin water sensitivity
        .quiescence_threshold = 150.0f,           // ~6 days — meristem goes dormant after sustained starvation, well before death (starvation_ticks_max = 2200)
        .starvation_ticks_max_stem = 8000,        // ~333 days — woody trunks carry significant parenchyma starch reserves
        .starvation_ticks_max_root = 6000,        // ~250 days — root tissue reserves somewhat less than trunk; fine roots die faster than woody roots
        .auxin_bias = -0.1f,                  // gentle basipetal shift (auxin accumulates toward root)
        .leaf_auxin_baseline = 0.15f,             // same scale as apical, but multiplier keeps it at 10%
        .leaf_growth_auxin_multiplier = 0.1f,     // single leaf at max growth = 10% of apical baseline

        .stem_auxin_max_boost = 0.5f,            // auxin promotes stem elongation by up to 50%
        .stem_auxin_half_saturation = 0.2f,
        .root_auxin_max_boost = -0.3f,            // auxin inhibits root elongation by up to 30%
        .root_auxin_half_saturation = 0.1f,       // roots are very sensitive
        .leaf_auxin_max_boost = 0.3f,             // auxin promotes leaf expansion by up to 30%
        .leaf_auxin_half_saturation = 0.2f,
        .apical_auxin_max_boost = 0.2f,           // mild promotion of tip extension
        .apical_auxin_half_saturation = 0.3f,     // high half-sat — apicals sit in high auxin
        .root_apical_auxin_max_boost = -0.2f,     // mild inhibition of root tip extension
        .root_apical_auxin_half_saturation = 0.1f,

        .cytokinin_production_rate = 5.0f,   // cytokinin per g sugar produced by leaves
        .cytokinin_diffusion_rate = 0.02f,         // low — same as water; long-distance via xylem bulk flow, not local diffusion
        .cytokinin_decay_rate = 0.05f,             // slower than auxin (0.12) — cytokinin travels via xylem bulk flow, less degradation in transit
        .cytokinin_threshold = 0.003f,        // just below xylem steady-state (~0.0047); activates when root signal is present
        .cytokinin_growth_threshold = 0.005f,       // low Km — trace cytokinin from roots is enough (xylem bulk flow in real plants means near-lossless delivery)
        .cytokinin_bias = 0.1f,               // gentle acropetal shift (cytokinin accumulates toward tips)

        .hormone_base_transport = 0.5f,      // generous floor — thin tips can still signal
        .hormone_transport_scale = 1.0f,     // moderate radius scaling for hormones
        .sugar_base_transport = 0.1f,        // higher floor — reduces diffusion wave artifacts
        .sugar_transport_scale = 5.0f,       // strong radius dependence for sugar

        // Corticular photosynthesis
        .stem_photosynthesis_rate = 0.005f,       // ~1/4 of leaf rate — stems are less efficient
        .stem_green_radius_threshold = 0.04f,     // 4mm radius — thicker than initial (1.5mm) but thin

        .growth_rate = 0.002f,              // ~5 mm/day = 0.2 mm/hr
        .shoot_plastochron = 24,            // 1 day between node creation (like real meristems)
        .branch_angle = 0.785f,             // ~45 degrees
        .cambium_responsiveness = 0.0002f, // dm/hr·bias — calibrated so main trunk (bias ~2.0) thickens at
                                            // ~0.00004 dm/hr, matching old thickening_rate. Lateral branches
                                            // with weaker canalization thicken proportionally less.
        .internode_elongation_rate = 0.004f, // dm/hr — intercalary stretch after creation
        .max_internode_length = 1.0f,       // 10 cm — elongation target
        .internode_maturation_ticks = 72,    // 3 days until elongation lockout (visual constraint)

        .root_growth_rate = 0.004f,         // ~1 cm/day = 0.4 mm/hr
        .root_plastochron = 24,             // 1 day between root node creation
        .root_branch_angle = 0.35f,         // ~20 degrees
        .root_internode_elongation_rate = 0.002f, // dm/hr
        .root_internode_maturation_ticks = 48,    // 2 days until elongation lockout
        .root_gravitropism_strength = .20f,
        .root_gravitropism_depth = 0.5f,
        .root_cytokinin_production_rate = 0.5f,    // raw rate (CK no longer gated on local auxin post-Task 1); set so ~20 active RAs feed xylem enough CK to drive SA growth
        .root_tip_auxin_production_rate = 0.002f,  // small PIN-recycling floor — self-equilibrium ~0.017 supports lateral root initiation without needing distant shoot auxin
        .root_auxin_growth_threshold = 0.10f,       // Km for auxin-gated root elongation (mirrors SAM's cytokinin_growth_threshold pattern)
        .root_ck_growth_floor = 0.001f,             // low Km — trace CK permits growth, zero stops it
        .root_auxin_activation_threshold = 0.01f,  // below RA self-eq (~0.017) so deep lateral RAs can activate without shoot-delivered auxin
        .root_cytokinin_inhibition_threshold = 0.075f, // lowered to reach into the realized xylem CK range (~0.09–0.10 at surface) so inhibition actually fires and keeps dormant laterals dormant
        .primary_root_lateral_delay_internodes = 5,    // primary tap root gets 5 internode head start before it starts laying down lateral buds

        .max_leaf_size = 1.5f,              // 15 cm side-length at maturity (realistic broad leaf)
        .leaf_growth_rate = 0.005f,         // ~0.5 mm/hr — full size (1.5dm) in ~300 hrs (~12 days)
        .leaf_bud_size = 0.02f,             // 2 mm bud
        .leaf_petiole_length = 0.5f,        // 5 cm petiole (realistic for broad leaves)
        .leaf_opacity = 0.85f,              // blocks 85%, transmits 15% — realistic for broadleaves
        .initial_radius = 0.015f,            // 1.5 mm — realistic seedling stem
        .root_initial_radius = 0.008f,       // 0.8 mm — realistic radicle
        .tip_offset = 0.01f,
        .growth_noise = 0.4f,               // ~23 degrees — larger than plagiotropism for organic shapes

        .sugar_diffusion_rate = 0.8f,        // high base rate — radius scaling is the real throttle
        .seed_sugar = 48.0f,                 // ~15 days heterotrophic growth

        .sugar_storage_density_wood = 500.0f, // g glucose max / dm³ — high cap so stems can pass sugar through
        .sugar_storage_density_leaf = 2.0f,   // g glucose max / dm² — enough buffer for export
        .sugar_cap_minimum = 0.1f,            // floor for tiny/new nodes
        .sugar_cap_meristem = 0.1f,           // meristem tips — microscopic tissue; 0.1g cap → 0.005g/tick demand

        // Water economy
        .water_absorption_rate = 0.5f,           // ml / (dm² · hr) — roots must keep up with full canopy transpiration
        .transpiration_rate = 0.04f,             // ml / (dm² · hr) — stomatal water loss
        .photosynthesis_water_ratio = 0.5f,      // 0.5 ml water per g sugar (small cost)
        .water_storage_density_stem = 800.0f,    // ml / dm³ — wood is ~80% water by volume
        .water_storage_density_leaf = 3.0f,      // ml / dm² — leaves hold water in vacuoles
        .water_cap_meristem = 1.0f,              // ml — small active reserve
        .water_diffusion_rate = 0.02f,           // local only (same range as auxin); xylem handles long-distance
        .water_bias = 0.05f,                     // slight upward bias (transpiration pull)
        .water_base_transport = 0.2f,            // throughput cap floor (rarely limiting at 0.02 diffusion)
        .water_transport_scale = 4.0f,           // radius scaling on cap (rarely limiting at 0.02 diffusion)

        .leaf_phototropism_rate = 0.02f,    // ~1.1 deg/hr — full correction in ~3 days
        .meristem_gravitropism_rate = 0.02f, // gentle always-on pull toward set-point angle
        .meristem_phototropism_rate = 0.1f,  // shade-scaled pull toward light

        // Gibberellin
        .ga_production_rate = 0.5f,
        .ga_leaf_age_max = 168,               // 7 days
        .ga_elongation_sensitivity = 2.0f,
        .ga_length_sensitivity = 1.5f,
        .ga_diffusion_rate = 0.05f,           // local only — GA signals the internode below young leaves, not the whole plant
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
        .stress_hormone_threshold = 0.3f,            // no hormone below 30% of breaking capacity
        .stress_hormone_production_rate = 0.5f,      // input is now 0-1 stress ratio excess
        .stress_hormone_diffusion_rate = 0.10f,   // local alarm signal — spreads a few nodes, not the whole plant
        .stress_hormone_decay_rate = 0.2f,        // fades quickly — local signal
        .stress_thickening_boost = 1.0f,          // 1:1 hormone-to-thickening boost
        .stress_elongation_inhibition = 1.0f,     // 1:1 hormone-to-elongation suppression
        .stress_gravitropism_boost = 0.5f,        // moderate upward correction
        .elastic_recovery_rate = 0.005f,          // slow spring-back (half of droop_rate)

        // Canalization
        .smoothing_rate = 0.1f,               // 10% per tick → ~20 tick response; natural decay when flux stops
        .canalization_weight = 1.0f,          // full effect by default

        // PIN transport
        .pin_capacity_per_area = 500.0f,      // AU/(dm²·tick) — primary calibration knob; adjust up to reduce sensitivity
        .pin_base_efficiency = 0.2f,          // cold-start: new stems transport 20% of capacity before PIN upregulates

        // Vascular transport
        .meristem_sink_fraction = 0.05f,      // max demand per tick = 5% of cap (0.005g at cap=0.1) — well
                                              // below leaf production; meristem can't starve the plant
        .vascular_radius_threshold = 0.01f,   // dm — below initial_radius (0.015), all new internodes qualify from birth

        // Compartmented vascular defaults — starting values from spec section 6.
        .base_radial_permeability_sugar = 1.0f,
        .radial_floor_fraction_sugar    = 0.1f,
        .radial_half_radius_sugar       = 0.3f,   // dm
        .base_radial_permeability_water = 1.0f,
        .radial_floor_fraction_water    = 0.1f,
        .radial_half_radius_water       = 0.3f,
        .phloem_fraction                = 0.3f,   // portion of π·r² allocated to sieve tubes.
                                                  //   Biological 5% is too small for our pipe
                                                  //   math: π·r²·L·0.05 at thin stems gives
                                                  //   sub-µg capacities and Jacobi max_move
                                                  //   clamps limit per-tick shoot delivery to
                                                  //   ~10⁻⁴ g — well below shoot demand.  Tuned
                                                  //   up to 30% for adequate flow rate.
        .xylem_fraction                 = 0.5f,   // same issue — bumped from 20% to 50%.
                                                  //   Real xylem can be 30–80% of stem cross-
                                                  //   section so 50% is within biological range.
        .leaf_reserve_fraction_sugar    = 0.3f,   // matches existing phloem_reserve_fraction value
        .meristem_sink_target_fraction  = 0.05f,  // matches existing meristem_sink_fraction
        .leaf_turgor_target_fraction    = 0.7f,   // leaves aim to fill to 70% of water_cap
        .root_water_reserve_fraction    = 0.3f,   // roots keep 30% of water_cap locally
                                                  // before pumping surplus into xylem
                                                  // (matches leaf_reserve_fraction_sugar).
    };
}

} // namespace botany
