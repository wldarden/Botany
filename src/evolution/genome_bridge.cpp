#include "evolution/genome_bridge.h"

namespace botany {

// Register a gene with strength computed as a percentage of its valid range.
// rate = probability of mutation per generation (0.1 = 10%)
// pct  = mutation stddev as fraction of (max - min)
static void reg(evolve::StructuredGenome& sg, const std::string& tag, float val,
                float rate, float min_val, float max_val, float pct) {
    float strength = (max_val - min_val) * pct;
    sg.add_gene({tag, {val}, {rate, strength, min_val, max_val}});
}

evolve::StructuredGenome build_genome_template(const Genome& g, float mutation_pct) {
    evolve::StructuredGenome sg;
    const float r = 0.1f;   // 10% chance of mutation per gene per generation
    const float p = mutation_pct;

    // --- Auxin group (7 genes) ---
    reg(sg, "apical_auxin_baseline",      g.apical_auxin_baseline,      r, 0.01f, 2.0f, p);
    reg(sg, "apical_growth_auxin_multiplier", g.apical_growth_auxin_multiplier, r, 0.0f, 10.0f, p);
    reg(sg, "auxin_diffusion_rate",       g.auxin_diffusion_rate,       r, 0.01f, 1.0f, p);
    reg(sg, "auxin_decay_rate",           g.auxin_decay_rate,           r, 0.001f, 0.5f, p);
    reg(sg, "auxin_threshold",            g.auxin_threshold,            r, 0.01f, 1.0f, p);
    reg(sg, "auxin_shade_boost",          g.auxin_shade_boost,          r, 0.0f, 2.0f, p);
    reg(sg, "auxin_sugar_half_saturation", g.auxin_sugar_half_saturation, r, 0.01f, 5.0f, p);
    reg(sg, "auxin_bias",                g.auxin_bias,                r, -1.0f, 0.0f, p);
    reg(sg, "leaf_auxin_baseline",            g.leaf_auxin_baseline,            r, 0.01f, 2.0f, p);
    reg(sg, "leaf_growth_auxin_multiplier",   g.leaf_growth_auxin_multiplier,   r, 0.0f, 1.0f, p);
    reg(sg, "stem_auxin_max_boost",              g.stem_auxin_max_boost,              r, -1.0f, 2.0f, p);
    reg(sg, "stem_auxin_half_saturation",        g.stem_auxin_half_saturation,        r, 0.01f, 1.0f, p);
    reg(sg, "root_auxin_max_boost",              g.root_auxin_max_boost,              r, -1.0f, 2.0f, p);
    reg(sg, "root_auxin_half_saturation",        g.root_auxin_half_saturation,        r, 0.01f, 1.0f, p);
    reg(sg, "leaf_auxin_max_boost",              g.leaf_auxin_max_boost,              r, -1.0f, 2.0f, p);
    reg(sg, "leaf_auxin_half_saturation",        g.leaf_auxin_half_saturation,        r, 0.01f, 1.0f, p);
    reg(sg, "apical_auxin_max_boost",            g.apical_auxin_max_boost,            r, -1.0f, 2.0f, p);
    reg(sg, "apical_auxin_half_saturation",      g.apical_auxin_half_saturation,      r, 0.01f, 1.0f, p);
    reg(sg, "root_apical_auxin_max_boost",       g.root_apical_auxin_max_boost,       r, -1.0f, 2.0f, p);
    reg(sg, "root_apical_auxin_half_saturation", g.root_apical_auxin_half_saturation, r, 0.01f, 1.0f, p);

    // --- Cytokinin group (6 genes) ---
    reg(sg, "cytokinin_production_rate",  g.cytokinin_production_rate,  r, 0.01f, 2.0f, p);
    reg(sg, "cytokinin_diffusion_rate",   g.cytokinin_diffusion_rate,   r, 0.01f, 1.0f, p);
    reg(sg, "cytokinin_decay_rate",       g.cytokinin_decay_rate,       r, 0.001f, 0.5f, p);
    reg(sg, "cytokinin_threshold",        g.cytokinin_threshold,        r, 0.01f, 1.0f, p);
    reg(sg, "cytokinin_growth_threshold", g.cytokinin_growth_threshold, r, 0.01f, 1.0f, p);
    reg(sg, "cytokinin_bias",            g.cytokinin_bias,            r, 0.0f, 1.0f, p);

    // --- Transport capacity group (4 genes) ---
    reg(sg, "hormone_base_transport",    g.hormone_base_transport,    r, 0.01f, 2.0f, p);
    reg(sg, "hormone_transport_scale",   g.hormone_transport_scale,   r, 0.1f, 5.0f, p);
    reg(sg, "sugar_base_transport",      g.sugar_base_transport,      r, 0.001f, 0.5f, p);
    reg(sg, "sugar_transport_scale",     g.sugar_transport_scale,     r, 0.5f, 20.0f, p);

    // --- Shoot growth group ---
    reg(sg, "growth_rate",                g.growth_rate,                r, 0.001f, 0.05f, p);
    reg(sg, "shoot_plastochron",          static_cast<float>(g.shoot_plastochron), r, 6.0f, 168.0f, p);
    reg(sg, "max_internode_length",       g.max_internode_length,       r, 0.05f,  3.0f, p);
    reg(sg, "branch_angle",              g.branch_angle,              r, 0.05f,  1.57f, p);
    reg(sg, "cambium_responsiveness",    g.cambium_responsiveness,    r, 0.0f, 0.001f, p);
    reg(sg, "internode_elongation_rate", g.internode_elongation_rate, r, 0.0005f, 0.02f, p);
    reg(sg, "internode_maturation_ticks", static_cast<float>(g.internode_maturation_ticks),
                                                                      r, 12.0f, 500.0f, p);

    // --- Root growth group ---
    reg(sg, "root_growth_rate",               g.root_growth_rate,               r, 0.001f, 0.05f, p);
    reg(sg, "root_plastochron",               static_cast<float>(g.root_plastochron), r, 6.0f, 168.0f, p);
    reg(sg, "root_branch_angle",             g.root_branch_angle,             r, 0.05f,  1.57f, p);
    reg(sg, "root_internode_elongation_rate", g.root_internode_elongation_rate, r, 0.0005f, 0.02f, p);
    reg(sg, "root_internode_maturation_ticks", static_cast<float>(g.root_internode_maturation_ticks),
                                                                               r, 12.0f, 500.0f, p);
    reg(sg, "root_gravitropism_strength",    g.root_gravitropism_strength,    r, 0.1f, 10.0f, p);
    reg(sg, "root_gravitropism_depth",       g.root_gravitropism_depth,       r, 0.1f, 5.0f, p);

    // --- Geometry group (8 genes) ---
    reg(sg, "max_leaf_size",          g.max_leaf_size,          r, 0.05f,   3.0f, p);
    reg(sg, "leaf_growth_rate",       g.leaf_growth_rate,       r, 0.0001f, 0.01f, p);
    reg(sg, "leaf_bud_size",          g.leaf_bud_size,          r, 0.005f,  0.1f, p);
    reg(sg, "leaf_petiole_length",   g.leaf_petiole_length,   r, 0.1f,    2.0f, p);
    reg(sg, "leaf_opacity",           g.leaf_opacity,           r, 0.3f,    1.0f, p);
    reg(sg, "initial_radius",         g.initial_radius,         r, 0.01f,   0.2f, p);
    reg(sg, "root_initial_radius",    g.root_initial_radius,    r, 0.005f,  0.1f, p);
    reg(sg, "tip_offset",             g.tip_offset,             r, 0.001f,  0.1f, p);
    reg(sg, "growth_noise",           g.growth_noise,           r, 0.01f,   0.8f, p);
    reg(sg, "leaf_phototropism_rate", g.leaf_phototropism_rate, r, 0.001f,  0.1f, p);
    reg(sg, "meristem_gravitropism_rate", g.meristem_gravitropism_rate, r, 0.001f, 0.2f, p);
    reg(sg, "meristem_phototropism_rate", g.meristem_phototropism_rate, r, 0.01f, 0.5f, p);
    reg(sg, "leaf_abscission_ticks", static_cast<float>(g.leaf_abscission_ticks),
                                                                      r, 48.0f, 2000.0f, p);
    reg(sg, "min_leaf_age_before_abscission", static_cast<float>(g.min_leaf_age_before_abscission),
                                                                      r, 24.0f, 1000.0f, p);

    // --- Sugar economy group ---
    reg(sg, "sugar_diffusion_rate",        g.sugar_diffusion_rate,        r, 0.01f,   1.0f, p);
    reg(sg, "seed_sugar",                  g.seed_sugar,                  r, 5.0f,    200.0f, p);
    reg(sg, "sugar_storage_density_wood",  g.sugar_storage_density_wood,  r, 50.0f,   2000.0f, p);
    reg(sg, "sugar_storage_density_leaf",  g.sugar_storage_density_leaf,  r, 0.05f,   5.0f, p);
    reg(sg, "sugar_cap_minimum",           g.sugar_cap_minimum,           r, 0.005f,  0.5f, p);
    reg(sg, "sugar_cap_meristem",          g.sugar_cap_meristem,          r, 0.1f,    10.0f, p);
    // --- Water economy group (10 genes) ---
    reg(sg, "water_absorption_rate",       g.water_absorption_rate,       r, 0.001f,  0.5f, p);
    reg(sg, "transpiration_rate",          g.transpiration_rate,          r, 0.001f,  0.5f, p);
    reg(sg, "photosynthesis_water_ratio",  g.photosynthesis_water_ratio,  r, 0.01f,   5.0f, p);
    reg(sg, "water_storage_density_stem",  g.water_storage_density_stem,  r, 100.0f,  2000.0f, p);
    reg(sg, "water_storage_density_leaf",  g.water_storage_density_leaf,  r, 0.1f,    10.0f, p);
    reg(sg, "water_cap_meristem",          g.water_cap_meristem,          r, 0.1f,    10.0f, p);
    reg(sg, "water_diffusion_rate",        g.water_diffusion_rate,        r, 0.1f,    1.0f, p);
    reg(sg, "water_bias",                  g.water_bias,                  r, 0.0f,    0.5f, p);
    reg(sg, "water_base_transport",        g.water_base_transport,        r, 0.01f,   1.0f, p);
    reg(sg, "water_transport_scale",       g.water_transport_scale,       r, 0.5f,    20.0f, p);
    // --- Gibberellin group (7 genes) ---
    reg(sg, "ga_production_rate",        g.ga_production_rate,        r, 0.01f,  2.0f, p);
    reg(sg, "ga_leaf_age_max",           static_cast<float>(g.ga_leaf_age_max),
                                                                      r, 24.0f,  2000.0f, p);
    reg(sg, "ga_elongation_sensitivity", g.ga_elongation_sensitivity, r, 0.1f,   5.0f, p);
    reg(sg, "ga_length_sensitivity",     g.ga_length_sensitivity,     r, 0.1f,   5.0f, p);
    reg(sg, "ga_diffusion_rate",         g.ga_diffusion_rate,         r, 0.01f,  1.0f, p);
    reg(sg, "ga_decay_rate",             g.ga_decay_rate,             r, 0.01f,  0.5f, p);

    // --- Ethylene group (11 genes) ---
    reg(sg, "ethylene_starvation_rate",       g.ethylene_starvation_rate,       r, 0.01f,  2.0f, p);
    reg(sg, "ethylene_shade_rate",            g.ethylene_shade_rate,            r, 0.01f,  2.0f, p);
    reg(sg, "ethylene_shade_threshold",       g.ethylene_shade_threshold,       r, 0.05f,  1.0f, p);
    reg(sg, "ethylene_age_rate",              g.ethylene_age_rate,              r, 0.001f, 0.5f, p);
    reg(sg, "ethylene_age_onset",             static_cast<float>(g.ethylene_age_onset),
                                                                                r, 100.0f, 5000.0f, p);
    reg(sg, "ethylene_crowding_rate",         g.ethylene_crowding_rate,         r, 0.01f,  1.0f, p);
    reg(sg, "ethylene_crowding_radius",       g.ethylene_crowding_radius,       r, 0.1f,   3.0f, p);
    reg(sg, "ethylene_diffusion_radius",      g.ethylene_diffusion_radius,      r, 0.1f,   5.0f, p);
    reg(sg, "ethylene_abscission_threshold",  g.ethylene_abscission_threshold,  r, 0.05f,  2.0f, p);
    reg(sg, "ethylene_elongation_inhibition", g.ethylene_elongation_inhibition, r, 0.1f,   5.0f, p);
    reg(sg, "senescence_duration",            static_cast<float>(g.senescence_duration),
                                                                                r, 12.0f,  500.0f, p);

    // --- Stress group (9 genes) ---
    reg(sg, "wood_density",                    g.wood_density,                    r, 10.0f, 200.0f, p);
    reg(sg, "wood_flexibility",                g.wood_flexibility,                r, 0.1f,  1.0f, p);
    reg(sg, "stress_hormone_threshold",        g.stress_hormone_threshold,        r, 0.0f,  0.8f, p);
    reg(sg, "stress_hormone_production_rate",  g.stress_hormone_production_rate,  r, 0.01f, 2.0f, p);
    reg(sg, "stress_hormone_diffusion_rate",   g.stress_hormone_diffusion_rate,   r, 0.01f, 0.5f, p);
    reg(sg, "stress_hormone_decay_rate",       g.stress_hormone_decay_rate,       r, 0.01f, 0.5f, p);
    reg(sg, "stress_thickening_boost",         g.stress_thickening_boost,         r, 0.0f,  5.0f, p);
    reg(sg, "stress_elongation_inhibition",    g.stress_elongation_inhibition,    r, 0.0f,  5.0f, p);
    reg(sg, "stress_gravitropism_boost",       g.stress_gravitropism_boost,       r, 0.0f,  5.0f, p);
    reg(sg, "elastic_recovery_rate",           g.elastic_recovery_rate,           r, 0.0f,  0.02f, p);

    // --- Canalization group (PIN transport params) ---
    reg(sg, "smoothing_rate",         g.smoothing_rate,         r, 0.01f, 0.5f,    p);
    reg(sg, "canalization_weight",    g.canalization_weight,    r, 0.0f,  5.0f,    p);
    reg(sg, "pin_capacity_per_area",  g.pin_capacity_per_area,  r, 50.0f, 2000.0f, p);
    reg(sg, "pin_base_efficiency",    g.pin_base_efficiency,    r, 0.05f, 0.5f,    p);

    // Vascular transport
    reg(sg, "xylem_conductance",           g.xylem_conductance,           r, 1.0f,  500.0f, p);
    reg(sg, "phloem_conductance",          g.phloem_conductance,          r, 1.0f,  50.0f,  p);
    reg(sg, "phloem_reserve_fraction",     g.phloem_reserve_fraction,     r, 0.05f, 0.8f,   p);
    reg(sg, "vascular_radius_threshold",   g.vascular_radius_threshold,   r, 0.001f, 0.05f, p);

    // --- Linkage groups ---
    sg.add_linkage_group({"auxin", {
        "apical_auxin_baseline", "apical_growth_auxin_multiplier",
        "auxin_diffusion_rate",
        "auxin_decay_rate", "auxin_threshold",
        "auxin_shade_boost", "auxin_sugar_half_saturation",
        "auxin_bias",
        "leaf_auxin_baseline", "leaf_growth_auxin_multiplier",
        "stem_auxin_max_boost", "stem_auxin_half_saturation",
        "root_auxin_max_boost", "root_auxin_half_saturation",
        "leaf_auxin_max_boost", "leaf_auxin_half_saturation",
        "apical_auxin_max_boost", "apical_auxin_half_saturation",
        "root_apical_auxin_max_boost", "root_apical_auxin_half_saturation"
    }});

    sg.add_linkage_group({"cytokinin", {
        "cytokinin_production_rate", "cytokinin_diffusion_rate",
        "cytokinin_decay_rate", "cytokinin_threshold", "cytokinin_growth_threshold",
        "cytokinin_bias"
    }});

    sg.add_linkage_group({"transport", {
        "hormone_base_transport", "hormone_transport_scale",
        "sugar_base_transport", "sugar_transport_scale"
    }});

    sg.add_linkage_group({"shoot_growth", {
        "growth_rate", "shoot_plastochron", "max_internode_length", "branch_angle",
        "cambium_responsiveness",
        "internode_elongation_rate", "internode_maturation_ticks"
    }});

    sg.add_linkage_group({"root_growth", {
        "root_growth_rate", "root_plastochron",
        "root_branch_angle", "root_internode_elongation_rate", "root_internode_maturation_ticks",
        "root_gravitropism_strength", "root_gravitropism_depth"
    }});

    sg.add_linkage_group({"geometry", {
        "max_leaf_size", "leaf_growth_rate", "leaf_bud_size", "leaf_petiole_length",
        "leaf_opacity",
        "initial_radius", "root_initial_radius", "tip_offset", "growth_noise",
        "leaf_phototropism_rate", "meristem_gravitropism_rate", "meristem_phototropism_rate",
        "leaf_abscission_ticks", "min_leaf_age_before_abscission"
    }});

    sg.add_linkage_group({"sugar_economy", {
        "sugar_diffusion_rate", "seed_sugar",
        "sugar_storage_density_wood", "sugar_storage_density_leaf",
        "sugar_cap_minimum", "sugar_cap_meristem",
        "xylem_conductance", "phloem_conductance", "phloem_reserve_fraction",
        "vascular_radius_threshold"
    }});

    sg.add_linkage_group({"water_economy", {
        "water_absorption_rate", "transpiration_rate", "photosynthesis_water_ratio",
        "water_storage_density_stem", "water_storage_density_leaf", "water_cap_meristem",
        "water_diffusion_rate", "water_bias", "water_base_transport", "water_transport_scale"
    }});

    sg.add_linkage_group({"gibberellin", {
        "ga_production_rate", "ga_leaf_age_max", "ga_elongation_sensitivity",
        "ga_length_sensitivity", "ga_diffusion_rate", "ga_decay_rate"
    }});

    sg.add_linkage_group({"ethylene", {
        "ethylene_starvation_rate", "ethylene_shade_rate", "ethylene_shade_threshold",
        "ethylene_age_rate", "ethylene_age_onset", "ethylene_crowding_rate",
        "ethylene_crowding_radius", "ethylene_diffusion_radius",
        "ethylene_abscission_threshold", "ethylene_elongation_inhibition",
        "senescence_duration"
    }});

    sg.add_linkage_group({"stress", {
        "wood_density", "wood_flexibility",
        "stress_hormone_threshold",
        "stress_hormone_production_rate", "stress_hormone_diffusion_rate",
        "stress_hormone_decay_rate",
        "stress_thickening_boost", "stress_elongation_inhibition", "stress_gravitropism_boost",
        "elastic_recovery_rate"
    }});

    sg.add_linkage_group({"canalization", {
        "smoothing_rate", "canalization_weight",
        "pin_capacity_per_area", "pin_base_efficiency"
    }});

    return sg;
}

evolve::StructuredGenome to_structured(const Genome& g, float mutation_pct) {
    return build_genome_template(g, mutation_pct);
}

Genome from_structured(const evolve::StructuredGenome& sg) {
    Genome g{};

    // Auxin
    g.apical_auxin_baseline   = sg.get("apical_auxin_baseline");
    g.apical_growth_auxin_multiplier = sg.get("apical_growth_auxin_multiplier");
    g.auxin_diffusion_rate    = sg.get("auxin_diffusion_rate");
    g.auxin_decay_rate        = sg.get("auxin_decay_rate");
    g.auxin_threshold              = sg.get("auxin_threshold");
    g.auxin_shade_boost            = sg.get("auxin_shade_boost");
    g.auxin_sugar_half_saturation  = sg.get("auxin_sugar_half_saturation");
    g.auxin_bias                   = sg.get("auxin_bias");
    g.leaf_auxin_baseline              = sg.get("leaf_auxin_baseline");
    g.leaf_growth_auxin_multiplier     = sg.get("leaf_growth_auxin_multiplier");
    g.stem_auxin_max_boost              = sg.get("stem_auxin_max_boost");
    g.stem_auxin_half_saturation        = sg.get("stem_auxin_half_saturation");
    g.root_auxin_max_boost              = sg.get("root_auxin_max_boost");
    g.root_auxin_half_saturation        = sg.get("root_auxin_half_saturation");
    g.leaf_auxin_max_boost              = sg.get("leaf_auxin_max_boost");
    g.leaf_auxin_half_saturation        = sg.get("leaf_auxin_half_saturation");
    g.apical_auxin_max_boost            = sg.get("apical_auxin_max_boost");
    g.apical_auxin_half_saturation      = sg.get("apical_auxin_half_saturation");
    g.root_apical_auxin_max_boost       = sg.get("root_apical_auxin_max_boost");
    g.root_apical_auxin_half_saturation = sg.get("root_apical_auxin_half_saturation");

    // Cytokinin
    g.cytokinin_production_rate   = sg.get("cytokinin_production_rate");
    g.cytokinin_diffusion_rate    = sg.get("cytokinin_diffusion_rate");
    g.cytokinin_decay_rate        = sg.get("cytokinin_decay_rate");
    g.cytokinin_threshold         = sg.get("cytokinin_threshold");
    g.cytokinin_growth_threshold  = sg.get("cytokinin_growth_threshold");
    g.cytokinin_bias              = sg.get("cytokinin_bias");

    // Transport capacity
    g.hormone_base_transport      = sg.get("hormone_base_transport");
    g.hormone_transport_scale     = sg.get("hormone_transport_scale");
    g.sugar_base_transport        = sg.get("sugar_base_transport");
    g.sugar_transport_scale       = sg.get("sugar_transport_scale");

    // Shoot growth
    g.growth_rate                = sg.get("growth_rate");
    g.shoot_plastochron          = static_cast<uint32_t>(sg.get("shoot_plastochron"));
    g.max_internode_length       = sg.get("max_internode_length");
    g.branch_angle               = sg.get("branch_angle");
    // Backwards compatibility: old genomes have thickening_rate but not cambium_responsiveness.
    // Graceful default if the gene is missing from a saved genome file.
    g.cambium_responsiveness     = sg.has_gene("cambium_responsiveness")
                                       ? sg.get("cambium_responsiveness")
                                       : 0.00002f;
    g.internode_elongation_rate  = sg.get("internode_elongation_rate");
    g.internode_maturation_ticks = static_cast<uint32_t>(sg.get("internode_maturation_ticks"));

    // Root growth
    g.root_growth_rate                  = sg.get("root_growth_rate");
    g.root_plastochron                  = static_cast<uint32_t>(sg.get("root_plastochron"));
    g.root_branch_angle                 = sg.get("root_branch_angle");
    g.root_internode_elongation_rate    = sg.get("root_internode_elongation_rate");
    g.root_internode_maturation_ticks   = static_cast<uint32_t>(sg.get("root_internode_maturation_ticks"));
    g.root_gravitropism_strength        = sg.get("root_gravitropism_strength");
    g.root_gravitropism_depth           = sg.get("root_gravitropism_depth");

    // Geometry
    g.max_leaf_size           = sg.get("max_leaf_size");
    g.leaf_growth_rate        = sg.get("leaf_growth_rate");
    g.leaf_bud_size           = sg.get("leaf_bud_size");
    g.leaf_petiole_length     = sg.get("leaf_petiole_length");
    g.leaf_opacity            = sg.get("leaf_opacity");
    g.initial_radius          = sg.get("initial_radius");
    g.root_initial_radius     = sg.get("root_initial_radius");
    g.tip_offset              = sg.get("tip_offset");
    g.growth_noise            = sg.get("growth_noise");
    g.leaf_phototropism_rate       = sg.get("leaf_phototropism_rate");
    g.meristem_gravitropism_rate   = sg.get("meristem_gravitropism_rate");
    g.meristem_phototropism_rate   = sg.get("meristem_phototropism_rate");
    g.leaf_abscission_ticks  = static_cast<uint32_t>(sg.get("leaf_abscission_ticks"));
    g.min_leaf_age_before_abscission = static_cast<uint32_t>(sg.get("min_leaf_age_before_abscission"));

    // Sugar economy
    g.sugar_diffusion_rate        = sg.get("sugar_diffusion_rate");
    g.seed_sugar                  = sg.get("seed_sugar");
    g.sugar_storage_density_wood  = sg.get("sugar_storage_density_wood");
    g.sugar_storage_density_leaf  = sg.get("sugar_storage_density_leaf");
    g.sugar_cap_minimum           = sg.get("sugar_cap_minimum");
    g.sugar_cap_meristem          = sg.get("sugar_cap_meristem");
    // Water economy
    g.water_absorption_rate       = sg.get("water_absorption_rate");
    g.transpiration_rate          = sg.get("transpiration_rate");
    g.photosynthesis_water_ratio  = sg.get("photosynthesis_water_ratio");
    g.water_storage_density_stem  = sg.get("water_storage_density_stem");
    g.water_storage_density_leaf  = sg.get("water_storage_density_leaf");
    g.water_cap_meristem          = sg.get("water_cap_meristem");
    g.water_diffusion_rate        = sg.get("water_diffusion_rate");
    g.water_bias                  = sg.get("water_bias");
    g.water_base_transport        = sg.get("water_base_transport");
    g.water_transport_scale       = sg.get("water_transport_scale");
    // Gibberellin
    g.ga_production_rate        = sg.get("ga_production_rate");
    g.ga_leaf_age_max           = static_cast<uint32_t>(sg.get("ga_leaf_age_max"));
    g.ga_elongation_sensitivity = sg.get("ga_elongation_sensitivity");
    g.ga_length_sensitivity     = sg.get("ga_length_sensitivity");
    g.ga_diffusion_rate         = sg.get("ga_diffusion_rate");
    g.ga_decay_rate             = sg.get("ga_decay_rate");

    // Ethylene
    g.ethylene_starvation_rate       = sg.get("ethylene_starvation_rate");
    g.ethylene_shade_rate            = sg.get("ethylene_shade_rate");
    g.ethylene_shade_threshold       = sg.get("ethylene_shade_threshold");
    g.ethylene_age_rate              = sg.get("ethylene_age_rate");
    g.ethylene_age_onset             = static_cast<uint32_t>(sg.get("ethylene_age_onset"));
    g.ethylene_crowding_rate         = sg.get("ethylene_crowding_rate");
    g.ethylene_crowding_radius       = sg.get("ethylene_crowding_radius");
    g.ethylene_diffusion_radius      = sg.get("ethylene_diffusion_radius");
    g.ethylene_abscission_threshold  = sg.get("ethylene_abscission_threshold");
    g.ethylene_elongation_inhibition = sg.get("ethylene_elongation_inhibition");
    g.senescence_duration            = static_cast<uint32_t>(sg.get("senescence_duration"));

    // Stress
    g.wood_density                    = sg.get("wood_density");
    g.wood_flexibility                = sg.get("wood_flexibility");
    g.stress_hormone_threshold        = sg.get("stress_hormone_threshold");
    g.stress_hormone_production_rate  = sg.get("stress_hormone_production_rate");
    g.stress_hormone_diffusion_rate   = sg.get("stress_hormone_diffusion_rate");
    g.stress_hormone_decay_rate       = sg.get("stress_hormone_decay_rate");
    g.stress_thickening_boost         = sg.get("stress_thickening_boost");
    g.stress_elongation_inhibition    = sg.get("stress_elongation_inhibition");
    g.stress_gravitropism_boost       = sg.get("stress_gravitropism_boost");
    g.elastic_recovery_rate           = sg.get("elastic_recovery_rate");

    // Canalization / PIN transport
    g.smoothing_rate        = sg.get("smoothing_rate");
    g.canalization_weight   = sg.get("canalization_weight");
    g.pin_capacity_per_area = sg.get("pin_capacity_per_area");
    g.pin_base_efficiency   = sg.get("pin_base_efficiency");

    // Vascular transport
    g.xylem_conductance          = sg.get("xylem_conductance");
    g.phloem_conductance         = sg.get("phloem_conductance");
    g.phloem_reserve_fraction    = sg.get("phloem_reserve_fraction");
    g.vascular_radius_threshold  = sg.has_gene("vascular_radius_threshold")
                                       ? sg.get("vascular_radius_threshold")
                                       : 0.01f;

    return g;
}

} // namespace botany
