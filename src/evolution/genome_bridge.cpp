#include "evolution/genome_bridge.h"

namespace botany {

static void reg(evolve::StructuredGenome& sg, const std::string& tag, float val,
                float rate, float strength, float min_val, float max_val) {
    sg.add_gene({tag, {val}, {rate, strength, min_val, max_val}});
}

evolve::StructuredGenome build_genome_template(const Genome& g) {
    evolve::StructuredGenome sg;

    // --- Auxin group (5 genes) ---
    reg(sg, "auxin_production_rate",   g.auxin_production_rate,   0.1f, 0.05f, 0.01f, 2.0f);
    reg(sg, "auxin_transport_rate",    g.auxin_transport_rate,    0.1f, 0.05f, 0.01f, 2.0f);
    reg(sg, "auxin_directional_bias",  g.auxin_directional_bias,  0.1f, 0.1f, -1.0f, 1.0f);
    reg(sg, "auxin_decay_rate",        g.auxin_decay_rate,        0.1f, 0.02f, 0.001f, 0.5f);
    reg(sg, "auxin_threshold",         g.auxin_threshold,         0.1f, 0.03f, 0.01f, 1.0f);

    // --- Cytokinin group (5 genes) ---
    reg(sg, "cytokinin_production_rate",  g.cytokinin_production_rate,  0.1f, 0.05f, 0.01f, 2.0f);
    reg(sg, "cytokinin_transport_rate",   g.cytokinin_transport_rate,   0.1f, 0.05f, 0.01f, 2.0f);
    reg(sg, "cytokinin_directional_bias", g.cytokinin_directional_bias, 0.1f, 0.1f, -1.0f, 1.0f);
    reg(sg, "cytokinin_decay_rate",       g.cytokinin_decay_rate,       0.1f, 0.02f, 0.001f, 0.5f);
    reg(sg, "cytokinin_threshold",        g.cytokinin_threshold,        0.1f, 0.03f, 0.01f, 1.0f);

    // --- Shoot growth group (7 genes) ---
    reg(sg, "growth_rate",                 g.growth_rate,                 0.1f, 0.001f,  0.001f, 0.05f);
    reg(sg, "max_internode_length",        g.max_internode_length,        0.1f, 0.05f,   0.05f,  3.0f);
    reg(sg, "min_internode_length",        g.min_internode_length,        0.1f, 0.03f,   0.01f,  1.0f);
    reg(sg, "branch_angle",                g.branch_angle,                0.1f, 0.1f,    0.05f,  1.57f);
    reg(sg, "thickening_rate",             g.thickening_rate,             0.1f, 0.00002f, 0.00001f, 0.001f);
    reg(sg, "internode_elongation_rate",   g.internode_elongation_rate,   0.1f, 0.001f,  0.0005f, 0.02f);
    reg(sg, "internode_maturation_ticks",  static_cast<float>(g.internode_maturation_ticks),
                                                                           0.1f, 10.0f,   12.0f,  500.0f);

    // --- Root growth group (8 genes) ---
    reg(sg, "root_growth_rate",                  g.root_growth_rate,                  0.1f, 0.001f, 0.001f, 0.05f);
    reg(sg, "root_max_internode_length",         g.root_max_internode_length,         0.1f, 0.05f,  0.05f,  3.0f);
    reg(sg, "root_min_internode_length",         g.root_min_internode_length,         0.1f, 0.03f,  0.01f,  1.0f);
    reg(sg, "root_branch_angle",                 g.root_branch_angle,                 0.1f, 0.05f,  0.05f,  1.57f);
    reg(sg, "root_internode_elongation_rate",    g.root_internode_elongation_rate,    0.1f, 0.001f, 0.0005f, 0.02f);
    reg(sg, "root_internode_maturation_ticks",   static_cast<float>(g.root_internode_maturation_ticks),
                                                                                       0.1f, 10.0f,  12.0f,  500.0f);
    reg(sg, "root_gravitropism_strength",        g.root_gravitropism_strength,        0.1f, 0.2f,   0.1f,   10.0f);
    reg(sg, "root_gravitropism_depth",           g.root_gravitropism_depth,           0.1f, 0.1f,   0.1f,   5.0f);

    // --- Geometry group (8 genes) ---
    reg(sg, "max_leaf_size",            g.max_leaf_size,            0.1f, 0.03f,   0.05f,  1.0f);
    reg(sg, "leaf_growth_rate",         g.leaf_growth_rate,         0.1f, 0.0005f, 0.0001f, 0.01f);
    reg(sg, "leaf_bud_size",            g.leaf_bud_size,            0.1f, 0.005f,  0.005f, 0.1f);
    reg(sg, "initial_radius",           g.initial_radius,           0.1f, 0.01f,   0.01f,  0.2f);
    reg(sg, "root_initial_radius",      g.root_initial_radius,      0.1f, 0.005f,  0.005f, 0.1f);
    reg(sg, "tip_offset",               g.tip_offset,               0.1f, 0.005f,  0.001f, 0.1f);
    reg(sg, "growth_noise",             g.growth_noise,             0.1f, 0.05f,   0.01f,  0.8f);
    reg(sg, "leaf_phototropism_rate",   g.leaf_phototropism_rate,   0.1f, 0.005f,  0.001f, 0.1f);

    // --- Sugar economy group (16 genes) ---
    reg(sg, "sugar_production_rate",       g.sugar_production_rate,       0.1f, 0.002f, 0.001f, 0.1f);
    reg(sg, "sugar_transport_conductance", g.sugar_transport_conductance, 0.1f, 2.0f,   1.0f,   100.0f);
    reg(sg, "sugar_maintenance_leaf",      g.sugar_maintenance_leaf,      0.1f, 0.002f, 0.001f, 0.1f);
    reg(sg, "sugar_maintenance_stem",      g.sugar_maintenance_stem,      0.1f, 0.005f, 0.001f, 0.2f);
    reg(sg, "sugar_maintenance_root",      g.sugar_maintenance_root,      0.1f, 0.02f,  0.01f,  0.5f);
    reg(sg, "sugar_maintenance_meristem",  g.sugar_maintenance_meristem,  0.1f, 0.001f, 0.0001f, 0.01f);
    reg(sg, "seed_sugar",                  g.seed_sugar,                  0.1f, 5.0f,   5.0f,   200.0f);
    reg(sg, "sugar_storage_density_wood",  g.sugar_storage_density_wood,  0.1f, 50.0f,  50.0f,  2000.0f);
    reg(sg, "sugar_storage_density_leaf",  g.sugar_storage_density_leaf,  0.1f, 0.1f,   0.05f,  5.0f);
    reg(sg, "sugar_cap_minimum",           g.sugar_cap_minimum,           0.1f, 0.01f,  0.005f, 0.5f);
    reg(sg, "sugar_cap_meristem",          g.sugar_cap_meristem,          0.1f, 0.5f,   0.1f,   10.0f);
    reg(sg, "sugar_save_shoot",            g.sugar_save_shoot,         0.1f, 0.005f, 0.001f, 0.2f);
    reg(sg, "sugar_save_root",             g.sugar_save_root,          0.1f, 0.003f, 0.001f, 0.1f);
    reg(sg, "sugar_save_stem",             g.sugar_save_stem,          0.1f, 0.005f, 0.001f, 0.2f);
    reg(sg, "sugar_activation_shoot",      g.sugar_activation_shoot,      0.1f, 0.1f,   0.05f,  5.0f);
    reg(sg, "sugar_activation_root",       g.sugar_activation_root,       0.1f, 0.05f,  0.05f,  3.0f);

    // --- Gibberellin group (7 genes) ---
    reg(sg, "ga_production_rate",        g.ga_production_rate,        0.1f, 0.05f, 0.01f,  2.0f);
    reg(sg, "ga_leaf_age_max",           static_cast<float>(g.ga_leaf_age_max),
                                                                      0.1f, 20.0f, 24.0f,  2000.0f);
    reg(sg, "ga_elongation_sensitivity", g.ga_elongation_sensitivity, 0.1f, 0.2f,  0.1f,   5.0f);
    reg(sg, "ga_length_sensitivity",     g.ga_length_sensitivity,     0.1f, 0.2f,  0.1f,   5.0f);
    reg(sg, "ga_transport_rate",         g.ga_transport_rate,         0.1f, 0.05f, 0.01f,  1.0f);
    reg(sg, "ga_directional_bias",       g.ga_directional_bias,       0.1f, 0.1f, -1.0f,   1.0f);
    reg(sg, "ga_decay_rate",             g.ga_decay_rate,             0.1f, 0.03f, 0.01f,  0.5f);

    // --- Ethylene group (11 genes) ---
    reg(sg, "ethylene_starvation_rate",       g.ethylene_starvation_rate,       0.1f, 0.05f, 0.01f,  2.0f);
    reg(sg, "ethylene_shade_rate",            g.ethylene_shade_rate,            0.1f, 0.05f, 0.01f,  2.0f);
    reg(sg, "ethylene_shade_threshold",       g.ethylene_shade_threshold,       0.1f, 0.05f, 0.05f,  1.0f);
    reg(sg, "ethylene_age_rate",              g.ethylene_age_rate,              0.1f, 0.01f, 0.001f, 0.5f);
    reg(sg, "ethylene_age_onset",             static_cast<float>(g.ethylene_age_onset),
                                                                                0.1f, 50.0f, 100.0f, 5000.0f);
    reg(sg, "ethylene_crowding_rate",         g.ethylene_crowding_rate,         0.1f, 0.02f, 0.01f,  1.0f);
    reg(sg, "ethylene_crowding_radius",       g.ethylene_crowding_radius,       0.1f, 0.1f,  0.1f,   3.0f);
    reg(sg, "ethylene_diffusion_radius",      g.ethylene_diffusion_radius,      0.1f, 0.2f,  0.1f,   5.0f);
    reg(sg, "ethylene_abscission_threshold",  g.ethylene_abscission_threshold,  0.1f, 0.1f,  0.05f,  2.0f);
    reg(sg, "ethylene_elongation_inhibition", g.ethylene_elongation_inhibition, 0.1f, 0.2f,  0.1f,   5.0f);
    reg(sg, "senescence_duration",            static_cast<float>(g.senescence_duration),
                                                                                0.1f, 10.0f, 12.0f,  500.0f);

    // --- Linkage groups ---
    sg.add_linkage_group({"auxin", {
        "auxin_production_rate", "auxin_transport_rate", "auxin_directional_bias",
        "auxin_decay_rate", "auxin_threshold"
    }});

    sg.add_linkage_group({"cytokinin", {
        "cytokinin_production_rate", "cytokinin_transport_rate", "cytokinin_directional_bias",
        "cytokinin_decay_rate", "cytokinin_threshold"
    }});

    sg.add_linkage_group({"shoot_growth", {
        "growth_rate", "max_internode_length", "min_internode_length", "branch_angle",
        "thickening_rate", "internode_elongation_rate", "internode_maturation_ticks"
    }});

    sg.add_linkage_group({"root_growth", {
        "root_growth_rate", "root_max_internode_length", "root_min_internode_length",
        "root_branch_angle", "root_internode_elongation_rate", "root_internode_maturation_ticks",
        "root_gravitropism_strength", "root_gravitropism_depth"
    }});

    sg.add_linkage_group({"geometry", {
        "max_leaf_size", "leaf_growth_rate", "leaf_bud_size", "initial_radius",
        "root_initial_radius", "tip_offset", "growth_noise", "leaf_phototropism_rate"
    }});

    sg.add_linkage_group({"sugar_economy", {
        "sugar_production_rate", "sugar_transport_conductance",
        "sugar_maintenance_leaf", "sugar_maintenance_stem",
        "sugar_maintenance_root", "sugar_maintenance_meristem",
        "seed_sugar", "sugar_storage_density_wood", "sugar_storage_density_leaf",
        "sugar_cap_minimum", "sugar_cap_meristem",
        "sugar_save_shoot", "sugar_save_root", "sugar_save_stem",
        "sugar_activation_shoot", "sugar_activation_root"
    }});

    sg.add_linkage_group({"gibberellin", {
        "ga_production_rate", "ga_leaf_age_max", "ga_elongation_sensitivity",
        "ga_length_sensitivity", "ga_transport_rate", "ga_directional_bias", "ga_decay_rate"
    }});

    sg.add_linkage_group({"ethylene", {
        "ethylene_starvation_rate", "ethylene_shade_rate", "ethylene_shade_threshold",
        "ethylene_age_rate", "ethylene_age_onset", "ethylene_crowding_rate",
        "ethylene_crowding_radius", "ethylene_diffusion_radius",
        "ethylene_abscission_threshold", "ethylene_elongation_inhibition",
        "senescence_duration"
    }});

    return sg;
}

evolve::StructuredGenome to_structured(const Genome& g) {
    return build_genome_template(g);
}

Genome from_structured(const evolve::StructuredGenome& sg) {
    Genome g{};

    // Auxin
    g.auxin_production_rate   = sg.get("auxin_production_rate");
    g.auxin_transport_rate    = sg.get("auxin_transport_rate");
    g.auxin_directional_bias  = sg.get("auxin_directional_bias");
    g.auxin_decay_rate        = sg.get("auxin_decay_rate");
    g.auxin_threshold         = sg.get("auxin_threshold");

    // Cytokinin
    g.cytokinin_production_rate   = sg.get("cytokinin_production_rate");
    g.cytokinin_transport_rate    = sg.get("cytokinin_transport_rate");
    g.cytokinin_directional_bias  = sg.get("cytokinin_directional_bias");
    g.cytokinin_decay_rate        = sg.get("cytokinin_decay_rate");
    g.cytokinin_threshold         = sg.get("cytokinin_threshold");

    // Shoot growth
    g.growth_rate                = sg.get("growth_rate");
    g.max_internode_length       = sg.get("max_internode_length");
    g.min_internode_length       = sg.get("min_internode_length");
    g.branch_angle               = sg.get("branch_angle");
    g.thickening_rate            = sg.get("thickening_rate");
    g.internode_elongation_rate  = sg.get("internode_elongation_rate");
    g.internode_maturation_ticks = static_cast<uint32_t>(sg.get("internode_maturation_ticks"));

    // Root growth
    g.root_growth_rate                  = sg.get("root_growth_rate");
    g.root_max_internode_length         = sg.get("root_max_internode_length");
    g.root_min_internode_length         = sg.get("root_min_internode_length");
    g.root_branch_angle                 = sg.get("root_branch_angle");
    g.root_internode_elongation_rate    = sg.get("root_internode_elongation_rate");
    g.root_internode_maturation_ticks   = static_cast<uint32_t>(sg.get("root_internode_maturation_ticks"));
    g.root_gravitropism_strength        = sg.get("root_gravitropism_strength");
    g.root_gravitropism_depth           = sg.get("root_gravitropism_depth");

    // Geometry
    g.max_leaf_size           = sg.get("max_leaf_size");
    g.leaf_growth_rate        = sg.get("leaf_growth_rate");
    g.leaf_bud_size           = sg.get("leaf_bud_size");
    g.initial_radius          = sg.get("initial_radius");
    g.root_initial_radius     = sg.get("root_initial_radius");
    g.tip_offset              = sg.get("tip_offset");
    g.growth_noise            = sg.get("growth_noise");
    g.leaf_phototropism_rate  = sg.get("leaf_phototropism_rate");

    // Sugar economy
    g.sugar_production_rate       = sg.get("sugar_production_rate");
    g.sugar_transport_conductance = sg.get("sugar_transport_conductance");
    g.sugar_maintenance_leaf      = sg.get("sugar_maintenance_leaf");
    g.sugar_maintenance_stem      = sg.get("sugar_maintenance_stem");
    g.sugar_maintenance_root      = sg.get("sugar_maintenance_root");
    g.sugar_maintenance_meristem  = sg.get("sugar_maintenance_meristem");
    g.seed_sugar                  = sg.get("seed_sugar");
    g.sugar_storage_density_wood  = sg.get("sugar_storage_density_wood");
    g.sugar_storage_density_leaf  = sg.get("sugar_storage_density_leaf");
    g.sugar_cap_minimum           = sg.get("sugar_cap_minimum");
    g.sugar_cap_meristem          = sg.get("sugar_cap_meristem");
    g.sugar_save_shoot         = sg.get("sugar_save_shoot");
    g.sugar_save_root          = sg.get("sugar_save_root");
    g.sugar_save_stem          = sg.get("sugar_save_stem");
    g.sugar_activation_shoot      = sg.get("sugar_activation_shoot");
    g.sugar_activation_root       = sg.get("sugar_activation_root");

    // Gibberellin
    g.ga_production_rate        = sg.get("ga_production_rate");
    g.ga_leaf_age_max           = static_cast<uint32_t>(sg.get("ga_leaf_age_max"));
    g.ga_elongation_sensitivity = sg.get("ga_elongation_sensitivity");
    g.ga_length_sensitivity     = sg.get("ga_length_sensitivity");
    g.ga_transport_rate         = sg.get("ga_transport_rate");
    g.ga_directional_bias       = sg.get("ga_directional_bias");
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

    return g;
}

} // namespace botany
