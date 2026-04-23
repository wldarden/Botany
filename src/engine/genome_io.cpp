#include "engine/genome_io.h"

#include <fstream>
#include <map>
#include <string>

namespace botany {

Genome load_genome_file(const std::string& path) {
    Genome g = default_genome();
    std::ifstream file(path);
    if (!file) return g;

    std::map<std::string, std::string> fields;
    std::string line;
    while (std::getline(file, line)) {
        auto eq = line.find('=');
        if (eq != std::string::npos) {
            fields[line.substr(0, eq)] = line.substr(eq + 1);
        }
    }

    auto get_f = [&](const char* name, float& out) {
        auto it = fields.find(name);
        if (it != fields.end()) out = std::stof(it->second);
    };
    auto get_u = [&](const char* name, uint32_t& out) {
        auto it = fields.find(name);
        if (it != fields.end()) out = static_cast<uint32_t>(std::stoul(it->second));
    };

    get_f("apical_auxin_baseline", g.apical_auxin_baseline);
    get_f("auxin_diffusion_rate", g.auxin_diffusion_rate);
    get_f("auxin_decay_rate", g.auxin_decay_rate);
    get_f("auxin_threshold", g.auxin_threshold);
    get_f("cytokinin_production_rate", g.cytokinin_production_rate);
    get_f("cytokinin_diffusion_rate", g.cytokinin_diffusion_rate);
    get_f("cytokinin_decay_rate", g.cytokinin_decay_rate);
    get_f("cytokinin_threshold", g.cytokinin_threshold);
    get_f("cytokinin_growth_threshold", g.cytokinin_growth_threshold);
    get_f("growth_rate", g.growth_rate);
    get_u("shoot_plastochron", g.shoot_plastochron);
    get_f("branch_angle", g.branch_angle);
    get_f("cambium_responsiveness", g.cambium_responsiveness);
    get_f("internode_elongation_rate", g.internode_elongation_rate);
    get_f("max_internode_length", g.max_internode_length);
    get_u("internode_maturation_ticks", g.internode_maturation_ticks);
    get_f("root_growth_rate", g.root_growth_rate);
    get_u("root_plastochron", g.root_plastochron);
    get_f("root_branch_angle", g.root_branch_angle);
    get_f("root_internode_elongation_rate", g.root_internode_elongation_rate);
    get_u("root_internode_maturation_ticks", g.root_internode_maturation_ticks);
    get_f("root_gravitropism_strength", g.root_gravitropism_strength);
    get_f("root_gravitropism_depth", g.root_gravitropism_depth);
    get_f("max_leaf_size", g.max_leaf_size);
    get_f("leaf_growth_rate", g.leaf_growth_rate);
    get_f("leaf_bud_size", g.leaf_bud_size);
    get_f("leaf_petiole_length", g.leaf_petiole_length);
    get_f("leaf_opacity", g.leaf_opacity);
    get_f("initial_radius", g.initial_radius);
    get_f("root_initial_radius", g.root_initial_radius);
    get_f("tip_offset", g.tip_offset);
    get_f("growth_noise", g.growth_noise);
    get_f("leaf_phototropism_rate", g.leaf_phototropism_rate);
    get_f("sugar_diffusion_rate", g.sugar_diffusion_rate);
    get_f("seed_sugar", g.seed_sugar);
    get_f("sugar_storage_density_wood", g.sugar_storage_density_wood);
    get_f("sugar_storage_density_leaf", g.sugar_storage_density_leaf);
    get_f("sugar_cap_minimum", g.sugar_cap_minimum);
    get_f("sugar_cap_meristem", g.sugar_cap_meristem);
    get_f("stem_photosynthesis_rate", g.stem_photosynthesis_rate);
    get_f("stem_green_radius_threshold", g.stem_green_radius_threshold);
    get_f("ga_production_rate", g.ga_production_rate);
    get_u("ga_leaf_age_max", g.ga_leaf_age_max);
    get_f("ga_elongation_sensitivity", g.ga_elongation_sensitivity);
    get_f("ga_length_sensitivity", g.ga_length_sensitivity);
    get_f("ga_diffusion_rate", g.ga_diffusion_rate);
    get_f("ga_decay_rate", g.ga_decay_rate);
    get_f("ethylene_starvation_rate", g.ethylene_starvation_rate);
    get_f("ethylene_shade_rate", g.ethylene_shade_rate);
    get_f("ethylene_shade_threshold", g.ethylene_shade_threshold);
    get_f("ethylene_age_rate", g.ethylene_age_rate);
    get_u("ethylene_age_onset", g.ethylene_age_onset);
    get_f("ethylene_crowding_rate", g.ethylene_crowding_rate);
    get_f("ethylene_crowding_radius", g.ethylene_crowding_radius);
    get_f("ethylene_diffusion_radius", g.ethylene_diffusion_radius);
    get_f("ethylene_abscission_threshold", g.ethylene_abscission_threshold);
    get_f("ethylene_elongation_inhibition", g.ethylene_elongation_inhibition);
    get_u("senescence_duration", g.senescence_duration);
    get_f("wood_density", g.wood_density);
    get_f("wood_flexibility", g.wood_flexibility);
    get_f("stress_hormone_threshold", g.stress_hormone_threshold);
    get_f("stress_hormone_production_rate", g.stress_hormone_production_rate);
    get_f("stress_hormone_diffusion_rate", g.stress_hormone_diffusion_rate);
    get_f("stress_hormone_decay_rate", g.stress_hormone_decay_rate);
    get_f("stress_thickening_boost", g.stress_thickening_boost);
    get_f("stress_elongation_inhibition", g.stress_elongation_inhibition);
    get_f("stress_gravitropism_boost", g.stress_gravitropism_boost);
    get_f("vascular_radius_threshold", g.vascular_radius_threshold);

    return g;
}

bool save_genome_file(const Genome& g, const std::string& path) {
    std::ofstream file(path);
    if (!file) return false;

    file << "apical_auxin_baseline=" << g.apical_auxin_baseline << "\n";
    file << "auxin_diffusion_rate=" << g.auxin_diffusion_rate << "\n";
    file << "auxin_decay_rate=" << g.auxin_decay_rate << "\n";
    file << "auxin_threshold=" << g.auxin_threshold << "\n";
    file << "cytokinin_production_rate=" << g.cytokinin_production_rate << "\n";
    file << "cytokinin_diffusion_rate=" << g.cytokinin_diffusion_rate << "\n";
    file << "cytokinin_decay_rate=" << g.cytokinin_decay_rate << "\n";
    file << "cytokinin_threshold=" << g.cytokinin_threshold << "\n";
    file << "cytokinin_growth_threshold=" << g.cytokinin_growth_threshold << "\n";
    file << "growth_rate=" << g.growth_rate << "\n";
    file << "shoot_plastochron=" << g.shoot_plastochron << "\n";
    file << "branch_angle=" << g.branch_angle << "\n";
    file << "cambium_responsiveness=" << g.cambium_responsiveness << "\n";
    file << "internode_elongation_rate=" << g.internode_elongation_rate << "\n";
    file << "max_internode_length=" << g.max_internode_length << "\n";
    file << "internode_maturation_ticks=" << g.internode_maturation_ticks << "\n";
    file << "root_growth_rate=" << g.root_growth_rate << "\n";
    file << "root_plastochron=" << g.root_plastochron << "\n";
    file << "root_branch_angle=" << g.root_branch_angle << "\n";
    file << "root_internode_elongation_rate=" << g.root_internode_elongation_rate << "\n";
    file << "root_internode_maturation_ticks=" << g.root_internode_maturation_ticks << "\n";
    file << "root_gravitropism_strength=" << g.root_gravitropism_strength << "\n";
    file << "root_gravitropism_depth=" << g.root_gravitropism_depth << "\n";
    file << "max_leaf_size=" << g.max_leaf_size << "\n";
    file << "leaf_growth_rate=" << g.leaf_growth_rate << "\n";
    file << "leaf_bud_size=" << g.leaf_bud_size << "\n";
    file << "leaf_petiole_length=" << g.leaf_petiole_length << "\n";
    file << "leaf_opacity=" << g.leaf_opacity << "\n";
    file << "initial_radius=" << g.initial_radius << "\n";
    file << "root_initial_radius=" << g.root_initial_radius << "\n";
    file << "tip_offset=" << g.tip_offset << "\n";
    file << "growth_noise=" << g.growth_noise << "\n";
    file << "leaf_phototropism_rate=" << g.leaf_phototropism_rate << "\n";
    file << "sugar_diffusion_rate=" << g.sugar_diffusion_rate << "\n";
    file << "seed_sugar=" << g.seed_sugar << "\n";
    file << "sugar_storage_density_wood=" << g.sugar_storage_density_wood << "\n";
    file << "sugar_storage_density_leaf=" << g.sugar_storage_density_leaf << "\n";
    file << "sugar_cap_minimum=" << g.sugar_cap_minimum << "\n";
    file << "sugar_cap_meristem=" << g.sugar_cap_meristem << "\n";
    file << "stem_photosynthesis_rate=" << g.stem_photosynthesis_rate << "\n";
    file << "stem_green_radius_threshold=" << g.stem_green_radius_threshold << "\n";
    file << "ga_production_rate=" << g.ga_production_rate << "\n";
    file << "ga_leaf_age_max=" << g.ga_leaf_age_max << "\n";
    file << "ga_elongation_sensitivity=" << g.ga_elongation_sensitivity << "\n";
    file << "ga_length_sensitivity=" << g.ga_length_sensitivity << "\n";
    file << "ga_diffusion_rate=" << g.ga_diffusion_rate << "\n";
    file << "ga_decay_rate=" << g.ga_decay_rate << "\n";
    file << "ethylene_starvation_rate=" << g.ethylene_starvation_rate << "\n";
    file << "ethylene_shade_rate=" << g.ethylene_shade_rate << "\n";
    file << "ethylene_shade_threshold=" << g.ethylene_shade_threshold << "\n";
    file << "ethylene_age_rate=" << g.ethylene_age_rate << "\n";
    file << "ethylene_age_onset=" << g.ethylene_age_onset << "\n";
    file << "ethylene_crowding_rate=" << g.ethylene_crowding_rate << "\n";
    file << "ethylene_crowding_radius=" << g.ethylene_crowding_radius << "\n";
    file << "ethylene_diffusion_radius=" << g.ethylene_diffusion_radius << "\n";
    file << "ethylene_abscission_threshold=" << g.ethylene_abscission_threshold << "\n";
    file << "ethylene_elongation_inhibition=" << g.ethylene_elongation_inhibition << "\n";
    file << "senescence_duration=" << g.senescence_duration << "\n";
    file << "wood_density=" << g.wood_density << "\n";
    file << "wood_flexibility=" << g.wood_flexibility << "\n";
    file << "stress_hormone_threshold=" << g.stress_hormone_threshold << "\n";
    file << "stress_hormone_production_rate=" << g.stress_hormone_production_rate << "\n";
    file << "stress_hormone_diffusion_rate=" << g.stress_hormone_diffusion_rate << "\n";
    file << "stress_hormone_decay_rate=" << g.stress_hormone_decay_rate << "\n";
    file << "stress_thickening_boost=" << g.stress_thickening_boost << "\n";
    file << "stress_elongation_inhibition=" << g.stress_elongation_inhibition << "\n";
    file << "stress_gravitropism_boost=" << g.stress_gravitropism_boost << "\n";
    file << "vascular_radius_threshold=" << g.vascular_radius_threshold << "\n";

    return static_cast<bool>(file);
}

} // namespace botany
