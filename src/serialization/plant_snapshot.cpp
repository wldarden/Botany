#include "serialization/plant_snapshot.h"
#include "engine/plant.h"
#include "engine/node/node.h"
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <climits>

namespace botany {

// -- Binary helpers (file-local). Values are written host-endian; this format
// is not expected to be portable across differently-endian machines, same as
// the existing serializer.cpp convention.
namespace {
template<typename T>
void write_val(std::ostream& out, const T& v) {
    out.write(reinterpret_cast<const char*>(&v), sizeof(T));
}

template<typename T>
T read_val(std::istream& in) {
    T v;
    in.read(reinterpret_cast<char*>(&v), sizeof(T));
    if (!in) throw std::runtime_error("plant_snapshot: unexpected EOF");
    return v;
}
} // namespace

void plant_snapshot_write_magic(std::ostream& out) {
    out.write(PLANT_SNAPSHOT_MAGIC, 4);
}

bool plant_snapshot_check_magic(std::istream& in) {
    char buf[4];
    in.read(buf, 4);
    if (!in) return false;
    return std::memcmp(buf, PLANT_SNAPSHOT_MAGIC, 4) == 0;
}

void write_genome_binary(std::ostream& out, const Genome& g) {
    write_val(out, g.apical_auxin_baseline);
    write_val(out, g.apical_growth_auxin_multiplier);
    write_val(out, g.auxin_diffusion_rate);
    write_val(out, g.auxin_decay_rate);
    write_val(out, g.auxin_threshold);
    write_val(out, g.auxin_shade_boost);
    write_val(out, g.auxin_sugar_half_saturation);
    write_val(out, g.auxin_water_half_saturation);
    write_val(out, g.cytokinin_sugar_half_saturation);
    write_val(out, g.cytokinin_water_half_saturation);
    write_val(out, g.quiescence_threshold);
    write_val(out, g.meristem_dormancy_death_ticks);
    write_val(out, g.starvation_ticks_max_stem);
    write_val(out, g.starvation_ticks_max_root);
    write_val(out, g.auxin_bias);
    write_val(out, g.leaf_auxin_baseline);
    write_val(out, g.leaf_growth_auxin_multiplier);
    write_val(out, g.stem_auxin_max_boost);
    write_val(out, g.stem_auxin_half_saturation);
    write_val(out, g.root_auxin_max_boost);
    write_val(out, g.root_auxin_half_saturation);
    write_val(out, g.leaf_auxin_max_boost);
    write_val(out, g.leaf_auxin_half_saturation);
    write_val(out, g.apical_auxin_max_boost);
    write_val(out, g.apical_auxin_half_saturation);
    write_val(out, g.root_apical_auxin_max_boost);
    write_val(out, g.root_apical_auxin_half_saturation);
    write_val(out, g.cytokinin_production_rate);
    write_val(out, g.cytokinin_diffusion_rate);
    write_val(out, g.cytokinin_decay_rate);
    write_val(out, g.cytokinin_threshold);
    write_val(out, g.cytokinin_growth_threshold);
    write_val(out, g.cytokinin_bias);
    write_val(out, g.hormone_base_transport);
    write_val(out, g.hormone_transport_scale);
    write_val(out, g.sugar_base_transport);
    write_val(out, g.sugar_transport_scale);
    write_val(out, g.stem_photosynthesis_rate);
    write_val(out, g.stem_green_radius_threshold);
    write_val(out, g.growth_rate);
    write_val(out, g.shoot_plastochron);
    write_val(out, g.branch_angle);
    write_val(out, g.cambium_responsiveness);
    write_val(out, g.internode_elongation_rate);
    write_val(out, g.max_internode_length);
    write_val(out, g.internode_maturation_ticks);
    write_val(out, g.root_growth_rate);
    write_val(out, g.root_plastochron);
    write_val(out, g.root_branch_angle);
    write_val(out, g.root_internode_elongation_rate);
    write_val(out, g.root_internode_maturation_ticks);
    write_val(out, g.root_gravitropism_strength);
    write_val(out, g.root_gravitropism_depth);
    write_val(out, g.root_cytokinin_production_rate);
    write_val(out, g.root_tip_auxin_production_rate);
    write_val(out, g.root_auxin_growth_threshold);
    write_val(out, g.root_ck_growth_floor);
    write_val(out, g.root_auxin_activation_threshold);
    write_val(out, g.root_cytokinin_inhibition_threshold);
    write_val(out, g.primary_root_lateral_delay_internodes);
    write_val(out, g.max_leaf_size);
    write_val(out, g.leaf_growth_rate);
    write_val(out, g.leaf_bud_size);
    write_val(out, g.leaf_petiole_length);
    write_val(out, g.leaf_opacity);
    write_val(out, g.initial_radius);
    write_val(out, g.root_initial_radius);
    write_val(out, g.tip_offset);
    write_val(out, g.growth_noise);
    write_val(out, g.sugar_diffusion_rate);
    write_val(out, g.seed_sugar);
    write_val(out, g.sugar_storage_density_wood);
    write_val(out, g.sugar_storage_density_leaf);
    write_val(out, g.sugar_cap_minimum);
    write_val(out, g.sugar_cap_meristem);
    write_val(out, g.water_absorption_rate);
    write_val(out, g.transpiration_rate);
    write_val(out, g.photosynthesis_water_ratio);
    write_val(out, g.water_storage_density_stem);
    write_val(out, g.water_storage_density_leaf);
    write_val(out, g.water_cap_meristem);
    write_val(out, g.water_diffusion_rate);
    write_val(out, g.water_bias);
    write_val(out, g.water_base_transport);
    write_val(out, g.water_transport_scale);
    write_val(out, g.leaf_phototropism_rate);
    write_val(out, g.meristem_gravitropism_rate);
    write_val(out, g.meristem_phototropism_rate);
    write_val(out, g.ga_production_rate);
    write_val(out, g.ga_leaf_age_max);
    write_val(out, g.ga_elongation_sensitivity);
    write_val(out, g.ga_length_sensitivity);
    write_val(out, g.ga_diffusion_rate);
    write_val(out, g.ga_decay_rate);
    write_val(out, g.leaf_abscission_ticks);
    write_val(out, g.min_leaf_age_before_abscission);
    write_val(out, g.ethylene_starvation_rate);
    write_val(out, g.ethylene_starvation_tick_threshold);
    write_val(out, g.ethylene_shade_rate);
    write_val(out, g.ethylene_shade_threshold);
    write_val(out, g.ethylene_age_rate);
    write_val(out, g.ethylene_age_onset);
    write_val(out, g.ethylene_crowding_rate);
    write_val(out, g.ethylene_crowding_radius);
    write_val(out, g.ethylene_diffusion_radius);
    write_val(out, g.ethylene_abscission_threshold);
    write_val(out, g.ethylene_elongation_inhibition);
    write_val(out, g.senescence_duration);
    write_val(out, g.wood_density);
    write_val(out, g.wood_flexibility);
    write_val(out, g.stress_hormone_threshold);
    write_val(out, g.stress_hormone_production_rate);
    write_val(out, g.stress_hormone_diffusion_rate);
    write_val(out, g.stress_hormone_decay_rate);
    write_val(out, g.stress_thickening_boost);
    write_val(out, g.stress_elongation_inhibition);
    write_val(out, g.stress_gravitropism_boost);
    write_val(out, g.elastic_recovery_rate);
    write_val(out, g.smoothing_rate);
    write_val(out, g.canalization_weight);
    write_val(out, g.pin_capacity_per_area);
    write_val(out, g.pin_base_efficiency);
    write_val(out, g.meristem_sink_fraction);
    write_val(out, g.vascular_radius_threshold);
    write_val(out, g.base_radial_permeability_sugar);
    write_val(out, g.radial_floor_fraction_sugar);
    write_val(out, g.radial_half_radius_sugar);
    write_val(out, g.base_radial_permeability_water);
    write_val(out, g.radial_floor_fraction_water);
    write_val(out, g.radial_half_radius_water);
    write_val(out, g.phloem_fraction);
    write_val(out, g.xylem_fraction);
    write_val(out, g.leaf_reserve_fraction_sugar);
    write_val(out, g.meristem_sink_target_fraction);
    write_val(out, g.leaf_turgor_target_fraction);
    write_val(out, g.root_water_reserve_fraction);
}

Genome read_genome_binary(std::istream& in) {
    Genome g{};
    g.apical_auxin_baseline = read_val<float>(in);
    g.apical_growth_auxin_multiplier = read_val<float>(in);
    g.auxin_diffusion_rate = read_val<float>(in);
    g.auxin_decay_rate = read_val<float>(in);
    g.auxin_threshold = read_val<float>(in);
    g.auxin_shade_boost = read_val<float>(in);
    g.auxin_sugar_half_saturation = read_val<float>(in);
    g.auxin_water_half_saturation = read_val<float>(in);
    g.cytokinin_sugar_half_saturation = read_val<float>(in);
    g.cytokinin_water_half_saturation = read_val<float>(in);
    g.quiescence_threshold = read_val<float>(in);
    g.meristem_dormancy_death_ticks = read_val<uint32_t>(in);
    g.starvation_ticks_max_stem = read_val<uint32_t>(in);
    g.starvation_ticks_max_root = read_val<uint32_t>(in);
    g.auxin_bias = read_val<float>(in);
    g.leaf_auxin_baseline = read_val<float>(in);
    g.leaf_growth_auxin_multiplier = read_val<float>(in);
    g.stem_auxin_max_boost = read_val<float>(in);
    g.stem_auxin_half_saturation = read_val<float>(in);
    g.root_auxin_max_boost = read_val<float>(in);
    g.root_auxin_half_saturation = read_val<float>(in);
    g.leaf_auxin_max_boost = read_val<float>(in);
    g.leaf_auxin_half_saturation = read_val<float>(in);
    g.apical_auxin_max_boost = read_val<float>(in);
    g.apical_auxin_half_saturation = read_val<float>(in);
    g.root_apical_auxin_max_boost = read_val<float>(in);
    g.root_apical_auxin_half_saturation = read_val<float>(in);
    g.cytokinin_production_rate = read_val<float>(in);
    g.cytokinin_diffusion_rate = read_val<float>(in);
    g.cytokinin_decay_rate = read_val<float>(in);
    g.cytokinin_threshold = read_val<float>(in);
    g.cytokinin_growth_threshold = read_val<float>(in);
    g.cytokinin_bias = read_val<float>(in);
    g.hormone_base_transport = read_val<float>(in);
    g.hormone_transport_scale = read_val<float>(in);
    g.sugar_base_transport = read_val<float>(in);
    g.sugar_transport_scale = read_val<float>(in);
    g.stem_photosynthesis_rate = read_val<float>(in);
    g.stem_green_radius_threshold = read_val<float>(in);
    g.growth_rate = read_val<float>(in);
    g.shoot_plastochron = read_val<uint32_t>(in);
    g.branch_angle = read_val<float>(in);
    g.cambium_responsiveness = read_val<float>(in);
    g.internode_elongation_rate = read_val<float>(in);
    g.max_internode_length = read_val<float>(in);
    g.internode_maturation_ticks = read_val<uint32_t>(in);
    g.root_growth_rate = read_val<float>(in);
    g.root_plastochron = read_val<uint32_t>(in);
    g.root_branch_angle = read_val<float>(in);
    g.root_internode_elongation_rate = read_val<float>(in);
    g.root_internode_maturation_ticks = read_val<uint32_t>(in);
    g.root_gravitropism_strength = read_val<float>(in);
    g.root_gravitropism_depth = read_val<float>(in);
    g.root_cytokinin_production_rate = read_val<float>(in);
    g.root_tip_auxin_production_rate = read_val<float>(in);
    g.root_auxin_growth_threshold = read_val<float>(in);
    g.root_ck_growth_floor = read_val<float>(in);
    g.root_auxin_activation_threshold = read_val<float>(in);
    g.root_cytokinin_inhibition_threshold = read_val<float>(in);
    g.primary_root_lateral_delay_internodes = read_val<uint32_t>(in);
    g.max_leaf_size = read_val<float>(in);
    g.leaf_growth_rate = read_val<float>(in);
    g.leaf_bud_size = read_val<float>(in);
    g.leaf_petiole_length = read_val<float>(in);
    g.leaf_opacity = read_val<float>(in);
    g.initial_radius = read_val<float>(in);
    g.root_initial_radius = read_val<float>(in);
    g.tip_offset = read_val<float>(in);
    g.growth_noise = read_val<float>(in);
    g.sugar_diffusion_rate = read_val<float>(in);
    g.seed_sugar = read_val<float>(in);
    g.sugar_storage_density_wood = read_val<float>(in);
    g.sugar_storage_density_leaf = read_val<float>(in);
    g.sugar_cap_minimum = read_val<float>(in);
    g.sugar_cap_meristem = read_val<float>(in);
    g.water_absorption_rate = read_val<float>(in);
    g.transpiration_rate = read_val<float>(in);
    g.photosynthesis_water_ratio = read_val<float>(in);
    g.water_storage_density_stem = read_val<float>(in);
    g.water_storage_density_leaf = read_val<float>(in);
    g.water_cap_meristem = read_val<float>(in);
    g.water_diffusion_rate = read_val<float>(in);
    g.water_bias = read_val<float>(in);
    g.water_base_transport = read_val<float>(in);
    g.water_transport_scale = read_val<float>(in);
    g.leaf_phototropism_rate = read_val<float>(in);
    g.meristem_gravitropism_rate = read_val<float>(in);
    g.meristem_phototropism_rate = read_val<float>(in);
    g.ga_production_rate = read_val<float>(in);
    g.ga_leaf_age_max = read_val<uint32_t>(in);
    g.ga_elongation_sensitivity = read_val<float>(in);
    g.ga_length_sensitivity = read_val<float>(in);
    g.ga_diffusion_rate = read_val<float>(in);
    g.ga_decay_rate = read_val<float>(in);
    g.leaf_abscission_ticks = read_val<uint32_t>(in);
    g.min_leaf_age_before_abscission = read_val<uint32_t>(in);
    g.ethylene_starvation_rate = read_val<float>(in);
    g.ethylene_starvation_tick_threshold = read_val<uint32_t>(in);
    g.ethylene_shade_rate = read_val<float>(in);
    g.ethylene_shade_threshold = read_val<float>(in);
    g.ethylene_age_rate = read_val<float>(in);
    g.ethylene_age_onset = read_val<uint32_t>(in);
    g.ethylene_crowding_rate = read_val<float>(in);
    g.ethylene_crowding_radius = read_val<float>(in);
    g.ethylene_diffusion_radius = read_val<float>(in);
    g.ethylene_abscission_threshold = read_val<float>(in);
    g.ethylene_elongation_inhibition = read_val<float>(in);
    g.senescence_duration = read_val<uint32_t>(in);
    g.wood_density = read_val<float>(in);
    g.wood_flexibility = read_val<float>(in);
    g.stress_hormone_threshold = read_val<float>(in);
    g.stress_hormone_production_rate = read_val<float>(in);
    g.stress_hormone_diffusion_rate = read_val<float>(in);
    g.stress_hormone_decay_rate = read_val<float>(in);
    g.stress_thickening_boost = read_val<float>(in);
    g.stress_elongation_inhibition = read_val<float>(in);
    g.stress_gravitropism_boost = read_val<float>(in);
    g.elastic_recovery_rate = read_val<float>(in);
    g.smoothing_rate = read_val<float>(in);
    g.canalization_weight = read_val<float>(in);
    g.pin_capacity_per_area = read_val<float>(in);
    g.pin_base_efficiency = read_val<float>(in);
    g.meristem_sink_fraction = read_val<float>(in);
    g.vascular_radius_threshold = read_val<float>(in);
    g.base_radial_permeability_sugar = read_val<float>(in);
    g.radial_floor_fraction_sugar = read_val<float>(in);
    g.radial_half_radius_sugar = read_val<float>(in);
    g.base_radial_permeability_water = read_val<float>(in);
    g.radial_floor_fraction_water = read_val<float>(in);
    g.radial_half_radius_water = read_val<float>(in);
    g.phloem_fraction = read_val<float>(in);
    g.xylem_fraction = read_val<float>(in);
    g.leaf_reserve_fraction_sugar = read_val<float>(in);
    g.meristem_sink_target_fraction = read_val<float>(in);
    g.leaf_turgor_target_fraction = read_val<float>(in);
    g.root_water_reserve_fraction = read_val<float>(in);
    return g;
}

namespace {

void write_chem_map(std::ostream& out, const std::unordered_map<ChemicalID, float>& m) {
    uint16_t count = static_cast<uint16_t>(m.size());
    write_val(out, count);
    for (const auto& kv : m) {
        write_val(out, static_cast<uint8_t>(kv.first));
        write_val(out, kv.second);
    }
}

std::unordered_map<ChemicalID, float> read_chem_map(std::istream& in) {
    uint16_t count = read_val<uint16_t>(in);
    std::unordered_map<ChemicalID, float> m;
    m.reserve(count);
    for (uint16_t i = 0; i < count; i++) {
        auto id = static_cast<ChemicalID>(read_val<uint8_t>(in));
        float v = read_val<float>(in);
        m[id] = v;
    }
    return m;
}

void write_bias_map(std::ostream& out, const std::unordered_map<Node*, float>& m) {
    uint16_t count = static_cast<uint16_t>(m.size());
    write_val(out, count);
    for (const auto& kv : m) {
        // Safe: by the time we serialize, every node has a valid id.
        uint32_t child_id = kv.first ? kv.first->id : UINT32_MAX;
        write_val(out, child_id);
        write_val(out, kv.second);
    }
}

std::unordered_map<uint32_t, float> read_bias_map(std::istream& in) {
    uint16_t count = read_val<uint16_t>(in);
    std::unordered_map<uint32_t, float> m;
    m.reserve(count);
    for (uint16_t i = 0; i < count; i++) {
        uint32_t id = read_val<uint32_t>(in);
        float v    = read_val<float>(in);
        m[id] = v;
    }
    return m;
}

} // namespace

void write_node_common(std::ostream& out, const Node& n, uint32_t parent_id) {
    write_val(out, n.id);
    write_val(out, parent_id);
    write_val(out, static_cast<uint8_t>(n.type));
    write_val(out, n.age);
    write_val(out, n.starvation_ticks);
    write_val(out, n.dormant_ticks);
    write_val(out, n.position);
    write_val(out, n.offset);
    write_val(out, n.rest_offset);
    write_val(out, n.radius);
    write_val(out, static_cast<uint8_t>(n.ever_active ? 1 : 0));
    write_chem_map(out, n.local().chemicals);
    write_bias_map(out, n.auxin_flow_bias);
}

NodeCommonRecord read_node_common(std::istream& in) {
    NodeCommonRecord r;
    r.id               = read_val<uint32_t>(in);
    r.parent_id        = read_val<uint32_t>(in);
    r.type             = static_cast<NodeType>(read_val<uint8_t>(in));
    r.age              = read_val<uint32_t>(in);
    r.starvation_ticks = read_val<uint32_t>(in);
    r.dormant_ticks    = read_val<uint32_t>(in);
    r.position         = read_val<glm::vec3>(in);
    r.offset           = read_val<glm::vec3>(in);
    r.rest_offset      = read_val<glm::vec3>(in);
    r.radius           = read_val<float>(in);
    r.ever_active      = read_val<uint8_t>(in) != 0;
    r.local_chemicals  = read_chem_map(in);
    r.auxin_flow_bias  = read_bias_map(in);
    return r;
}

// Stubs — filled in by later tasks.
SaveResult save_plant_snapshot(const Plant&, uint64_t, const std::string&) {
    return SaveResult{false, "", "save_plant_snapshot not implemented yet"};
}

LoadedPlant load_plant_snapshot(const std::string&, const std::optional<Genome>&) {
    throw std::runtime_error("load_plant_snapshot not implemented yet");
}

} // namespace botany
