#include "serialization/plant_snapshot.h"
#include "engine/plant.h"
#include "engine/node/node.h"
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"
#include "engine/compartments.h"
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <climits>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

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

void write_leaf_trailer(std::ostream& out, const LeafNode& l) {
    write_val(out, l.leaf_size);
    write_val(out, l.light_exposure);
    write_val(out, l.senescence_ticks);
    write_val(out, l.deficit_ticks);
    write_val(out, l.facing);
}

LeafTrailer read_leaf_trailer(std::istream& in) {
    LeafTrailer t;
    t.leaf_size        = read_val<float>(in);
    t.light_exposure   = read_val<float>(in);
    t.senescence_ticks = read_val<uint32_t>(in);
    t.deficit_ticks    = read_val<uint32_t>(in);
    t.facing           = read_val<glm::vec3>(in);
    return t;
}

void write_apical_trailer(std::ostream& out, const ApicalNode& a) {
    write_val(out, static_cast<uint8_t>(a.active ? 1 : 0));
    write_val(out, static_cast<uint8_t>(a.is_primary ? 1 : 0));
    write_val(out, a.growth_dir);
    write_val(out, a.ticks_since_last_node);
}

ApicalTrailer read_apical_trailer(std::istream& in) {
    ApicalTrailer t;
    t.active                = read_val<uint8_t>(in) != 0;
    t.is_primary            = read_val<uint8_t>(in) != 0;
    t.growth_dir            = read_val<glm::vec3>(in);
    t.ticks_since_last_node = read_val<uint32_t>(in);
    return t;
}

void write_root_apical_trailer(std::ostream& out, const RootApicalNode& r) {
    write_val(out, static_cast<uint8_t>(r.active ? 1 : 0));
    write_val(out, static_cast<uint8_t>(r.is_primary ? 1 : 0));
    write_val(out, r.growth_dir);
    write_val(out, r.ticks_since_last_node);
    write_val(out, r.internodes_spawned);
}

RootApicalTrailer read_root_apical_trailer(std::istream& in) {
    RootApicalTrailer t;
    t.active                = read_val<uint8_t>(in) != 0;
    t.is_primary            = read_val<uint8_t>(in) != 0;
    t.growth_dir            = read_val<glm::vec3>(in);
    t.ticks_since_last_node = read_val<uint32_t>(in);
    t.internodes_spawned    = read_val<uint32_t>(in);
    return t;
}

void write_conduit_pools(std::ostream& out, const Node& n) {
    const TransportPool* p = n.phloem();
    const TransportPool* x = n.xylem();
    static const std::unordered_map<ChemicalID, float> kEmpty;
    write_chem_map(out, p ? p->chemicals : kEmpty);
    write_chem_map(out, x ? x->chemicals : kEmpty);
}

ConduitPools read_conduit_pools(std::istream& in) {
    ConduitPools p;
    p.phloem = read_chem_map(in);
    p.xylem  = read_chem_map(in);
    return p;
}

namespace {

std::string make_timestamp_now() {
    auto t = std::chrono::system_clock::now();
    std::time_t tt = std::chrono::system_clock::to_time_t(t);
    std::tm lt{};
#if defined(_WIN32)
    localtime_s(&lt, &tt);
#else
    localtime_r(&tt, &lt);
#endif
    std::ostringstream os;
    os << std::put_time(&lt, "%Y%m%d_%H%M%S");
    return os.str();
}

// DFS from seed so parent always precedes children in the file.
void dfs_collect(const Node* n, std::vector<const Node*>& out) {
    if (!n) return;
    out.push_back(n);
    for (const Node* c : n->children) dfs_collect(c, out);
}

void write_node_full(std::ostream& out, const Node& n) {
    uint32_t parent_id = n.parent ? n.parent->id : UINT32_MAX;
    write_node_common(out, n, parent_id);
    if (n.type == NodeType::STEM || n.type == NodeType::ROOT) {
        write_conduit_pools(out, n);
    }
    switch (n.type) {
        case NodeType::LEAF:        write_leaf_trailer(out, *n.as_leaf()); break;
        case NodeType::APICAL:      write_apical_trailer(out, *n.as_apical()); break;
        case NodeType::ROOT_APICAL: write_root_apical_trailer(out, *n.as_root_apical()); break;
        case NodeType::STEM:
        case NodeType::ROOT:        /* no trailer beyond conduit pools */ break;
    }
}

} // namespace

SaveResult save_plant_snapshot(const Plant& p, uint64_t engine_tick, const std::string& dir) {
    SaveResult r;
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) { r.error = "cannot create dir: " + ec.message(); return r; }

    std::ostringstream name;
    name << "plant_" << make_timestamp_now() << "_tick" << engine_tick << ".tree";
    auto full = std::filesystem::path(dir) / name.str();

    std::ofstream out(full, std::ios::binary);
    if (!out) { r.error = "cannot open output: " + full.string(); return r; }

    // Header
    plant_snapshot_write_magic(out);
    write_val(out, PLANT_SNAPSHOT_VERSION);
    write_val(out, engine_tick);
    write_genome_binary(out, p.genome());

    // Collect nodes in DFS order (parents precede children).
    std::vector<const Node*> order;
    dfs_collect(p.seed(), order);

    // next_node_id: one past the max id we saw, consistent with create_node contract.
    uint32_t max_id = 0;
    for (const Node* n : order) if (n->id >= max_id) max_id = n->id;
    uint32_t next_id = max_id + 1;

    uint32_t count = static_cast<uint32_t>(order.size());
    write_val(out, count);
    write_val(out, next_id);

    for (const Node* n : order) write_node_full(out, *n);

    out.flush();
    if (!out) { r.error = "write failed"; return r; }

    r.ok = true;
    r.path = full.string();
    return r;
}

LoadedPlant load_plant_snapshot(const std::string& path,
                                const std::optional<Genome>& genome_override) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("plant_snapshot: cannot open " + path);

    if (!plant_snapshot_check_magic(in))
        throw std::runtime_error("plant_snapshot: bad magic in " + path);

    uint32_t version = read_val<uint32_t>(in);
    if (version != PLANT_SNAPSHOT_VERSION)
        throw std::runtime_error("plant_snapshot: unsupported version " + std::to_string(version));

    uint64_t engine_tick = read_val<uint64_t>(in);
    Genome file_genome   = read_genome_binary(in);
    uint32_t node_count  = read_val<uint32_t>(in);
    uint32_t next_id     = read_val<uint32_t>(in);

    Genome active = genome_override ? *genome_override : file_genome;

    auto plant = Plant::from_empty(active);

    struct LoadedRecord {
        NodeCommonRecord common;
        std::optional<ConduitPools> pools;
        std::optional<LeafTrailer>  leaf;
        std::optional<ApicalTrailer> apical;
        std::optional<RootApicalTrailer> root_apical;
    };

    std::vector<LoadedRecord> records;
    records.reserve(node_count);

    // Pass 1a: read records into memory.
    for (uint32_t i = 0; i < node_count; i++) {
        LoadedRecord rec;
        rec.common = read_node_common(in);
        if (rec.common.type == NodeType::STEM || rec.common.type == NodeType::ROOT)
            rec.pools = read_conduit_pools(in);
        switch (rec.common.type) {
            case NodeType::LEAF:        rec.leaf = read_leaf_trailer(in); break;
            case NodeType::APICAL:      rec.apical = read_apical_trailer(in); break;
            case NodeType::ROOT_APICAL: rec.root_apical = read_root_apical_trailer(in); break;
            case NodeType::STEM:
            case NodeType::ROOT:        break;
        }
        records.push_back(std::move(rec));
    }

    // Pass 1b: allocate subclass + fill common state + trailers + pools.
    std::unordered_map<uint32_t, Node*> by_id;
    for (const LoadedRecord& rec : records) {
        std::unique_ptr<Node> n;
        switch (rec.common.type) {
            case NodeType::STEM:        n = std::make_unique<StemNode>(rec.common.id, rec.common.position, rec.common.radius); break;
            case NodeType::ROOT:        n = std::make_unique<RootNode>(rec.common.id, rec.common.position, rec.common.radius); break;
            case NodeType::LEAF:        n = std::make_unique<LeafNode>(rec.common.id, rec.common.position, rec.common.radius); break;
            case NodeType::APICAL:      n = std::make_unique<ApicalNode>(rec.common.id, rec.common.position, rec.common.radius); break;
            case NodeType::ROOT_APICAL: n = std::make_unique<RootApicalNode>(rec.common.id, rec.common.position, rec.common.radius); break;
        }
        n->age              = rec.common.age;
        n->starvation_ticks = rec.common.starvation_ticks;
        n->dormant_ticks    = rec.common.dormant_ticks;
        n->offset           = rec.common.offset;
        n->rest_offset      = rec.common.rest_offset;
        n->position         = rec.common.position;
        n->ever_active      = rec.common.ever_active;
        n->local().chemicals = rec.common.local_chemicals;

        if (rec.pools) {
            if (auto* p = n->phloem()) p->chemicals = rec.pools->phloem;
            if (auto* x = n->xylem())  x->chemicals = rec.pools->xylem;
        }
        if (rec.leaf) {
            auto* l = n->as_leaf();
            l->leaf_size        = rec.leaf->leaf_size;
            l->light_exposure   = rec.leaf->light_exposure;
            l->senescence_ticks = rec.leaf->senescence_ticks;
            l->deficit_ticks    = rec.leaf->deficit_ticks;
            l->facing           = rec.leaf->facing;
        }
        if (rec.apical) {
            auto* a = n->as_apical();
            a->active                = rec.apical->active;
            a->is_primary            = rec.apical->is_primary;
            a->growth_dir            = rec.apical->growth_dir;
            a->ticks_since_last_node = rec.apical->ticks_since_last_node;
        }
        if (rec.root_apical) {
            auto* r = n->as_root_apical();
            r->active                = rec.root_apical->active;
            r->is_primary            = rec.root_apical->is_primary;
            r->growth_dir            = rec.root_apical->growth_dir;
            r->ticks_since_last_node = rec.root_apical->ticks_since_last_node;
            r->internodes_spawned    = rec.root_apical->internodes_spawned;
        }

        by_id[rec.common.id] = n.get();
        plant->install_node(std::move(n));
    }

    // Pass 2: wire parent/children pointers.
    for (const LoadedRecord& rec : records) {
        Node* self = by_id.at(rec.common.id);
        if (rec.common.parent_id == UINT32_MAX) continue;
        auto pit = by_id.find(rec.common.parent_id);
        if (pit == by_id.end())
            throw std::runtime_error("plant_snapshot: node " + std::to_string(rec.common.id)
                                   + " references unknown parent " + std::to_string(rec.common.parent_id));
        self->parent = pit->second;
        pit->second->children.push_back(self);
    }

    // Pass 3: re-key auxin_flow_bias from child_id → child ptr.
    for (const LoadedRecord& rec : records) {
        Node* self = by_id.at(rec.common.id);
        for (const auto& kv : rec.common.auxin_flow_bias) {
            auto cit = by_id.find(kv.first);
            if (cit != by_id.end()) self->auxin_flow_bias[cit->second] = kv.second;
        }
    }

    plant->set_next_id(next_id);

    LoadedPlant out;
    out.plant = std::move(plant);
    out.genome = active;
    out.engine_tick = engine_tick;
    return out;
}

} // namespace botany
