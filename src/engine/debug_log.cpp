#include "engine/debug_log.h"
#include "engine/plant.h"
#include "engine/sugar.h"
#include "engine/world_params.h"
#include "engine/node/node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"
#include "engine/chemical/chemical.h"
#include "engine/vascular_sub_stepped.h"
#include <iomanip>
#include <filesystem>
#include <glm/geometric.hpp>

namespace botany {

void DebugLog::open(const std::string& path) {
    std::filesystem::path p(path);
    if (p.has_parent_path()) {
        std::filesystem::create_directories(p.parent_path());
    }
    file_.open(path);
    header_written_ = false;
}

void DebugLog::close() {
    file_.close();
}

void DebugLog::log_tick(uint32_t tick, const Plant& plant, const WorldParams& world) {
    if (!file_.is_open()) return;

    if (!header_written_) {
        file_ << "tick,node_id,parent_id,type,age,children,"
              << "sugar,sugar_cap,maintenance,production,"
              << "auxin,cytokinin,gibberellin,ethylene,stress_hormone,"
              << "auxin_produced,cytokinin_produced,"
              << "light_exposure,leaf_size,angle_eff,senescence,"
              << "starvation_ticks,pos_x,pos_y,pos_z,radius,length,"
              << "total_mass,stress,"
              // Compartmented vascular diagnostics (2026-04-20):
              << "water,water_cap,"
              << "phloem_sugar,phloem_cap,"
              << "xylem_water,xylem_cytokinin,xylem_cap,"
              << "radial_perm_sugar,radial_perm_water,"
              << "active,is_primary\n";
        header_written_ = true;
    }

    const Genome& g = plant.genome();

    plant.for_each_node([&](const Node& node) {
        const char* type_str = "?";
        switch (node.type) {
            case NodeType::STEM: type_str = "STEM"; break;
            case NodeType::ROOT: type_str = "ROOT"; break;
            case NodeType::LEAF: type_str = "LEAF"; break;
            case NodeType::APICAL: type_str = "SA"; break;
            case NodeType::ROOT_APICAL: type_str = "RA"; break;
        }

        float sugar = node.local().chemical(ChemicalID::Sugar);
        float cap = sugar_cap(node, g);
        float maint = node.maintenance_cost(world);
        float parent_sugar = node.parent ? node.parent->local().chemical(ChemicalID::Sugar) : -1.0f;

        float light_exp = 0.0f;
        float leaf_sz = 0.0f;
        float angle_eff = 0.0f;
        float production = 0.0f;
        uint32_t senescence = 0;

        if (auto* leaf = node.as_leaf()) {
            light_exp = leaf->light_exposure;
            leaf_sz = leaf->leaf_size;
            senescence = leaf->senescence_ticks;

            // Compute what production would be (using facing direction, same as photosynthesize)
            float facing_len = glm::length(leaf->facing);
            if (facing_len > 1e-4f) {
                glm::vec3 leaf_normal = leaf->facing / facing_len;
                angle_eff = std::max(0.0f, leaf_normal.y);
            } else {
                angle_eff = 1.0f;
            }
            production = light_exp * angle_eff * world.light_level * leaf_sz * world.sugar_production_rate;
        }

        int32_t parent_id = node.parent ? static_cast<int32_t>(node.parent->id) : -1;

        // Compartmented vascular: phloem/xylem pool contents and capacities.
        // Columns are 0 when the node doesn't own the pool (leaves, meristems).
        const float water_loc       = node.local().chemical(ChemicalID::Water);
        const float water_cap_val   = water_cap(node, g);
        const auto* phl             = node.phloem();
        const auto* xyl             = node.xylem();
        const float phloem_sugar    = phl ? phl->chemical(ChemicalID::Sugar) : 0.0f;
        const float phloem_cap_val  = phl ? phloem_capacity(node, g)         : 0.0f;
        const float xylem_water     = xyl ? xyl->chemical(ChemicalID::Water)     : 0.0f;
        const float xylem_cyto      = xyl ? xyl->chemical(ChemicalID::Cytokinin) : 0.0f;
        const float xylem_cap_val   = xyl ? xylem_capacity(node, g)          : 0.0f;
        const float perm_sugar      = radial_permeability_sugar(node.radius, g);
        const float perm_water      = radial_permeability_water(node.radius, g);

        // Meristem activation state — 0 if not a meristem.
        int active_flag = 0;
        int primary_flag = 0;
        if (const auto* sa = node.as_apical()) {
            active_flag  = sa->active ? 1 : 0;
            primary_flag = sa->is_primary ? 1 : 0;
        } else if (const auto* ra = node.as_root_apical()) {
            active_flag  = ra->active ? 1 : 0;
            primary_flag = ra->is_primary ? 1 : 0;
        }

        file_ << tick << "," << node.id << "," << parent_id << ","
              << type_str << "," << node.age << "," << node.children.size() << ","
              << std::fixed << std::setprecision(6)
              << sugar << "," << cap << "," << maint << "," << production << ","
              << node.local().chemical(ChemicalID::Auxin) << ","
              << node.local().chemical(ChemicalID::Cytokinin) << ","
              << node.local().chemical(ChemicalID::Gibberellin) << ","
              << node.local().chemical(ChemicalID::Ethylene) << ","
              << node.local().chemical(ChemicalID::Stress) << ","
              << node.tick_auxin_produced << ","
              << node.tick_cytokinin_produced << ","
              << light_exp << "," << leaf_sz << "," << angle_eff << ","
              << senescence << ","
              << node.starvation_ticks << ","
              << node.position.x << "," << node.position.y << "," << node.position.z << ","
              << node.radius << "," << glm::length(node.offset) << ","
              << node.total_mass << "," << node.stress << ","
              << water_loc << "," << water_cap_val << ","
              << phloem_sugar << "," << phloem_cap_val << ","
              << xylem_water << "," << xylem_cyto << "," << xylem_cap_val << ","
              << perm_sugar << "," << perm_water << ","
              << active_flag << "," << primary_flag << "\n";
    });
}

} // namespace botany
