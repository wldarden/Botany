#include "engine/debug_log.h"
#include "engine/plant.h"
#include "engine/sugar.h"
#include "engine/world_params.h"
#include "engine/node/node.h"
#include "engine/node/leaf_node.h"
#include "engine/chemical/chemical.h"
#include <iomanip>
#include <glm/geometric.hpp>

namespace botany {

void DebugLog::open(const std::string& path) {
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
              << "light_exposure,leaf_size,angle_eff,senescence,"
              << "starvation_ticks,pos_x,pos_y,pos_z,radius,"
              << "total_mass,stress\n";
        header_written_ = true;
    }

    const Genome& g = plant.genome();

    plant.for_each_node([&](const Node& node) {
        const char* type_str = "?";
        switch (node.type) {
            case NodeType::STEM: type_str = "STEM"; break;
            case NodeType::ROOT: type_str = "ROOT"; break;
            case NodeType::LEAF: type_str = "LEAF"; break;
            case NodeType::SHOOT_APICAL: type_str = "SA"; break;
            case NodeType::SHOOT_AXILLARY: type_str = "SX"; break;
            case NodeType::ROOT_APICAL: type_str = "RA"; break;
            case NodeType::ROOT_AXILLARY: type_str = "RX"; break;
        }

        float sugar = node.chemical(ChemicalID::Sugar);
        float cap = sugar_cap(node, g);
        float maint = node.maintenance_cost(world);
        float parent_sugar = node.parent ? node.parent->chemical(ChemicalID::Sugar) : -1.0f;

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

        file_ << tick << "," << node.id << "," << parent_id << ","
              << type_str << "," << node.age << "," << node.children.size() << ","
              << std::fixed << std::setprecision(6)
              << sugar << "," << cap << "," << maint << "," << production << ","
              << node.chemical(ChemicalID::Auxin) << ","
              << node.chemical(ChemicalID::Cytokinin) << ","
              << node.chemical(ChemicalID::Gibberellin) << ","
              << node.chemical(ChemicalID::Ethylene) << ","
              << node.chemical(ChemicalID::Stress) << ","
              << light_exp << "," << leaf_sz << "," << angle_eff << ","
              << senescence << ","
              << node.starvation_ticks << ","
              << node.position.x << "," << node.position.y << "," << node.position.z << ","
              << node.radius << ","
              << node.total_mass << "," << node.stress << "\n";
    });
}

} // namespace botany
