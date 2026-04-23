// src/engine/engine.cpp
#include "engine/engine.h"
#include "engine/light.h"
#include "engine/node/node.h"
#include "engine/genome.h"
#include <filesystem>
#include <fstream>

namespace botany {

namespace {

const char* node_type_str(NodeType t) {
    switch (t) {
        case NodeType::STEM:        return "STEM";
        case NodeType::ROOT:        return "ROOT";
        case NodeType::LEAF:        return "LEAF";
        case NodeType::APICAL:      return "SA";
        case NodeType::ROOT_APICAL: return "RA";
    }
    return "?";
}

// Per-tick canalization logger: one row per parent-child edge.
// Reads last_auxin_flux (preserved through tick by Plant::tick_tree) and
// auxin_flow_bias (persistent). File is truncated on tick 0, appended after.
void write_canalization_log(uint32_t tick, const Plant& plant, const Genome& g) {
    std::filesystem::create_directories("debug");
    auto mode = (tick == 0)
        ? (std::ios::out | std::ios::trunc)
        : (std::ios::out | std::ios::app);
    std::ofstream csv("debug/canalization_log.csv", mode);
    if (!csv.is_open()) return;
    if (tick == 0) {
        csv << "tick,parent_id,child_id,child_type,child_radius,"
               "auxin_flux,auxin_flow_bias,pin_capacity,pin_saturation\n";
    }
    plant.for_each_node([&](const Node& parent) {
        for (const Node* child : parent.children) {
            float flux = 0.0f;
            auto fit = parent.last_auxin_flux.find(const_cast<Node*>(child));
            if (fit != parent.last_auxin_flux.end()) flux = fit->second;
            float bias = 0.0f;
            auto bit = parent.auxin_flow_bias.find(const_cast<Node*>(child));
            if (bit != parent.auxin_flow_bias.end()) bias = bit->second;
            float r2 = child->radius * child->radius;
            float capacity = r2 * g.pin_capacity_per_area;
            float saturation = (capacity > 1e-8f)
                ? std::min(flux / capacity, 1.0f) : 0.0f;
            csv << tick << ','
                << parent.id << ',' << child->id << ','
                << node_type_str(child->type) << ','
                << child->radius << ','
                << flux << ','
                << bias << ','
                << capacity << ','
                << saturation << '\n';
        }
    });
}

} // namespace

PlantID Engine::create_plant(const Genome& genome, glm::vec3 position) {
    PlantID id = static_cast<PlantID>(plants_.size());
    plants_.push_back(std::make_unique<Plant>(genome, position));
    return id;
}

void Engine::tick() {
    bool perf = perf_log_.is_open();

    // World-level light computation — skip on non-update ticks for performance.
    // Light/shadow changes slowly; recomputing every tick is wasteful.
    // Skipped entirely when the GPU LightSystem is driving light_exposure instead.
    if (world_params_.cpu_light_enabled) {
        uint32_t interval = std::max(1u, world_params_.light_update_interval);
        if (tick_ % interval == 0) {
            ScopedTimer t(perf_log_.stats().light_ms);
            compute_light_exposure(plants_, world_params_, &shadow_map_);
        }
    }

    world_params_.current_tick = tick_;

    {
        ScopedTimer t(perf_log_.stats().plant_tick_ms);
        for (auto& plant : plants_) {
            plant->tick(world_params_, perf ? &perf_log_.stats() : nullptr);
        }
    }

    if (debug_log_.is_open()) {
        ScopedTimer t(perf_log_.stats().debug_log_ms);
        for (const auto& plant : plants_) {
            debug_log_.log_tick(tick_, *plant, world_params_);
        }
    }

    if (global_economy_log_.is_open()) {
        for (const auto& plant : plants_) {
            global_economy_log_.log_tick(tick_, *plant, world_params_);
        }
    }

    if (world_params_.canalization_debug_log) {
        for (const auto& plant : plants_) {
            write_canalization_log(tick_, *plant, plant->genome());
        }
    }

    if (perf) {
        // Count nodes and aggregate per-node sugar accounting for this tick.
        auto& s = perf_log_.stats();
        for (const auto& plant : plants_) {
            s.node_count += plant->node_count();
            plant->for_each_node([&](const Node& n) {
                s.sugar_spent_maintenance += n.tick_sugar_maintenance;
                if (n.tick_sugar_activity  < 0.0f) s.sugar_spent_growth    += -n.tick_sugar_activity;
                if (n.tick_sugar_transport < 0.0f) s.sugar_spent_transport += -n.tick_sugar_transport;
            });
        }
        perf_log_.flush(tick_);
    }

    // Autocompress: runs at tick-interval boundaries, between full ticks.
    if (compression_enabled_
        && tick_ > 0
        && (tick_ % compression_interval_) == 0
        && !plants_.empty()) {
        last_compression_ = compress_plant(*plants_[0], compression_params_);
    }

    tick_++;
}

void Engine::reset() {
    plants_.clear();
    tick_ = 0;
}

const Plant& Engine::get_plant(PlantID id) const {
    return *plants_.at(id);
}

Plant& Engine::get_plant_mut(PlantID id) {
    return *plants_.at(id);
}

PlantID Engine::adopt_plant(std::unique_ptr<Plant> plant) {
    PlantID id = static_cast<PlantID>(plants_.size());
    plants_.push_back(std::move(plant));
    return id;
}

void Engine::set_tick(uint32_t tick) {
    tick_ = tick;
}

void Engine::enable_autocompress(bool enabled) {
    compression_enabled_ = enabled;
}

void Engine::set_compression_interval(uint32_t ticks) {
    compression_interval_ = ticks == 0 ? 1 : ticks;
}

void Engine::set_compression_params(const CompressionParams& params) {
    compression_params_ = params;
}

CompressionResult Engine::trigger_compression() {
    if (plants_.empty()) {
        last_compression_ = CompressionResult{};
        return last_compression_;
    }
    last_compression_ = compress_plant(*plants_[0], compression_params_);
    return last_compression_;
}

const CompressionResult& Engine::last_compression() const {
    return last_compression_;
}

} // namespace botany
