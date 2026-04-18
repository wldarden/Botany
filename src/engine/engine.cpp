// src/engine/engine.cpp
#include "engine/engine.h"
#include "engine/light.h"

namespace botany {

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

} // namespace botany
