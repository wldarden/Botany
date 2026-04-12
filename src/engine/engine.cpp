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
    // World-level light computation (shared shadow map across all plants)
    compute_light_exposure(plants_, world_params_, &shadow_map_);

    for (auto& plant : plants_) {
        plant->tick(world_params_);
    }
    if (debug_log_.is_open()) {
        for (const auto& plant : plants_) {
            debug_log_.log_tick(tick_, *plant, world_params_);
        }
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
