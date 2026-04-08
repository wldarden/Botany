// src/engine/engine.h
#pragma once

#include <cstdint>
#include <memory>
#include <vector>
#include <glm/vec3.hpp>
#include "engine/genome.h"
#include "engine/plant.h"

namespace botany {

using PlantID = uint32_t;

class Engine {
public:
    Engine() = default;

    PlantID create_plant(const Genome& genome, glm::vec3 position);
    void tick();

    const Plant& get_plant(PlantID id) const;
    Plant& get_plant_mut(PlantID id);

    uint32_t get_tick() const { return tick_; }
    uint32_t plant_count() const { return static_cast<uint32_t>(plants_.size()); }

private:
    uint32_t tick_ = 0;
    std::vector<std::unique_ptr<Plant>> plants_;
};

} // namespace botany
