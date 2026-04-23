// src/engine/engine.h
#pragma once

#include <cstdint>
#include <memory>
#include <vector>
#include <glm/vec3.hpp>
#include "engine/genome.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include "engine/light.h"
#include "engine/debug_log.h"
#include "engine/global_economy_log.h"
#include "engine/perf_log.h"

namespace botany {

using PlantID = uint32_t;

class Engine {
public:
    Engine() = default;

    PlantID create_plant(const Genome& genome, glm::vec3 position);
    void tick();
    void reset();

    // Install a plant constructed outside the engine (e.g., from a snapshot).
    // Returns the assigned PlantID.  Takes ownership.
    PlantID adopt_plant(std::unique_ptr<Plant> plant);

    // Override the current tick counter (used when loading a snapshot so HUDs
    // and age-dependent signals see the correct elapsed sim time).
    void set_tick(uint32_t tick);

    const Plant& get_plant(PlantID id) const;
    Plant& get_plant_mut(PlantID id);

    uint32_t get_tick() const { return tick_; }
    uint32_t plant_count() const { return static_cast<uint32_t>(plants_.size()); }
    const std::vector<std::unique_ptr<Plant>>& all_plants() const { return plants_; }

    const WorldParams& world_params() const { return world_params_; }
    WorldParams& world_params_mut() { return world_params_; }

    const ShadowMapViz& shadow_map() const { return shadow_map_; }
    DebugLog& debug_log() { return debug_log_; }
    GlobalEconomyLog& global_economy_log() { return global_economy_log_; }
    PerfLog& perf_log() { return perf_log_; }

private:
    uint32_t tick_ = 0;
    std::vector<std::unique_ptr<Plant>> plants_;
    WorldParams world_params_;
    ShadowMapViz shadow_map_;
    DebugLog debug_log_;
    GlobalEconomyLog global_economy_log_;
    PerfLog perf_log_;
};

} // namespace botany
