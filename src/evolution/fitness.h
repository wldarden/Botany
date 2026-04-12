#pragma once

#include <cstdint>
#include "engine/genome.h"
#include "engine/world_params.h"

namespace botany {

struct PlantStats {
    uint32_t survival_ticks = 0;
    uint32_t node_count = 0;
    uint32_t leaf_count = 0;
    float total_sugar_produced = 0.0f;
    float height = 0.0f;
    float crown_ratio = 0.0f;
    uint32_t branch_depth = 0;
    float leaf_height_spread = 0.0f;
};

struct FitnessWeights {
    float survival     = 1.0f;
    float biomass      = 1.0f;
    float sugar        = 1.0f;
    float leaves       = 1.0f;
    float height       = 1.0f;
    float crown_ratio  = 1.0f;
    float branch_depth = 1.0f;
    float leaf_spread  = 1.0f;
};

PlantStats evaluate_plant(const Genome& genome, const WorldParams& world, uint32_t max_ticks);
float compute_fitness(const PlantStats& stats, const PlantStats& gen_max, const FitnessWeights& weights);

} // namespace botany
