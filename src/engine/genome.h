// src/engine/genome.h
#pragma once

#include <cstdint>

namespace botany {

struct Genome {
    // Hormone production & sensitivity
    float auxin_production_rate;
    float auxin_transport_rate;
    float auxin_decay_rate;
    float auxin_threshold;

    float cytokinin_production_rate;
    float cytokinin_transport_rate;
    float cytokinin_decay_rate;
    float cytokinin_threshold;

    // Shoot growth
    float growth_rate;
    float max_internode_length;
    uint32_t internode_spacing;
    float branch_angle;
    float thickening_rate;

    // Root growth
    float root_growth_rate;
    float root_max_internode_length;
    uint32_t root_internode_spacing;
    float root_branch_angle;

    // Geometry
    float leaf_size;
    float initial_radius;
};

inline Genome default_genome() {
    return Genome{
        .auxin_production_rate = 1.0f,
        .auxin_transport_rate = 0.3f,
        .auxin_decay_rate = 0.05f,
        .auxin_threshold = 0.5f,

        .cytokinin_production_rate = 1.0f,
        .cytokinin_transport_rate = 0.3f,
        .cytokinin_decay_rate = 0.05f,
        .cytokinin_threshold = 0.5f,

        .growth_rate = 0.1f,
        .max_internode_length = 1.0f,
        .internode_spacing = 5,
        .branch_angle = 0.785f,  // ~45 degrees
        .thickening_rate = 0.001f,

        .root_growth_rate = 0.08f,
        .root_max_internode_length = 0.8f,
        .root_internode_spacing = 4,
        .root_branch_angle = 0.925f,  // ~53 degrees

        .leaf_size = 0.3f,
        .initial_radius = 0.05f,
    };
}

} // namespace botany
