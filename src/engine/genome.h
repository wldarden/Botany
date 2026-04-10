// src/engine/genome.h
#pragma once

#include <cstdint>

namespace botany {

struct Genome {
    // Hormone production & sensitivity
    float auxin_production_rate;
    float auxin_transport_rate;
    float auxin_spillback_rate;   // fraction of junction auxin that spills back into branches
    float auxin_decay_rate;
    float auxin_threshold;

    float cytokinin_production_rate;
    float cytokinin_transport_rate;
    float cytokinin_decay_rate;
    float cytokinin_threshold;

    // Shoot growth
    float growth_rate;
    float max_internode_length;
    float branch_angle;
    float thickening_rate;

    // Root growth
    float root_growth_rate;
    float root_max_internode_length;
    float root_branch_angle;
    float root_gravitropism_strength;  // how strongly roots turn downward near surface
    float root_gravitropism_depth;     // depth at which gravitropism correction begins

    // Geometry
    float leaf_size;
    float initial_radius;
    float root_initial_radius;
    float tip_offset;     // small forward offset when chaining new tip/axillary nodes
    float growth_noise;   // max angular perturbation per segment (radians)

    // Sugar / photosynthesis
    float sugar_production_rate;
    float sugar_transport_conductance;
    float sugar_maintenance_leaf;
    float sugar_maintenance_stem;
    float sugar_maintenance_root;
    float sugar_maintenance_meristem;

    float seed_sugar;                   // initial sugar stored in the seed node

    // Sugar costs — per activity type
    float sugar_cost_shoot_growth;      // cost per unit of shoot extension
    float sugar_cost_root_growth;       // cost per unit of root extension
    float sugar_cost_thickening;        // cost per unit of radius increase
    float sugar_cost_activation;        // cost to activate an axillary meristem

    // Sugar save thresholds — per node type
    float sugar_save_shoot;             // reserve for shoot apical meristems
    float sugar_save_root;              // reserve for root apical meristems
    float sugar_save_stem;              // reserve for stem thickening
};

inline Genome default_genome() {
    return Genome{
        .auxin_production_rate = 1.0f,
        .auxin_transport_rate = 0.3f,
        .auxin_spillback_rate = 0.1f,
        .auxin_decay_rate = 0.05f,
        .auxin_threshold = 0.15f,

        .cytokinin_production_rate = 1.0f,
        .cytokinin_transport_rate = 0.3f,
        .cytokinin_decay_rate = 0.05f,
        .cytokinin_threshold = 0.15f,

        .growth_rate = 0.1f,
        .max_internode_length = 1.0f,
        .branch_angle = 0.785f,  // ~45 degrees
        .thickening_rate = 0.001f,

        .root_growth_rate = 0.04f,
        .root_max_internode_length = 0.8f,
        .root_branch_angle = 0.35f,   // ~20 degrees — tight spread, pushing through soil
        .root_gravitropism_strength = 2.0f,
        .root_gravitropism_depth = 0.5f,

        .leaf_size = 0.3f,
        .initial_radius = 0.05f,
        .root_initial_radius = 0.025f,
        .tip_offset = 0.01f,
        .growth_noise = 0.03f,  // ~1.7 degrees

        .sugar_production_rate = 0.5f,
        .sugar_transport_conductance = 10.0f,
        .sugar_maintenance_leaf = 0.02f,
        .sugar_maintenance_stem = 0.01f,
        .sugar_maintenance_root = 0.01f,
        .sugar_maintenance_meristem = 0.005f,

        .seed_sugar = 50.0f,
        .sugar_cost_shoot_growth = 0.3f,
        .sugar_cost_root_growth = 0.2f,
        .sugar_cost_thickening = 0.1f,
        .sugar_cost_activation = 1.0f,

        .sugar_save_shoot = 0.05f,
        .sugar_save_root = 0.05f,
        .sugar_save_stem = 0.02f,
    };
}

} // namespace botany
