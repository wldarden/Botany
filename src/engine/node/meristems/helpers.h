// src/engine/meristems/helpers.h
// Shared helper functions used by meristem tick() implementations.
#pragma once

#include "engine/node/node.h"
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>
#include <cmath>
#include <random>

namespace botany {
namespace meristem_helpers {

inline std::mt19937& rng() {
    static std::mt19937 instance{std::random_device{}()};
    return instance;
}

// Compute growth direction for a node: normalized vector from parent to this node
inline glm::vec3 growth_direction(const Node& node) {
    if (node.parent) {
        float len = glm::length(node.offset);
        if (len > 0.0001f) {
            return node.offset / len;
        }
    }
    bool root_type = node.type == NodeType::ROOT
                  || node.type == NodeType::ROOT_APICAL;
    return root_type ? glm::vec3(0.0f, -1.0f, 0.0f) : glm::vec3(0.0f, 1.0f, 0.0f);
}

// Compute a branch direction in the node's local reference frame.
inline glm::vec3 branch_direction(const glm::vec3& main_dir, float angle, uint32_t seed) {
    glm::vec3 gravity(0.0f, -1.0f, 0.0f);

    glm::vec3 perp;
    float alignment = std::abs(glm::dot(main_dir, gravity));
    if (alignment < 0.99f) {
        perp = glm::normalize(glm::cross(main_dir, gravity));
    } else {
        perp = glm::normalize(glm::cross(main_dir, glm::vec3(1.0f, 0.0f, 0.0f)));
    }

    glm::vec3 perp2 = glm::cross(main_dir, perp);

    float rotate = static_cast<float>(seed) * 2.399f; // golden angle in radians
    float cos_r = std::cos(rotate);
    float sin_r = std::sin(rotate);
    glm::vec3 radial = perp * cos_r + perp2 * sin_r;

    float cos_a = std::cos(angle);
    float sin_a = std::sin(angle);
    return glm::normalize(main_dir * cos_a + radial * sin_a);
}

// Roll a target internode length based on sugar availability + randomness.
inline float roll_internode_length(float min_len, float max_len, float sugar_fraction) {
    float base = min_len + (max_len - min_len) * sugar_fraction;
    std::uniform_real_distribution<float> jitter(0.8f, 1.2f);
    float result = base * jitter(rng());
    return std::clamp(result, min_len, max_len);
}

// Perturb a direction by a random angle up to max_angle radians
inline glm::vec3 perturb(const glm::vec3& dir, float max_angle) {
    if (max_angle < 1e-6f) return dir;

    std::uniform_real_distribution<float> angle_dist(0.0f, max_angle);
    std::uniform_real_distribution<float> rot_dist(0.0f, 2.0f * 3.14159f);

    float tilt = angle_dist(rng());
    float rot = rot_dist(rng());

    glm::vec3 p;
    if (std::abs(dir.y) < 0.9f) {
        p = glm::normalize(glm::cross(dir, glm::vec3(0.0f, 1.0f, 0.0f)));
    } else {
        p = glm::normalize(glm::cross(dir, glm::vec3(1.0f, 0.0f, 0.0f)));
    }
    glm::vec3 p2 = glm::cross(dir, p);

    glm::vec3 radial = p * std::cos(rot) + p2 * std::sin(rot);
    return glm::normalize(dir * std::cos(tilt) + radial * std::sin(tilt));
}

// Compute growth fraction: sugar funds growth; a permissive-signal chemical
// gates the rate.  Called by two sites:
//   (a) ApicalNode::elongate — sugar + cytokinin delivered via xylem from roots
//   (b) RootApicalNode::elongate — sugar + cytokinin produced locally at the RA
// Both use Michaelis-Menten: signal / (signal + Km).  No hard cap — more
// signal always helps, with diminishing returns.  The parameter name
// "cytokinin" reflects the dominant caller pattern; semantically it's any
// chemical you want to graded-gate the rate on.
inline float growth_fraction(float sugar, float max_cost,
                             float cytokinin, float cyt_threshold) {
    if (max_cost < 1e-6f) return 1.0f;
    float sugar_gf = std::min(sugar / max_cost, 1.0f);
    if (sugar_gf < 1e-6f) return 0.0f;
    float cyt_gf = cytokinin / (cytokinin + std::max(cyt_threshold, 1e-6f));
    return sugar_gf * cyt_gf;
}

// Sugar-only growth fraction for root apicals.
// Real root tips maintain their own auxin via PIN recycling — elongation
// is limited by photosynthate supply, not exogenous hormone signal.
inline float sugar_growth_fraction(float sugar, float max_cost) {
    if (max_cost < 1e-6f) return 1.0f;
    return std::min(sugar / max_cost, 1.0f);
}

// Turgor pressure fraction: water availability gates cell expansion.
// Linear — turgor scales directly with water content relative to capacity.
inline float turgor_fraction(float water, float water_cap) {
    if (water_cap < 1e-6f) return 1.0f;
    return std::clamp(water / water_cap, 0.0f, 1.0f);
}

// Saturating auxin growth multiplier (Michaelis-Menten).
// Returns 1.0 at zero auxin, asymptotes to 1.0 + max_boost.
// Positive max_boost = promotion, negative = inhibition.
inline float auxin_growth_factor(float auxin, float max_boost, float half_sat) {
    if (std::abs(max_boost) < 1e-8f) return 1.0f;
    float saturation = auxin / (auxin + std::max(half_sat, 1e-6f));
    return 1.0f + max_boost * saturation;
}

// Multiplicative metabolic gating — captures short-timescale biochemistry +
// transcription response to sugar and water at a meristem.  Both stresses
// compound.  Floors encode hormone-specific buffering (conjugate pools for
// auxin, etc.); caller supplies floor per hormone type.
//
// At sugar=0 and water=0, returns floor_s * floor_w (minimum production).
// At saturation (sugar >> K_s, water >> K_w), returns 1.0 (full production).
//
// Used by SA auxin, RA auxin, and RA cytokinin production.
inline float metabolic_factor(float sugar, float K_sugar, float floor_sugar,
                              float water, float K_water, float floor_water) {
    float sf = floor_sugar + (1.0f - floor_sugar) * sugar / (sugar + K_sugar);
    float wf = floor_water + (1.0f - floor_water) * water / (water + K_water);
    return sf * wf;
}

} // namespace meristem_helpers
} // namespace botany
