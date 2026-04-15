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

// Compute growth fraction: sugar funds growth, cytokinin gates the rate.
// Cytokinin is produced by photosynthesizing leaves — no producing leaves = no growth.
// Cytokinin uses Michaelis-Menten kinetics: cyt / (cyt + Km). No hard cap —
// more cytokinin always helps, but with diminishing returns. Meristems far
// from roots get less cytokinin and grow proportionally slower.
inline float growth_fraction(float sugar, float max_cost,
                             float cytokinin, float cyt_threshold) {
    if (max_cost < 1e-6f) return 1.0f;
    float sugar_gf = std::min(sugar / max_cost, 1.0f);
    if (sugar_gf < 1e-6f) return 0.0f;
    float cyt_gf = cytokinin / (cytokinin + std::max(cyt_threshold, 1e-6f));
    return sugar_gf * cyt_gf;
}

} // namespace meristem_helpers
} // namespace botany
