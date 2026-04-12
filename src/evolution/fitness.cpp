#include "evolution/fitness.h"
#include "engine/engine.h"
#include "engine/node/node.h"
#include "engine/node/leaf_node.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace botany {

// Count branch depth: how many branch-from-branch generations.
// First child continues same branch; other children start new branches.
static uint32_t compute_branch_depth(const Plant& plant) {
    uint32_t max_depth = 0;
    struct Entry { const Node* node; uint32_t depth; };
    std::vector<Entry> stack;
    stack.push_back({plant.seed(), 0});

    while (!stack.empty()) {
        auto [node, depth] = stack.back();
        stack.pop_back();
        max_depth = std::max(max_depth, depth);

        if (node->children.empty()) continue;

        bool first = true;
        for (const Node* child : node->children) {
            if (child->type == NodeType::LEAF) continue;
            if (child->is_meristem()) continue;
            if (first) {
                stack.push_back({child, depth});
                first = false;
            } else {
                stack.push_back({child, depth + 1});
            }
        }
    }
    return max_depth;
}

static bool has_active_meristems(const Plant& plant) {
    bool found = false;
    plant.for_each_node([&](const Node& n) {
        if (n.type == NodeType::SHOOT_APICAL || n.type == NodeType::ROOT_APICAL) {
            found = true;
        }
    });
    return found;
}

// Collect stats from a single plant after simulation.
static PlantStats collect_stats(const Plant& plant, uint32_t survival_ticks) {
    PlantStats stats;
    stats.survival_ticks = survival_ticks;
    stats.node_count = plant.node_count();
    stats.total_sugar_produced = plant.total_sugar_produced();

    float min_y = 1e9f, max_y = -1e9f;
    float min_x = 1e9f, max_x = -1e9f;
    float min_z = 1e9f, max_z = -1e9f;
    std::vector<float> leaf_ys;

    plant.for_each_node([&](const Node& n) {
        if (!std::isfinite(n.position.x) || !std::isfinite(n.position.y) || !std::isfinite(n.position.z))
            return;  // skip nodes with NaN/Inf positions (degenerate genomes)

        min_y = std::min(min_y, n.position.y);
        max_y = std::max(max_y, n.position.y);
        min_x = std::min(min_x, n.position.x);
        max_x = std::max(max_x, n.position.x);
        min_z = std::min(min_z, n.position.z);
        max_z = std::max(max_z, n.position.z);

        if (n.type == NodeType::LEAF) {
            auto* leaf = n.as_leaf();
            if (leaf && leaf->senescence_ticks == 0) {
                stats.leaf_count++;
                if (std::isfinite(n.position.y)) {
                    leaf_ys.push_back(n.position.y);
                }
            }
        }
    });

    stats.height = (std::isfinite(max_y) && max_y > 0.0f) ? max_y : 0.0f;

    float width_x = max_x - min_x;
    float width_z = max_z - min_z;
    float canopy_width = std::max(width_x, width_z);
    float cr = stats.height > 0.01f ? canopy_width / stats.height : 0.0f;
    stats.crown_ratio = std::isfinite(cr) ? cr : 0.0f;

    stats.branch_depth = compute_branch_depth(plant);

    if (leaf_ys.size() > 1) {
        float mean = 0.0f;
        for (float y : leaf_ys) mean += y;
        mean /= static_cast<float>(leaf_ys.size());

        float variance = 0.0f;
        for (float y : leaf_ys) {
            float d = y - mean;
            variance += d * d;
        }
        variance /= static_cast<float>(leaf_ys.size());
        float spread = std::sqrt(variance);
        stats.leaf_height_spread = std::isfinite(spread) ? spread : 0.0f;
    }

    return stats;
}

PlantStats evaluate_plant(const Genome& genome, const WorldParams& world, uint32_t max_ticks) {
    Engine engine;
    engine.world_params_mut() = world;
    PlantID pid = engine.create_plant(genome, glm::vec3(0.0f));

    uint32_t ticks = 0;
    for (; ticks < max_ticks; ticks++) {
        engine.tick();
        if (!has_active_meristems(engine.get_plant(pid))) break;
    }

    return collect_stats(engine.get_plant(pid), ticks);
}

std::vector<PlantStats> evaluate_group(const std::vector<Genome>& genomes,
                                       const WorldParams& world, uint32_t max_ticks,
                                       float spacing) {
    Engine engine;
    engine.world_params_mut() = world;

    // Place plants spaced along X axis, centered around origin
    uint32_t n = static_cast<uint32_t>(genomes.size());
    float total_width = (n > 1) ? spacing * static_cast<float>(n - 1) : 0.0f;
    float start_x = -total_width * 0.5f;

    std::vector<PlantID> pids;
    pids.reserve(n);
    for (uint32_t i = 0; i < n; i++) {
        float x = start_x + spacing * static_cast<float>(i);
        pids.push_back(engine.create_plant(genomes[i], glm::vec3(x, 0.0f, 0.0f)));
    }

    // Run sim — stop when ALL plants have lost their meristems
    uint32_t ticks = 0;
    for (; ticks < max_ticks; ticks++) {
        engine.tick();

        bool any_active = false;
        for (uint32_t i = 0; i < n; i++) {
            if (has_active_meristems(engine.get_plant(pids[i]))) {
                any_active = true;
                break;
            }
        }
        if (!any_active) break;
    }

    // Collect stats per plant
    std::vector<PlantStats> results;
    results.reserve(n);
    for (uint32_t i = 0; i < n; i++) {
        results.push_back(collect_stats(engine.get_plant(pids[i]), ticks));
    }
    return results;
}

float compute_fitness(const PlantStats& stats, const PlantStats& gen_max, const FitnessWeights& w) {
    auto norm = [](float val, float max_val) -> float {
        return max_val > 1e-9f ? val / max_val : 0.0f;
    };

    float score = 0.0f;
    score += w.survival     * norm(static_cast<float>(stats.survival_ticks), static_cast<float>(gen_max.survival_ticks));
    score += w.biomass      * norm(static_cast<float>(stats.node_count), static_cast<float>(gen_max.node_count));
    score += w.sugar        * norm(stats.total_sugar_produced, gen_max.total_sugar_produced);
    score += w.leaves       * norm(static_cast<float>(stats.leaf_count), static_cast<float>(gen_max.leaf_count));
    score += w.height       * norm(stats.height, gen_max.height);
    score += w.crown_ratio  * norm(stats.crown_ratio, gen_max.crown_ratio);
    score += w.branch_depth * norm(static_cast<float>(stats.branch_depth), static_cast<float>(gen_max.branch_depth));
    score += w.leaf_spread  * norm(stats.leaf_height_spread, gen_max.leaf_height_spread);
    return score;
}

} // namespace botany
