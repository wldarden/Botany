#include "engine/ethylene.h"
#include "engine/plant.h"
#include "engine/node/node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/world_params.h"
#include <cmath>
#include <vector>
#include <glm/geometric.hpp>

namespace botany {

void compute_ethylene(Plant& plant, const WorldParams& /*world*/) {
    const Genome& g = plant.genome();

    // Collect node pointers and positions for spatial queries
    struct NodeInfo {
        Node* node;
        glm::vec3 pos;
    };
    std::vector<NodeInfo> nodes;
    plant.for_each_node_mut([&](Node& node) {
        nodes.push_back({&node, node.position});
    });

    // Phase 1: Reset (record prior tick's ethylene as consumed) and compute local production
    for (auto& info : nodes) {
        Node& node = *info.node;
        // Record last tick's value as consumed before zeroing (reset-each-tick signal model).
        float eth_prev = node.local().chemical(ChemicalID::Ethylene);
        if (eth_prev > 0.0f) {
            node.tick_chem_consumed[static_cast<size_t>(ChemicalID::Ethylene)] += eth_prev;
        }
        node.local().chemical(ChemicalID::Ethylene) = 0.0f;

        float eth_produced = 0.0f;

        // Trigger 1: sustained sugar starvation.  Gated on starvation_ticks
        // (consecutive ticks with sugar <= 0) rather than instantaneous sugar
        // level — transient dips during active growth shouldn't flood the
        // plant with a signal meant to flag genuine carbon crisis.
        if (node.starvation_ticks >= g.ethylene_starvation_tick_threshold) {
            eth_produced += g.ethylene_starvation_rate;
        }

        // Trigger 2: Low light (LEAF only)
        if (auto* leaf = node.as_leaf()) {
            if (leaf->light_exposure < g.ethylene_shade_threshold) {
                eth_produced += g.ethylene_shade_rate * (1.0f - leaf->light_exposure);
            }
        }

        // Trigger 3: Old age (LEAF only)
        if (node.type == NodeType::LEAF &&
            node.age > g.ethylene_age_onset) {
            float age_past = static_cast<float>(node.age - g.ethylene_age_onset);
            eth_produced += g.ethylene_age_rate * age_past
                           / static_cast<float>(g.ethylene_age_onset);
        }

        // Trigger 4: Crowding — count nearby nodes within crowding_radius
        float cr2 = g.ethylene_crowding_radius * g.ethylene_crowding_radius;
        int nearby_count = 0;
        for (const auto& other : nodes) {
            if (other.node == &node) continue;
            glm::vec3 diff = other.pos - info.pos;
            if (glm::dot(diff, diff) < cr2) {
                nearby_count++;
            }
        }
        eth_produced += g.ethylene_crowding_rate * static_cast<float>(nearby_count);

        node.local().chemical(ChemicalID::Ethylene) += eth_produced;
        if (eth_produced > 0.0f) {
            node.tick_chem_produced[static_cast<size_t>(ChemicalID::Ethylene)] += eth_produced;
        }
    }

    // Phase 2: Spatial gas diffusion — mass-conserving.  Each source's
    // emitted ethylene is distributed across itself and in-range neighbors
    // by falloff weights normalized to sum to 1.  Total gas is preserved
    // regardless of neighborhood density.  Previously the formula was
    // non-conservative (each source cloned its full emission onto every
    // receiver), so a dense cluster of 10 stressed neighbors delivered 10×
    // the concentration to a central receiver — flooding whole seedlings
    // with ethylene from a single stressed leaf.
    const float dr = g.ethylene_diffusion_radius;
    const float dr2 = dr * dr;

    // Snapshot the per-node emission before we overwrite.  Phase 1 left
    // each node's local().Ethylene = its own produced amount; we
    // redistribute these values and then write the redistributed totals
    // back.  Without the snapshot we'd double-count as sources become
    // receivers in the same pass.
    std::vector<float> emitted(nodes.size());
    for (size_t i = 0; i < nodes.size(); i++) {
        emitted[i] = nodes[i].node->local().chemical(ChemicalID::Ethylene);
    }

    std::vector<float> redistributed(nodes.size(), 0.0f);

    for (size_t i = 0; i < nodes.size(); i++) {
        const float src = emitted[i];
        if (src <= 0.0f) continue;

        // Pass 1 over neighbors: compute total weight (self = 1.0 + neighbors)
        float sum_weights = 1.0f;
        for (size_t j = 0; j < nodes.size(); j++) {
            if (i == j) continue;
            glm::vec3 diff = nodes[j].pos - nodes[i].pos;
            float dist2 = glm::dot(diff, diff);
            if (dist2 >= dr2) continue;
            float dist = std::sqrt(dist2);
            sum_weights += 1.0f - dist / dr;
        }

        // Pass 2: distribute src proportionally
        const float per_unit = src / sum_weights;
        redistributed[i] += per_unit; // self weight = 1.0
        for (size_t j = 0; j < nodes.size(); j++) {
            if (i == j) continue;
            glm::vec3 diff = nodes[j].pos - nodes[i].pos;
            float dist2 = glm::dot(diff, diff);
            if (dist2 >= dr2) continue;
            float dist = std::sqrt(dist2);
            const float w = 1.0f - dist / dr;
            redistributed[j] += per_unit * w;
        }
    }

    for (size_t i = 0; i < nodes.size(); i++) {
        nodes[i].node->local().chemical(ChemicalID::Ethylene) = redistributed[i];
    }
}

} // namespace botany
