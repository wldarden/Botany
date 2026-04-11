#include "engine/ethylene.h"
#include "engine/plant.h"
#include "engine/node.h"
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

    // Phase 1: Reset and compute local production
    for (auto& info : nodes) {
        Node& node = *info.node;
        node.ethylene = 0.0f;

        // Trigger 1: Sugar starvation
        if (node.sugar <= 0.0f) {
            node.ethylene += g.ethylene_starvation_rate;
        }

        // Trigger 2: Low light (LEAF only)
        if (node.type == NodeType::LEAF &&
            node.light_exposure < g.ethylene_shade_threshold) {
            node.ethylene += g.ethylene_shade_rate * (1.0f - node.light_exposure);
        }

        // Trigger 3: Old age (LEAF only)
        if (node.type == NodeType::LEAF &&
            node.age > g.ethylene_age_onset) {
            float age_past = static_cast<float>(node.age - g.ethylene_age_onset);
            node.ethylene += g.ethylene_age_rate * age_past
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
        node.ethylene += g.ethylene_crowding_rate * static_cast<float>(nearby_count);
    }

    // Phase 2: Spatial gas diffusion (compute-then-apply)
    std::vector<float> received(nodes.size(), 0.0f);
    float dr = g.ethylene_diffusion_radius;
    float dr2 = dr * dr;

    for (size_t i = 0; i < nodes.size(); i++) {
        if (nodes[i].node->ethylene <= 0.0f) continue;
        for (size_t j = 0; j < nodes.size(); j++) {
            if (i == j) continue;
            glm::vec3 diff = nodes[j].pos - nodes[i].pos;
            float dist2 = glm::dot(diff, diff);
            if (dist2 >= dr2) continue;
            float dist = std::sqrt(dist2);
            float falloff = 1.0f - dist / dr;
            received[j] += nodes[i].node->ethylene * falloff;
        }
    }

    for (size_t i = 0; i < nodes.size(); i++) {
        nodes[i].node->ethylene += received[i];
    }
}

void process_abscission(Plant& plant) {
    const Genome& g = plant.genome();

    // Increment senescence on senescing leaves, collect leaves to remove
    std::vector<Node*> to_remove;
    plant.for_each_node_mut([&](Node& node) {
        if (node.type != NodeType::LEAF) return;

        // Start senescence if ethylene exceeds threshold and not yet senescing
        if (node.senescence_ticks == 0 &&
            node.ethylene > g.ethylene_abscission_threshold) {
            node.senescence_ticks = 1;
        }

        // Advance senescence
        if (node.senescence_ticks > 0) {
            node.senescence_ticks++;
            if (node.senescence_ticks >= g.senescence_duration) {
                to_remove.push_back(&node);
            }
        }
    });

    for (Node* n : to_remove) {
        plant.remove_subtree(n);
    }
}

} // namespace botany
