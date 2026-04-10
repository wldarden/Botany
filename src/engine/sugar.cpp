#include "engine/sugar.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include "engine/node.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace botany {

void produce_sugar(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();
    plant.for_each_node_mut([&](Node& node) {
        if (node.type == NodeType::LEAF) {
            node.sugar += world.light_level * node.leaf_size * g.sugar_production_rate;
        }
    });
}

void consume_sugar(Plant& plant) {
    const Genome& g = plant.genome();
    plant.for_each_node_mut([&](Node& node) {
        float cost = 0.0f;
        switch (node.type) {
            case NodeType::LEAF:
                cost = g.sugar_maintenance_leaf * node.leaf_size;
                break;
            case NodeType::STEM:
                cost = g.sugar_maintenance_stem * node.radius;
                break;
            case NodeType::ROOT:
                cost = g.sugar_maintenance_root * node.radius;
                break;
        }
        if (node.meristem && node.meristem->is_tip() && node.meristem->active) {
            cost += g.sugar_maintenance_meristem;
        }
        node.sugar = std::max(0.0f, node.sugar - cost);

        // Track starvation
        if (node.sugar <= 0.0f) {
            node.starvation_ticks++;
        } else {
            node.starvation_ticks = 0;
        }
    });
}

void diffuse_sugar(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    // Build edge list once — all parent-child connections
    struct Edge {
        Node* a;
        Node* b;
        float capacity;
    };
    std::vector<Edge> edges;
    plant.for_each_node_mut([&](Node& node) {
        if (node.parent) {
            float min_radius = std::min(node.radius, node.parent->radius);
            // LEAF nodes have radius 0, use a small baseline for leaf connections
            if (node.type == NodeType::LEAF || node.parent->type == NodeType::LEAF) {
                min_radius = std::max(min_radius, 0.01f);
            }
            float capacity = min_radius * min_radius * 3.14159f * g.sugar_transport_conductance;
            edges.push_back({&node, node.parent, capacity});
        }
    });

    // Run multiple diffusion iterations
    for (int iter = 0; iter < world.sugar_diffusion_iterations; iter++) {
        // Compute all flows first, then apply (avoids order-dependent artifacts)
        struct Flow {
            Node* from;
            Node* to;
            float amount;
        };
        std::vector<Flow> flows;

        for (const Edge& e : edges) {
            float gradient = e.a->sugar - e.b->sugar;
            if (std::abs(gradient) < 1e-6f) continue;

            float flow = gradient * e.capacity;
            if (flow > 0.0f) {
                flow = std::min(flow, e.a->sugar);
                flows.push_back({e.a, e.b, flow});
            } else {
                flow = std::min(-flow, e.b->sugar);
                flows.push_back({e.b, e.a, flow});
            }
        }

        for (const Flow& f : flows) {
            f.from->sugar -= f.amount;
            f.to->sugar += f.amount;
        }
    }
}

void prune_starved_nodes(Plant& plant, const WorldParams& world) {
    std::vector<Node*> to_prune;
    plant.for_each_node_mut([&](Node& node) {
        if (node.starvation_ticks >= world.starvation_ticks_max) {
            // Only prune if parent is NOT also starved (prune from the top)
            bool parent_starved = node.parent &&
                node.parent->starvation_ticks >= world.starvation_ticks_max;
            if (!parent_starved && node.parent != nullptr) {  // never prune seed (parent == nullptr)
                to_prune.push_back(&node);
            }
        }
    });

    for (Node* n : to_prune) {
        plant.remove_subtree(n);
    }
}

void transport_sugar(Plant& plant, const WorldParams& world) {
    produce_sugar(plant, world);
    diffuse_sugar(plant, world);
    consume_sugar(plant);
    prune_starved_nodes(plant, world);
}

} // namespace botany
