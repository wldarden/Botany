#include "engine/sugar.h"
#include "engine/node/node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/genome.h"
#include <algorithm>
#include <glm/geometric.hpp>

namespace botany {

// Seed node (no parent) is a root-shoot junction. Its cap should match its
// children's scale so concentration gradients work naturally — not so small
// that it bottlenecks, not so large that it becomes a concentration sink.
// We sum its children's caps to make it proportional to the connected network.
static float seed_resource_cap(const Node& seed, float per_node_cap_fn(const Node&, const Genome&), const Genome& g) {
    float total = 0.0f;
    for (const Node* child : seed.children) {
        total += per_node_cap_fn(*child, g);
    }
    return std::max(total, g.sugar_cap_minimum);
}


float sugar_cap(const Node& node, const Genome& g) {
    // Seed: cap = sum of children's caps (scales with connected network),
    // but always at least seed_sugar so the initial reserve fits.
    if (!node.parent) return std::max(seed_resource_cap(node, sugar_cap, g), g.seed_sugar);

    float length = std::max(glm::length(node.offset), 0.01f);

    switch (node.type) {
        case NodeType::LEAF: {
            float ls = node.as_leaf()->leaf_size;
            float area = ls * ls;
            return std::max(g.sugar_cap_minimum, area * g.sugar_storage_density_leaf);
        }
        case NodeType::STEM:
        case NodeType::ROOT: {
            float volume = 3.14159f * node.radius * node.radius * length;
            return std::max(g.sugar_cap_minimum, volume * g.sugar_storage_density_wood);
        }
        case NodeType::APICAL:
        case NodeType::ROOT_APICAL:
            return g.sugar_cap_meristem;
    }
    return g.sugar_cap_minimum;
}

float water_cap(const Node& node, const Genome& g) {
    // Seed: cap = sum of children's caps (scales with connected network)
    if (!node.parent) return seed_resource_cap(node, water_cap, g);

    float length = std::max(glm::length(node.offset), 0.01f);

    switch (node.type) {
        case NodeType::LEAF: {
            float ls = node.as_leaf()->leaf_size;
            float area = ls * ls;
            return std::max(g.sugar_cap_minimum, area * g.water_storage_density_leaf);
        }
        case NodeType::STEM:
        case NodeType::ROOT: {
            float volume = 3.14159f * node.radius * node.radius * length;
            return std::max(g.sugar_cap_minimum, volume * g.water_storage_density_stem);
        }
        case NodeType::APICAL:
        case NodeType::ROOT_APICAL:
            return g.water_cap_meristem;
    }
    return g.sugar_cap_minimum;
}

} // namespace botany
