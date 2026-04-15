#include "engine/sugar.h"
#include "engine/node/node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/genome.h"
#include <algorithm>
#include <glm/geometric.hpp>

namespace botany {

float sugar_cap(const Node& node, const Genome& g) {
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
            float cap = volume * g.sugar_storage_density_wood;
            if (!node.parent) {
                cap = std::max(cap, g.seed_sugar);
            }
            return std::max(g.sugar_cap_minimum, cap);
        }
        case NodeType::APICAL:
        case NodeType::ROOT_APICAL:
            return g.sugar_cap_meristem;
    }
    return g.sugar_cap_minimum;
}

} // namespace botany
