// src/engine/meristems/meristem.cpp
#include "engine/meristems/meristem.h"
#include "engine/plant.h"
#include "engine/node.h"
#include "engine/world_params.h"
#include <vector>

namespace botany {

void tick_meristems(Plant& plant, const WorldParams& world) {
    // Collect meristem nodes first, since ticking may add new nodes
    std::vector<Node*> meristem_nodes;
    plant.for_each_node_mut([&](Node& n) {
        // Each node ticks itself (age, growth, etc.)
        n.tick(plant, world);

        if (n.meristem) {
            meristem_nodes.push_back(&n);
        }
    });

    // Meristem dispatch: apical/axillary growth
    for (Node* node : meristem_nodes) {
        if (node->meristem) {
            node->meristem->tick(*node, plant, world);
        }
    }
}

} // namespace botany
