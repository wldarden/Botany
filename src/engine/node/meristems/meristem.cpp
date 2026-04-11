// src/engine/meristems/meristem.cpp
#include "engine/node/meristems/meristem.h"
#include "engine/plant.h"
#include "engine/node/node.h"
#include "engine/world_params.h"

namespace botany {

static void tick_recursive(Node& node, Plant& plant, const WorldParams& world) {
    node.tick(plant, world);

    // Snapshot children: meristem ticks may reparent or add siblings
    auto children = node.children;
    for (Node* child : children) {
        tick_recursive(*child, plant, world);
    }
}

void tick_meristems(Plant& plant, const WorldParams& world) {
    tick_recursive(*plant.seed_mut(), plant, world);
    plant.flush_removals();
}

} // namespace botany
