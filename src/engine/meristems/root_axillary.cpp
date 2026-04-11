// src/engine/meristems/root_axillary.cpp
#include "engine/meristems/root_axillary.h"
#include "engine/meristems/root_apical.h"
#include "engine/plant.h"
#include "engine/world_params.h"

namespace botany {

void RootAxillaryMeristem::tick(Node& node, Plant& plant, const WorldParams& world) {
    if (active) return;

    const Genome& g = plant.genome();
    // Sense cytokinin on the parent root node, not this side-branch node
    float stem_cytokinin = node.parent ? node.parent->cytokinin : node.cytokinin;
    if (stem_cytokinin < g.cytokinin_threshold) {
        // Parent root must have sufficient sugar supply
        float parent_sugar = node.parent ? node.parent->sugar : node.sugar;
        if (parent_sugar < g.sugar_activation_root) return;

        // Activation costs sugar
        if (node.sugar < world.sugar_cost_activation) return;
        node.sugar -= world.sugar_cost_activation;

        // Replace this dormant root axillary with an active root apical
        auto* apical = plant.create_meristem<RootApicalMeristem>();
        apical->ticks_since_last_node = 0;
        node.meristem = apical;
        node.radius = g.root_initial_radius;
    }
}

} // namespace botany
