// src/engine/meristems/shoot_axillary.cpp
#include "engine/meristems/shoot_axillary.h"
#include "engine/meristems/shoot_apical.h"
#include "engine/plant.h"
#include "engine/world_params.h"

namespace botany {

void ShootAxillaryMeristem::tick(Node& node, Plant& plant, const WorldParams& world) {
    if (active) return; // already activated — nothing to do

    const Genome& g = plant.genome();
    // Sense auxin on the parent stem node (where auxin actually flows),
    // not on this side-branch node which never receives basipetal auxin.
    // Shoot axillaries activate when auxin drops low enough —
    // meaning they're far from any active shoot tip.
    float stem_auxin = node.parent ? node.parent->auxin : node.auxin;
    if (stem_auxin < g.auxin_threshold) {
        // Parent stem must be well-fed (proxy for light reaching this branch)
        float parent_sugar = node.parent ? node.parent->sugar : node.sugar;
        if (parent_sugar < g.sugar_activation_shoot) return;

        // Activation costs sugar
        if (node.sugar < world.sugar_cost_activation) return;
        node.sugar -= world.sugar_cost_activation;

        // Replace this dormant axillary with an active shoot apical
        auto* apical = plant.create_meristem<ShootApicalMeristem>();
        apical->ticks_since_last_node = 0;
        node.meristem = apical;
        node.radius = g.initial_radius;
    }
}

} // namespace botany
