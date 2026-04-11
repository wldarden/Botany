#include "engine/gibberellin.h"
#include "engine/plant.h"
#include "engine/node.h"

namespace botany {

void compute_gibberellin(Plant& plant) {
    const Genome& g = plant.genome();

    // Phase 1: Reset all GA to zero
    plant.for_each_node_mut([](Node& node) {
        node.gibberellin = 0.0f;
    });

    // Phase 2: Young leaves produce GA on their parent (and grandparent)
    plant.for_each_node_mut([&](Node& node) {
        if (node.type != NodeType::LEAF) return;
        if (node.age >= g.ga_leaf_age_max) return;
        if (node.leaf_size < 1e-6f) return;

        float production = node.leaf_size * g.ga_production_rate;

        // Apply to parent
        if (node.parent) {
            node.parent->gibberellin += production;

            // Apply reduced fraction to grandparent
            if (node.parent->parent) {
                node.parent->parent->gibberellin += production * 0.3f;
            }
        }
    });
}

} // namespace botany
