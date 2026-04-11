#include "engine/gibberellin.h"
#include "engine/plant.h"
#include "engine/node/node.h"
#include "engine/node/leaf_node.h"

namespace botany {

void compute_gibberellin(Plant& plant) {
    const Genome& g = plant.genome();

    // Phase 1: Reset all GA to zero
    plant.for_each_node_mut([](Node& node) {
        node.chemical(ChemicalID::Gibberellin) = 0.0f;
        node.gibberellin = 0.0f;
    });

    // Phase 2: Young leaves produce GA on their parent (and grandparent)
    plant.for_each_node_mut([&](Node& node) {
        auto* leaf = node.as_leaf();
        if (!leaf) return;
        if (node.age >= g.ga_leaf_age_max) return;
        if (leaf->leaf_size < 1e-6f) return;

        float production = leaf->leaf_size * g.ga_production_rate;

        // Apply to parent
        if (node.parent) {
            node.parent->chemical(ChemicalID::Gibberellin) += production;
            node.parent->gibberellin = node.parent->chemical(ChemicalID::Gibberellin);

            // Apply reduced fraction to grandparent
            if (node.parent->parent) {
                node.parent->parent->chemical(ChemicalID::Gibberellin) += production * 0.3f;
                node.parent->parent->gibberellin = node.parent->parent->chemical(ChemicalID::Gibberellin);
            }
        }
    });
}

} // namespace botany
