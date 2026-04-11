#include "engine/node/meristems/root_axillary.h"
#include "engine/node/meristems/root_apical.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include <algorithm>

namespace botany {

RootAxillaryNode::RootAxillaryNode(uint32_t id, glm::vec3 position, float radius)
    : MeristemNode(id, NodeType::ROOT_AXILLARY, position, radius, false)
{}

void RootAxillaryNode::tick(Plant& plant, const WorldParams& world) {
    MeristemNode::tick(plant, world);
    if (active) return;

    const Genome& g = plant.genome();

    // Sense cytokinin on the parent root node
    float stem_cytokinin = parent ? parent->cytokinin : cytokinin;
    if (stem_cytokinin < g.cytokinin_threshold) {
        float parent_sugar = parent ? parent->sugar : sugar;
        if (parent_sugar < g.sugar_activation_root) return;
        if (sugar < world.sugar_cost_activation) return;

        // Create replacement root apical node
        Node* apical = plant.create_node(NodeType::ROOT_APICAL, offset, g.root_initial_radius);
        apical->sugar = sugar - world.sugar_cost_activation;

        // Replace self in parent's children
        if (parent) {
            auto& siblings = parent->children;
            auto it = std::find(siblings.begin(), siblings.end(), static_cast<Node*>(this));
            if (it != siblings.end()) *it = apical;
            apical->parent = parent;
        }

        parent = nullptr;
        plant.queue_removal(this);
    }
}

} // namespace botany
