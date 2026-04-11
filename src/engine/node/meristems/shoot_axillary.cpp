#include "engine/node/meristems/shoot_axillary.h"
#include "engine/node/meristems/shoot_apical.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include <algorithm>

namespace botany {

ShootAxillaryNode::ShootAxillaryNode(uint32_t id, glm::vec3 position, float radius)
    : MeristemNode(id, NodeType::SHOOT_AXILLARY, position, radius, false)
{}

void ShootAxillaryNode::tick(Plant& plant, const WorldParams& world) {
    MeristemNode::tick(plant, world);
    if (active) return; // already activated

    const Genome& g = plant.genome();

    // Sense auxin on the parent stem node
    float stem_auxin = parent ? parent->auxin : auxin;
    if (stem_auxin < g.auxin_threshold) {
        float parent_sugar = parent ? parent->sugar : sugar;
        if (parent_sugar < g.sugar_activation_shoot) return;
        if (sugar < world.sugar_cost_activation) return;

        // Create replacement shoot apical node
        Node* apical = plant.create_node(NodeType::SHOOT_APICAL, offset, g.initial_radius);
        apical->sugar = sugar - world.sugar_cost_activation;

        // Replace self in parent's children
        if (parent) {
            auto& siblings = parent->children;
            auto it = std::find(siblings.begin(), siblings.end(), static_cast<Node*>(this));
            if (it != siblings.end()) *it = apical;
            apical->parent = parent;
        }

        // Queue self for deferred removal
        parent = nullptr;
        plant.queue_removal(this);
    }
}

} // namespace botany
