#include "engine/node/meristems/root_axillary.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include <algorithm>

namespace botany {

RootAxillaryNode::RootAxillaryNode(uint32_t id, glm::vec3 position, float radius)
    : Node(id, NodeType::ROOT_AXILLARY, position, radius)
{}

void RootAxillaryNode::tick(Plant& plant, const WorldParams& world) {
    Node::tick(plant, world);
    if (active) return;

    const Genome& g = plant.genome();
    if (!can_activate(g, world)) return;
    activate(plant, g, world);
}

bool RootAxillaryNode::can_activate(const Genome& g, const WorldParams& world) const {
    float stem_cytokinin = parent ? parent->cytokinin : cytokinin;
    if (stem_cytokinin >= g.cytokinin_threshold) return false;

    float parent_sugar = parent ? parent->sugar : sugar;
    if (parent_sugar < g.sugar_activation_root) return false;
    if (sugar < world.sugar_cost_activation) return false;

    return true;
}

void RootAxillaryNode::activate(Plant& plant, const Genome& g, const WorldParams& world) {
    Node* apical = plant.create_node(NodeType::ROOT_APICAL, offset, g.root_initial_radius);
    apical->sugar = sugar - world.sugar_cost_activation;

    if (parent) {
        auto& siblings = parent->children;
        auto it = std::find(siblings.begin(), siblings.end(), static_cast<Node*>(this));
        if (it != siblings.end()) *it = apical;
        apical->parent = parent;
    }

    parent = nullptr;
    plant.queue_removal(this);
}

} // namespace botany
