#include "engine/node/meristems/root_axillary.h"
#include "engine/plant.h"
#include "engine/world_params.h"

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
    float stem_cytokinin = parent ? parent->chemical(ChemicalID::Cytokinin) : chemical(ChemicalID::Cytokinin);
    if (stem_cytokinin >= g.cytokinin_threshold) return false;

    float parent_sugar_val = parent ? parent->chemical(ChemicalID::Sugar) : chemical(ChemicalID::Sugar);
    if (parent_sugar_val < g.sugar_activation_root) return false;
    if (chemical(ChemicalID::Sugar) < world.sugar_cost_activation) return false;

    return true;
}

void RootAxillaryNode::activate(Plant& plant, const Genome& g, const WorldParams& world) {
    Node* apical = plant.create_node(NodeType::ROOT_APICAL, offset, g.root_initial_radius);
    apical->chemical(ChemicalID::Sugar) = chemical(ChemicalID::Sugar) - world.sugar_cost_activation;

    if (parent) {
        parent->replace_child(this, apical);
    }
    plant.queue_removal(this);
}

} // namespace botany
