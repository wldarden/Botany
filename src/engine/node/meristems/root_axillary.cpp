#include "engine/node/meristems/root_axillary.h"
#include "engine/plant.h"
#include "engine/sugar.h"
#include "engine/world_params.h"

namespace botany {

RootAxillaryNode::RootAxillaryNode(uint32_t id, glm::vec3 position, float radius)
    : Node(id, NodeType::ROOT_AXILLARY, position, radius)
{}

void RootAxillaryNode::tick(Plant& plant, const WorldParams& world) {
    if (active) {
        Node::tick(plant, world);
        return;
    }

    // Dormant fast path: skip mass/stress/droop (no children, tiny node).
    // Just update position, age, transport (for sugar accumulation), and check activation.
    position = parent ? parent->position + offset : offset;
    age++;
    const Genome& g = plant.genome();
    transport_chemicals(g);

    // Mass bookkeeping (constant, no children)
    total_mass = world.meristem_mass;
    mass_moment = total_mass * position;

    if (!can_activate(g, world)) return;
    activate(plant, g, world);
}

bool RootAxillaryNode::can_activate(const Genome& g, const WorldParams& world) const {
    // Cytokinin from producing leaves signals "the plant can support new roots"
    float local_cyt = parent ? parent->chemical(ChemicalID::Cytokinin) : chemical(ChemicalID::Cytokinin);
    if (local_cyt < g.cytokinin_threshold) return false;

    if (chemical(ChemicalID::Sugar) < world.sugar_cost_activation) return false;

    return true;
}

void RootAxillaryNode::activate(Plant& plant, const Genome& g, const WorldParams& world) {
    Node* apical = plant.create_node(NodeType::ROOT_APICAL, offset, g.root_initial_radius);
    apical->chemical(ChemicalID::Sugar) = chemical(ChemicalID::Sugar) - world.sugar_cost_activation;

    if (parent) {
        Node* stem = parent;  // save before replace_child nulls this->parent
        stem->replace_child(this, apical);
        apical->position = stem->position + apical->offset;
    }
    plant.queue_removal(this);
}

} // namespace botany
