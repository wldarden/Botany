#include "engine/node/meristems/shoot_axillary.h"
#include "engine/node/meristems/shoot_apical.h"
#include "engine/plant.h"
#include "engine/sugar.h"
#include "engine/world_params.h"
#include <glm/geometric.hpp>

namespace botany {

ShootAxillaryNode::ShootAxillaryNode(uint32_t id, glm::vec3 position, float radius)
    : Node(id, NodeType::SHOOT_AXILLARY, position, radius)
{}

void ShootAxillaryNode::tick(Plant& plant, const WorldParams& world) {
    Node::tick(plant, world);
    if (active) return;

    const Genome& g = plant.genome();
    if (!can_activate(g, world)) return;
    activate(plant, g, world);
}

bool ShootAxillaryNode::can_activate(const Genome& g, const WorldParams& world) const {
    // Low auxin removes inhibition (apical dominance weakened)
    float stem_auxin = parent ? parent->chemical(ChemicalID::Auxin) : chemical(ChemicalID::Auxin);
    if (stem_auxin >= g.auxin_threshold) return false;

    // Cytokinin from producing leaves signals "the plant can support a new branch."
    // Sugar can't gate this — it's everywhere during seed-funded growth.
    float local_cyt = parent ? parent->chemical(ChemicalID::Cytokinin) : chemical(ChemicalID::Cytokinin);
    if (local_cyt < g.cytokinin_threshold) return false;

    if (chemical(ChemicalID::Sugar) < world.sugar_cost_activation) return false;

    return true;
}

void ShootAxillaryNode::activate(Plant& plant, const Genome& g, const WorldParams& world) {
    Node* apical = plant.create_node(NodeType::SHOOT_APICAL, offset, g.initial_radius);
    apical->chemical(ChemicalID::Sugar) = chemical(ChemicalID::Sugar) - world.sugar_cost_activation;

    // Set the branch's initial growth direction from the offset angle.
    // Without this, the SA re-rolls direction from its tiny offset and
    // quickly converges to vertical, making all branches parallel the trunk.
    float olen = glm::length(offset);
    if (olen > 1e-4f) {
        apical->as_shoot_apical()->growth_dir = offset / olen;
    }

    if (parent) {
        Node* stem = parent;  // save before replace_child nulls this->parent
        stem->replace_child(this, apical);
        apical->position = stem->position + apical->offset;
    }
    plant.queue_removal(this);
}

} // namespace botany
