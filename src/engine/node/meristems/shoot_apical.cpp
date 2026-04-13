#include "engine/node/meristems/shoot_apical.h"
#include "engine/node/meristems/helpers.h"
#include "engine/plant.h"
#include "engine/node/leaf_node.h"
#include "engine/world_params.h"

namespace botany {

using namespace meristem_helpers;

ShootApicalNode::ShootApicalNode(uint32_t id, glm::vec3 position, float radius)
    : Node(id, NodeType::SHOOT_APICAL, position, radius)
{}

void ShootApicalNode::tick(Plant& plant, const WorldParams& world) {
    chemical(ChemicalID::Auxin) += plant.genome().auxin_production_rate;
    Node::tick(plant, world);
    ticks_since_last_node++;
}

void ShootApicalNode::grow(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    grow_tip(g, world);

    // Time-based spawning (plastochron): create nodes at regular intervals,
    // like real meristems. Internode length comes from elongation afterward.
    // Don't spawn if starving — a meristem with no sugar shouldn't create nodes it can't feed.
    if (parent && ticks_since_last_node >= g.shoot_plastochron && starvation_ticks == 0) {
        spawn_internode(plant, g);
    }
}

void ShootApicalNode::roll_direction(const Genome& g) {
    growth_dir = perturb(growth_direction(*this), g.growth_noise);
}

void ShootApicalNode::grow_tip(const Genome& g, const WorldParams& world) {
    float max_cost = g.growth_rate * world.sugar_cost_shoot_growth;
    float gf = growth_fraction(chemical(ChemicalID::Sugar), max_cost,
                               chemical(ChemicalID::Cytokinin), g.cytokinin_growth_threshold);
    if (gf < 1e-6f) return;

    // Roll a fresh direction each tick (small perturbation from parent direction)
    if (glm::length(growth_dir) < 1e-4f) roll_direction(g);

    // Stress hormone pulls growth direction toward vertical
    float stress_grav = chemical(ChemicalID::Stress) * g.stress_gravitropism_boost;
    if (stress_grav > 1e-6f) {
        glm::vec3 up(0.0f, 1.0f, 0.0f);
        float blend = std::min(stress_grav, 0.5f);  // cap at 50% pull toward vertical
        growth_dir = glm::normalize(glm::mix(growth_dir, up, blend));
    }

    float actual_rate = g.growth_rate * gf;
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_shoot_growth;
    offset += growth_dir * actual_rate;
}

void ShootApicalNode::spawn_internode(Plant& plant, const Genome& g) {

    // Create new interior stem node and insert it between us and our parent
    Node* internode = plant.create_node(NodeType::STEM, offset, radius);
    parent->replace_child(this, internode);
    offset = growth_dir * g.tip_offset;
    internode->add_child(this);

    // Phyllotaxis: compute lateral offset at golden-angle rotation
    glm::vec3 branch_dir = branch_direction(growth_dir, g.branch_angle, phyllotaxis_index);
    glm::vec3 radial = branch_dir - growth_dir * glm::dot(branch_dir, growth_dir);
    float rl = glm::length(radial);
    if (rl > 1e-4f) radial /= rl;
    glm::vec3 lateral_offset = radial * internode->radius + branch_dir * g.tip_offset;

    spawn_axillary(plant, internode, g, lateral_offset);
    spawn_leaf(plant, internode, g, lateral_offset);

    phyllotaxis_index++;
    ticks_since_last_node = 0;
    growth_dir = glm::vec3(0.0f); // re-roll direction for next internode
}

void ShootApicalNode::spawn_axillary(Plant& plant, Node* internode, const Genome& g, const glm::vec3& lateral_offset) {
    Node* axillary = plant.create_node(NodeType::SHOOT_AXILLARY, lateral_offset, g.initial_radius * 0.5f);
    internode->add_child(axillary);
}

void ShootApicalNode::spawn_leaf(Plant& plant, Node* internode, const Genome& g, const glm::vec3& lateral_offset) {
    // Start leaf close to stem — petiole grows with leaf size
    Node* leaf = plant.create_node(NodeType::LEAF, lateral_offset, 0.0f);
    leaf->as_leaf()->leaf_size = g.leaf_bud_size;
    float len = glm::length(lateral_offset);
    if (len > 1e-4f) {
        leaf->as_leaf()->facing = lateral_offset / len;
    }
    // Meristem gives the new leaf some sugar to bootstrap it
    float gift = std::min(chemical(ChemicalID::Sugar) * 0.1f, 0.5f);
    chemical(ChemicalID::Sugar) -= gift;
    leaf->chemical(ChemicalID::Sugar) = gift;
    internode->add_child(leaf);
}

float ShootApicalNode::maintenance_cost(const Genome& g) const {
    return g.sugar_maintenance_meristem;
}

} // namespace botany
