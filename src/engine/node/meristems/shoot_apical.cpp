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

    // Roll direction once at the start of each internode
    if (target_internode_length < 1e-6f) roll_direction(g);

    if (!grow_tip(g, world)) return;
    if (parent) spawn_internode(plant, g);
}

void ShootApicalNode::roll_direction(const Genome& g) {
    growth_dir = perturb(growth_direction(*this), g.growth_noise);
}

bool ShootApicalNode::grow_tip(const Genome& g, const WorldParams& world) {
    float max_cost = g.growth_rate * world.sugar_cost_shoot_growth;
    float gf = sugar_growth_fraction(chemical(ChemicalID::Sugar), g.sugar_save_shoot, max_cost);
    if (gf < 1e-6f) return false;

    if (target_internode_length < 1e-6f) {
        target_internode_length = roll_internode_length(
            g.min_internode_length, g.max_internode_length, gf);
    }

    float actual_rate = g.growth_rate * gf;
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_shoot_growth;
    offset += growth_dir * actual_rate;
    return true;
}

void ShootApicalNode::spawn_internode(Plant& plant, const Genome& g) {
    float dist = glm::length(offset);
    if (dist <= target_internode_length) return;

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
    target_internode_length = 0.0f;
    ticks_since_last_node = 0;
    // Next call to grow() will roll a fresh direction
}

void ShootApicalNode::spawn_axillary(Plant& plant, Node* internode, const Genome& g, const glm::vec3& lateral_offset) {
    Node* axillary = plant.create_node(NodeType::SHOOT_AXILLARY, lateral_offset, g.initial_radius * 0.5f);
    internode->add_child(axillary);
}

void ShootApicalNode::spawn_leaf(Plant& plant, Node* internode, const Genome& g, const glm::vec3& lateral_offset) {
    Node* leaf = plant.create_node(NodeType::LEAF, lateral_offset, 0.0f);
    leaf->as_leaf()->leaf_size = g.leaf_bud_size;
    internode->add_child(leaf);
}

float ShootApicalNode::maintenance_cost(const Genome& g) const {
    return g.sugar_maintenance_meristem;
}

} // namespace botany
