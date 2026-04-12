#include "engine/node/meristems/root_apical.h"
#include "engine/node/meristems/helpers.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include <glm/geometric.hpp>

namespace botany {

using namespace meristem_helpers;

RootApicalNode::RootApicalNode(uint32_t id, glm::vec3 position, float radius)
    : Node(id, NodeType::ROOT_APICAL, position, radius)
{}

void RootApicalNode::tick(Plant& plant, const WorldParams& world) {
    // Cytokinin production moved to leaves (proportional to sugar production)
    Node::tick(plant, world);
    ticks_since_last_node++;
}

void RootApicalNode::grow(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    // Roll direction once at the start of each internode
    if (target_internode_length < 1e-6f) roll_direction(g);

    if (!grow_tip(g, world)) return;
    if (parent) spawn_internode(plant, g);
}

void RootApicalNode::roll_direction(const Genome& g) {
    growth_dir = apply_gravitropism(
        perturb(growth_direction(*this), g.growth_noise), g);
}

glm::vec3 RootApicalNode::apply_gravitropism(const glm::vec3& dir, const Genome& g) const {
    if (position.y <= -g.root_gravitropism_depth) return dir;

    float exposure = (position.y + g.root_gravitropism_depth)
                   / g.root_gravitropism_depth;
    exposure = glm::clamp(exposure, 0.0f, 1.0f);
    float strength = exposure * g.root_gravitropism_strength;
    glm::vec3 down(0.0f, -1.0f, 0.0f);
    return glm::normalize(dir + down * strength);
}

bool RootApicalNode::grow_tip(const Genome& g, const WorldParams& world) {
    float max_cost = g.root_growth_rate * world.sugar_cost_root_growth;
    float gf = growth_fraction(chemical(ChemicalID::Sugar), max_cost,
                               chemical(ChemicalID::Cytokinin), g.cytokinin_growth_threshold);
    if (gf < 1e-6f) return false;

    if (target_internode_length < 1e-6f) {
        target_internode_length = roll_internode_length(
            g.root_min_internode_length, g.root_max_internode_length, gf);
    }

    float actual_rate = g.root_growth_rate * gf;
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_root_growth;
    offset += growth_dir * actual_rate;
    return true;
}

void RootApicalNode::spawn_internode(Plant& plant, const Genome& g) {
    float dist = glm::length(offset);
    if (dist <= target_internode_length) return;

    // Create new interior root node and insert it between us and our parent
    Node* internode = plant.create_node(NodeType::ROOT, offset, radius);
    parent->replace_child(this, internode);
    offset = growth_dir * g.tip_offset;
    internode->add_child(this);

    // Lateral branching: compute offset and spawn axillary bud
    if (!plant.root_meristems_at_cap()) {
        glm::vec3 branch_dir_val = branch_direction(growth_dir, g.root_branch_angle, id);
        glm::vec3 ax_radial = branch_dir_val - growth_dir * glm::dot(branch_dir_val, growth_dir);
        float ax_rl = glm::length(ax_radial);
        if (ax_rl > 1e-4f) ax_radial /= ax_rl;
        glm::vec3 lateral_offset = ax_radial * internode->radius + branch_dir_val * g.tip_offset;
        spawn_axillary(plant, internode, g, lateral_offset);
    }

    target_internode_length = 0.0f;
    ticks_since_last_node = 0;
    // Next call to grow() will roll a fresh direction
}

void RootApicalNode::spawn_axillary(Plant& plant, Node* internode, const Genome& g, const glm::vec3& lateral_offset) {
    Node* axillary = plant.create_node(NodeType::ROOT_AXILLARY, lateral_offset, g.root_initial_radius * 0.5f);
    internode->add_child(axillary);
}

float RootApicalNode::maintenance_cost(const Genome& g) const {
    return g.sugar_maintenance_meristem;
}

} // namespace botany
