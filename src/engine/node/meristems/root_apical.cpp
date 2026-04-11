#include "engine/node/meristems/root_apical.h"
#include "engine/node/meristems/helpers.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include <algorithm>
#include <glm/geometric.hpp>

namespace botany {

using namespace meristem_helpers;

RootApicalNode::RootApicalNode(uint32_t id, glm::vec3 position, float radius)
    : Node(id, NodeType::ROOT_APICAL, position, radius)
{}

void RootApicalNode::tick(Plant& plant, const WorldParams& world) {
    cytokinin += plant.genome().cytokinin_production_rate;

    Node::tick(plant, world);
    ticks_since_last_node++;

    const Genome& g = plant.genome();
    glm::vec3 dir = apply_gravitropism(
        perturb(growth_direction(*this), g.growth_noise), g);

    if (!grow(g, world, dir)) return;
    if (parent) split_internode(plant, g, dir);
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

bool RootApicalNode::grow(const Genome& g, const WorldParams& world, const glm::vec3& dir) {
    float max_cost = g.root_growth_rate * world.sugar_cost_root_growth;
    float gf = sugar_growth_fraction(sugar, g.sugar_save_root, max_cost);
    if (gf < 1e-6f) return false;

    if (target_internode_length < 1e-6f) {
        target_internode_length = roll_internode_length(
            g.root_min_internode_length, g.root_max_internode_length, gf);
    }

    float actual_rate = g.root_growth_rate * gf;
    sugar -= actual_rate * world.sugar_cost_root_growth;
    offset += dir * actual_rate;
    return true;
}

void RootApicalNode::split_internode(Plant& plant, const Genome& g, const glm::vec3& dir) {
    float dist = glm::length(offset);
    if (dist <= target_internode_length) return;

    Node* my_parent = parent;

    Node* internode = plant.create_node(NodeType::ROOT, offset, radius);

    auto& siblings = my_parent->children;
    auto it = std::find(siblings.begin(), siblings.end(), static_cast<Node*>(this));
    if (it != siblings.end()) *it = internode;
    internode->parent = my_parent;

    parent = nullptr;
    offset = dir * g.tip_offset;
    internode->add_child(this);

    if (!plant.root_meristems_at_cap()) {
        glm::vec3 branch_dir_val = branch_direction(dir, g.root_branch_angle, id);
        glm::vec3 ax_radial = branch_dir_val - dir * glm::dot(branch_dir_val, dir);
        float ax_rl = glm::length(ax_radial);
        if (ax_rl > 1e-4f) ax_radial /= ax_rl;
        glm::vec3 ax_offset = ax_radial * internode->radius + branch_dir_val * g.tip_offset;
        Node* axillary = plant.create_node(NodeType::ROOT_AXILLARY, ax_offset, g.root_initial_radius * 0.5f);
        internode->add_child(axillary);
    }

    target_internode_length = 0.0f;
    ticks_since_last_node = 0;
}

float RootApicalNode::maintenance_cost(const Genome& g) const {
    return g.sugar_maintenance_meristem;
}

} // namespace botany
