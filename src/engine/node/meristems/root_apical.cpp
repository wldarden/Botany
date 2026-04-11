#include "engine/node/meristems/root_apical.h"
#include "engine/node/meristems/helpers.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include <algorithm>
#include <glm/geometric.hpp>

namespace botany {

using namespace meristem_helpers;

RootApicalNode::RootApicalNode(uint32_t id, glm::vec3 position, float radius)
    : MeristemNode(id, NodeType::ROOT_APICAL, position, radius, true)
{}

void RootApicalNode::tick(Plant& plant, const WorldParams& world) {
    // Produce cytokinin before base tick (so transport moves it this frame)
    if (active) cytokinin += plant.genome().cytokinin_production_rate;

    MeristemNode::tick(plant, world);
    if (!active) return;

    const Genome& g = plant.genome();
    glm::vec3 dir = perturb(growth_direction(*this), g.growth_noise);

    // Gravitropism: roots turn downward near the surface
    if (position.y > -g.root_gravitropism_depth) {
        float exposure = (position.y + g.root_gravitropism_depth)
                       / g.root_gravitropism_depth;
        exposure = glm::clamp(exposure, 0.0f, 1.0f);
        float strength = exposure * g.root_gravitropism_strength;
        glm::vec3 down(0.0f, -1.0f, 0.0f);
        dir = glm::normalize(dir + down * strength);
    }

    // Sugar-scaled growth
    float max_cost = g.root_growth_rate * world.sugar_cost_root_growth;
    float gf = sugar_growth_fraction(sugar, g.sugar_save_root, max_cost);
    if (gf < 1e-6f) return;

    if (target_internode_length < 1e-6f) {
        target_internode_length = roll_internode_length(
            g.root_min_internode_length, g.root_max_internode_length, gf);
    }

    float actual_rate = g.root_growth_rate * gf;
    sugar -= actual_rate * world.sugar_cost_root_growth;
    offset += dir * actual_rate;

    // Chain growth: insert internode when long enough
    if (!parent) return;
    float dist = glm::length(offset);
    if (dist <= target_internode_length) return;

    Node* my_parent = parent;

    // Create new interior root node
    Node* internode = plant.create_node(NodeType::ROOT, offset, radius);

    // Replace this meristem in parent's children
    auto& siblings = my_parent->children;
    auto it = std::find(siblings.begin(), siblings.end(), static_cast<Node*>(this));
    if (it != siblings.end()) *it = internode;
    internode->parent = my_parent;

    // Re-attach this meristem as child of internode
    parent = nullptr;
    offset = dir * g.tip_offset;
    internode->add_child(this);

    // Root axillary bud (if not at cap)
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

} // namespace botany
