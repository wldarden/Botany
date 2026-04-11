// src/engine/meristems/root_apical.cpp
#include "engine/meristems/root_apical.h"
#include "engine/meristems/root_axillary.h"
#include "engine/meristems/helpers.h"
#include "engine/plant.h"
#include "engine/world_params.h"

namespace botany {

using namespace meristem_helpers;

void RootApicalMeristem::tick(Node& node, Plant& plant, const WorldParams& world) {
    if (!active) return;

    const Genome& g = plant.genome();

    glm::vec3 dir = perturb(growth_direction(node), g.growth_noise);

    // Gravitropism: roots sense gravity and turn downward when near or above
    // the soil surface. The closer to the surface, the stronger the correction.
    if (node.position.y > -g.root_gravitropism_depth) {
        float exposure = (node.position.y + g.root_gravitropism_depth)
                       / g.root_gravitropism_depth; // 0 at depth, 1+ at surface
        exposure = glm::clamp(exposure, 0.0f, 1.0f);
        float strength = exposure * g.root_gravitropism_strength;
        glm::vec3 down(0.0f, -1.0f, 0.0f);
        dir = glm::normalize(dir + down * strength);
    }

    // Sugar-scaled growth
    float max_cost = g.root_growth_rate * world.sugar_cost_root_growth;
    float gf = sugar_growth_fraction(node.sugar, g.sugar_save_root, max_cost);
    if (gf < 1e-6f) return;

    // Roll target internode length on first growth tick of this internode
    if (target_internode_length < 1e-6f) {
        target_internode_length = roll_internode_length(
            g.root_min_internode_length, g.root_max_internode_length, gf);
    }

    float actual_rate = g.root_growth_rate * gf;
    float actual_cost = actual_rate * world.sugar_cost_root_growth;
    node.sugar -= actual_cost;

    node.offset += dir * actual_rate;

    // Chain growth: same logic as shoot — interior node gets one root axillary
    if (node.parent) {
        float dist = glm::length(node.offset);
        if (dist > target_internode_length) {
            // This node becomes interior — give it a root axillary meristem
            // (only if we haven't hit the root meristem cap)
            if (!plant.root_meristems_at_cap()) {
                glm::vec3 branch_dir_val = branch_direction(dir, g.root_branch_angle, node.id);
                glm::vec3 ax_radial = branch_dir_val - dir * glm::dot(branch_dir_val, dir);
                float ax_rl = glm::length(ax_radial);
                if (ax_rl > 1e-4f) ax_radial /= ax_rl;
                glm::vec3 ax_offset = ax_radial * node.radius + branch_dir_val * g.tip_offset;
                Node* axillary = plant.create_node(NodeType::ROOT, ax_offset, g.root_initial_radius * 0.5f);
                axillary->meristem = plant.create_meristem<RootAxillaryMeristem>();
                node.add_child(axillary);
            }

            // Create new tip slightly ahead so it inherits the root's growth direction
            Node* new_tip = plant.create_node(NodeType::ROOT, dir * g.tip_offset, g.root_initial_radius);
            new_tip->meristem = node.meristem;
            target_internode_length = 0.0f;  // re-roll for next internode
            node.meristem = nullptr;
            node.add_child(new_tip);
        }
    }
}

} // namespace botany
