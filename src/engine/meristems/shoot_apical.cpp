// src/engine/meristems/shoot_apical.cpp
#include "engine/meristems/shoot_apical.h"
#include "engine/meristems/shoot_axillary.h"
#include "engine/meristems/helpers.h"
#include "engine/plant.h"
#include "engine/world_params.h"

namespace botany {

using namespace meristem_helpers;

void ShootApicalMeristem::tick(Node& node, Plant& plant, const WorldParams& world) {
    if (!active) return;

    const Genome& g = plant.genome();
    glm::vec3 dir = perturb(growth_direction(node), g.growth_noise);

    if (!grow(node, g, world, dir)) return;
    if (node.parent) split_internode(node, plant, g, dir);
}

bool ShootApicalMeristem::grow(Node& node, const Genome& g, const WorldParams& world, const glm::vec3& dir) {
    float max_cost = g.growth_rate * world.sugar_cost_shoot_growth;
    float gf = sugar_growth_fraction(node.sugar, g.sugar_save_shoot, max_cost);
    if (gf < 1e-6f) return false;

    // Roll target internode length on first growth tick of this internode
    if (target_internode_length < 1e-6f) {
        target_internode_length = roll_internode_length(
            g.min_internode_length, g.max_internode_length, gf);
    }

    float actual_rate = g.growth_rate * gf;
    node.sugar -= actual_rate * world.sugar_cost_shoot_growth;
    node.offset += dir * actual_rate;
    return true;
}

void ShootApicalMeristem::split_internode(Node& node, Plant& plant, const Genome& g, const glm::vec3& dir) {
    float dist = glm::length(node.offset);
    if (dist <= target_internode_length) return;

    // Phyllotaxis: axillary bud and leaf placed at golden-angle offset
    glm::vec3 branch_dir = branch_direction(dir, g.branch_angle, phyllotaxis_index);
    glm::vec3 radial = branch_dir - dir * glm::dot(branch_dir, dir);
    float rl = glm::length(radial);
    if (rl > 1e-4f) radial /= rl;
    glm::vec3 lateral_offset = radial * node.radius + branch_dir * g.tip_offset;

    // Axillary bud
    Node* axillary = plant.create_node(NodeType::STEM, lateral_offset, g.initial_radius * 0.5f);
    axillary->meristem = plant.create_meristem<ShootAxillaryMeristem>();
    node.add_child(axillary);

    // Leaf at same angular position (leaf axil)
    Node* leaf_node = plant.create_node(NodeType::LEAF, lateral_offset, 0.0f);
    leaf_node->as_leaf()->leaf_size = g.leaf_bud_size;
    node.add_child(leaf_node);

    phyllotaxis_index++;

    // New tip carries meristem forward
    Node* new_tip = plant.create_node(NodeType::STEM, dir * g.tip_offset, g.initial_radius);
    new_tip->meristem = node.meristem;
    target_internode_length = 0.0f;
    node.meristem = nullptr;
    node.add_child(new_tip);
}

} // namespace botany
