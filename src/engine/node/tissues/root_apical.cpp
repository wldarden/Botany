#include "engine/node/tissues/root_apical.h"
#include "engine/node/meristems/helpers.h"
#include "engine/plant.h"
#include "engine/sugar.h"
#include "engine/world_params.h"
#include <glm/geometric.hpp>

namespace botany {

using namespace meristem_helpers;

RootApicalNode::RootApicalNode(uint32_t id, glm::vec3 position, float radius)
    : Node(id, NodeType::ROOT_APICAL, position, radius)
{}

void RootApicalNode::update_tissue(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    absorb_water(g, world);

    if (!active) {
        if (can_activate(g, world)) activate(g, world);
        return;
    }

    ticks_since_last_node++;
    elongate(g, world);
    check_spawn(plant, g);
}

void RootApicalNode::absorb_water(const Genome& g, const WorldParams& world) {
    // Hemisphere approximation for root tip surface area
    float surface_area = 2.0f * 3.14159f * radius * radius;
    float absorbed = g.water_absorption_rate * surface_area * world.soil_moisture;
    float cap = water_cap(*this, g);
    chemical(ChemicalID::Water) = std::min(chemical(ChemicalID::Water) + absorbed, cap);
}

void RootApicalNode::check_spawn(Plant& plant, const Genome& g) {
    if (parent && ticks_since_last_node >= g.root_plastochron && starvation_ticks == 0) {
        spawn_internode(plant, g);
    }
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

void RootApicalNode::elongate(const Genome& g, const WorldParams& world) {
    float max_cost = g.root_growth_rate * world.sugar_cost_root_growth;
    float gf = growth_fraction(chemical(ChemicalID::Sugar), max_cost,
                               chemical(ChemicalID::Cytokinin), g.cytokinin_growth_threshold);
    if (gf < 1e-6f) return;

    if (glm::length(growth_dir) < 1e-4f) roll_direction(g);

    float auxin_boost = meristem_helpers::auxin_growth_factor(
        chemical(ChemicalID::Auxin), g.root_apical_auxin_max_boost, g.root_apical_auxin_half_saturation);
    float actual_rate = g.root_growth_rate * gf * auxin_boost;
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_root_growth;
    offset += growth_dir * actual_rate;
}

void RootApicalNode::spawn_internode(Plant& plant, const Genome& g) {

    // Create new interior root node and insert it between us and our parent
    Node* internode = plant.create_node(NodeType::ROOT, offset, radius);
    internode->rest_offset = internode->offset;  // remember stress-free direction
    parent->replace_child(this, internode);
    internode->position = internode->parent->position + internode->offset;
    offset = growth_dir * g.tip_offset;
    internode->add_child(this);
    position = internode->position + offset;

    // Lateral branching: compute offset and spawn axillary bud
    if (!plant.root_meristems_at_cap()) {
        glm::vec3 branch_dir_val = branch_direction(growth_dir, g.root_branch_angle, id);
        glm::vec3 ax_radial = branch_dir_val - growth_dir * glm::dot(branch_dir_val, growth_dir);
        float ax_rl = glm::length(ax_radial);
        if (ax_rl > 1e-4f) ax_radial /= ax_rl;
        glm::vec3 lateral_offset = ax_radial * internode->radius + branch_dir_val * g.tip_offset;
        spawn_axillary(plant, internode, g, lateral_offset);
    }

    ticks_since_last_node = 0;
    growth_dir = glm::vec3(0.0f); // re-roll direction for next internode
}

void RootApicalNode::spawn_axillary(Plant& plant, Node* internode, const Genome& g, const glm::vec3& lateral_offset) {
    Node* bud = plant.create_node(NodeType::ROOT_APICAL, lateral_offset, g.root_initial_radius * 0.5f);
    bud->as_root_apical()->active = false;
    internode->add_child(bud);
    bud->position = internode->position + bud->offset;
}

float RootApicalNode::maintenance_cost(const WorldParams& world) const {
    return world.sugar_maintenance_meristem;
}

bool RootApicalNode::can_activate(const Genome& g, const WorldParams& world) const {
    // Cytokinin from producing leaves signals "the plant can support new roots"
    float local_cyt = parent ? parent->chemical(ChemicalID::Cytokinin) : chemical(ChemicalID::Cytokinin);
    if (local_cyt < g.cytokinin_threshold) return false;

    if (chemical(ChemicalID::Sugar) < world.sugar_cost_activation) return false;

    return true;
}

void RootApicalNode::activate(const Genome& g, const WorldParams& world) {
    active = true;
    chemical(ChemicalID::Sugar) -= world.sugar_cost_activation;
    radius = g.root_initial_radius;
}

} // namespace botany
