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

    // Quiescence: active meristem reverts to dormant after sustained sugar
    // starvation rather than dying.  Matches real root-tip biology — tips
    // survive weeks of stress by going quiescent, reactivate later via
    // normal can_activate() when sugar returns.  Runs before the active
    // check so the newly-dormant RA is handled by the dormant branch below.
    if (active && starvation_ticks >= static_cast<uint32_t>(g.quiescence_threshold)) {
        active = false;
        starvation_ticks = 0;  // dormant meristem pays no maintenance; clean slate for reactivation
    }

    if (!active) {
        // Dormant RAs do not synthesize auxin or cytokinin — matches shoot apical
        // behavior and real plant biology (hormone biosynthesis is an active process
        // gated on metabolism).  Dormant buds depend on PIN-transported auxin from
        // upstream to eventually cross the activation threshold, which is the
        // canalization-driven branching pattern real plants exhibit: laterals
        // activate along well-canalized paths, not uniformly everywhere.
        //
        // Quiescent meristems cannot starve to death — they have no active metabolism.
        // Reset starvation_ticks each tick so check_starvation() (which runs after
        // update_tissue) cannot accumulate to starvation_ticks_max.
        starvation_ticks = 0;
        if (can_activate(g, world)) activate(g, world);
        return;
    }

    // Active root tips: auxin self-maintenance (local PIN-recycling maximum) and
    // cytokinin production gated by local auxin ("more auxin → stronger cyto
    // signal").  Both productions match the shoot apical pattern.
    // Metabolic gating: sugar + water each contribute an MM factor; floor 0.1
    // matches SA convention (auxin conjugate pools buffer short-term stress).
    float mf_auxin = metabolic_factor(
        chemical(ChemicalID::Sugar), g.auxin_sugar_half_saturation, 0.1f,
        chemical(ChemicalID::Water), g.auxin_water_half_saturation, 0.1f);
    float produced_auxin = g.root_tip_auxin_production_rate * mf_auxin;
    chemical(ChemicalID::Auxin) += produced_auxin;
    tick_auxin_produced += produced_auxin;

    // Cytokinin: floor 0.05 (smaller than auxin's 0.1) — CK has less
    // conjugate-pool buffering than auxin, so its response to sugar is sharper.
    // This is the primary root-to-shoot feedback brake: low root sugar → low
    // CK → less CK delivered to shoot via xylem → fewer SA activations.
    float mf_cyto = metabolic_factor(
        chemical(ChemicalID::Sugar), g.cytokinin_sugar_half_saturation, 0.05f,
        chemical(ChemicalID::Water), g.cytokinin_water_half_saturation, 0.05f);
    float cyto_produced = g.root_cytokinin_production_rate * chemical(ChemicalID::Auxin) * mf_cyto;
    chemical(ChemicalID::Cytokinin) += cyto_produced;
    tick_cytokinin_produced += cyto_produced;

    ticks_since_last_node++;
    elongate(g, world);
    check_spawn(plant, g);
}

void RootApicalNode::absorb_water(const Genome& g, const WorldParams& world) {
    // Hemisphere surface area with root hair multiplier.
    // Real root tips have dense root hairs that increase effective
    // absorption area by 10-100x over bare geometry.
    constexpr float root_hair_multiplier = 20.0f;
    float surface_area = 2.0f * 3.14159f * radius * radius * root_hair_multiplier;
    float cap = water_cap(*this, g);
    float fill_fraction = (cap > 1e-6f) ? chemical(ChemicalID::Water) / cap : 1.0f;
    float gradient = std::max(0.0f, world.soil_moisture - fill_fraction);
    float absorbed = g.water_absorption_rate * surface_area * gradient;
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
    // Root elongation is sugar-gated only — real root tips maintain their own
    // auxin maximum via PIN recycling, so we don't gate on exogenous auxin.
    float gf = sugar_growth_fraction(chemical(ChemicalID::Sugar), max_cost);
    if (gf < 1e-6f) return;
    float water_gf = turgor_fraction(chemical(ChemicalID::Water), water_cap(*this, g));
    if (water_gf < 1e-6f) return;
    gf *= water_gf;

    if (glm::length(growth_dir) < 1e-4f) roll_direction(g);

    float actual_rate = g.root_growth_rate * gf;
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

void RootApicalNode::compute_growth_reserve(const Genome& g, const WorldParams& world) {
    sugar_reserved_for_growth = 0.0f;
    if (!active) return;

    // Mirror elongate(): max sugar cost at full growth rate.
    float max_cost = g.root_growth_rate * world.sugar_cost_root_growth;
    sugar_reserved_for_growth = std::min(max_cost, chemical(ChemicalID::Sugar));
}

float RootApicalNode::maintenance_cost(const WorldParams& world) const {
    return active ? world.sugar_maintenance_meristem : 0.0f;
}

bool RootApicalNode::can_activate(const Genome& g, const WorldParams& world) const {
    // Auxin from the shoot promotes root activation — "the shoot wants more roots".
    // Existing cytokinin (produced by already-active roots) inhibits activation
    // of additional roots, preventing runaway branching.
    if (chemical(ChemicalID::Auxin) < g.root_auxin_activation_threshold) return false;

    if (chemical(ChemicalID::Cytokinin) > g.root_cytokinin_inhibition_threshold) return false;

    if (chemical(ChemicalID::Sugar) < world.sugar_cost_activation) return false;

    return true;
}

void RootApicalNode::activate(const Genome& g, const WorldParams& world) {
    active = true;
    chemical(ChemicalID::Sugar) -= world.sugar_cost_activation;
    radius = g.root_initial_radius;
}

} // namespace botany
