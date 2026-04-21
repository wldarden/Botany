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

    // Record presence of primary root meristem for this tick's tracker.
    // Piggyback on the DFS walk (no extra traversal).  If the tracker is
    // still -1 at end of tick_tree, the primary died and promotion runs.
    if (is_primary) plant.primary_ra_id_this_tick = static_cast<int32_t>(id);

    absorb_water(g, world);

    // Quiescence: active meristem reverts to dormant after sustained sugar
    // starvation rather than dying.  Matches real root-tip biology — tips
    // survive weeks of stress by going quiescent, reactivate later via
    // normal can_activate() when sugar returns.  Runs before the active
    // check so the newly-dormant RA is handled by the dormant branch below.
    // EXCEPT: primary meristems never quiesce (same rule as ApicalNode) —
    // they sustain themselves at metabolic_factor floor or die normally.
    if (active && !is_primary && starvation_ticks >= static_cast<uint32_t>(g.quiescence_threshold)) {
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
    // cytokinin production from local metabolic state (sugar + water).  Both
    // productions use Michaelis-Menten metabolic gating.
    // Metabolic gating: sugar + water each contribute an MM factor; floor 0.1
    // matches SA convention (auxin conjugate pools buffer short-term stress).
    float mf_auxin = metabolic_factor(
        local().chemical(ChemicalID::Sugar), g.auxin_sugar_half_saturation, 0.1f,
        local().chemical(ChemicalID::Water), g.auxin_water_half_saturation, 0.1f);
    float produced_auxin = g.root_tip_auxin_production_rate * mf_auxin;
    local().chemical(ChemicalID::Auxin) += produced_auxin;
    tick_auxin_produced += produced_auxin;

    // Cytokinin: floor 0.05 (smaller than auxin's 0.1) — CK has less
    // conjugate-pool buffering than auxin, so its response to sugar is sharper.
    // CK production is driven purely by local metabolic state (sugar + water).
    // It is NOT gated on local auxin — in real root tips, CK synthesis reflects
    // "this tip is metabolically healthy," independent of the distant shoot.
    // (Prior versions gated on local_auxin, which broke for tall plants where
    // shoot-derived auxin never reaches the RA; see hormone-biology spec.)
    float mf_cyto = metabolic_factor(
        local().chemical(ChemicalID::Sugar), g.cytokinin_sugar_half_saturation, 0.05f,
        local().chemical(ChemicalID::Water), g.cytokinin_water_half_saturation, 0.05f);
    float cyto_produced = g.root_cytokinin_production_rate * mf_cyto;
    local().chemical(ChemicalID::Cytokinin) += cyto_produced;
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
    float fill_fraction = (cap > 1e-6f) ? local().chemical(ChemicalID::Water) / cap : 1.0f;
    float gradient = std::max(0.0f, world.soil_moisture - fill_fraction);
    float absorbed = g.water_absorption_rate * surface_area * gradient;
    float water_before = local().chemical(ChemicalID::Water);
    local().chemical(ChemicalID::Water) = std::min(water_before + absorbed, cap);
    tick_chem_produced[static_cast<size_t>(ChemicalID::Water)] += local().chemical(ChemicalID::Water) - water_before;
}

void RootApicalNode::check_spawn(Plant& plant, const Genome& g) {
    // dist_from_parent guard mirrors ApicalNode::check_spawn.  Without it,
    // a starved root apical that stops elongating still spawns one internode
    // per plastochron tick — and each such spawn creates a zero-length root
    // node (because growth_dir was reset to 0 by the previous spawn and
    // never re-rolled due to elongate()'s early-return).  These zero-length
    // nodes have zero xylem/phloem capacity and permanently trap any
    // chemical injected into them (observed: 37 mg of cytokinin stranded
    // in a chain of zero-cap roots at y≈-0.76 after 484 ticks).
    const float dist_from_parent = glm::length(offset);
    if (parent
        && ticks_since_last_node >= g.root_plastochron
        && dist_from_parent >= g.tip_offset
        && starvation_ticks == 0) {
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
    // Root elongation is gated by local cytokinin (not auxin).  CK at the RA
    // represents "this tip is metabolically healthy" — it's produced from the
    // RA's own sugar+water (see update_tissue).  Using CK as the modulator
    // keeps the gate scale-free: a 10m tree's deep RAs still produce their
    // own CK and gate their own growth, with no dependency on distant-shoot
    // signals that don't scale (see hormone-biology spec).  Sugar remains the
    // actual rate limiter — CK only permits growth if sugar is also available.
    float gf = growth_fraction(local().chemical(ChemicalID::Sugar), max_cost,
                               local().chemical(ChemicalID::Cytokinin), g.root_ck_growth_floor);
    if (gf < 1e-6f) return;
    float water_gf = turgor_fraction(local().chemical(ChemicalID::Water), water_cap(*this, g));
    if (water_gf < 1e-6f) return;
    gf *= water_gf;

    if (glm::length(growth_dir) < 1e-4f) roll_direction(g);

    float actual_rate = g.root_growth_rate * gf;
    local().chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_root_growth;
    tick_chem_consumed[static_cast<size_t>(ChemicalID::Sugar)] += actual_rate * world.sugar_cost_root_growth;
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
    bud->ever_active = false;  // dormant lateral — renderer skips until activated
    internode->add_child(bud);
    bud->position = internode->position + bud->offset;
}

float RootApicalNode::maintenance_cost(const WorldParams& world) const {
    return active ? world.sugar_maintenance_meristem : 0.0f;
}

bool RootApicalNode::can_activate(const Genome& g, const WorldParams& world) const {
    // Auxin from the shoot promotes root activation — "the shoot wants more roots".
    // Auxin reaches dormant RAs via local diffusion (short-range) and via PIN
    // cascade through root_forwarded, both of which populate local() for a
    // dormant bud.  So local auxin is the right sensor here.
    if (local().chemical(ChemicalID::Auxin) < g.root_auxin_activation_threshold) return false;

    // Cytokinin (produced by already-active root tips) inhibits activation of
    // additional roots, preventing runaway branching.  Sense both the own
    // local pool AND the parent root's xylem pool — take the max.  A dormant
    // bud never runs vascular extract, so its own local() CK is never
    // populated by the xylem stream; reading upstream xylem lets the bud
    // perceive passing sap.  Without this, the inhibition check was
    // effectively dead code and lateral RAs over-proliferated in CK-rich
    // regions of the root system.  Honoring own local too keeps fixture
    // tests that deposit CK directly at the bud behaving correctly.
    const auto* xyl = nearest_xylem_upstream();
    float ambient_ck = std::max(
        local().chemical(ChemicalID::Cytokinin),
        xyl ? xyl->chemical(ChemicalID::Cytokinin) : 0.0f);
    if (ambient_ck > g.root_cytokinin_inhibition_threshold) return false;

    if (local().chemical(ChemicalID::Sugar) < world.sugar_cost_activation) return false;

    return true;
}

void RootApicalNode::activate(const Genome& g, const WorldParams& world) {
    ever_active = true;
    active = true;
    local().chemical(ChemicalID::Sugar) -= world.sugar_cost_activation;
    tick_chem_consumed[static_cast<size_t>(ChemicalID::Sugar)] += world.sugar_cost_activation;
    radius = g.root_initial_radius;
}

} // namespace botany
