#include "engine/node/tissues/apical.h"
#include "engine/node/meristems/helpers.h"
#include "engine/plant.h"
#include "engine/node/tissues/leaf.h"
#include "engine/sugar.h"
#include "engine/world_params.h"
#include <unordered_set>

namespace botany {

using namespace meristem_helpers;

ApicalNode::ApicalNode(uint32_t id, glm::vec3 position, float radius)
    : Node(id, NodeType::APICAL, position, radius)
{}

void ApicalNode::update_tissue(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    // Record presence of primary shoot meristem for this tick's tracker.
    // Piggyback on the DFS walk — no extra traversal.  If this tracker is
    // still -1 at end of tick_tree, Plant::promote_primary_meristems runs
    // and picks a new primary from the surviving laterals.
    if (is_primary) plant.primary_sa_id_this_tick = static_cast<int32_t>(id);

    // Quiescence: active meristem reverts to dormant after sustained sugar
    // starvation rather than dying.  Symmetric with RA quiescence.
    // EXCEPT: primary meristems never quiesce — they sustain themselves at
    // metabolic_factor floor production (10% auxin, 5% CK at full starvation)
    // or die via normal starvation.  Only lateral (axillary) buds can quiesce.
    // Prevents the hormone deadlock where all meristems dormant simultaneously
    // = no hormone production anywhere = nothing can reactivate.
    if (active && !is_primary && starvation_ticks >= static_cast<uint32_t>(g.quiescence_threshold)) {
        active = false;
        starvation_ticks = 0;
    }

    if (!active) {
        // Dormant buds accumulate nothing — zero maintenance cost, and we
        // keep starvation_ticks at 0 so check_starvation can't drive them
        // to death while dormant.  They can always reactivate via can_activate().
        starvation_ticks = 0;
        if (can_activate(g, world)) activate(g, world);
        return;
    }

    produce_auxin(plant, g, world);
    photosynthesize(plant, g, world);
    ticks_since_last_node++;
    elongate(plant, g, world);
    check_spawn(plant, g);
}

void ApicalNode::produce_auxin(Plant& /*plant*/, const Genome& g, const WorldParams& world) {
    float max_cost = g.growth_rate * world.sugar_cost_meristem_growth;
    float growth_gf = meristem_helpers::growth_fraction(
        local().chemical(ChemicalID::Sugar), max_cost,
        local().chemical(ChemicalID::Cytokinin), g.cytokinin_growth_threshold);

    float base = g.apical_auxin_baseline;
    float local_light = estimate_local_light();
    float light_factor = 1.0f + g.auxin_shade_boost * (1.0f - local_light);
    // Metabolic gating: sugar and water both required for full auxin synthesis.
    // Floor 0.1 matches conjugate-pool buffering — zero metabolism still yields ~1%.
    float mf = metabolic_factor(
        local().chemical(ChemicalID::Sugar), g.auxin_sugar_half_saturation, 0.1f,
        local().chemical(ChemicalID::Water), g.auxin_water_half_saturation, 0.1f);
    float modulated_baseline = base * light_factor * mf;
    float produced = modulated_baseline * (1.0f + g.apical_growth_auxin_multiplier * growth_gf);
    local().chemical(ChemicalID::Auxin) += produced;
    tick_auxin_produced += produced;
}

void ApicalNode::photosynthesize(Plant& plant, const Genome& g, const WorldParams& world) {
    // Green shoot tips photosynthesize at low efficiency — roughly enough
    // to self-sustain in full light, but not enough to be net exporters.
    float light = estimate_local_light() * world.light_level;
    float production = light * world.sugar_maintenance_meristem * world.sugar_meristem_photosynthesis;
    if (production > 0.0f) {
        local().chemical(ChemicalID::Sugar) += production;
        plant.add_sugar_produced(production);

    }
}

void ApicalNode::check_spawn(Plant& plant, const Genome& g) {
    float dist_from_parent = glm::length(offset);
    // Minimum internode length 2 × tip_offset (= 2mm at default genome).
    // After each spawn the SAM resets to ~1mm (tip_offset), so requiring
    // 2mm before the next spawn means every spawned internode has
    // geometric length ≥ 1mm of actual accumulated growth plus the
    // 1mm tip, giving conduit pipes a minimum physical size that
    // matters for transport throughput.  Was g.tip_offset (1mm), which
    // let the SAM spawn internodes with barely any grown length.
    if (parent && ticks_since_last_node >= g.shoot_plastochron
        && dist_from_parent >= 2.0f * g.tip_offset
        && starvation_ticks == 0) {
        spawn_internode(plant, g);
    }
}

float ApicalNode::estimate_local_light() const {
    // Average light_exposure of leaves on nearby internodes (parent chain).
    // Gives the meristem a local sense of canopy shade.
    float total_light = 0.0f;
    int leaf_count = 0;

    const Node* stem = parent;
    for (int i = 0; i < 3 && stem; ++i, stem = stem->parent) {
        for (const Node* child : stem->children) {
            if (child->type == NodeType::LEAF) {
                total_light += child->as_leaf()->light_exposure;
                leaf_count++;
            }
        }
    }

    if (leaf_count == 0) return 1.0f;  // no leaves yet — assume full sun
    return total_light / static_cast<float>(leaf_count);
}


void ApicalNode::roll_direction(const Genome& g, const WorldParams& world) {
    growth_dir = perturb(growth_direction(*this), g.growth_noise);

    // Plagiotropism: pull toward set-point angle (vertical for trunk, branch angle for branches).
    // Stress hormone boosts the correction.
    {
        float base_pull = g.meristem_gravitropism_rate;
        float stress_pull = local().chemical(ChemicalID::Stress) * g.stress_gravitropism_boost;
        float blend = std::min(base_pull + stress_pull, 0.5f);
        if (blend > 1e-6f) {
            growth_dir = glm::normalize(glm::mix(growth_dir, set_point_dir, blend));
        }
    }

    // Phototropism: in shade, grow toward the light source.
    {
        float light = estimate_local_light() * world.light_level;
        float shade = 1.0f - light;
        float photo_blend = std::min(shade * g.meristem_phototropism_rate, 0.3f);
        if (photo_blend > 1e-6f) {
            growth_dir = glm::normalize(glm::mix(growth_dir, world.light_direction, photo_blend));
        }
    }
}

void ApicalNode::elongate(Plant& plant, const Genome& g, const WorldParams& world) {
    float max_cost = g.growth_rate * world.sugar_cost_meristem_growth;
    float gf = growth_fraction(local().chemical(ChemicalID::Sugar), max_cost,
                               local().chemical(ChemicalID::Cytokinin), g.cytokinin_growth_threshold);
    if (gf < 1e-6f) return;
    float water_gf = turgor_fraction(local().chemical(ChemicalID::Water), water_cap(*this, g));
    if (water_gf < 1e-6f) return;
    gf *= water_gf;

    // Roll a fresh direction on internode spawn (noise + plagiotropism + phototropism)
    if (glm::length(growth_dir) < 1e-4f) roll_direction(g, world);

    // Deflect away from nearby non-ancestor stems to avoid growing through them.
    // Only meristems do this check — there are very few active at any time.
    {
        // Build ancestor set (don't deflect from our own trunk)
        std::unordered_set<const Node*> ancestors;
        for (const Node* p = parent; p; p = p->parent) ancestors.insert(p);

        glm::vec3 push(0.0f);
        plant.for_each_node([&](const Node& other) {
            if (other.type != NodeType::STEM) return;
            if (ancestors.count(&other)) return;

            glm::vec3 diff = position - other.position;
            float dist = glm::length(diff);
            float min_dist = other.radius + radius + 0.05f;  // 5mm margin
            if (dist < min_dist && dist > 1e-6f) {
                // Push proportional to overlap
                push += (diff / dist) * (min_dist - dist);
            }
        });

        float push_len = glm::length(push);
        if (push_len > 1e-6f) {
            // Blend push direction into growth_dir
            float blend = std::min(push_len * 2.0f, 0.5f);
            growth_dir = glm::normalize(glm::mix(growth_dir, push / push_len, blend));
        }
    }

    // Ground collision: don't grow downward when at or below ground level
    if (position.y <= 0.0f && growth_dir.y < 0.0f) {
        growth_dir.y = 0.0f;
        float len = glm::length(growth_dir);
        if (len > 1e-6f) {
            growth_dir /= len;  // re-normalize, keeps horizontal component
        } else {
            growth_dir = glm::vec3(0.0f, 1.0f, 0.0f);  // straight up if no horizontal component
        }
    }

    float auxin_boost = meristem_helpers::auxin_growth_factor(
        local().chemical(ChemicalID::Auxin), g.apical_auxin_max_boost, g.apical_auxin_half_saturation);
    float actual_rate = g.growth_rate * gf * auxin_boost;
    local().chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_meristem_growth;
    offset += growth_dir * actual_rate;
}

void ApicalNode::spawn_internode(Plant& plant, const Genome& g) {

    // Create new interior stem node and insert it between us and our parent
    Node* internode = plant.create_node(NodeType::STEM, offset, radius);
    internode->rest_offset = internode->offset;  // remember stress-free direction
    parent->replace_child(this, internode);
    internode->position = internode->parent->position + internode->offset;

    // After spawning, reset the SAM to a 1mm tip ahead of the new internode.
    // Safety floor: if growth_dir is zero (happens if elongate returned early
    // due to low sugar/cytokinin and the previous spawn zeroed growth_dir),
    // fall back to "up" so the SAM is never at length 0.  Guarantees SAM
    // offset >= tip_offset = 1mm at all times.
    glm::vec3 dir = growth_dir;
    if (glm::length(dir) < 1e-4f) {
        dir = glm::vec3(0.0f, 1.0f, 0.0f);
    }
    offset = dir * g.tip_offset;
    internode->add_child(this);
    position = internode->position + offset;

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

void ApicalNode::spawn_axillary(Plant& plant, Node* internode, const Genome& g, const glm::vec3& lateral_offset) {
    Node* bud = plant.create_node(NodeType::APICAL, lateral_offset, g.initial_radius * 0.5f);
    bud->as_apical()->active = false;
    bud->ever_active = false;  // dormant lateral — renderer skips until activated
    internode->add_child(bud);
    bud->position = internode->position + bud->offset;
}

void ApicalNode::spawn_leaf(Plant& plant, Node* internode, const Genome& g, const glm::vec3& lateral_offset) {
    // Start leaf close to stem — petiole grows with leaf size
    Node* leaf = plant.create_node(NodeType::LEAF, lateral_offset, 0.0f);
    leaf->as_leaf()->leaf_size = g.leaf_bud_size;
    float len = glm::length(lateral_offset);
    if (len > 1e-4f) {
        leaf->as_leaf()->facing = lateral_offset / len;
    }
    // Meristem gives the new leaf some sugar to bootstrap it
    float gift = std::min(local().chemical(ChemicalID::Sugar) * 0.1f, 0.5f);
    local().chemical(ChemicalID::Sugar) -= gift;
    leaf->local().chemical(ChemicalID::Sugar) = gift;
    internode->add_child(leaf);
    leaf->position = internode->position + leaf->offset;
}

float ApicalNode::maintenance_cost(const WorldParams& world) const {
    return active ? world.sugar_maintenance_meristem : 0.0f;
}

bool ApicalNode::can_activate(const Genome& g, const WorldParams& world) const {
    // Low auxin removes inhibition (apical dominance weakened)
    float stem_auxin = parent ? parent->local().chemical(ChemicalID::Auxin) : local().chemical(ChemicalID::Auxin);
    if (stem_auxin >= g.auxin_threshold) return false;

    // Cytokinin from roots (xylem-delivered) signals "the root system is healthy".
    // Sense both the own local pool AND the parent stem's xylem pool — take
    // the max.  A dormant bud doesn't run vascular extract (no water demand),
    // so its own local() CK is never populated by the xylem stream; reading
    // the upstream xylem directly lets the bud perceive sap flowing past it,
    // analogous to how real dormant buds sense adjacent xylem without active
    // uptake.  We still honor local CK because any mechanism that deposits
    // CK directly at the bud (test fixtures, future radial flow) should also
    // be able to activate it.
    const auto* xyl = nearest_xylem_upstream();
    float ambient_ck = std::max(
        local().chemical(ChemicalID::Cytokinin),
        xyl ? xyl->chemical(ChemicalID::Cytokinin) : 0.0f);
    if (ambient_ck < g.cytokinin_threshold) return false;

    if (local().chemical(ChemicalID::Sugar) < world.sugar_cost_activation) return false;

    return true;
}

void ApicalNode::activate(const Genome& g, const WorldParams& world) {
    ever_active = true;
    active = true;
    local().chemical(ChemicalID::Sugar) -= world.sugar_cost_activation;
    radius = g.initial_radius;

    // Set growth direction from branch offset
    float olen = glm::length(offset);
    if (olen > 1e-4f) {
        glm::vec3 branch_dir = offset / olen;
        growth_dir = branch_dir;
        set_point_dir = branch_dir;
    }
}

} // namespace botany
