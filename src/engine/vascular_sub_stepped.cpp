#include "engine/vascular_sub_stepped.h"
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/node/node.h"
#include "engine/sugar.h"
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"
#include <glm/geometric.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace botany {

namespace {

// Flat view of the plant for efficient sub-step iteration.  Populated once
// at the start of vascular_sub_stepped(), then reused across all N iterations.
struct FlatNodes {
    std::vector<Node*> all;             // every node, DFS order
    std::vector<Node*> conduits;        // only nodes with phloem/xylem
    std::vector<Node*> leaves;          // LeafNode only
    std::vector<Node*> roots;           // RootNode only
    std::vector<Node*> apicals;         // ApicalNode only (includes SA)
    std::vector<Node*> root_apicals;    // RootApicalNode only
};

void collect_dfs(Node* n, FlatNodes& flat) {
    if (!n) return;
    flat.all.push_back(n);
    if (n->phloem() || n->xylem()) flat.conduits.push_back(n);
    switch (n->type) {
        case NodeType::LEAF:        flat.leaves.push_back(n); break;
        case NodeType::ROOT:        flat.roots.push_back(n); break;
        case NodeType::APICAL:      flat.apicals.push_back(n); break;
        case NodeType::ROOT_APICAL: flat.root_apicals.push_back(n); break;
        default: break;
    }
    for (Node* child : n->children) collect_dfs(child, flat);
}

FlatNodes flatten(Plant& plant) {
    FlatNodes flat;
    collect_dfs(plant.seed_mut(), flat);
    return flat;
}

} // anonymous namespace

void vascular_sub_stepped(Plant& /*plant*/, const Genome& /*g*/, const WorldParams& /*world*/) {
    // Empty stub — per-sub-step logic is added in Tasks 10-18.
}

VascularBudget compute_budget(Node& n, const Genome& g, const WorldParams& /*world*/) {
    VascularBudget b;
    const float sugar = n.local().chemical(ChemicalID::Sugar);
    const float water = n.local().chemical(ChemicalID::Water);
    const float cap_s = sugar_cap(n, g);
    const float cap_w = water_cap(n, g);

    switch (n.type) {
        case NodeType::LEAF: {
            const float reserve = g.leaf_reserve_fraction_sugar * cap_s;
            b.sugar_supply = std::max(0.0f, sugar - reserve);
            const float turgor_target = g.leaf_turgor_target_fraction * cap_w;
            b.water_demand = std::max(0.0f, turgor_target - water);
            break;
        }
        case NodeType::APICAL: {
            if (n.as_apical() && n.as_apical()->active) {
                const float target = g.meristem_sink_target_fraction * cap_s;
                b.sugar_demand = std::max(0.0f, target - sugar);
                const float turgor_target = g.leaf_turgor_target_fraction * cap_w;
                b.water_demand = std::max(0.0f, turgor_target - water);
            }
            break;
        }
        case NodeType::ROOT_APICAL: {
            if (n.as_root_apical() && n.as_root_apical()->active) {
                const float target = g.meristem_sink_target_fraction * cap_s;
                b.sugar_demand = std::max(0.0f, target - sugar);
                const float turgor_target = g.leaf_turgor_target_fraction * cap_w;
                b.water_demand = std::max(0.0f, turgor_target - water);
            }
            // Root apicals produce cytokinin.  Treat all local cytokinin as
            // surplus to inject into parent xylem.
            b.cytokinin_supply = n.local().chemical(ChemicalID::Cytokinin);
            break;
        }
        // STEM and ROOT get sugar via radial flow from their own phloem.
        // No direct inject/extract budget for them.
        default: break;
    }
    return b;
}

void inject_step(Node& source, const VascularBudget& b, uint32_t N, const Genome& /*g*/) {
    if (N == 0) return;

    // Sugar injection: leaves push into parent phloem (active pump).
    if (b.sugar_supply > 0.0f) {
        if (auto* target_pool = source.nearest_phloem_upstream()) {
            const float slice  = b.sugar_supply / static_cast<float>(N);
            const float actual = std::min(slice, source.local().chemical(ChemicalID::Sugar));
            source.local().chemical(ChemicalID::Sugar) -= actual;
            target_pool->chemical(ChemicalID::Sugar)   += actual;
        }
    }

    // Cytokinin injection: root apicals push into parent xylem.
    if (b.cytokinin_supply > 0.0f) {
        if (auto* target_pool = source.nearest_xylem_upstream()) {
            const float slice  = b.cytokinin_supply / static_cast<float>(N);
            const float actual = std::min(slice, source.local().chemical(ChemicalID::Cytokinin));
            source.local().chemical(ChemicalID::Cytokinin) -= actual;
            target_pool->chemical(ChemicalID::Cytokinin)   += actual;
        }
    }
}

void extract_step(Node& sink, const VascularBudget& b, uint32_t N, const Genome& /*g*/) {
    if (N == 0) return;

    if (b.sugar_demand > 0.0f) {
        if (auto* source_pool = sink.nearest_phloem_upstream()) {
            const float slice  = b.sugar_demand / static_cast<float>(N);
            const float actual = std::min(slice, source_pool->chemical(ChemicalID::Sugar));
            source_pool->chemical(ChemicalID::Sugar) -= actual;
            sink.local().chemical(ChemicalID::Sugar) += actual;
        }
    }

    if (b.water_demand > 0.0f) {
        if (auto* source_pool = sink.nearest_xylem_upstream()) {
            const float slice  = b.water_demand / static_cast<float>(N);
            const float actual = std::min(slice, source_pool->chemical(ChemicalID::Water));
            source_pool->chemical(ChemicalID::Water) -= actual;
            sink.local().chemical(ChemicalID::Water) += actual;
            // Cytokinin rides along passively — move cytokinin proportional
            // to water drawn from this pool.
            const float water_after = source_pool->chemical(ChemicalID::Water);
            if (water_after + actual > 1e-8f) {
                const float cyto_ratio = actual / (water_after + actual);
                const float cyto_move  = source_pool->chemical(ChemicalID::Cytokinin) * cyto_ratio;
                source_pool->chemical(ChemicalID::Cytokinin) -= cyto_move;
                sink.local().chemical(ChemicalID::Cytokinin) += cyto_move;
            }
        }
    }
}

float radial_permeability_sugar(float radius, const Genome& g) {
    const float base   = g.base_radial_permeability_sugar;
    const float floor  = g.radial_floor_fraction_sugar;
    const float r_half = g.radial_half_radius_sugar;
    const float ratio  = radius / r_half;
    return base * (floor + (1.0f - floor) / (1.0f + ratio * ratio));
}

float radial_permeability_water(float radius, const Genome& g) {
    const float base   = g.base_radial_permeability_water;
    const float floor  = g.radial_floor_fraction_water;
    const float r_half = g.radial_half_radius_water;
    const float ratio  = radius / r_half;
    return base * (floor + (1.0f - floor) / (1.0f + ratio * ratio));
}

void radial_flow_step(Node& n, uint32_t N, const Genome& g) {
    if (N == 0) return;
    const float inv_N = 1.0f / static_cast<float>(N);

    // --- Sugar: phloem <-> local ---
    if (auto* phl = n.phloem()) {
        const float perm    = radial_permeability_sugar(n.radius, g);
        const float cap_phl = phloem_capacity(n, g);
        const float cap_loc = sugar_cap(n, g);
        if (cap_phl > 1e-8f && cap_loc > 1e-8f) {
            const float conc_phl = phl->chemical(ChemicalID::Sugar) / cap_phl;
            const float conc_loc = n.local().chemical(ChemicalID::Sugar) / cap_loc;
            // Flow toward equal concentration; positive = phloem -> local.
            const float dconc = conc_phl - conc_loc;
            float flow = perm * dconc * inv_N;
            // Cap the flow to avoid overshoot in one step.
            const float max_equalize_volume = std::min(cap_phl, cap_loc) * 0.5f;
            flow = std::clamp(flow, -max_equalize_volume, max_equalize_volume);
            // Don't draw more than what's present in either side.
            if (flow > 0.0f) flow = std::min(flow, phl->chemical(ChemicalID::Sugar));
            else             flow = -std::min(-flow, n.local().chemical(ChemicalID::Sugar));

            phl->chemical(ChemicalID::Sugar)      -= flow;
            n.local().chemical(ChemicalID::Sugar) += flow;
        }
    }

    // --- Water + cytokinin: xylem <-> local ---
    if (auto* xyl = n.xylem()) {
        const float perm    = radial_permeability_water(n.radius, g);
        const float cap_xyl = xylem_capacity(n, g);
        const float cap_loc = water_cap(n, g);
        if (cap_xyl > 1e-8f && cap_loc > 1e-8f) {
            const float conc_xyl = xyl->chemical(ChemicalID::Water) / cap_xyl;
            const float conc_loc = n.local().chemical(ChemicalID::Water) / cap_loc;
            const float dconc = conc_xyl - conc_loc;
            float flow = perm * dconc * inv_N;
            const float max_equalize_volume = std::min(cap_xyl, cap_loc) * 0.5f;
            flow = std::clamp(flow, -max_equalize_volume, max_equalize_volume);
            if (flow > 0.0f) flow = std::min(flow, xyl->chemical(ChemicalID::Water));
            else             flow = -std::min(-flow, n.local().chemical(ChemicalID::Water));

            xyl->chemical(ChemicalID::Water)      -= flow;
            n.local().chemical(ChemicalID::Water) += flow;

            // Cytokinin rides with water proportionally.
            const float water_source_before = (flow > 0.0f)
                ? xyl->chemical(ChemicalID::Water) + flow
                : n.local().chemical(ChemicalID::Water) - flow;
            if (water_source_before > 1e-8f && std::abs(flow) > 1e-8f) {
                const float cyto_ratio = std::abs(flow) / water_source_before;
                if (flow > 0.0f) {
                    // Water moved xylem -> local, cytokinin follows.
                    const float cyto_move = xyl->chemical(ChemicalID::Cytokinin) * cyto_ratio;
                    xyl->chemical(ChemicalID::Cytokinin)      -= cyto_move;
                    n.local().chemical(ChemicalID::Cytokinin) += cyto_move;
                } else {
                    // Water moved local -> xylem, cytokinin follows.
                    const float cyto_move = n.local().chemical(ChemicalID::Cytokinin) * cyto_ratio;
                    n.local().chemical(ChemicalID::Cytokinin) -= cyto_move;
                    xyl->chemical(ChemicalID::Cytokinin)      += cyto_move;
                }
            }
        }
    }
}

void jacobi_step(Node& parent, Node& child, const Genome& g) {
    // Phloem sugar.
    if (auto* pp = parent.phloem()) {
        if (auto* cp = child.phloem()) {
            const float cap_p = phloem_capacity(parent, g);
            const float cap_c = phloem_capacity(child,  g);
            if (cap_p > 1e-8f && cap_c > 1e-8f) {
                const float pressure_p = pp->chemical(ChemicalID::Sugar) / cap_p;
                const float pressure_c = cp->chemical(ChemicalID::Sugar) / cap_c;
                const float base_cond  = std::min(cap_p, cap_c);
                const float bias       = parent.get_bias_multiplier(&child, g);
                const float conductance = base_cond * bias;
                float flow = conductance * (pressure_p - pressure_c);
                const float max_move = 0.5f * std::min(cap_p, cap_c);
                flow = std::clamp(flow, -max_move, max_move);
                if (flow > 0.0f) flow = std::min(flow, pp->chemical(ChemicalID::Sugar));
                else             flow = -std::min(-flow, cp->chemical(ChemicalID::Sugar));

                pp->chemical(ChemicalID::Sugar) -= flow;
                cp->chemical(ChemicalID::Sugar) += flow;
            }
        }
    }

    // Xylem water + cytokinin.
    if (auto* pp = parent.xylem()) {
        if (auto* cp = child.xylem()) {
            const float cap_p = xylem_capacity(parent, g);
            const float cap_c = xylem_capacity(child,  g);
            if (cap_p > 1e-8f && cap_c > 1e-8f) {
                const float pressure_p = pp->chemical(ChemicalID::Water) / cap_p;
                const float pressure_c = cp->chemical(ChemicalID::Water) / cap_c;
                const float base_cond  = std::min(cap_p, cap_c);
                const float bias       = parent.get_bias_multiplier(&child, g);
                const float conductance = base_cond * bias;
                float flow = conductance * (pressure_p - pressure_c);
                const float max_move = 0.5f * std::min(cap_p, cap_c);
                flow = std::clamp(flow, -max_move, max_move);
                if (flow > 0.0f) flow = std::min(flow, pp->chemical(ChemicalID::Water));
                else             flow = -std::min(-flow, cp->chemical(ChemicalID::Water));

                // Water move.
                pp->chemical(ChemicalID::Water) -= flow;
                cp->chemical(ChemicalID::Water) += flow;

                // Cytokinin rides with water.
                if (std::abs(flow) > 1e-8f) {
                    TransportPool* src = (flow > 0.0f) ? pp : cp;
                    TransportPool* dst = (flow > 0.0f) ? cp : pp;
                    // Water in source BEFORE the transfer.
                    const float water_source_before = (flow > 0.0f)
                        ? pp->chemical(ChemicalID::Water) + flow
                        : cp->chemical(ChemicalID::Water) - flow;
                    if (water_source_before > 1e-8f) {
                        const float cyto_ratio = std::abs(flow) / water_source_before;
                        const float cyto_move  = src->chemical(ChemicalID::Cytokinin) * cyto_ratio;
                        src->chemical(ChemicalID::Cytokinin) -= cyto_move;
                        dst->chemical(ChemicalID::Cytokinin) += cyto_move;
                    }
                }
            }
        }
    }
}

// Returns the length of the internode (distance from parent to this node).
// For the seed (no parent), returns 1.0 so the seed has a meaningful capacity.
static float node_length(const Node& n) {
    if (!n.parent) return 1.0f;
    return glm::length(n.offset);
}

float phloem_capacity(const Node& n, const Genome& g) {
    if (!n.phloem()) return 0.0f;
    const float L = node_length(n);
    return 3.14159265358979f * n.radius * n.radius * L * g.phloem_fraction;
}

float xylem_capacity(const Node& n, const Genome& g) {
    if (!n.xylem()) return 0.0f;
    const float L = node_length(n);
    return 3.14159265358979f * n.radius * n.radius * L * g.xylem_fraction;
}

} // namespace botany
