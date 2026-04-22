#include "engine/vascular_sub_stepped.h"
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/node/node.h"
#include "engine/sugar.h"
#include "engine/perf_log.h"
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
    uint32_t max_conduit_depth = 0;     // number of conduit edges from seed to deepest conduit node
};

// DFS collection + max-conduit-chain-depth tracking.  `depth` = number of
// conduit edges traversed from the seed to n.  Conduit edge := both endpoints
// have phloem or xylem (i.e., both are in flat.conduits).  Only conduit edges
// lengthen the Jacobi pressure-propagation chain, so non-conduit descents
// (through leaves / meristem tips) don't count.
void collect_dfs(Node* n, FlatNodes& flat, uint32_t depth = 0) {
    if (!n) return;
    flat.all.push_back(n);
    const bool n_is_conduit = (n->phloem() || n->xylem());
    if (n_is_conduit) {
        flat.conduits.push_back(n);
        flat.max_conduit_depth = std::max(flat.max_conduit_depth, depth);
    }
    switch (n->type) {
        case NodeType::LEAF:        flat.leaves.push_back(n); break;
        case NodeType::ROOT:        flat.roots.push_back(n); break;
        case NodeType::APICAL:      flat.apicals.push_back(n); break;
        case NodeType::ROOT_APICAL: flat.root_apicals.push_back(n); break;
        default: break;
    }
    for (Node* child : n->children) {
        const bool child_is_conduit = (child->phloem() || child->xylem());
        const uint32_t next_depth = (n_is_conduit && child_is_conduit) ? depth + 1 : depth;
        collect_dfs(child, flat, next_depth);
    }
}

FlatNodes flatten(Plant& plant) {
    FlatNodes flat;
    collect_dfs(plant.seed_mut(), flat);
    return flat;
}

// --- Precomputed per-edge Jacobi coefficients --------------------------------
// All values derived from plant geometry (radius, length) and canalization bias,
// neither of which changes inside vascular_sub_stepped's N-sub-step loop.  Compute
// once per tick, reuse N times.
struct EdgeJacobiCoeffs {
    bool  has_phloem_pair   = false;
    float cap_p_phloem      = 0.0f;
    float cap_c_phloem      = 0.0f;
    float base_cond_phloem  = 0.0f;  // min(cap_p, cap_c) — used for instrumentation
    float conductance_phloem = 0.0f; // base_cond * bias — used for flow

    bool  has_xylem_pair    = false;
    float cap_p_xylem       = 0.0f;
    float cap_c_xylem       = 0.0f;
    float base_cond_xylem   = 0.0f;
    float conductance_xylem = 0.0f;
};

EdgeJacobiCoeffs compute_edge_jacobi_coeffs(Node& parent, Node& child, const Genome& g) {
    EdgeJacobiCoeffs c;
    if (parent.phloem() && child.phloem()) {
        c.cap_p_phloem = phloem_capacity(parent, g);
        c.cap_c_phloem = phloem_capacity(child,  g);
        if (c.cap_p_phloem > 1e-8f && c.cap_c_phloem > 1e-8f) {
            c.base_cond_phloem  = std::min(c.cap_p_phloem, c.cap_c_phloem);
            const float bias    = parent.get_bias_multiplier(&child, g);
            c.conductance_phloem = c.base_cond_phloem * bias;
            c.has_phloem_pair = true;
        }
    }
    if (parent.xylem() && child.xylem()) {
        c.cap_p_xylem = xylem_capacity(parent, g);
        c.cap_c_xylem = xylem_capacity(child,  g);
        if (c.cap_p_xylem > 1e-8f && c.cap_c_xylem > 1e-8f) {
            c.base_cond_xylem  = std::min(c.cap_p_xylem, c.cap_c_xylem);
            const float bias   = parent.get_bias_multiplier(&child, g);
            c.conductance_xylem = c.base_cond_xylem * bias;
            c.has_xylem_pair = true;
        }
    }
    return c;
}

// Inner Jacobi body — takes precomputed coeffs.  Behavior identical to the
// original jacobi_step(); see jacobi_step() body for comments on the physics.
void jacobi_step_impl(Node& parent, Node& child, const EdgeJacobiCoeffs& c) {
    if (c.has_phloem_pair) {
        auto* pp = parent.phloem();
        auto* cp = child.phloem();
        const float pressure_p = pp->chemical(ChemicalID::Sugar) / c.cap_p_phloem;
        const float pressure_c = cp->chemical(ChemicalID::Sugar) / c.cap_c_phloem;
        float flow = c.conductance_phloem * (pressure_p - pressure_c);
        if (flow > 0.0f) flow = std::min(flow, pp->chemical(ChemicalID::Sugar));
        else             flow = -std::min(-flow, cp->chemical(ChemicalID::Sugar));

        pp->chemical(ChemicalID::Sugar) -= flow;
        cp->chemical(ChemicalID::Sugar) += flow;

        const size_t SI = static_cast<size_t>(ChemicalID::Sugar);
        parent.tick_edge_flux[SI][&child] += flow;
        parent.tick_edge_cap[SI][&child]  += c.base_cond_phloem;
    }

    if (c.has_xylem_pair) {
        auto* pp = parent.xylem();
        auto* cp = child.xylem();
        const float pressure_p = pp->chemical(ChemicalID::Water) / c.cap_p_xylem;
        const float pressure_c = cp->chemical(ChemicalID::Water) / c.cap_c_xylem;
        float flow = c.conductance_xylem * (pressure_p - pressure_c);
        if (flow > 0.0f) flow = std::min(flow, pp->chemical(ChemicalID::Water));
        else             flow = -std::min(-flow, cp->chemical(ChemicalID::Water));

        pp->chemical(ChemicalID::Water) -= flow;
        cp->chemical(ChemicalID::Water) += flow;

        // Cytokinin has its own pressure gradient in xylem; see jacobi_step() for rationale.
        const float cyto_pressure_p = pp->chemical(ChemicalID::Cytokinin) / c.cap_p_xylem;
        const float cyto_pressure_c = cp->chemical(ChemicalID::Cytokinin) / c.cap_c_xylem;
        float cyto_flow = c.conductance_xylem * (cyto_pressure_p - cyto_pressure_c);
        if (cyto_flow > 0.0f) cyto_flow = std::min(cyto_flow, pp->chemical(ChemicalID::Cytokinin));
        else                  cyto_flow = -std::min(-cyto_flow, cp->chemical(ChemicalID::Cytokinin));

        pp->chemical(ChemicalID::Cytokinin) -= cyto_flow;
        cp->chemical(ChemicalID::Cytokinin) += cyto_flow;

        const size_t WI = static_cast<size_t>(ChemicalID::Water);
        const size_t CI = static_cast<size_t>(ChemicalID::Cytokinin);
        parent.tick_edge_flux[WI][&child] += flow;
        parent.tick_edge_cap[WI][&child]  += c.base_cond_xylem;
        parent.tick_edge_flux[CI][&child] += cyto_flow;
        parent.tick_edge_cap[CI][&child]  += c.base_cond_xylem;
    }
}

// --- Precomputed per-node radial coefficients --------------------------------
// Constant across the N-sub-step loop (depend only on node radius).
struct NodeRadialCoeffs {
    bool  has_phloem    = false;
    float perm_sugar    = 0.0f;
    float cap_phl       = 0.0f;
    float cap_loc_sugar = 0.0f;

    bool  has_xylem     = false;
    float perm_water    = 0.0f;
    float cap_xyl       = 0.0f;
    float cap_loc_water = 0.0f;
};

NodeRadialCoeffs compute_node_radial_coeffs(const Node& n, const Genome& g) {
    NodeRadialCoeffs c;
    if (n.phloem()) {
        c.cap_phl       = phloem_capacity(n, g);
        c.cap_loc_sugar = sugar_cap(n, g);
        if (c.cap_phl > 1e-8f && c.cap_loc_sugar > 1e-8f) {
            c.perm_sugar = radial_permeability_sugar(n.radius, g);
            c.has_phloem = true;
        }
    }
    if (n.xylem()) {
        c.cap_xyl       = xylem_capacity(n, g);
        c.cap_loc_water = water_cap(n, g);
        if (c.cap_xyl > 1e-8f && c.cap_loc_water > 1e-8f) {
            c.perm_water = radial_permeability_water(n.radius, g);
            c.has_xylem  = true;
        }
    }
    return c;
}

// Inner radial body — takes precomputed coeffs.  Behavior identical to the
// original radial_flow_step(); see that function for comments on the physics.
void radial_flow_step_impl(Node& n, uint32_t N, const NodeRadialCoeffs& c) {
    if (N == 0) return;
    const float inv_N = 1.0f / static_cast<float>(N);

    if (c.has_phloem) {
        auto* phl = n.phloem();
        const float conc_phl = phl->chemical(ChemicalID::Sugar) / c.cap_phl;
        const float conc_loc = n.local().chemical(ChemicalID::Sugar) / c.cap_loc_sugar;
        const float dconc = conc_phl - conc_loc;
        float flow = c.perm_sugar * dconc * inv_N;
        const float max_equalize_volume = std::min(c.cap_phl, c.cap_loc_sugar) * 0.5f;
        flow = std::clamp(flow, -max_equalize_volume, max_equalize_volume);
        if (flow > 0.0f) flow = std::min(flow, phl->chemical(ChemicalID::Sugar));
        else             flow = -std::min(-flow, n.local().chemical(ChemicalID::Sugar));

        phl->chemical(ChemicalID::Sugar)      -= flow;
        n.local().chemical(ChemicalID::Sugar) += flow;
    }

    if (c.has_xylem) {
        auto* xyl = n.xylem();
        const float conc_xyl = xyl->chemical(ChemicalID::Water) / c.cap_xyl;
        const float conc_loc = n.local().chemical(ChemicalID::Water) / c.cap_loc_water;
        const float dconc = conc_xyl - conc_loc;
        float flow = c.perm_water * dconc * inv_N;
        const float max_equalize_volume = std::min(c.cap_xyl, c.cap_loc_water) * 0.5f;
        flow = std::clamp(flow, -max_equalize_volume, max_equalize_volume);
        if (flow > 0.0f) flow = std::min(flow, xyl->chemical(ChemicalID::Water));
        else             flow = -std::min(-flow, n.local().chemical(ChemicalID::Water));

        xyl->chemical(ChemicalID::Water)      -= flow;
        n.local().chemical(ChemicalID::Water) += flow;
        // Cytokinin is intentionally NOT radially exchanged; see radial_flow_step() for why.
    }
}

} // anonymous namespace

void vascular_sub_stepped(Plant& plant, const Genome& g, const WorldParams& world) {
    FlatNodes flat = flatten(plant);

    // --- Adaptive N ---
    // Jacobi propagates pressure ~1 hop per sub-step, so N needs to be at least
    // the longest conduit chain for pressure to reach the deepest meristem.
    // Short plants don't need N=25 — they equilibrate in a handful of passes.
    // world.vascular_substeps is now treated as a CEILING, not a fixed count.
    const uint32_t N_ceiling    = std::max<uint32_t>(1, world.vascular_substeps);
    // Floor matters for branching: at a junction, canalization-bias-weighted
    // Jacobi distributes flow unequally between sibling sub-trees, and extract
    // runs in DFS order.  Those asymmetries only fully wash out after enough
    // passes.  Empirically, very shallow plants (chain depth 1) still need
    // ~10 iterations for sink concentrations to equalize within 1e-5.
    const uint32_t N_min         = 10;
    const uint32_t safety_margin = 5;
    const uint32_t N = std::min(N_ceiling,
                                std::max(N_min, flat.max_conduit_depth + safety_margin));

    // --- Part A: Budget snapshot (computed once) ---
    std::vector<VascularBudget> budgets;
    budgets.reserve(flat.all.size());
    for (Node* n : flat.all) {
        budgets.push_back(compute_budget(*n, g, world));
    }

    // --- Part A': Jacobi edge coefficients (computed once) ---
    // Every value inside EdgeJacobiCoeffs is constant across the N sub-steps
    // (geometry + canalization bias don't change inside vascular).  Hoisting
    // these out of the inner loop is pure redundancy elimination — the older
    // code recomputed phloem_capacity, xylem_capacity, and get_bias_multiplier
    // on every sub-step for every edge.
    struct JacobiEdge {
        Node* parent;
        Node* child;
        EdgeJacobiCoeffs coeffs;
    };
    std::vector<JacobiEdge> jacobi_edges;
    for (Node* n : flat.conduits) {
        for (Node* child : n->children) {
            if (child->phloem() || child->xylem()) {
                jacobi_edges.push_back({n, child, compute_edge_jacobi_coeffs(*n, *child, g)});
            }
        }
    }

    // --- Part A'': Radial per-node coefficients (computed once) ---
    // radial_permeability_* and capacities are functions of node radius only,
    // also constant across the N-loop.  Stored in parallel with flat.conduits.
    std::vector<NodeRadialCoeffs> radial_coeffs;
    radial_coeffs.reserve(flat.conduits.size());
    for (Node* n : flat.conduits) {
        radial_coeffs.push_back(compute_node_radial_coeffs(*n, g));
    }

    // --- Part B: Sub-step loop (N iterations) ---
    //
    // Order: inject → Jacobi → radial → extract.
    //
    // Earlier the order was inject → radial → extract → Jacobi, which had a
    // pathological failure mode for cytokinin: root apicals injected cyto
    // into their parent root's xylem, then radial flow immediately absorbed
    // most of it into that root's local_env (because radial equilibrates
    // between xylem and local every sub-step), and Jacobi only ran AFTER
    // that absorption.  The result was ~37 mg of cytokinin stranded in
    // deep root local pools and 0 reaching any stem xylem.
    //
    // Running Jacobi before radial gives the pressure wave a chance to
    // propagate longitudinally (root.xylem → seed.xylem → stem.xylem → ...)
    // before radial absorption at each node pulls chemical out into local.
    // Extract runs last so sinks (meristems, leaves) see the post-Jacobi,
    // post-radial conduit state when pulling their demand.
    PerfStats* perf = plant.perf();  // may be null when perf logging is off

    for (uint32_t iter = 0; iter < N; ++iter) {
        // Step 1: Inject at sources (leaves for sugar, roots for water via
        // root pressure, root apicals for cytokinin).
        {
            ScopedTimer t(perf ? perf->vascular_inject_ms : ScopedTimer::dummy);
            for (size_t i = 0; i < flat.all.size(); ++i) {
                const VascularBudget& b = budgets[i];
                if (b.sugar_supply > 0.0f
                    || b.water_supply > 0.0f
                    || b.cytokinin_supply > 0.0f) {
                    inject_step(*flat.all[i], b, N, g);
                }
            }
        }

        // Step 2: Longitudinal Jacobi across every conduit edge.
        // Uses precomputed per-edge coeffs (see Part A' above).
        {
            ScopedTimer t(perf ? perf->vascular_jacobi_ms : ScopedTimer::dummy);
            for (const JacobiEdge& e : jacobi_edges) {
                jacobi_step_impl(*e.parent, *e.child, e.coeffs);
            }
        }

        // Step 3: Radial flow on every conduit (stem, root).
        // Uses precomputed per-node coeffs (see Part A'' above).
        {
            ScopedTimer t(perf ? perf->vascular_radial_ms : ScopedTimer::dummy);
            for (size_t i = 0; i < flat.conduits.size(); ++i) {
                radial_flow_step_impl(*flat.conduits[i], N, radial_coeffs[i]);
            }
        }

        // Step 4: Extract at sinks (meristems for sugar, leaves/meristems for water).
        {
            ScopedTimer t(perf ? perf->vascular_extract_ms : ScopedTimer::dummy);
            for (size_t i = 0; i < flat.all.size(); ++i) {
                const VascularBudget& b = budgets[i];
                if (b.sugar_demand > 0.0f || b.water_demand > 0.0f) {
                    extract_step(*flat.all[i], b, N, g);
                }
            }
        }
    }
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
            // Root pressure: root apicals absorb water from soil, keep a local
            // reserve, and actively pump surplus into the nearest upstream
            // xylem.  This models osmotic root pressure — the main mechanism
            // driving xylem flow upward in young plants before they have
            // enough leaves to generate transpiration pull.
            const float water_reserve = g.root_water_reserve_fraction * cap_w;
            b.water_supply = std::max(0.0f, water - water_reserve);
            // Root apicals produce cytokinin to signal the rest of the plant —
            // the cytokinin's job is to travel up xylem to the shoot, not to
            // accumulate locally.  Inject all of it into parent xylem.
            //
            // Note: the older implementation kept a reserve equal to
            // root_cytokinin_inhibition_threshold (0.15 g) in an attempt to
            // preserve the axillary-activation inhibition check in
            // RootApicalNode::can_activate.  That turned out to be wrong:
            // per-tick production is ~1.3e-4 g, three orders of magnitude
            // below the inhibition threshold, so local cytokinin never
            // reaches anywhere near the threshold under any policy.  The
            // reserve effectively pinned b.cytokinin_supply = 0, leaving
            // every xylem pool in the plant empty of cytokinin and starving
            // shoot meristems (observed at tick 483: all 91 nodes had
            // xylem_cytokinin = 0, primary SAM had local cyto = 0 exactly,
            // SAM growth_fraction hard-gated to 0).
            b.cytokinin_supply = n.local().chemical(ChemicalID::Cytokinin);
            break;
        }
        case NodeType::ROOT: {
            // Established root nodes also generate root pressure: they absorb
            // water from soil (via RootNode::absorb_water → local_env) and
            // actively load surplus into their own xylem.  Unlike root apicals
            // which inject into their parent's xylem (they have no xylem of
            // their own), root nodes load their own xylem directly — this is
            // handled in inject_step by checking self.xylem() before falling
            // back to walk-up.
            const float water_reserve = g.root_water_reserve_fraction * cap_w;
            b.water_supply = std::max(0.0f, water - water_reserve);
            break;
        }
        case NodeType::STEM: {
            // The seed is a STEM with no parent.  It's initialized with
            // g.seed_sugar (48 g by default) — the heterotrophic reserve
            // the seedling lives on before leaves produce enough to feed
            // the plant themselves.  Real seeds actively mobilize these
            // reserves into the phloem (cotyledon → phloem transport,
            // driven by enzymes and active sucrose loaders), not just
            // slowly leak via diffusion.
            //
            // Without explicit mobilization, the seed's 48 g sits in
            // local_env and dribbles into phloem only via passive radial
            // flow (capacity ~0.3 mg for a 1.7cm-radius seed stem), which
            // massively throttles how fast the plant can access its
            // reserves.  Observed: after 484 ticks a plant had still only
            // consumed ~2.5 g — most of the sugar trapped in seed.local.
            //
            // Fix: seed actively loads sugar into its own phloem, above
            // a reserve fraction, analogous to how roots pump water (root
            // pressure) and leaves pump sugar.  Non-seed stems get no
            // sugar_supply — they're conduits, not sugar sources.
            if (!n.parent) {
                const float reserve = g.leaf_reserve_fraction_sugar * cap_s;
                b.sugar_supply = std::max(0.0f, sugar - reserve);
            }
            break;
        }
        default: break;
    }
    return b;
}

void inject_step(Node& source, const VascularBudget& b, uint32_t N, const Genome& /*g*/) {
    if (N == 0) return;

    // Sugar injection: leaves push into parent phloem (active pump).  For
    // the seed (which has its own phloem), load its own phloem directly
    // so the seedling's reserve mobilization feeds both shoot and root
    // chains via Jacobi.  Same pattern as root nodes pumping water into
    // their own xylem.
    if (b.sugar_supply > 0.0f) {
        TransportPool* target_pool = source.phloem();
        if (!target_pool) target_pool = source.nearest_phloem_upstream();
        if (target_pool) {
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

    // Water injection (root pressure): roots and root apicals actively pump
    // water from local_env into xylem.  Root nodes load their own xylem;
    // root apicals (no xylem of their own) fall back to walk-up parent's
    // xylem.  This creates positive pressure at the bottom of the xylem
    // stream, pushing water upward through Jacobi even when leaves aren't
    // transpiring enough to generate cohesion-tension pull from above.
    if (b.water_supply > 0.0f) {
        TransportPool* target_pool = source.xylem();
        if (!target_pool) target_pool = source.nearest_xylem_upstream();
        if (target_pool) {
            const float slice  = b.water_supply / static_cast<float>(N);
            const float actual = std::min(slice, source.local().chemical(ChemicalID::Water));
            source.local().chemical(ChemicalID::Water) -= actual;
            target_pool->chemical(ChemicalID::Water)   += actual;
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

// Public radial_flow_step — thin wrapper around the precomputed _impl.  Tests
// and any external caller that doesn't have a NodeRadialCoeffs handy go through
// this path, which computes the coeffs on the fly.  The hot path in
// vascular_sub_stepped() precomputes coeffs once and calls _impl N times.
//
// Design note on radial flow: cytokinin is deliberately NOT radially exchanged
// between xylem and local_env.  Attempted a cytokinin own-gradient radial flow
// first; it still drained cytokinin into root local pools because local_cap
// (≈ water_cap of the root volume) is 3+ orders of magnitude larger than
// xylem_cap, so the equilibrium concentration the radial flow chases lands
// almost all the chemical in the local side regardless of transport physics.
// Real xylem sap does exchange with surrounding parenchyma (pit membrane
// permeability), but at a tiny rate that our per-sub-step radial formula can't
// capture faithfully.  Keeping cytokinin confined to xylem + explicit extract
// is a reasonable simplification that gets the transport signal through to the
// shoot.
void radial_flow_step(Node& n, uint32_t N, const Genome& g) {
    const NodeRadialCoeffs c = compute_node_radial_coeffs(n, g);
    radial_flow_step_impl(n, N, c);
}

// Public jacobi_step — thin wrapper around the precomputed _impl.  See
// jacobi_step_impl() for the per-sub-step physics.  Design notes preserved
// from the earlier inline implementation:
//
// * Removed max_move = 0.5 × min(cap) clamp.  It was throttling flow by up to
//   5 orders of magnitude when pressure gradients were large (e.g.,
//   over-pressurized root xylem trying to drain to seed xylem).  The
//   source-chemical clamp is still sufficient to prevent pools going negative.
//
// * Cytokinin in xylem uses its OWN pressure-driven Jacobi, independent of
//   water flow.  Previously it "rode water" proportionally, which created a
//   pathological trap: root pressure pumped water into root xylem, radial
//   flow dragged water (and cytokinin) back into root local, and since radial
//   flow is water-coupled, the cytokinin never escaped the root zone (~90% of
//   plant cytokinin ended up stranded in root.local_env, 0.09% reaching
//   SAMs).  An independent cytokinin pressure-Jacobi lets the signal
//   propagate upward from root apicals (source) to SAMs (sinks) along its own
//   concentration gradient — biologically sensible for xylem sap diffusion.
//
// * Per-edge instrumentation (tick_edge_flux / tick_edge_cap) is accumulated
//   inside _impl for all N sub-steps so the "Transport Capacity Used"
//   overlay reads the full-tick flux ÷ full-tick cap correctly.
void jacobi_step(Node& parent, Node& child, const Genome& g) {
    const EdgeJacobiCoeffs c = compute_edge_jacobi_coeffs(parent, child, g);
    jacobi_step_impl(parent, child, c);
}

// Returns the length of the internode (distance from parent to this node)
// FOR CAPACITY CALCULATIONS.  For the seed (no parent), returns 1.0 so the
// seed has a meaningful capacity.
//
// Floor at 5 × tip_offset (5mm at default geometry).  Freshly-spawned
// internodes have geometric length = tip_offset (1mm), which gives them
// sub-µg transport capacity and makes them flow bottlenecks in the chain.
// This is an artifact of our discretization — a real plant's xylem is a
// continuous pipe whose flow rate doesn't care how we chop it into
// segments.  By flooring the effective length, newly-spawned short
// internodes get "plant segment" capacity rather than their literal
// geometric tiny length.  Geometry (for rendering, physics, droop,
// phototropism etc) is unaffected.
static float node_length(const Node& n, const Genome& g) {
    if (!n.parent) return 1.0f;
    const float geometric = glm::length(n.offset);
    // Floor at 2 × tip_offset.  5× floored too aggressively — with bigger
    // caps, radial_flow_step's max_equalize_volume = 0.5 × min(cap) also
    // scales up, which makes root xylems absorb cytokinin into local faster
    // than Jacobi can propagate it upward, starving shoot meristems.
    // 2× gives a modest capacity boost without over-amplifying radial.
    const float floor = 2.0f * g.tip_offset;
    return std::max(floor, geometric);
}

float phloem_capacity(const Node& n, const Genome& g) {
    if (!n.phloem()) return 0.0f;
    const float L = node_length(n, g);
    return 3.14159265358979f * n.radius * n.radius * L * g.phloem_fraction;
}

float xylem_capacity(const Node& n, const Genome& g) {
    if (!n.xylem()) return 0.0f;
    const float L = node_length(n, g);
    return 3.14159265358979f * n.radius * n.radius * L * g.xylem_fraction;
}

} // namespace botany
