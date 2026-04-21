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

void vascular_sub_stepped(Plant& plant, const Genome& g, const WorldParams& world) {
    const uint32_t N = std::max<uint32_t>(1, world.vascular_substeps);

    FlatNodes flat = flatten(plant);

    // --- Part A: Budget snapshot (computed once) ---
    std::vector<VascularBudget> budgets;
    budgets.reserve(flat.all.size());
    for (Node* n : flat.all) {
        budgets.push_back(compute_budget(*n, g, world));
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
    for (uint32_t iter = 0; iter < N; ++iter) {
        // Step 1: Inject at sources (leaves for sugar, roots for water via
        // root pressure, root apicals for cytokinin).
        for (size_t i = 0; i < flat.all.size(); ++i) {
            const VascularBudget& b = budgets[i];
            if (b.sugar_supply > 0.0f
                || b.water_supply > 0.0f
                || b.cytokinin_supply > 0.0f) {
                inject_step(*flat.all[i], b, N, g);
            }
        }

        // Step 2: Longitudinal Jacobi across every conduit edge.
        for (Node* n : flat.conduits) {
            for (Node* child : n->children) {
                if (child->phloem() || child->xylem()) {
                    jacobi_step(*n, *child, g);
                }
            }
        }

        // Step 3: Radial flow on every conduit (stem, root).
        for (Node* n : flat.conduits) {
            radial_flow_step(*n, N, g);
        }

        // Step 4: Extract at sinks (meristems for sugar, leaves/meristems for water).
        for (size_t i = 0; i < flat.all.size(); ++i) {
            const VascularBudget& b = budgets[i];
            if (b.sugar_demand > 0.0f || b.water_demand > 0.0f) {
                extract_step(*flat.all[i], b, N, g);
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

            // Cytokinin is NOT subject to radial flow between xylem and
            // local_env.  It stays in xylem until something extracts it
            // (sinks via extract_step) or propagates it longitudinally
            // (Jacobi).
            //
            // Attempted an own-gradient radial flow first; it still
            // drained cytokinin into root local pools because
            // local_cap (≈ water_cap of the root volume) is 3+ orders
            // of magnitude larger than xylem_cap (just the vessel
            // cross-section × length), so the equilibrium concentration
            // the radial flow chases lands almost all the chemical in
            // the local side regardless of transport physics.
            //
            // Real xylem sap does exchange with surrounding parenchyma
            // (pit membrane permeability), but at a tiny rate that our
            // per-sub-step radial formula can't capture faithfully.
            // Keeping cytokinin confined to xylem + explicit extract is
            // a reasonable simplification that gets the transport signal
            // through to the shoot.
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
                // Removed max_move = 0.5 * min(cap) clamp.  It was throttling
                // flow by up to 5 orders of magnitude when pressure gradients
                // were large (e.g., over-pressurized root xylem trying to
                // drain to seed xylem).  The source-chemical clamp below
                // is still sufficient to prevent pools going negative.
                if (flow > 0.0f) flow = std::min(flow, pp->chemical(ChemicalID::Sugar));
                else             flow = -std::min(-flow, cp->chemical(ChemicalID::Sugar));

                pp->chemical(ChemicalID::Sugar) -= flow;
                cp->chemical(ChemicalID::Sugar) += flow;

                // --- Per-edge instrumentation (phloem sugar) ---
                // Accumulate signed flux and theoretical cap across all N sub-steps.
                // flow is positive when sugar moves parent→child, negative when reversed.
                // cap_this_substep = base_cond (bottleneck capacity, before bias).
                // Over the full tick this sums to base_cond * N, the total theoretical
                // throughput for "Transport Capacity Used" overlay (flux / cap).
                const size_t SI = static_cast<size_t>(ChemicalID::Sugar);
                parent.tick_edge_flux[SI][&child] += flow;
                parent.tick_edge_cap[SI][&child]  += base_cond;
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
                // Removed max_move clamp here too — same reasoning as phloem.
                // Source-chemical clamp below still prevents negative pools.
                if (flow > 0.0f) flow = std::min(flow, pp->chemical(ChemicalID::Water));
                else             flow = -std::min(-flow, cp->chemical(ChemicalID::Water));

                // Water move.
                pp->chemical(ChemicalID::Water) -= flow;
                cp->chemical(ChemicalID::Water) += flow;

                // Cytokinin in xylem has its OWN concentration gradient —
                // not tied to water flow.  Previously it "rode water"
                // proportionally, which created a pathological trap:
                // root pressure pumped water into root xylem, radial flow
                // then dragged water (and cytokinin) back into root local,
                // and since radial flow is also water-coupled, the cytokinin
                // never escaped the root zone.  ~90% of plant cytokinin
                // ended up stranded in root.local_env, 0.09% reaching SAMs.
                //
                // Giving cytokinin its own pressure-driven Jacobi (same
                // formula shape as phloem sugar) lets it propagate upward
                // from its source (root apicals) to its sinks (SAMs) along
                // its own concentration gradient.  Biologically sensible —
                // cytokinin in xylem sap diffuses along its own potential
                // independent of bulk water flow.
                const float cyto_pressure_p = pp->chemical(ChemicalID::Cytokinin) / cap_p;
                const float cyto_pressure_c = cp->chemical(ChemicalID::Cytokinin) / cap_c;
                float cyto_flow = conductance * (cyto_pressure_p - cyto_pressure_c);
                if (cyto_flow > 0.0f) cyto_flow = std::min(cyto_flow, pp->chemical(ChemicalID::Cytokinin));
                else                  cyto_flow = -std::min(-cyto_flow, cp->chemical(ChemicalID::Cytokinin));

                pp->chemical(ChemicalID::Cytokinin) -= cyto_flow;
                cp->chemical(ChemicalID::Cytokinin) += cyto_flow;
            }
        }
    }
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
