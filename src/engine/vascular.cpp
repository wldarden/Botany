// src/engine/vascular.cpp — Xylem & phloem vascular transport.
// Two-phase global pass: post-order aggregates supply/demand,
// pre-order distributes actual flow. Runs before the DFS tree walk.
#include "engine/vascular.h"
#include "engine/plant.h"
#include "engine/node/node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"
#include "engine/sugar.h"
#include "engine/genome.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <filesystem>
#include <glm/geometric.hpp>

namespace botany {

bool has_vasculature(const Node& n, const Genome& g) {
    if (!n.parent) return true;  // seed is always a vascular junction
    // Vascular admission is radius-based. All stem/root nodes with radius >=
    // vascular_radius_threshold qualify. initial_radius (0.015 dm) > threshold (0.01 dm)
    // so newly spawned internodes join immediately. Leaves and meristems never join
    // (they rely on last-mile local diffusion from the nearest vascular node).
    if (n.type == NodeType::STEM || n.type == NodeType::ROOT)
        return n.radius >= g.vascular_radius_threshold;
    return false;
}

// Vascular pipe capacity based on cross-sectional area.
static float pipe_capacity(const Node& n, float conductance) {
    return 3.14159f * n.radius * n.radius * conductance;
}

// Per-node temporary data for the two-phase pass.
struct VascNodeInfo {
    Node* node = nullptr;
    int parent_idx = -1;          // index into flat array
    std::vector<int> child_idxs;  // indices of children in flat array

    float supply = 0.0f;          // what this subtree can provide
    float demand = 0.0f;          // what this subtree wants
    float capacity = 1e30f;       // pipe capacity at this node (bottleneck)
    bool is_conduit = false;      // mature stem/root — passes through
};

// Build a flat array of VascNodeInfo in pre-order (DFS).
static void build_flat(Node* node, int parent_idx,
                       std::vector<VascNodeInfo>& flat) {
    int my_idx = static_cast<int>(flat.size());
    flat.push_back({node, parent_idx, {}, 0, 0, 1e30f, false});
    if (parent_idx >= 0) {
        flat[parent_idx].child_idxs.push_back(my_idx);
    }
    for (Node* child : node->children) {
        build_flat(child, my_idx, flat);
    }
}

// Run one vascular system (phloem or xylem) for one chemical.
static const char* chem_name(ChemicalID id) {
    switch (id) {
        case ChemicalID::Sugar:     return "Sugar";
        case ChemicalID::Water:     return "Water";
        case ChemicalID::Cytokinin: return "Cytokinin";
        default:                    return "Other";
    }
}

static const char* node_type_name(NodeType t) {
    switch (t) {
        case NodeType::STEM:        return "STEM";
        case NodeType::ROOT:        return "ROOT";
        case NodeType::LEAF:        return "LEAF";
        case NodeType::APICAL:      return "SA";
        case NodeType::ROOT_APICAL: return "RA";
        default:                    return "?";
    }
}


// ── Münch pressure-flow phloem helpers ────────────────────────────────────────

// Phloem ring cross-section area — used ONLY for pipe capacity (flow_vol = velocity × ring_area).
// Represents the full active conducting layer: sieve tubes + companion cells + inner
// bark + cambium. Heartwood excluded. Scales linearly with r.
// t is clamped to r so tiny young nodes (small radius) cannot produce negative ring area
// when the fixed ring thickness temporarily exceeds their outer radius.
static float phloem_ring_area(float r, float t) {
    constexpr float kPi = 3.14159265f;
    t = std::min(t, r);
    return kPi * (2.0f * r * t - t * t);
}

// Effective pipe radius for an edge between two nodes.  STEM/ROOT nodes carry
// real vasculature; LEAF nodes intentionally have radius == 0 (they track
// leaf_size instead), and meristems have only tiny tip radii.  For edges
// involving leaves/meristems the pipe width is determined by the conduit side
// — the leaf or meristem taps into the parent stem's vascular bundle rather
// than carrying its own.  For conduit↔conduit edges the bottleneck is min(both).
static float edge_pipe_radius(const Node& a, const Node& b) {
    auto is_conduit = [](const Node& n) {
        return n.type == NodeType::STEM || n.type == NodeType::ROOT;
    };
    bool ac = is_conduit(a), bc = is_conduit(b);
    if (ac && bc) return std::min(a.radius, b.radius);
    if (ac)       return a.radius;
    if (bc)       return b.radius;
    return std::max({a.radius, b.radius, 1e-4f});  // degenerate fallback
}

// Total node volume — used for concentration (pressure) and sugar cap.
// Using full cylinder (not ring) prevents the ~85× concentration asymmetry that would
// make stem conduits appear as higher-pressure than leaf sources at realistic sugar levels.
// Seed has zero-length offset (it IS the trunk base); give it a minimum generous volume.
static float node_volume(const Node& n, const WorldParams& world) {
    constexpr float kPi = 3.14159265f;
    constexpr float kSeedMinLength = 0.2f;  // dm — seed is a storage organ, not a pipe
    switch (n.type) {
        case NodeType::STEM:
        case NodeType::ROOT: {
            float length = glm::length(n.offset);
            if (!n.parent) length = std::max(length, kSeedMinLength);  // seed fix
            length = std::max(length, 1e-4f);  // general safety floor
            return kPi * n.radius * n.radius * length;
        }
        case NodeType::LEAF: {
            float ls = n.as_leaf()->leaf_size;
            return ls * ls * world.leaf_thickness;
        }
        case NodeType::APICAL:
        case NodeType::ROOT_APICAL: {
            return (4.0f / 3.0f) * kPi * n.radius * n.radius * n.radius;
        }
        default:
            return 0.0f;
    }
}

// Phloem osmotic pressure at a vascular node.
// pressure = (sugar / node_vol) × osmotic_coeff × water_frac
// Using node_volume (full cylinder) gives realistic concentration gradients.
// water_frac couples drought to phloem: a dehydrated source cannot build pressure.
static float compute_phloem_pressure(const Node& n, const Genome& g, const WorldParams& world) {
    float vol = node_volume(n, world);
    if (vol <= 0.0f) return 0.0f;
    float sugar_conc = n.chemical(ChemicalID::Sugar) / vol;  // g/dm³
    float wc = water_cap(n, g);
    float water_frac = (wc > 0.0f)
        ? std::clamp(n.chemical(ChemicalID::Water) / wc, 0.0f, 1.0f)
        : 0.0f;
    return sugar_conc * g.phloem_osmotic_coefficient * water_frac;
}

// Permeability of a node's sieve tube membrane — controls how much sugar
// crosses from pipe lumen into local tissue per unit concentration gradient.
static float unloading_permeability(NodeType t, const Genome& g) {
    switch (t) {
        case NodeType::APICAL:      return g.phloem_unloading_meristem;
        case NodeType::ROOT_APICAL: return g.phloem_unloading_meristem;
        case NodeType::LEAF:        return g.phloem_unloading_leaf;
        case NodeType::ROOT:        return g.phloem_unloading_root;
        case NodeType::STEM:        return g.phloem_unloading_stem;
        default:                    return g.phloem_unloading_stem;
    }
}

// ── phloem_resolve: Münch pressure-flow with Jacobi simultaneous resolution ──
//
// Three phases run in sequence:
//   1. Leaf loading pass (pre-Jacobi) — concentration-gradient transfer from each
//      leaf into its parent vascular conduit.  Velocity-capped + equalization-capped.
//   2. Jacobi pipe-network iterations — all edges resolved simultaneously for
//      phloem_iterations rounds.  No ordering bias; fairness scaling prevents
//      over-draft.
//   3. Meristem unloading pass (post-Jacobi) — concentration-gradient transfer
//      from each vascular conduit into its meristem child.
//
// Leaves and meristems are NOT in the Jacobi pipe graph; they interact with the
// network only via the endpoint passes.  All other chemicals (auxin, water, etc.)
// are handled by separate passes and are unaffected here.
void phloem_resolve(Plant& plant, const Genome& g, const WorldParams& world) {
    Node* seed = plant.seed_mut();
    if (!seed) return;

    // Build flat pre-order traversal (same helper used by xylem_resolve).
    std::vector<VascNodeInfo> flat;
    flat.reserve(plant.node_count());
    build_flat(seed, -1, flat);
    int N = static_cast<int>(flat.size());
    if (N == 0) return;

    // ── Logging setup ─────────────────────────────────────────────────────────
    // Allocate tracking arrays only when phloem_debug_log is enabled so the
    // normal (non-logging) path has zero overhead.
    const bool logging = world.phloem_debug_log;
    std::vector<float> log_pre_sugar;       // sugar at each node before phloem runs
    std::vector<float> log_loaded;          // sugar loaded from leaf children (indexed by parent)
    std::vector<float> log_unloaded;        // sugar unloaded to meristem children (indexed by parent)
    std::vector<float> log_flow_in;         // sugar received via Jacobi (accumulated across iterations)
    std::vector<float> log_flow_out;        // sugar sent via Jacobi (accumulated across iterations)

    // Per-edge Jacobi totals for phloem_edge_log.csv.
    // Edge key = child index (each non-seed node has exactly one parent_idx).
    // desired = pressure-driven flow out of sender; unload = what crossed the
    // membrane into receiver. Both summed across phloem_iterations.
    std::vector<float> edge_desired_to_parent;
    std::vector<float> edge_desired_from_parent;
    std::vector<float> edge_unload_to_parent;
    std::vector<float> edge_unload_from_parent;
    std::vector<float> edge_final_pressure;  // last-iter pressure at each node
    if (logging) {
        log_pre_sugar.resize(N, 0.0f);
        log_loaded.resize(N, 0.0f);
        log_unloaded.resize(N, 0.0f);
        log_flow_in.resize(N, 0.0f);
        log_flow_out.resize(N, 0.0f);
        edge_desired_to_parent.resize(N, 0.0f);
        edge_desired_from_parent.resize(N, 0.0f);
        edge_unload_to_parent.resize(N, 0.0f);
        edge_unload_from_parent.resize(N, 0.0f);
        edge_final_pressure.resize(N, 0.0f);
        for (int i = 0; i < N; ++i)
            log_pre_sugar[i] = flat[i].node->chemical(ChemicalID::Sugar);
    }

    // ── Phase 1: Leaf loading pass ────────────────────────────────────────────
    // Each leaf transfers sugar into its parent conduit proportional to the
    // concentration gradient × permeability.  Both velocity cap and equalization
    // cap apply.  Leaves are NOT in the Jacobi network — this is their only path
    // to export.
    for (int i = 0; i < N; ++i) {
        Node& leaf = *flat[i].node;
        if (leaf.type != NodeType::LEAF) continue;
        if (flat[i].parent_idx < 0) continue;
        Node& parent = *flat[flat[i].parent_idx].node;
        if (!has_vasculature(parent, g)) continue;

        float leaf_vol   = node_volume(leaf, world);
        float parent_vol = node_volume(parent, world);
        if (leaf_vol <= 0.0f || parent_vol <= 0.0f) continue;

        float leaf_conc   = leaf.chemical(ChemicalID::Sugar) / leaf_vol;
        float parent_conc = parent.chemical(ChemicalID::Sugar) / parent_vol;
        float gradient    = leaf_conc - parent_conc;
        if (gradient <= 0.0f) continue;  // leaf not above parent — no loading

        float r_eff      = edge_pipe_radius(leaf, parent);
        float ring_area  = phloem_ring_area(r_eff, world.phloem_ring_thickness);
        float velocity   = std::min(gradient * g.conductance_per_pressure,
                                    world.max_phloem_velocity);
        float flow_vol   = velocity * ring_area;               // dm³/tick
        float max_vel    = flow_vol * leaf_conc;               // g — velocity-limited

        float equil_conc = (leaf.chemical(ChemicalID::Sugar) + parent.chemical(ChemicalID::Sugar))
                           / std::max(leaf_vol + parent_vol, 1e-8f);
        float max_equil  = std::max(0.0f,
                               leaf.chemical(ChemicalID::Sugar) - equil_conc * leaf_vol);

        float load = std::min(max_vel, max_equil);
        load = std::min(load, leaf.chemical(ChemicalID::Sugar));  // safety

        // Storage cap for the parent conduit uses the per-type sugar_cap() model
        // (stem = volume × sugar_storage_density_wood, leaf = area × leaf-density,
        // seed = max of children-sum and seed_sugar, meristem = fixed g.sugar_cap_meristem).
        // Keeps phloem caps consistent with the rest of the sim instead of using a
        // volumetric-only cap that made young internodes choke at ~0.002 g.
        float parent_cap  = sugar_cap(parent, g);
        float parent_room = parent_cap - parent.chemical(ChemicalID::Sugar);
        load = std::min(load, std::max(0.0f, parent_room));

        if (load > 1e-8f) {
            leaf.chemical(ChemicalID::Sugar)   -= load;
            parent.chemical(ChemicalID::Sugar) += load;
            if (logging) log_loaded[flat[i].parent_idx] += load;
        }
    }

    // ── Phase 2: Demand-driven source-to-sink allocation ──────────────────────
    // Replaces the previous pairwise-Jacobi Münch resolution.  Pairwise
    // equilibration attenuates sugar signals exponentially along long shoot
    // chains — observed in a 2596-tick run where the primary SA at the tip
    // of a 30-stem chain had sugar=0 for 1900+ ticks despite the seed holding
    // 0.23 g.  The demand-driven pattern (mirroring xylem_resolve) eliminates
    // hop-by-hop attenuation: sinks pull sugar directly from the plant-wide
    // supply pool in proportion to their demand.
    //
    // Classification (computed once per tick):
    //   SOURCE (supply):
    //     - STEM / ROOT / seed with sugar > phloem_reserve — surplus above the
    //       parenchyma reserve is available to flow to sinks.
    //   SINK (demand):
    //     - Active APICAL / ROOT_APICAL wanting to fill to meristem_sink_fraction
    //       × cap per tick.  Dormant meristems have zero demand.
    //
    // STEM / ROOT / seed are SOURCES ONLY — they supply surplus above their
    // reserve.  They do NOT actively pull sugar to refill their parenchyma.
    // Biologically this matches real phloem: mature stem/root parenchyma
    // doesn't actively unload from the sieve tube stream.  Sugar gets
    // deposited into stems/roots by leaf loading (Phase 1) and by transit
    // through neighboring nodes.  Maintenance consumes what's there; if a
    // conduit empties, it lives on its 8000-tick (stem) / 6000-tick (root)
    // starvation threshold — stand-in for starch reserves we'll model later.
    //
    // Earlier iteration of this code had conduits demand `reserve - sugar`.
    // That starved meristems because ~1000 conduit nodes each demanding tiny
    // refill dominated the ~16 meristems each demanding 0.005g.  Active
    // meristems got ~0.00002g/tick — far below their 0.0005g maintenance —
    // so SAs stopped growing despite ample leaf production.  Removing
    // conduit demand restores the semantic: meristems are the pulling force,
    // conduits are pipes.
    //
    // LEAF loading was already handled in Phase 1 above.  Leaves keep their
    // remaining sugar for their own use (photosynthesis/respiration).
    //
    // Mass conservation: Σ source_deduct = delivered = Σ sink_add, by
    // construction of proportional scales.  Verified every tick by the SUMMARY
    // row's conservation_error column in phloem_log.csv.
    std::vector<float> supply(N, 0.0f);
    std::vector<float> demand(N, 0.0f);

    for (int i = 0; i < N; ++i) {
        Node& n = *flat[i].node;
        float cap = sugar_cap(n, g);
        if (cap <= 1e-8f) continue;
        float reserve = cap * g.phloem_reserve_fraction;
        float sug = n.chemical(ChemicalID::Sugar);

        if (n.type == NodeType::APICAL) {
            // Active shoot apicals are sinks; dormant ones are neither.
            if (n.as_apical()->active) {
                float target = cap * g.meristem_sink_fraction;
                demand[i] = std::max(0.0f, target - sug);
            }
        } else if (n.type == NodeType::ROOT_APICAL) {
            if (n.as_root_apical()->active) {
                float target = cap * g.meristem_sink_fraction;
                demand[i] = std::max(0.0f, target - sug);
            }
        } else if (n.type == NodeType::STEM || n.type == NodeType::ROOT) {
            // Conduit / storage tissue — SOURCE ONLY (no demand, see comment
            // block above).  Seed (STEM with no parent): sugar_cap() inflates
            // the seed cap to max(child_sum, seed_sugar) so the seed can hold
            // its initial 48g reserve.  But phloem_reserve_fraction × 48g =
            // 14.4g, which would lock the seed from donating below 14.4g
            // even though it's the only source in the plant.  For the seed
            // we use the network-proportional cap (sum of children's
            // sugar_cap) instead, so the reserve scales with the connected
            // network, not the initialization constant.
            if (n.parent == nullptr) {
                float network_cap = 0.0f;
                for (const Node* child : n.children)
                    network_cap += sugar_cap(*child, g);
                network_cap = std::max(network_cap, g.sugar_cap_minimum);
                reserve = network_cap * g.phloem_reserve_fraction;
            }
            if (sug > reserve) supply[i] = sug - reserve;
            // No else branch — conduits don't pull.
        }
        // LEAF nodes: participated in Phase 1.  They keep their remaining
        // sugar for their own use (photosynthesis-triggered growth, respiration).
    }

    float total_supply = 0.0f, total_demand = 0.0f;
    for (int i = 0; i < N; ++i) {
        total_supply += supply[i];
        total_demand += demand[i];
    }
    float delivered = std::min(total_supply, total_demand);

    if (delivered > 1e-12f) {
        float supply_scale = (total_supply > 1e-12f) ? delivered / total_supply : 0.0f;
        float demand_scale = (total_demand > 1e-12f) ? delivered / total_demand : 0.0f;

        for (int i = 0; i < N; ++i) {
            Node& n = *flat[i].node;
            if (supply[i] > 0.0f) {
                float deducted = supply[i] * supply_scale;
                n.chemical(ChemicalID::Sugar) -= deducted;
                if (logging) log_flow_out[i] += deducted;
            }
            if (demand[i] > 0.0f) {
                float added = demand[i] * demand_scale;
                n.chemical(ChemicalID::Sugar) += added;
                if (logging) log_flow_in[i] += added;
                // Meristems track "unloaded" separately so the SUMMARY log
                // still distinguishes meristem-sink delivery from conduit
                // refill.  Store under the meristem's parent id, matching the
                // old Phase 3 log convention.
                if ((n.type == NodeType::APICAL || n.type == NodeType::ROOT_APICAL)
                    && flat[i].parent_idx >= 0) {
                    log_unloaded[flat[i].parent_idx] += added;
                }
            }
        }
    }

    // ── Phloem debug log ──────────────────────────────────────────────────────
    // Writes debug/phloem_log.csv with one row per node and a SUMMARY row.
    // File is truncated on tick 0 and appended thereafter.
    // Conservation check: SUMMARY.conservation_error = total_sugar_after - total_sugar_before.
    if (logging) {
        std::filesystem::create_directories("debug");
        auto mode = (world.current_tick == 0)
            ? (std::ios::out | std::ios::trunc)
            : (std::ios::out | std::ios::app);
        std::ofstream csv("debug/phloem_log.csv", mode);
        if (csv.is_open()) {
            if (world.current_tick == 0) {
                csv << "tick,node_id,node_type,parent_id,"
                       "sugar,volume,concentration,pressure,water_fraction,"
                       "sugar_loaded_from_leaf,sugar_unloaded_to_meristem,"
                       "sugar_flow_in,sugar_flow_out,net_flow\n";
            }

            float total_before = 0.0f;
            float total_after  = 0.0f;
            float total_loaded   = 0.0f;
            float total_unloaded = 0.0f;
            float total_flow_in  = 0.0f;
            float total_flow_out = 0.0f;

            for (int i = 0; i < N; ++i) {
                const Node& n      = *flat[i].node;
                float pre           = log_pre_sugar[i];
                float post          = n.chemical(ChemicalID::Sugar);
                float vol           = node_volume(n, world);
                float conc          = (vol > 0.0f) ? pre / vol : 0.0f;
                float pressure      = compute_phloem_pressure(n, g, world);  // uses current sugar
                float wc            = water_cap(n, g);
                float wfrac         = (wc > 0.0f)
                    ? std::clamp(n.chemical(ChemicalID::Water) / wc, 0.0f, 1.0f)
                    : 0.0f;
                float loaded        = log_loaded[i];
                float unloaded      = log_unloaded[i];
                float fin           = log_flow_in[i];
                float fout          = log_flow_out[i];
                float net           = fin - fout;
                int parent_id = (flat[i].parent_idx >= 0)
                    ? static_cast<int>(flat[flat[i].parent_idx].node->id)
                    : -1;

                csv << world.current_tick << ','
                    << n.id << ','
                    << node_type_name(n.type) << ','
                    << parent_id << ','
                    << pre << ','
                    << vol << ','
                    << conc << ','
                    << pressure << ','
                    << wfrac << ','
                    << loaded << ','
                    << unloaded << ','
                    << fin << ','
                    << fout << ','
                    << net << '\n';

                total_before   += pre;
                total_after    += post;
                total_loaded   += loaded;
                total_unloaded += unloaded;
                total_flow_in  += fin;
                total_flow_out += fout;
            }

            float conservation_error = total_after - total_before;
            csv << world.current_tick << ",SUMMARY,"
                << N << ','
                << total_before << ','
                << total_loaded << ','
                << total_unloaded << ','
                << total_flow_in << ','
                << conservation_error << '\n';
        }

        // ── Per-edge phloem log ───────────────────────────────────────────
        // One row per non-seed node: flow between it and its parent, summed
        // across all Jacobi iterations. edge_desired_to_parent[i] = flow from
        // i toward its parent (i was sender). edge_desired_from_parent[i] =
        // flow from parent into i (parent was sender).
        auto edge_mode = (world.current_tick == 0)
            ? (std::ios::out | std::ios::trunc)
            : (std::ios::out | std::ios::app);
        std::ofstream ecsv("debug/phloem_edge_log.csv", edge_mode);
        if (ecsv.is_open()) {
            if (world.current_tick == 0) {
                ecsv << "tick,parent_id,child_id,child_type,"
                        "pressure_parent,pressure_child,"
                        "conc_parent,conc_child,"
                        "desired_to_parent,unload_to_parent,"
                        "desired_from_parent,unload_from_parent,"
                        "net_child_gain\n";
            }
            for (int i = 1; i < N; ++i) {  // skip seed (no parent edge)
                int pi = flat[i].parent_idx;
                if (pi < 0) continue;
                const Node& parent = *flat[pi].node;
                const Node& child  = *flat[i].node;
                float p_vol = node_volume(parent, world);
                float c_vol = node_volume(child,  world);
                float p_conc = (p_vol > 0.0f) ? parent.chemical(ChemicalID::Sugar) / p_vol : 0.0f;
                float c_conc = (c_vol > 0.0f) ? child.chemical(ChemicalID::Sugar)  / c_vol : 0.0f;
                float net_child_gain = edge_unload_from_parent[i] - edge_desired_to_parent[i];
                ecsv << world.current_tick << ','
                     << parent.id << ',' << child.id << ','
                     << node_type_name(child.type) << ','
                     << edge_final_pressure[pi] << ','
                     << edge_final_pressure[i]  << ','
                     << p_conc << ',' << c_conc << ','
                     << edge_desired_to_parent[i] << ','
                     << edge_unload_to_parent[i]  << ','
                     << edge_desired_from_parent[i] << ','
                     << edge_unload_from_parent[i]  << ','
                     << net_child_gain << '\n';
            }
        }
    }
}

// ── xylem_resolve: demand-driven water + cytokinin transport ──────────────
//
// Real xylem is NOT local diffusion — it's bulk flow under transpiration tension.
// A transpiring leaf creates a negative pressure that propagates through the
// cohesion-tension water column and pulls water from roots in the same instant
// (classical plant physiology).  Our previous wfrac-gradient Jacobi model could
// not capture this: adjacent-pair equilibration bottlenecks flow at every edge,
// and water pools at roots while shoot tips stay drought-stunted no matter how
// many iterations we run.
//
// Replacing it with a source-to-sink flow network, mirroring how phloem_resolve
// handles sugar sources (leaves) and sinks (meristems):
//
//   1. Classify every node as SOURCE (has water to give), SINK (needs water),
//      or NEUTRAL (conduit/no participation).
//   2. Sum total supply and total demand plant-wide.
//   3. Deliver flow = min(total_supply, total_demand) — mass-conservative.
//   4. Each source contributes proportionally: source.water -= supply[i] × scale_s.
//      Each sink receives proportionally: sink.water += demand[j] × scale_d.
//   5. Cytokinin rides the water: each source also loses cyto proportional to the
//      water it sent; sinks receive cyto pooled across the xylem stream.
//
// Classifications:
//   SOURCE: ROOT, ROOT_APICAL — absorbed from soil, offer water above a reserve
//   SINK:   LEAF — transpiration + photosynthesis water cost + fill-to-cap demand
//           APICAL, STEM, ROOT (buffer demand to maintain growth turgor)
//   NEUTRAL: seed — it's a junction; its water replenishes via root subtree and
//           is drained by shoot subtree through the source/sink accounting.
//
// Mass conservation is trivially exact: Σ deducted = delivered = Σ added.
// Pipe capacity is NOT modeled explicitly; throughput is demand-limited and the
// xylem conductance is high enough that real tree-scale flows stay below capacity.
// If that becomes false we can add a pipe-throughput cap as a secondary limit.
void xylem_resolve(Plant& plant, const Genome& g, const WorldParams& world) {
    Node* seed = plant.seed_mut();
    if (!seed) return;

    std::vector<VascNodeInfo> flat;
    flat.reserve(plant.node_count());
    build_flat(seed, -1, flat);
    const int N = static_cast<int>(flat.size());
    if (N == 0) return;

    // Logging arrays — capture pre-amounts so SUMMARY rows show conservation.
    const bool logging = world.xylem_debug_log;
    std::vector<float> pre_water(N, 0.0f), pre_cyto(N, 0.0f);
    std::vector<float> edge_water_up(N, 0.0f);   // per-node net water gained (for logging only)
    std::vector<float> edge_water_down(N, 0.0f); // per-node water lost (for logging only)
    std::vector<float> edge_cyto_up(N, 0.0f);
    std::vector<float> edge_cyto_down(N, 0.0f);
    if (logging) {
        for (int i = 0; i < N; ++i) {
            pre_water[i] = flat[i].node->chemical(ChemicalID::Water);
            pre_cyto[i]  = flat[i].node->chemical(ChemicalID::Cytokinin);
        }
    }

    // ── Pre-compute water caps ─────────────────────────────────────────────
    std::vector<float> wcap(N, 0.0f);
    for (int i = 0; i < N; ++i) wcap[i] = water_cap(*flat[i].node, g);

    // ── Phase 1: compute per-node supply (sources) and demand (sinks) ──────
    //
    // SOURCE: ROOT, ROOT_APICAL offer water above a reserve fraction of their cap
    // (keeps some for their own turgor/growth).  Reserve fraction tuned low (0.3)
    // because roots are never short on water in moist soil.
    //
    // SINK: LEAF demand combines three terms:
    //   - transpiration this tick (area × transp_rate × light × stomatal)
    //   - photosynthesis water cost (production this tick × photo_water_ratio)
    //   - fill-to-cap buffer so wfrac doesn't sag between consumption events
    //
    // SAs, STEMs, and ROOTs have a turgor-buffer demand to reach a target wfrac
    // (0.5 for stems/SAs).  This keeps the whole vascular chain moist enough for
    // growth; real plants don't run stems bone dry either.
    //
    // Seed is neutral — no supply or demand of its own; water flows through it
    // from root subtree to shoot subtree as part of the source/sink accounting.
    constexpr float kRootReserveFrac   = 0.3f;  // roots keep 30% of cap as reserve
    constexpr float kStemTargetTurgor  = 0.5f;  // stems aim for 50% wfrac
    constexpr float kShootMeristemTurgor = 0.7f; // SAs want 70% wfrac for growth
    constexpr float kRootMeristemTurgor  = 0.5f; // RAs want 50% (less critical)

    std::vector<float> supply(N, 0.0f);
    std::vector<float> demand(N, 0.0f);

    for (int i = 0; i < N; ++i) {
        Node& n = *flat[i].node;
        float w = n.chemical(ChemicalID::Water);
        float cap = wcap[i];
        if (cap <= 1e-8f) continue;

        if (n.type == NodeType::ROOT || n.type == NodeType::ROOT_APICAL) {
            // Source side.  ROOT_APICAL is both absorbing tissue and xylem
            // stream — offer whatever's above the reserve floor.
            float reserve = cap * kRootReserveFrac;
            supply[i] = std::max(0.0f, w - reserve);
            // Still give RAs a small turgor-buffer demand so they don't drain
            // below the reserve during big transpiration pulls.
            if (n.type == NodeType::ROOT_APICAL) {
                float target = cap * kRootMeristemTurgor;
                demand[i] = std::max(0.0f, target - w);
            }
        } else if (n.type == NodeType::LEAF) {
            // Sink.  Anticipate this tick's consumption (transpiration + photo)
            // plus fill any deficit to cap.  Using last-tick's observed
            // light_exposure and wfrac-derived stomatal conductance.
            const LeafNode* leaf = n.as_leaf();
            float area = leaf->leaf_size * leaf->leaf_size;
            float light = leaf->light_exposure * world.light_level;
            float wfrac_now = w / cap;
            float stomatal = std::clamp(wfrac_now, 0.2f, 1.0f);
            float transp = area * g.transpiration_rate * light * stomatal;
            float production = area * world.sugar_production_rate * light * stomatal;
            float photo_water = production * g.photosynthesis_water_ratio;
            float buffer = std::max(0.0f, cap - w);
            demand[i] = transp + photo_water + buffer;
        } else if (n.type == NodeType::APICAL) {
            float target = cap * kShootMeristemTurgor;
            demand[i] = std::max(0.0f, target - w);
        } else if (n.type == NodeType::STEM) {
            // Seed (no parent) is neutral — stays out of the accounting so it
            // behaves as a passive junction.
            if (n.parent == nullptr) continue;
            float target = cap * kStemTargetTurgor;
            demand[i] = std::max(0.0f, target - w);
        }
    }

    float total_supply = 0.0f, total_demand = 0.0f;
    for (int i = 0; i < N; ++i) {
        total_supply += supply[i];
        total_demand += demand[i];
    }
    float delivered = std::min(total_supply, total_demand);

    // ── Phase 2: apply transfers proportionally ────────────────────────────
    // Mass conservation: Σ source_deduct = delivered = Σ sink_add.
    // The proportional scaling means every source sends the same FRACTION of its
    // offered supply, and every sink receives the same FRACTION of its demand.
    // This avoids prioritization decisions and matches bulk-flow physics where a
    // single-junction tree redistributes flow uniformly under full mixing.
    if (delivered > 1e-12f) {
        float supply_scale = (total_supply > 1e-12f) ? delivered / total_supply : 0.0f;
        float demand_scale = (total_demand > 1e-12f) ? delivered / total_demand : 0.0f;

        // ── Water transfers ────────────────────────────────────────────────
        for (int i = 0; i < N; ++i) {
            Node& n = *flat[i].node;
            if (supply[i] > 0.0f) {
                float deducted = supply[i] * supply_scale;
                n.chemical(ChemicalID::Water) -= deducted;
                if (logging) edge_water_down[i] += deducted;
            }
            if (demand[i] > 0.0f) {
                float added = demand[i] * demand_scale;
                n.chemical(ChemicalID::Water) += added;
                if (logging) edge_water_up[i] += added;
            }
        }

        // ── Cytokinin rides the water stream ──────────────────────────────
        // Each source loses cyto proportional to the water it sent, using its
        // pre-transfer cyto-per-water concentration.  Pool all the transported
        // cyto at the seed "junction" and redistribute to sinks in the same
        // proportion the sinks received water.  Mass-conservative.
        float total_cyto_transported = 0.0f;
        for (int i = 0; i < N; ++i) {
            Node& n = *flat[i].node;
            if (supply[i] <= 0.0f) continue;
            float water_sent = supply[i] * supply_scale;
            float water_before = pre_water[i];  // use pre-tick water (matches cyto basis)
            if (water_before <= 1e-8f) continue;
            float cyto_per_water = pre_cyto[i] / water_before;
            float cyto_sent = std::min(water_sent * cyto_per_water, n.chemical(ChemicalID::Cytokinin));
            if (cyto_sent <= 0.0f) continue;
            n.chemical(ChemicalID::Cytokinin) -= cyto_sent;
            total_cyto_transported += cyto_sent;
            if (logging) edge_cyto_down[i] += cyto_sent;
        }

        if (total_cyto_transported > 1e-12f && total_demand > 1e-12f) {
            // Distribute to sinks proportional to their share of received water.
            for (int i = 0; i < N; ++i) {
                Node& n = *flat[i].node;
                if (demand[i] <= 0.0f) continue;
                float water_received = demand[i] * demand_scale;
                float cyto_share = total_cyto_transported * (water_received / delivered);
                n.chemical(ChemicalID::Cytokinin) += cyto_share;
                if (logging) edge_cyto_up[i] += cyto_share;
            }
        }
    }

    // ── Conservation log ──────────────────────────────────────────────────────
    // One row per (node, chemical) with before/after and net change.  SUMMARY
    // row per (tick, chemical) with conservation_error = total_after − total_before.
    // With the running-budget caps above, conservation_error should stay within
    // float rounding (~1e-6 ml) every tick.
    if (logging) {
        std::filesystem::create_directories("debug");
        bool fresh_file = (world.current_tick == 0);
        auto mode = fresh_file
            ? (std::ios::out | std::ios::trunc)
            : (std::ios::out | std::ios::app);
        std::ofstream csv("debug/xylem_log.csv", mode);
        if (csv.is_open()) {
            if (fresh_file) {
                csv << "tick,chemical,node_id,node_type,parent_id,"
                       "amount_before,amount_after,edge_flow_up,edge_flow_down,"
                       "net_change\n";
            }

            auto write_rows = [&](ChemicalID chem, const std::vector<float>& pre,
                                  const std::vector<float>& eup,
                                  const std::vector<float>& edown) {
                float total_before = 0.0f, total_after = 0.0f;
                for (int i = 0; i < N; ++i) {
                    const Node& n = *flat[i].node;
                    float before = pre[i];
                    float after  = n.chemical(chem);
                    int parent_id = (flat[i].parent_idx >= 0)
                        ? static_cast<int>(flat[flat[i].parent_idx].node->id)
                        : -1;
                    csv << world.current_tick << ','
                        << chem_name(chem) << ','
                        << n.id << ','
                        << node_type_name(n.type) << ','
                        << parent_id << ','
                        << before << ',' << after << ','
                        << eup[i] << ',' << edown[i] << ','
                        << (after - before) << '\n';
                    total_before += before;
                    total_after  += after;
                }
                csv << world.current_tick << ',' << chem_name(chem) << ",SUMMARY,,"
                    << N << ','
                    << total_before << ',' << total_after << ",,,"
                    << (total_after - total_before) << '\n';
            };

            write_rows(ChemicalID::Water,     pre_water, edge_water_up, edge_water_down);
            write_rows(ChemicalID::Cytokinin, pre_cyto,  edge_cyto_up,  edge_cyto_down);
        }
    }
}

void vascular_transport(Plant& plant, const Genome& g, const WorldParams& world) {
    Node* seed = plant.seed_mut();
    if (!seed) return;

    // Phloem: Münch pressure-flow (Jacobi resolve with velocity + equalization caps).
    phloem_resolve(plant, g, world);

    // Xylem: water + cytokinin pressure-flow (mass-conserving Jacobi).
    xylem_resolve(plant, g, world);
}

} // namespace botany
