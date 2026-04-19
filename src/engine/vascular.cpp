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

    // ── Phase 2: Jacobi simultaneous pipe-network resolution ─────────────────
    // phloem_iterations inner loops.  Each iteration:
    //   a) Compute pressure at every vascular node from current sugar state.
    //   b) Compute desired flow for every edge processed in pre-order.  Pipe
    //      throughput (velocity × ring_area) and equalization both cap the
    //      transfer.  Unload permeability throttles how much of `desired`
    //      crosses the sieve-tube membrane into receiver tissue per tick.
    //   c) Apply iteration delta.
    //
    // Transit vs. storage — biologically correct interpretation:
    // The node's sugar chemical represents BOTH parenchymatic storage AND
    // sieve-tube transit sugar.  Real phloem sap runs at 150–250 g/dm³ which
    // far exceeds parenchyma-tissue storage concentration — transit sugar
    // can legitimately push a node's total above sugar_cap().  No receiver_
    // room cap here: a high-pressure push can fill an intermediate stem
    // above cap, and the next iteration's gradient forwards the excess.
    // Only the sender_remaining cap (sender can't export more than it has)
    // and the max_equil cap (no overshoot between adjacent nodes) constrain
    // per-edge flow.  Mass is conserved exactly every edge.
    for (uint32_t iter = 0; iter < world.phloem_iterations; ++iter) {

        // a) Pressure at every vascular node
        std::vector<float> pressure(N, 0.0f);
        for (int i = 0; i < N; ++i) {
            if (has_vasculature(*flat[i].node, g))
                pressure[i] = compute_phloem_pressure(*flat[i].node, g, world);
        }
        if (logging && iter == world.phloem_iterations - 1) {
            for (int i = 0; i < N; ++i) edge_final_pressure[i] = pressure[i];
        }

        std::vector<float> iter_delta(N, 0.0f);

        for (int i = 0; i < N; ++i) {
            if (!has_vasculature(*flat[i].node, g)) continue;
            Node& cur = *flat[i].node;
            float cur_vol = node_volume(cur, world);
            if (cur_vol <= 0.0f) continue;

            // Process each adjacent vascular edge where cur is the high-pressure sender
            auto process_edge = [&](int j) {
                if (!has_vasculature(*flat[j].node, g)) return;
                float dp = pressure[i] - pressure[j];
                if (dp <= 0.0f) return;  // cur is the high-pressure side

                Node& next     = *flat[j].node;
                float next_vol = node_volume(next, world);
                if (next_vol <= 0.0f) return;

                float r_eff     = edge_pipe_radius(cur, next);
                float ring_area = phloem_ring_area(r_eff, world.phloem_ring_thickness);
                float velocity  = std::min(dp * g.conductance_per_pressure,
                                           world.max_phloem_velocity);
                float flow_vol  = velocity * ring_area;  // dm³/tick

                float cur_conc = cur.chemical(ChemicalID::Sugar) / cur_vol;  // g/dm³
                float max_vel  = flow_vol * cur_conc;  // g — velocity-limited sugar

                // Equalization: cur cannot drop below the concentration that would
                // result if both nodes pooled their sugar.  Same principle as the
                // diffusion equalization clamp in Node::compute_transport_flow().
                float equil_conc = (cur.chemical(ChemicalID::Sugar) + next.chemical(ChemicalID::Sugar))
                                   / std::max(cur_vol + next_vol, 1e-8f);
                float max_equil  = std::max(0.0f,
                                       cur.chemical(ChemicalID::Sugar) - equil_conc * cur_vol);

                float desired = std::min(max_vel, max_equil);
                desired = std::min(desired, cur.chemical(ChemicalID::Sugar));  // safety

                // MASS-CONSERVATION: cap desired by sender's RUNNING sugar budget
                // this iteration — what cur_sugar has left after prior outflows
                // scheduled earlier in this iteration (iter_delta[i] is negative
                // when cur is a net sender).
                float sender_remaining = cur.chemical(ChemicalID::Sugar) + iter_delta[i];
                desired = std::min(desired, std::max(0.0f, sender_remaining));
                if (desired <= 0.0f) return;

                // Unloading: fraction of desired that crosses the membrane into next's tissue.
                // Bidirectional formula: gradient × permeability × transfer_volume.
                // Membrane permeability is a THROTTLE on the edge — it caps how much
                // sugar actually transfers this iteration. It MUST NOT be a leak.
                float next_conc = next.chemical(ChemicalID::Sugar) / next_vol;
                float grad_unload = std::max(0.0f, cur_conc - next_conc);
                float perm    = unloading_permeability(next.type, g);
                float unload  = std::min(grad_unload * perm * desired, desired);
                if (unload <= 0.0f) return;

                // NO receiver_room cap here — in real phloem the sieve tube
                // carries sugar through nodes at concentrations far above the
                // parenchymatic storage cap (sieve sap is ~150–250 g/dm³ sucrose
                // in a stem whose bulk-tissue concentration is much lower).  The
                // node's sugar pool represents both parenchymatic storage AND
                // transit sugar in the sieve tube — it is allowed to temporarily
                // exceed sugar_cap() when a high-pressure push is in transit.
                // The next Jacobi iteration will see the elevated concentration
                // and forward it along the pressure gradient.  Equilibration
                // (max_equil clamp above) still prevents overshoot between any
                // two adjacent nodes; the sender_remaining cap still prevents a
                // sender from pushing more than it holds.  Conservation is the
                // exact per-edge invariant sender -= X; receiver += X.

                // Mass-conserving edge transfer: sender loses exactly what
                // receiver gains.  Verified every tick by SUMMARY conservation
                // row in phloem_log.csv.
                iter_delta[i] -= unload;
                iter_delta[j] += unload;
                if (logging) {
                    log_flow_out[i] += unload;
                    log_flow_in[j]  += unload;
                    // Per-edge tracking: identify edge by whether j is parent of i
                    // or a child. Sums across all Jacobi iterations.  `desired` is
                    // logged for diagnostics — it is the pressure-demanded flow
                    // before membrane throttling; `unload` is what actually moved.
                    if (flat[i].parent_idx == j) {
                        edge_desired_to_parent[i] += desired;
                        edge_unload_to_parent[i]  += unload;
                    } else {
                        // j is a child of i. Store totals on the child's row
                        // so each edge owns one row.
                        edge_desired_from_parent[j] += desired;
                        edge_unload_from_parent[j]  += unload;
                    }
                }
            };

            if (flat[i].parent_idx >= 0) process_edge(flat[i].parent_idx);
            for (int ci : flat[i].child_idxs) process_edge(ci);
        }

        // c) Apply iteration delta.  Sender cap (sender_remaining in process_edge)
        // keeps each node's sugar ≥ 0; the max(,0) is a float-rounding safety rail.
        // No upper cap: intermediate stems may temporarily hold sugar above their
        // parenchymatic sugar_cap() when sieve-tube transit is in progress.  That
        // elevated concentration feeds the pressure gradient for the next
        // iteration, which forwards the excess toward sinks.
        for (int i = 0; i < N; ++i) {
            if (std::abs(iter_delta[i]) < 1e-10f) continue;
            Node& n   = *flat[i].node;
            float new_val = n.chemical(ChemicalID::Sugar) + iter_delta[i];
            n.chemical(ChemicalID::Sugar) = std::max(new_val, 0.0f);
        }
    }

    // ── Phase 3: Meristem unloading pass ─────────────────────────────────────
    // Each meristem draws sugar from its parent vascular conduit.  Velocity-capped
    // + equalization-capped; meristem room cap prevents over-saturation.
    for (int i = 0; i < N; ++i) {
        Node& mer = *flat[i].node;
        if (mer.type != NodeType::APICAL && mer.type != NodeType::ROOT_APICAL) continue;
        if (flat[i].parent_idx < 0) continue;
        Node& parent = *flat[flat[i].parent_idx].node;
        if (!has_vasculature(parent, g)) continue;

        float parent_vol = node_volume(parent, world);
        float mer_vol    = node_volume(mer, world);
        if (parent_vol <= 0.0f || mer_vol <= 0.0f) continue;

        float parent_conc = parent.chemical(ChemicalID::Sugar) / parent_vol;
        float mer_conc    = mer.chemical(ChemicalID::Sugar) / mer_vol;
        float gradient    = parent_conc - mer_conc;
        if (gradient <= 0.0f) continue;  // parent not above meristem — no unloading

        float r_eff     = edge_pipe_radius(parent, mer);
        float ring_area = phloem_ring_area(r_eff, world.phloem_ring_thickness);
        float velocity  = std::min(gradient * g.conductance_per_pressure,
                                   world.max_phloem_velocity);
        float flow_vol  = velocity * ring_area;
        float max_vel   = flow_vol * parent_conc;

        float equil_conc = (parent.chemical(ChemicalID::Sugar) + mer.chemical(ChemicalID::Sugar))
                           / std::max(parent_vol + mer_vol, 1e-8f);
        float max_equil  = std::max(0.0f,
                               parent.chemical(ChemicalID::Sugar) - equil_conc * parent_vol);

        float unload = std::min(max_vel, max_equil);
        unload = std::min(unload, parent.chemical(ChemicalID::Sugar));  // safety

        // Meristems and leaves don't have their own distributed vasculature —
        // they tap the parent conduit's phloem.  Storage cap is the per-type
        // sugar_cap() (g.sugar_cap_meristem = 0.1 g for meristems).  Draining
        // ~0.005 g/tick from the parent is enough to fuel normal elongation.
        float mer_cap  = sugar_cap(mer, g);
        float mer_room = mer_cap - mer.chemical(ChemicalID::Sugar);
        unload = std::min(unload, std::max(0.0f, mer_room));

        if (unload > 1e-8f) {
            parent.chemical(ChemicalID::Sugar) -= unload;
            mer.chemical(ChemicalID::Sugar)    += unload;
            if (logging) log_unloaded[flat[i].parent_idx] += unload;
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

// ── xylem_resolve: mass-conserving Jacobi pressure-flow for water + cytokinin ──
//
// Mirrors phloem_resolve in structure: flat pre-order traversal, Jacobi outer
// iterations, per-edge mass-conserving transfer with running-budget caps on both
// sender and receiver.  Key differences from phloem_resolve:
//
//   • Driving "pressure" is water_fraction = water / water_cap (dimensionless, 0–1).
//     No osmotic-coefficient multiplier — dp = frac_sender − frac_receiver directly
//     drives flow.  This is a pragmatic proxy for xylem water potential: the
//     thirstiest node pulls the hardest.  Leaves transpire → frac drops → they
//     pull from stems → stems pull from roots → roots absorb from soil.
//
//   • Pipe area = full node cross-section (π·r²) rather than a thin ring.
//     In real stems the xylem (dead wood vessels) occupies the bulk of the cross
//     section, whereas phloem is a thin living layer.  So xylem throughput
//     scales with r² instead of r.
//
//   • Two separate passes: water first (driven by its own fraction gradient),
//     then cytokinin (driven by the same water gradient; amount moved per edge
//     = flow_vol × sender_cyto_per_water).  Cyto rides the xylem stream.
//
//   • Leaves, meristems, and all other nodes participate — there is no Phase-1
//     loading or Phase-3 unloading like phloem.  The whole network is one pipe.
//     Transpiration at leaves is modeled elsewhere (LeafNode::transpire) as a
//     direct subtraction from leaf water; that creates the gradient this pass
//     resolves.
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
    std::vector<float> edge_water_up(N, 0.0f);   // water flowing from i to parent(i)
    std::vector<float> edge_water_down(N, 0.0f); // water flowing from i to children (+ sum over children)
    std::vector<float> edge_cyto_up(N, 0.0f);
    std::vector<float> edge_cyto_down(N, 0.0f);
    if (logging) {
        for (int i = 0; i < N; ++i) {
            pre_water[i] = flat[i].node->chemical(ChemicalID::Water);
            pre_cyto[i]  = flat[i].node->chemical(ChemicalID::Cytokinin);
        }
    }

    // ── Pre-compute water caps (same for both passes; unchanged across iters) ──
    std::vector<float> wcap(N, 0.0f);
    for (int i = 0; i < N; ++i) wcap[i] = water_cap(*flat[i].node, g);

    // ── Water pass ─────────────────────────────────────────────────────────────
    for (uint32_t iter = 0; iter < world.xylem_iterations; ++iter) {
        // Compute water_fraction at every node (the "pressure" proxy).
        std::vector<float> wfrac(N, 0.0f);
        for (int i = 0; i < N; ++i) {
            if (wcap[i] > 0.0f)
                wfrac[i] = flat[i].node->chemical(ChemicalID::Water) / wcap[i];
        }

        std::vector<float> iter_delta(N, 0.0f);

        for (int i = 0; i < N; ++i) {
            Node& cur = *flat[i].node;

            auto process_edge = [&](int j) {
                Node& next = *flat[j].node;
                float dp = wfrac[i] - wfrac[j];
                if (dp <= 1e-6f) return;  // cur must be the fuller (higher-frac) side

                float r_eff    = edge_pipe_radius(cur, next);
                float pipe_area = 3.14159265f * r_eff * r_eff;  // full cross-section
                float velocity  = std::min(dp * g.xylem_conductance_per_pressure,
                                           world.max_xylem_velocity);
                float flow_vol  = velocity * pipe_area;  // dm³/tick

                // Equalization cap — cur cannot drop below the fraction that
                // would result if both nodes pooled their water.  Prevents
                // oscillation from overshoot.
                float total_water = cur.chemical(ChemicalID::Water) + next.chemical(ChemicalID::Water);
                float total_cap   = wcap[i] + wcap[j];
                float equil_frac  = (total_cap > 1e-8f) ? total_water / total_cap : 0.0f;
                float max_equil   = std::max(0.0f,
                    cur.chemical(ChemicalID::Water) - equil_frac * wcap[i]);

                float desired = std::min(flow_vol, max_equil);

                // Sender cap only — the xylem equivalent of the "sieve-tube
                // transit" relaxation in phloem.  Real xylem vessels carry
                // water through a node at high throughput; the node's water
                // chemical represents bulk tissue water AND vessel transit,
                // and is allowed to temporarily exceed water_cap in transit.
                // The next iteration sees the elevated water_frac and pushes
                // the excess forward along the gradient.  Equalization cap
                // (max_equil) prevents overshoot between adjacent nodes.
                float sender_remaining = cur.chemical(ChemicalID::Water) + iter_delta[i];
                desired = std::min(desired, std::max(0.0f, sender_remaining));
                if (desired <= 0.0f) return;

                iter_delta[i] -= desired;
                iter_delta[j] += desired;
                if (logging) {
                    if (flat[i].parent_idx == j) edge_water_up[i] += desired;
                    else                         edge_water_down[i] += desired;
                }
            };

            if (flat[i].parent_idx >= 0) process_edge(flat[i].parent_idx);
            for (int ci : flat[i].child_idxs) process_edge(ci);
        }

        // Apply — sender cap still guarantees ≥ 0 so the max(,0) is a
        // float-rounding safety rail (cannot fire for >~1e-6 ml).  No upper
        // cap: water_cap bounds physical storage, but the chemical pool may
        // temporarily hold extra water in transit through the xylem.  That
        // transit excess drives the gradient for the next iteration.
        for (int i = 0; i < N; ++i) {
            if (std::abs(iter_delta[i]) < 1e-10f) continue;
            Node& n = *flat[i].node;
            n.chemical(ChemicalID::Water) = std::max(n.chemical(ChemicalID::Water) + iter_delta[i], 0.0f);
        }
    }

    // ── Cytokinin pass ────────────────────────────────────────────────────────
    // Cytokinin rides in the water stream.  Drive = same water_fraction gradient
    // (cyto goes wherever water goes).  Mass per edge = flow_vol × sender's
    // cyto-per-water concentration.  Mass-conserving via the same running caps.
    for (uint32_t iter = 0; iter < world.xylem_iterations; ++iter) {
        std::vector<float> wfrac(N, 0.0f);
        for (int i = 0; i < N; ++i) {
            if (wcap[i] > 0.0f)
                wfrac[i] = flat[i].node->chemical(ChemicalID::Water) / wcap[i];
        }

        std::vector<float> iter_delta(N, 0.0f);

        for (int i = 0; i < N; ++i) {
            Node& cur = *flat[i].node;

            auto process_edge = [&](int j) {
                Node& next = *flat[j].node;
                float dp = wfrac[i] - wfrac[j];
                if (dp <= 1e-6f) return;

                float cur_water = cur.chemical(ChemicalID::Water);
                if (cur_water <= 1e-8f) return;
                float cyto_per_water = cur.chemical(ChemicalID::Cytokinin) / cur_water;
                if (cyto_per_water <= 0.0f) return;

                float r_eff     = edge_pipe_radius(cur, next);
                float pipe_area = 3.14159265f * r_eff * r_eff;
                float velocity  = std::min(dp * g.xylem_conductance_per_pressure,
                                           world.max_xylem_velocity);
                float flow_vol  = velocity * pipe_area;  // dm³/tick

                float desired = flow_vol * cyto_per_water;  // mass of cyto carried

                float sender_remaining = cur.chemical(ChemicalID::Cytokinin) + iter_delta[i];
                desired = std::min(desired, std::max(0.0f, sender_remaining));
                if (desired <= 0.0f) return;

                // Cytokinin has no storage cap (signal molecule) — receiver accepts
                // whatever arrives.  Decay (Node::decay_chemicals) prevents
                // unbounded accumulation.
                iter_delta[i] -= desired;
                iter_delta[j] += desired;
                if (logging) {
                    if (flat[i].parent_idx == j) edge_cyto_up[i] += desired;
                    else                         edge_cyto_down[i] += desired;
                }
            };

            if (flat[i].parent_idx >= 0) process_edge(flat[i].parent_idx);
            for (int ci : flat[i].child_idxs) process_edge(ci);
        }

        for (int i = 0; i < N; ++i) {
            if (std::abs(iter_delta[i]) < 1e-10f) continue;
            Node& n = *flat[i].node;
            n.chemical(ChemicalID::Cytokinin) = std::max(
                n.chemical(ChemicalID::Cytokinin) + iter_delta[i], 0.0f);
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
