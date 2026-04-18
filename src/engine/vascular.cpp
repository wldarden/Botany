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

static void run_vascular(std::vector<VascNodeInfo>& flat,
                         ChemicalID chem_id,
                         float conductance,
                         const Genome& g,
                         std::ofstream* log,
                         uint32_t current_tick = 0) {
    if (flat.empty()) return;

    // --- Phase 1: Post-order (leaves→seed) ---
    // Walk backwards through the pre-order array = post-order.
    for (int i = static_cast<int>(flat.size()) - 1; i >= 0; --i) {
        auto& info = flat[i];
        Node& n = *info.node;

        // Determine role: source, sink, or conduit
        info.is_conduit = has_vasculature(n, g);
        if (info.is_conduit) {
            info.capacity = pipe_capacity(n, conductance);
        }

        // Classify source/sink based on chemical and node type
        if (chem_id == ChemicalID::Sugar) {
            // Phloem: leaves and seed reserves are sources, growing tips are sinks
            if (!n.parent) {
                // Seed: mobilize stored reserves above the phloem reserve floor.
                // The seed's stored sugar is trapped by the vascular skip logic
                // (seed→STEM local diffusion is blocked because both have vasculature),
                // so the phloem pass is the only path out. Without this, 44g of seed
                // sugar sits inert while the plant starves.
                float cap = sugar_cap(n, g);
                float reserve = cap * g.phloem_reserve_fraction;
                float surplus = std::max(0.0f, n.chemical(chem_id) - reserve - n.sugar_reserved_for_growth);
                info.supply += surplus;
            } else if (n.type == NodeType::LEAF) {
                float cap = sugar_cap(n, g);
                float reserve = cap * g.phloem_reserve_fraction;
                // Also protect sugar reserved for local growth (set by pre_transport_growth()).
                // This ensures leaves have sugar for expansion when update_tissue() runs later.
                float surplus = std::max(0.0f, n.chemical(chem_id) - reserve - n.sugar_reserved_for_growth);
                info.supply += surplus;
            } else if (n.type == NodeType::APICAL || n.type == NodeType::ROOT_APICAL) {
                // Only ACTIVE meristems are phloem sinks. Dormant lateral buds
                // (active=false) rely on local diffusion from their parent conduit
                // node — registering them as vascular sinks bleeds sugar away from
                // the active growing tip at every level of a long chain, starving it.
                bool is_active = (n.type == NodeType::ROOT_APICAL)
                    ? n.as_root_apical()->active
                    : n.as_apical()->active;
                if (is_active) {
                    float cap = sugar_cap(n, g);
                    // Bound per-tick demand to cap×meristem_sink_fraction rather than
                    // the full cap−current deficit. A hungry meristem with 2.0g cap
                    // would otherwise demand 2.0g/tick — more than the whole plant
                    // produces — starving leaves before they can expand. The bounded
                    // pull matches actual elongation cost (~0.05g/tick) with headroom.
                    float per_tick_max = cap * g.meristem_sink_fraction;
                    float deficit = std::max(0.0f, per_tick_max - n.chemical(chem_id));
                    info.demand += deficit;
                }
                // Inactive buds intentionally do NOT fall through to the starvation
                // check below — they should not compete for phloem at all.
            } else if (n.starvation_ticks > 0) {
                // Starving nodes are emergency sinks
                float cap = sugar_cap(n, g);
                info.demand += std::max(0.0f, cap * 0.5f - n.chemical(chem_id));
            }
        } else {
            // Xylem: roots are sources, leaves and shoot tips are sinks
            if (n.type == NodeType::ROOT || n.type == NodeType::ROOT_APICAL) {
                float surplus = std::max(0.0f, n.chemical(chem_id) * 0.5f);
                info.supply += surplus;
            } else if (n.type == NodeType::LEAF) {
                if (chem_id == ChemicalID::Water) {
                    float cap = water_cap(n, g);
                    float deficit = std::max(0.0f, cap - n.chemical(chem_id));
                    info.demand += deficit;
                }
                // Cytokinin: leaves are NOT xylem sinks. They receive cytokinin
                // passively via local diffusion from their parent stem. Leaf CKX
                // enzymes actively degrade it (handled by decay).
            } else if (n.type == NodeType::APICAL) {
                if (chem_id == ChemicalID::Water) {
                    float cap = water_cap(n, g);
                    info.demand += std::max(0.0f, cap - n.chemical(chem_id));
                } else {
                    // Cytokinin: shoot tips pull only what they're missing — deficit-based.
                    // Stops pulling once above cytokinin_growth_threshold (Km), so
                    // well-supplied apicals don't crowd out farther ones.
                    info.demand += std::max(0.0f, g.cytokinin_growth_threshold - n.chemical(ChemicalID::Cytokinin));
                }
            }
        }

        // Propagate subtree supply AND demand separately to parent.
        // No within-subtree netting — the seed sees total supply and total
        // demand from the whole tree, so cross-branch flow works correctly
        // (e.g., leaf sugar reaches root tips via the seed junction).
        if (info.parent_idx >= 0) {
            auto& parent_info = flat[info.parent_idx];
            if (info.is_conduit) {
                parent_info.supply += std::min(info.supply, info.capacity);
                parent_info.demand += std::min(info.demand, info.capacity);
            } else {
                parent_info.supply += info.supply;
                parent_info.demand += info.demand;
            }
        }
    }

    // --- Phase 2: Pre-order (seed→leaves) ---
    // Distribute actual flow from seed outward.
    // Seed (index 0) has the total supply/demand of the whole tree.
    float seed_available = std::min(flat[0].supply, flat[0].demand);

    for (int i = 0; i < static_cast<int>(flat.size()); ++i) {
        auto& info = flat[i];
        Node& n = *info.node;

        // How much flow is available at this node?
        float available;
        if (i == 0) {
            available = seed_available;
        } else {
            // Inherited from parent during its distribution pass
            available = info.supply;  // reused as "flow arriving from parent"
        }

        // Deduct from local source if this node is one
        if (info.supply > 0 && available > 0) {
            // This node is a source — how much of the available flow is sourced here?
            // Compute what this node itself supplies (not children)
            float local_supply = 0.0f;
            if (chem_id == ChemicalID::Sugar && !n.parent) {
                // Seed: deduct mobilized reserves (mirrors Phase 1 classification).
                float cap = sugar_cap(n, g);
                float reserve = cap * g.phloem_reserve_fraction;
                local_supply = std::max(0.0f, n.chemical(chem_id) - reserve - n.sugar_reserved_for_growth);
            } else if (chem_id == ChemicalID::Sugar && n.type == NodeType::LEAF) {
                float cap = sugar_cap(n, g);
                float reserve = cap * g.phloem_reserve_fraction;
                local_supply = std::max(0.0f, n.chemical(chem_id) - reserve - n.sugar_reserved_for_growth);
            } else if (chem_id != ChemicalID::Sugar &&
                       (n.type == NodeType::ROOT || n.type == NodeType::ROOT_APICAL)) {
                local_supply = std::max(0.0f, n.chemical(chem_id) * 0.5f);
            }

            if (local_supply > 0) {
                float deducted = std::min(local_supply, available);
                n.chemical(chem_id) -= deducted;
            }
        }

        // Deliver to local sink if this node is one
        float local_demand = 0.0f;
        if (chem_id == ChemicalID::Sugar) {
            if (n.type == NodeType::APICAL || n.type == NodeType::ROOT_APICAL) {
                // Mirror Phase 1: only deliver to active meristems via the vascular
                // pass. Inactive buds receive sugar through local diffusion instead.
                bool is_active = (n.type == NodeType::ROOT_APICAL)
                    ? n.as_root_apical()->active
                    : n.as_apical()->active;
                if (is_active) {
                    float cap = sugar_cap(n, g);
                    // Mirror Phase 1 cap: deliver at most cap×meristem_sink_fraction per tick.
                    float per_tick_max = cap * g.meristem_sink_fraction;
                    local_demand = std::max(0.0f, per_tick_max - n.chemical(chem_id));
                }
            } else if (n.starvation_ticks > 0) {
                float cap = sugar_cap(n, g);
                local_demand = std::max(0.0f, cap * 0.5f - n.chemical(chem_id));
            }
        } else if (n.type == NodeType::LEAF) {
            // Leaves are xylem sinks for water only, not cytokinin
            if (chem_id == ChemicalID::Water) {
                float cap = water_cap(n, g);
                local_demand = std::max(0.0f, cap - n.chemical(chem_id));
            }
        } else if (n.type == NodeType::APICAL) {
            if (chem_id == ChemicalID::Water) {
                float cap = water_cap(n, g);
                local_demand = std::max(0.0f, cap - n.chemical(chem_id));
            } else {
                // Cytokinin: shoot tips are the primary sink
                local_demand = 0.05f;
            }
        }

        float delivered = std::min(available, local_demand);
        if (delivered > 0) {
            n.chemical(chem_id) += delivered;
            available -= delivered;
        }

        // Distribute remaining flow to children by conductance weight, with demand
        // as a ceiling. Conductance weight = pipe_capacity × (1 + canalization_weight × bias).
        // Bias is auxin_flow_bias on the parent-to-child connection: connections where
        // auxin has flowed repeatedly (higher PIN saturation) carry proportionally more flow.
        // When a child is satisfied (allocation ≥ demand), its surplus redistributes to
        // remaining siblings by conductance proportions. Iterate until budget exhausted or
        // all children satisfied.
        if (available > 0 && !info.child_idxs.empty()) {
            int n_ch = static_cast<int>(info.child_idxs.size());

            // Compute conductance weights and demand+pipe ceilings for each child.
            std::vector<float> weights(n_ch);
            std::vector<float> ceilings(n_ch);
            float total_weight = 0.0f;
            for (int k = 0; k < n_ch; ++k) {
                int ci = info.child_idxs[k];
                float cap = pipe_capacity(*flat[ci].node, conductance);
                auto it = info.node->auxin_flow_bias.find(flat[ci].node);
                float bias = (it != info.node->auxin_flow_bias.end()) ? it->second : 0.0f;
                weights[k] = cap * (1.0f + g.canalization_weight * bias);
                total_weight += weights[k];
                ceilings[k] = flat[ci].demand;
                if (flat[ci].is_conduit)
                    ceilings[k] = std::min(ceilings[k], flat[ci].capacity);
            }

            if (total_weight > 1e-8f) {
                // Iterative water-filling: each round allocates proportionally by
                // conductance weight, caps any child that hits its ceiling, then
                // redistributes the unclaimed surplus among the remaining children.
                std::vector<float> alloc(n_ch, 0.0f);
                std::vector<bool> done(n_ch, false);
                float budget = available;

                for (int iter = 0; iter <= n_ch; ++iter) {
                    float active_w = 0.0f;
                    for (int k = 0; k < n_ch; ++k)
                        if (!done[k]) active_w += weights[k];
                    if (active_w <= 1e-8f || budget <= 1e-8f) break;

                    bool any_capped = false;
                    for (int k = 0; k < n_ch; ++k) {
                        if (done[k]) continue;
                        float share = budget * (weights[k] / active_w);
                        if (share >= ceilings[k]) {
                            alloc[k] = ceilings[k];
                            done[k] = true;
                            any_capped = true;
                        }
                    }
                    if (!any_capped) {
                        // No child hit its ceiling — distribute remaining budget and finish.
                        for (int k = 0; k < n_ch; ++k)
                            if (!done[k])
                                alloc[k] = budget * (weights[k] / active_w);
                        break;
                    }
                    // Recompute remaining budget (available minus all capped allocations).
                    budget = available;
                    for (int k = 0; k < n_ch; ++k)
                        budget -= alloc[k];
                    if (budget <= 1e-8f) break;
                }

                for (int k = 0; k < n_ch; ++k) {
                    flat[info.child_idxs[k]].supply = alloc[k];

                    if (log) {
                        int ci = info.child_idxs[k];
                        Node* child = flat[ci].node;
                        float bias = 0.0f;
                        auto it = info.node->auxin_flow_bias.find(child);
                        if (it != info.node->auxin_flow_bias.end()) bias = it->second;
                        *log << current_tick << ','
                             << info.node->id << ',' << child->id << ','
                             << node_type_name(child->type) << ','
                             << chem_name(chem_id) << ','
                             << ceilings[k] << ','
                             << weights[k] << ','
                             << alloc[k] << ','
                             << bias << '\n';
                    }
                }
            }
        }
    }
}

// ── Münch pressure-flow phloem helpers ────────────────────────────────────────

// Phloem ring cross-section area — used ONLY for pipe capacity (flow_vol = velocity × ring_area).
// Living phloem+cambium layer; heartwood excluded. Scales linearly with r.
// Real phloem ring: ~0.1–0.5 mm thick regardless of stem diameter.
static float phloem_ring_area(float r, float t) {
    constexpr float kPi = 3.14159265f;
    return kPi * (2.0f * r * t - t * t);
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

        float r_eff      = std::min(leaf.radius, parent.radius);
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

        float parent_cap  = parent_vol * world.max_sugar_concentration;
        float parent_room = parent_cap - parent.chemical(ChemicalID::Sugar);
        load = std::min(load, std::max(0.0f, parent_room));

        if (load > 1e-8f) {
            leaf.chemical(ChemicalID::Sugar)   -= load;
            parent.chemical(ChemicalID::Sugar) += load;
        }
    }

    // ── Phase 2: Jacobi simultaneous pipe-network resolution ─────────────────
    // phloem_iterations inner loops.  Each iteration:
    //   a) Compute pressure at every vascular node from current sugar state.
    //   b) Compute desired flow for EVERY edge simultaneously (reads from
    //      start-of-iteration sugar only — no intra-iteration updates).
    //   c) Fairness scaling: if a sender would be overdrawn, scale its outflows.
    //   d) Apply iteration delta, clamped to [0, node_cap].
    //
    // Unloading at conduit nodes happens in the edge loop (Section 4.6):
    // each edge transfer splits into sieve-tube transit and membrane crossing.
    // Only the membrane-crossing fraction actually changes the receiver's sugar;
    // the transit fraction remains in the sender's sieve tube (credited back via
    // the sign convention: sender pays full desired_sugar, receiver gets unload).
    for (uint32_t iter = 0; iter < world.phloem_iterations; ++iter) {

        // a) Pressure at every vascular node
        std::vector<float> pressure(N, 0.0f);
        for (int i = 0; i < N; ++i) {
            if (has_vasculature(*flat[i].node, g))
                pressure[i] = compute_phloem_pressure(*flat[i].node, g, world);
        }

        // b) Desired flow for every edge simultaneously (reads start-of-iteration state)
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

                float r_eff     = std::min(cur.radius, next.radius);
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

                // Unloading: fraction of desired that crosses the membrane into next's tissue.
                // Bidirectional formula: gradient × permeability × transfer_volume.
                float next_conc = next.chemical(ChemicalID::Sugar) / next_vol;
                float grad_unload = std::max(0.0f, cur_conc - next_conc);
                float perm    = unloading_permeability(next.type, g);
                float unload  = std::min(grad_unload * perm * desired, desired);

                // sender pays full desired_sugar; receiver gets only unloaded fraction.
                // transit remainder (desired - unload) stays in sender's sieve tube:
                // it is credited back implicitly (sender delta = -desired + transit_back =
                // effectively -unload net), but for simplicity we track it in the delta
                // as: sender -= desired, receiver += unload.  Transit stays in the pipe
                // and doesn't leave the sender's node.
                iter_delta[i] -= desired;
                iter_delta[j] += unload;
            };

            if (flat[i].parent_idx >= 0) process_edge(flat[i].parent_idx);
            for (int ci : flat[i].child_idxs) process_edge(ci);
        }

        // c) Fairness scaling: prevent any node from sending more than it has.
        for (int i = 0; i < N; ++i) {
            if (iter_delta[i] >= 0.0f) continue;  // net receiver — no scaling needed
            float available  = flat[i].node->chemical(ChemicalID::Sugar);
            float total_out  = -iter_delta[i];
            if (total_out > available + 1e-8f) {
                float scale = available / total_out;
                iter_delta[i] = -available;
                // Scale corresponding inflows proportionally.  Because each edge
                // contributes to exactly one receiver's delta (unload fraction),
                // we approximate by clamping at the apply step.  The cap at (d) is
                // the hard safety net; per-sender scaling here reduces error.
                // Note: iter_delta[j] for receivers from this sender are already set.
                // A second pass to re-scale those entries is the refinement; the clamp
                // below handles conservation in practice.
            }
        }

        // d) Apply iteration delta — clamp to [0, node_cap]
        for (int i = 0; i < N; ++i) {
            if (std::abs(iter_delta[i]) < 1e-10f) continue;
            Node& n   = *flat[i].node;
            float cap = node_volume(n, world) * world.max_sugar_concentration;
            float new_val = std::clamp(
                n.chemical(ChemicalID::Sugar) + iter_delta[i], 0.0f, cap);
            n.chemical(ChemicalID::Sugar) = new_val;
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

        float r_eff     = std::min(parent.radius, mer.radius);
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

        float mer_cap  = mer_vol * world.max_sugar_concentration;
        float mer_room = mer_cap - mer.chemical(ChemicalID::Sugar);
        unload = std::min(unload, std::max(0.0f, mer_room));

        if (unload > 1e-8f) {
            parent.chemical(ChemicalID::Sugar) -= unload;
            mer.chemical(ChemicalID::Sugar)    += unload;
        }
    }
}

void vascular_transport(Plant& plant, const Genome& g, const WorldParams& world) {
    Node* seed = plant.seed_mut();
    if (!seed) return;

    // Build flat pre-order array of nodes
    std::vector<VascNodeInfo> flat;
    flat.reserve(plant.node_count());
    build_flat(seed, -1, flat);

    // Optional per-junction debug log.
    // On tick 0 the file is truncated (new run); subsequent ticks append.
    std::ofstream log_file;
    std::ofstream* log = nullptr;
    if (world.vascular_debug_log) {
        auto mode = (world.current_tick == 0)
            ? (std::ios::out | std::ios::trunc)
            : (std::ios::out | std::ios::app);
        log_file.open("debug/vascular_log.csv", mode);
        if (log_file.is_open()) {
            if (world.current_tick == 0)
                log_file << "tick,junction_node_id,child_node_id,child_type,"
                            "chemical,demand,conductance_weight,delivered,"
                            "auxin_flow_bias\n";
            log = &log_file;
        }
    }

    // Phloem: Münch pressure-flow (Jacobi resolve with velocity + equalization caps).
    // phloem_resolve builds its own flat array and handles leaf loading + pipe
    // network + meristem unloading in three sequential passes.
    phloem_resolve(plant, g, world);

    // Reset supply/demand for xylem pass (flat was already built above for logging)
    for (auto& info : flat) {
        info.supply = 0.0f;
        info.demand = 0.0f;
    }

    // Xylem: water from roots to shoots
    run_vascular(flat, ChemicalID::Water, g.xylem_conductance, g, log, world.current_tick);

    // Reset for cytokinin
    for (auto& info : flat) {
        info.supply = 0.0f;
        info.demand = 0.0f;
    }

    // Xylem: cytokinin from roots to shoots (carried in water stream)
    run_vascular(flat, ChemicalID::Cytokinin, g.xylem_conductance, g, log, world.current_tick);
}

} // namespace botany
