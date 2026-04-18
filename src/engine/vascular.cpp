// src/engine/vascular.cpp — Xylem & phloem vascular transport.
//
// phloem_resolve():  Münch pressure-driven sugar flow (see §4 of design doc).
//   Three sub-steps after DFS tick: leaf loading, BFS on conduit network,
//   meristem unloading. Sugar does NOT diffuse locally.
//
// xylem_resolve():  Phase 1/Phase 2 water+cytokinin transport (unchanged).
//   Post-order aggregates supply/demand; pre-order distributes via water-filling.
//   Runs after phloem_resolve so it sees final post-consumption water state.
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

        // Xylem: roots are sources, leaves and shoot tips are sinks.
        // (Sugar is handled by phloem_resolve, not run_vascular.)
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
            // passively via local diffusion from their parent stem.
        } else if (n.type == NodeType::APICAL) {
            if (chem_id == ChemicalID::Water) {
                float cap = water_cap(n, g);
                info.demand += std::max(0.0f, cap - n.chemical(chem_id));
            } else {
                // Cytokinin: shoot tips pull only what they're missing — deficit-based.
                info.demand += std::max(0.0f, g.cytokinin_growth_threshold - n.chemical(ChemicalID::Cytokinin));
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
            // Xylem only: roots are sources (sugar handled by phloem_resolve).
            if (n.type == NodeType::ROOT || n.type == NodeType::ROOT_APICAL) {
                local_supply = std::max(0.0f, n.chemical(chem_id) * 0.5f);
            }

            if (local_supply > 0) {
                float deducted = std::min(local_supply, available);
                n.chemical(chem_id) -= deducted;
            }
        }

        // Deliver to local sink if this node is one (xylem: water + cytokinin only)
        float local_demand = 0.0f;
        if (n.type == NodeType::LEAF) {
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

// ---------------------------------------------------------------------------
// phloem_resolve — Münch pressure-driven sugar transport
// ---------------------------------------------------------------------------
// Helpers — phloem geometry

static float phloem_ring_area(float r, float t) {
    constexpr float pi = 3.14159265f;
    return pi * (2.0f * r * t - t * t);
}

// Volume of the phloem compartment for a node (dm³):
//   STEM/ROOT: living phloem ring × length (heartwood excluded)
//   SEED (no parent): total cylinder volume — mobilisable starch storage organ
//   LEAF:      total mesophyll volume (all living tissue)
//   APICAL/ROOT_APICAL: sphere volume (entirely meristematic tissue)
static float phloem_volume(const Node& n, const WorldParams& world) {
    constexpr float pi = 3.14159265f;
    float t = world.phloem_ring_thickness;
    if (!n.parent) {
        // Seed — use total cylinder volume so stored reserves have meaningful concentration
        float length = glm::length(n.offset);
        if (length < 1e-6f) length = 0.01f;
        return pi * n.radius * n.radius * length;
    }
    switch (n.type) {
        case NodeType::STEM:
        case NodeType::ROOT: {
            float length = glm::length(n.offset);
            if (length < 1e-6f) return 0.0f;
            float area = phloem_ring_area(n.radius, t);
            return area * length;
        }
        case NodeType::LEAF: {
            if (const auto* leaf = n.as_leaf()) {
                float ls = leaf->leaf_size;
                return ls * ls * world.leaf_thickness;
            }
            return 0.0f;
        }
        case NodeType::APICAL:
        case NodeType::ROOT_APICAL: {
            return (4.0f / 3.0f) * pi * n.radius * n.radius * n.radius;
        }
        default:
            return 0.0f;
    }
}

static float compute_phloem_pressure(const Node& n, const Genome& g, const WorldParams& world) {
    float vol = phloem_volume(n, world);
    if (vol <= 0.0f) return 0.0f;
    float sugar_conc = n.chemical(ChemicalID::Sugar) / vol;
    float wc = water_cap(n, g);
    float water_frac = (wc > 0.0f)
        ? std::clamp(n.chemical(ChemicalID::Water) / wc, 0.0f, 1.0f)
        : 0.0f;
    return sugar_conc * g.phloem_osmotic_coefficient * water_frac;
}

static float phloem_local_speed(const Node& a, const Node& b, const WorldParams& w) {
    float r_eff = std::min(a.radius, b.radius);
    float r_ref = w.phloem_reference_radius;
    if (r_ref < 1e-8f) return 0.0f;
    return w.base_phloem_speed * (r_eff * r_eff) / (r_ref * r_ref);
}

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

void phloem_resolve(Plant& plant, const Genome& g, const WorldParams& world) {
    Node* seed = plant.seed_mut();
    if (!seed) return;

    // Build flat pre-order array (shared with xylem_resolve helpers)
    std::vector<VascNodeInfo> flat;
    flat.reserve(plant.node_count());
    build_flat(seed, -1, flat);
    int N = static_cast<int>(flat.size());
    if (N == 0) return;

    // --- Step 3a: Leaf loading pass (pre-BFS endpoint sources) ---
    // Each leaf loads sugar into its parent conduit ring proportional to the
    // concentration gradient × permeability × pipe capacity.
    // Capped at (a) available leaf sugar and (b) room in parent ring.
    for (int i = 0; i < N; ++i) {
        Node& n = *flat[i].node;
        if (n.type != NodeType::LEAF) continue;
        if (flat[i].parent_idx < 0) continue;

        Node& par = *flat[flat[i].parent_idx].node;
        float par_vol = phloem_volume(par, world);
        if (par_vol <= 0.0f) continue;

        float leaf_vol = phloem_volume(n, world);
        float leaf_conc = (leaf_vol > 0.0f) ? n.chemical(ChemicalID::Sugar) / leaf_vol : 0.0f;
        float par_conc = par.chemical(ChemicalID::Sugar) / par_vol;

        float gradient = std::max(0.0f, leaf_conc - par_conc);
        if (gradient <= 0.0f) continue;

        float r_eff = std::min(n.radius, par.radius);
        float pipe_cap = phloem_ring_area(r_eff, world.phloem_ring_thickness) * g.phloem_conductance;
        float load_raw = gradient * g.phloem_unloading_leaf * pipe_cap;

        // Cap 1: can't load more than the leaf has
        float load = std::min(load_raw, n.chemical(ChemicalID::Sugar));
        // Cap 2: can't overflow parent ring
        float par_room = par_vol * world.max_sugar_concentration - par.chemical(ChemicalID::Sugar);
        load = std::min(load, std::max(0.0f, par_room));

        if (load > 0.0f) {
            n.chemical(ChemicalID::Sugar) -= load;
            par.chemical(ChemicalID::Sugar) += load;
        }
    }

    // --- Step 3b: BFS on conduit network (stems, roots, seed) ---
    // Compute phloem pressure at every vascular conduit node.
    std::vector<float> pressure(N, 0.0f);
    for (int i = 0; i < N; ++i) {
        if (has_vasculature(*flat[i].node, g))
            pressure[i] = compute_phloem_pressure(*flat[i].node, g, world);
    }

    // Accumulate deltas. Shared across all source BFS walks so later walks see
    // earlier fills when computing available room (prevents over-cap).
    std::vector<float> delta(N, 0.0f);

    struct BFSEntry { int idx; float budget; float stream_conc; };
    std::vector<BFSEntry> sources;
    sources.reserve(N);

    for (int i = 0; i < N; ++i) {
        if (!has_vasculature(*flat[i].node, g)) continue;
        bool is_source = false;
        if (flat[i].parent_idx >= 0 && has_vasculature(*flat[flat[i].parent_idx].node, g))
            if (pressure[i] > pressure[flat[i].parent_idx]) is_source = true;
        for (int ci : flat[i].child_idxs)
            if (has_vasculature(*flat[ci].node, g) && pressure[i] > pressure[ci])
                is_source = true;
        if (is_source) {
            float vol = phloem_volume(*flat[i].node, world);
            float src_conc = (vol > 0.0f) ? flat[i].node->chemical(ChemicalID::Sugar) / vol : 0.0f;
            sources.push_back({i, 1.0f, src_conc});
        }
    }

    // Walk from each source outward (DFS stack simulates BFS budget walk).
    for (auto& entry : sources) {
        struct WalkState { int from_idx; int to_idx; float budget; float stream_conc; };
        std::vector<WalkState> walk;
        walk.push_back({-1, entry.idx, entry.budget, entry.stream_conc});

        while (!walk.empty()) {
            auto [from_i, cur_i, budget, stream_conc] = walk.back();
            walk.pop_back();

            Node& cur = *flat[cur_i].node;

            // Collect downhill vascular neighbours
            std::vector<int> nexts;
            if (flat[cur_i].parent_idx >= 0 && flat[cur_i].parent_idx != from_i
                    && has_vasculature(*flat[flat[cur_i].parent_idx].node, g)
                    && pressure[flat[cur_i].parent_idx] < pressure[cur_i])
                nexts.push_back(flat[cur_i].parent_idx);
            for (int ci : flat[cur_i].child_idxs)
                if (ci != from_i && has_vasculature(*flat[ci].node, g)
                        && pressure[ci] < pressure[cur_i])
                    nexts.push_back(ci);

            for (int nxt_i : nexts) {
                Node& nxt = *flat[nxt_i].node;

                float edge_len = glm::length(nxt.offset);
                float speed = phloem_local_speed(cur, nxt, world);
                float time_cost = (speed > 1e-8f) ? edge_len / speed : 1e30f;

                float time_fraction = std::min(1.0f, budget / time_cost);
                float remaining = budget - time_cost * time_fraction;

                float r_eff = std::min(cur.radius, nxt.radius);
                float pipe_cap = phloem_ring_area(r_eff, world.phloem_ring_thickness) * g.phloem_conductance;

                // Source capacity guard: flow from cur limited by its available sugar
                float src_available = cur.chemical(ChemicalID::Sugar) + delta[cur_i];
                float flow_vol = stream_conc * pipe_cap * time_fraction;
                flow_vol = std::min(flow_vol, std::max(0.0f, src_available));

                // Destination ring cap guard: don't overfill the ring
                float nxt_vol = phloem_volume(nxt, world);
                float nxt_cap = nxt_vol * world.max_sugar_concentration;
                float nxt_room = nxt_cap - (nxt.chemical(ChemicalID::Sugar) + delta[nxt_i]);
                flow_vol = std::min(flow_vol, std::max(0.0f, nxt_room));

                if (flow_vol <= 0.0f) continue;

                // Passive unloading
                float nxt_local_conc = (nxt_vol > 0.0f)
                    ? (nxt.chemical(ChemicalID::Sugar) + delta[nxt_i]) / nxt_vol : 0.0f;
                float gradient = std::max(0.0f, stream_conc - nxt_local_conc);
                float perm = unloading_permeability(nxt.type, g);
                float unload = std::min(gradient * perm * flow_vol, flow_vol);

                delta[nxt_i] += unload;
                delta[cur_i] -= flow_vol;

                float stream_after = (flow_vol > 1e-6f)
                    ? stream_conc * (1.0f - unload / flow_vol) : 0.0f;

                // Transit dead-end credit: when BFS cannot continue, credit
                // the in-pipe transit back to dest so conservation holds.
                bool will_continue = (remaining > 1e-4f) && (stream_after > 1e-8f);
                if (will_continue && !nexts.empty()) {
                    walk.push_back({cur_i, nxt_i, remaining, stream_after});
                } else {
                    delta[nxt_i] += (flow_vol - unload);
                }
            }
        }
    }

    // Apply deltas atomically (clamped to [0, cap])
    for (int i = 0; i < N; ++i) {
        if (delta[i] == 0.0f) continue;
        Node& n = *flat[i].node;
        float vol = phloem_volume(n, world);
        float cap = vol * world.max_sugar_concentration;
        float new_val = std::clamp(n.chemical(ChemicalID::Sugar) + delta[i], 0.0f, cap);
        n.chemical(ChemicalID::Sugar) = new_val;
    }

    // --- Step 3c: Meristem unloading pass (post-BFS terminal sinks) ---
    for (int i = 0; i < N; ++i) {
        Node& n = *flat[i].node;
        if (n.type != NodeType::APICAL && n.type != NodeType::ROOT_APICAL) continue;
        if (flat[i].parent_idx < 0) continue;

        Node& par = *flat[flat[i].parent_idx].node;
        float par_vol = phloem_volume(par, world);
        if (par_vol <= 0.0f) continue;

        float par_conc = par.chemical(ChemicalID::Sugar) / par_vol;
        float mer_vol = phloem_volume(n, world);
        float mer_conc = (mer_vol > 0.0f) ? n.chemical(ChemicalID::Sugar) / mer_vol : 0.0f;

        float gradient = std::max(0.0f, par_conc - mer_conc);
        if (gradient <= 0.0f) continue;

        float r_eff = std::min(n.radius, par.radius);
        float pipe_cap = phloem_ring_area(r_eff, world.phloem_ring_thickness) * g.phloem_conductance;
        float unload_raw = gradient * g.phloem_unloading_meristem * pipe_cap;

        // Cap: can't exceed parent sugar or meristem room
        float mer_room = std::max(0.0f, mer_vol * world.max_sugar_concentration - n.chemical(ChemicalID::Sugar));
        float unload = std::min({unload_raw, par.chemical(ChemicalID::Sugar), mer_room});

        if (unload > 0.0f) {
            par.chemical(ChemicalID::Sugar) -= unload;
            n.chemical(ChemicalID::Sugar) += unload;
        }
    }
}

// ---------------------------------------------------------------------------
// xylem_resolve — Phase 1/Phase 2 water + cytokinin transport (unchanged)
// ---------------------------------------------------------------------------
void xylem_resolve(Plant& plant, const Genome& g, const WorldParams& world) {
    Node* seed = plant.seed_mut();
    if (!seed) return;

    // Build flat pre-order array of nodes
    std::vector<VascNodeInfo> flat;
    flat.reserve(plant.node_count());
    build_flat(seed, -1, flat);

    // Optional per-junction debug log.
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
