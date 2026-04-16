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

namespace botany {

bool has_vasculature(const Node& n, const Genome& g) {
    if (!n.parent) return true;  // seed is always a vascular junction
    if (n.type == NodeType::STEM)
        return n.age >= g.cambium_maturation_ticks;
    if (n.type == NodeType::ROOT)
        return n.age >= g.root_cambium_maturation_ticks;
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
static void run_vascular(std::vector<VascNodeInfo>& flat,
                         ChemicalID chem_id,
                         float conductance,
                         const Genome& g) {
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
            // Phloem: leaves are sources, growing tips are sinks
            if (n.type == NodeType::LEAF) {
                float cap = sugar_cap(n, g);
                float reserve = cap * g.phloem_reserve_fraction;
                float surplus = std::max(0.0f, n.chemical(chem_id) - reserve);
                info.supply += surplus;
            } else if (n.type == NodeType::APICAL || n.type == NodeType::ROOT_APICAL) {
                // Meristems want sugar for growth — demand is their capacity deficit
                float cap = sugar_cap(n, g);
                float deficit = std::max(0.0f, cap - n.chemical(chem_id));
                info.demand += deficit;
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
                } else {
                    // Cytokinin: leaves don't strongly demand it
                    info.demand += 0.01f;
                }
            } else if (n.type == NodeType::APICAL) {
                if (chem_id == ChemicalID::Water) {
                    float cap = water_cap(n, g);
                    info.demand += std::max(0.0f, cap - n.chemical(chem_id));
                } else {
                    // Cytokinin: shoot tips need it for growth
                    info.demand += 0.05f;
                }
            }
        }

        // Propagate subtree supply/demand up to parent
        // Conduit nodes pass through; non-conduit nodes only contribute
        // their own supply/demand (children can't pipe through them).
        if (info.parent_idx >= 0) {
            auto& parent_info = flat[info.parent_idx];
            if (info.is_conduit) {
                // Conduit: pass subtree totals, bottlenecked by capacity
                float net_supply = std::max(0.0f, info.supply - info.demand);
                float net_demand = std::max(0.0f, info.demand - info.supply);
                parent_info.supply += std::min(net_supply, info.capacity);
                parent_info.demand += std::min(net_demand, info.capacity);
            } else {
                // Non-conduit (leaf, meristem, young node): contribute own supply/demand
                // but don't carry children's flow (they use local diffusion)
                float net_supply = std::max(0.0f, info.supply - info.demand);
                float net_demand = std::max(0.0f, info.demand - info.supply);
                parent_info.supply += net_supply;
                parent_info.demand += net_demand;
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
            if (chem_id == ChemicalID::Sugar && n.type == NodeType::LEAF) {
                float cap = sugar_cap(n, g);
                float reserve = cap * g.phloem_reserve_fraction;
                local_supply = std::max(0.0f, n.chemical(chem_id) - reserve);
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
                float cap = sugar_cap(n, g);
                local_demand = std::max(0.0f, cap - n.chemical(chem_id));
            } else if (n.starvation_ticks > 0) {
                float cap = sugar_cap(n, g);
                local_demand = std::max(0.0f, cap * 0.5f - n.chemical(chem_id));
            }
        } else if (n.type == NodeType::LEAF || n.type == NodeType::APICAL) {
            if (chem_id == ChemicalID::Water) {
                float cap = water_cap(n, g);
                local_demand = std::max(0.0f, cap - n.chemical(chem_id));
            } else {
                local_demand = (n.type == NodeType::APICAL) ? 0.05f : 0.01f;
            }
        }

        float delivered = std::min(available, local_demand);
        if (delivered > 0) {
            n.chemical(chem_id) += delivered;
            available -= delivered;
        }

        // Distribute remaining flow to children proportionally to their demand
        if (available > 0 && !info.child_idxs.empty()) {
            float total_child_demand = 0.0f;
            for (int ci : info.child_idxs) {
                total_child_demand += flat[ci].demand;
            }
            if (total_child_demand > 1e-8f) {
                for (int ci : info.child_idxs) {
                    float share = available * (flat[ci].demand / total_child_demand);
                    // Cap by child's pipe capacity
                    if (flat[ci].is_conduit) {
                        share = std::min(share, flat[ci].capacity);
                    }
                    // Reuse supply field to pass flow down to child
                    flat[ci].supply = share;
                }
            }
        }
    }
}

void vascular_transport(Plant& plant, const Genome& g) {
    Node* seed = plant.seed_mut();
    if (!seed) return;

    // Build flat pre-order array of nodes
    std::vector<VascNodeInfo> flat;
    flat.reserve(plant.node_count());
    build_flat(seed, -1, flat);

    // Phloem: sugar from leaves to sinks
    run_vascular(flat, ChemicalID::Sugar, g.phloem_conductance, g);

    // Reset supply/demand for xylem pass
    for (auto& info : flat) {
        info.supply = 0.0f;
        info.demand = 0.0f;
    }

    // Xylem: water from roots to shoots
    run_vascular(flat, ChemicalID::Water, g.xylem_conductance, g);

    // Reset for cytokinin
    for (auto& info : flat) {
        info.supply = 0.0f;
        info.demand = 0.0f;
    }

    // Xylem: cytokinin from roots to shoots (carried in water stream)
    run_vascular(flat, ChemicalID::Cytokinin, g.xylem_conductance, g);
}

} // namespace botany
