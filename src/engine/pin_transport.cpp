// src/engine/pin_transport.cpp — PIN-mediated polar auxin transport.
// Three phases per tick:
//   A. Shoot post-order: each shoot node pumps auxin toward its parent (basipetal).
//   B. Seed junction: collects from shoot children, distributes to root children by r².
//   C. Root pre-order: distributes auxin from seed toward root tips (acropetal).
// Runs after vascular_transport(), before the DFS tree walk.
#include "engine/pin_transport.h"
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/node/node.h"
#include "engine/chemical/chemical.h"
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>

namespace botany {

// --- Tree traversal helpers ---

static void collect_shoot_postorder(Node* n, std::vector<Node*>& out) {
    for (Node* child : n->children) {
        if (child->type != NodeType::ROOT && child->type != NodeType::ROOT_APICAL)
            collect_shoot_postorder(child, out);
    }
    out.push_back(n);
}

static void collect_root_preorder(Node* n, std::vector<Node*>& out) {
    out.push_back(n);
    for (Node* child : n->children) {
        if (child->type == NodeType::ROOT || child->type == NodeType::ROOT_APICAL)
            collect_root_preorder(child, out);
    }
}

// --- Main pass ---

void pin_transport(Plant& plant, const Genome& g) {
    Node* seed = plant.seed_mut();
    if (!seed) return;

    // Partition seed's children into shoot vs root sides.
    std::vector<Node*> shoot_children_of_seed;
    std::vector<Node*> root_children_of_seed;
    for (Node* c : seed->children) {
        if (c->type == NodeType::ROOT || c->type == NodeType::ROOT_APICAL)
            root_children_of_seed.push_back(c);
        else
            shoot_children_of_seed.push_back(c);
    }

    // ----------------------------------------------------------------
    // Phase A: Shoot post-order — leaves → seed.
    //   Each shoot node pumps auxin toward its parent.
    //   Anti-teleportation: use parent->transport_received for non-seed parents.
    //   Seed children use local accumulator (bug 2 fix — transport_received on seed
    //   isn't flushed until end of seed->tick(), so seed would distribute stale data).
    // ----------------------------------------------------------------
    std::vector<Node*> shoot_nodes;
    for (Node* c : shoot_children_of_seed)
        collect_shoot_postorder(c, shoot_nodes);

    float seed_collected = 0.0f;

    for (Node* n : shoot_nodes) {
        Node* par = n->parent;
        if (!par) continue;

        float r2 = n->radius * n->radius;
        float max_cap = r2 * g.pin_capacity_per_area;

        // Efficiency: cold-start floor + upregulation from canalization history.
        // Read last tick's auxin_flow_bias — PIN hasn't updated it yet this tick.
        float flow_bias = 0.0f;
        auto it = par->auxin_flow_bias.find(n);
        if (it != par->auxin_flow_bias.end()) flow_bias = it->second;
        float efficiency = g.pin_base_efficiency + flow_bias * (1.0f - g.pin_base_efficiency);

        float available = n->chemical(ChemicalID::Auxin);
        float moved = std::min(available, max_cap * efficiency);
        if (moved < 1e-8f) continue;

        n->chemical(ChemicalID::Auxin) -= moved;
        par->last_auxin_flux[n] += moved;  // record for update_canalization

        if (par == seed) {
            seed_collected += moved;
        } else {
            par->transport_received[ChemicalID::Auxin] += moved;
        }
    }

    // ----------------------------------------------------------------
    // Phase B: Seed junction — distribute this tick's shoot collection
    //   to root children weighted by radius (thicker root = larger PIN area).
    // ----------------------------------------------------------------
    std::unordered_map<Node*, float> root_forwarded;

    if (seed_collected > 1e-8f && !root_children_of_seed.empty()) {
        float total_r = 0.0f;
        for (Node* rc : root_children_of_seed) total_r += rc->radius;

        for (Node* rc : root_children_of_seed) {
            float share = (total_r > 1e-8f)
                ? rc->radius / total_r
                : 1.0f / static_cast<float>(root_children_of_seed.size());
            float max_cap = rc->radius * rc->radius * g.pin_capacity_per_area;
            float to_send = std::min(seed_collected * share, max_cap);
            root_forwarded[rc] = to_send;
            seed->last_auxin_flux[rc] += to_send;
            seed_collected -= to_send;
        }
    }
    seed->chemical(ChemicalID::Auxin) += seed_collected;  // remainder stays at seed

    // ----------------------------------------------------------------
    // Phase C: Root pre-order — seed → root tips.
    //   Each root node receives its forwarded auxin and passes the rest
    //   deeper. Writes directly to chemical(Auxin) so root update_tissue()
    //   sees PIN-delivered auxin this tick (PIN runs before the DFS walk).
    // ----------------------------------------------------------------
    std::vector<Node*> root_nodes;
    for (Node* c : root_children_of_seed)
        collect_root_preorder(c, root_nodes);

    for (Node* n : root_nodes) {
        auto it = root_forwarded.find(n);
        if (it == root_forwarded.end()) continue;
        float incoming = it->second;
        if (incoming < 1e-8f) continue;

        // Collect root children of this node.
        std::vector<Node*> root_children;
        float total_r = 0.0f;
        for (Node* child : n->children) {
            if (child->type == NodeType::ROOT || child->type == NodeType::ROOT_APICAL) {
                root_children.push_back(child);
                total_r += child->radius;
            }
        }

        // Distribute to root children; remainder stays at this node.
        float distributed = 0.0f;
        for (Node* rc : root_children) {
            float share = (total_r > 1e-8f)
                ? rc->radius / total_r
                : 1.0f / static_cast<float>(root_children.size());
            float max_cap = rc->radius * rc->radius * g.pin_capacity_per_area;
            float to_send = std::min(incoming * share, max_cap);
            root_forwarded[rc] += to_send;
            n->last_auxin_flux[rc] += to_send;
            distributed += to_send;
        }

        // incoming − distributed stays at this node.
        n->chemical(ChemicalID::Auxin) += (incoming - distributed);
    }
}

} // namespace botany
