// src/engine/plant.cpp
#include "engine/plant.h"
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"
#include "engine/ethylene.h"
#include "engine/vascular_sub_stepped.h"
#include "engine/pin_transport.h"
#include "engine/perf_log.h"
#include "engine/world_params.h"
#include "engine/sugar.h"
#include <algorithm>
#include <unordered_set>

namespace botany {

Plant::Plant(const Genome& genome, glm::vec3 position)
    : genome_(genome)
{
    // Seed node
    Node* seed = create_node(NodeType::STEM, position, genome.initial_radius);
    seed->local().chemical(ChemicalID::Sugar) = genome.seed_sugar;

    // Bootstrap cytokinin: the embryo contains hormones throughout,
    // enough to start growth until roots produce their own.
    float seed_cyt = genome.cytokinin_growth_threshold * 5.0f;
    seed->local().chemical(ChemicalID::Cytokinin) = seed_cyt;

    // Bootstrap auxin: the embryo also contains stored IAA. This gives the
    // primary root its first growth impulse before shoot auxin can diffuse down,
    // which takes at least one tick (seed ticks before shoot in DFS order).
    float seed_auxin = genome.root_auxin_growth_threshold * 2.0f;
    seed->local().chemical(ChemicalID::Auxin) = seed_auxin;

    // Shoot apical meristem node (child of seed) — flagged as primary so it
    // cannot quiesce.  Real plants maintain their main apex continuously; only
    // lateral buds dormant.  If this meristem dies, Plant::tick_tree re-promotes
    // the strongest lateral (highest canalization path) to primary.
    Node* shoot = create_node(NodeType::APICAL, glm::vec3(0.0f, 0.01f, 0.0f), genome.initial_radius);
    shoot->local().chemical(ChemicalID::Cytokinin) = seed_cyt;
    shoot->local().chemical(ChemicalID::Auxin) = seed_auxin;
    shoot->as_apical()->is_primary = true;
    seed->add_child(shoot);

    // Root apical meristem node (child of seed) — same primary flag.
    Node* root = create_node(NodeType::ROOT_APICAL, glm::vec3(0.0f, -0.01f, 0.0f), genome.root_initial_radius);
    root->local().chemical(ChemicalID::Cytokinin) = seed_cyt;
    root->as_root_apical()->is_primary = true;
    seed->add_child(root);

    // Seed water: initialize at capacity (computed after children are known).
    // Using seed_sugar here would be wrong — it's 48 g-glucose, far over the ~2 ml
    // water capacity of a seedling, which causes the seed to flood children and
    // masks root-to-leaf flow in the visualizer.
    seed->local().chemical(ChemicalID::Water) = water_cap(*seed, genome);

    // Set initial world positions (subsequent ticks compute inline during DFS)
    seed->position = seed->offset;
    shoot->position = seed->position + shoot->offset;
    root->position = seed->position + root->offset;
}

void Plant::tick(const WorldParams& world, PerfStats* perf) {
    perf_ = perf;
    tick_tree(world, perf);
    perf_ = nullptr;
}

void Plant::remove_subtree(Node* node) {
    if (!node) return;

    // Remove from parent's children list
    if (node->parent) {
        auto& siblings = node->parent->children;
        siblings.erase(
            std::remove(siblings.begin(), siblings.end(), node),
            siblings.end()
        );
    }

    // Collect all nodes in the subtree (BFS)
    std::vector<Node*> to_remove;
    std::vector<Node*> queue = {node};
    while (!queue.empty()) {
        Node* n = queue.back();
        queue.pop_back();
        to_remove.push_back(n);
        for (Node* child : n->children) {
            queue.push_back(child);
        }
    }

    // Update root meristem count
    for (Node* n : to_remove) {
        if (n->type == NodeType::ROOT_APICAL) {
            if (root_meristem_count_ > 0) root_meristem_count_--;
        }
    }

    // Build a set for fast lookup
    std::unordered_set<Node*> remove_set(to_remove.begin(), to_remove.end());

    // Remove nodes
    nodes_.erase(
        std::remove_if(nodes_.begin(), nodes_.end(),
            [&](const std::unique_ptr<Node>& n) {
                return remove_set.count(n.get()) > 0;
            }),
        nodes_.end()
    );
}

void Plant::queue_removal(Node* node) {
    pending_removals_.push_back(node);
}

void Plant::flush_removals() {
    if (pending_removals_.empty()) return;

    std::unordered_set<Node*> remove_set(pending_removals_.begin(), pending_removals_.end());
    for (Node* n : pending_removals_) {
        if (n->type == NodeType::ROOT_APICAL) {
            if (root_meristem_count_ > 0) root_meristem_count_--;
        }
    }

    nodes_.erase(
        std::remove_if(nodes_.begin(), nodes_.end(),
            [&](const std::unique_ptr<Node>& n) {
                return remove_set.count(n.get()) > 0;
            }),
        nodes_.end()
    );
    pending_removals_.clear();
}

static void tick_recursive(Node& node, Plant& plant, const WorldParams& world) {
    node.tick(plant, world);

    // Snapshot children: meristem ticks may reparent or add siblings
    auto children = node.children;
    for (Node* child : children) {
        tick_recursive(*child, plant, world);
    }
}

void Plant::tick_tree(const WorldParams& world, PerfStats* /*perf*/) {
    // Clear last_auxin_flux on every node before anything writes to it this tick.
    // PIN writes first, then per-node diffusion accumulates during the DFS walk,
    // and update_canalization reads. The post-tick logger in engine.cpp reads the
    // same map, so we must NOT clear in update_canalization — clear once here.
    // Also reset per-tick hormone production counters (diagnostic only).
    for (auto& n : nodes_) {
        n->last_auxin_flux.clear();
        n->tick_auxin_produced = 0.0f;
        n->tick_cytokinin_produced = 0.0f;
    }

    // Reset the per-tick primary-meristem trackers.  Meristems set these
    // during their own update_tissue (piggyback on the DFS walk below).
    primary_sa_id_this_tick = -1;
    primary_ra_id_this_tick = -1;

    pin_transport(*this, genome_);              // PIN-mediated polar auxin transport (before metabolism)
    tick_recursive(*nodes_[0], *this, world);   // Phase 1: per-node metabolism
    vascular_sub_stepped(*this, genome_, world); // Phase 2: vascular transport (after metabolism)
    flush_removals();

    // After the DFS walk + removals: if either primary meristem's tracker is
    // still -1, the primary died (or never existed).  Promote the best
    // surviving lateral to primary based on canalization dominance.
    if (primary_sa_id_this_tick < 0 || primary_ra_id_this_tick < 0) {
        promote_primary_meristems();
    }
}

// Walk the chain of parent→child edges from node back to root, summing the
// auxin_flow_bias values (each stored on the parent, keyed by child pointer).
// Represents how dominantly-canalized the path to this node has been —
// higher sum = stronger historical auxin flux = more deserving of becoming
// the new primary meristem.  A zero-flux bud gets zero; the dominant axis
// accumulates over each step of its ancestry.
static float canalization_score(const Node* n) {
    float sum = 0.0f;
    for (const Node* cur = n; cur && cur->parent; cur = cur->parent) {
        auto it = cur->parent->auxin_flow_bias.find(const_cast<Node*>(cur));
        if (it != cur->parent->auxin_flow_bias.end()) sum += it->second;
    }
    return sum;
}

// Called at end of tick_tree when a primary tracker is -1.  Walks all
// nodes once (rare event), picks the surviving non-primary meristem with
// the highest canalization_score, and promotes it: sets is_primary=true
// and ensures it's active.  If no candidate of that type exists, the plant
// has genuinely lost that subsystem (shoot or root entirely dead) — a real
// game-over state for that half of the plant.
void Plant::promote_primary_meristems() {
    if (primary_sa_id_this_tick < 0) {
        ApicalNode* best = nullptr;
        float best_score = -1.0f;
        for (auto& n : nodes_) {
            ApicalNode* sa = n->as_apical();
            if (!sa) continue;
            float score = canalization_score(n.get());
            if (score > best_score) {
                best_score = score;
                best = sa;
            }
        }
        if (best) {
            best->is_primary = true;
            best->active = true;
            best->starvation_ticks = 0;
        }
    }
    if (primary_ra_id_this_tick < 0) {
        RootApicalNode* best = nullptr;
        float best_score = -1.0f;
        for (auto& n : nodes_) {
            RootApicalNode* ra = n->as_root_apical();
            if (!ra) continue;
            float score = canalization_score(n.get());
            if (score > best_score) {
                best_score = score;
                best = ra;
            }
        }
        if (best) {
            best->is_primary = true;
            best->active = true;
            best->starvation_ticks = 0;
        }
    }
}


Node* Plant::create_node(NodeType type, glm::vec3 position, float radius) {
    std::unique_ptr<Node> node;
    uint32_t id = next_id();
    switch (type) {
        case NodeType::STEM: node = std::make_unique<StemNode>(id, position, radius); break;
        case NodeType::ROOT: node = std::make_unique<RootNode>(id, position, radius); break;
        case NodeType::LEAF: node = std::make_unique<LeafNode>(id, position, radius); break;
        case NodeType::APICAL: node = std::make_unique<ApicalNode>(id, position, radius); break;
        case NodeType::ROOT_APICAL:
            node = std::make_unique<RootApicalNode>(id, position, radius);
            root_meristem_count_++;
            break;
    }
    Node* ptr = node.get();
    nodes_.push_back(std::move(node));
    return ptr;
}

void Plant::for_each_node(std::function<void(const Node&)> fn) const {
    for (const auto& node : nodes_) {
        fn(*node);
    }
}

void Plant::for_each_node_mut(std::function<void(Node&)> fn) {
    for (auto& node : nodes_) {
        fn(*node);
    }
}

} // namespace botany
