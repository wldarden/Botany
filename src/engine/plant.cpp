// src/engine/plant.cpp
#include "engine/plant.h"
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"
#include "engine/ethylene.h"
#include "engine/vascular.h"
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
    seed->chemical(ChemicalID::Sugar) = genome.seed_sugar;

    // Bootstrap cytokinin: the embryo contains hormones throughout,
    // enough to start growth until roots produce their own.
    float seed_cyt = genome.cytokinin_growth_threshold * 5.0f;
    seed->chemical(ChemicalID::Cytokinin) = seed_cyt;

    // Bootstrap auxin: the embryo also contains stored IAA. This gives the
    // primary root its first growth impulse before shoot auxin can diffuse down,
    // which takes at least one tick (seed ticks before shoot in DFS order).
    float seed_auxin = genome.root_auxin_growth_threshold * 2.0f;
    seed->chemical(ChemicalID::Auxin) = seed_auxin;

    // Shoot apical meristem node (child of seed)
    Node* shoot = create_node(NodeType::APICAL, glm::vec3(0.0f, 0.01f, 0.0f), genome.initial_radius);
    shoot->chemical(ChemicalID::Cytokinin) = seed_cyt;
    shoot->chemical(ChemicalID::Auxin) = seed_auxin;
    seed->add_child(shoot);

    // Root apical meristem node (child of seed)
    Node* root = create_node(NodeType::ROOT_APICAL, glm::vec3(0.0f, -0.01f, 0.0f), genome.root_initial_radius);
    root->chemical(ChemicalID::Cytokinin) = seed_cyt;
    seed->add_child(root);

    // Seed water: initialize at capacity (computed after children are known).
    // Using seed_sugar here would be wrong — it's 48 g-glucose, far over the ~2 ml
    // water capacity of a seedling, which causes the seed to flood children and
    // masks root-to-leaf flow in the visualizer.
    seed->chemical(ChemicalID::Water) = water_cap(*seed, genome);

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

void Plant::pre_transport_growth(const WorldParams& world) {
    // Each node computes how much sugar it will spend on growth this tick.
    // This runs before vascular_transport() so the phloem pass sees
    // (sugar - sugar_reserved_for_growth) as available leaf supply, ensuring
    // nodes have sugar for their own growth when update_tissue() runs later.
    for (auto& node : nodes_) {
        node->compute_growth_reserve(genome_, world);
    }
}

void Plant::tick_tree(const WorldParams& world, PerfStats* /*perf*/) {
    pre_transport_growth(world);                // nodes claim sugar for growth before vascular
    vascular_transport(*this, genome_, world);  // bulk flow for sugar/water/cytokinin
    pin_transport(*this, genome_);              // PIN-mediated polar auxin transport
    tick_recursive(*nodes_[0], *this, world);
    flush_removals();
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
