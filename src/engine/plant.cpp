// src/engine/plant.cpp
#include "engine/plant.h"
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/node/leaf_node.h"
#include "engine/node/meristems/shoot_apical.h"
#include "engine/node/meristems/shoot_axillary.h"
#include "engine/node/meristems/root_apical.h"
#include "engine/node/meristems/root_axillary.h"
#include "engine/sugar.h"
#include "engine/gibberellin.h"
#include "engine/ethylene.h"
#include "engine/world_params.h"
#include <algorithm>
#include <unordered_set>

namespace botany {

Plant::Plant(const Genome& genome, glm::vec3 position)
    : genome_(genome)
{
    // Seed node
    Node* seed = create_node(NodeType::STEM, position, genome.initial_radius);
    seed->chemical(ChemicalID::Sugar) = genome.seed_sugar;

    // Shoot apical meristem node (child of seed)
    Node* shoot = create_node(NodeType::SHOOT_APICAL, glm::vec3(0.0f, 0.01f, 0.0f), genome.initial_radius);
    seed->add_child(shoot);

    // Root apical meristem node (child of seed)
    Node* root = create_node(NodeType::ROOT_APICAL, glm::vec3(0.0f, -0.01f, 0.0f), genome.root_initial_radius);
    seed->add_child(root);

    recompute_world_positions();
}

void Plant::tick(const WorldParams& world) {
    // Auxin, cytokinin, sugar diffusion: now transported locally
    // by each node during Node::transport_chemicals()
    compute_gibberellin(*this);
    transport_sugar(*this, world);
    compute_ethylene(*this, world);
    process_abscission(*this);
    tick_tree(world);
    recompute_world_positions();
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
        if (n->type == NodeType::ROOT_APICAL || n->type == NodeType::ROOT_AXILLARY) {
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
        if (n->type == NodeType::ROOT_APICAL || n->type == NodeType::ROOT_AXILLARY) {
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

void Plant::tick_tree(const WorldParams& world) {
    tick_recursive(*nodes_[0], *this, world);
    flush_removals();
}

static void recompute_recursive(Node& node) {
    for (Node* child : node.children) {
        child->position = node.position + child->offset;
        recompute_recursive(*child);
    }
}

void Plant::recompute_world_positions() {
    Node& seed = *nodes_[0];
    seed.position = seed.offset;
    recompute_recursive(seed);
}

Node* Plant::create_node(NodeType type, glm::vec3 position, float radius) {
    std::unique_ptr<Node> node;
    uint32_t id = next_id();
    switch (type) {
        case NodeType::STEM: node = std::make_unique<StemNode>(id, position, radius); break;
        case NodeType::ROOT: node = std::make_unique<RootNode>(id, position, radius); break;
        case NodeType::LEAF: node = std::make_unique<LeafNode>(id, position, radius); break;
        case NodeType::SHOOT_APICAL: node = std::make_unique<ShootApicalNode>(id, position, radius); break;
        case NodeType::SHOOT_AXILLARY: node = std::make_unique<ShootAxillaryNode>(id, position, radius); break;
        case NodeType::ROOT_APICAL:
            node = std::make_unique<RootApicalNode>(id, position, radius);
            root_meristem_count_++;
            break;
        case NodeType::ROOT_AXILLARY:
            node = std::make_unique<RootAxillaryNode>(id, position, radius);
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
