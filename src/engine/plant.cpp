// src/engine/plant.cpp
#include "engine/plant.h"
#include "engine/hormone.h"
#include "engine/sugar.h"
#include "engine/gibberellin.h"
#include "engine/ethylene.h"
#include "engine/meristems/meristem.h"
#include "engine/meristems/meristem_types.h"
#include "engine/world_params.h"
#include <algorithm>
#include <unordered_set>

namespace botany {

Plant::Plant(const Genome& genome, glm::vec3 position)
    : genome_(genome)
{
    // Seed node
    Node* seed = create_node(NodeType::STEM, position, genome.initial_radius);
    seed->sugar = genome.seed_sugar;

    // Shoot apical meristem node (offset slightly above seed)
    Node* shoot = create_node(NodeType::STEM, glm::vec3(0.0f, 0.01f, 0.0f), genome.initial_radius);
    shoot->meristem = create_meristem<ShootApicalMeristem>();
    seed->add_child(shoot);

    // Root apical meristem node (offset slightly below seed)
    Node* root = create_node(NodeType::ROOT, glm::vec3(0.0f, -0.01f, 0.0f), genome.root_initial_radius);
    root->meristem = create_meristem<RootApicalMeristem>();
    seed->add_child(root);

    recompute_world_positions();
}

void Plant::tick(const WorldParams& world) {
    transport_auxin(*this);
    transport_cytokinin(*this);
    compute_gibberellin(*this);
    transport_sugar(*this, world);
    compute_ethylene(*this, world);
    process_abscission(*this);
    tick_meristems(*this, world);
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

    // Build a set of raw pointers to remove for fast lookup
    std::unordered_set<Node*> remove_set(to_remove.begin(), to_remove.end());

    // Collect meristems that belong to removed nodes and update root meristem count
    std::unordered_set<Meristem*> meristems_to_remove;
    for (Node* n : to_remove) {
        if (n->meristem) {
            if (n->meristem->type() == MeristemType::ROOT_APICAL ||
                n->meristem->type() == MeristemType::ROOT_AXILLARY) {
                if (root_meristem_count_ > 0) root_meristem_count_--;
            }
            meristems_to_remove.insert(n->meristem);
            n->meristem = nullptr;
        }
    }

    // Remove meristems
    meristems_.erase(
        std::remove_if(meristems_.begin(), meristems_.end(),
            [&](const std::unique_ptr<Meristem>& m) {
                return meristems_to_remove.count(m.get()) > 0;
            }),
        meristems_.end()
    );

    // Remove nodes
    nodes_.erase(
        std::remove_if(nodes_.begin(), nodes_.end(),
            [&](const std::unique_ptr<Node>& n) {
                return remove_set.count(n.get()) > 0;
            }),
        nodes_.end()
    );
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
