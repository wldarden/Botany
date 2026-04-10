// src/engine/plant.cpp
#include "engine/plant.h"
#include "engine/meristem_types.h"

namespace botany {

Plant::Plant(const Genome& genome, glm::vec3 position)
    : genome_(genome)
{
    // Seed node
    Node* seed = create_node(NodeType::STEM, position, genome.initial_radius);

    // Shoot apical meristem node (just above seed)
    Node* shoot = create_node(NodeType::STEM, position + glm::vec3(0.0f, 0.01f, 0.0f), genome.initial_radius);
    shoot->meristem = create_meristem<ShootApicalMeristem>();
    seed->add_child(shoot);

    // Root apical meristem node (just below seed)
    Node* root = create_node(NodeType::ROOT, position - glm::vec3(0.0f, 0.01f, 0.0f), genome.root_initial_radius);
    root->meristem = create_meristem<RootApicalMeristem>();
    seed->add_child(root);
}

Node* Plant::create_node(NodeType type, glm::vec3 position, float radius) {
    auto node = std::make_unique<Node>(next_id(), type, position, radius);
    Node* ptr = node.get();
    nodes_.push_back(std::move(node));
    return ptr;
}

Leaf* Plant::create_leaf(float size) {
    auto l = std::make_unique<Leaf>(Leaf{size});
    Leaf* ptr = l.get();
    leaves_.push_back(std::move(l));
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
