// src/engine/plant.cpp
#include "engine/plant.h"

namespace botany {

Plant::Plant(const Genome& genome, glm::vec3 position)
    : genome_(genome)
{
    // Seed node
    Node* seed = create_node(NodeType::STEM, position, genome.initial_radius);

    // Shoot apical meristem node (just above seed)
    Node* shoot = create_node(NodeType::STEM, position + glm::vec3(0.0f, 0.01f, 0.0f), genome.initial_radius);
    Meristem* shoot_m = create_meristem(MeristemType::APICAL, true);
    shoot->meristem = shoot_m;
    seed->add_child(shoot);

    // Root apical meristem node (just below seed)
    Node* root = create_node(NodeType::ROOT, position - glm::vec3(0.0f, 0.01f, 0.0f), genome.initial_radius);
    Meristem* root_m = create_meristem(MeristemType::ROOT_APICAL, true);
    root->meristem = root_m;
    seed->add_child(root);
}

Node* Plant::create_node(NodeType type, glm::vec3 position, float radius) {
    auto node = std::make_unique<Node>(next_id(), type, position, radius);
    Node* ptr = node.get();
    nodes_.push_back(std::move(node));
    return ptr;
}

Meristem* Plant::create_meristem(MeristemType type, bool active) {
    auto m = std::make_unique<Meristem>(Meristem{type, active, 0});
    Meristem* ptr = m.get();
    meristems_.push_back(std::move(m));
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
