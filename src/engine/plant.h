// src/engine/plant.h
#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <vector>
#include <glm/vec3.hpp>
#include "engine/genome.h"
#include "engine/node.h"

namespace botany {

class Plant {
public:
    Plant(const Genome& genome, glm::vec3 position);

    const Genome& genome() const { return genome_; }

    const Node* seed() const { return nodes_[0].get(); }
    Node* seed_mut() { return nodes_[0].get(); }

    uint32_t node_count() const { return static_cast<uint32_t>(nodes_.size()); }

    Node* create_node(NodeType type, glm::vec3 position, float radius);
    Meristem* create_meristem(MeristemType type, bool active);
    Leaf* create_leaf(float size);

    void for_each_node(std::function<void(const Node&)> fn) const;
    void for_each_node_mut(std::function<void(Node&)> fn);

    uint32_t next_id() { return next_id_++; }

private:
    Genome genome_;
    uint32_t next_id_ = 0;
    std::vector<std::unique_ptr<Node>> nodes_;
    std::vector<std::unique_ptr<Meristem>> meristems_;
    std::vector<std::unique_ptr<Leaf>> leaves_;
};

} // namespace botany
