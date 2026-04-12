// src/engine/plant.h
#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <vector>
#include <glm/vec3.hpp>
#include "engine/genome.h"
#include "engine/node/node.h"

namespace botany {

class Plant {
public:
    Plant(const Genome& genome, glm::vec3 position);

    const Genome& genome() const { return genome_; }

    const Node* seed() const { return nodes_[0].get(); }
    Node* seed_mut() { return nodes_[0].get(); }

    uint32_t node_count() const { return static_cast<uint32_t>(nodes_.size()); }

    uint32_t root_meristem_count() const { return root_meristem_count_; }
    static constexpr uint32_t max_root_meristems = 100;

    Node* create_node(NodeType type, glm::vec3 position, float radius);

    bool root_meristems_at_cap() const {
        return root_meristem_count_ >= max_root_meristems;
    }

    float total_sugar_produced() const { return total_sugar_produced_; }
    void add_sugar_produced(float amount) { total_sugar_produced_ += amount; }

    void tick(const struct WorldParams& world);
    void remove_subtree(Node* node);
    void recompute_world_positions();

    void queue_removal(Node* node);
    void flush_removals();

    void for_each_node(std::function<void(const Node&)> fn) const;
    void for_each_node_mut(std::function<void(Node&)> fn);

    uint32_t next_id() { return next_id_++; }

private:
    void tick_tree(const WorldParams& world);
    Genome genome_;
    uint32_t next_id_ = 0;
    uint32_t root_meristem_count_ = 0;
    std::vector<std::unique_ptr<Node>> nodes_;
    std::vector<Node*> pending_removals_;
    float total_sugar_produced_ = 0.0f;
};

} // namespace botany
