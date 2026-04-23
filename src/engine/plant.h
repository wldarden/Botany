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

    // Factory for loading from a snapshot: constructs a Plant with the given
    // genome but no default seed/SA/RA.  Caller must install_node() for every
    // node from the snapshot and then set_next_id() to match the snapshot's
    // next-id counter before the plant is ticked.
    static std::unique_ptr<Plant> from_empty(const Genome& genome);

    // Install a pre-constructed node.  Takes ownership.  Caller is responsible
    // for setting parent / children links, chemicals, and any subclass state
    // before this plant is ticked.  Used by the snapshot loader.
    void install_node(std::unique_ptr<Node> node);

    // Override the id counter so future create_node() calls don't collide with
    // ids already loaded from a snapshot.
    void set_next_id(uint32_t next);

    const Genome& genome() const { return genome_; }

    const Node* seed() const { return nodes_[0].get(); }
    Node* seed_mut() { return nodes_[0].get(); }

    uint32_t node_count() const { return static_cast<uint32_t>(nodes_.size()); }

    uint32_t root_meristem_count() const { return root_meristem_count_; }
    static constexpr uint32_t max_root_meristems = 10000;

    Node* create_node(NodeType type, glm::vec3 position, float radius);

    bool root_meristems_at_cap() const {
        return root_meristem_count_ >= max_root_meristems;
    }

    float total_sugar_produced() const { return total_sugar_produced_; }
    void add_sugar_produced(float amount) { total_sugar_produced_ += amount; }

    void tick(const struct WorldParams& world, struct PerfStats* perf = nullptr);
    void remove_subtree(Node* node);

    void queue_removal(Node* node);
    void flush_removals();

    void for_each_node(std::function<void(const Node&)> fn) const;
    void for_each_node_mut(std::function<void(Node&)> fn);

    uint32_t next_id() { return next_id_++; }

    // Per-tick primary-meristem trackers.  Reset to -1 at the start of each
    // tick; set to the observing meristem's id when a primary SA/RA is seen
    // during tick_tree's DFS walk (cheap, piggybacked on existing traversal).
    // After the walk, if either is still -1 the original primary died this
    // tick (or earlier) — triggers a one-time re-promotion pass that picks
    // the surviving non-primary meristem with the strongest canalization
    // path back to the seed.  Public because the meristem tick sets them.
    int32_t primary_sa_id_this_tick = -1;
    int32_t primary_ra_id_this_tick = -1;

private:
    // Snapshot-load constructor: tag-dispatched so it doesn't collide with the
    // normal (Genome, vec3) seeding constructor.  Leaves nodes_ empty.  The
    // loader populates nodes via install_node() and then calls set_next_id().
    Plant(const Genome& genome, bool /*empty_tag*/);

    void tick_tree(const WorldParams& world, PerfStats* perf);
    void promote_primary_meristems();  // called after tick_tree walk if a primary is missing
    Genome genome_;
    uint32_t next_id_ = 0;
    uint32_t root_meristem_count_ = 0;
    std::vector<std::unique_ptr<Node>> nodes_;
    std::vector<Node*> pending_removals_;
    float total_sugar_produced_ = 0.0f;
    PerfStats* perf_ = nullptr;  // set during tick, null when no perf logging

public:
    PerfStats* perf() const { return perf_; }
};

} // namespace botany
