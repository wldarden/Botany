#pragma once

#include <cstdint>
#include <iostream>
#include <vector>
#include <glm/vec3.hpp>
#include "engine/genome.h"
#include "engine/node/node.h"

namespace botany {

class Engine;

struct NodeSnapshot {
    uint32_t id;
    uint32_t parent_id;  // UINT32_MAX if no parent
    NodeType type;
    glm::vec3 position;
    float radius;
    float auxin;
    float cytokinin;
    float sugar;
    float leaf_size;
    glm::vec3 facing;
};

struct TickSnapshot {
    uint32_t tick_number;
    std::vector<NodeSnapshot> nodes;
};

struct RecordingHeader {
    uint32_t num_ticks;
    Genome genome;
};

struct Recording {
    Genome genome;
};

void save_recording_header(std::ostream& out, const Genome& genome, uint32_t num_ticks);
RecordingHeader load_recording_header(std::istream& in);

void save_tick(std::ostream& out, const Engine& engine, uint32_t plant_id);
TickSnapshot load_tick(std::istream& in);

} // namespace botany
