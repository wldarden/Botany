#pragma once

#include <cstdint>
#include <vector>
#include <glm/vec3.hpp>

namespace botany {

enum class NodeType { STEM, ROOT };

enum class MeristemType { APICAL, AXILLARY, ROOT_APICAL, ROOT_AXILLARY };

struct Meristem {
    MeristemType type;
    bool active;
    uint32_t ticks_since_last_node;
};

struct Leaf {
    float size;
};

struct Node {
    uint32_t id;
    Node* parent;
    std::vector<Node*> children;

    glm::vec3 position;
    float radius;

    NodeType type;
    uint32_t age;
    float auxin;
    float cytokinin;

    Meristem* meristem;
    Leaf* leaf;

    Node(uint32_t id, NodeType type, glm::vec3 position, float radius);

    void add_child(Node* child);
};

} // namespace botany
