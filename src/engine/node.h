#pragma once

#include <cstdint>
#include <vector>
#include <glm/vec3.hpp>

namespace botany {

enum class NodeType { STEM, ROOT, LEAF };

enum class MeristemType { APICAL, AXILLARY, ROOT_APICAL, ROOT_AXILLARY };

class Plant;

class Meristem {
public:
    virtual ~Meristem() = default;

    virtual MeristemType type() const = 0;
    virtual bool is_tip() const = 0;
    virtual void tick(struct Node& node, Plant& plant) = 0;

    bool active = false;
    uint32_t ticks_since_last_node = 0;

protected:
    Meristem(bool active) : active(active) {}
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
    float sugar = 0.0f;
    float leaf_size = 0.0f;

    Meristem* meristem;

    Node(uint32_t id, NodeType type, glm::vec3 position, float radius);

    void add_child(Node* child);
};

} // namespace botany
