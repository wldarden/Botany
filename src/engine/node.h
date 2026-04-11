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
    virtual void tick(struct Node& node, Plant& plant, const struct WorldParams& world) = 0;

    bool active = false;
    uint32_t ticks_since_last_node = 0;

protected:
    Meristem(bool active) : active(active) {}
};

struct Node {
    uint32_t id;
    Node* parent;
    std::vector<Node*> children;

    glm::vec3 offset;       // vector from parent to this node (seed: world position)
    glm::vec3 position;     // cached world position (recomputed from offset chain)
    float radius;

    NodeType type;
    uint32_t age;
    float auxin;
    float cytokinin;
    float sugar = 0.0f;
    float leaf_size = 0.0f;
    float light_exposure = 1.0f;   // 0..1, computed per tick from shading
    uint32_t starvation_ticks = 0;
    float gibberellin = 0.0f;            // GA concentration (reset each tick)
    float ethylene = 0.0f;               // ethylene concentration (reset each tick)
    uint32_t senescence_ticks = 0;       // 0 = healthy, >0 = senescing (irreversible)

    Meristem* meristem;

    Node(uint32_t id, NodeType type, glm::vec3 position, float radius);

    void add_child(Node* child);
};

} // namespace botany
