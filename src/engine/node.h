#pragma once

#include <cstdint>
#include <vector>
#include <glm/vec3.hpp>

namespace botany {

enum class NodeType { STEM, ROOT, LEAF };

enum class MeristemType { APICAL, AXILLARY, ROOT_APICAL, ROOT_AXILLARY };

class Plant;
struct WorldParams;
class StemNode;
class RootNode;
class LeafNode;

class Meristem {
public:
    virtual ~Meristem() = default;

    virtual MeristemType type() const = 0;
    virtual bool is_tip() const = 0;
    virtual void tick(class Node& node, Plant& plant, const WorldParams& world) = 0;

    bool active = false;
    uint32_t ticks_since_last_node = 0;

protected:
    Meristem(bool active) : active(active) {}
};

class Node {
public:
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
    uint32_t starvation_ticks = 0;
    float gibberellin = 0.0f;
    float ethylene = 0.0f;

    Meristem* meristem;

    Node(uint32_t id, NodeType type, glm::vec3 position, float radius);
    virtual ~Node() = default;

    void add_child(Node* child);
    virtual void tick(Plant& plant, const WorldParams& world);

    // Fast downcasting (gated on node_type enum, no RTTI)
    StemNode*       as_stem();
    const StemNode* as_stem() const;
    RootNode*       as_root();
    const RootNode* as_root() const;
    LeafNode*       as_leaf();
    const LeafNode* as_leaf() const;
};

class StemNode : public Node {
public:
    StemNode(uint32_t id, glm::vec3 position, float radius);
    void tick(Plant& plant, const WorldParams& world) override;
};

class RootNode : public Node {
public:
    RootNode(uint32_t id, glm::vec3 position, float radius);
    void tick(Plant& plant, const WorldParams& world) override;
};

class LeafNode : public Node {
public:
    float leaf_size = 0.0f;
    float light_exposure = 1.0f;
    uint32_t senescence_ticks = 0;

    LeafNode(uint32_t id, glm::vec3 position, float radius);
    void tick(Plant& plant, const WorldParams& world) override;
};

} // namespace botany
