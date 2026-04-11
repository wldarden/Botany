#pragma once

#include <cstdint>
#include <vector>
#include <glm/vec3.hpp>

namespace botany {

enum class NodeType {
    STEM, ROOT, LEAF,
    SHOOT_APICAL, SHOOT_AXILLARY, ROOT_APICAL, ROOT_AXILLARY
};

class Plant;
struct WorldParams;
class StemNode;
class RootNode;
class LeafNode;
class MeristemNode;
class ShootApicalNode;
class ShootAxillaryNode;
class RootApicalNode;
class RootAxillaryNode;

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

    Node(uint32_t id, NodeType type, glm::vec3 position, float radius);
    virtual ~Node() = default;

    void add_child(Node* child);
    virtual void tick(Plant& plant, const WorldParams& world);
    void transport_chemicals(const struct Genome& g);

    bool is_meristem() const;

    // Fast downcasting (gated on node_type enum, no RTTI)
    StemNode*       as_stem();
    const StemNode* as_stem() const;
    RootNode*       as_root();
    const RootNode* as_root() const;
    LeafNode*       as_leaf();
    const LeafNode* as_leaf() const;
    MeristemNode*       as_meristem();
    const MeristemNode* as_meristem() const;
    ShootApicalNode*       as_shoot_apical();
    const ShootApicalNode* as_shoot_apical() const;
    ShootAxillaryNode*       as_shoot_axillary();
    const ShootAxillaryNode* as_shoot_axillary() const;
    RootApicalNode*       as_root_apical();
    const RootApicalNode* as_root_apical() const;
    RootAxillaryNode*       as_root_axillary();
    const RootAxillaryNode* as_root_axillary() const;
};

} // namespace botany
