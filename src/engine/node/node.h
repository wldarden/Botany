#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>
#include <glm/vec3.hpp>
#include "engine/chemical/chemical.h"

// Hash for ChemicalID so it works as unordered_map key
template<>
struct std::hash<botany::ChemicalID> {
    std::size_t operator()(botany::ChemicalID id) const noexcept {
        return static_cast<std::size_t>(id);
    }
};

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
class ShootApicalNode;
class ShootAxillaryNode;
class RootApicalNode;
class RootAxillaryNode;

class Node {
public:
    // --- Graph structure ---
    uint32_t id;
    Node* parent;
    std::vector<Node*> children;

    // --- Spatial ---
    glm::vec3 offset;       // vector from parent to this node (seed: world position)
    glm::vec3 position;     // cached world position (recomputed from offset chain)
    float radius;

    // --- Identity ---
    NodeType type;
    uint32_t age;

    // --- Chemicals ---
    float auxin;
    float cytokinin;
    float sugar = 0.0f;
    uint32_t starvation_ticks = 0;
    float gibberellin = 0.0f;
    float ethylene = 0.0f;

    // Chemical storage — map-based, replaces individual fields.
    // During migration, both map and old fields exist.
    std::unordered_map<ChemicalID, float> chemicals;

    float& chemical(ChemicalID id) { return chemicals[id]; }
    float chemical(ChemicalID id) const {
        auto it = chemicals.find(id);
        return it != chemicals.end() ? it->second : 0.0f;
    }

    // --- Lifecycle ---
    Node(uint32_t id, NodeType type, glm::vec3 position, float radius);
    virtual ~Node() = default;

    // --- Tree operations ---
    void add_child(Node* child);
    void replace_child(Node* old_child, Node* new_child);

    // --- Tick pipeline ---
    virtual void tick(Plant& plant, const WorldParams& world);
    virtual float maintenance_cost(const struct Genome& g) const;
    virtual void grow(Plant& plant, const WorldParams& world);
    void transport_chemicals(const struct Genome& g);
    void die(Plant& plant);

    // --- Type queries ---
    bool is_meristem() const;

    // --- Downcasting (gated on NodeType enum, no RTTI) ---
    StemNode*       as_stem();
    const StemNode* as_stem() const;
    RootNode*       as_root();
    const RootNode* as_root() const;
    LeafNode*       as_leaf();
    const LeafNode* as_leaf() const;
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
