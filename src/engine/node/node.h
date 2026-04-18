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
    APICAL, ROOT_APICAL
};

class Plant;
struct WorldParams;
struct Genome;
class StemNode;
class RootNode;
class LeafNode;
class ApicalNode;
class RootApicalNode;

class Node {
public:
    // --- Graph structure ---
    uint32_t id;
    Node* parent;
    std::vector<Node*> children;

    // --- Spatial ---
    glm::vec3 offset;       // vector from parent to this node (seed: world position)
    glm::vec3 rest_offset{0.0f};  // stress-free direction (for elastic recovery after droop)
    glm::vec3 position;     // cached world position (recomputed from offset chain)
    float radius;

    // --- Identity ---
    NodeType type;
    uint32_t age;

    // --- Chemicals ---
    uint32_t starvation_ticks = 0;

    // --- Mass / stress (computed each tick, children's values are one tick stale) ---
    float total_mass = 0.0f;       // self mass + all children's total_mass
    glm::vec3 mass_moment{0.0f};   // self_mass * position + Σ child.mass_moment
    float stress = 0.0f;           // torque / cross-section (structural load)

    // --- Per-tick sugar accounting (updated every tick, one tick stale in inspector) ---
    float tick_sugar_maintenance = 0.0f; // consumed by maintenance (positive = consumed)
    float tick_sugar_activity    = 0.0f; // net change from tissue work: positive = produced, negative = spent on growth
    float tick_sugar_transport   = 0.0f; // net change from transport: negative = exported, positive = imported

    // Chemical storage — map-based, sole storage for all chemical values.
    std::unordered_map<ChemicalID, float> chemicals;

    float& chemical(ChemicalID id) { return chemicals[id]; }
    float chemical(ChemicalID id) const {
        auto it = chemicals.find(id);
        return it != chemicals.end() ? it->second : 0.0f;
    }

    // Transport received buffer — chemicals received from parent's Phase 2
    // this tick. NOT visible to this node's own transport (anti-teleportation).
    // Flushed into chemical() after this node's transport completes.
    std::unordered_map<ChemicalID, float> transport_received;

    // Canalization — per-child transport bias (stored on parent, keyed by child pointer)
    std::unordered_map<Node*, float> auxin_flow_bias;       // transient — fast, decays toward PIN saturation
    std::unordered_map<Node*, float> last_auxin_flux;       // transient per-tick: auxin moved per child

    float get_bias_multiplier(Node* child, const Genome& g) const;

    // Returns the auxin_flow_bias this node's parent has recorded for it.
    // Stored on the parent keyed by child pointer — zero if no parent or no entry yet.
    // For the seed (no parent) returns the max of all its children's biases.
    float get_parent_auxin_flow_bias() const;

    // --- Lifecycle ---
    Node(uint32_t id, NodeType type, glm::vec3 position, float radius);
    virtual ~Node() = default;

    // --- Tree operations ---
    void add_child(Node* child);
    void replace_child(Node* old_child, Node* new_child);

    // --- Tick pipeline ---
    void tick(Plant& plant, const WorldParams& world);       // non-virtual: all universal processes
    virtual float maintenance_cost(const WorldParams& world) const;

    void transport_with_children(const Genome& g);
    void update_canalization(const Genome& g);
    void decay_chemicals(const Genome& g);
    void die(Plant& plant);

private:
    // --- Tick helpers (called by tick in order) ---
    void sync_world_position();
    bool pay_maintenance(Plant& plant, const WorldParams& world);
    bool update_physics(Plant& plant, const Genome& g, const WorldParams& world);
    void transport_chemicals(const Genome& g);
    void deduct_maintenance_sugar(const WorldParams& world);
    bool check_starvation(Plant& plant, const WorldParams& world);
    void compute_mass(const Genome& g, const WorldParams& world);
    void compute_stress(const Genome& g, const WorldParams& world);
    bool apply_droop_and_break(Plant& plant, const Genome& g, const WorldParams& world);

public:

    // --- Type queries ---
    bool is_meristem() const;

    // --- Downcasting (gated on NodeType enum, no RTTI) ---
    StemNode*       as_stem();
    const StemNode* as_stem() const;
    RootNode*       as_root();
    const RootNode* as_root() const;
    LeafNode*       as_leaf();
    const LeafNode* as_leaf() const;
    ApicalNode*       as_apical();
    const ApicalNode* as_apical() const;
    RootApicalNode*       as_root_apical();
    const RootApicalNode* as_root_apical() const;

protected:
    virtual void update_tissue(Plant& plant, const WorldParams& world);  // subclass entry point
};

} // namespace botany
