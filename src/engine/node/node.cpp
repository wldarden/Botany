#include "engine/node/node.h"
#include "engine/plant.h"
#include "engine/sugar.h"
#include "engine/world_params.h"
#include "engine/chemical/chemical_registry.h"
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/node/leaf_node.h"
#include "engine/node/meristems/shoot_apical.h"
#include "engine/node/meristems/shoot_axillary.h"
#include "engine/node/meristems/root_apical.h"
#include "engine/node/meristems/root_axillary.h"

namespace botany {

// --- Node (base) ---

Node::Node(uint32_t id, NodeType type, glm::vec3 position, float radius)
    : id(id)
    , parent(nullptr)
    , offset(position)
    , position(position)
    , radius(radius)
    , type(type)
    , age(0)
{
    // Initialize all chemical map entries to zero
    chemicals[ChemicalID::Auxin] = 0.0f;
    chemicals[ChemicalID::Cytokinin] = 0.0f;
    chemicals[ChemicalID::Gibberellin] = 0.0f;
    chemicals[ChemicalID::Sugar] = 0.0f;
    chemicals[ChemicalID::Ethylene] = 0.0f;
}

void Node::add_child(Node* child) {
    children.push_back(child);
    child->parent = this;
}

void Node::replace_child(Node* old_child, Node* new_child) {
    for (auto& c : children) {
        if (c == old_child) {
            c = new_child;
            new_child->parent = this;
            old_child->parent = nullptr;
            return;
        }
    }
}

void Node::tick(Plant& plant, const WorldParams& world) {
    age++;
    const Genome& g = plant.genome();

    // Maintenance sugar consumption
    float cost = maintenance_cost(g);
    chemical(ChemicalID::Sugar) = std::max(0.0f, chemical(ChemicalID::Sugar) - cost);

    // Cap clamp
    float cap = sugar_cap(*this, g);
    chemical(ChemicalID::Sugar) = std::min(chemical(ChemicalID::Sugar), cap);

    // Starvation tracking + death
    if (chemical(ChemicalID::Sugar) <= 0.0f) starvation_ticks++;
    else starvation_ticks = 0;

    if (starvation_ticks >= world.starvation_ticks_max && parent) {
        die(plant);
        return;
    }

    transport_chemicals(g);
    grow(plant, world);
}

void Node::grow(Plant& /*plant*/, const WorldParams& /*world*/) {}

void Node::die(Plant& plant) {
    // Detach from parent
    if (parent) {
        auto& siblings = parent->children;
        siblings.erase(std::remove(siblings.begin(), siblings.end(), this), siblings.end());
        parent = nullptr;
    }

    // Collect self + all descendants, queue for deferred removal
    std::vector<Node*> to_die = {this};
    size_t i = 0;
    while (i < to_die.size()) {
        for (Node* child : to_die[i]->children) {
            to_die.push_back(child);
        }
        i++;
    }

    // Clear children so recursive tick snapshot is empty
    children.clear();

    for (Node* n : to_die) {
        plant.queue_removal(n);
    }
}

float Node::maintenance_cost(const Genome& /*g*/) const {
    return 0.0f;
}

// Generic biased transport: blend of gradient diffusion and directional push.
// bias < 0: basipetal (push to parent). bias > 0: acropetal (pull from parent). 0: pure gradient.
static void transport_chemical(float& my_val, float& parent_val,
                               float rate, float bias, float decay) {
    float abs_b = std::abs(bias);
    float gw = 1.0f - abs_b;
    float dw = abs_b;

    float gradient_flow = (my_val - parent_val) * gw * rate;
    float directional_flow = (bias < 0)
        ?  my_val * dw * rate
        : -parent_val * dw * rate;

    float flow = gradient_flow + directional_flow;
    if (flow > 0.0f) flow = std::min(flow, my_val);
    else             flow = std::max(flow, -parent_val);

    my_val -= flow;
    parent_val += flow;
    my_val *= (1.0f - decay);
}

void Node::transport_chemicals(const Genome& g) {
    if (parent) {
        // Hormones: biased local transport via registry
        for (const auto& hp : hormone_params(g)) {
            transport_chemical(chemical(hp.id), parent->chemical(hp.id),
                hp.transport_rate, hp.directional_bias, hp.decay_rate);
        }

        // Sugar: cap-aware gradient diffusion (unchanged)
        float my_cap = sugar_cap(*this, g);
        float parent_cap = sugar_cap(*parent, g);
        float diff = chemical(ChemicalID::Sugar) - parent->chemical(ChemicalID::Sugar);
        float min_r = std::min(radius, parent->radius);
        if (type == NodeType::LEAF || is_meristem()
            || parent->type == NodeType::LEAF || parent->is_meristem())
            min_r = std::max(min_r, 0.01f);
        float conductance = std::min(min_r * min_r * 3.14159f * g.sugar_transport_conductance, 0.25f);
        float flow = diff * conductance;
        if (flow > 0.0f) {
            float headroom = std::max(0.0f, parent_cap - parent->chemical(ChemicalID::Sugar));
            flow = std::min({flow, chemical(ChemicalID::Sugar), headroom});
        } else {
            float headroom = std::max(0.0f, my_cap - chemical(ChemicalID::Sugar));
            flow = std::max({flow, -parent->chemical(ChemicalID::Sugar), -headroom});
        }
        chemical(ChemicalID::Sugar) -= flow;
        parent->chemical(ChemicalID::Sugar) += flow;
    } else {
        // Seed: just decay hormones
        for (const auto& hp : hormone_params(g)) {
            chemical(hp.id) *= (1.0f - hp.decay_rate);
        }
    }
}

bool Node::is_meristem() const {
    return type == NodeType::SHOOT_APICAL || type == NodeType::SHOOT_AXILLARY
        || type == NodeType::ROOT_APICAL  || type == NodeType::ROOT_AXILLARY;
}

// Downcasting helpers
StemNode*       Node::as_stem()       { return type == NodeType::STEM ? static_cast<StemNode*>(this) : nullptr; }
const StemNode* Node::as_stem() const { return type == NodeType::STEM ? static_cast<const StemNode*>(this) : nullptr; }
RootNode*       Node::as_root()       { return type == NodeType::ROOT ? static_cast<RootNode*>(this) : nullptr; }
const RootNode* Node::as_root() const { return type == NodeType::ROOT ? static_cast<const RootNode*>(this) : nullptr; }
LeafNode*       Node::as_leaf()       { return type == NodeType::LEAF ? static_cast<LeafNode*>(this) : nullptr; }
const LeafNode* Node::as_leaf() const { return type == NodeType::LEAF ? static_cast<const LeafNode*>(this) : nullptr; }

ShootApicalNode*       Node::as_shoot_apical()       { return type == NodeType::SHOOT_APICAL ? static_cast<ShootApicalNode*>(this) : nullptr; }
const ShootApicalNode* Node::as_shoot_apical() const { return type == NodeType::SHOOT_APICAL ? static_cast<const ShootApicalNode*>(this) : nullptr; }
ShootAxillaryNode*       Node::as_shoot_axillary()       { return type == NodeType::SHOOT_AXILLARY ? static_cast<ShootAxillaryNode*>(this) : nullptr; }
const ShootAxillaryNode* Node::as_shoot_axillary() const { return type == NodeType::SHOOT_AXILLARY ? static_cast<const ShootAxillaryNode*>(this) : nullptr; }
RootApicalNode*       Node::as_root_apical()       { return type == NodeType::ROOT_APICAL ? static_cast<RootApicalNode*>(this) : nullptr; }
const RootApicalNode* Node::as_root_apical() const { return type == NodeType::ROOT_APICAL ? static_cast<const RootApicalNode*>(this) : nullptr; }
RootAxillaryNode*       Node::as_root_axillary()       { return type == NodeType::ROOT_AXILLARY ? static_cast<RootAxillaryNode*>(this) : nullptr; }
const RootAxillaryNode* Node::as_root_axillary() const { return type == NodeType::ROOT_AXILLARY ? static_cast<const RootAxillaryNode*>(this) : nullptr; }

} // namespace botany
