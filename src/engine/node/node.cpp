#include "engine/node/node.h"
#include "engine/plant.h"
#include "engine/sugar.h"
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/node/leaf_node.h"
#include "engine/node/meristem_node.h"
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
    , auxin(0.0f)
    , cytokinin(0.0f)
{}

void Node::add_child(Node* child) {
    children.push_back(child);
    child->parent = this;
}

void Node::tick(Plant& plant, const WorldParams& /*world*/) {
    age++;
    transport_chemicals(plant.genome());
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
        // Auxin: basipetal
        transport_chemical(auxin, parent->auxin,
            g.auxin_transport_rate, g.auxin_directional_bias, g.auxin_decay_rate);
        // Cytokinin: acropetal
        transport_chemical(cytokinin, parent->cytokinin,
            g.cytokinin_transport_rate, g.cytokinin_directional_bias, g.cytokinin_decay_rate);
        // Sugar: gradient-based, cap-aware, conductance scales with radius
        float my_cap = sugar_cap(*this, g);
        float parent_cap = sugar_cap(*parent, g);
        float diff = sugar - parent->sugar;
        float min_r = std::min(radius, parent->radius);
        if (type == NodeType::LEAF || is_meristem()
            || parent->type == NodeType::LEAF || parent->is_meristem())
            min_r = std::max(min_r, 0.01f);
        float conductance = std::min(min_r * min_r * 3.14159f * g.sugar_transport_conductance, 0.25f);
        float flow = diff * conductance;
        if (flow > 0.0f) {
            float headroom = std::max(0.0f, parent_cap - parent->sugar);
            flow = std::min({flow, sugar, headroom});
        } else {
            float headroom = std::max(0.0f, my_cap - sugar);
            flow = std::max({flow, -parent->sugar, -headroom});
        }
        sugar -= flow;
        parent->sugar += flow;
    } else {
        // Seed: no parent, just decay hormones
        auxin *= (1.0f - g.auxin_decay_rate);
        cytokinin *= (1.0f - g.cytokinin_decay_rate);
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

MeristemNode*       Node::as_meristem()       { return is_meristem() ? static_cast<MeristemNode*>(this) : nullptr; }
const MeristemNode* Node::as_meristem() const { return is_meristem() ? static_cast<const MeristemNode*>(this) : nullptr; }

ShootApicalNode*       Node::as_shoot_apical()       { return type == NodeType::SHOOT_APICAL ? static_cast<ShootApicalNode*>(this) : nullptr; }
const ShootApicalNode* Node::as_shoot_apical() const { return type == NodeType::SHOOT_APICAL ? static_cast<const ShootApicalNode*>(this) : nullptr; }
ShootAxillaryNode*       Node::as_shoot_axillary()       { return type == NodeType::SHOOT_AXILLARY ? static_cast<ShootAxillaryNode*>(this) : nullptr; }
const ShootAxillaryNode* Node::as_shoot_axillary() const { return type == NodeType::SHOOT_AXILLARY ? static_cast<const ShootAxillaryNode*>(this) : nullptr; }
RootApicalNode*       Node::as_root_apical()       { return type == NodeType::ROOT_APICAL ? static_cast<RootApicalNode*>(this) : nullptr; }
const RootApicalNode* Node::as_root_apical() const { return type == NodeType::ROOT_APICAL ? static_cast<const RootApicalNode*>(this) : nullptr; }
RootAxillaryNode*       Node::as_root_axillary()       { return type == NodeType::ROOT_AXILLARY ? static_cast<RootAxillaryNode*>(this) : nullptr; }
const RootAxillaryNode* Node::as_root_axillary() const { return type == NodeType::ROOT_AXILLARY ? static_cast<const RootAxillaryNode*>(this) : nullptr; }

} // namespace botany
