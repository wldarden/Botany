#include "engine/compression.h"
#include "engine/plant.h"
#include "engine/node/node.h"
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/genome.h"
#include <cmath>
#include <glm/geometric.hpp>

namespace botany {
namespace {

float effective_max_combined_length(const Genome& g, const CompressionParams& p) {
    return p.max_combined_length > 0.0f ? p.max_combined_length
                                        : 2.0f * g.max_internode_length;
}

// Count children of N that participate in the structural chain.  Leaves,
// apicals, and root_apicals are not structural — they decorate an internode
// but do not branch the trunk/root chain.
uint32_t count_structural_children(const Node& n) {
    uint32_t k = 0;
    for (const Node* c : n.children) {
        if (c->type == NodeType::STEM || c->type == NodeType::ROOT) ++k;
    }
    return k;
}

bool is_conduit_type(NodeType t) {
    return t == NodeType::STEM || t == NodeType::ROOT;
}

} // namespace

bool can_merge(const Node& parent, const Node& child,
               const Genome& g, const CompressionParams& params) {
    // Rule 1: same type, and must be a conduit type.
    if (parent.type != child.type) return false;
    if (!is_conduit_type(parent.type)) return false;

    // Rule 2: parent must have a grandparent (parent != seed).
    if (!parent.parent) return false;

    // Rule 3: parent has exactly ONE structural child and it's this child.
    if (count_structural_children(parent) != 1) return false;
    if (child.parent != &parent) return false;

    // Rule 4: angle between parent.offset and child.offset is small.
    float lenP = glm::length(parent.offset);
    float lenC = glm::length(child.offset);
    if (lenP < 1e-6f || lenC < 1e-6f) return false;
    float cos_gate = std::cos(params.max_angle_rad);
    float cos_pair = glm::dot(parent.offset, child.offset) / (lenP * lenC);
    if (cos_pair < cos_gate) return false;

    // Rule 5: radius ratio within tolerance.
    float rmax = std::max(parent.radius, child.radius);
    if (rmax < 1e-6f) return false;
    if (std::fabs(parent.radius - child.radius) / rmax > params.max_radius_ratio) return false;

    // Rule 6: both past elongation maturation.
    if (parent.age < g.internode_maturation_ticks) return false;
    if (child.age  < g.internode_maturation_ticks) return false;

    // Rule 7: combined length below threshold.
    if (lenP + lenC >= effective_max_combined_length(g, params)) return false;

    return true;
}

CompressionResult compress_plant(Plant& /*plant*/, const CompressionParams& /*params*/) {
    CompressionResult r;
    return r; // Filled in by Tasks 3-4.
}

} // namespace botany
