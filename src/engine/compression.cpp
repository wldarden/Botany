#include "engine/compression.h"
#include "engine/plant.h"
#include "engine/node/node.h"
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/genome.h"
#include "engine/sugar.h"              // sugar_cap, water_cap
#include "engine/vascular_sub_stepped.h" // phloem_capacity, xylem_capacity
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <utility>
#include <vector>
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

void merge_pair(Plant& plant, Node& parent, Node& child, CompressionResult& result) {
    const Genome& g = plant.genome();

    const float lenP = glm::length(parent.offset);
    const float lenC = glm::length(child.offset);
    const float total_len = lenP + lenC;

    // Volume-preserving radius: sqrt((rP²·LP + rC²·LC) / (LP + LC))
    const float r_merged = (total_len > 1e-8f)
        ? std::sqrt((parent.radius * parent.radius * lenP +
                     child.radius  * child.radius  * lenC) / total_len)
        : parent.radius;

    // Extend P's geometry so it now occupies C's old world position.
    const glm::vec3 new_offset      = parent.offset      + child.offset;
    const glm::vec3 new_rest_offset = parent.rest_offset + child.rest_offset;

    // Sum local chemicals.
    for (const auto& kv : child.local().chemicals) {
        parent.local().chemical(kv.first) += kv.second;
    }

    // Sum conduit pool chemicals.  Both are STEM or both ROOT (can_merge gate).
    if (auto* pp = parent.phloem()) {
        if (auto* cp = child.phloem()) {
            for (const auto& kv : cp->chemicals) {
                pp->chemical(kv.first) += kv.second;
            }
        }
    }
    if (auto* px = parent.xylem()) {
        if (auto* cx = child.xylem()) {
            for (const auto& kv : cx->chemicals) {
                px->chemical(kv.first) += kv.second;
            }
        }
    }

    // Remove C from P's children list (we are absorbing C).
    parent.children.erase(
        std::remove(parent.children.begin(), parent.children.end(), &child),
        parent.children.end());

    // Reparent C's children onto P.  Their world-space positions stay put
    // (P now sits where C did), so their offsets (= child_position - parent_position)
    // are unchanged.
    for (Node* gc : child.children) {
        gc->parent = &parent;
        parent.children.push_back(gc);
    }
    child.children.clear(); // prevent die() cascade

    // Canalization bias transfer:
    // 1. Every C→grandchild entry migrates into parent.auxin_flow_bias.
    for (const auto& kv : child.auxin_flow_bias) {
        parent.auxin_flow_bias[kv.first] = kv.second;
    }
    // 2. The parent→child edge no longer exists — drop its entry.
    parent.auxin_flow_bias.erase(&child);

    // Apply new geometry BEFORE clamping (caps depend on radius + length).
    parent.offset      = new_offset;
    parent.rest_offset = new_rest_offset;
    parent.radius      = r_merged;

    // Clamp local + pool chemicals to new caps; log clamp losses.
    auto clamp_local = [&](ChemicalID id, float cap, float& running_delta) {
        float& v = parent.local().chemical(id);
        if (v > cap) {
            running_delta -= (v - cap);
            v = cap;
        }
    };
    const float sc = sugar_cap(parent, g);
    const float wc = water_cap(parent, g);
    clamp_local(ChemicalID::Sugar, sc, result.delta_sugar);
    clamp_local(ChemicalID::Water, wc, result.delta_water);
    // Auxin + cytokinin: no explicit per-node cap in the current model;
    // record nothing but keep the field wired for future clamps.
    (void)result.delta_auxin;
    (void)result.delta_cytokinin;

    // Pool caps: phloem_capacity / xylem_capacity are concentration
    // denominators (volume of pool in dm³), not hard ceilings on stored mass.
    // Vascular operation self-limits pool contents via Jacobi + radial, so we
    // do not clamp here — doing so would silently destroy mass on merges
    // where the summed pool content happens to exceed the tiny cross-section
    // capacity (e.g., 0.04 dm radius × 2 dm length → 0.002 dm³ xylem, easily
    // exceeded by realistic water values).

    // Queue child for deferred removal.
    plant.queue_removal(&child);
    result.merges_performed++;
}

namespace {

// Single pass.  For each candidate (parent, child) pair where parent has
// exactly one structural child, invoke merge_pair if can_merge accepts it.
// Collects pointers first to avoid mutating the tree mid-traversal.
uint32_t compress_plant_single_pass(Plant& plant,
                                    const Genome& g,
                                    const CompressionParams& params,
                                    CompressionResult& result) {
    // Snapshot candidate pairs before mutation.
    std::vector<std::pair<Node*, Node*>> candidates;
    plant.for_each_node_mut([&](Node& parent) {
        if (!is_conduit_type(parent.type)) return;
        if (count_structural_children(parent) != 1) return;
        for (Node* child : parent.children) {
            if (child->type == parent.type && is_conduit_type(child->type)) {
                candidates.emplace_back(&parent, child);
                break;
            }
        }
    });

    uint32_t merges_this_pass = 0;
    // Track nodes that are already consumed this pass so we don't try to
    // merge a node that just got absorbed.
    std::unordered_set<Node*> consumed;
    for (auto& pair : candidates) {
        Node* p = pair.first;
        Node* c = pair.second;
        if (consumed.count(p) || consumed.count(c)) continue;
        if (!can_merge(*p, *c, g, params)) continue;
        merge_pair(plant, *p, *c, result);
        consumed.insert(c);
        merges_this_pass++;
    }
    return merges_this_pass;
}

} // namespace

CompressionResult compress_plant(Plant& plant, const CompressionParams& params) {
    CompressionResult r;
    const Genome& g = plant.genome();
    for (uint32_t pass = 0; pass < params.max_passes; ++pass) {
        uint32_t merged = compress_plant_single_pass(plant, g, params, r);
        r.passes_run++;
        plant.flush_removals();
        if (merged == 0) break;
    }
    return r;
}

} // namespace botany
