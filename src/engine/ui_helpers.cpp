#include "engine/ui_helpers.h"
#include "engine/node/node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/vascular_sub_stepped.h"
#include <glm/geometric.hpp>
#include <cmath>

namespace botany {

const TransportPool* vascular_scope(const Node& n, ChemicalID chem) {
    switch (chem) {
        case ChemicalID::Sugar:
            return n.nearest_phloem_upstream();
        case ChemicalID::Water:
        case ChemicalID::Cytokinin:
            return n.nearest_xylem_upstream();
        case ChemicalID::Auxin:
        case ChemicalID::Gibberellin:
        case ChemicalID::Ethylene:
        case ChemicalID::Stress:
        case ChemicalID::Count:
            return nullptr;
    }
    return nullptr;
}

float compute_maintenance_cost(const Node& n, const Genome& /*g*/, const WorldParams& w) {
    // Preview of what pay_maintenance() (i.e., maintenance_cost() virtual call) would
    // return for this node.  Must match each subclass's maintenance_cost() exactly.
    // Update here whenever those overrides change.

    if (n.as_stem()) {
        // StemNode::maintenance_cost — living ring scaled by lateral surface area (πrL)
        float length = std::max(glm::length(n.offset), 0.01f);
        return w.sugar_maintenance_stem * 3.14159f * n.radius * length;
    }
    if (n.as_root()) {
        // RootNode::maintenance_cost — same model as stem, root-tissue rate
        float length = std::max(glm::length(n.offset), 0.01f);
        return w.sugar_maintenance_root * 3.14159f * n.radius * length;
    }
    if (auto* lf = n.as_leaf()) {
        // LeafNode::maintenance_cost — scales with leaf area (leaf_size²)
        return w.sugar_maintenance_leaf * lf->leaf_size * lf->leaf_size;
    }
    if (auto* ap = n.as_apical()) {
        // ApicalNode::maintenance_cost — active meristems pay flat rate; dormant pay 0
        return ap->active ? w.sugar_maintenance_meristem : 0.0f;
    }
    if (auto* ra = n.as_root_apical()) {
        // RootApicalNode::maintenance_cost — same as apical
        return ra->active ? w.sugar_maintenance_meristem : 0.0f;
    }
    // Base Node::maintenance_cost returns 0 — seed node falls here
    return 0.0f;
}

float hydraulic_maturity(const Node& n, const Genome& g) {
    if (!n.as_stem() && !n.as_root()) return 0.0f;
    float base = g.base_radial_permeability_sugar;
    if (base <= 1e-9f) return 0.0f;
    float perm = radial_permeability_sugar(n.radius, g);
    float closed = 1.0f - (perm / base);
    if (closed < 0.0f) closed = 0.0f;
    if (closed > 1.0f) closed = 1.0f;
    return closed;
}

int nodes_to_seed(const Node& n) {
    int count = 0;
    const Node* p = n.parent;
    while (p) { ++count; p = p->parent; }
    return count;
}

} // namespace botany
