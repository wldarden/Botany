#include "engine/ui_helpers.h"
#include "engine/node/node.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/vascular_sub_stepped.h"

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

float compute_maintenance_cost(const Node& n, const Genome& g, const WorldParams& w) {
    // Phase 1 stub — Task 18 replaces this with exact per-type formulas cross-checked
    // against each subclass's pay_maintenance() override. For now return 0 so the
    // starvation overlay gets a sensible "no demand" reading until then.
    (void)n; (void)g; (void)w;
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
