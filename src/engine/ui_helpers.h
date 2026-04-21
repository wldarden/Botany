#pragma once

#include "engine/chemical/chemical.h"
#include "engine/compartments.h"

namespace botany {

class Node;
class Genome;
struct WorldParams;

// Returns the conduit pool this chemical lives in for this node.
// - Stem/Root with vascular radius: own phloem (Sugar) or own xylem (Water, Cytokinin).
// - Leaf/Apical/RootApical/seed specialty: walks up to the nearest ancestor stem/root.
// - Signaling chems (Auxin, Gibberellin, Ethylene, Stress): returns nullptr (no conduit).
// - If no upstream conduit exists (e.g., seedling before vascular admission): nullptr.
const TransportPool* vascular_scope(const Node& n, ChemicalID chem);

// Preview of what pay_maintenance() would deduct this tick, without side effects.
// Used by the starvation overlay to evaluate sugar coverage.
// NOTE: this is stubbed for Phase 1; Task 18 completes with full per-type formulas.
float compute_maintenance_cost(const Node& n, const Genome& g, const WorldParams& w);

// 1 - radial_permeability_sugar(r) / base_radial_permeability_sugar.
// 0% = young leaky stem. ~90% = fully closed trunk (asymptote at the floor fraction).
// Defined for stems/roots; returns 0 for other types.
float hydraulic_maturity(const Node& n, const Genome& g);

// Walk-up count of parent hops until parent == nullptr (seed has 0).
int nodes_to_seed(const Node& n);

} // namespace botany
