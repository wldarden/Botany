#pragma once
#include <cstdint>

namespace botany {

class Plant;
class Node;
struct Genome;
struct WorldParams;

// Per-node supply/demand snapshot for one tick.  Computed once at the start
// of vascular_sub_stepped() and then amortized across N sub-steps.
struct VascularBudget {
    float sugar_supply     = 0.0f;  // leaf: above reserve; 0 otherwise
    float sugar_demand     = 0.0f;  // active meristem: to fill toward sink_target
    float water_demand     = 0.0f;  // leaf/meristem: to fill toward turgor_target
    float cytokinin_supply = 0.0f;  // root apical: produced this tick
};

VascularBudget compute_budget(Node& n, const Genome& g, const WorldParams& world);

// Source injection: transfer budget/N from the node's own local_env into its
// nearest upstream conduit.  Only does anything when the node has a
// sugar_supply or cytokinin_supply budget (leaves, root apicals).
void inject_step(Node& source, const VascularBudget& b, uint32_t N, const Genome& g);

// Sink extraction: pull budget/N from nearest upstream conduit into sink's
// local_env.  Capped by actual pool content so the pool never goes negative.
// Water extraction also carries cytokinin proportionally (cytokinin rides in
// xylem solution).
void extract_step(Node& sink, const VascularBudget& b, uint32_t N, const Genome& g);

// Sub-stepped vascular transport.  Replaces the pairwise-Jacobi
// vascular_transport in Phase E.  For the duration of Phase D this function
// coexists with the old one but is not wired into Plant::tick() yet.
//
// Algorithm: N = world.vascular_substeps iterations of
//   1. Inject — sources push budget/N into parent's conduit
//   2. Radial flow — stem/root local_env ⇄ own phloem/xylem, radius-dependent
//      permeability
//   3. Extract — sinks pull budget/N from parent's conduit
//   4. Longitudinal Jacobi — one pass of neighbor pressure equalization
// See spec: docs/superpowers/specs/2026-04-19-compartmented-vascular-model-design.md
void vascular_sub_stepped(Plant& plant, const Genome& g, const WorldParams& world);

// Radial permeability between a stem/root's own local_env and its own
// phloem/xylem.  Curve: perm(r) = base × (floor + (1 - floor) / (1 + (r/r_half)²))
// Young thin stems (r=0): perm ≈ base.  Mature thick trunks (r→∞): perm ≈ base × floor.
// See spec section 6.
float radial_permeability_sugar(float radius, const Genome& g);
float radial_permeability_water(float radius, const Genome& g);

// Phloem/xylem capacity for a stem or root node: π · r² · length · fraction.
// Returns 0 for nodes without the matching pool (leaves, meristems).
float phloem_capacity(const Node& n, const Genome& g);
float xylem_capacity(const Node& n, const Genome& g);

} // namespace botany
