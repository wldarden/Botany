#pragma once

namespace botany {

class Plant;
struct Genome;
struct WorldParams;

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

} // namespace botany
