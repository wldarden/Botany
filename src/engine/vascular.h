// src/engine/vascular.h — Xylem & phloem vascular transport.
// Global pass that moves sugar (phloem) and water+cytokinin (xylem)
// via bulk flow, bypassing intermediate conduit nodes.
// Local diffusion (Node::transport_with_children) handles everything else.
#pragma once

#include "engine/chemical/chemical.h"
#include "engine/world_params.h"

namespace botany {

class Plant;
struct Genome;
class Node;

// Run vascular transport for the entire plant.
// Call once per tick, before the DFS tree walk.
// Pass world so the call can gate optional debug logging via world.vascular_debug_log.
void vascular_transport(Plant& plant, const Genome& g, const WorldParams& world);

// Münch pressure-flow phloem resolve.
// Three passes: leaf loading (pre-Jacobi), Jacobi pipe-network resolution,
// meristem unloading (post-Jacobi).  Called from vascular_transport().
void phloem_resolve(Plant& plant, const Genome& g, const WorldParams& world);

// Is this chemical transported via vasculature (bulk flow)?
inline bool is_vascular_chemical(ChemicalID id) {
    return id == ChemicalID::Sugar || id == ChemicalID::Water || id == ChemicalID::Cytokinin;
}

// Does this node have mature vasculature (xylem/phloem)?
// Seed always has vasculature (trunk base / junction).
bool has_vasculature(const Node& n, const Genome& g);

} // namespace botany
