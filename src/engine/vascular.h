// src/engine/vascular.h — Xylem & phloem vascular transport.
//
// Two separate resolve passes run after the DFS tick each hour:
//
//   phloem_resolve()  — Münch pressure-driven sugar transport.
//     Three sub-steps: (3a) pre-BFS leaf loading, (3b) distance-limited BFS
//     on the stem/root/seed conduit network, (3c) post-BFS meristem unloading.
//     Sugar does NOT diffuse locally — it moves only through phloem_resolve.
//
//   xylem_resolve()   — Phase 1/Phase 2 water+cytokinin transport.
//     Unchanged from previous vascular_transport for xylem chemicals.
//     Water does NOT diffuse locally — it moves only through xylem_resolve.
//
// Local diffusion (Node::transport_with_children) handles auxin, GA, and stress only.
#pragma once

#include "engine/chemical/chemical.h"
#include "engine/world_params.h"

namespace botany {

class Plant;
struct Genome;
class Node;

// Münch pressure-flow phloem resolve — runs after the DFS tick.
void phloem_resolve(Plant& plant, const Genome& g, const WorldParams& world);

// Xylem Phase 1/Phase 2 resolve for Water and Cytokinin — runs after phloem_resolve.
void xylem_resolve(Plant& plant, const Genome& g, const WorldParams& world);

// Is this chemical a vascular chemical (skipped in local diffusion)?
// Sugar and Water are handled exclusively by phloem_resolve / xylem_resolve.
// Cytokinin is carried in xylem and also skipped on mature conduit edges.
inline bool is_vascular_chemical(ChemicalID id) {
    return id == ChemicalID::Sugar || id == ChemicalID::Water || id == ChemicalID::Cytokinin;
}

// Does this node have mature vasculature (xylem/phloem)?
// Seed always has vasculature (trunk base / junction).
bool has_vasculature(const Node& n, const Genome& g);

} // namespace botany
