#pragma once

#include <unordered_map>
#include "engine/chemical/chemical.h"

namespace botany {

// Per-node local compartment.  Every node owns one.  Holds the chemicals the
// node itself uses for metabolism, growth, and signaling — i.e. everything
// that is NOT in the long-distance transport stream.
struct LocalEnv {
    std::unordered_map<ChemicalID, float> chemicals;

    float& chemical(ChemicalID id) { return chemicals[id]; }
    float chemical(ChemicalID id) const {
        auto it = chemicals.find(id);
        return it != chemicals.end() ? it->second : 0.0f;
    }
};

// Per-stem/root vascular conduit.  A StemNode or RootNode owns one phloem and
// one xylem TransportPool — representing the sieve tubes (phloem) and vessel
// elements (xylem) in that segment.  Specialty nodes (leaves, meristems) do
// not own any TransportPool.
struct TransportPool {
    std::unordered_map<ChemicalID, float> chemicals;

    float& chemical(ChemicalID id) { return chemicals[id]; }
    float chemical(ChemicalID id) const {
        auto it = chemicals.find(id);
        return it != chemicals.end() ? it->second : 0.0f;
    }
};

} // namespace botany
