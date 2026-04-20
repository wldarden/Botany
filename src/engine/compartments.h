#pragma once

#include <unordered_map>
#include "engine/chemical/chemical.h"

// Hash for ChemicalID (duplicate of node.h's definition, but compartments.h
// needs to be includable standalone by headers/tests that don't pull in node.h)
namespace std {
    template<>
    struct hash<botany::ChemicalID> {
        std::size_t operator()(botany::ChemicalID id) const noexcept {
            return static_cast<std::size_t>(id);
        }
    };
}

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
