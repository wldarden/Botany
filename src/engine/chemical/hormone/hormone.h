// src/engine/chemical/hormone/hormone.h
#pragma once

#include "engine/chemical/chemical.h"
#include <algorithm>

namespace botany {

// Per-chemical diffusion parameters extracted from a plant's genome.
struct ChemicalDiffusionParams {
    ChemicalID id;
    float diffusion_rate;
    float decay_rate;
};

// Single diffusion function used by all chemicals.
// High flows to low. Rate controls speed. Clamp so nothing goes negative.
inline void diffuse(float& my_val, float& neighbor_val, float rate) {
    float flow = (my_val - neighbor_val) * rate;
    if (flow > 0.0f) flow = std::min(flow, my_val);
    else             flow = std::max(flow, -neighbor_val);
    my_val -= flow;
    neighbor_val += flow;
}

} // namespace botany
