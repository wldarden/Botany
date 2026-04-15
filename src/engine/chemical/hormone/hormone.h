// src/engine/chemical/hormone/hormone.h
#pragma once

#include "engine/chemical/chemical.h"
#include <algorithm>

namespace botany {

// Per-chemical diffusion parameters extracted from a plant's genome.
struct ChemicalDiffusionParams {
    ChemicalID id;
    float diffusion_rate;   // gradient responsiveness
    float decay_rate;       // per-tick exponential decay
    float bias;             // signed equilibrium shift (neg = toward root, pos = toward tips)
    float base_transport;   // throughput floor (units per tick)
    float transport_scale;  // radius amplification factor for throughput
};

// Compute the desired transport flow between a parent and one child.
// Returns a weight: positive = child wants to receive from parent,
//                   negative = child wants to give to parent.
// Pure function — does not mutate any values.
inline float compute_transport_flow(
    float child_val, float parent_val,
    float child_cap, float parent_cap,
    float child_radius, float parent_radius,
    float reference_radius,
    const ChemicalDiffusionParams& params)
{
    bool has_cap = (child_cap > 0.0f && parent_cap > 0.0f);

    // Concentrations
    float child_conc  = has_cap ? child_val / child_cap   : child_val;
    float parent_conc = has_cap ? parent_val / parent_cap : parent_val;

    // Shifted equilibrium: bias > 0 pushes toward tips, < 0 toward root.
    // (parent_conc - child_conc) is positive when parent has more → child should receive.
    // Subtracting bias: positive bias makes it easier for child to receive (acropetal).
    float diff = (parent_conc - child_conc) - params.bias;
    float desired = diff * params.diffusion_rate;

    // Scale back to absolute units for capped chemicals
    if (has_cap) {
        desired *= (child_cap + parent_cap) * 0.5f;
    }

    // Radius scaling: thicker connections transport more
    float min_r = reference_radius * 0.2f;
    float conn_r = std::max(std::min(child_radius, parent_radius), min_r);
    float radius_factor = conn_r / std::max(reference_radius, 1e-6f);
    desired *= radius_factor;

    // Throughput cap
    float max_transport = params.base_transport + radius_factor * params.transport_scale;
    return std::clamp(desired, -max_transport, max_transport);
}

} // namespace botany
