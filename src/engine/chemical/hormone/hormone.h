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

// Unified transport: concentration-based with shifted equilibrium and throughput cap.
//
// my_cap / parent_cap: capacity for this chemical.
//   > 0 = real capacity (sugar): concentration = level/cap, destination headroom clamped.
//   <= 0 = no capacity (hormones): concentration = raw level, no headroom clamp.
//
// Flow sign: positive = toward parent (root-ward), negative = toward child (tip-ward).
inline void transport_chemical(
    float& my_val, float& parent_val,
    float my_cap, float parent_cap,
    float node_radius, float parent_radius,
    float reference_radius,
    const ChemicalDiffusionParams& params)
{
    bool has_cap = (my_cap > 0.0f && parent_cap > 0.0f);

    // Concentrations
    float my_conc, parent_conc;
    if (has_cap) {
        my_conc    = my_val / my_cap;
        parent_conc = parent_val / parent_cap;
    } else {
        my_conc    = my_val;
        parent_conc = parent_val;
    }

    // Shifted equilibrium: bias offsets where "equal" is
    float effective_diff = (my_conc - parent_conc) - params.bias;

    // Desired flow (positive = toward parent/root)
    float desired_flow = effective_diff * params.diffusion_rate;

    // Scale concentration gradient back to absolute units
    if (has_cap) {
        float avg_cap = (my_cap + parent_cap) * 0.5f;
        desired_flow *= avg_cap;
    }

    // Radius scaling: thicker connections transport more
    float min_r = reference_radius * 0.2f;
    float conn_r = std::max(std::min(node_radius, parent_radius), min_r);
    float radius_factor = conn_r / std::max(reference_radius, 1e-6f);
    desired_flow *= radius_factor;

    // Throughput cap: bottlenecked by thinner pipe
    float max_transport = params.base_transport + radius_factor * params.transport_scale;
    float flow = std::clamp(desired_flow, -max_transport, max_transport);

    // Source clamp: can't send more than you have
    if (flow > 0.0f) {
        flow = std::min(flow, my_val);
        if (has_cap) flow = std::min(flow, std::max(0.0f, parent_cap - parent_val));
    } else {
        flow = std::max(flow, -parent_val);
        if (has_cap) flow = std::max(flow, -std::max(0.0f, my_cap - my_val));
    }

    my_val    -= flow;
    parent_val += flow;
}

} // namespace botany
