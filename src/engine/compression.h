#pragma once

#include <cstdint>

namespace botany {

class Plant;
class Node;
struct Genome;

// Tunable thresholds for compression.  See design doc
// docs/superpowers/specs/2026-04-23-plant-compression-design.md.
struct CompressionParams {
    float    max_angle_rad       = 0.175f;  // ≈ 10°
    float    max_radius_ratio    = 0.20f;   // |r1-r2|/max(r1,r2)
    float    max_combined_length = 0.0f;    // 0 → auto = 2 × g.max_internode_length
    uint32_t max_passes          = 5;       // multi-pass convergence cap
};

// Summary of a compress_plant run.  delta_* are negative when cap clamping
// discards chemicals; zero if the merge was mass-conservative.
struct CompressionResult {
    uint32_t merges_performed = 0;
    uint32_t passes_run       = 0;
    float    delta_sugar      = 0.0f;
    float    delta_water      = 0.0f;
    float    delta_auxin      = 0.0f;
    float    delta_cytokinin  = 0.0f;
};

// Predicate: can parent P absorb its single structural child C under params?
// Used by compress_plant's scan and exposed here for fine-grained testing.
bool can_merge(const Node& parent, const Node& child,
               const Genome& g, const CompressionParams& params);

// Main entry.  Scans plant and merges every accepted adjacent (parent, child)
// stem/root pair in up to params.max_passes iterations, flushing deferred
// removals between passes.  Runs mutations directly on the plant; must only
// be called between full ticks (never from inside Plant::tick / Engine::tick).
CompressionResult compress_plant(Plant& plant, const CompressionParams& params);

} // namespace botany
