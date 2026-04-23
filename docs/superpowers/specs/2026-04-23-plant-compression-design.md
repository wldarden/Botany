# Plant Compression (Node Fusion) — Design

**Date:** 2026-04-23
**Status:** Draft (approved in brainstorming)

## Problem

A 4-year sim produces ~4500 nodes for a mature tree, but a lot of those are structurally redundant: long trunks consist of dozens of adjacent, near-collinear, similar-radius STEM nodes that contribute metabolic bookkeeping without meaningful architectural information. Per-tick vascular sub-stepping scales with conduit depth, so redundant nodes cost CPU time every tick for no benefit.

## Goal

Periodically (auto or on-demand) merge adjacent stem/root nodes that pass a structural similarity test, reducing node count while preserving plant topology, chemical mass, and canalization history within acceptable tolerance.

## Non-goals

- Reversibility (once merged, two nodes cannot be unmerged).
- Perfect conservation (small cap-clamping loss is acceptable; see "Decisions locked").
- Merging across node types (stem↔root, or stem↔leaf) — only same-type pairs.
- Merging inside `Plant::tick` or any sub-phase (vascular, metabolism). Compression runs between full ticks only.

## Decisions locked in brainstorming

1. **Timing:** runs entirely outside `Engine::tick()`. Either between full ticks (auto mode) or on demand (manual button). Never mid-tick.
2. **Lossiness:** lossy is acceptable. When the volume-preserving merged cap is smaller than summed contents, clamp and log the delta. A per-merge loss of < 1% is the expected worst case.
3. **UI:** realtime viewer gets a Compression panel with (a) auto-compress checkbox, (b) interval input (default 1000 ticks) visible when auto is on, (c) manual "Compress Now" button visible when auto is off, (d) last-run stats (merges count, sugar/water Δ).
4. **Scope:** STEM-STEM and ROOT-ROOT merges only. Seed is never merged.

## Modules

### `src/engine/compression.{h,cpp}`

Public API:

```cpp
struct CompressionParams {
    float max_angle_rad       = 0.175f;  // 10°
    float max_radius_ratio    = 0.20f;   // |r1-r2|/max(r1,r2)
    float max_combined_length = 0.0f;    // 0 = auto (2 × g.max_internode_length)
    uint32_t max_passes       = 5;       // iterate until no more merges or cap reached
};

struct CompressionResult {
    uint32_t merges_performed = 0;
    uint32_t passes_run       = 0;
    float    delta_sugar      = 0.0f;    // negative = clamped loss
    float    delta_water      = 0.0f;
    float    delta_auxin      = 0.0f;
    float    delta_cytokinin  = 0.0f;
};

CompressionResult compress_plant(Plant& plant, const CompressionParams& params);
```

Deliberately NOT a member of `Plant` — it's an off-tick meta-operation, not a per-tick plant behavior. Keeps `plant.cpp` focused on per-tick work.

### `src/engine/engine.h` / `engine.cpp` additions

```cpp
// Enable/disable auto-compression triggered inside Engine::tick().
void enable_autocompress(bool enabled);
void set_compression_interval(uint32_t ticks);
void set_compression_params(const CompressionParams& params);

// Trigger a one-off compression pass (manual button, or external CLI tool).
CompressionResult trigger_compression();

// Last run's result, for UI display.
const CompressionResult& last_compression() const;
```

Internal state on `Engine`:
```cpp
bool compression_enabled_ = false;
uint32_t compression_interval_ = 1000;
CompressionParams compression_params_;
CompressionResult last_compression_;
```

In `Engine::tick()` (end of the method, after `plants_[0]->tick(...)`, before incrementing `tick_`):
```cpp
if (compression_enabled_ && tick_ > 0 && (tick_ % compression_interval_) == 0) {
    if (!plants_.empty()) {
        last_compression_ = compress_plant(*plants_[0], compression_params_);
    }
}
```

### `src/app_realtime.cpp` additions

New ImGui collapsing header "Compression" inside the Controls window. Contains the checkbox, conditional controls, and last-run stats as described in Design §3.

## Merge algorithm

### Candidate pair selection

DFS from the seed. For each parent P visited, examine its structural children list (children whose type is STEM or ROOT — leaves, apicals, root apicals are ignored for the "single structural child" test). Attempt merge of P with its child C if P has exactly one structural child.

### Merge acceptance criteria

Merge P + C is permitted if and only if ALL hold:

1. `P.type == C.type` and type is STEM or ROOT.
2. P has a parent (P is not the seed).
3. P has exactly one structural (STEM/ROOT) child, namely C. Leaves, apicals, and root apicals on P do not count against this rule.
4. Angle between P's `offset` and C's `offset` (both normalized) satisfies `dot(normalize(P.offset), normalize(C.offset)) >= cos(params.max_angle_rad)`.
5. `abs(P.radius - C.radius) / max(P.radius, C.radius) < params.max_radius_ratio`.
6. Both P and C have `age >= g.internode_maturation_ticks` (no active elongation).
7. `length(P) + length(C) < effective_max_length`, where `effective_max_length = params.max_combined_length > 0 ? params.max_combined_length : 2.0f * plant.genome().max_internode_length`.

Where `length(N) = glm::length(N.offset)`.

The sentinel `max_combined_length = 0` means "auto from genome" — compress_plant reads the plant's genome to pick a sensible default. Explicit positive values override.

### Merge execution (P absorbs C, C is queued for removal)

Given an accepted pair (P, C):

1. **Compute new geometry.**
   - `r_merged = sqrt((P.radius² · length(P) + C.radius² · length(C)) / (length(P) + length(C)))` — volume-preserving.
   - `new_offset = P.offset + C.offset` (grandparent → P's new position now matches what was grandparent → C's old position; i.e., P absorbs C's geometry).
   - `new_rest_offset = P.rest_offset + C.rest_offset`.

2. **Sum chemicals.**
   - For each (id, value) in `C.local().chemicals`: `P.local().chemical(id) += value`.
   - If both are STEM or both are ROOT (they are), sum phloem: `P.phloem()->chemical(id) += C.phloem()->chemical(id)` for every id in C's phloem. Same for xylem.

3. **Update P's geometry** — after summing chemicals, assign `P.offset = new_offset`, `P.rest_offset = new_rest_offset`, `P.radius = r_merged`. Position syncs automatically on next tick via `sync_world_position`.

4. **Reparent C's children.**
   - For each `child` in `C.children`:
     - `child.parent = &P`.
     - Recompute `child.offset` so `child.position` stays at the same world location: `child.offset = child.position - P.new_world_position`. Since P's new world position equals C's old world position, `child.offset = child.position - C.old_position = child.old_offset`. So offsets for C's children DO NOT change (they were already expressed relative to C's position, which is where P now sits).
   - `P.children.push_back(child)` for each child.
   - Clear `C.children` to prevent `die()` from cascading destruction.

5. **Move canalization bias entries.**
   - For each `(child_ptr, bias)` in `C.auxin_flow_bias`: `P.auxin_flow_bias[child_ptr] = bias`. (The child_ptr is still valid — children were reparented, not reallocated.)
   - `P.auxin_flow_bias.erase(&C)` — the P→C edge no longer exists.

6. **Clamp to new caps, log deltas.**
   - Compute `sugar_cap(P_post_merge)`, `water_cap(P_post_merge)`, and the vascular caps.
   - For each local chemical with a cap: if value > cap, record `(value - cap)` as a loss delta for that chemical and clamp to cap.
   - Same for phloem/xylem pool entries where caps apply.

7. **Remove C.** `plant.queue_removal(&C)`. Flushed after the full scan.

### Multi-pass convergence

A single DFS scan catches only non-overlapping pairs. After one pass, new adjacencies emerge (e.g., merged-P and its grandchild may now be mergeable). Loop up to `params.max_passes` (default 5) calling a single-pass scan until `merges_performed == 0` in a pass or the cap is hit.

Each pass calls `plant.flush_removals()` before the next pass to keep pointer state consistent.

## UI (realtime viewer)

Inside the existing Controls window, new `ImGui::CollapsingHeader("Compression")` containing:

- `ImGui::Checkbox("Auto-compress", &autocompress)` — toggles `Engine::enable_autocompress`.
- If autocompress is true: `ImGui::InputInt("Every N ticks", &interval)` — calls `set_compression_interval` on change.
- If autocompress is false: `ImGui::Button("Compress Now")` — calls `trigger_compression`.
- `ImGui::Separator();`
- `ImGui::Text("Last run: %u merges across %u passes", result.merges_performed, result.passes_run);`
- If any `delta_*` is nonzero: show them as formatted values. If all zero: "no loss".

## Testing

`tests/test_compression.cpp`:

1. **Sugar conservation within tolerance** — grow a 40-tick plant, record total sugar, call compress_plant, record total sugar. Assert `(pre - post) >= 0` and `(pre - post) < 0.01 * pre`.
2. **Angle gate refuses bent pair** — construct a plant where P's offset = (0,1,0) and C's offset = (0.9, 0.4, 0) (large angle). Verify no merge happens at that pair.
3. **Branching gate refuses mergeable-looking parent with 2 kids** — P has 2 structural STEM children. Verify no merge on P.
4. **Leaf attached to intermediate stem survives with same world position** — stem A → stem B (with leaf L attached via B.add_child(L)) → stem C. Merge B into A. Assert L.parent == A and glm::distance(L.position_after_next_tick, L.position_before_merge) < epsilon. (Position sync is deferred to next tick's sync_world_position; test either invokes that or compares L.offset reconstruction.)
5. **Axillary bud survives** — same pattern with a dormant ApicalNode.
6. **Canalization bias carries** — construct a chain A → B → C. Set `B.auxin_flow_bias[&C] = 0.73` before merge. Merge B into A (A absorbs B). After merge, assert `A.auxin_flow_bias[&C]` exists with value 0.73. (Rationale: the B→C edge becomes the A→C edge post-merge, and its canalization history should transfer with it.)
7. **Idempotence** — compress_plant twice in a row. Second call returns `merges_performed == 0`.
8. **Volume-preserving radius** — construct a specific pair with r=0.05, L=1 and r=0.03, L=1. After merge, assert radius within 1e-4 of `sqrt((0.0025 + 0.0009) / 2) = sqrt(0.0017) ≈ 0.04123`.
9. **Save/load round-trip post-compression** — grow, compress, save, load, verify node count matches post-compression count.

Manual smoke:
- 4000+ tick run.
- Hit "Compress Now".
- Observe node count drop.
- Continue ticking; no visual popping, no crashes, plant continues growing sensibly.
- Save snapshot, quit, load, continue.

## Open risks

- **Edge conductance shifts.** Merging two vascular edges in series changes effective conductance. In practice the merged node's cap approximates the union of the original two nodes' caps, so the Jacobi equalizer should behave similarly, but if whole-plant vascular throughput shifts noticeably after a big compression pass, that's a flag. Manual smoke + tick counter logs should catch this.
- **Canalization cliff on the retired P→C edge.** The flux history encoded in `P.auxin_flow_bias[&C]` is destroyed. The merged node has no equivalent of that bias (it was the edge's history, not a node's). This is architecturally unavoidable when fusing an edge into a node. The main trunk's "strongest path" invariant should be unaffected since the adjacent edges that mattered for trunk canalization (grandparent→P and C→grandchild) survive, with their biases intact.
- **Lossiness accumulates.** Repeated compressions over many sim-years could drop total plant mass a few percent from cap clamping. If this becomes observable, revisit the volume-preserving radius formula or switch to strict mode.
- **Tests 2 and 3 rely on hand-constructed plants.** Those plant setups have historically been fragile in this codebase (vascular structure needs specific field init). Writing them may require a helper or extending `Plant::install_node` usage.
