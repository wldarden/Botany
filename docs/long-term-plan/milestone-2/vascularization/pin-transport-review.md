# PIN Transport Design — Cold Code Review

**Reviewer:** Claude Sonnet 4.6  
**Date:** 2026-04-17  
**Doc reviewed:** `pin-transport-design.md` (commit `3af71c5`)  
**Source files read:** `node.h`, `node.cpp`, `genome.h`, `vascular.cpp`

---

## Overall Assessment

The design is biologically sound and architecturally clean. The necessary infrastructure
(`auxin_flow_bias`, `last_auxin_flux`, `transport_received`, `update_canalization`) all exists
in the codebase. The structural refactoring described (removing `structural_flow_bias`, switching
cambium driver) is coherent and accounts for all 30+ occurrences in source.

Two critical bugs and several implementer-blocking ambiguities need resolution before coding
begins. The worked examples use a radius that doesn't match the actual genome default — the
parameter calibration section needs correction.

---

## What Was Verified Correct

| Claim in design doc | Code location | Status |
|---------------------|---------------|--------|
| `auxin_flow_bias` is a parent-owned map keyed by child pointer | `node.h:78` | ✓ |
| `last_auxin_flux` is a parent-owned map, records per-connection flux | `node.h:80` | ✓ |
| `transport_received` buffer exists, anti-teleports cross-node writes | `node.h:75`, flushed at `node.cpp:92-95` | ✓ |
| `update_canalization()` reads `last_auxin_flux` to update biases | `node.cpp:560-579` | ✓ |
| `initial_radius` genome param | `genome.h:83`, default `0.015 dm` | ✓ |
| Seed node is the root of the tree, identified by `parent == nullptr` | `node.h:41`, `plant.cpp:22` | ✓ |
| Tick order: vascular_transport → DFS walk | `plant.cpp:143-146` | ✓ |
| Seed type is STEM | `plant.cpp:22` | ✓ |
| Biases are cleaned up on `die()`, transferred on `replace_child()` | `node.cpp:52-60`, `node.cpp:312-313` | ✓ |

---

## `structural_flow_bias` Audit

The design deletes `structural_flow_bias` entirely. All usages found in source, against the
doc's removal plan:

| File | What it does | Removal plan in doc |
|------|-------------|---------------------|
| `node.h:79` | Map declaration on Node | Delete field |
| `node.cpp:52-55` | Transfer during `replace_child()` | Delete (structural bias doesn't transfer anymore) |
| `node.cpp:291-297` | `get_parent_structural_bias()` reads it | Delete function |
| `node.cpp:304-306` | `get_bias_multiplier()` sums flow + structural | Remove structural term; keep flow term |
| `node.cpp:312-313` | Erasure during `die()` | Delete |
| `node.cpp:572-577` | `update_canalization()` ratchet update | Delete — replace with saturation formula |
| `stem_node.cpp:33-35` | **Cambium driver** — `bias = get_parent_structural_bias()` | Switch to `auxin_flow_bias` |
| `root_node.cpp:34` | Same as stem | Same fix |
| `vascular.cpp:22,28-30` | `has_vasculature()` gate | Replace with `radius >= vascular_radius_threshold` |
| `vascular.cpp:269-270,323-324` | Conductance weight multiplier | Remove — vascular weighting by `r²` only |
| `vascular.cpp:362` | CSV debug log header | Remove or replace |
| `tissues/apical.cpp:181,189` | Initial stamp at internode creation | Remove stamp (initial radius already set) |
| `tissues/root_apical.cpp:103,109` | Same | Same |
| `genome.h` | `structural_threshold`, `structural_growth_rate`, `structural_max`, `vascular_conductance_threshold`, comment on `cambium_responsiveness` | Delete params, update comment |

**Verdict:** Design correctly accounts for all usages. No unaddressed occurrences found.

---

## Critical Issues

### 1. `last_auxin_flux` Is Cleared Before `update_canalization` Can Read PIN Flux

**Location:** `node.cpp:346` — `last_auxin_flux.clear()` is the first thing
`transport_with_children()` does.

**The bug:** PIN runs before the DFS walk. It writes per-connection flux into
`parent->last_auxin_flux[child]`. Then the DFS walk starts. When each parent node
reaches `transport_chemicals()`:

```
transport_with_children():          ← LINE 346: clears last_auxin_flux ← PIN flux erased
  ... diffusion writes last_auxin_flux ...
update_canalization():              ← reads last_auxin_flux: sees diffusion only, not PIN
```

`update_canalization()` never sees the PIN flux. `auxin_flow_bias` is driven only by
the tiny diffusion-based flux, defeating the entire canalization mechanism.

**Fix:** Move the `last_auxin_flux` clear to happen once per tick, BEFORE the PIN pass.
The cleanest approach: in `update_canalization()`, clear `last_auxin_flux` AFTER reading
it (not before writing). Remove the `last_auxin_flux.clear()` from `transport_with_children()`.
Then in `Plant::tick()`:

```
1. pre-tick: all nodes' last_auxin_flux are already clear (cleared at end of last tick)
2. vascular_transport()
3. pin_transport()        → writes last_auxin_flux
4. DFS walk:
   → transport_with_children()  → adds diffusion flux to last_auxin_flux (no clear)
   → update_canalization()      → reads combined (PIN + diffusion) flux, then clears
```

This makes PIN and diffusion flux additive in `last_auxin_flux`, which is the correct
biological model: both move auxin through the same connection, both reinforce canalization.

**Doc must specify:** where the clear moves to, and that both passes accumulate into the
same map.

---

### 2. `transient_gain` / `transient_rate` Are Not Addressed

**Location:** `node.cpp:567-570` — `update_canalization()` currently computes:

```cpp
float target = flux * g.transient_gain;          // line 568 — raw multiplier, not saturation
float& flow_bias = auxin_flow_bias[child];
flow_bias += (target - flow_bias) * g.transient_rate;   // line 570 — lerp rate
```

The design specifies a different target formula:

```
current_saturation = auxin_flux / (radius² × pin_capacity_per_area)
auxin_flow_bias = lerp(previous, current_saturation, smoothing_rate)
```

The design adds `smoothing_rate` (new param, default 0.1) and `pin_capacity_per_area`.
But it never says what to do with the **existing** `transient_gain` (default 2.0) and
`transient_rate` (default 0.2). An implementer encounters three questions:

- Is `smoothing_rate` a new param, or is it a rename of `transient_rate`?
- Is `transient_gain` deleted?
- If `transient_rate` is renamed and its default changes from 0.2 to 0.1, that's a behavioral change in any running plant.

**What the implementation must do:**
- **Delete** `transient_gain` — the saturation formula replaces the raw multiplier
- **Rename** `transient_rate` → `smoothing_rate`, change default from 0.2 → 0.1
- Update `update_canalization()` to use the saturation formula with `pin_capacity_per_area`
  and the child node's `radius` (accessed via the child pointer)

**Add to Implementation Step 2:** `Delete transient_gain. Rename transient_rate → smoothing_rate,
change default from 0.2 → 0.1.`

---

## Major Issues

### 3. Worked Examples Use Wrong Radius

**Location:** "Effective Reach" section in New Genome Parameters.

The doc says `r = initial_radius ≈ 0.05 dm`, but `genome.h:83` sets `initial_radius = 0.015 dm`.
The actual numbers for a fresh internode:

```
r  = 0.015 dm
r² = 0.000225 dm²
max_capacity = 0.000225 × 100 = 0.0225 AU/tick  (at full efficiency)
base capacity = 0.0225 × 0.2  = 0.0045 AU/tick  (at pin_base_efficiency = 0.2)
```

The doc says `max_capacity = 0.25` and `base = 0.05` — both off by 11×. The calibration
target (r² ≈ 0.005 dm²) corresponds to `r ≈ 0.071 dm`, which is ~5× initial_radius — a
moderately thickened stem, not "a thin stem."

**Actual steady-state reach with initial_radius = 0.015:**

With max_capacity = 0.0225 and decay = 0.12 per tick:
```
Node 0 (apex):   level ≈ 1.06 AU  (production 0.15 – outflow 0.0225 = decay-balanced)
Node 1:          level ≈ 0.020
Node 5:          level ≈ 0.012
Node 10:         level ≈ 0.006
```

This reaches node 10 with ~0.006 AU — below the `auxin_threshold = 0.15` for apical
dominance. The gradient forms but is too shallow to suppress lateral buds beyond a few nodes.

**Implication:** `pin_capacity_per_area = 100` may be too low for a plant at initial
radius. Consider starting at 500–1000 to give thin stems sufficient throughput before
canalization builds. Alternatively, raise `pin_base_efficiency` (0.2 → 0.5). The
calibration note in the doc should use actual initial_radius numbers, not ≈ 0.05.

---

### 4. `vascular_radius_threshold` Default Not Specified

The doc correctly says to replace the `structural_flow_bias`-based `vascular_conductance_threshold`
with a radius-based `vascular_radius_threshold` in `has_vasculature()`. But the default
value is never given.

**Suggested default:** `0.012 dm` — slightly below `initial_radius = 0.015` so fresh
internodes qualify for vasculature immediately, while degenerate zero-radius nodes (if any)
are excluded. Could also simply be `initial_radius × 0.9` to track the initial_radius
setting.

**Specify the default in Implementation Step 2** and add it to the doc's genome params
table (or note it's derived from initial_radius).

---

## Moderate Issues

### 5. Seed Junction Timing Ambiguity

The design says the seed distributes auxin to root children from `chemical(ChemicalID::Auxin)`
after "flushing its transport_received buffer from shoot-side children." But `transport_received`
is normally flushed inside `Node::tick()` during the DFS walk — which hasn't happened yet
during the PIN pass.

**What actually happens:** During the shoot-side PIN pass, shoot children write into
`seed->transport_received[Auxin]`. This is not yet in `seed->chemical(Auxin)`. When the
root-distribution step runs during the same PIN pass, `chemical(Auxin)` contains only the
**previous tick's** seed auxin, not the freshly received shoot flow.

**Impact:** One-tick pipeline delay in the root-side signal. Shoot auxin arrives at the
seed's `transport_received` this tick; flows to root children next tick (after DFS flushing).
This is a natural one-tick lag, not a bug — but the doc's phrasing ("after flushing its
transport_received buffer") implies an immediate in-pass flush that doesn't happen through
the normal mechanism.

**Clarification needed:** Either:
a) Acknowledge the one-tick delay (probably fine — at 10–30 tick response timescales it
   doesn't matter), OR
b) Specify that `pin_transport()` explicitly flushes `seed->transport_received` mid-pass
   before distributing to root children.

Option (a) is simpler. The doc should say: root distribution uses the seed's auxin from
the previous tick. The pipeline delay is one tick and is negligible.

---

## Minor Issues

### 6. Feedback Chain Shows Non-Lerped Formula

In the "full feedback chain" diagram:

```
→ auxin_flow_bias = flux / (radius² × pin_capacity_per_area)   [PIN saturation, 0→1]
```

This shows `auxin_flow_bias` as the raw saturation, but the actual update is:
```
auxin_flow_bias = lerp(previous, saturation, smoothing_rate)
```

This is the steady-state value (when lerp has converged), so it's correct at equilibrium
but may confuse an implementer expecting the exact update formula. Consider labeling it
`→ current_saturation ≈` or adding a note.

---

## Parameter Sanity Check

With `initial_radius = 0.015 dm`, `pin_capacity_per_area = 100`, `pin_base_efficiency = 0.2`:

| Scenario | r (dm) | r² (dm²) | max_capacity (AU/tick) | base capacity (AU/tick) |
|----------|--------|----------|----------------------|------------------------|
| Fresh internode (actual initial_radius) | 0.015 | 0.000225 | **0.023** | **0.0045** |
| Doc example ("thin stem") | 0.05 | 0.0025 | 0.25 | 0.050 |
| Doc calibration target | 0.071 | 0.005 | 0.50 | 0.100 |
| Established trunk (5× initial) | 0.075 | 0.005625 | **0.56** | **0.11** |

The doc examples are accurate for a moderately thickened stem but do not represent the
cold-start case. The established trunk (r = 0.075) can trivially pass all available auxin
— this is correct. The bootstrapping period (r ≤ 0.03) is where the numbers need tuning.

**Recommendation:** Raise `pin_capacity_per_area` to 500 (or `pin_base_efficiency` to 0.5)
as the initial default, then tune down empirically. With cap = 500 and initial_radius = 0.015:
```
base capacity = 0.000225 × 500 × 0.2 = 0.0225 AU/tick
```
— 15% of apical production (0.15), enough to establish the initial gradient.

---

## Implementation Additions Required

Additions to the existing steps based on this review:

**Step 1 (Files):** No change.

**Step 2 (Genome params):** Add:
- Delete `transient_gain`
- Rename `transient_rate` → `smoothing_rate`, change default from 0.2 → 0.1
- Add `vascular_radius_threshold` (suggested default: `0.012 dm`)

**Step 3 (Wiring):** No change.

**Step 4 (Regression test):** Add:
- Assert `auxin_flow_bias` is driven by PIN flux (not just diffusion): build a plant
  where PIN is active but diffusion is effectively blocked (e.g., very short chain,
  measure `auxin_flow_bias` at the connection between node 0 and node 5 vs expected
  saturation from PIN flux)
- Assert `last_auxin_flux` contains non-zero values after both PIN and diffusion (verify
  additive accumulation, not overwrite)

**Step 5 (CLAUDE.md):** Also update:
- Note `transient_gain` deleted, `transient_rate` renamed
- Note `vascular_radius_threshold` added

**New: Step 0 — Pre-implementation:** Before writing any PIN code, change the
`last_auxin_flux` clear point: move from `transport_with_children()` line 346 to
the end of `update_canalization()` (after reading, before returning). This is a
prerequisite for PIN flux to survive until `update_canalization` reads it.

---

## Verdict

The design is implementable. Fix the two critical issues (last_auxin_flux clear timing,
transient_gain/transient_rate resolution) and the three major issues (example radius,
vascular_radius_threshold default, seed timing clarification) before writing code. Everything
else (structural_flow_bias removal, cambium driver switch, tick order) is correctly specified.
