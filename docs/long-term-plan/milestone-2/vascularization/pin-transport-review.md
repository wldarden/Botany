# PIN Transport Design — Cold Review

**Doc reviewed:** `pin-transport-design.md`
**Source files read:** `node.h`, `node.cpp`, `genome.h`, `vascular.cpp`, `stem_node.cpp`,
`root_node.cpp`, `apical.cpp`, `root_apical.cpp`; grepped tests for `structural_flow_bias`

---

## Overall Assessment

The design is biologically sound and architecturally clean. The necessary infrastructure
(`auxin_flow_bias`, `last_auxin_flux`, `transport_received`, `update_canalization`) all exists
in the codebase. The structural refactoring (removing `structural_flow_bias`, switching cambium
driver to `auxin_flow_bias`) is coherent.

**Not ready to implement as written.** Three issues would produce silent wrong behavior:

1. `last_auxin_flux.clear()` erases PIN flux before `update_canalization` can read it.
2. Seed junction distributes from `chemical(Auxin)` but shoot-side auxin is in
   `transport_received` during the PIN pass — one tick stale.
3. Parameter calibration uses `r ≈ 0.05 dm` but actual `initial_radius = 0.015 dm`
   (off by 11× in throughput examples).

Fixing these three and filling six completeness gaps (described below) would make it
ready to implement.

---

## What Was Verified Correct

| Claim in design doc | Code location | Status |
|---------------------|---------------|--------|
| `auxin_flow_bias` is parent-owned, keyed by child pointer | `node.h:78` | ✓ |
| `last_auxin_flux` is parent-owned, records per-connection flux | `node.h:80` | ✓ |
| `transport_received` buffer exists, anti-teleports cross-node writes | `node.h:75`, flushed at `node.cpp:92–95` | ✓ |
| `update_canalization()` reads `last_auxin_flux` to update biases | `node.cpp:560–579` | ✓ |
| `initial_radius` genome param | `genome.h:83`, default `0.015 dm` | ✓ (contradicts doc examples — see §2.1) |
| Seed identified by `parent == nullptr` | `node.h:41` | ✓ |
| Tick order: vascular_transport → DFS walk | `plant.cpp` | ✓ |
| Biases cleaned up on `die()`, transferred on `replace_child()` | `node.cpp:52–60`, `node.cpp:312–313` | ✓ |
| `pipe_capacity()` uses pure `π × r² × conductance` (no bias term) | `vascular.cpp:36–38` | ✓ |

---

## Critical Bugs

### Bug 1: `last_auxin_flux.clear()` Wipes PIN-Recorded Flux

**Severity: BLOCKING**

`node.cpp:346` — `last_auxin_flux.clear()` is the first thing `transport_with_children()`
does. Tick order is:

```
1. vascular_transport()
2. pin_transport()              ← writes last_auxin_flux on parent nodes
3. DFS walk
   └── transport_chemicals()
         └── transport_with_children()  ← line 346: clears last_auxin_flux  ← BUG
               update_canalization()    ← reads empty map; sees diffusion only
```

`auxin_flow_bias` would be driven entirely by 5%-per-tick diffusion flux, not PIN — the
entire canalization mechanism silently becomes the old one.

**Fix:** Move the clear so that PIN and diffusion both accumulate into the same map before
`update_canalization` reads it:

```
pre-tick: last_auxin_flux already clear (cleared at end of prior tick)
pin_transport()          → writes last_auxin_flux
DFS walk:
  transport_with_children()  → accumulates diffusion flux (no clear here)
  update_canalization()      → reads combined flux, then clears last_auxin_flux
```

Remove `last_auxin_flux.clear()` from `transport_with_children()`. Add it to the end of
`update_canalization()` instead. Add this as an explicit step in the Implementation Plan.

---

### Bug 2: Seed Junction Distributes Stale Auxin

**Severity: BLOCKING**

The seed junction pseudocode (design doc lines 155–165) distributes from
`chemical(ChemicalID::Auxin)`:

```cpp
float to_send = std::min(chemical(ChemicalID::Auxin) * share, max_cap);
```

But during the PIN pass, shoot-side children post their auxin into
`seed->transport_received[Auxin]`. Per anti-teleportation rules, `transport_received` is
flushed into `chemical()` only at the end of `Node::tick()` — after the DFS walk,
AFTER the PIN pass. At the moment the seed distributes to root children, the
newly-collected shoot auxin is in `transport_received`, not yet in `chemical(Auxin)`.
The seed would distribute last tick's auxin level, not the current one.

The design says "After flushing its `transport_received` buffer from shoot-side children"
but provides no mechanism for that flush to occur during the PIN pass.

**Fix:** The PIN pass should accumulate incoming shoot auxin into a local variable,
not into `transport_received`:

```cpp
// Collect from shoot children into a local accumulator
float collected = 0.0f;
for (Node* child : shoot_children) {
    float moved = std::min(child->chemical(Auxin), max_cap * efficiency);
    child->chemical(Auxin) -= moved;
    seed->last_auxin_flux[child] += moved;
    collected += moved;
}
// Distribute to root children via transport_received (normal anti-teleport path)
for (Node* root_child : root_children) {
    float share = root_child->radius / total_root_radius;
    float to_send = std::min(collected * share, root_max_cap);
    root_child->transport_received[Auxin] += to_send;
    collected -= to_send;
}
seed->chemical(Auxin) += collected;  // remainder stays at seed
```

The seed acts as a synchronous relay within the PIN pass — no buffer needed for the
transit node itself. Root children get `transport_received` normally (anti-teleport).

---

## Parameter Calibration Error

### Wrong `initial_radius` in Calibration Examples

**Location:** "Effective Reach" section, lines ~377–395 in the design doc.

The doc says "newly-created thin stem (`r = initial_radius ≈ 0.05 dm`)" but
`genome.h:237` sets `initial_radius = 0.015 dm`. The actual numbers:

```
r  = 0.015 dm
r² = 0.000225 dm²
max_capacity = 0.000225 × 100 = 0.0225 AU/tick  (at full efficiency)
cold-start   = 0.0225 × 0.2  = 0.0045 AU/tick   (at pin_base_efficiency)
```

The doc's examples show `0.25` and `0.05` — off by 11×. The calibration target (r² ≈
0.005 dm²) corresponds to `r ≈ 0.07 dm`, which is 4.7× initial_radius — not a "thin
stem." The established trunk example (r = 0.25 dm, max_capacity = 6.25 AU/tick) is
correct.

**Impact:** `pin_capacity_per_area = 100` may need post-implementation tuning. The
cold-start saturation formula `saturation = efficiency = pin_base_efficiency = 0.2` is
self-consistent regardless of radius, so the feedback structure is sound. Only the
absolute throughput numbers in the prose are wrong.

**Required fix:** Correct the thin-stem example to use `r = 0.015 dm`. Note that
`pin_capacity_per_area = 100` was calibrated for a different radius and will need
empirical tuning.

---

## Consistency Issues

### `stress_boost` Missing from Thickening Formula

Design doc (line ~262):
```
delta_radius = cambium_responsiveness × auxin_flow_bias × sugar_gf
```

Current `stem_node.cpp:50–51`:
```cpp
float stress_boost = 1.0f + chemical(ChemicalID::Stress) * g.stress_thickening_boost;
float actual_rate = g.cambium_responsiveness * bias * sugar_gf * stress_boost;
```

`stress_thickening_boost` is not on any removal list and remains in the genome. The term
appears intentionally retained (thigmomorphogenesis — mechanical stress accelerates
cambium). Update the formula in the design doc to include `× stress_boost`.

### `transient_gain` and `transient_rate` Not Listed for Removal

The new `update_canalization` formula replaces the existing one completely:

**Old (`node.cpp:567–570`):**
```cpp
float target = flux * g.transient_gain;
flow_bias += (target - flow_bias) * g.transient_rate;
```

**New (design doc):**
```
current_saturation = auxin_flux / (radius² × pin_capacity_per_area)
auxin_flow_bias = lerp(previous, current_saturation, smoothing_rate)
```

`transient_gain` and `transient_rate` become unused but are not listed for removal in
Implementation Plan Step 2. The correct disposition:

- **Delete `transient_gain`** — the saturation denominator (`pin_capacity_per_area`)
  replaces the raw flux multiplier.
- **Rename `transient_rate` → `smoothing_rate`**, change default from `0.2` → `0.1`.
  The lerp role is identical; this is a rename with recalibration, not a new param.

Add to Implementation Plan Step 2: "Delete `transient_gain`. Rename `transient_rate` →
`smoothing_rate`, update default from 0.2 → 0.1."

### `vascular_radius_threshold` Has No Default Value

Mentioned in Implementation Plan Step 2 as a new genome param that replaces
`vascular_conductance_threshold`. Does not appear in the New Genome Parameters table.
No default value specified anywhere.

With `initial_radius = 0.015 dm`, set `vascular_radius_threshold = 0.01 dm` so all
newly-created internodes qualify from birth (initial_radius > threshold). Add to the
parameter table with this value and rationale.

### Plan.md Pseudocode Has a Null-Pointer UB

`plan.md` Step 4 pseudocode, line 235:
```cpp
auto it = parent ? parent->auxin_flow_bias.find(this) : parent->auxin_flow_bias.end();
```

If `parent` is null, `parent->auxin_flow_bias.end()` dereferences a null pointer —
undefined behavior. The correct form:
```cpp
float flow_bias = 0.0f;
if (parent) {
    auto it = parent->auxin_flow_bias.find(this);
    if (it != parent->auxin_flow_bias.end()) flow_bias = it->second;
}
```

---

## Completeness Gaps

### Test Suite Updates Completely Unmentioned

The Implementation Plan adds `tests/test_pin_transport.cpp` but says nothing about
the existing tests that will break. Grep results:

- `tests/test_node.cpp` — 8 uses of `structural_flow_bias` in assertions
- `tests/test_meristem.cpp` — 14 uses, including a named test "Thickening proportional
  to `structural_flow_bias`" and tests that pre-populate the map to bootstrap thickening
- `tests/test_vascularization.cpp` — 15 uses across all 4 integration tests, including
  "Zero structural_flow_bias → zero thickening" and "Canalization ratchet: auxin flux
  builds structural_flow_bias"

All 4 vascularization tests need rewriting around the new `auxin_flow_bias` / radius
model. Many meristem tests need their setup updated from `structural_flow_bias` injection
to auxin injection.

Add to Implementation Plan: "Update `test_node.cpp`, `test_meristem.cpp`, and
`test_vascularization.cpp` to remove `structural_flow_bias` setup and assertions;
rewrite vascularization tests around `auxin_flow_bias` and radius."

### `update_canalization()` Rewrite Not in Implementation Plan

The changed formula is described in the design body, but no implementation step names
`node.cpp`'s `update_canalization()` function as a target. An implementer following
only the steps could miss it. Add as an explicit sub-step to Step 2 or Step 4.

### `get_bias_multiplier()` Not Mentioned

`node.cpp:300–307`:
```cpp
return 1.0f + g.canalization_weight * (flow + structural);
```

After removing `structural_flow_bias`, the `structural` term disappears: becomes
`1.0f + g.canalization_weight * flow`. This function gates local diffusion sibling
weighting. It is not listed in any implementation step. Add to the cleanup list.

### Root Thickening Comment Contradicts New Design

`root_node.cpp:34–36`:
```cpp
// structural_flow_bias from sugar and cytokinin transport (not auxin —
// real root polar auxin transport governs patterning/gravitropism, not
// cambial signaling).
```

With the new design, root thickening reads `auxin_flow_bias`, which is built from PIN
auxin flux flowing acropetally through root internodes. This comment explicitly says
"not auxin" for root cambial signaling, which becomes incorrect. Update the comment.

---

## `structural_flow_bias` Audit — All Sites

The design says to delete `structural_flow_bias` entirely. Every source site and its
required action:

| File | Lines | Site | Action |
|------|-------|------|--------|
| `node.h` | 79, 84 | Declaration + `get_parent_structural_bias` comment | Delete field and method |
| `node.cpp` | 57–60 | `replace_child()` transfers structural bias | Delete block |
| `node.cpp` | 284–296 | `get_parent_structural_bias()` | Delete method |
| `node.cpp` | 304–305 | `get_bias_multiplier()` reads structural bias | Remove structural term |
| `node.cpp` | 313 | `die()` erases from parent structural bias | Delete line |
| `node.cpp` | 573–577 | `update_canalization()` ratchet | Replace with new lerp + saturation formula |
| `stem_node.cpp` | 22, 33, 38 | `thicken()` reads structural bias | Replace with `auxin_flow_bias` |
| `root_node.cpp` | 34, 38 | `thicken()` via `get_parent_structural_bias()` | Replace; update comment |
| `apical.cpp` | 181–190 | `spawn_internode()` stamps initial structural bias | Delete stamping block |
| `root_apical.cpp` | 103–109 | Same stamping | Delete stamping block |
| `vascular.cpp` | 22, 28–29 | `has_vasculature()` checks structural bias | Replace with `radius >= vascular_radius_threshold` |
| `vascular.cpp` | 254, 269–271 | Distribution weights use structural bias | Remove augmentation — weight by `r²` only |
| `vascular.cpp` | 323–324, 362 | Debug log reads + CSV column header | Remove or replace |
| `genome.h` | 56–57 | Comment references structural bias | Update |
| `genome.h` | 157–165 | `structural_threshold`, `structural_growth_rate`, `structural_max`, `vascular_conductance_threshold` | Delete all four |

The design doc explicitly addresses the `genome.h` removals and core `node.cpp` deletions.
Sites not mentioned in the Implementation Plan: `apical.cpp`, `root_apical.cpp`,
`get_bias_multiplier()`, the vascular debug log, and the test files (~37 references).

---

## Resolved Open Questions

**Root tip auxin maximum / PIN2**

The sim has `root_tip_auxin_production_rate = 0.03` (genome.h:226) — a local auxin
maximum at root tips via direct production, a reasonable engineering patch. PIN2's
columella recirculation loop is biologically accurate but architecturally complex.
**Decision: keep the current `root_tip_auxin_production_rate` as the root-tip auxin
source. PIN2 is out of scope for this milestone.** Note `root_tip_auxin_production_rate`
in the new genome parameters section so implementers know it is intentionally retained.

**Phototropic PIN redistribution**

Phototropism is currently handled geometrically. Lateral auxin redistribution via PIN3/PIN4
would require lateral PIN polarity modeling — a separate system.
**Decision: out of scope. Note as a future milestone direction.** The PIN infrastructure
would support this later without major restructuring.

**Wound response**

Stem removal and vascular regeneration are not modeled. Rethinking wound-stress interaction
is also out of scope.
**Decision: out of scope. Note as future direction.** The `last_auxin_flux` map in PIN
would naturally support wound vascular regeneration when eventually modeled.
