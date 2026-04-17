# PIN Transport Design — Cold Review

Reviewer reads the design from scratch with no prior context. All line numbers refer to
`pin-transport-design.md` unless otherwise noted.

---

## 1. Critical Bugs

These would produce silent wrong behavior if not caught before implementation.

### 1.1 `last_auxin_flux.clear()` Wipes PIN-Recorded Flux

**Severity: BLOCKING**

The design says (line 143):

> PIN records flux into `last_auxin_flux` on the parent, the same map that
> `update_canalization` already reads to update `auxin_flow_bias` each tick.

But in `node.cpp:346`:

```cpp
void Node::transport_with_children(const Genome& g) {
    if (children.empty()) return;
    last_auxin_flux.clear();   // ← wipes the map at the start of each call
    ...
```

Tick order is:
```
1. vascular_transport()
2. pin_transport()              ← PIN records into last_auxin_flux
3. DFS walk
   └── per node: transport_chemicals()
         └── transport_with_children()  ← clears last_auxin_flux  ← BUG
               update_canalization()    ← now reads an empty or diffusion-only map
```

When the DFS walk reaches each node, `transport_with_children` immediately clears
`last_auxin_flux`, erasing everything PIN recorded. `update_canalization` then sees only
local diffusion flux, not PIN flux. `auxin_flow_bias` will be driven by 5%-per-tick
diffusion, not by PIN — completely defeating the purpose of the new system.

**Fix options:**
- (A) Move `last_auxin_flux.clear()` to after `update_canalization()` completes, and have
  PIN *accumulate* into it during the PIN pass. Then `transport_with_children` accumulates
  local diffusion on top, and `update_canalization` sees the combined flux.
- (B) Have PIN write to a separate `pin_auxin_flux` map and update `update_canalization`
  to sum both maps. Cleaner separation but more state.
- Option A is simpler and consistent with the design's intent.

This fix must be part of the Implementation Plan.

---

### 1.2 Seed Junction: Collected Auxin Is in `transport_received`, Not `chemical()`

**Severity: BLOCKING**

The seed junction pseudocode (lines 155–165) distributes from `chemical(ChemicalID::Auxin)`:

```cpp
float to_send = std::min(chemical(ChemicalID::Auxin) * share, max_cap);
root_child->transport_received[ChemicalID::Auxin] += to_send;
chemical(ChemicalID::Auxin) -= to_send;
```

But during the PIN pass (before the DFS tick), shoot-side children have just posted auxin
into `seed->transport_received[Auxin]`. Per anti-teleportation rules, `transport_received`
is flushed into `chemical()` only at the end of `Node::tick()` — after the DFS walk, not
during the PIN pass. The PIN pass runs *before* the DFS walk.

So at the moment the seed junction code runs, the newly-collected shoot auxin is in
`transport_received`, NOT in `chemical(Auxin)`. The distribution to root children would
send the seed's **previous tick** auxin level, not the current tick's incoming signal.

The design says "After flushing its `transport_received` buffer from shoot-side children"
but provides no mechanism for this flush to actually occur during the PIN pass.

**Fix:** The PIN pass implementation needs to explicitly fold the seed's shoot-side
collected auxin into a local accumulator variable rather than routing through
`transport_received`. Concretely:

```cpp
// Collect from shoot children
float collected = 0.0f;
for (Node* child : shoot_children) {
    float moved = ...;
    child->chemical(Auxin) -= moved;
    collected += moved;
    seed->last_auxin_flux[child] += moved;
}
// Distribute directly (not via transport_received)
for (Node* root_child : root_children) {
    float share = root_child->radius / total_root_radius;
    float to_send = std::min(collected * share, max_cap);
    root_child->transport_received[Auxin] += to_send;
    collected -= to_send;
}
seed->chemical(Auxin) += collected;  // remainder stays at seed
```

The anti-teleportation rule applies to root children (they get `transport_received`), but
the seed itself acts as a synchronous relay within the PIN pass — no buffer needed for the
transit node.

---

## 2. Parameter Calibration Error

### 2.1 Wrong `initial_radius` in Calibration Examples

**Severity: HIGH — misleads implementer on tuning**

The design (lines 377–378) states:

> On a **newly-created thin stem** (`r = initial_radius ≈ 0.05 dm`, efficiency = 0.2):
> ```
> max_capacity = 0.05² × 100 = 0.25 AU/tick
> moves        = min(available, 0.25 × 0.2) = min(available, 0.05 AU/tick)
> ```

The actual `initial_radius` in `genome.h:237` is **`0.015 dm`**, not `0.05 dm`. With the
real value:

```
max_capacity = 0.015² × 100 = 0.000225 × 100 = 0.0225 AU/tick
moves        = min(available, 0.0225 × 0.2)  = min(available, 0.0045 AU/tick)
```

The doc overestimates thin-stem throughput by **11×**. The calibration target (lines
230–233) also uses r² = 0.005 dm² (r ≈ 0.07 dm, nearly 5× initial_radius) for what it
calls "a thin active stem."

**Actual behavior with `pin_capacity_per_area = 100`:**

| Stem | radius (dm) | r² (dm²) | max_capacity (AU/tick) | cold-start moves |
|------|------------|----------|----------------------|-----------------|
| New internode | 0.015 | 0.000225 | 0.0225 | 0.0045 |
| Mature lateral | 0.05 | 0.0025 | 0.25 | 0.05 (or 0.25 at full eff.) |
| Main trunk | 0.25 | 0.0625 | 6.25 | 6.25 (at full eff.) |

The doc's established trunk example (r = 0.25 dm) is correct. The self-limiting property
holds once the trunk is thick enough, but new internodes are far more bottlenecked than
described.

**Impact on the saturation formula:** At cold-start, a new stem moves
`r² × capacity × efficiency = r² × capacity × pin_base_efficiency` AU/tick. Saturation =
`flux / (r² × capacity) = (r² × capacity × efficiency) / (r² × capacity) = efficiency =
pin_base_efficiency = 0.2` by definition. So cold-start saturation is always exactly
`pin_base_efficiency` — this is self-consistent and correct, regardless of radius. The
calibration table is fine as a reference for mid-life stems; just remove the claim that
`initial_radius ≈ 0.05 dm`.

**Required fix:** Correct the "newly-created thin stem" example to use `r = 0.015 dm`.
Note that `pin_capacity_per_area = 100` may need post-implementation tuning given the
actual initial radius is 3× smaller than the doc's example.

---

## 3. Consistency Issues

### 3.1 `stress_boost` Missing from Thickening Formula

**Design doc (line 262):**
```
delta_radius = cambium_responsiveness × auxin_flow_bias × sugar_gf
```

**Current `stem_node.cpp:50–51`:**
```cpp
float stress_boost = 1.0f + chemical(ChemicalID::Stress) * g.stress_thickening_boost;
float actual_rate = g.cambium_responsiveness * bias * sugar_gf * stress_boost;
```

The existing code has a `stress_boost` term (thigmomorphogenesis — mechanical stress
accelerates cambium). The design doc drops it silently. This is either an intentional
removal (the doc should say so) or an omission.

Given `stress_thickening_boost` stays in the genome after this change (it's not on the
removal list), the intent is almost certainly to keep `stress_boost`. Update the formula
in the design doc to match.

### 3.2 `transient_gain` and `transient_rate` Not Listed for Removal

The new `update_canalization` formula (lines 192–197) completely replaces the existing
one:

**Old (node.cpp:567–569):**
```cpp
float target = flux * g.transient_gain;
flow_bias += (target - flow_bias) * g.transient_rate;
```

**New (design doc):**
```
current_saturation = auxin_flux / (radius² × pin_capacity_per_area)
auxin_flow_bias    = lerp(auxin_flow_bias, current_saturation, smoothing_rate)
```

The new formula makes `transient_gain` and `transient_rate` completely unused. But
Implementation Plan Step 2 only lists removing `structural_growth_rate`, `structural_threshold`,
`structural_max`, `vascular_conductance_threshold`. It does not say to remove `transient_gain`
and `transient_rate`, which remain in `genome.h:154–155` and `default_genome()`.

The new `smoothing_rate` serves the same role as `transient_rate` (lerp speed), and
`pin_capacity_per_area` serves the same role as `transient_gain` (flux → bias scaling
denominator). Remove `transient_gain` and `transient_rate` from the genome and add to
Step 2 of the Implementation Plan.

### 3.3 `vascular_radius_threshold` Has No Default Value

The design mentions `vascular_radius_threshold` in Implementation Plan Step 2 as a new
genome param that replaces `vascular_conductance_threshold`. It does not appear in the
"New Genome Parameters" table (which only lists the three PIN params), and no default value
is specified anywhere in the document.

The current `vascular_conductance_threshold = 0.005` with initial structural bias stamp
`0.01 + length * 0.1` means all new internodes qualify immediately. With the radius-based
gate, `vascular_radius_threshold` should be set just below `initial_radius = 0.015 dm` —
e.g., `0.01 dm` — so all newly-created internodes qualify from birth. Specify this in the
parameter table.

### 3.4 Plan.md Pseudocode Has a Null-Pointer UB

`plan.md` line 235 (Step 4 pseudocode):
```cpp
auto it = parent ? parent->auxin_flow_bias.find(this) : parent->auxin_flow_bias.end();
```

If `parent` is null, `parent->auxin_flow_bias.end()` dereferences a null pointer —
undefined behavior. The correct form is:

```cpp
float flow_bias = 0.0f;
if (parent) {
    auto it = parent->auxin_flow_bias.find(this);
    if (it != parent->auxin_flow_bias.end()) flow_bias = it->second;
}
```

---

## 4. Completeness Gaps

### 4.1 Test Suite Updates Completely Unmentioned

The Implementation Plan says to add `tests/test_pin_transport.cpp` (Step 4) but is silent
on the extensive test suite changes that removal of `structural_flow_bias` requires. Grep
results from the test directory:

- `tests/test_node.cpp` — 8 direct uses of `structural_flow_bias` in test assertions
- `tests/test_meristem.cpp` — 14 uses, including a named test "Thickening proportional to
  structural_flow_bias" and multiple tests that manually pre-populate `structural_flow_bias`
  to bootstrap thickening
- `tests/test_vascularization.cpp` — 15 uses across all 4 integration tests, including
  "Zero structural_flow_bias → zero thickening" and "Canalization ratchet: auxin flux builds
  structural_flow_bias"

All 4 tests in `test_vascularization.cpp` would need rewriting to use the new auxin_flow_bias
/ radius-based model. Many tests in `test_meristem.cpp` pre-populate `structural_flow_bias`
as setup to get thickening to occur — these would need to inject auxin and PIN flux instead.

Add a Step to the Implementation Plan: "Update test_node.cpp, test_meristem.cpp, and
test_vascularization.cpp to remove structural_flow_bias setup and assertions; rewrite
vascularization tests around auxin_flow_bias and radius."

### 4.2 Initial Bias Stamping in `apical.cpp` and `root_apical.cpp` Not Addressed

`apical.cpp:181–190`:
```cpp
// Stamp initial structural_flow_bias: the meristem pre-specifies procambium...
float initial_bias = 0.01f + internode_length * 0.1f;
float& entry = internode->parent->structural_flow_bias[internode];
entry = std::max(entry, initial_bias);
```
Same pattern in `root_apical.cpp:103–109`.

The design says `structural_flow_bias` is deleted entirely, but neither `apical.cpp` nor
`root_apical.cpp` appears anywhere in the Implementation Plan. With the radius-based
`has_vasculature()` gate, the initial bias stamp is no longer needed (new internodes
qualify via `initial_radius >= vascular_radius_threshold`). These two sites should simply
delete the stamping block.

Add to Step 2 or Step 4: "Delete structural bias stamping in `ApicalNode::spawn_internode`
and `RootApicalNode`'s equivalent."

### 4.3 `get_bias_multiplier()` Not Mentioned

`node.cpp:300–307`:
```cpp
float Node::get_bias_multiplier(Node* child, const Genome& g) const {
    float flow = 0.0f, structural = 0.0f;
    ...
    return 1.0f + g.canalization_weight * (flow + structural);
}
```

After removing `structural_flow_bias`, this should return
`1.0f + g.canalization_weight * flow`. This function is called during local diffusion
sibling weighting. The design's table (line 311) says the result for diffusion is
`radius_factor × (1 + canalization_weight × auxin_flow_bias)`, which is correct after
this change. But `get_bias_multiplier` is not listed anywhere in the Implementation Plan.

### 4.4 `update_canalization()` Code Change Not in Implementation Plan

The new `auxin_flow_bias` formula is described in the design body but does not appear as
an explicit step in the Implementation Plan. The plan says to update `genome.h` (Step 2)
and the thickening function (Step 4), but `update_canalization()` in `node.cpp:560–578`
needs a full rewrite of its bias update logic. Add this as an explicit sub-step.

### 4.5 Vascular Distribution Weights Use `structural_flow_bias`

`vascular.cpp:269–271`:
```cpp
auto it = info.node->structural_flow_bias.find(flat[ci].node);
float bias = (it != info.node->structural_flow_bias.end()) ? it->second : 0.0f;
weights[k] = cap * (1.0f + g.canalization_weight * bias);
```

The design table (line 309) documents this correctly ("r² (pipe_capacity)" after removal).
But the Implementation Plan Step 2 says only "Add `vascular_radius_threshold`" under
genome changes; it doesn't explicitly say "remove bias augmentation from vascular
distribution weights in `run_vascular()`." The plan.md Step 5 says "In `pipe_capacity()`,
remove any bias augmentation" — but `pipe_capacity()` itself has no bias augmentation; the
augmentation is in the weights computation *calling* `pipe_capacity()`. The cite is
slightly off. Make it explicit: the weights computation in `run_vascular()` must change.

### 4.6 Vascular Debug Log Column Header References `structural_flow_bias`

`vascular.cpp:362`:
```cpp
"tick,junction_node_id,child_node_id,child_type,"
"chemical,demand,conductance_weight,delivered,"
"structural_flow_bias\n";
```

After removal, this column no longer has meaning. Either drop the column or replace it with
something useful (e.g., `auxin_flow_bias`). Minor but would produce a misleading CSV.

### 4.7 Root Thickening Comment Contradicts New Design

`root_node.cpp:34–36`:
```cpp
// structural_flow_bias from sugar and cytokinin transport (not auxin —
// real root polar auxin transport governs patterning/gravitropism, not
// cambial signaling).
```

With the new design, root thickening is driven by `auxin_flow_bias`, which is built by
auxin flux — including the shoot-derived auxin that PIN delivers acropetally through root
internodes. This comment explicitly says "not auxin" for root cambial signaling, which will
be wrong. The comment needs updating to reflect that root cambium is now activated by
PIN-transported auxin flux, consistent with the shoot model.

---

## 5. `structural_flow_bias` Grep vs. Doc Claims

The design says to delete `structural_flow_bias` entirely. Every current site and its
disposition:

| File | Lines | Site | Action |
|------|-------|------|--------|
| `node.h` | 79, 84 | Declaration + `get_parent_structural_bias` comment | Delete member + method |
| `node.cpp` | 57–60 | `replace_child()` transfers structural bias | Delete block |
| `node.cpp` | 284–296 | `get_parent_structural_bias()` implementation | Delete method |
| `node.cpp` | 304–305 | `get_bias_multiplier()` reads structural bias | Remove structural term |
| `node.cpp` | 313 | `die()` erases from parent structural bias | Delete line |
| `node.cpp` | 573–577 | `update_canalization()` ratchet update | Replace with new lerp formula |
| `stem_node.cpp` | 22, 33, 38 | `thicken()` reads structural bias | Replace with `auxin_flow_bias` |
| `root_node.cpp` | 34, 38 | `thicken()` reads via `get_parent_structural_bias()` | Replace; update comment |
| `apical.cpp` | 181–190 | `spawn_internode()` stamps initial structural bias | Delete stamping block |
| `root_apical.cpp` | 103–109 | Same stamping | Delete stamping block |
| `vascular.cpp` | 22, 28–29 | `has_vasculature()` checks structural bias | Replace with radius check |
| `vascular.cpp` | 254, 269–271 | Distribution weights use structural bias | Remove augmentation |
| `vascular.cpp` | 323–324, 362 | Debug log reads and labels structural bias | Update log |
| `genome.h` | 56–57 | Comment references structural bias | Update comment |
| `genome.h` | 157–165 | `structural_threshold`, `structural_growth_rate`, `structural_max`, `vascular_conductance_threshold` | Delete all four |

The design doc explicitly accounts for all **genome.h** removals and the main **node.cpp**
deletions, but does not name `apical.cpp`, `root_apical.cpp`, `get_bias_multiplier()`, or
the vascular debug log. Test files are also unmentioned (see §4.1).

---

## 6. Resolved Open Questions

The three open questions at the bottom of the design document are resolved as follows:

**Root tip auxin maximum / PIN2**

The sim has `root_tip_auxin_production_rate = 0.03` (genome.h:226), which produces a local
auxin maximum at root tips via direct production — a reasonable engineering patch. PIN2's
recirculation loop in the root columella is biologically elegant but architecturally
complex. **Decision: keep the current `root_tip_auxin_production_rate` as the root-tip
auxin source. PIN2 is out of scope for this milestone.** Note `root_tip_auxin_production_rate`
in the new genome parameters section so implementers know it is intentionally retained.

**Phototropic PIN redistribution**

Phototropism is currently handled geometrically (`meristem_phototropism_rate`,
`leaf_phototropism_rate`). Lateral auxin redistribution through PIN3/PIN4 is biologically
accurate but would require lateral PIN polarity modeling, which is a separate system.
**Decision: out of scope. Note as a future milestone direction.** The PIN transport pass
infrastructure would support adding this later without major restructuring.

**Wound response**

Stem removal and vascular regeneration (Sachs wound response) are not modeled. The plant
currently has gravity-based stress (`chemical(ChemicalID::Stress)`) but no cutting/damage
events. Rethinking wound-stress interaction is also out of scope.
**Decision: out of scope. Note as future direction.** The `last_auxin_flux` infrastructure
in PIN would naturally support wound vascular regeneration when it is eventually modeled —
blocked connections would lose auxin flux, causing re-routing.

---

## 7. Minor Notes

- The description of PIN direction for LEAF nodes (line 105) says "node pumps auxin toward
  its parent" — correct for the petiole direction toward the stem. Unambiguous. ✓
- The `root_cytokinin_production_rate` behavior is correctly left unchanged. ✓
- Local diffusion parameters are explicitly stated as unchanged (line 406). ✓
- `pipe_capacity()` in `vascular.cpp` currently has no bias augmentation itself — the
  augmentation is in the calling code. The phrase "pure radius: `π × r² × conductance`"
  accurately describes the function; the calling code is what changes. ✓
- The three new genome params (`pin_capacity_per_area`, `pin_base_efficiency`,
  `smoothing_rate`) should be added to a named linkage group in `genome_bridge.cpp`. The
  `canalization` group is the natural home. Add this to Step 2.

---

## 8. Overall Assessment

**Not ready to implement as written.** Three issues would produce silent wrong behavior:

1. **PIN flux gets erased by `last_auxin_flux.clear()`** — `auxin_flow_bias` would not
   respond to PIN transport at all.
2. **Seed junction distributes stale auxin** — shoot-derived signal would be one tick
   delayed at the seed-to-root handoff, breaking the through-plant gradient.
3. **Parameter calibration uses wrong `initial_radius`** — implementer would tune
   `pin_capacity_per_area` against examples that are wrong by 11×.

Additionally, the Implementation Plan is incomplete for someone who hasn't already read
every file in the codebase:
- No step for updating `update_canalization()`
- No step for deleting the initial bias stamping in `apical.cpp`/`root_apical.cpp`
- No step for updating `get_bias_multiplier()`
- No step for updating the test suite (which has ~37 uses of `structural_flow_bias`)
- No default value for `vascular_radius_threshold`
- `transient_gain` / `transient_rate` not listed for removal

The **design intent is sound and the biology is well-motivated**. The architecture
(post-order collect → seed junction → pre-order distribute, using `transport_received` and
`last_auxin_flux`) is consistent with existing infrastructure. Fixing the three critical
issues above and filling in the six completeness gaps would make this ready to implement.

Estimated scope of corrections: one pass over the design body to fix the seed junction
description and the thin-stem example, and one pass over the Implementation Plan to
enumerate all changed sites explicitly.
