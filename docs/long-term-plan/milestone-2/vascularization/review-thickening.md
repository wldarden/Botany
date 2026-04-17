# Thickening Code Review

## What the plan claims

The plan argues that the current thickening model is structurally wrong in two ways:

1. **Auxin is a gate, not a driver.** Any stem above `auxin_thickening_threshold` thickens at the same fixed `thickening_rate`. A main trunk with heavy auxin flux and a barely-active lateral with just-threshold auxin thicken identically given equal sugar. The parameter that should vary — how much auxin has been flowing through this connection over time — has no effect once the gate is open.

2. **Age controls vascular admission.** `has_vasculature()` uses `cambium_maturation_ticks` (336 ticks, 14 days for stems) as the gate into the vascular bulk transport network. This makes vascular development purely time-dependent rather than flux-dependent.

**The proposed fix:**

- Replace `thickening_rate` + `auxin_thickening_threshold` gate with `cambium_responsiveness × structural_flow_bias × sugar_available_fraction`.
- Replace age-based `has_vasculature()` with a `structural_flow_bias > vascular_conductance_threshold` check.
- Augment `pipe_capacity()` to include `structural_flow_bias` in the conductance calculation.
- The plan assumes `structural_flow_bias` is accessible as `parent->structural_flow_bias[this]`.

---

## What the code actually does

### `StemNode::thicken()` — [src/engine/node/stem_node.cpp:21–43](src/engine/node/stem_node.cpp)

```cpp
if (age < g.cambium_maturation_ticks) return;          // age gate (default 336 ticks / 14 days)
float effective_rate = g.thickening_rate;               // base rate (default 0.00004 dm/hr)
float auxin_gf = min(auxin / auxin_thickening_threshold, 1.0f);  // saturating gate, NOT linear driver
effective_rate *= auxin_gf;
float sugar_gf = min(sugar / max_cost, 1.0f);
float stress_boost = 1.0f + stress * stress_thickening_boost;   // stress hormone adds extra
float actual_rate = effective_rate * sugar_gf * stress_boost;
radius += actual_rate;
```

**Key point:** `auxin_gf` is clamped to 1.0. A trunk with auxin=0.06 and a trunk with auxin=0.6 both produce `auxin_gf = 1.0` (both above `auxin_thickening_threshold = 0.03`). The plan's description of this as a "gate not a driver" is accurate. There is no proportionality to sustained flux — only a saturation threshold.

**Stress boost is present and unmentioned in the plan.** `stress_thickening_boost = 1.0f` (default). This can double the thickening rate under high mechanical stress. The plan's proposed formula `cambium_responsiveness × structural_flow_bias × sugar_fraction` omits it entirely.

### `RootNode::thicken()` — [src/engine/node/root_node.cpp:30–47](src/engine/node/root_node.cpp)

```cpp
if (age < g.root_cambium_maturation_ticks) return;   // age gate (default 168 ticks / 7 days)
float effective_rate = g.thickening_rate;             // same parameter as StemNode
float sugar_gf = min(sugar / max_cost, 1.0f);
radius += effective_rate * sugar_gf;
```

**Root thickening deliberately excludes auxin.** The comment reads: "Root secondary growth is sugar-gated only — not auxin-driven. Real root tips maintain their own auxin via PIN recycling." No `auxin_gf`, no `stress_boost`. It uses the same `thickening_rate` as stems.

**The plan does not acknowledge this deliberate difference.** The plan speaks of `StemNode::thicken() (and the equivalent in RootNode)` as if they work the same way. They don't.

### `has_vasculature()` — [src/engine/vascular.cpp:18–25](src/engine/vascular.cpp)

```cpp
if (!n.parent) return true;                                    // seed always vascular
if (n.type == NodeType::STEM) return n.age >= g.cambium_maturation_ticks;
if (n.type == NodeType::ROOT) return n.age >= g.root_cambium_maturation_ticks;
return false;
```

Separate age thresholds for STEM (336 ticks) and ROOT (168 ticks). Leaves, meristems, and axillaries are never vascular (`return false` in all other cases).

### `pipe_capacity()` — [src/engine/vascular.cpp:28–30](src/engine/vascular.cpp)

```cpp
return 3.14159f * n.radius * n.radius * conductance;
```

Pure radius-squared area. No `structural_flow_bias` involvement.

### Genome parameters for thickening — [src/engine/genome.h:56–69, 149, 207–219, 294](src/engine/genome.h)

| Parameter | Default | Notes |
|---|---|---|
| `thickening_rate` | `0.00004f` dm/hr | ~3.5 mm radius/year; shared by stems and roots |
| `auxin_thickening_threshold` | `0.03f` | Saturation point; stems only |
| `cambium_maturation_ticks` | `336` (14 days) | Stem vascular admission and thickening unlock |
| `root_cambium_maturation_ticks` | `168` (7 days) | Root equivalents |
| `stress_thickening_boost` | `1.0f` | Stem only; absent from root thickening |

All five are registered as evolvable parameters in `genome_bridge.cpp`. `thickening_rate` and `auxin_thickening_threshold` are in the `shoot_growth` linkage group ([genome_bridge.cpp:187–191](src/evolution/genome_bridge.cpp)). `stress_thickening_boost` is in the `stress` linkage group ([genome_bridge.cpp:238](src/evolution/genome_bridge.cpp)).

### `structural_flow_bias` location — [src/engine/node/node.h:79](src/engine/node/node.h)

```cpp
std::unordered_map<Node*, float> structural_flow_bias;  // persistent — slow, permanent
```

Stored on the parent node, keyed by child pointer. Transferred on `replace_child` ([node.cpp:57–60](src/engine/node/node.cpp)). Cleaned up on `die()` ([node.cpp:291](src/engine/node/node.cpp)). **The plan's access pattern `parent->structural_flow_bias[this]` is correct.**

### Radius feedback into transport

Yes, radius feeds back into transport via `pipe_capacity()`. Thickening → larger radius → larger `pipe_capacity` → more sugar and water can move through that connection per tick. The feedback loop exists, but driven by radius alone, not by `structural_flow_bias` independently.

---

## Gaps and corrections

### 1. Plan correctly identifies the gate-vs-driver problem

The plan's core diagnosis is accurate. `auxin_gf` is `min(auxin / threshold, 1.0f)` — it saturates at 1.0 and provides no distinction between barely-above-threshold and heavily-above-threshold. Two stems above 0.03 thicken identically regardless of auxin level.

### 2. Root thickening biology is not what the plan assumes

The plan treats roots symmetrically with stems — "and the equivalent in RootNode." But the current root code explicitly avoids the auxin gate because root cambium responds to different signals than shoot cambium. The plan proposes `structural_flow_bias` as the driver for both, but `structural_flow_bias` is built by auxin flux. If root auxin flux is low (by design — roots maintain their own auxin via PIN recycling and the vascular system doesn't push much shoot-auxin into root tissue), then root structural bias will be low and roots will barely thicken. This may be biologically correct or not, but it's an implicit change in root biology that the plan doesn't call out.

### 3. `stress_thickening_boost` is unmentioned

The stem thickening formula includes a stress hormone multiplier (`1 + stress * stress_thickening_boost`). The plan's new formula omits it. The plan should either explicitly drop this (and explain why mechanical stress no longer accelerates cambium activity) or add a `× stress_factor` term to the proposed formula.

### 4. No delay gate after removing age check

Currently, no thickening occurs during the first 14 days (336 ticks) for stems. The plan replaces this with a `structural_flow_bias > vascular_conductance_threshold` gate. Since bias starts at 0 and accumulates slowly (`structural_growth_rate = 0.005` per tick when flux exceeds `structural_threshold = 0.05`), the effective delay is variable: a main stem carrying heavy auxin from day 1 will cross the threshold faster than a lateral. This is the biologically desirable behavior, but it means early-plant behavior changes significantly — currently zero thickening for 2 weeks; new model has thickening begin whenever flux accumulates enough bias, potentially within 20–40 ticks on the main axis.

### 5. `cambium_maturation_ticks` used in two separate places

The age threshold governs both `thicken()` (the growth step) and `has_vasculature()` (vascular transport admission). The plan mentions removing it from `has_vasculature()` (step 4) and from `thicken()` (step 3) but treats them as one change. In implementation they need to be disconnected — you might want to test `thicken()` and `has_vasculature()` independently.

### 6. `pipe_capacity()` augmentation is optional in the plan but structurally important

The plan suggests augmenting `pipe_capacity()` with `structural_flow_bias` but frames it as "or start as a fixed constant." If not done, the self-reinforcing loop is incomplete: wider stem carries more auxin → more bias → more thickening → wider stem → more capacity. But if capacity scales only by radius² (not bias), then bias and radius both grow but only radius drives capacity. The loop still works but bias is doing indirect work rather than direct work.

---

## What would break

### Tests that fail immediately after changing `thicken()` to use `structural_flow_bias`

**`test_meristem.cpp:42` — "Secondary growth thickens interior nodes, not tips"**
Sets `auxin = 1.0f` and `age = cambium_maturation_ticks` for all nodes, then checks `seed->radius` increases. With the new model, the seed has `structural_flow_bias = {}` (empty map — no parent) and `structural_flow_bias[this] = 0` for all children. Zero bias = zero `delta_radius`. Test fails: seed doesn't thicken.

**`test_meristem.cpp:364` — "Thickening deducts sugar"**
Same issue: sets `auxin = 1.0f`, expects sugar to decrease and radius to increase after one tick. New model: zero bias → no thickening → sugar unchanged → test fails.

**`test_meristem.cpp:519` — "Thickening scales with sugar level"**
Sets `auxin = 1.0f`, compares low-sugar vs high-sugar thickening. New model: structural bias is zero in a fresh plant regardless of auxin, so both thicken equally (at zero). Test fails.

**`test_meristem.cpp:352` — "Thickening does not occur without sugar"**
Sets `sugar = 0`, checks radius doesn't change. New model: `delta_radius = cambium_responsiveness × 0 × 0 = 0` (bias is zero in a fresh plant anyway). This test would still pass, but for the wrong reason — it now tests nothing meaningful since the failure mode changed.

### Tests that fail after changing `has_vasculature()`

No tests in `tests/` explicitly check vascular admission (no `has_vasculature` calls in test files). But tests that depend on sugar transport working correctly through a plant could be indirectly affected if bias-based admission delays vascular connectivity in ways the tests don't anticipate.

### Downstream code that reads radius and would behave differently

- `pipe_capacity()` ([vascular.cpp:28](src/engine/vascular.cpp)) — radius increase from thickening still feeds capacity; no structural change needed here unless you also augment with bias.
- `maintenance_cost()` in `StemNode` and `RootNode` — volume-based; proportional to `radius²`. Thicker stems cost more to maintain. If bias-driven thickening grows radius differently (slower early, faster on dominant branches), maintenance costs change accordingly.
- Droop/break physics — `stress_max_weight` in `Node` involves radius. Structural changes to radius growth profile affect how quickly a plant can support its own weight.

### Evolution parameters that need updating

`genome_bridge.cpp` registers `thickening_rate` and `auxin_thickening_threshold` in the `shoot_growth` linkage group. Evolved genomes (including `best_genome_monday.txt` and `best_genome_night.txt` in `src/data/`) have tuned values for these. After the refactor, saved genome files will have stale/extra fields and will be missing `cambium_responsiveness` and `vascular_conductance_threshold`.

---

## Implementation feasibility

Moderate. The infrastructure is already there — `structural_flow_bias` is maintained, transferred on chain growth, and accessible from parent nodes. The actual code changes are small (~10–15 lines in `stem_node.cpp`, `root_node.cpp`, and `vascular.cpp`). The harder work is:

1. **Genome parameter surgery** — remove `thickening_rate`, `auxin_thickening_threshold`, `cambium_maturation_ticks`, `root_cambium_maturation_ticks`; add `cambium_responsiveness` and `vascular_conductance_threshold`. These are referenced in `genome.h`, `genome_bridge.cpp`, `app_realtime.cpp`, `app_evolve.cpp`, `app_dump.cpp`, and `app_sugar_test.cpp`. Six files minimum.

2. **Test rewrites** — the three thickening tests (lines 42, 364, 519 in `test_meristem.cpp`) need to be rewritten. The new tests need to actually build up `structural_flow_bias` by running a plant for enough ticks, then check thickening. That requires longer integration-style test scenarios rather than single-tick micro-tests.

3. **Tuning `cambium_responsiveness`** — the default value needs to produce similar thickening to the current `thickening_rate = 0.00004`. The relationship is: `actual_rate ≈ cambium_responsiveness × structural_flow_bias`. A mature main stem with well-developed structural bias (say ~0.5–1.0) hitting `cambium_responsiveness = 0.00008` would give similar rates to the current default. But bias takes time to build up — early in plant life, the new model thickens less aggressively than the old one, which might need compensation.

4. **Root biology decision** — decide explicitly whether root `structural_flow_bias` is expected to accumulate (roots do carry some auxin via xylem transport) and at what levels. If root structural bias stays near zero, roots won't thicken, which may or may not be the desired behavior.

**Risks:**
- Very young plants may exhibit no thickening at all (zero bias at tick 0), where previously thickening started at tick 336. This changes visual growth trajectory noticeably.
- If `structural_threshold = 0.05` is too high relative to typical root auxin levels, root bias never accumulates and root thickening goes to zero permanently.
- Removing the age-based vascular gate means very young internodes may be admitted to the vascular pass immediately if they have any bias — even one tick of auxin flow above threshold starts the bias clock. This is probably fine but changes the vascular pass behavior in ways that need verification.

---

## Suggestions

1. **Address the `stress_thickening_boost` question explicitly** before coding. Either add a `× (1 + stress * stress_sensitivity)` term to the new formula (keeping the biological behavior) or justify removing it. The current 1.0 default means mechanical stress doubles thickening rate on stressed stems — that's a significant effect to quietly drop.

2. **Keep a `structural_flow_bias` floor for testing** (e.g., for `cambium_responsiveness = 0` plants, explicitly verify zero thickening). The "monocot" behavior the plan describes is compelling and worth a test.

3. **Split step 3 and step 4.** Refactoring `thicken()` and refactoring `has_vasculature()` are independent changes. The thickening refactor can be tested with the age gate still in place (using age-gated stems but bias-driven rates). Once that's verified, separately replace the age gate.

4. **Write a small helper to prime `structural_flow_bias`** for tests. Since all thickening tests need non-zero bias, a utility like `prime_structural_bias(plant, N_ticks)` would avoid copy-paste setup across multiple test cases.

5. **Root thickening is worth a comment in the plan.** The current code has a principled reason for excluding auxin from root thickening. If the plan intends `structural_flow_bias` to be the driver for roots, it should explicitly address what level of root structural bias is expected to accumulate, or carve out a separate `root_cambium_responsiveness` parameter.

6. **Consider `vascular_conductance_threshold` carefully.** The plan says a node joins the vascular network when bias exceeds this threshold. But the `structural_growth_rate = 0.005` per tick means at full flux you reach a bias of 0.5 after 100 ticks (about 4 days). If `vascular_conductance_threshold = 0.1`, main-axis internodes get vascular connectivity in ~20 ticks. Is that reasonable? Currently it's 14 days (336 ticks). The gap is large and the biological rationale for the new timing should be confirmed before committing to a default.
