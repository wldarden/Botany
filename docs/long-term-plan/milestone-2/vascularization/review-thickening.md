# Thickening Code Review

Review of `StemNode::thicken()` and `RootNode::thicken()` against `plan.md`.

---

## What the Plan Claims

The plan argues that secondary thickening is vascular development — wood is old xylem, and the cambium divides in response to auxin flowing *through* it. The current model inverts this causality: it uses auxin as a gate (open/closed) rather than a driver (rate proportional to flux). The proposed fix makes `structural_flow_bias` — the canalization ratchet already maintained in `node.cpp` — the direct driver of thickening:

```
delta_radius = cambium_responsiveness × structural_flow_bias × sugar_available_fraction
```

This eliminates three genome parameters (`thickening_rate`, `auxin_thickening_threshold`, `cambium_maturation_ticks` / `root_cambium_maturation_ticks`) and replaces them with one (`cambium_responsiveness`). It also replaces the age-based `has_vasculature()` gate in `vascular.cpp` with a bias threshold, making vascular admission path-dependent rather than time-dependent.

The plan presents this as a five-step refactor, with Steps 1–2 (vascular Phase 2 bias weighting, local diffusion rate reductions) feeding better canalization data into Steps 3–4 (thicken refactor, vascular maturation gate). Step 5 is documentation.

---

## What the Code Does

### `StemNode::thicken()` — `src/engine/node/stem_node.cpp:21–43`

```cpp
void StemNode::thicken(const Genome& g, const WorldParams& world) {
    if (age < g.cambium_maturation_ticks) return;            // L24: age gate

    float density_scale = g.wood_density / world.reference_wood_density;
    float effective_rate = g.thickening_rate;                // L27: fixed rate

    float auxin_gf = std::min(                               // L31: auxin as gate, not driver
        chemical(ChemicalID::Auxin) / std::max(g.auxin_thickening_threshold, 1e-6f), 1.0f);
    effective_rate *= auxin_gf;                              // L32: clamps at 1.0 — caps not scales

    float max_cost = effective_rate * world.sugar_cost_stem_growth * density_scale;
    float sugar_gf = (max_cost > 1e-6f) ?
        std::min(chemical(ChemicalID::Sugar) / max_cost, 1.0f) : 1.0f;
    if (sugar_gf <= 1e-6f) return;                          // L37: sugar gate

    float stress_boost = 1.0f + chemical(ChemicalID::Stress) * g.stress_thickening_boost; // L39
    float actual_rate = effective_rate * sugar_gf * stress_boost;
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_stem_growth * density_scale;
    radius += actual_rate;                                   // L42
}
```

**Inputs:** `age`, `g.cambium_maturation_ticks`, `g.thickening_rate`, `g.auxin_thickening_threshold`, node's auxin level, node's sugar level, `g.wood_density`, `world.reference_wood_density`, `world.sugar_cost_stem_growth`, node's stress hormone, `g.stress_thickening_boost`.

**`structural_flow_bias` is not read here.** The function is self-contained and never touches the canalization maps.

**Called from** `StemNode::update_tissue()` at `stem_node.cpp:17`, which is called from `Node::tick()` at `node.cpp:79` (inside `update_tissue(plant, world)`), after the maintenance/starvation check at `node.cpp:75` but before physics and transport.

### `RootNode::thicken()` — `src/engine/node/root_node.cpp:30–47`

```cpp
void RootNode::thicken(const Genome& g, const WorldParams& world) {
    if (age < g.root_cambium_maturation_ticks) return;       // L33: age gate

    // Root secondary growth is sugar-gated only — not auxin-driven.
    float effective_rate = g.thickening_rate;                // L38: same fixed rate as stem

    float max_cost = effective_rate * world.sugar_cost_stem_growth;
    float sugar_gf = (max_cost > 1e-6f) ?
        std::min(chemical(ChemicalID::Sugar) / max_cost, 1.0f) : 1.0f;
    if (sugar_gf <= 1e-6f) return;                          // L42: sugar gate

    float actual_rate = effective_rate * sugar_gf;
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_stem_growth;
    radius += actual_rate;                                   // L46
}
```

Root thickening is simpler: no auxin involvement, no stress boost, no wood density scale. Just age gate + sugar gate + fixed rate.

### `has_vasculature()` — `src/engine/vascular.cpp:18–25`

```cpp
bool has_vasculature(const Node& n, const Genome& g) {
    if (!n.parent) return true;
    if (n.type == NodeType::STEM)  return n.age >= g.cambium_maturation_ticks;
    if (n.type == NodeType::ROOT)  return n.age >= g.root_cambium_maturation_ticks;
    return false;
}
```

Age-based, as described. `structural_flow_bias` is not referenced in `vascular.cpp` at all.

### Where `structural_flow_bias` is actually used

`structural_flow_bias` is maintained per-parent (keyed by child pointer) in `node.cpp`:

- **Accumulated** at `node.cpp:551–553`: flux above `g.structural_threshold` increments `structural_flow_bias[child]` by `g.structural_growth_rate` per tick.
- **Used for transport weighting** at `node.cpp:282–284`: `get_bias_multiplier()` returns `1 + canalization_weight * (transient + structural)` and is applied in `transport_with_children()` to redistribute local diffusion among siblings.
- **Transferred on reparenting** at `node.cpp:57–60`: chain growth preserves bias history.
- **Erased on death** at `node.cpp:291`.

It is used only for local diffusion weighting. It does not feed into thickening or vascular Phase 2.

---

## Gaps and Corrections

### Gap 1: Plan Steps 1–4 are entirely unimplemented

None of the five implementation steps have been started. The `structural_flow_bias` → thickening connection does not exist. The vascular Phase 2 does not weight by bias. The maturation gate is still age-based. The genome still has `thickening_rate` and `auxin_thickening_threshold`. This document reviews the *current* state; all of plan.md is future work.

### Gap 2: Auxin is a saturating gate, not a proportional driver — but the plan overstates how wrong this is

`auxin_gf = min(auxin / threshold, 1.0)` — this is a ramp from 0 to 1 that saturates at the threshold. Below threshold, thickening is proportionally reduced; above threshold, thickening runs at full `thickening_rate` regardless of how much more auxin is present. The plan calls this "a gate, not a driver," which is accurate once auxin exceeds `auxin_thickening_threshold = 0.03`. In practice, the main trunk reliably exceeds this threshold and the lateral branches may or may not — so the distinction does matter for hierarchy, but it is a partial driver in the sub-threshold range.

### Gap 3: Root thickening uses no auxin at all — consistent with biology, inconsistent with plan's formula

The plan says the new formula applies to both `StemNode::thicken()` and `RootNode::thicken()`. But root thickening already has no auxin involvement — the comment at `root_node.cpp:35–37` explains why (root polar auxin transport drives gravitropism/patterning, not cambial signaling). This is biologically defensible. The new formula `cambium_responsiveness × structural_flow_bias × sugar_available_fraction` would still apply to roots correctly because `structural_flow_bias` accumulates on root connections from sugar/cytokinin/water transport flux — not solely from auxin. But the plan's biological rationale (cambium activated by auxin flux through it) applies less cleanly to roots. This is a minor inconsistency worth flagging when implementing Step 3.

### Gap 4: `stress_thickening_boost` has no place in the new model

`StemNode::thicken()` includes `stress_boost = 1.0 + stress * g.stress_thickening_boost` (line 39). This models thigmomorphogenesis — mechanical stress induces cambial activity. The plan does not mention this term and its proposed formula has no stress component. Under Step 3, this either needs to be preserved as an additive term or explicitly dropped. If dropped, plants under heavy mechanical load will no longer thicken faster at loaded nodes, which is a real biological phenomenon. Worth deciding intentionally.

### Gap 5: `wood_density` scaling is absent from the plan's formula

Current stem thickening applies `density_scale = g.wood_density / world.reference_wood_density` to both the sugar cost and the thickening rate (lines 26, 35, 41). Dense-wood plants thicken more slowly per unit sugar. The plan's formula `cambium_responsiveness × structural_flow_bias × sugar_available_fraction` uses `sugar_available_fraction` (sugar/capacity) rather than raw sugar, which implicitly normalizes for capacity. But it does not include density. The sugar cost side should still scale with density to maintain the energy accounting.

### Gap 6: `sugar_available_fraction` is a different formulation than current sugar gating

Current code: `sugar_gf = min(sugar / max_cost, 1.0)` — this is a fraction-of-affordability gate. Plan's formula: `sugar_available_fraction = sugar / sugar_cap` — this is a fullness fraction. These behave differently. Plan's version means a node at 50% sugar capacity thickens at half the bias-driven rate even if it has plenty of sugar for the actual cost. The current version means a node thickens at full rate as long as it can afford the full cost. Plan's version ties thickening rate to relative storage fullness, which is more biological (cambium division is sensitive to osmotic/turgor state, not just whether sugar is available). This distinction is intentional in the plan but deserves a comment.

---

## What Would Break if We Change It

### Tests that assume current behavior

**`test_meristem.cpp:42–64` — "Secondary growth thickens interior nodes, not tips"**

Sets `n.age = g.cambium_maturation_ticks` and `n.chemical(ChemicalID::Auxin) = 1.0f`, then checks `seed->radius > seed_r_before`. This test hardcodes the age gate and auxin gate. Under the new model: (a) the age gate disappears — the check `n.age = g.cambium_maturation_ticks` becomes meaningless, and (b) auxin 1.0 would not directly drive thickening — `structural_flow_bias` would need to be pre-populated for anything to happen. The test would fail silently (seed doesn't thicken because no bias has accumulated yet).

**`test_meristem.cpp:525–558` — Sugar reserve fraction test**

Seeds two plants, gives one very low sugar (`0.00002f`), confirms partial thickening. This tests the sugar-gated rate scaling. Under the new model the sugar term changes from `min(sugar / max_cost, 1.0)` to `sugar / sugar_cap` — still a sugar fraction, so this test concept survives, but the threshold values change. The exact numeric assertions would break.

**`test_evolution.cpp:46–47` — Genome serialization**

Round-trips `stress_thickening_boost` through the evolution bridge. If we keep `stress_thickening_boost`, this test is unaffected. If we drop it, it needs updating.

### Code that references removed genome params

`genome_bridge.cpp` (evolution serialization) maps `thickening_rate`, `auxin_thickening_threshold`, `cambium_maturation_ticks`, and `root_cambium_maturation_ticks` as named genes with mutation configs. These would all need to be removed and `cambium_responsiveness` (and `vascular_conductance_threshold`) added. Evolution runs against `best_genome.txt` files — existing saved genomes would have stale fields; the bridge should degrade gracefully (fields not found = use genome default).

### Vascular maturation gate

If `has_vasculature()` switches from age to bias threshold, every node that was previously admitted to vascular transport based on age alone is now conditionally excluded until its bias exceeds `vascular_conductance_threshold`. Young plants in particular would be fully local-diffusion-only until the main axis builds enough bias — this might make the first few days of plant development significantly slower for sugar transport if the threshold is too high. Needs careful calibration.

---

## Implementation Feasibility

**Steps are well-scoped.** Each step has a single entry point: Step 1 touches only `vascular.cpp`'s Phase 2 distribution loop. Step 3 touches only `StemNode::thicken()`, `RootNode::thicken()`, and `genome.h`. Step 4 touches only `has_vasculature()` in `vascular.cpp` and `pipe_capacity()`.

**The main risk is calibration, not architecture.** `structural_flow_bias` accumulates slowly (`structural_growth_rate = 0.005` per tick, saturates at `structural_max = 2.0` after ~400 ticks of sustained flux above `structural_threshold = 0.15`). If `cambium_responsiveness` is too low, plants will barely thicken for hundreds of hours. If too high, they'll thicken faster than currently. The plan suggests 0.001–0.01 per tick — needs bracketing against the current `thickening_rate = 0.00004f` to match observable trunk width.

**Step 1 has no risk** — adding bias weighting to vascular distribution is additive. If all biases are zero (new plant), `(1 + canalization_weight * 0) = 1`, behavior is identical to current.

**Step 2 (diffusion rate reduction) is the highest-risk step** — changing auxin diffusion rates affects branching patterns across the entire plant. Should be done with `botany_realtime` running and a before/after visual comparison. Do not combine with Step 3 in one commit.

**Steps 3 and 4 should go in one commit** — they are the same decision (move from age/rate model to bias/responsiveness model), and splitting them creates an inconsistent intermediate state where thickening reads bias but the vascular gate still uses age.

**Accessing `structural_flow_bias` from child to parent** requires `parent->structural_flow_bias[this]` inside `thicken()`. The `parent` pointer is available on `Node` but accessing a map keyed by `this` pointer from the child is slightly unusual. The plan explicitly calls this out (`node.cpp:133`). A safer access pattern: `(parent && parent->structural_flow_bias.count(this)) ? parent->structural_flow_bias.at(this) : 0.0f`.

---

## Suggestions

**1. Add a pre-Step 3 regression test now.** Before touching `thicken()`, write a test that: (a) grows a multi-branch plant for 500+ ticks, (b) reads `structural_flow_bias` values on the seed node for its shoot and root children, (c) asserts they are non-zero, and (d) asserts the main axis has higher bias than a lateral branch. This test will pass against the current code (bias already accumulates), and it becomes the Step 6 integration test from the plan. Having it in place before the refactor means any mis-wiring of the new thickening logic will show up immediately.

**2. Preserve `stress_thickening_boost` as an additive term.** The plan formula is `cambium_responsiveness × bias × sugar_fraction`. Add stress: `cambium_responsiveness × bias × sugar_fraction × (1 + stress_hormone * stress_thickening_boost)`. This preserves thigmomorphogenesis without complicating the core formula. Drop it only if empirical observation shows it doesn't produce useful behavior.

**3. Root thickening needs a different bias source.** Roots don't carry much auxin (the gradient is basipetal, so root nodes have low auxin). Their `structural_flow_bias` will be driven by sugar transport (which does run through the vascular pass) and cytokinin (acropetal, so roots have higher cytokinin). These still accumulate bias on well-used connections, so the formula still works — but the comment in `root_node.cpp` saying "not auxin-driven" should be updated to explain that roots use flow bias from sugar/cytokinin rather than auxin.

**4. Don't implement Step 2 (diffusion rates) in the same commit as Step 3.** Step 2 changes how much `structural_flow_bias` accumulates (slower diffusion = more localized auxin = more selective bias accumulation). If done simultaneously with the thickening refactor, it's impossible to tell which change produced any visual difference. Two separate commits with separate `botany_realtime` verification sessions.

**5. Consider a `cambium_bias_floor` parameter.** Under the new model, a stem with zero bias will never thicken, even if it's old and has been carrying auxin at below-`structural_threshold` levels for years. A small bias floor (`structural_threshold = 0.15` currently, so any flux below 0.15 produces zero bias) means that a lightly-used lateral branch stays perpetually thin. This is biologically reasonable for weak branches, but a `cambium_bias_floor` genome parameter (default 0) would allow monocot-adjacent plants to control exactly how thin their lateral stems can become.
