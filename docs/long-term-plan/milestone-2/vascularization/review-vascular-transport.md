# Vascular Transport Code Review

Reviewed against `docs/long-term-plan/milestone-2/vascularization/plan.md`.

Files examined: `src/engine/vascular.cpp`, `src/engine/vascular.h`, `src/engine/node/stem_node.cpp`, `src/engine/node/root_node.cpp`, `src/engine/node/node.h`, `src/engine/node/node.cpp`, `src/engine/genome.h`.

---

## What the Plan Claims

The plan describes five changes to bring the sim in line with the two-system transport architecture:

1. **Phase 2 bias weights** — junction distribution in the vascular pass should multiply each child's demand-proportional share by `(1 + canalization_weight * structural_flow_bias)`, using the same bias map that local diffusion already reads.
2. **Slow local diffusion** — reduce auxin, GA, and stress diffusion rates to make signals genuinely local.
3. **Refactor `thicken()`** — replace `thickening_rate × auxin_gf × sugar_gf` with `cambium_responsiveness × structural_flow_bias × sugar_available_fraction`. Remove `thickening_rate`, `auxin_thickening_threshold`, and the age-based cambium gate.
4. **Bias-based vascular maturation** — replace `has_vasculature()`'s age check with a `structural_flow_bias > vascular_conductance_threshold` test. Augment `pipe_capacity()` with a bias term: `π × r² × conductance × (1 + bias_conductance_scale × bias)`.
5. **Code annotation** — make the two-system architecture legible in comments.

---

## What the Code Does

### `has_vasculature()` — `vascular.cpp:18–25`

```cpp
bool has_vasculature(const Node& n, const Genome& g) {
    if (!n.parent) return true;  // seed always vascular
    if (n.type == NodeType::STEM)
        return n.age >= g.cambium_maturation_ticks;
    if (n.type == NodeType::ROOT)
        return n.age >= g.root_cambium_maturation_ticks;
    return false;
}
```

Pure age gate. No connection to `structural_flow_bias`. `cambium_maturation_ticks` (default 336 h = 14 days) and `root_cambium_maturation_ticks` (168 h = 7 days) are the only criteria.

### `pipe_capacity()` — `vascular.cpp:28–30`

```cpp
static float pipe_capacity(const Node& n, float conductance) {
    return 3.14159f * n.radius * n.radius * conductance;
}
```

Cross-sectional area times conductance. No bias augmentation.

### Phase 1 — `vascular.cpp:64–133`

Post-order walk. Classifies each node as source, sink, or conduit and propagates aggregated subtree supply and demand to parent. Conduits cap propagated flow at `pipe_capacity`. Supply and demand are kept separate (no within-subtree netting) so cross-branch phloem flow works — leaf sugar can reach root tips via the seed junction.

**Source/sink classification:**
- Sugar phloem: leaf nodes above `phloem_reserve_fraction × cap` are sources; `APICAL` and `ROOT_APICAL` are sinks (capacity deficit); any node with `starvation_ticks > 0` is an emergency sink (wants 50% of cap).
- Water/cytokinin xylem: root nodes (`ROOT`, `ROOT_APICAL`) are sources (50% of current level); `LEAF` and `APICAL` are water sinks; `APICAL` is a cytokinin sink with a **fixed demand of 0.05f** (not deficit-based, no capacity check).

### Phase 2 — `vascular.cpp:135–223`

Pre-order walk. Distributes `seed_available = min(total_supply, total_demand)` from seed outward.

**Junction distribution (lines 206–222):**
```cpp
float total_child_demand = 0.0f;
for (int ci : info.child_idxs) {
    total_child_demand += flat[ci].demand;
}
if (total_child_demand > 1e-8f) {
    for (int ci : info.child_idxs) {
        float share = available * (flat[ci].demand / total_child_demand);
        if (flat[ci].is_conduit) {
            share = std::min(share, flat[ci].capacity);
        }
        flat[ci].supply = share;
    }
}
```

Purely demand-proportional. No `structural_flow_bias` weighting. The `structural_flow_bias` map exists on every `Node` (public field, `unordered_map<Node*, float>`, `node.h:79`) and is accessible here as `flat[i].node->structural_flow_bias.find(flat[ci].node)` — it is not read.

### `thicken()` — `stem_node.cpp:21–42`

```cpp
void StemNode::thicken(const Genome& g, const WorldParams& world) {
    if (age < g.cambium_maturation_ticks) return;

    float density_scale = g.wood_density / world.reference_wood_density;
    float effective_rate = g.thickening_rate;

    float auxin_gf = std::min(chemical(ChemicalID::Auxin) /
                              std::max(g.auxin_thickening_threshold, 1e-6f), 1.0f);
    effective_rate *= auxin_gf;

    float max_cost = effective_rate * world.sugar_cost_stem_growth * density_scale;
    float sugar_gf = (max_cost > 1e-6f) ? std::min(chemical(ChemicalID::Sugar) / max_cost, 1.0f) : 1.0f;
    if (sugar_gf <= 1e-6f) return;

    float stress_boost = 1.0f + chemical(ChemicalID::Stress) * g.stress_thickening_boost;
    float actual_rate = effective_rate * sugar_gf * stress_boost;
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_stem_growth * density_scale;
    radius += actual_rate;
}
```

Driven by `thickening_rate × auxin_gf × sugar_gf × stress_boost`. `stress_boost` is applied multiplicatively (currently additive factor). `auxin_gf` is a linear ramp from 0 to 1 over the range `[0, auxin_thickening_threshold]` — it saturates at threshold, so anything above threshold thickens at full rate. The age gate (`cambium_maturation_ticks`) applies here too.

`RootNode::thicken()` (`root_node.cpp:30–`) mirrors this structure with the same genome parameters.

### Genome parameters — `genome.h`

Present and unchanged from current state:
- `thickening_rate = 0.00004f` (dm/hr)
- `auxin_thickening_threshold = 0.03f`
- `cambium_maturation_ticks = 336`
- `root_cambium_maturation_ticks = 168`
- `canalization_weight = 1.0f` (exists and used in local diffusion)
- `phloem_conductance = 8.0f`, `xylem_conductance = 10.0f`
- `phloem_reserve_fraction = 0.3f`

Missing (plan adds these):
- `cambium_responsiveness` (replaces `thickening_rate`)
- `vascular_conductance_threshold` (replaces age gate)
- `bias_conductance_scale` (augments `pipe_capacity`)

---

## Gaps and Corrections

### Gap 1 — Phase 2 ignores `structural_flow_bias` (plan Step 1, most important)

The distribution loop at `vascular.cpp:212–213` does a plain demand-proportional split. The bias weight `(1 + canalization_weight * structural_flow_bias[child])` described in the plan is not applied. `structural_flow_bias` is publicly accessible via `flat[i].node->structural_flow_bias`, so there is no structural barrier to adding it — it just has not been done.

### Gap 2 — `pipe_capacity()` does not include bias (plan Step 4)

`pipe_capacity` uses radius only. The plan's formula `π × r² × conductance × (1 + bias_conductance_scale × bias)` has not been added. This is a one-line change once `bias_conductance_scale` exists in the genome.

### Gap 3 — `has_vasculature()` uses age, not bias threshold (plan Step 4)

The age-based maturation gate is still in place. The plan wants `structural_flow_bias > vascular_conductance_threshold` (read from `n.parent->structural_flow_bias[&n]`). There is a small implementation wrinkle: the function signature takes `const Node& n`, while the bias map is keyed by `Node*` (non-const). Looking up `n.parent->structural_flow_bias.find(const_cast<Node*>(&n))` works but is inelegant. Cleaner to change the signature to `Node& n` (non-const) or use a helper on Node that accepts `const Node*` for lookup.

### Gap 4 — `thicken()` not refactored to read `structural_flow_bias` (plan Step 3)

`thicken()` still uses `thickening_rate`, `auxin_gf`, and the age gate. The plan's formula (`cambium_responsiveness × structural_flow_bias × sugar_available_fraction`) has not been implemented. Access pattern `parent->structural_flow_bias[this]` works directly since `parent` is a public `Node*` on every node.

One difference to flag: the plan's new formula drops the `stress_boost` multiplier that exists in the current code (`1.0f + chemical(Stress) * stress_thickening_boost`). The plan's formula is `cambium_responsiveness × bias × sugar_fraction` with no stress term. This is presumably intentional — stress drives cambium activity indirectly via auxin flux rather than directly — but worth confirming before removing it, since `stress_thickening_boost` has a non-zero default.

### Gap 5 — Cytokinin sink demand is hardcoded (minor correctness issue, not in plan)

In Phase 1 (`vascular.cpp:114`) and Phase 2 (`vascular.cpp:195`), `APICAL` nodes demand exactly `0.05f` of cytokinin regardless of how much they already have. There is no capacity or deficit check. This means a fully-saturated shoot tip makes the same demand as a depleted one. For water and sugar the code correctly uses `cap - current` deficit — cytokinin should follow the same pattern. Not mentioned in the plan, but worth fixing for consistency.

### Gap 6 — Phase 2 capacity cap does not redistribute uncapped surplus (minor, not in plan)

At `vascular.cpp:215–217`, when a child's share is capped by `pipe_capacity`, the excess is silently discarded rather than redistributed to siblings with remaining headroom. For most plants this is a small error, but in a heavily-branched plant with bottlenecked laterals it could cause meaningful flow loss. A second redistribution pass (similar to the logic in `transport_with_children`) would fix it.

---

## Implementation Feasibility

All four plan changes are straightforward given the current code structure.

**Step 1 (bias weights in Phase 2)** is the easiest and highest leverage. Replace the Phase 2 distribution loop body (approx. `vascular.cpp:207–219`) with:

```cpp
float total_biased_demand = 0.0f;
for (int ci : info.child_idxs) {
    auto it = info.node->structural_flow_bias.find(flat[ci].node);
    float bias = (it != info.node->structural_flow_bias.end()) ? it->second : 0.0f;
    flat[ci].capacity = /* existing capacity */;
    // store biased weight temporarily in supply field, overwritten below
    flat[ci].supply = flat[ci].demand * (1.0f + g.canalization_weight * bias);
    total_biased_demand += flat[ci].supply;
}
if (total_biased_demand > 1e-8f) {
    for (int ci : info.child_idxs) {
        float share = available * (flat[ci].supply / total_biased_demand);
        if (flat[ci].is_conduit)
            share = std::min(share, flat[ci].capacity);
        flat[ci].supply = share;
    }
}
```

This does not require any new genome parameters — `canalization_weight` already exists.

**Step 3 (refactor `thicken()`)** requires adding `cambium_responsiveness` to genome and reading `parent->structural_flow_bias[this]`. If the map has no entry for `this` (e.g., the node was just created), default to 0.0f — no bias, no thickening. That matches the desired behavior for brand-new internodes. The `density_scale` and sugar-cost accounting should be preserved since those physics are correct.

**Step 4 (bias maturation gate)** is a two-line change to `has_vasculature()` once `vascular_conductance_threshold` is in the genome. Recommend changing the function signature from `const Node&` to `Node&` to avoid `const_cast`, or adding a const-safe lookup method on `Node`.

**Step 4 (`pipe_capacity` augmentation)** is a one-line change once `bias_conductance_scale` is a genome parameter. Whether to use `0.0` (disabled) or a tuned default is the main uncertainty — start at 0.0 to keep the existing behavior and tune upward.

---

## Suggestions

1. **Implement Step 1 first, in isolation.** It requires no genome changes, no structural refactoring, and no test breakage risk. Verification is visual: `--color sugar` on a branched plant should show the main axis carrying more sugar than equally-sized laterals after a few hundred ticks of canalization accumulation.

2. **Fix the cytokinin demand calculation** (Gap 5) as part of Step 1. It's a two-line change and makes all three vascular chemicals use consistent logic.

3. **Add a redistribution pass** after Phase 2 capacity-capping (Gap 6). The pattern exists in `transport_with_children` — extract it into a shared helper if desired.

4. **Keep `stress_boost` in the new `thicken()` formula**, or explicitly decide to drop it. The plan's formula omits it without comment. If mechanical stress is expected to drive cambium activity beyond what auxin flux captures (e.g., a bent branch that thickens to resist load), the term belongs. If the plan intends stress-thickening to emerge purely from increased auxin flux in stressed regions, drop it and remove `stress_thickening_boost` from the genome.

5. **For `has_vasculature()` with bias threshold**: changing the signature to `Node&` is cleaner than `const_cast`. The function is only called from `vascular.cpp` where the nodes are mutable anyway.

6. **The `phloem_reserve_fraction` logic** (leaves keep 30% of sugar cap to themselves) is not mentioned in the plan but is a meaningful design choice. It prevents vascular bulk-flow from completely draining photosynthesizing leaves. Worth documenting explicitly as a design constant, not just a parameter.
