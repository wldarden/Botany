# Vascular Transport Code Review

## What the plan claims

The plan assumes the following about the existing vascular transport implementation:

1. **Phase 2 junction splitting** is demand-proportional but doesn't yet use `structural_flow_bias`. Step 1 of the implementation plan is to multiply each child's demand share by `(1 + canalization_weight * structural_flow_bias[child])` — the same formula local diffusion already uses via `get_bias_multiplier()`.

2. **`structural_flow_bias` is accessible** at junctions via `node->parent->structural_flow_bias[node]` — the existing canalization map on each parent node, keyed by child pointer.

3. **`has_vasculature()` uses an age-based gate** (`cambium_maturation_ticks`, `root_cambium_maturation_ticks`) that the plan wants to replace with a bias threshold (`vascular_conductance_threshold`).

4. **`pipe_capacity()` uses only radius** (`π × r² × conductance`), and the plan wants to augment it with bias: `capacity = π × r² × conductance × (1 + bias_conductance_scale × structural_flow_bias)`.

5. **`StemNode::thicken()` is an independent rate-driven step** that the plan wants to replace with `cambium_responsiveness × structural_flow_bias × sugar_available_fraction`.

6. The plan assumes the `canalization_weight` genome parameter already exists and that no new genome machinery is needed for Step 1.

---

## What the code actually does

### `has_vasculature()` — `vascular.cpp:18–25`

Age-based gate, exactly as the plan describes:

```cpp
bool has_vasculature(const Node& n, const Genome& g) {
    if (!n.parent) return true;  // seed is always a vascular junction
    if (n.type == NodeType::STEM)
        return n.age >= g.cambium_maturation_ticks;
    if (n.type == NodeType::ROOT)
        return n.age >= g.root_cambium_maturation_ticks;
    return false;
}
```

LEAF, APICAL, and ROOT_APICAL nodes are never conduits regardless of age. Non-STEM/ROOT node types always return false.

### `pipe_capacity()` — `vascular.cpp:28–30`

Pure area formula, no bias:

```cpp
static float pipe_capacity(const Node& n, float conductance) {
    return 3.14159f * n.radius * n.radius * conductance;
}
```

### Phase 1: Post-order aggregation — `vascular.cpp:64–139`

Walks the flat pre-order array backwards (giving post-order). For each node:

- Classifies source/sink based on `chem_id` and `NodeType`
- Conduit nodes propagate subtree totals capped at pipe capacity; non-conduits only contribute their own local supply/demand
- Supply and demand are aggregated up to the seed

**Sugar (phloem) source/sink rules (lines 77–92):**
- Source: `LEAF` with sugar above `phloem_reserve_fraction × sugar_cap`
- Sink: `APICAL` and `ROOT_APICAL` (capacity deficit)
- Emergency sink: any node with `starvation_ticks > 0` (50% of capacity)

**Xylem (water + cytokinin) source/sink rules (lines 96–116):**
- Source: `ROOT` and `ROOT_APICAL` — 50% of their current chemical level
- Water sink: `LEAF` (capacity deficit) and `APICAL` (capacity deficit)
- Cytokinin sink: `LEAF` (fixed 0.01f) and `APICAL` (fixed 0.05f)

Note: STEM nodes that are conduits pass subtree flow but are neither sources nor sinks in any pass. The seed node (STEM type, no parent) accumulates its subtree's net supply and demand at `flat[0]`.

### Phase 2: Pre-order distribution — `vascular.cpp:141–222`

Seed (index 0) starts with `min(supply, demand)` as the flow budget. Walks outward. At each node:

1. Deducts from local source if applicable
2. Delivers to local sink (up to `min(available, local_demand)`)
3. Distributes remainder to children **purely proportionally to Phase 1 demand** — lines 205–221:

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
        flat[ci].supply = share;  // reuse supply field to pass flow down
    }
}
```

**No bias weighting anywhere in this loop.** Siblings of equal demand receive equal flow regardless of how much auxin has historically flowed through one vs. the other.

### Seed node behavior

- Always a conduit (`has_vasculature()` returns true when `!n.parent`)
- In Phase 1, accumulates the tree's net supply/demand
- In Phase 2, starts with `min(supply, demand)` — distributes to children proportionally
- Seed can itself be a source or sink (it's STEM type, but the source/sink code only checks LEAF/ROOT/APICAL — seed is none of those, so it's neither; it's purely a junction/conduit)

### `StemNode::thicken()` — `stem_node.cpp:21–43`

Current model has three gates:

1. Age gate: `age < g.cambium_maturation_ticks` → return early (same cutoff as vascular maturation)
2. Auxin gate: `auxin_gf = min(auxin / auxin_thickening_threshold, 1.0f)` — scales rate, doesn't block entirely
3. Sugar gate: exits if sugar can't fund any growth

Then: `radius += thickening_rate * auxin_gf * sugar_gf * stress_boost`

Stress hormone additionally boosts the rate via `1.0f + stress * g.stress_thickening_boost`.

`structural_flow_bias` is not read here at all.

### How `structural_flow_bias` is stored and accessed

- Declared on `Node` as `std::unordered_map<Node*, float> structural_flow_bias` (node.h:79)
- Keyed by child pointer; stored on the parent
- Updated by `Node::update_canalization()` (node.cpp:538–555): increments by `structural_growth_rate` per tick when `last_auxin_flux[child] > structural_threshold`, capped at `structural_max`
- Used in local diffusion via `get_bias_multiplier(child, g)` (node.cpp:278–285): returns `1.0f + canalization_weight * (flow_bias + structural_bias)`
- Preserved across `replace_child()` (node.cpp:57–60): the bias entry is rekeyed to the new child pointer when chain growth reparents
- Cleaned up via `die()` (node.cpp:291): entry erased from parent's map

In Phase 2 of the vascular pass, when distributing to children at node index `i`, the parent's map is directly accessible as `flat[i].node->structural_flow_bias`. Looking up a specific child: `flat[i].node->structural_flow_bias.find(flat[ci].node)`. **This access requires no struct changes to `VascNodeInfo` or any changes to how the flat array is built.**

### NodeType discrepancy

The actual `NodeType` enum in `node.h:19–22` contains **5 types**: `STEM, ROOT, LEAF, APICAL, ROOT_APICAL`.

CLAUDE.md documents **7 types** (`STEM, ROOT, LEAF, SHOOT_APICAL, SHOOT_AXILLARY, ROOT_APICAL, ROOT_AXILLARY`). The vascular code uses `APICAL` (not `SHOOT_APICAL`), and there are no axillary node types in the codebase. The plan and CLAUDE.md appear to describe an aspirational or historical architecture that is partially implemented. The vascular code is consistent with the 5-type enum that actually exists.

---

## Gaps and corrections

### G1: Phase 2 does not weight by `structural_flow_bias`

**The most important gap.** The plan's central premise — that vascular junctions use canalization history to privilege well-developed branches — is not implemented. Line 212: `float share = available * (flat[ci].demand / total_child_demand)` uses raw demand with no bias multiplier. A main trunk and a lateral branch of equal demand receive exactly equal vascular flow, regardless of how many hundreds of ticks the main trunk has been carrying auxin.

### G2: `has_vasculature()` uses age, not bias

`cambium_maturation_ticks` (336 hours / 14 days) and `root_cambium_maturation_ticks` (168 hours / 7 days) gate vascular participation. A young internode on the main axis with heavy auxin flux is excluded; an old internode on a dormant lateral that activated recently is admitted when its age clock expires. Plan's Step 4 is not implemented.

### G3: `pipe_capacity()` does not include bias

Only `π × r² × conductance`. The plan's `bias_conductance_scale` augmentation does not exist, and `bias_conductance_scale` is not in `genome.h` at all.

### G4: `thicken()` does not read `structural_flow_bias`

Still uses `thickening_rate × auxin_gf × sugar_gf × stress_boost`. The `cambium_responsiveness` parameter proposed in the plan does not exist in `genome.h`. `thickening_rate` and `auxin_thickening_threshold` still exist and are used. Plan's Step 3 is not implemented.

### G5: `const`-correctness wrinkle in `has_vasculature()`

`has_vasculature()` signature is `bool has_vasculature(const Node& n, const Genome& g)`. To do the bias-based check, it would need `n.parent->structural_flow_bias.find(const Node*)` — but the map key type is `Node*` (non-const). `&n` yields `const Node*`, which doesn't match. Fixes: change the map key to `const Node*`, or change the function signature to take `Node*` (non-const), or use `const_cast<Node*>(&n)` for the lookup. The plan doesn't mention this.

### G6: Phase 1 demand aggregation doesn't account for bias

The plan proposes multiplying demand by the bias weight in Phase 2, but Phase 1 aggregates unweighted demand upward. This means the seed's total demand budget is computed without knowing how bias will redistribute it. In practice this is fine — redistribution among siblings doesn't change the total flow budget — but it means the supply/demand balance computed in Phase 1 is slightly "wrong" relative to what Phase 2 will actually deliver. Low-demand children with high bias will receive more than their Phase 1 demand suggested; high-demand children with low bias receive less. The underallocated high-demand children may appear starved from the vascular pass's perspective even though total supply is adequate.

A cleaner approach: also weight Phase 1 demand aggregation. Then Phase 1 and Phase 2 use consistent demand figures. The plan doesn't discuss this.

### G7: Cytokinin sourcing is identical to water sourcing

For cytokinin, `surplus = n.chemical(chem_id) * 0.5f` for ROOT and ROOT_APICAL nodes. Cytokinin is produced by root apicals during their `tick()`, which runs after the vascular pass. So the cytokinin available to the vascular pass this tick is last tick's production. This is correct (one-tick lag is expected and consistent with the anti-teleportation model), but the plan doesn't discuss it. It's worth noting in comments.

---

## Implementation feasibility

### Step 1: Add `structural_flow_bias` weighting to Phase 2

**Easy.** No struct changes, no new genome parameters (both `canalization_weight` and `structural_flow_bias` maps already exist). The change is isolated to the child-distribution loop in `run_vascular()`.

Exact change at `vascular.cpp:205–221`:

```cpp
// Current:
float total_child_demand = 0.0f;
for (int ci : info.child_idxs) {
    total_child_demand += flat[ci].demand;
}

// Proposed:
float total_weighted_demand = 0.0f;
for (int ci : info.child_idxs) {
    float bias = 0.0f;
    auto it = info.node->structural_flow_bias.find(flat[ci].node);
    if (it != info.node->structural_flow_bias.end()) bias = it->second;
    total_weighted_demand += flat[ci].demand * (1.0f + g.canalization_weight * bias);
}

// In distribution loop, replace:
// float share = available * (flat[ci].demand / total_child_demand);
// with:
float bias = 0.0f;
auto it = info.node->structural_flow_bias.find(flat[ci].node);
if (it != info.node->structural_flow_bias.end()) bias = it->second;
float weighted_demand = flat[ci].demand * (1.0f + g.canalization_weight * bias);
float share = (total_weighted_demand > 1e-8f)
    ? available * (weighted_demand / total_weighted_demand)
    : 0.0f;
```

The `Genome& g` parameter is already present in `run_vascular()`. No other changes required.

**One hidden risk**: if a child has zero Phase 1 demand but nonzero structural bias, it still gets zero share (because `0.0 * (1 + bias) = 0`). This is correct behavior — bias amplifies demand, it doesn't create demand from nothing.

### Step 4: Bias-based maturation gate

**Moderate.** Requires fixing the `const`-correctness issue (see G5). Best approach: change `has_vasculature()` to accept `Node*` (non-const). This affects all call sites (inside `run_vascular()`), which already have non-const access via `info.node`.

New signature: `bool has_vasculature(Node* n, const Genome& g)`. Inside:
```cpp
bool has_vasculature(Node* n, const Genome& g) {
    if (!n->parent) return true;  // seed
    if (n->type != NodeType::STEM && n->type != NodeType::ROOT) return false;
    auto it = n->parent->structural_flow_bias.find(n);
    float bias = (it != n->parent->structural_flow_bias.end()) ? it->second : 0.0f;
    return bias >= g.vascular_conductance_threshold;
}
```

Add `vascular_conductance_threshold` to `genome.h` and `default_genome()`. Remove `cambium_maturation_ticks` and `root_cambium_maturation_ticks` (also used in `stem_node.cpp:22` and presumably `root_node.cpp` — audit both before removing).

**Risk**: For the first N ticks of a seedling's life, `structural_flow_bias` on the seed→first-stem connection may be below threshold, meaning the first stem is not a conduit and the whole vascular system does nothing. Ensure the seed node is exempt (it already is via `!n->parent`), and consider a low default threshold so young plants still have functional vasculature.

### Step 3: Refactor `thicken()`

**Straightforward** once `structural_flow_bias` access is established. Add `cambium_responsiveness` to genome, remove `thickening_rate` and `auxin_thickening_threshold` (also audit `genome_bridge.cpp` in the evolution subsystem — these are likely named genes there and would need to be swapped).

---

## Suggestions

**1. Do Step 1 first, exactly as scoped.** It's three lines of change, uses no new genome parameters, and requires zero interaction with the thickening or maturation refactors. It's verifiable immediately by running `--color sugar` and observing unequal distribution at junctions after canalization builds up. The other steps (3 and 4) are bigger and can be done independently.

**2. Consider weighted Phase 1 aggregation.** The plan doesn't mention it, but propagating bias-weighted demand up in Phase 1 would make supply/demand accounting consistent end-to-end. Without it, a subtree with many high-bias connections will receive more vascular flow than its Phase 1 demand predicted — which is fine for behavior but makes the accounting harder to reason about.

**3. The age-to-bias transition for `has_vasculature()` needs a seedling safety check.** A brand-new plant has no structural bias anywhere. If `vascular_conductance_threshold > 0`, the entire vascular system is inactive until auxin has flowed for enough ticks to accumulate bias. During this window, leaves can't deliver sugar via phloem and roots can't deliver water via xylem. This may or may not be acceptable depending on whether the seed's initial sugar reserves can fund the gap. Consider either a low default threshold (0.01), an explicit age fallback, or ensuring `seed_sugar` is large enough for the first ~100 ticks without vascular sugar delivery.

**4. `structural_threshold = 0.15f` in the default genome is high.** This means only connections carrying auxin flux above 0.15 units/tick accumulate structural bias. Given `apical_auxin_baseline = 0.15f` and decay, many connections may never cross this threshold in a small plant. Step 1 will have no visible effect until canalization has actually built up. Verify this isn't zero in practice by logging `structural_flow_bias` values on seed's children after 200+ ticks before assuming Step 1 is working.

**5. Genome bridge audit before removing parameters.** `src/evolution/genome_bridge.cpp` translates between `botany::Genome` and `evolve::StructuredGenome`. Before removing `thickening_rate`, `auxin_thickening_threshold`, `cambium_maturation_ticks`, and `root_cambium_maturation_ticks`, check whether they appear as named genes in the bridge. If the evolution system is actively running, removing them without updating the bridge will break serialization of saved genomes.

**6. Document the one-tick lag for cytokinin vascular transport.** Root apicals produce cytokinin during `tick()`, which runs after `vascular_transport()`. The cytokinin available to the vascular pass each tick is always from the previous tick's production. This is correct and consistent with the anti-teleportation model, but it's easy to forget and debug. A comment in `vascular.cpp` near the cytokinin pass noting this would prevent future confusion.
