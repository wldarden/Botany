# Canalization & Maturation Gate Code Review

## What the plan claims

### Canalization data
The plan assumes `structural_flow_bias` is a per-connection map stored on the parent node, keyed by child pointer. It proposes reading it in `thicken()` as `parent->structural_flow_bias[this]` and in `has_vasculature()` as `parent->structural_flow_bias[node]`. It also assumes `canalization_weight` and `structural_threshold` already exist in the genome.

### Maturation gate
The plan claims the current gate is purely age-based (`has_vasculature()` checks `node.age >= cambium_maturation_ticks`) and proposes replacing it with a bias threshold: a node participates in vascular transport when its parent's `structural_flow_bias` entry for it exceeds a new genome parameter `vascular_conductance_threshold`. It also proposes removing `cambium_maturation_ticks` and `root_cambium_maturation_ticks` from the genome.

### Vascular Phase 2 weighting
The plan assumes Phase 2 currently distributes flow purely by demand, and proposes multiplying each child's share by `(1 + canalization_weight * structural_flow_bias)` — the same formula already used in local diffusion.

### Thickening formula
The plan proposes removing `thickening_rate` and `auxin_thickening_threshold` and replacing `StemNode::thicken()` and `RootNode::thicken()` with `delta_radius = cambium_responsiveness × structural_flow_bias × sugar_available_fraction`.

---

## Canalization: what the code actually does

### Data structure
All three canalization maps live directly on `Node` as `std::unordered_map<Node*, float>`:

```
// node.h:78-80
std::unordered_map<Node*, float> auxin_flow_bias;       // transient — fast, decays
std::unordered_map<Node*, float> structural_flow_bias;  // persistent — slow, permanent
std::unordered_map<Node*, float> last_auxin_flux;       // per-tick: auxin moved per child
```

These are **per-parent-node maps keyed by child pointer**. The plan's assumption is exactly correct: `parent->structural_flow_bias[child]` is the data access pattern.

### Update logic
`node.cpp:538–557` — `Node::update_canalization()` runs after every `transport_with_children()` call:

```cpp
// node.cpp:538
void Node::update_canalization(const Genome& g) {
    for (Node* child : children) {
        float flux = 0.0f;
        auto it = last_auxin_flux.find(child);
        if (it != last_auxin_flux.end()) flux = it->second;

        // Transient: exponential chase toward flux-derived target
        float target = flux * g.transient_gain;
        float& flow_bias = auxin_flow_bias[child];
        flow_bias += (target - flow_bias) * g.transient_rate;

        // Structural: slow ratchet, only grows above threshold, never decays
        float& struct_bias = structural_flow_bias[child];
        if (flux > g.structural_threshold) {
            struct_bias += g.structural_growth_rate;
        }
        struct_bias = std::min(struct_bias, g.structural_max);
    }
}
```

`last_auxin_flux` is cleared at the start of `transport_with_children()` (`node.cpp:324`) and populated in both Phase 1 (`node.cpp:424`) and Phase 2 (`node.cpp:521`). It measures actual auxin moved per child per tick — the raw signal driving canalization.

### Where biases are read
`node.cpp:278–285` — `get_bias_multiplier()`:

```cpp
float Node::get_bias_multiplier(Node* child, const Genome& g) const {
    float flow = 0.0f, structural = 0.0f;
    auto it_f = auxin_flow_bias.find(child);
    if (it_f != auxin_flow_bias.end()) flow = it_f->second;
    auto it_s = structural_flow_bias.find(child);
    if (it_s != structural_flow_bias.end()) structural = it_s->second;
    return 1.0f + g.canalization_weight * (flow + structural);
}
```

This is called at two points in `transport_with_children()`:
- `node.cpp:384` — Phase 1 bias_mult computed for each child
- `node.cpp:462` — Phase 2 `weight = desired * bias_mult` for receiver weighting

So canalization already biases **local diffusion** (both phases). The combined `flow + structural` formula is used — both transient and structural bias contribute together.

### Lifecycle
- Biases transfer on `replace_child()` (`node.cpp:52–64`) — chain growth preserves branch history
- Biases cleaned up on `die()` (`node.cpp:290–292`) — erased from parent's maps when a node is removed

---

## Maturation gate: what the code actually does

### The function
`vascular.cpp:18–25`:

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

Exact conditions:
- **Seed** (no parent): always true
- **STEM**: `age >= cambium_maturation_ticks` (default 336 hours = 14 days)
- **ROOT**: `age >= root_cambium_maturation_ticks` (default 168 hours = 7 days)
- **All other types** (LEAF, APICAL, ROOT_APICAL): always false

This is a hard binary gate. There is no gradient or partial vasculature — a node is either fully in the vascular network or fully excluded.

### Where it is checked
Two call sites:

1. **`vascular.cpp:71`** — during Phase 1 post-order walk: `info.is_conduit = has_vasculature(n, g)`. Determines whether a node passes subtree supply/demand up through its pipe (conduit) or contributes only its own supply/demand.

2. **`node.cpp:353–359`** — during `transport_with_children()` to skip local diffusion of vascular chemicals on mature-to-mature edges:
   ```cpp
   bool parent_vascular = has_vasculature(*this, g);
   // ...
   if (is_vascular_chemical(dp.id) && parent_vascular && has_vasculature(*child, g)) {
       continue;
   }
   ```
   This prevents double-transport: if both parent and child are in the vascular network, local diffusion skips Sugar/Water/Cytokinin (the vascular pass already handled them).

### What excluded nodes do
Non-conduit nodes in the vascular pass still have their supply/demand tallied (they contribute to the totals), but they don't pipe subtree flow. Their children's supply/demand doesn't propagate through them — only their own local amount does (`vascular.cpp:132–136`). For delivery, they receive sugar/water/cytokinin as sinks, but the flow reaches them via the local diffusion pass (`transport_with_children`) rather than the vascular bulk pass.

---

## Gaps and corrections

### 1. Vascular Phase 2 does not use structural_flow_bias — confirmed gap
The plan is correct: `run_vascular()` Phase 2 at `vascular.cpp:210–219` distributes flow purely by demand:

```cpp
float share = available * (flat[ci].demand / total_child_demand);
```

No bias multiplier applied. `VascNodeInfo` doesn't carry any canalization data. This is the gap Step 1 of the plan addresses.

**Access pattern for the fix:** During `build_flat()` or Phase 2, each `VascNodeInfo` needs `parent->structural_flow_bias[this_node]`. The flat array has `parent_idx`, so at Phase 2 distribution time: `flat[flat[ci].parent_idx].node->structural_flow_bias[flat[ci].node]` gives the bias. Alternatively, add a `float struct_bias` field to `VascNodeInfo` and populate it during `build_flat()`. The latter is cleaner.

### 2. The plan's thicken() access pattern is exactly correct
`parent->structural_flow_bias[this]` works as written. The map is `unordered_map<Node*, float>`, parent is public, `this` is `Node*`. In `StemNode::thicken()`, `this` is the child node, `parent` is the parent. No issues.

**One subtlety**: If the node has no bias entry yet (new node, first tick), `operator[]` on `unordered_map` inserts a default 0.0f. That's fine for thickening (zero bias = zero delta) but means the map grows by one entry on first access. Use `.find()` or `.count()` if zero-initialization of the map entry is undesirable.

### 3. has_vasculature() const-correctness issue with bias lookup
`has_vasculature(const Node& n, const Genome& g)` takes a const ref. The proposed replacement accesses `n.parent->structural_flow_bias.find((Node*)&n)`. This requires a `const_cast` or a non-const lookup because:
- `structural_flow_bias` is `unordered_map<Node*, float>` with non-const key type
- `&n` in a const context is `const Node*`, not `Node*`
- `find(const Node* key)` won't compile without a cast

The fix is straightforward: cast `&n` to `Node*` for the map lookup, since we're only reading, not modifying. Or change the function signature to `bool has_vasculature(Node& n, const Genome& g)`. This is a minor implementation detail but must be handled.

### 4. Root thickening comment vs. plan
`RootNode::thicken()` (`root_node.cpp:30–47`) already does NOT use auxin gating (comment: "Root secondary growth is sugar-gated only"). The plan says to replace it with `structural_flow_bias`-driven thickening. This is consistent — the plan replaces the sugar-only root approach with a bias-driven one, same as the shoot.

However, note the plan's rationale for root thickening is less biologically obvious than for shoot: real root cambium activity isn't directly driven by polar auxin transport the way shoot cambium is. The plan acknowledges this only by saying "same formula applies." Worth flagging as a design choice, not a bug.

### 5. Parameters that don't exist yet
The plan references two new genome parameters that are **not yet in `genome.h`**:
- `cambium_responsiveness` — proposed replacement for `thickening_rate`
- `vascular_conductance_threshold` — minimum structural bias for vascular participation

Both `thickening_rate` (`genome.h:57`) and `auxin_thickening_threshold` (`genome.h:58`) are still present. `cambium_maturation_ticks` (`genome.h:62`) and `root_cambium_maturation_ticks` (`genome.h:69`) are also still present.

### 6. plan mentions `pipe_capacity()` augmentation — not yet addressed
The plan proposes: `capacity = π × r² × conductance × (1 + bias_conductance_scale × structural_flow_bias)`. The current `pipe_capacity()` at `vascular.cpp:28–30` uses only radius:

```cpp
static float pipe_capacity(const Node& n, float conductance) {
    return 3.14159f * n.radius * n.radius * conductance;
}
```

This augmentation would require passing structural_flow_bias into `pipe_capacity()` — either as a parameter or by restructuring `VascNodeInfo`. If `VascNodeInfo` carries `float struct_bias`, it's a clean single-field addition.

---

## Implementation feasibility

### Step 1: Bias vascular Phase 2
**Easy.** Add `float struct_bias = 0.0f` to `VascNodeInfo`. Populate it in `build_flat()`: if the node has a parent, look up `flat[parent_idx].node->structural_flow_bias.find(node)`. In Phase 2, when computing each child's share: `share *= (1 + g.canalization_weight * flat[ci].struct_bias)`. Renormalize across siblings so total distributed doesn't change. No circular dependencies.

### Step 3: Replace thicken() with bias-driven formula
**Moderate.** Steps:
1. Add `cambium_responsiveness` to `Genome`, remove `thickening_rate` + `auxin_thickening_threshold`
2. In `StemNode::thicken()`: remove age gate and auxin gate, compute `bias = parent ? parent->structural_flow_bias.count(this) ? parent->structural_flow_bias.at(this) : 0.0f : 0.0f`, then `delta = g.cambium_responsiveness * bias * sugar_fraction`
3. In `RootNode::thicken()`: same, but root already has no auxin gate — just remove age gate and swap formula
4. Update `genome_bridge.cpp` for evolution system

The access pattern is straightforward. No circular dependency — `thicken()` is called from `update_tissue()` inside `tick()`, which runs before `transport_chemicals()`, which is where `update_canalization()` updates the bias. This means thickening always runs on the **previous tick's** bias — one tick lag, which is biologically reasonable (cambium responds to accumulated history, not instantaneous flux).

### Step 4: Replace age gate with bias gate
**Easy in code, needs care on semantics.** Change `has_vasculature()` from:
```cpp
return n.age >= g.cambium_maturation_ticks;
```
to:
```cpp
if (!n.parent) return false; // non-seed check already handled above
auto it = n.parent->structural_flow_bias.find(const_cast<Node*>(&n));
return it != n.parent->structural_flow_bias.end()
    && it->second >= g.vascular_conductance_threshold;
```

**Potential issue:** New nodes have no entry in `structural_flow_bias` at all for their first tick (the parent hasn't run `update_canalization()` for them yet). The `find()` correctly returns `end()` → false → no vasculature. This is correct behavior.

**Bootstrap concern:** With the age gate, every node eventually matures regardless of auxin history. With the bias gate, a node on a path that never carries enough auxin to exceed `structural_threshold` will *never* develop vasculature. This is the intended behavior per the plan, but it means `vascular_conductance_threshold` must be set carefully — if it's too high relative to `structural_growth_rate * ticks`, even the main trunk won't develop vasculature in reasonable simulation time.

**No circular dependencies.** `has_vasculature()` is called from:
1. `vascular.cpp` — runs before `tick_tree()`, so bias values are from previous tick. Fine.
2. `node.cpp` `transport_with_children()` — runs during DFS walk, bias values from current tick's `update_canalization()` haven't run yet for this node. Fine.

---

## Suggestions

### 1. Add struct_bias to VascNodeInfo now, even before Step 1
Even if Step 1 (bias weighting in Phase 2) isn't implemented yet, adding `float struct_bias` to `VascNodeInfo` and populating it in `build_flat()` is a trivial addition that makes Step 1, the pipe_capacity augmentation, and future debugging much easier. Zero cost to add it, significant convenience benefit.

### 2. Don't combine transient and structural in the vascular pass
The current local diffusion uses `get_bias_multiplier()` which sums `auxin_flow_bias + structural_flow_bias`. The plan correctly says only `structural_flow_bias` should gate vascular participation — transient bias is PIN protein redistribution, not physical xylem/phloem. The vascular phase should use `structural_flow_bias` only. This differs from local diffusion, which uses both.

### 3. Consider keeping cambium_maturation_ticks as a floor
The plan says to remove the age gate entirely and replace with the bias gate. An alternative: keep a very short age floor (e.g., 24–48 ticks = 1–2 days) as a minimum before any node can develop vasculature, combined with the bias gate. This prevents newly created nodes from being classified as vascular on their first tick even if they inherit bias through `replace_child()`. The current `replace_child()` transfers bias to the new node — a freshly inserted internode above a meristem could immediately have nonzero bias if the meristem had accumulated some. Whether this is the desired behavior is a design question.

### 4. Bootstrap for the very first tick
The seed starts with `structural_flow_bias` empty. On tick 1, the shoot apical produces auxin. `update_canalization()` runs and records flux for the seed→apical connection. `structural_threshold = 0.15` means we need flux > 0.15 per tick before structural bias accumulates. At the default `structural_growth_rate = 0.005`, it takes ~200 ticks above threshold to reach 1.0 bias. With the proposed `vascular_conductance_threshold` (value TBD), the main stem won't have vasculature for the first N ticks. This is biologically correct (primary phloem/xylem from the procambium comes first, secondary vasculature takes time), but the threshold value matters a lot for gameplay/simulation timing.

### 5. The `last_auxin_flux` measurement undercounts root paths
`last_auxin_flux` tracks auxin moved in `transport_with_children()` — local diffusion only. Since the vascular pass doesn't currently route auxin (auxin is explicitly excluded from `is_vascular_chemical()`), all auxin movement is through local diffusion. This is correct. But it means `structural_flow_bias` on root paths is driven by whatever auxin makes it down the shoot-to-root local diffusion chain. Given auxin's basipetal bias and the diffusion rate, the seed→root axis should accumulate some bias, but likely less than the seed→shoot axis. This asymmetry may be intended or may need adjustment.

### 6. File path for chemical_registry
The plan (Step 2) refers to `chemical_registry.h` for diffusion rate tuning. The actual file is `src/engine/chemical/chemical_registry.h:24`, not a standalone config. The diffusion rates come from the genome, so tuning them means changing `default_genome()` in `genome.h:168`. That's where the change actually lands.
