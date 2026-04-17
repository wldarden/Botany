# Canalization & Maturation Gate — Code Review

Review of how the existing canalization system works relative to the vascularization plan's claims, with gap analysis and implementation notes.

---

## What the Plan Claims

1. **Canalization storage** — `structural_flow_bias` is per-connection, stored on the parent node keyed by child pointer. `auxin_flow_bias` is the fast transient layer; `structural_flow_bias` is the slow permanent layer.

2. **Maturation gate replacement** — `has_vasculature()` currently uses `node.age >= cambium_maturation_ticks`. The plan proposes replacing it with a bias threshold: a node is vascular when `parent->structural_flow_bias[this] >= vascular_conductance_threshold`. Seed (no parent) stays always-vascular.

3. **Thickening formula refactor** — Replace the current age-gate + `auxin_gf` formula in `StemNode::thicken()` and `RootNode::thicken()` with:
   ```
   delta_radius = cambium_responsiveness × structural_flow_bias × sugar_available_fraction
   ```
   Remove `thickening_rate`, `auxin_thickening_threshold`, `cambium_maturation_ticks`, `root_cambium_maturation_ticks`.

4. **Vascular pipe capacity** — Augment `pipe_capacity()` in `vascular.cpp` to scale conductance by `structural_flow_bias`, so well-canalized connections carry more bulk flow.

---

## Canalization: What the Code Does

### Storage (per-connection, on parent)

`node.h:78-80`:
```cpp
std::unordered_map<Node*, float> auxin_flow_bias;       // transient — fast, decays
std::unordered_map<Node*, float> structural_flow_bias;  // persistent — slow, permanent
std::unordered_map<Node*, float> last_auxin_flux;       // per-tick auxin moved per child
```
Both bias maps and the flux record are on the **parent** node, keyed by child pointer. **This matches the plan exactly.**

### Update logic

`node.cpp:538-556` — `Node::update_canalization(const Genome& g)` called at the end of `transport_chemicals()` after `transport_with_children()`:

```cpp
// node.cpp:538
void Node::update_canalization(const Genome& g) {
    for (Node* child : children) {
        float flux = 0.0f;
        auto it = last_auxin_flux.find(child);
        if (it != last_auxin_flux.end()) flux = it->second;

        // Transient: exponential chase toward flux * transient_gain
        float target = flux * g.transient_gain;
        float& flow_bias = auxin_flow_bias[child];
        flow_bias += (target - flow_bias) * g.transient_rate;  // line 548

        // Structural: ratchet, never decays
        float& struct_bias = structural_flow_bias[child];
        if (flux > g.structural_threshold) {
            struct_bias += g.structural_growth_rate;           // line 553
        }
        struct_bias = std::min(struct_bias, g.structural_max); // line 555
    }
}
```

`last_auxin_flux` is cleared at the start of each `transport_with_children()` call (`node.cpp:324`), then accumulated during both inflow and outflow phases (`node.cpp:425`, `522`).

### How bias is read

`node.cpp:278-284` — `Node::get_bias_multiplier(Node* child, const Genome& g)`:
```cpp
return 1.0f + g.canalization_weight * (flow_bias + structural_bias);
```

This multiplier is applied in Phase 2 of `transport_with_children()` — the parent-giving-to-children pass — at `node.cpp:462` and `511`. It biases the **proportional share** a child receives from its parent. It does **not** amplify the total budget; it only redistributes what the parent was already giving.

### Lifecycle

- `replace_child()` (`node.cpp:51-60`): transfers both bias entries from old child to new child pointer — chain growth preserves branch history.
- `die()` (`node.cpp:289-291`): erases `this` from parent's both bias maps on node removal.
- New children start with no entry (treated as 0.0 by `find()` guard).

---

## Maturation Gate: What the Code Does

`vascular.cpp:18-25`:
```cpp
bool has_vasculature(const Node& n, const Genome& g) {
    if (!n.parent) return true;          // seed always vascular
    if (n.type == NodeType::STEM)
        return n.age >= g.cambium_maturation_ticks;
    if (n.type == NodeType::ROOT)
        return n.age >= g.root_cambium_maturation_ticks;
    return false;
}
```

Declared in `vascular.h:26`; also called inside `node.cpp:353-359` during `transport_with_children()` to skip local diffusion for vascular chemicals on mature-to-mature edges.

`genome.h:61,69` — defaults:
```
cambium_maturation_ticks      = 336   // 14 days until shoot thickening
root_cambium_maturation_ticks = 168   // 7 days until root thickening
```

Used in two places:
1. `vascular.cpp:71` — `info.is_conduit` flag for the bulk transport pass
2. `stem_node.cpp:24` and `root_node.cpp:33` — guard at the top of `thicken()`

### Current thickening formula (`stem_node.cpp:21-43`):
```cpp
if (age < g.cambium_maturation_ticks) return;
float auxin_gf = std::min(auxin / auxin_thickening_threshold, 1.0f);
float effective_rate = g.thickening_rate * auxin_gf;
float sugar_gf = std::min(sugar / max_cost, 1.0f);
radius += effective_rate * sugar_gf * stress_boost;
```

---

## Gaps and Corrections

### Gap 1: Thickening has zero coupling to `structural_flow_bias` (critical)

The plan describes `structural_flow_bias` as the future driver of thickening, but in the current code `StemNode::thicken()` uses local auxin concentration and an age gate. `structural_flow_bias` is used **only** in `transport_with_children()` to redistribute local diffusion shares. There is no connection between the canalization system and secondary growth at all. This is the core gap the plan addresses.

### Gap 2: `vascular.cpp::pipe_capacity()` ignores structural bias

`vascular.cpp:28-30`:
```cpp
static float pipe_capacity(const Node& n, float conductance) {
    return 3.14159f * n.radius * n.radius * conductance;
}
```

The plan calls for augmenting this with `structural_flow_bias`, but the vascular pass currently has no access to bias data. `VascNodeInfo` doesn't carry it; the flat-array walk would need the parent's `structural_flow_bias[child]` for each node.

### Gap 3: Thickening reads `auxin_flow_bias` indirectly via concentration, not via the stored bias

Current `auxin_gf` uses local `chemical(ChemicalID::Auxin)` at the node — which is what diffused down. The plan replaces this with the accumulated `structural_flow_bias`, which is a cumulative history of how much auxin *flowed through* the connection, not how much is present. These are different signals; the structural bias is more suitable as it's noise-resistant and permanent.

### Gap 4: `thicken()` is called on the child node but the bias is on the parent

`StemNode::thicken()` runs on the node itself. To read `structural_flow_bias` it needs `parent->structural_flow_bias[this]`. There is no current helper for this. Adding `float structural_bias_from_parent() const` to `Node` would be the cleanest approach.

### Correction 1: Plan's storage claim is accurate

The plan says biases are "per-connection, stored on parent, keyed by child pointer." This is exactly what the code does. `node.h:78-79` confirms.

### Correction 2: Plan omits `last_auxin_flux` as the intermediate

The update pipeline is: `transport_with_children()` → accumulates `last_auxin_flux[child]` → `update_canalization()` reads it → updates both bias maps. The plan describes the end-state correctly but doesn't mention this intermediate recording map.

---

## Implementation Feasibility

The plan is implementable in roughly three independent steps:

### Step A: Replace age gate in `has_vasculature()` with bias threshold

In `vascular.cpp:18-25`, change:
```cpp
return n.age >= g.cambium_maturation_ticks;
```
to:
```cpp
if (!n.parent) return true;
auto it = n.parent->structural_flow_bias.find(const_cast<Node*>(&n));
return it != n.parent->structural_flow_bias.end() &&
       it->second >= g.vascular_conductance_threshold;
```
Add `vascular_conductance_threshold` to genome (suggested default: ~0.1, fires after ~20 ticks of sustained flow at `structural_growth_rate = 0.005`).

The same change applies to the guards in `stem_node.cpp:24` and `root_node.cpp:33`.

### Step B: Replace thickening formula

In `StemNode::thicken()` and `RootNode::thicken()`:
1. Remove the age guard
2. Add a `structural_bias_from_parent()` helper on `Node` that returns `parent->structural_flow_bias.count(this) ? parent->structural_flow_bias.at(this) : 0.0f`
3. Replace `effective_rate` computation:
   ```cpp
   float struct_bias = structural_bias_from_parent();
   if (struct_bias < g.vascular_conductance_threshold) return;
   float effective_rate = g.cambium_responsiveness * struct_bias;
   ```
4. Keep the sugar gate (`sugar_gf`) and stress boost.

Remove from genome: `thickening_rate`, `auxin_thickening_threshold`, `cambium_maturation_ticks`, `root_cambium_maturation_ticks`.

### Step C: Augment `pipe_capacity()` with structural bias (optional, lower priority)

`VascNodeInfo` would need a `structural_bias` field. In `build_flat()`, look up `parent->structural_flow_bias[node]` for each non-root node and store it. Then:
```cpp
static float pipe_capacity(const VascNodeInfo& info, float conductance) {
    float base = 3.14159f * info.node->radius * info.node->radius * conductance;
    return base * (1.0f + info.structural_bias);
}
```
This is the most invasive change (modifies the flat-array build) but is well-isolated.

### Threshold calibration note

With `structural_growth_rate = 0.005/tick` and `structural_threshold = 0.15` (minimum flux to count):
- At 1 tick/hour, threshold 0.1 fires at ~20 hours
- Threshold 0.25 fires at ~50 hours (~2 days)
- Threshold 1.0 fires at ~200 hours (~8 days)

Current `cambium_maturation_ticks = 336 hours`. To preserve similar timing, set `vascular_conductance_threshold ≈ 0.5–1.0`. A lower value (0.1–0.2) would make vascular tissue develop faster on actively-growing branches, which might actually be more biologically realistic.

---

## Suggestions

1. **Add `structural_bias_from_parent()`** to `Node` before implementing Step B — it's a one-liner helper that avoids exposing the map internals in `stem_node.cpp` and `root_node.cpp`.

2. **Do Step A before Step B** — the maturation gate change can be deployed and tested independently (thickening will still use old formula temporarily, but vascular bulk flow will fire on bias instead of age). This makes bisecting easier.

3. **Keep `structural_max = 2.0` as-is** — the cap means the bias-based thickening rate will plateau at a finite maximum, which is correct biological behavior (trunk can't thicken infinitely fast).

4. **Step C is optional for correctness** — vascular bulk flow already routes proportionally to demand, which is mostly correct. Bias-scaling of `pipe_capacity` is a refinement that rewards established connections, but the first two steps are the load-bearing changes.
