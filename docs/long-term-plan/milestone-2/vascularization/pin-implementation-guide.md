# PIN Transport — Implementation Guide

**For:** An implementing agent with no prior context on this codebase.  
**Goal:** Replace the `structural_flow_bias` canalization ratchet with a PIN-protein-based
auxin transport pass. Three separate systems result: vascular (sugar/water/cytokinin),
PIN (auxin long-range), and local diffusion (auxin/GA/stress short-range).

**Design sources (read before touching code):**
- `docs/long-term-plan/milestone-2/vascularization/pin-transport-design.md`
- `docs/long-term-plan/milestone-2/vascularization/pin-transport-review.md`

**Build command:** `/usr/local/bin/cmake --build build` — run after every step.

---

## Critical Bug Fixes (from review doc)

Two bugs in the design doc pseudocode must be corrected before implementing:

**Bug 1 (`node.cpp:346`):** `last_auxin_flux.clear()` runs at the START of
`transport_with_children()`. If PIN writes flux into `last_auxin_flux` first and diffusion
runs after, the clear wipes PIN's data. Fix: remove the clear from `transport_with_children()`
and add it to the END of `update_canalization()` instead.

**Bug 2 (seed junction):** The design doc says the seed "distributes from `chemical(Auxin)`",
but shoot-side auxin arrives into `seed->transport_received` and isn't flushed until end of
`seed->tick()`. The seed would distribute last tick's auxin. Fix: accumulate incoming shoot
auxin into a local `seed_collected` variable during Phase A, then distribute from that in
Phase B. See exact code in Step 1 below.

---

## Step 1 — Add `pin_transport()` (new file)

### 1a. Create `src/engine/pin_transport.h`

```cpp
#pragma once
namespace botany { class Plant; struct Genome; }
namespace botany {
    void pin_transport(Plant& plant, const Genome& g);
}
```

### 1b. Create `src/engine/pin_transport.cpp`

```cpp
// src/engine/pin_transport.cpp — PIN-mediated polar auxin transport.
// Three phases per tick:
//   A. Shoot post-order: each shoot node pumps auxin toward its parent (basipetal).
//   B. Seed junction: collects from shoot children, distributes to root children by r².
//   C. Root pre-order: distributes auxin from seed toward root tips (acropetal).
// Runs after vascular_transport(), before the DFS tree walk.
#include "engine/pin_transport.h"
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/node/node.h"
#include "engine/chemical/chemical.h"
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>

namespace botany {

// --- Tree traversal helpers ---

static void collect_shoot_postorder(Node* n, std::vector<Node*>& out) {
    for (Node* child : n->children) {
        if (child->type != NodeType::ROOT && child->type != NodeType::ROOT_APICAL)
            collect_shoot_postorder(child, out);
    }
    out.push_back(n);
}

static void collect_root_preorder(Node* n, std::vector<Node*>& out) {
    out.push_back(n);
    for (Node* child : n->children) {
        if (child->type == NodeType::ROOT || child->type == NodeType::ROOT_APICAL)
            collect_root_preorder(child, out);
    }
}

// --- Main pass ---

void pin_transport(Plant& plant, const Genome& g) {
    Node* seed = plant.seed_mut();
    if (!seed) return;

    // Partition seed's children into shoot vs root sides.
    std::vector<Node*> shoot_children_of_seed;
    std::vector<Node*> root_children_of_seed;
    for (Node* c : seed->children) {
        if (c->type == NodeType::ROOT || c->type == NodeType::ROOT_APICAL)
            root_children_of_seed.push_back(c);
        else
            shoot_children_of_seed.push_back(c);
    }

    // ----------------------------------------------------------------
    // Phase A: Shoot post-order — leaves → seed.
    //   Each shoot node pumps auxin toward its parent.
    //   Anti-teleportation: use parent->transport_received for non-seed parents.
    //   Seed children use local accumulator (bug 2 fix — see review doc).
    // ----------------------------------------------------------------
    std::vector<Node*> shoot_nodes;
    for (Node* c : shoot_children_of_seed)
        collect_shoot_postorder(c, shoot_nodes);

    float seed_collected = 0.0f;

    for (Node* n : shoot_nodes) {
        Node* par = n->parent;
        if (!par) continue;  // should not happen for shoot nodes

        float r2 = n->radius * n->radius;
        float max_cap = r2 * g.pin_capacity_per_area;

        // Efficiency: cold-start floor + upregulation from canalization history.
        // Read last tick's auxin_flow_bias — PIN hasn't updated it yet this tick.
        float flow_bias = 0.0f;
        auto it = par->auxin_flow_bias.find(n);
        if (it != par->auxin_flow_bias.end()) flow_bias = it->second;
        float efficiency = g.pin_base_efficiency + flow_bias * (1.0f - g.pin_base_efficiency);

        float available = n->chemical(ChemicalID::Auxin);
        float moved = std::min(available, max_cap * efficiency);
        if (moved < 1e-8f) continue;

        n->chemical(ChemicalID::Auxin) -= moved;
        par->last_auxin_flux[n] += moved;  // record for update_canalization

        if (par == seed) {
            // Bug 2 fix: accumulate locally, not into transport_received.
            // Seed's transport_received isn't flushed until end of seed->tick().
            seed_collected += moved;
        } else {
            par->transport_received[ChemicalID::Auxin] += moved;
        }
    }

    // ----------------------------------------------------------------
    // Phase B: Seed junction — distribute this tick's shoot collection
    //   to root children weighted by radius (thicker root = larger PIN area).
    //   Seed is always fully canalized — no efficiency floor applies here.
    // ----------------------------------------------------------------
    // root_forwarded: how much auxin is heading to each root node this tick.
    std::unordered_map<Node*, float> root_forwarded;

    if (seed_collected > 1e-8f && !root_children_of_seed.empty()) {
        float total_r = 0.0f;
        for (Node* rc : root_children_of_seed) total_r += rc->radius;

        for (Node* rc : root_children_of_seed) {
            float share = (total_r > 1e-8f)
                ? rc->radius / total_r
                : 1.0f / static_cast<float>(root_children_of_seed.size());
            float max_cap = rc->radius * rc->radius * g.pin_capacity_per_area;
            float to_send = std::min(seed_collected * share, max_cap);
            root_forwarded[rc] = to_send;
            seed->last_auxin_flux[rc] += to_send;  // record for canalization
            seed_collected -= to_send;
        }
    }
    seed->chemical(ChemicalID::Auxin) += seed_collected;  // remainder stays at seed

    // ----------------------------------------------------------------
    // Phase C: Root pre-order — seed → root tips.
    //   Each root node receives its forwarded auxin and passes the rest
    //   deeper. Writes directly to chemical(Auxin) so root update_tissue()
    //   sees PIN-delivered auxin this tick (PIN runs before the DFS walk).
    // ----------------------------------------------------------------
    std::vector<Node*> root_nodes;
    for (Node* c : root_children_of_seed)
        collect_root_preorder(c, root_nodes);

    for (Node* n : root_nodes) {
        auto it = root_forwarded.find(n);
        if (it == root_forwarded.end()) continue;
        float incoming = it->second;
        if (incoming < 1e-8f) continue;

        // Collect root children of this node.
        std::vector<Node*> root_children;
        float total_r = 0.0f;
        for (Node* child : n->children) {
            if (child->type == NodeType::ROOT || child->type == NodeType::ROOT_APICAL) {
                root_children.push_back(child);
                total_r += child->radius;
            }
        }

        // Distribute to root children; remainder stays at this node.
        float distributed = 0.0f;
        for (Node* rc : root_children) {
            float share = (total_r > 1e-8f)
                ? rc->radius / total_r
                : 1.0f / static_cast<float>(root_children.size());
            float max_cap = rc->radius * rc->radius * g.pin_capacity_per_area;
            float to_send = std::min(incoming * share, max_cap);
            root_forwarded[rc] += to_send;
            n->last_auxin_flux[rc] += to_send;  // record for canalization
            distributed += to_send;
        }

        // Deposit: remainder (incoming − distributed) stays at this node.
        n->chemical(ChemicalID::Auxin) += (incoming - distributed);
    }
}

} // namespace botany
```

---

## Step 2 — Rewrite `update_canalization()` and move `last_auxin_flux.clear()`

### 2a. Remove `last_auxin_flux.clear()` from `transport_with_children()`

**File:** `src/engine/node/node.cpp`  
**Line:** 346 — first line inside `transport_with_children()`.

Delete this line:
```cpp
last_auxin_flux.clear();
```

### 2b. Rewrite `update_canalization()` (lines 560–578)

Replace the entire body with the saturation-based lerp. The new formula computes
instantaneous PIN saturation (flux / capacity) and smooths it exponentially.

**Old body (`node.cpp:560–578`):**
```cpp
void Node::update_canalization(const Genome& g) {
    for (Node* child : children) {
        float flux = 0.0f;
        auto it = last_auxin_flux.find(child);
        if (it != last_auxin_flux.end()) flux = it->second;

        float target = flux * g.transient_gain;
        float& flow_bias = auxin_flow_bias[child];
        flow_bias += (target - flow_bias) * g.transient_rate;

        float& struct_bias = structural_flow_bias[child];
        if (flux > g.structural_threshold) {
            struct_bias += g.structural_growth_rate;
        }
        struct_bias = std::min(struct_bias, g.structural_max);
    }
}
```

**New body:**
```cpp
void Node::update_canalization(const Genome& g) {
    for (Node* child : children) {
        float flux = 0.0f;
        auto it = last_auxin_flux.find(child);
        if (it != last_auxin_flux.end()) flux = it->second;

        // Saturation = fraction of this connection's PIN capacity currently used.
        // r² is the child's cross-sectional area; pin_capacity_per_area is the max
        // flux density. Clamp to [0, 1] — physical PIN can't exceed full saturation.
        float r2 = child->radius * child->radius;
        float capacity = r2 * g.pin_capacity_per_area;
        float saturation = (capacity > 1e-8f)
            ? std::min(flux / capacity, 1.0f) : 0.0f;

        // Exponential smoothing: bias chases saturation at smoothing_rate.
        // Natural decay: when flux stops, saturation→0 and bias lerps toward 0
        // over ~1/smoothing_rate ticks. No separate decay param needed.
        float& flow_bias = auxin_flow_bias[child];
        flow_bias += (saturation - flow_bias) * g.smoothing_rate;
    }

    // Bug 1 fix: clear AFTER update_canalization reads, not at start of transport.
    // PIN fills last_auxin_flux before diffusion runs; both must accumulate
    // into the same map before canalization reads it.
    last_auxin_flux.clear();
}
```

---

## Step 3 — Delete `structural_flow_bias` (all 14 sites)

Work through these in order. Build after each sub-step.

### 3a. `src/engine/node/node.h`

**Line 79** — delete the field declaration:
```cpp
std::unordered_map<Node*, float> structural_flow_bias;  // persistent — slow, permanent
```

**Lines 84–86** — delete the accessor comment and declaration:
```cpp
// Returns the structural_flow_bias this node's parent has recorded for it.
// Stored on the parent keyed by child pointer — zero if no parent or no entry yet.
float get_parent_structural_bias() const;
```

### 3b. `src/engine/node/node.cpp`

**Lines 57–60** — delete the structural bias transfer block inside `replace_child()`:
```cpp
auto it_struct = structural_flow_bias.find(old_child);
if (it_struct != structural_flow_bias.end()) {
    structural_flow_bias[new_child] = it_struct->second;
    structural_flow_bias.erase(it_struct);
}
```

**Lines 284–298** — delete the entire `get_parent_structural_bias()` method:
```cpp
float Node::get_parent_structural_bias() const {
    if (!parent) {
        float max_bias = 0.0f;
        for (const auto& [child, bias] : structural_flow_bias)
            max_bias = std::max(max_bias, bias);
        return max_bias;
    }
    auto it = parent->structural_flow_bias.find(const_cast<Node*>(this));
    if (it == parent->structural_flow_bias.end()) return 0.0f;
    return it->second;
}
```

**Lines 300–307** — update `get_bias_multiplier()` to remove the structural term:
```cpp
// Old:
float Node::get_bias_multiplier(Node* child, const Genome& g) const {
    float flow = 0.0f, structural = 0.0f;
    auto it_f = auxin_flow_bias.find(child);
    if (it_f != auxin_flow_bias.end()) flow = it_f->second;
    auto it_s = structural_flow_bias.find(child);
    if (it_s != structural_flow_bias.end()) structural = it_s->second;
    return 1.0f + g.canalization_weight * (flow + structural);
}

// New:
float Node::get_bias_multiplier(Node* child, const Genome& g) const {
    float flow = 0.0f;
    auto it_f = auxin_flow_bias.find(child);
    if (it_f != auxin_flow_bias.end()) flow = it_f->second;
    return 1.0f + g.canalization_weight * flow;
}
```

**Line 313** — delete inside `die()` (the structural bias erase):
```cpp
parent->structural_flow_bias.erase(this);
```

### 3c. `src/engine/node/stem_node.cpp`

**Lines 22–38** — replace the bias-reading block in `StemNode::thicken()`:

```cpp
// Old — reads structural_flow_bias via get_parent_structural_bias():
float bias;
if (!parent) {
    bias = 0.0f;
    for (auto& [child, b] : structural_flow_bias) {
        bias = std::max(bias, b);
    }
} else {
    bias = get_parent_structural_bias();
}

// New — reads auxin_flow_bias (PIN saturation, transient):
float bias = 0.0f;
if (!parent) {
    // Seed node: no parent records its bias. Use max of recorded child biases —
    // the seed is the transit junction, its thickening reflects the dominant path.
    for (auto& [child, b] : auxin_flow_bias)
        bias = std::max(bias, b);
} else {
    auto it = parent->auxin_flow_bias.find(const_cast<Node*>(this));
    if (it != parent->auxin_flow_bias.end()) bias = it->second;
}
```

Also update the comment at line 22 to read:
```cpp
// Cambium activity is driven by PIN auxin flux saturation (auxin_flow_bias), not age.
// Connections carrying high auxin flux have high saturation → fast thickening.
// Connections with low or zero flux stay thin (auxin_flow_bias decays toward 0).
// Stress hormone (thigmomorphogenesis) provides an additional boost.
```

### 3d. `src/engine/node/root_node.cpp`

**Lines 32–38** — replace in `RootNode::thicken()`:

```cpp
// Old:
// Same vascular-driven model as StemNode. Root connections accumulate
// structural_flow_bias from sugar and cytokinin transport (not auxin —
// real root polar auxin transport governs patterning/gravitropism, not
// cambial signaling). Well-used root connections thicken; lateral roots
// with low flux stay thin.
float bias = get_parent_structural_bias();

// New:
// Same PIN-driven model as StemNode. Shoot-derived auxin flows acropetally
// through root internodes via the PIN pass, building auxin_flow_bias on
// root connections. Well-used root connections thicken; lateral roots
// with low PIN flux stay thin.
float bias = 0.0f;
if (parent) {
    auto it = parent->auxin_flow_bias.find(const_cast<Node*>(this));
    if (it != parent->auxin_flow_bias.end()) bias = it->second;
}
```

### 3e. `src/engine/node/tissues/apical.cpp`

**Lines 185–191** — delete the structural bias stamping block inside `spawn_internode()`:
```cpp
if (internode->parent) {
    float internode_length = glm::length(internode->offset);
    float initial_bias = 0.01f + internode_length * 0.1f;
    float& entry = internode->parent->structural_flow_bias[internode];
    entry = std::max(entry, initial_bias);
}
```

No replacement. The first PIN flux through this connection builds `auxin_flow_bias`
immediately via `update_canalization`. No bootstrap stamp is needed.

### 3f. `src/engine/node/tissues/root_apical.cpp`

**Lines 103–111** — delete the identical stamping block inside `spawn_internode()`:
```cpp
if (internode->parent) {
    float internode_length = glm::length(internode->offset);
    float initial_bias = 0.01f + internode_length * 0.1f;
    float& entry = internode->parent->structural_flow_bias[internode];
    entry = std::max(entry, initial_bias);
}
```

No replacement for the same reason.

### 3g. `src/engine/vascular.cpp`

**Lines 19–32** — replace `has_vasculature()` entirely:

```cpp
// Old:
bool has_vasculature(const Node& n, const Genome& g) {
    if (!n.parent) return true;  // seed is always a vascular junction
    if (n.type == NodeType::STEM || n.type == NodeType::ROOT) {
        auto it = n.parent->structural_flow_bias.find(const_cast<Node*>(&n));
        if (it == n.parent->structural_flow_bias.end()) return false;
        return it->second >= g.vascular_conductance_threshold;
    }
    return false;
}

// New — admission by radius; all stems ≥ vascular_radius_threshold qualify.
// initial_radius = 0.015 dm > vascular_radius_threshold = 0.01 dm, so newly
// spawned internodes join immediately. Leaves and meristems never join (they
// rely on last-mile local diffusion). No structural bias needed.
bool has_vasculature(const Node& n, const Genome& g) {
    if (!n.parent) return true;  // seed is always a vascular junction
    if (n.type == NodeType::STEM || n.type == NodeType::ROOT)
        return n.radius >= g.vascular_radius_threshold;
    return false;
}
```

**Lines 269–271** — remove structural bias augmentation from vascular distribution weights:

```cpp
// Old:
auto it = info.node->structural_flow_bias.find(flat[ci].node);
float bias = (it != info.node->structural_flow_bias.end()) ? it->second : 0.0f;
weights[k] = cap * (1.0f + g.canalization_weight * bias);

// New — weight by pipe cross-section only (r² × conductance).
// Canalization is preserved: thicker stems still get proportionally more flow.
// The double-weighting (radius AND structural bias) is removed by design.
weights[k] = cap;
```

**Lines 323–324** — update debug log to remove the structural bias lookup and CSV column:

```cpp
// Old (inside the if (log) block after alloc):
float bias = 0.0f;
auto it = info.node->structural_flow_bias.find(child);
if (it != info.node->structural_flow_bias.end()) bias = it->second;
*log << current_tick << ','
     << info.node->id << ',' << child->id << ','
     << node_type_name(child->type) << ','
     << chem_name(chem_id) << ','
     << ceilings[k] << ','
     << weights[k] << ','
     << alloc[k] << ','
     << bias << '\n';

// New — drop structural_flow_bias column:
*log << current_tick << ','
     << info.node->id << ',' << child->id << ','
     << node_type_name(child->type) << ','
     << chem_name(chem_id) << ','
     << ceilings[k] << ','
     << weights[k] << ','
     << alloc[k] << '\n';
```

Also update the CSV header at line 362 — remove `structural_flow_bias` from the column list:
```cpp
// Old:
log_file << "tick,junction_node_id,child_node_id,child_type,"
            "chemical,demand,conductance_weight,delivered,"
            "structural_flow_bias\n";

// New:
log_file << "tick,junction_node_id,child_node_id,child_type,"
            "chemical,demand,conductance_weight,delivered\n";
```

---

## Step 4 — Update `genome.h`

### 4a. Update the `cambium_responsiveness` comment (line 56–57)

```cpp
// Old:
float cambium_responsiveness;     // dm/hr·bias — thickening rate per unit structural_flow_bias.
                                  // delta_radius = cambium_responsiveness × structural_flow_bias × sugar_gf × stress_boost.

// New:
float cambium_responsiveness;     // dm/hr·bias — thickening rate per unit auxin_flow_bias (PIN saturation 0→1).
                                  // delta_radius = cambium_responsiveness × auxin_flow_bias × sugar_gf × stress_boost.
```

### 4b. Replace the canalization params block (lines 153–165)

```cpp
// Old:
// Canalization — auxin flow history biases transport
float transient_gain;                 // target bias per unit of auxin flux
float transient_rate;                 // how fast transient bias chases its target (0-1)
float structural_threshold;           // minimum auxin flux to grow structural bias
float structural_growth_rate;         // structural bias increment per tick above threshold
float structural_max;                 // cap on structural bias
float canalization_weight;            // global scaling on combined bias effect (0 = disabled)

// Vascular transport
float xylem_conductance;
float phloem_conductance;
float phloem_reserve_fraction;
float vascular_conductance_threshold; // minimum structural_flow_bias for bulk vascular admission.
                                      // ...

// New:
// Canalization — PIN transport history biases auxin flow and cambium
float smoothing_rate;                 // lerp rate for auxin_flow_bias toward current PIN saturation (~20 tick response)
float canalization_weight;            // global scaling on auxin_flow_bias effect in local diffusion (0 = disabled)

// PIN transport
float pin_capacity_per_area;          // AU/(dm²·tick) — max auxin transport per unit cross-section at full efficiency.
                                      // Also the denominator in saturation = flux / (r² × pin_capacity_per_area).
float pin_base_efficiency;            // [0–1] — cold-start PIN efficiency when auxin_flow_bias = 0 (constitutively active PINs)

// Vascular transport
float xylem_conductance;
float phloem_conductance;
float phloem_reserve_fraction;
float vascular_radius_threshold;      // dm — minimum radius for bulk vascular admission.
                                      // Set below initial_radius (0.015 dm) so newly spawned
                                      // internodes qualify from birth. Only pre-initial-radius
                                      // manually-created nodes stay excluded.
```

### 4c. Update `default_genome()` initializers

Remove:
```cpp
.transient_gain = 2.0f,
.transient_rate = 0.2f,
.structural_threshold = 0.15f,
.structural_growth_rate = 0.005f,
.structural_max = 2.0f,
// ...
.vascular_conductance_threshold = 0.005f,
```

Add:
```cpp
.smoothing_rate = 0.1f,            // 10% per tick → ~20 tick response time; natural decay when flux stops
.canalization_weight = 1.0f,
.pin_capacity_per_area = 500.0f,   // AU/(dm²·tick) — tuning knob; adjust up to reduce sensitivity
.pin_base_efficiency = 0.2f,       // cold-start: new stems transport 20% of capacity immediately
.vascular_radius_threshold = 0.01f, // dm — below initial_radius (0.015), all new internodes qualify
```

**Calibration note:** With `initial_radius = 0.015 dm`, a new stem's cold-start capacity is
`0.015² × 500 × 0.2 = 0.0225 AU/tick` — about 15% of the apical's baseline production
(`0.15 AU/tick`). Auxin will pool at the tip briefly until `auxin_flow_bias` builds.
Once fully upregulated: `0.015² × 500 × 1.0 = 0.1125 AU/tick` — 75% of apical output.
Tune `pin_capacity_per_area` upward to reduce cambium sensitivity, downward to increase it.

---

## Step 5 — Wire into `plant.cpp`

**File:** `src/engine/plant.cpp`

Add include at the top of the file:
```cpp
#include "engine/pin_transport.h"
```

**Line 144** — `tick_tree()` currently reads:
```cpp
void Plant::tick_tree(const WorldParams& world, PerfStats* /*perf*/) {
    vascular_transport(*this, genome_, world);  // bulk flow for sugar/water/cytokinin
    tick_recursive(*nodes_[0], *this, world);
    flush_removals();
}
```

Add the `pin_transport` call between vascular and DFS:
```cpp
void Plant::tick_tree(const WorldParams& world, PerfStats* /*perf*/) {
    vascular_transport(*this, genome_, world);  // bulk flow for sugar/water/cytokinin
    pin_transport(*this, genome_);              // PIN-mediated polar auxin transport
    tick_recursive(*nodes_[0], *this, world);
    flush_removals();
}
```

---

## Step 6 — Register `pin_transport.cpp` in CMakeLists.txt

**File:** `CMakeLists.txt`  
**Lines 69–84** — add `src/engine/pin_transport.cpp` to `botany_engine`:

```cmake
add_library(botany_engine STATIC
    src/engine/node/node.cpp
    src/engine/node/stem_node.cpp
    src/engine/node/root_node.cpp
    src/engine/node/tissues/leaf.cpp
    src/engine/node/tissues/apical.cpp
    src/engine/node/tissues/root_apical.cpp
    src/engine/plant.cpp
    src/engine/pin_transport.cpp      # ← ADD
    src/engine/sugar.cpp
    src/engine/light.cpp
    src/engine/debug_log.cpp
    src/engine/ethylene.cpp
    src/engine/vascular.cpp
    src/engine/engine.cpp
    src/engine/world_params.cpp
)
```

Build now and fix any compilation errors before proceeding to tests.

---

## Step 7 — Update tests

Three existing test files use `structural_flow_bias` and will fail to compile.
One new test file is needed.

### 7a. `tests/test_node.cpp` (8 uses)

Grep for `structural_flow_bias` in this file. For each:
- If the test sets `parent->structural_flow_bias[child] = X` to bootstrap thickening:
  replace with `parent->auxin_flow_bias[child] = X` (same map semantics, different type).
- If the test asserts `structural_flow_bias` grew: assert `auxin_flow_bias` instead.
- Remove any assertion that `structural_flow_bias` grows monotonically — it no longer does
  (it decays when flux stops). Assert that `auxin_flow_bias > 0` after auxin flux, instead.

### 7b. `tests/test_meristem.cpp` (14 uses)

Same pattern. Specifically:
- Tests that pre-populate `structural_flow_bias` to bootstrap thickening: switch to
  `auxin_flow_bias`.
- Test named "Thickening proportional to `structural_flow_bias`": rename to
  "Thickening proportional to `auxin_flow_bias`" and verify it still passes by injecting
  auxin into the meristem and letting the PIN pass build `auxin_flow_bias` organically,
  OR inject `auxin_flow_bias` directly if the genome has PIN disabled.
- Remove any test that checks the ratchet behavior (bias only increases). The new bias
  lerps both up and down.

### 7c. `tests/test_vascularization.cpp` (15 uses — all 4 tests need rewriting)

The 4 tests were written around `structural_flow_bias`. Rewrite their intent:

| Old test | New intent |
|----------|-----------|
| "Zero structural_flow_bias → zero thickening" | "Zero auxin_flow_bias → zero thickening": create stem, don't inject auxin, verify radius stays at initial_radius after N ticks |
| "Bias-proportional growth rate" | "auxin_flow_bias-proportional thickening": inject known auxin_flow_bias values, verify `delta_radius = cambium_responsiveness × bias × sugar_gf × stress_boost` |
| "Canalization ratchet: auxin flux builds structural_flow_bias" | "PIN flux builds auxin_flow_bias": run plant for N ticks with active apical; verify auxin_flow_bias > 0 on main axis connections |
| "Conductance-weighted phloem distribution" | Keep intent; remove `structural_flow_bias` from setup; distribution is now pure-r² (pipe capacity only) |

### 7d. Create `tests/test_pin_transport.cpp`

```cpp
// tests/test_pin_transport.cpp — Integration tests for PIN polar auxin transport.
#include <catch2/catch_test_macros.hpp>
#include "engine/plant.h"
#include "engine/world_params.h"
#include "engine/pin_transport.h"

using namespace botany;

static WorldParams static_world() {
    WorldParams w = default_world_params();
    w.starvation_ticks_max = 1000000u;
    w.light_level = 0.0f;
    return w;
}

TEST_CASE("PIN transport: auxin crosses the seed junction to root side") {
    Genome g = default_genome();
    g.shoot_plastochron = 1000000u;  // no node spawning
    g.root_plastochron  = 1000000u;
    WorldParams w = static_world();

    Plant plant(g, glm::vec3(0));

    // Run for 50 ticks to establish flow
    for (int i = 0; i < 50; ++i) plant.tick(w);

    // Find root apical
    float root_apical_auxin = 0.0f;
    float seed_auxin = 0.0f;
    float apical_auxin = 0.0f;
    plant.for_each_node([&](const Node& n) {
        if (!n.parent) seed_auxin = n.chemical(ChemicalID::Auxin);
        if (n.type == NodeType::APICAL) apical_auxin = n.chemical(ChemicalID::Auxin);
        if (n.type == NodeType::ROOT_APICAL) root_apical_auxin = n.chemical(ChemicalID::Auxin);
    });

    // Gradient direction: apex > seed >= root apical (basipetal)
    REQUIRE(apical_auxin > seed_auxin);
    // Root side receives PIN-transported auxin (was zero before PIN)
    REQUIRE(root_apical_auxin > 0.0f);
}

TEST_CASE("PIN transport: auxin_flow_bias builds on main axis") {
    Genome g = default_genome();
    g.shoot_plastochron = 1000000u;
    g.root_plastochron  = 1000000u;
    WorldParams w = static_world();

    Plant plant(g, glm::vec3(0));
    for (int i = 0; i < 100; ++i) plant.tick(w);

    // Seed should have non-zero auxin_flow_bias for its shoot child
    bool found_bias = false;
    plant.for_each_node([&](const Node& n) {
        if (!n.parent) {  // seed
            for (const auto& [child, bias] : n.auxin_flow_bias) {
                if (bias > 0.0f) found_bias = true;
            }
        }
    });
    REQUIRE(found_bias);
}
```

Register it in `CMakeLists.txt` alongside the existing test files.

---

## Step 8 — Also update `src/evolution/genome_bridge.cpp`

Search for all `structural_flow_bias`-adjacent genome field names in `genome_bridge.cpp`:
`transient_gain`, `transient_rate`, `structural_threshold`, `structural_growth_rate`,
`structural_max`, `vascular_conductance_threshold`.

For each:
- Remove from `build_genome_template()` / `to_structured()` / `from_structured()`.
- Add entries for `smoothing_rate`, `pin_capacity_per_area`, `pin_base_efficiency`,
  `vascular_radius_threshold` with sensible mutation ranges:
  - `smoothing_rate`: range [0.01, 0.5], strength ~20%
  - `pin_capacity_per_area`: range [50, 2000], strength ~20%
  - `pin_base_efficiency`: range [0.05, 0.5], strength ~20%
  - `vascular_radius_threshold`: range [0.001, 0.05], strength ~20%

---

## Commit Strategy

### Commit 1 — New PIN pass + genome params

Stage: `src/engine/pin_transport.h`, `src/engine/pin_transport.cpp`,
`src/engine/genome.h`, `src/engine/plant.cpp`, `src/engine/node/node.cpp`
(both the `last_auxin_flux.clear()` move and the `update_canalization` rewrite),
`CMakeLists.txt`.

Message: `feat: add pin_transport pass, rewrite update_canalization to saturation lerp`

The code compiles at this point; `structural_flow_bias` still exists, so thickening
still uses the old driver. Tests pass (old tests still run; new ones not yet added).

### Commit 2 — Delete `structural_flow_bias` everywhere

Stage all modified files from Step 3 (`node.h`, `node.cpp`, `stem_node.cpp`,
`root_node.cpp`, `apical.cpp`, `root_apical.cpp`, `vascular.cpp`).
Also stage `src/evolution/genome_bridge.cpp` changes from Step 8.

Message: `refactor: delete structural_flow_bias, switch thickening to auxin_flow_bias`

Build and verify all tests still pass before committing.

### Commit 3 — Test updates

Stage: updated `test_node.cpp`, `test_meristem.cpp`, `test_vascularization.cpp`,
new `test_pin_transport.cpp`, and `CMakeLists.txt` test registration.

Message: `test: update vascularization/meristem tests for PIN model, add test_pin_transport`

---

## Verification Checklist

After all commits, run `./build/botany_realtime --color auxin` and confirm:

- [ ] Bright auxin at shoot tips, smooth gradient declining over 15–20 nodes
- [ ] Non-zero (dim) auxin visible at seed and upper root internodes
- [ ] Root apical meristems show non-zero auxin (full loop connected)
- [ ] Lateral branches show lower auxin than main axis at same height
- [ ] Run `./build/botany_tests` — all tests pass

Run `grep -r "structural_flow_bias" src/ tests/` — should return zero matches.
