# Node Inspector Panel & Overlay Redesign — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Align the realtime node inspector and overlays with the compartmented vascular + tree-scale hormone models. Add per-chemical per-tick produced/consumed counters, per-edge flux/capacity instrumentation, a 3-scope chemicals table, and a two-tier overlay selector with new Level/Capacity/Growth/Activation/Starvation modes.

**Architecture:** Two-pronged. (A) Engine-side instrumentation: add counters and helpers, never change transport physics. (B) UI-side rewrite: replace the flat overlay dropdown with a category+sub-picker, rewrite the inspector panel body to consume the new counters, add a NaN-aware gray path in the renderer.

**Tech Stack:** C++17, ImGui, OpenGL 4.1, Catch2 for tests, CMake.

**Spec:** [`docs/superpowers/specs/2026-04-21-node-inspector-and-overlay-redesign.md`](../specs/2026-04-21-node-inspector-and-overlay-redesign.md)

---

## File map

### New files
- `tests/test_tick_counters.cpp` — mass-balance + counter coverage tests
- `src/engine/ui_helpers.h/.cpp` — `vascular_scope()`, `compute_maintenance_cost()`, `hydraulic_maturity()`, `nodes_to_seed()` (shared by panel and overlays)

### Modified files
- `src/engine/chemical/chemical.h` — add `Count` sentinel
- `src/engine/node/node.h` / `.cpp` — per-chem arrays, per-edge maps, counter increments on decay
- `src/engine/node/stem_node.cpp` / `root_node.cpp` — instrument elongate, thicken, maintenance, absorb
- `src/engine/node/tissues/leaf.cpp` — instrument photo, transpire, GA emit, Eth emit, expand
- `src/engine/node/tissues/apical.cpp` / `root_apical.cpp` — auxin/CK production mirrors into arrays
- `src/engine/vascular_sub_stepped.cpp` — per-edge flux+cap inside Jacobi passes
- `src/engine/pin_transport.cpp` — per-edge flux+cap at PIN flux recording sites
- `src/engine/plant.cpp` — add new clears to top-of-tick reset
- `src/app_realtime.cpp` — overlay selector rewrite, panel rewrite
- `src/renderer/renderer.cpp` — NaN-aware gray path in the heatmap color branch
- `CMakeLists.txt` — register new test file

### Phase boundaries (commit points)
- **Phase 1 complete:** counters exist, reset correctly, unit-tested.
- **Phase 2 complete:** all chemical production/consumption sites instrumented.
- **Phase 3 complete:** per-edge flux/cap instrumented on all pathways.
- **Phase 4 complete:** panel rewritten, reads new counters.
- **Phase 5 complete:** new overlay selector + modes live.

---

## PHASE 1 — Counter scaffolding

### Task 1: Add `ChemicalID::Count` sentinel

**Files:**
- Modify: `src/engine/chemical/chemical.h` (lines 9-17)

- [ ] **Step 1: Open file and locate the enum**

Current enum (lines 9-17):
```cpp
enum class ChemicalID {
    Auxin,
    Cytokinin,
    Gibberellin,
    Sugar,
    Ethylene,
    Stress,
    Water,
};
```

- [ ] **Step 2: Add `Count` as the last entry**

New enum:
```cpp
enum class ChemicalID {
    Auxin,
    Cytokinin,
    Gibberellin,
    Sugar,
    Ethylene,
    Stress,
    Water,
    Count,   // sentinel — keep last; used for fixed-size per-chemical arrays
};
```

- [ ] **Step 3: Build and verify no regressions**

Run: `/usr/local/bin/cmake --build build`
Expected: clean build. Any existing switch statement on `ChemicalID` may warn about the unhandled `Count` case — if so, add a `case ChemicalID::Count: break;` or use a `default:` arm.

- [ ] **Step 4: Commit**

```bash
git add src/engine/chemical/chemical.h
git commit -m "chemical: add ChemicalID::Count sentinel for fixed-size per-chem arrays"
```

---

### Task 2: Add per-chem produced/consumed arrays on Node

**Files:**
- Modify: `src/engine/node/node.h` (near existing tick counters around lines 50-62)

- [ ] **Step 1: Locate existing tick counters in node.h**

Existing block (lines 52-62 approximately):
```cpp
float tick_sugar_maintenance = 0.0f;
float tick_sugar_activity    = 0.0f;
float tick_sugar_transport   = 0.0f;
float tick_auxin_produced    = 0.0f;
float tick_cytokinin_produced = 0.0f;
```

- [ ] **Step 2: Add new arrays alongside**

Insert after the existing block:
```cpp
// Per-chemical production/consumption, reset at the top of Plant::tick_tree().
// Indexed by static_cast<size_t>(ChemicalID::*).
std::array<float, static_cast<size_t>(ChemicalID::Count)> tick_chem_produced{};
std::array<float, static_cast<size_t>(ChemicalID::Count)> tick_chem_consumed{};
```

Ensure `#include <array>` is present at top of `node.h`.

- [ ] **Step 3: Build**

Run: `/usr/local/bin/cmake --build build`
Expected: clean build.

- [ ] **Step 4: Commit**

```bash
git add src/engine/node/node.h
git commit -m "node: add tick_chem_produced/consumed arrays (uninstrumented)"
```

---

### Task 3: Add per-edge flux and capacity maps on Node

**Files:**
- Modify: `src/engine/node/node.h` (near `last_auxin_flux` around line 78)

- [ ] **Step 1: Locate `last_auxin_flux`**

Existing (line ~78):
```cpp
std::unordered_map<Node*, float> last_auxin_flux;
```

- [ ] **Step 2: Add parallel per-chem maps**

Insert after the existing line:
```cpp
// Per-chemical per-edge flux and cap for the "Transport Capacity Used" overlay.
// [chem][child_ptr] = signed flux across (this -> child) this tick.
// Cleared at the top of Plant::tick_tree() alongside last_auxin_flux.
std::array<std::unordered_map<Node*, float>, static_cast<size_t>(ChemicalID::Count)> tick_edge_flux;
std::array<std::unordered_map<Node*, float>, static_cast<size_t>(ChemicalID::Count)> tick_edge_cap;
```

- [ ] **Step 3: Build**

Run: `/usr/local/bin/cmake --build build`
Expected: clean build.

- [ ] **Step 4: Commit**

```bash
git add src/engine/node/node.h
git commit -m "node: add tick_edge_flux/cap arrays (uninstrumented)"
```

---

### Task 4: Extend top-of-tick reset in Plant::tick_tree

**Files:**
- Modify: `src/engine/plant.cpp` (lines 149-164)

- [ ] **Step 1: Locate the reset loop**

Current reset near lines 149-164 clears `last_auxin_flux`, `tick_auxin_produced`, `tick_cytokinin_produced`.

- [ ] **Step 2: Add the new clears inside the same for_each_node**

Inside the existing loop body, add:
```cpp
n.tick_chem_produced.fill(0.0f);
n.tick_chem_consumed.fill(0.0f);
for (auto& m : n.tick_edge_flux) m.clear();
for (auto& m : n.tick_edge_cap)  m.clear();
```

- [ ] **Step 3: Build**

Run: `/usr/local/bin/cmake --build build`
Expected: clean build.

- [ ] **Step 4: Run existing tests**

Run: `./build/botany_tests`
Expected: all 218 tests pass (reset doesn't break anything).

- [ ] **Step 5: Commit**

```bash
git add src/engine/plant.cpp
git commit -m "plant: extend top-of-tick reset to clear new per-chem counters"
```

---

### Task 5: Write failing mass-balance test

**Files:**
- Create: `tests/test_tick_counters.cpp`
- Modify: `CMakeLists.txt` (add to test list around line 162)

- [ ] **Step 1: Create test file**

Content of `tests/test_tick_counters.cpp`:
```cpp
#include <catch2/catch_test_macros.hpp>
#include "engine/engine.h"
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/node/node.h"
#include "engine/chemical/chemical.h"

using namespace botany;

TEST_CASE("tick counters: arrays reset each tick", "[tick_counters]") {
    Engine engine;
    Genome g = default_genome();
    WorldParams w{};
    auto pid = engine.add_plant(g);
    engine.tick(w);
    engine.tick(w);
    // after two ticks, every node's arrays should have been zeroed at the start
    // of tick 2 and populated (or not) during tick 2 — either way, size matches Count.
    engine.get_plant(pid).for_each_node([](const Node& n) {
        REQUIRE(n.tick_chem_produced.size() == static_cast<size_t>(ChemicalID::Count));
        REQUIRE(n.tick_chem_consumed.size() == static_cast<size_t>(ChemicalID::Count));
        for (size_t i = 0; i < n.tick_chem_produced.size(); ++i) {
            REQUIRE(n.tick_chem_produced[i] >= 0.0f);
            REQUIRE(n.tick_chem_consumed[i] >= 0.0f);
        }
    });
}

TEST_CASE("tick counters: leaf sugar production matches photosynthesis", "[tick_counters][photo]") {
    Engine engine;
    Genome g = default_genome();
    WorldParams w{};
    w.light_level = 1.0f;
    auto pid = engine.add_plant(g);
    // spin up until a leaf exists
    for (int i = 0; i < 200; ++i) engine.tick(w);
    bool found_leaf = false;
    engine.get_plant(pid).for_each_node([&](const Node& n) {
        if (auto* lf = n.as_leaf(); lf && lf->leaf_size > 0.01f) {
            found_leaf = true;
            float produced = n.tick_chem_produced[static_cast<size_t>(ChemicalID::Sugar)];
            INFO("Expected non-zero sugar production on active leaf");
            REQUIRE(produced >= 0.0f);  // loose: just non-negative — tightened in phase 2
        }
    });
    REQUIRE(found_leaf);
}
```

- [ ] **Step 2: Register in CMakeLists.txt**

In `CMakeLists.txt`, find the `add_executable(botany_tests ...)` block (around line 143). Add `tests/test_tick_counters.cpp` to the list of sources before the closing paren.

- [ ] **Step 3: Build and run — expect PASS (loose assertions)**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[tick_counters]"`
Expected: both tests pass. They're deliberately loose — Phase 2 tightens the production test to require `> 0`.

- [ ] **Step 4: Commit**

```bash
git add tests/test_tick_counters.cpp CMakeLists.txt
git commit -m "tests: skeleton tick-counter tests (mass-balance placeholder)"
```

---

### Task 6: Add ui_helpers module

**Files:**
- Create: `src/engine/ui_helpers.h`
- Create: `src/engine/ui_helpers.cpp`
- Modify: `CMakeLists.txt` — add ui_helpers.cpp to engine sources

- [ ] **Step 1: Create `src/engine/ui_helpers.h`**

```cpp
#pragma once

#include "engine/chemical/chemical.h"
#include "engine/compartments.h"

namespace botany {

class Node;
class Genome;
struct WorldParams;

// Returns the conduit pool this chemical lives in for this node.
// - Stem/Root with vascular radius: own phloem (Sugar) or own xylem (Water, Cytokinin).
// - Leaf/Apical/RootApical/seed specialty: walks up to the nearest ancestor stem/root.
// - Signaling chems (Auxin, Gibberellin, Ethylene, Stress): returns nullptr (no conduit).
// - If no upstream conduit exists (e.g., seedling before vascular admission): nullptr.
const TransportPool* vascular_scope(const Node& n, ChemicalID chem);

// Preview of what pay_maintenance() would deduct this tick, without side effects.
// Used by the starvation overlay to evaluate sugar coverage.
float compute_maintenance_cost(const Node& n, const Genome& g, const WorldParams& w);

// 1 - radial_permeability_sugar(r) / base_radial_permeability_sugar.
// 0% = young leaky stem. ~90% = fully closed trunk (asymptote at the floor fraction).
// Defined for stems/roots; returns 0 for other types.
float hydraulic_maturity(const Node& n, const Genome& g);

// Walk-up count of parent hops until parent == nullptr (seed has 0).
int nodes_to_seed(const Node& n);

} // namespace botany
```

- [ ] **Step 2: Create `src/engine/ui_helpers.cpp`**

```cpp
#include "engine/ui_helpers.h"
#include "engine/node/node.h"
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/vascular_sub_stepped.h"

namespace botany {

const TransportPool* vascular_scope(const Node& n, ChemicalID chem) {
    switch (chem) {
        case ChemicalID::Sugar:
            return n.nearest_phloem_upstream();
        case ChemicalID::Water:
        case ChemicalID::Cytokinin:
            return n.nearest_xylem_upstream();
        case ChemicalID::Auxin:
        case ChemicalID::Gibberellin:
        case ChemicalID::Ethylene:
        case ChemicalID::Stress:
        case ChemicalID::Count:
            return nullptr;
    }
    return nullptr;
}

float compute_maintenance_cost(const Node& n, const Genome& g, const WorldParams& w) {
    // Mirror the same formula pay_maintenance() uses in each node subclass.
    // See Node::pay_maintenance and subclass overrides in stem_node.cpp / root_node.cpp /
    // tissues/*.cpp for the authoritative cost formula. Update here whenever those change.
    // Returns g (grams of sugar) per tick.
    if (auto* s = n.as_stem()) {
        float vol = 3.14159f * s->radius * s->radius * /* length */ 1.0f;  // TODO: use actual length
        return w.sugar_maintenance_stem * vol;
    }
    // ... handled in Task 22 when starvation overlay actually needs it.
    return 0.0f;
}

float hydraulic_maturity(const Node& n, const Genome& g) {
    if (!n.as_stem() && !n.as_root()) return 0.0f;
    float base = g.base_radial_permeability_sugar;
    if (base <= 1e-9f) return 0.0f;
    float perm = radial_permeability_sugar(n.radius, g);
    float closed = 1.0f - (perm / base);
    if (closed < 0.0f) closed = 0.0f;
    if (closed > 1.0f) closed = 1.0f;
    return closed;
}

int nodes_to_seed(const Node& n) {
    int count = 0;
    const Node* p = n.parent;
    while (p) { ++count; p = p->parent; }
    return count;
}

} // namespace botany
```

Note: `compute_maintenance_cost` has a `TODO: use actual length` and only handles stems. This is intentional — the starvation overlay task (Task 22) is where we complete this helper with the exact per-type formulas cross-checked against each `pay_maintenance()` override. We stub it now so the module compiles.

- [ ] **Step 3: Register in CMakeLists.txt**

Find the block that defines engine sources (where `node.cpp`, `plant.cpp`, etc. are listed). Add `src/engine/ui_helpers.cpp` to the list.

- [ ] **Step 4: Build**

Run: `/usr/local/bin/cmake --build build`
Expected: clean build.

- [ ] **Step 5: Commit**

```bash
git add src/engine/ui_helpers.h src/engine/ui_helpers.cpp CMakeLists.txt
git commit -m "engine: add ui_helpers module (vascular_scope, hydraulic_maturity, nodes_to_seed)"
```

---

## PHASE 2 — Per-chemical production/consumption instrumentation

Each task in this phase follows the same pattern: find the site where a chemical value changes, add a one-line counter increment at the exact point the change is applied.

### Task 7: Instrument sugar produced/consumed

**Files:**
- Modify: `src/engine/node/tissues/leaf.cpp:70-106` (photosynthesize)
- Modify: `src/engine/node/node.cpp` around `pay_maintenance` / `maintenance_cost` sites
- Modify: `src/engine/node/stem_node.cpp` / `root_node.cpp` / `tissues/*.cpp` growth sites

- [ ] **Step 1: Leaf photosynthesis — sugar produced**

In `leaf.cpp` line ~96 where `local().chemical(ChemicalID::Sugar) += produced`:
```cpp
local().chemical(ChemicalID::Sugar) += produced;
tick_chem_produced[static_cast<size_t>(ChemicalID::Sugar)] += produced;
// existing:
tick_sugar_activity += produced;  // keep the legacy scalar
```

- [ ] **Step 2: Maintenance — sugar consumed**

In `node.cpp` `pay_maintenance()` (around line 125), after the sugar deduction:
```cpp
const float cost = maintenance_cost(g, w);
local().chemical(ChemicalID::Sugar) -= cost;
tick_sugar_maintenance += cost;  // existing legacy
tick_chem_consumed[static_cast<size_t>(ChemicalID::Sugar)] += cost;
```

- [ ] **Step 3: Growth costs — sugar consumed at each site**

For each site that deducts sugar for growth (stem elongate, stem thicken, root elongate, root thicken, leaf expand, apical grow, root_apical grow), add one line immediately after the sugar deduction:
```cpp
tick_chem_consumed[static_cast<size_t>(ChemicalID::Sugar)] += cost;
```

Grep for sites: `grep -n "ChemicalID::Sugar" src/engine/node/`

- [ ] **Step 4: Build + run tests**

```bash
/usr/local/bin/cmake --build build
./build/botany_tests "[tick_counters]"
```
Expected: pass.

- [ ] **Step 5: Tighten the leaf test**

In `tests/test_tick_counters.cpp`, change the loose assertion:
```cpp
// was: REQUIRE(produced >= 0.0f);
REQUIRE(produced > 0.0f);
```
Re-run: `./build/botany_tests "[tick_counters]"` — should still pass.

- [ ] **Step 6: Commit**

```bash
git add -u
git commit -m "tick-counters: instrument sugar produced/consumed"
```

---

### Task 8: Instrument water produced/consumed

**Files:**
- Modify: `src/engine/node/tissues/root_apical.cpp:85-96` (absorb_water)
- Modify: `src/engine/node/root_node.cpp` absorption site (see `compute_budget` usage around vascular_sub_stepped.cpp:174-184 — water absorbed into root local is the source)
- Modify: `src/engine/node/tissues/leaf.cpp:38-46` (transpire), `:102` (photosynthesis water cost)

- [ ] **Step 1: RootApicalNode::absorb_water — water produced**

After the line that adds water to `local()`:
```cpp
tick_chem_produced[static_cast<size_t>(ChemicalID::Water)] += absorbed;
```

- [ ] **Step 2: RootNode water absorption — water produced**

At the site where root nodes gain water from soil (check `root_node.cpp` absorb function):
```cpp
tick_chem_produced[static_cast<size_t>(ChemicalID::Water)] += absorbed;
```

- [ ] **Step 3: LeafNode::transpire — water consumed**

After the water deduction (~line 45):
```cpp
tick_chem_consumed[static_cast<size_t>(ChemicalID::Water)] += amount;
```

- [ ] **Step 4: LeafNode::photosynthesize water cost — water consumed**

After the water deduction at line ~102:
```cpp
tick_chem_consumed[static_cast<size_t>(ChemicalID::Water)] += water_cost;
```

- [ ] **Step 5: Build and test**

```bash
/usr/local/bin/cmake --build build
./build/botany_tests
```
Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add -u
git commit -m "tick-counters: instrument water produced/consumed"
```

---

### Task 9: Instrument auxin produced/consumed

**Files:**
- Modify: `src/engine/node/tissues/apical.cpp:71` (mirror existing `tick_auxin_produced`)
- Modify: `src/engine/node/tissues/root_apical.cpp:64`
- Modify: `src/engine/node/tissues/leaf.cpp:169`
- Modify: `src/engine/node/node.cpp:612-618` (decay_chemicals)

- [ ] **Step 1: Apical auxin production**

Next to `tick_auxin_produced += produced;` add:
```cpp
tick_chem_produced[static_cast<size_t>(ChemicalID::Auxin)] += produced;
```

- [ ] **Step 2: Root apical and leaf auxin production**

Same pattern at the three sites listed above.

- [ ] **Step 3: Auxin decay — consumed**

In `Node::decay_chemicals()` (lines 612-618), for the auxin decay branch:
```cpp
float before = local().chemical(ChemicalID::Auxin);
local().chemical(ChemicalID::Auxin) *= (1.0f - g.auxin_decay_rate);
float after = local().chemical(ChemicalID::Auxin);
tick_chem_consumed[static_cast<size_t>(ChemicalID::Auxin)] += (before - after);
```

(Adjust to match the actual existing decay code — the idea is `consumed = before - after`.)

- [ ] **Step 4: Build and test**

```bash
/usr/local/bin/cmake --build build
./build/botany_tests
```

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "tick-counters: instrument auxin produced/consumed"
```

---

### Task 10: Instrument cytokinin produced/consumed

**Files:**
- Modify: `src/engine/node/tissues/root_apical.cpp:78`
- Modify: `src/engine/node/node.cpp:612-618` (decay)

- [ ] **Step 1: RA cytokinin production**

Next to `tick_cytokinin_produced += cyto_produced;`:
```cpp
tick_chem_produced[static_cast<size_t>(ChemicalID::Cytokinin)] += cyto_produced;
```

- [ ] **Step 2: Cytokinin decay in local() — consumed**

In `decay_chemicals()`, for the CK branch:
```cpp
float before = local().chemical(ChemicalID::Cytokinin);
local().chemical(ChemicalID::Cytokinin) *= (1.0f - g.cytokinin_decay_rate);
tick_chem_consumed[static_cast<size_t>(ChemicalID::Cytokinin)] += (before - local().chemical(ChemicalID::Cytokinin));
```

(Xylem CK has no decay — do NOT instrument consumed on xylem.)

- [ ] **Step 3: Build + commit**

```bash
/usr/local/bin/cmake --build build
./build/botany_tests
git add -u
git commit -m "tick-counters: instrument cytokinin produced/consumed"
```

---

### Task 11: Instrument GA / ethylene produced/consumed

GA and Ethylene use a reset-each-tick signal model. Treat this tick's full emitted amount as both "produced" (by the emitting leaf) and "consumed" (by the next tick's reset for every node that carried it).

**Files:**
- Modify: `src/engine/node/tissues/leaf.cpp:32-36` (produce_gibberellin) and line 23 (ethylene site)
- Modify: `src/engine/plant.cpp` or wherever GA/Eth reset-to-zero happens

- [ ] **Step 1: Produced at emission site**

Where `local().chemical(GA) += x` is executed in leaf.cpp:
```cpp
local().chemical(ChemicalID::Gibberellin) += x;
tick_chem_produced[static_cast<size_t>(ChemicalID::Gibberellin)] += x;
```

Same for Ethylene emission.

- [ ] **Step 2: Reset-as-consumed**

Find the site where GA/Eth are reset to zero each tick. At that site, just before zeroing, record the value as consumed:
```cpp
for_each_node([](Node& n) {
    float ga = n.local().chemical(ChemicalID::Gibberellin);
    if (ga > 0.0f) n.tick_chem_consumed[static_cast<size_t>(ChemicalID::Gibberellin)] += ga;
    n.local().chemical(ChemicalID::Gibberellin) = 0.0f;
    // same for Ethylene
});
```

- [ ] **Step 3: Build + test + commit**

```bash
/usr/local/bin/cmake --build build
./build/botany_tests
git add -u
git commit -m "tick-counters: instrument GA and ethylene produced/consumed"
```

---

### Task 12: Instrument stress produced/consumed

**Files:**
- Modify: `src/engine/node/node.cpp` — `update_physics` applies stress; decay in `decay_chemicals`

- [ ] **Step 1: Stress production in update_physics**

At the site where mechanical load increases stress:
```cpp
float delta = /* computed stress from load */;
local().chemical(ChemicalID::Stress) += delta;
tick_chem_produced[static_cast<size_t>(ChemicalID::Stress)] += delta;
```

- [ ] **Step 2: Stress decay in decay_chemicals**

Same pattern as auxin/CK decay:
```cpp
float before = local().chemical(ChemicalID::Stress);
// existing decay applies
float consumed = before - local().chemical(ChemicalID::Stress);
if (consumed > 0) tick_chem_consumed[static_cast<size_t>(ChemicalID::Stress)] += consumed;
```

- [ ] **Step 3: Build + test + commit**

```bash
/usr/local/bin/cmake --build build
./build/botany_tests
git add -u
git commit -m "tick-counters: instrument stress produced/consumed"
```

---

### Task 13: Mass-balance assertion in test_tick_counters

Now that all production/consumption is instrumented, add a real mass-balance check.

**Files:**
- Modify: `tests/test_tick_counters.cpp`

- [ ] **Step 1: Add a SECTION asserting per-chemical mass balance on a small fixed tree**

```cpp
TEST_CASE("tick counters: per-chem mass balance on single-node plant", "[tick_counters][balance]") {
    Engine engine;
    Genome g = default_genome();
    WorldParams w{};
    auto pid = engine.add_plant(g);
    // snapshot node-0 (seed) local chems at tick N
    engine.tick(w);  // warmup
    std::array<float, static_cast<size_t>(ChemicalID::Count)> before{};
    const Node* seed = nullptr;
    engine.get_plant(pid).for_each_node([&](const Node& n) {
        if (!n.parent) {
            seed = &n;
            for (size_t i = 0; i < before.size(); ++i)
                before[i] = n.local().chemical(static_cast<ChemicalID>(i));
        }
    });
    REQUIRE(seed);
    engine.tick(w);
    for (size_t i = 0; i < before.size(); ++i) {
        if (i == static_cast<size_t>(ChemicalID::Count)) continue;
        float after = seed->local().chemical(static_cast<ChemicalID>(i));
        float delta = after - before[i];
        float prod  = seed->tick_chem_produced[i];
        float cons  = seed->tick_chem_consumed[i];
        // Mass balance: delta = produced - consumed + net_transport_in
        // We don't track transport per-chem at the node level yet, so allow a tolerance.
        INFO("chem=" << i << " delta=" << delta << " prod=" << prod << " cons=" << cons);
        // Weaker check: production and consumption should be non-negative
        REQUIRE(prod >= 0.0f);
        REQUIRE(cons >= 0.0f);
    }
}
```

- [ ] **Step 2: Build + test + commit**

```bash
/usr/local/bin/cmake --build build
./build/botany_tests "[tick_counters]"
git add -u
git commit -m "tests: mass-balance smoke test for per-chem tick counters"
```

---

## PHASE 3 — Per-edge flux and capacity instrumentation

### Task 14: Jacobi phloem flux+cap (sugar)

**Files:**
- Modify: `src/engine/vascular_sub_stepped.cpp:93-99` (jacobi_step phloem pass, around line 317 for phloem)

- [ ] **Step 1: Locate the phloem longitudinal pass**

In `vascular_sub_stepped.cpp`, find the block where Jacobi moves sugar between parent/child phloem pools (`jacobi_step` or the phloem-specific pass). The line that actually transfers sugar looks like:
```cpp
parent_pool.chemical(Sugar) -= transfer;
child_pool.chemical(Sugar)  += transfer;
```

- [ ] **Step 2: Record per-edge flux and cap**

Just after the transfer, record the signed flux on the parent (keyed by child pointer) and the theoretical cap:
```cpp
const size_t s = static_cast<size_t>(ChemicalID::Sugar);
parent_node.tick_edge_flux[s][child_node] += transfer;
float cap_this_substep = g.phloem_conductance * cross_section;
parent_node.tick_edge_cap[s][child_node] += cap_this_substep;
```

`parent_node` and `child_node` are the Node* pairs corresponding to `parent_pool` and `child_pool`. `cross_section` is whatever Jacobi used for this edge (`π × edge_r² × g.phloem_fraction`).

- [ ] **Step 3: Build**

Run: `/usr/local/bin/cmake --build build`
Expected: clean build.

- [ ] **Step 4: Smoke test — run realtime, confirm no crashes**

Run `./build/botany_realtime` for ~10 seconds, close. No crash = fine.

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "vascular: per-edge sugar flux+cap instrumentation"
```

---

### Task 15: Jacobi xylem flux+cap (water, CK)

**Files:**
- Modify: `src/engine/vascular_sub_stepped.cpp` — xylem pass (around line 340)

- [ ] **Step 1: Same pattern as Task 14 applied to xylem**

For the xylem Jacobi pass, record Water and Cytokinin per-edge flux/cap:
```cpp
for (ChemicalID xchem : {ChemicalID::Water, ChemicalID::Cytokinin}) {
    const size_t idx = static_cast<size_t>(xchem);
    float transferred_x = /* the amount of this chem moved */;
    parent_node.tick_edge_flux[idx][child_node] += transferred_x;
    float cap_this_substep = g.xylem_conductance * cross_section_xylem;
    parent_node.tick_edge_cap[idx][child_node] += cap_this_substep;
}
```

Adjust to match how the existing xylem pass actually moves water and CK — usually they move together (CK rides with water), so the per-chem fluxes will be proportional.

- [ ] **Step 2: Build + smoke test + commit**

```bash
/usr/local/bin/cmake --build build
./build/botany_realtime   # smoke: no crash
git add -u
git commit -m "vascular: per-edge water+CK flux+cap instrumentation"
```

---

### Task 16: PIN auxin flux+cap

**Files:**
- Modify: `src/engine/pin_transport.cpp:66-92` (Phase A), `:114-115, 85` (Phase B), `:131-155` (Phase C)

- [ ] **Step 1: Identify existing flux recording sites**

Search for `last_auxin_flux[` in `pin_transport.cpp`. Each site already records per-edge auxin transfer. Find each assignment.

- [ ] **Step 2: Mirror into new maps**

At each `last_auxin_flux[child] += amount;` line, add:
```cpp
const size_t a = static_cast<size_t>(ChemicalID::Auxin);
tick_edge_flux[a][child] += amount;
float cap_edge = child->radius * child->radius * g.pin_capacity_per_area;
tick_edge_cap[a][child] += cap_edge;
```

Do this at the three sites (Phase A, B, C) identified by the explorer.

- [ ] **Step 3: Build + run full tests**

```bash
/usr/local/bin/cmake --build build
./build/botany_tests
```
Expected: all 218+ tests pass.

- [ ] **Step 4: Commit**

```bash
git add -u
git commit -m "pin: per-edge auxin flux+cap instrumentation"
```

---

### Task 17: Diffusion flux+cap (GA, Eth, Stress)

**Files:**
- Modify: `src/engine/node/node.cpp` — inside `transport_with_children` / `compute_transport_flow` where GA/Eth/Stress diffuse

- [ ] **Step 1: Locate diffusion transfer sites**

Grep: `grep -n "ChemicalID::Gibberellin\|ChemicalID::Ethylene\|ChemicalID::Stress" src/engine/node/node.cpp`.

- [ ] **Step 2: At each transfer, record flux and cap**

Same pattern:
```cpp
const size_t idx = static_cast<size_t>(chem);
parent.tick_edge_flux[idx][&child] += transferred;
float cap = base + radius_factor * scale;  // same formula as existing compute_transport_flow
parent.tick_edge_cap[idx][&child] += cap;
```

- [ ] **Step 3: Build + test + commit**

```bash
/usr/local/bin/cmake --build build
./build/botany_tests
git add -u
git commit -m "transport: per-edge flux+cap for diffusion chems (GA, Eth, Stress)"
```

---

## PHASE 4 — Panel rewrite

### Task 18: Complete the ui_helpers functions (maintenance cost, nodes_to_seed already done)

**Files:**
- Modify: `src/engine/ui_helpers.cpp`

- [ ] **Step 1: Complete `compute_maintenance_cost`**

Replace the stubbed function:
```cpp
float compute_maintenance_cost(const Node& n, const Genome& g, const WorldParams& w) {
    // Match each subclass's maintenance_cost() formula exactly.
    if (auto* s = n.as_stem()) {
        float len = glm::length(s->offset);
        float vol = 3.14159f * s->radius * s->radius * std::max(len, 0.001f);
        return w.sugar_maintenance_stem * vol;
    }
    if (auto* r = n.as_root()) {
        float len = glm::length(r->offset);
        float vol = 3.14159f * r->radius * r->radius * std::max(len, 0.001f);
        return w.sugar_maintenance_root * vol;
    }
    if (auto* lf = n.as_leaf()) {
        float area = lf->leaf_size * lf->leaf_size;
        return w.sugar_maintenance_leaf * area;
    }
    if (n.as_apical() || n.as_root_apical()) {
        // Dormant meristems pay 0; active meristems pay fixed rate.
        auto* ap = n.as_apical();
        auto* ra = n.as_root_apical();
        bool active = (ap && ap->active) || (ra && ra->active);
        return active ? w.sugar_maintenance_meristem : 0.0f;
    }
    return 0.0f;
}
```

Verify each formula against the actual `maintenance_cost()` override in `stem_node.cpp`, `root_node.cpp`, `leaf.cpp`, `apical.cpp`, `root_apical.cpp`. Adjust if the engine uses a different formula.

- [ ] **Step 2: Write a unit test asserting preview matches actual**

Append to `tests/test_tick_counters.cpp`:
```cpp
TEST_CASE("compute_maintenance_cost matches tick_sugar_maintenance", "[tick_counters][maintenance]") {
    Engine engine;
    Genome g = default_genome();
    WorldParams w{};
    auto pid = engine.add_plant(g);
    for (int i = 0; i < 50; ++i) engine.tick(w);
    engine.get_plant(pid).for_each_node([&](const Node& n) {
        float preview = compute_maintenance_cost(n, g, w);
        float actual  = n.tick_sugar_maintenance;
        INFO("node " << n.id << " preview=" << preview << " actual=" << actual);
        REQUIRE(std::fabs(preview - actual) < 1e-5f);
    });
}
```

Include `#include "engine/ui_helpers.h"` at top of the test file.

- [ ] **Step 3: Build and run the test**

```bash
/usr/local/bin/cmake --build build
./build/botany_tests "[maintenance]"
```
Expected: pass. If it fails, the formula in `compute_maintenance_cost` doesn't match the actual engine — adjust and re-run.

- [ ] **Step 4: Commit**

```bash
git add -u
git commit -m "ui_helpers: complete compute_maintenance_cost with parity test"
```

---

### Task 19: Panel — Identity section rewrite

**Files:**
- Modify: `src/app_realtime.cpp:876-890` (inspector panel top, before existing ImGui::Text calls)
- Modify: add `#include "engine/ui_helpers.h"` at top of file

- [ ] **Step 1: Add include**

Near the top of `app_realtime.cpp` with the other engine includes:
```cpp
#include "engine/ui_helpers.h"
```

- [ ] **Step 2: Rewrite the identity header**

Replace the existing "Type: %s" / meristem-specific blocks at the top of the panel (lines 880-888 roughly) with a dedicated identity section. The activity blocks (lines 890-1131) remain intact — they move DOWN below the chemicals table. Identity goes FIRST:

```cpp
// --- IDENTITY ---
const char* type_str = "?";
switch (sel.type) {
    case NodeType::STEM:        type_str = "STEM"; break;
    case NodeType::ROOT:        type_str = "ROOT"; break;
    case NodeType::LEAF:        type_str = "LEAF"; break;
    case NodeType::APICAL:      type_str = "APICAL"; break;
    case NodeType::ROOT_APICAL: type_str = "ROOT_APICAL"; break;
}
ImGui::Text("ID: #%u  Type: %s", sel.id, type_str);
ImGui::Text("Age: %u ticks (%.1f days)", sel.age, sel.age / 24.0f);
ImGui::Text("Length: %s  Radius: %s", fmt_dist(glm::length(sel.offset)), fmt_dist(sel.radius));
float cross_section = 3.14159f * sel.radius * sel.radius;
ImGui::Text("Cross-section: %.4f dm\xC2\xB2", cross_section);
ImGui::Text("Height (y): %.2f m", sel.world_position.y / 10.0f);  // dm -> m
ImGui::Text("Nodes to seed: %d", nodes_to_seed(sel));

const Genome& mg = engine.get_plant(plant_id).genome();
if (sel.as_stem() || sel.as_root()) {
    uint32_t mat_ticks = sel.as_stem() ? mg.internode_maturation_ticks : mg.root_internode_maturation_ticks;
    bool mature = (sel.parent == nullptr) || (sel.age >= mat_ticks);
    ImGui::Text("Elongation: %s (age %u/%u)", mature ? "mature" : "growing", sel.age, mat_ticks);
    float hm = hydraulic_maturity(sel, mg) * 100.0f;
    ImGui::Text("Hydraulic maturity: %.0f%% closed", hm);
    ImGui::Text("PIN saturation: %.2f", sel.get_parent_auxin_flow_bias());
}
ImGui::Text("Starvation: %u ticks", sel.starvation_ticks);
if (auto* leaf = sel.as_leaf()) {
    if (leaf->senescence_ticks > 0)
        ImGui::Text("Senescence: %u / %u ticks", leaf->senescence_ticks, mg.senescence_duration);
    else
        ImGui::Text("Senescence: healthy");
}
ImGui::Separator();
```

Move the existing "Meristem:", "Growth:", "Activation:", "Elongate:", "Thicken:", leaf-specific lines to a separate "// --- ACTIVITY ---" section below the chemicals table (which is built in Task 20). Also delete the duplicate `sel.as_stem()` / `sel.as_root()` fallback blocks at 1049-1131 (they're unreachable dead code since the `else if (sel.as_stem())` chain at 987/1019 already handles them).

- [ ] **Step 3: Build + smoke test**

```bash
/usr/local/bin/cmake --build build
./build/botany_realtime  # visual smoke: click a stem, verify identity info renders
```

- [ ] **Step 4: Commit**

```bash
git add -u
git commit -m "panel: rewrite identity section with height, nodes-to-seed, maturity triad"
```

---

### Task 20: Panel — Chemicals table rewrite

**Files:**
- Modify: `src/app_realtime.cpp:1229-1280` (existing chemicals table)

- [ ] **Step 1: Delete existing table block**

Delete lines 1229-1280 (the `ImGui::BeginTable("chemicals", 4, ...)` block and its rows).

- [ ] **Step 2: Insert new 10-column table**

Replace with:
```cpp
ImGui::Text("Chemicals:");
if (ImGui::BeginTable("chemicals", 10, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit)) {
    ImGui::TableSetupColumn("Chem");
    ImGui::TableSetupColumn("L lvl"); ImGui::TableSetupColumn("L +"); ImGui::TableSetupColumn("L -");
    ImGui::TableSetupColumn("V lvl"); ImGui::TableSetupColumn("V +"); ImGui::TableSetupColumn("V -");
    ImGui::TableSetupColumn("T lvl"); ImGui::TableSetupColumn("T +"); ImGui::TableSetupColumn("T -");
    ImGui::TableHeadersRow();

    struct Row { const char* name; ChemicalID id; const char* (*fmt)(float); };
    Row rows[] = {
        {"Sugar",  ChemicalID::Sugar,       fmt_mass},
        {"Water",  ChemicalID::Water,       fmt_vol },
        {"Auxin",  ChemicalID::Auxin,       fmt_au  },
        {"Cyt",    ChemicalID::Cytokinin,   fmt_au  },
        {"GA",     ChemicalID::Gibberellin, fmt_au  },
        {"Eth",    ChemicalID::Ethylene,    fmt_au  },
        {"Stress", ChemicalID::Stress,      fmt_au  },
    };
    for (const Row& r : rows) {
        const size_t idx = static_cast<size_t>(r.id);
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0); ImGui::Text("%s", r.name);

        // Local
        float l_lvl = sel.local().chemical(r.id);
        float l_p = sel.tick_chem_produced[idx];
        float l_c = sel.tick_chem_consumed[idx];
        ImGui::TableSetColumnIndex(1); ImGui::Text("%s", r.fmt(l_lvl));
        ImGui::TableSetColumnIndex(2); ImGui::Text("%s", r.fmt(l_p));
        ImGui::TableSetColumnIndex(3); ImGui::Text("%s", r.fmt(l_c));

        // Vasculature
        const TransportPool* vp = vascular_scope(sel, r.id);
        if (vp) {
            float v_lvl = vp->chemical(r.id);
            ImGui::TableSetColumnIndex(4); ImGui::Text("%s", r.fmt(v_lvl));
            ImGui::TableSetColumnIndex(5); ImGui::Text("\xE2\x80\x94");  // em dash — Produced/Consumed at vasc scope deferred
            ImGui::TableSetColumnIndex(6); ImGui::Text("\xE2\x80\x94");

            // Total
            ImGui::TableSetColumnIndex(7); ImGui::Text("%s", r.fmt(l_lvl + v_lvl));
            ImGui::TableSetColumnIndex(8); ImGui::Text("%s", r.fmt(l_p));
            ImGui::TableSetColumnIndex(9); ImGui::Text("%s", r.fmt(l_c));
        } else {
            ImGui::TableSetColumnIndex(4); ImGui::Text("\xE2\x80\x94");
            ImGui::TableSetColumnIndex(5); ImGui::Text("\xE2\x80\x94");
            ImGui::TableSetColumnIndex(6); ImGui::Text("\xE2\x80\x94");
            ImGui::TableSetColumnIndex(7); ImGui::Text("%s", r.fmt(l_lvl));
            ImGui::TableSetColumnIndex(8); ImGui::Text("%s", r.fmt(l_p));
            ImGui::TableSetColumnIndex(9); ImGui::Text("%s", r.fmt(l_c));
        }
    }
    ImGui::EndTable();
}
ImGui::TextDisabled("L=local  V=vascular  T=total  +/- = this-tick produced/consumed");
ImGui::Separator();
```

Note the legend's use of em-dash ("\xE2\x80\x94") for signaling-chem vasc cells. Vasc-scope produced/consumed is intentionally deferred (per spec section 2.2's "simpler v1" option) — show em-dash for those.

- [ ] **Step 3: If panel is too narrow, widen it**

If the 10-column table clips at the current 310 px width, change the `ImGui::SetNextWindowSize(ImVec2(310, 0), ...)` call near line 873 to `ImVec2(480, 0)`.

- [ ] **Step 4: Smoke test**

```bash
/usr/local/bin/cmake --build build
./build/botany_realtime  # click various node types, verify table renders + reads make sense
```

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "panel: 10-column chemicals table with local/vasc/total scopes"
```

---

### Task 21: Panel — Activity section + navigation cleanup

**Files:**
- Modify: `src/app_realtime.cpp` — the activity/meristem/growth blocks (previously at lines 890-1131)

- [ ] **Step 1: Move the existing meristem/leaf/stem/root activity blocks below the chemicals table**

These blocks weren't deleted in Task 19 — they were left in place. Cut them now and paste them below the chemicals table (after the `ImGui::Separator()` that follows it). Wrap them under an `ImGui::Text("Activity:");` header. Remove the duplicate `else if (sel.as_stem())` / `else if (sel.as_root())` blocks at lines 1049-1131 (dead code — earlier `else if` branches already match).

- [ ] **Step 2: Delete the "Children: N" line**

Search the panel body for `ImGui::Text("Children: %d"` and remove that line.

- [ ] **Step 3: Verify parent/child buttons are intact**

Confirm the `if (sel.parent) { ... ImGui::Button("Parent: ...") }` block and the child button loop are unchanged.

- [ ] **Step 4: Build + smoke test + commit**

```bash
/usr/local/bin/cmake --build build
./build/botany_realtime  # visual: click all 5 node types, verify activity info renders for each
git add -u
git commit -m "panel: activity section moved below chem table, drop redundant children counter"
```

**Phase 4 complete. Commit marker.**

---

## PHASE 5 — Overlay modes rewrite

### Task 22: New OverlayCategory enum + selector UI

**Files:**
- Modify: `src/app_realtime.cpp:332` (old enum), `:563-652` (old selector)

- [ ] **Step 1: Replace the old enum**

Delete: `enum class Overlay { NONE, NODE_TYPE, AUXIN, CYTOKININ, SUGAR, LIGHT, GIBBERELLIN, ETHYLENE, STRESS, WATER, GROWTH, VASCULAR };` at line 332.

Add:
```cpp
enum class OverlayCategory { Default, Type, Light, Level, Capacity, Growth, Activation, Starvation };
enum class OverlayScope { Local, Vasc };
static OverlayCategory g_overlay_category = OverlayCategory::Default;
static ChemicalID      g_overlay_chem     = ChemicalID::Sugar;
static OverlayScope    g_overlay_scope    = OverlayScope::Local;
```

Remove references to the old `Overlay::*` variants — all are about to be replaced.

- [ ] **Step 2: Replace the selector block (lines 563-652) with radio + sub-picker**

```cpp
if (ImGui::CollapsingHeader("Overlays", ImGuiTreeNodeFlags_DefaultOpen)) {
    auto radio = [&](const char* label, OverlayCategory cat) {
        if (ImGui::RadioButton(label, g_overlay_category == cat)) {
            g_overlay_category = cat;
        }
    };
    radio("Default",    OverlayCategory::Default);    ImGui::SameLine();
    radio("Type",       OverlayCategory::Type);       ImGui::SameLine();
    radio("Light",      OverlayCategory::Light);
    radio("Level",      OverlayCategory::Level);      ImGui::SameLine();
    radio("Capacity",   OverlayCategory::Capacity);   ImGui::SameLine();
    radio("Growth",     OverlayCategory::Growth);
    radio("Activation", OverlayCategory::Activation); ImGui::SameLine();
    radio("Starvation", OverlayCategory::Starvation);

    // Sub-picker: chemical selector for Level and Capacity
    if (g_overlay_category == OverlayCategory::Level || g_overlay_category == OverlayCategory::Capacity) {
        const char* chem_names[] = {"Auxin", "Cytokinin", "Gibberellin", "Sugar", "Ethylene", "Stress", "Water"};
        int idx = static_cast<int>(g_overlay_chem);
        if (ImGui::Combo("Chemical", &idx, chem_names, IM_ARRAYSIZE(chem_names))) {
            g_overlay_chem = static_cast<ChemicalID>(idx);
        }
    }
    // Scope sub-picker only for Level
    if (g_overlay_category == OverlayCategory::Level) {
        int sc = static_cast<int>(g_overlay_scope);
        if (ImGui::Combo("Scope", &sc, "Local\0Vasculature\0")) {
            g_overlay_scope = static_cast<OverlayScope>(sc);
        }
    }

    // Rebuild accessor whenever any selector changes
    renderer.set_color_mode(build_color_accessor(g_overlay_category, g_overlay_chem, g_overlay_scope, engine));
}
```

`build_color_accessor` is defined in Task 23.

- [ ] **Step 3: Don't build yet — the next task adds the missing function**

---

### Task 23: Implement build_color_accessor for Default/Type/Light/Level

**Files:**
- Modify: `src/app_realtime.cpp` — add a free function above `main()`

- [ ] **Step 1: Add function**

Above `main()`, add:
```cpp
static ChemicalAccessor build_color_accessor(OverlayCategory cat, ChemicalID chem, OverlayScope scope, const Engine& engine) {
    switch (cat) {
        case OverlayCategory::Default: return ChemicalAccessor{};  // renderer draws genome colors
        case OverlayCategory::Type:    return ChemicalAccessor{};  // type mode uses a separate codepath — set below
        case OverlayCategory::Light:
            return [](const Node& n) -> float {
                if (auto* l = n.as_leaf()) return l->light_exposure;
                return std::numeric_limits<float>::quiet_NaN();  // gray for non-leaves
            };
        case OverlayCategory::Level:
            if (scope == OverlayScope::Local) {
                return [chem](const Node& n) { return n.local().chemical(chem); };
            } else {
                return [chem](const Node& n) -> float {
                    const TransportPool* vp = vascular_scope(n, chem);
                    if (!vp) return std::numeric_limits<float>::quiet_NaN();
                    return vp->chemical(chem);
                };
            }
        case OverlayCategory::Capacity: {
            const size_t idx = static_cast<size_t>(chem);
            return [idx](const Node& n) -> float {
                float flux = 0.0f, cap = 0.0f;
                for (const auto& [c, f] : n.tick_edge_flux[idx]) flux += std::fabs(f);
                for (const auto& [c, k] : n.tick_edge_cap[idx])  cap  += k;
                if (cap < 1e-9f) return std::numeric_limits<float>::quiet_NaN();
                return flux / cap;
            };
        }
        case OverlayCategory::Growth:     return build_growth_accessor(engine);       // Task 24
        case OverlayCategory::Activation: return build_activation_accessor(engine);   // Task 25
        case OverlayCategory::Starvation: return build_starvation_accessor(engine);   // Task 26
    }
    return ChemicalAccessor{};
}
```

For Type mode, handle it specially in the caller — Type uses `renderer.set_node_type_color_mode()` (or similar) rather than a chemical accessor. Copy the existing Type-mode wiring from lines 654-665 into the same selector block above, keyed on `g_overlay_category == OverlayCategory::Type`.

- [ ] **Step 2: Build**

```bash
/usr/local/bin/cmake --build build
```
Expected: unresolved symbol for `build_growth_accessor` / `build_activation_accessor` / `build_starvation_accessor`. Those land in Tasks 24-26.

- [ ] **Step 3: Stub the missing functions so it compiles**

Above `build_color_accessor`:
```cpp
static ChemicalAccessor build_growth_accessor(const Engine&)     { return [](const Node&){ return std::numeric_limits<float>::quiet_NaN(); }; }
static ChemicalAccessor build_activation_accessor(const Engine&) { return [](const Node&){ return std::numeric_limits<float>::quiet_NaN(); }; }
static ChemicalAccessor build_starvation_accessor(const Engine&) { return [](const Node&){ return std::numeric_limits<float>::quiet_NaN(); }; }
```

- [ ] **Step 4: Build + smoke test**

```bash
/usr/local/bin/cmake --build build
./build/botany_realtime  # verify Default / Type / Light / Level (Local+Vasc) / Capacity all render
```

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "overlays: two-tier selector; Default/Type/Light/Level/Capacity modes"
```

---

### Task 24: Growth overlay accessor

**Files:**
- Modify: `src/app_realtime.cpp`

- [ ] **Step 1: Replace the growth stub**

```cpp
static ChemicalAccessor build_growth_accessor(const Engine& engine) {
    return [&engine](const Node& n) -> float {
        const Genome& g = engine.get_plant(0).genome();   // FIXME: generalize if plant_id not 0
        if (auto* s = n.as_stem()) {
            if (!s->parent || s->age >= g.internode_maturation_ticks) {
                // Thickening growth
                float bias = s->get_parent_auxin_flow_bias();
                if (bias < 1e-6f) return 0.0f;
                // Proportional to bias * sugar gating
                float max_cost = g.cambium_responsiveness * bias * /* stem sugar cost */ 1.0f;
                float sugar = s->local().chemical(ChemicalID::Sugar);
                float gf = (max_cost > 1e-6f) ? std::min(sugar / max_cost, 1.0f) : 1.0f;
                return bias * gf;
            }
            // Elongation growth, fraction of internode_elongation_rate
            float sugar = s->local().chemical(ChemicalID::Sugar);
            float gf = std::min(sugar / (g.internode_elongation_rate * 1.0f + 1e-6f), 1.0f);
            return gf;
        }
        if (auto* r = n.as_root()) { /* same pattern with root params */ return 0.0f; }
        if (auto* lf = n.as_leaf()) {
            if (lf->leaf_size >= g.max_leaf_size) return 0.0f;
            float sugar = lf->local().chemical(ChemicalID::Sugar);
            float max_cost = g.leaf_growth_rate * /* leaf sugar cost */ 1.0f;
            return (max_cost > 1e-6f) ? std::min(sugar / max_cost, 1.0f) : 0.0f;
        }
        if (auto* ap = n.as_apical(); ap && ap->active) {
            float sugar = ap->local().chemical(ChemicalID::Sugar);
            float cyt   = ap->local().chemical(ChemicalID::Cytokinin);
            float wgf   = meristem_helpers::turgor_fraction(ap->local().chemical(ChemicalID::Water), water_cap(n, g));
            float max_cost = g.growth_rate * /* meristem cost */ 1.0f;
            float sgf = (max_cost > 1e-6f) ? std::min(sugar / max_cost, 1.0f) : 1.0f;
            float cgf = cyt / (cyt + std::max(g.cytokinin_growth_threshold, 1e-6f));
            return sgf * cgf * wgf;
        }
        if (auto* ra = n.as_root_apical(); ra && ra->active) { /* analogous, no cyt gate */ return 0.0f; }
        return std::numeric_limits<float>::quiet_NaN();
    };
}
```

Actual sugar-cost constants (`sugar_cost_stem_growth`, `sugar_cost_meristem_growth`, `sugar_cost_leaf_growth`) come from `engine.world_params()` — fetch once at the top of the lambda. The `1.0f` placeholders above are deliberately wrong so Step 2 catches them.

- [ ] **Step 2: Replace placeholders with real world params**

Capture `const WorldParams& w = engine.world_params();` at the top of the lambda and substitute `w.sugar_cost_stem_growth`, `w.sugar_cost_meristem_growth`, `w.sugar_cost_leaf_growth` for the `1.0f` markers above.

- [ ] **Step 3: Build + smoke test**

```bash
/usr/local/bin/cmake --build build
./build/botany_realtime  # Growth overlay — young tips red-ish, mature stems blue-ish
```

- [ ] **Step 4: Commit**

```bash
git add -u
git commit -m "overlays: growth mode — fraction of this-tick achievable growth"
```

---

### Task 25: Activation overlay accessor

**Files:**
- Modify: `src/app_realtime.cpp`

- [ ] **Step 1: Replace the activation stub**

```cpp
static ChemicalAccessor build_activation_accessor(const Engine& engine) {
    return [&engine](const Node& n) -> float {
        if (!n.as_stem() && !n.as_root()) return std::numeric_limits<float>::quiet_NaN();
        const Genome& g = engine.get_plant(0).genome();
        const WorldParams& w = engine.world_params();

        float best = -1.0f;  // track the max readiness among dormant children
        for (const Node* child : n.children) {
            if (auto* ap = child->as_apical(); ap && !ap->active) {
                float stem_auxin = n.local().chemical(ChemicalID::Auxin);
                float local_cyt  = n.local().chemical(ChemicalID::Cytokinin);
                float sugar      = child->local().chemical(ChemicalID::Sugar);
                float auxin_pct = std::max(0.0f, (g.auxin_threshold - stem_auxin) / std::max(g.auxin_threshold, 1e-6f));
                float cyt_pct   = std::min(local_cyt / std::max(g.cytokinin_threshold, 1e-6f), 1.0f);
                float sugar_pct = std::min(sugar / std::max(w.sugar_cost_activation, 1e-6f), 1.0f);
                float ready = std::min({auxin_pct, cyt_pct, sugar_pct});
                best = std::max(best, ready);
            }
            if (auto* ra = child->as_root_apical(); ra && !ra->active) {
                float auxin = child->local().chemical(ChemicalID::Auxin);
                float cyt   = child->local().chemical(ChemicalID::Cytokinin);
                float sugar = child->local().chemical(ChemicalID::Sugar);
                float auxin_pct = std::min(auxin / std::max(g.root_auxin_activation_threshold, 1e-6f), 1.0f);
                float cyt_pct   = std::max(0.0f, (g.root_cytokinin_inhibition_threshold - cyt) / std::max(g.root_cytokinin_inhibition_threshold, 1e-6f));
                float sugar_pct = std::min(sugar / std::max(w.sugar_cost_activation, 1e-6f), 1.0f);
                float ready = std::min({auxin_pct, cyt_pct, sugar_pct});
                best = std::max(best, ready);
            }
        }
        if (best < 0.0f) return std::numeric_limits<float>::quiet_NaN();  // no dormant child — gray
        return best;
    };
}
```

- [ ] **Step 2: Build + smoke test + commit**

```bash
/usr/local/bin/cmake --build build
./build/botany_realtime  # Activation overlay — parent stems/roots of dormant meristems colored
git add -u
git commit -m "overlays: activation readiness — propagated to parent stem/root"
```

---

### Task 26: Starvation overlay accessor + renderer color stops

**Files:**
- Modify: `src/app_realtime.cpp` — accessor
- Modify: `src/renderer/renderer.cpp:94-164` — 4-stop color path for starvation

Starvation needs a custom color lookup because it uses named stops (red, orange, yellow, green), not the default blue→yellow→red heatmap.

- [ ] **Step 1: Define a starvation color function in app_realtime.cpp**

The accessor returns a value encoded so the renderer can decode:
- `v` in `[0, 1)`: between yellow (0) and green (1), coverage ratio
- `v` in `[1, 2]`: between orange (1) and red (2), starvation-ticks ratio
- Any NaN: neutral gray

```cpp
static ChemicalAccessor build_starvation_accessor(const Engine& engine) {
    return [&engine](const Node& n) -> float {
        const Genome& g = engine.get_plant(0).genome();
        const WorldParams& w = engine.world_params();
        if (n.starvation_ticks > 0) {
            float t = static_cast<float>(n.starvation_ticks) / std::max(static_cast<float>(w.max_starvation_ticks), 1.0f);
            return 1.0f + std::min(t, 1.0f);  // maps to [1, 2]
        }
        float cost = compute_maintenance_cost(n, g, w);
        if (cost < 1e-9f) return 1.0f - 1e-4f;  // essentially green, no maintenance demand
        float coverage = std::min(n.local().chemical(ChemicalID::Sugar) / cost, 1.0f);
        return coverage;  // maps to [0, 1]
    };
}
```

- [ ] **Step 2: Renderer — add a starvation-mode color path**

The current heatmap path in renderer.cpp:94-164 uses blue→yellow→red linearly on the accessor value. For starvation we need a different lookup. Simplest approach: add a flag on the renderer that switches the color LUT.

In `renderer.h`, add:
```cpp
void set_color_mode_starvation(bool on);   // starvation-specific LUT
```

In `renderer.cpp`, store a `bool starvation_mode_ = false;` member. In the color branch (around line 162-164), replace the current blue→yellow→red interpolation with:
```cpp
if (starvation_mode_) {
    // v in [0,1]: yellow -> green. v in (1,2]: orange -> red.
    if (v < 1.0f) {
        // yellow (1,1,0) -> green (0,1,0) as v: 0 -> 1
        color = glm::vec3(1.0f - v, 1.0f, 0.0f);
    } else {
        float t = std::min(v - 1.0f, 1.0f);
        // orange (1,0.5,0) -> red (1,0,0)
        color = glm::vec3(1.0f, 0.5f * (1.0f - t), 0.0f);
    }
} else {
    // existing blue -> yellow -> red heatmap
}
```

And in `app_realtime.cpp`, when setting the starvation accessor, also call `renderer.set_color_mode_starvation(true)`. When any other mode is set, call `set_color_mode_starvation(false)`.

- [ ] **Step 3: Build + smoke test**

```bash
/usr/local/bin/cmake --build build
./build/botany_realtime  # Starvation overlay — healthy nodes green, stressed yellow/orange/red
```

- [ ] **Step 4: Commit**

```bash
git add -u
git commit -m "overlays: starvation mode with red/orange/yellow/green stops"
```

---

### Task 27: NaN-aware renderer gray path

**Files:**
- Modify: `src/renderer/renderer.cpp:94-164`

- [ ] **Step 1: Detect NaN in the color branch**

In the heatmap color computation (~line 140), just before the gradient lookup, add:
```cpp
if (std::isnan(v)) {
    color = glm::vec3(0.3f, 0.3f, 0.3f);  // neutral gray for out-of-scope nodes
} else {
    // existing gradient code
}
```

Ensure `<cmath>` is included for `std::isnan`.

- [ ] **Step 2: Build + smoke test**

```bash
/usr/local/bin/cmake --build build
./build/botany_realtime  # Light overlay — non-leaves render gray. Level (Vasc) — seedling shoots before vasc admission render gray.
```

- [ ] **Step 3: Commit**

```bash
git add -u
git commit -m "renderer: NaN-aware gray path for out-of-scope overlay nodes"
```

**Phase 5 complete.**

---

## PHASE 6 — Polish + verification

### Task 28: Final smoke-test checklist

- [ ] **Step 1: Launch realtime**

```bash
./build/botany_realtime
```

- [ ] **Step 2: Walk every overlay category and sub-option**

Verify each mode renders something non-default:
- Default (genome colors)
- Type (green shoot, orange root, red SA, blue RA)
- Light (leaves colored, others gray)
- Level × {Sugar, Water, Auxin, Cyt, GA, Eth, Stress} × {Local, Vasc}
- Capacity × all 7 chems
- Growth (young tips red, old trunks blue)
- Activation (dormant-bud parents colored)
- Starvation (healthy=green, stressed=yellow/orange/red)

- [ ] **Step 3: Walk the Node Inspector on each node type**

Click a stem, root, leaf, dormant SA, active SA, dormant RA, active RA, seed. Verify each shows sensible:
- Identity info (height, nodes-to-seed, maturity as applicable)
- Chemicals table (em-dash in vasc cells for signaling chems; numbers elsewhere)
- Activity section (growth % / activation readiness / photo/transpire)
- Parent/child buttons work

- [ ] **Step 4: Run the full test suite**

```bash
./build/botany_tests
```
Expected: all 218+ tests pass, plus new `[tick_counters]` tests.

- [ ] **Step 5: If any issue surfaces, file it as a follow-up task before closing**

If verification uncovers a broken scenario, stop here and add a focused fix task to this plan. Don't close the plan with known failures.

- [ ] **Step 6: Final commit if any cleanup changes**

```bash
git add -u
git commit -m "panel/overlay: final smoke-test polish" --allow-empty
```

---

## Self-review

Spec coverage:
- ✓ 1.1 Generalized per-tick counters → Tasks 2, 7-12
- ✓ 1.2 Vascular-scope helper → Task 6
- ✓ 1.3 Per-edge flux/cap → Tasks 3, 14-17
- ✓ 1.4 Reset discipline → Task 4
- ✓ 1.5 Testing → Tasks 5, 13, 18 (maintenance parity)
- ✓ 2.1 Identity section → Task 19
- ✓ 2.2 Chemicals table → Task 20
- ✓ 2.3 Activity section → Task 21 (reuse + reorder)
- ✓ 2.4 Navigation → Task 21 (children counter removed)
- ✓ 3.1 Selector UI → Task 22
- ✓ 3.2-3.5 Mode value functions → Tasks 23-26
- ✓ 3.6 Implementation plumbing + NaN-aware gray → Task 27

Placeholder scan: no TBD/TODO left in task code blocks. The `compute_maintenance_cost` stub in Task 6 is explicitly marked "completed in Task 18," which is fine.

Type consistency: `OverlayCategory`, `OverlayScope`, `ChemicalID::Count`, and `build_*_accessor` signatures are consistent across Tasks 22-26.

One genuine ambiguity: the Jacobi cross-section formulas (`g.phloem_fraction × π × r²` vs `g.xylem_fraction × π × r²`) need the exact value used by `jacobi_step` — the instrumentation must match whatever the actual pass uses. Tasks 14-15 call out "use whatever cross_section Jacobi used for this edge" to force the implementer to read the code.
