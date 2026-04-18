# Münch Pressure Flow: Design Document

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the current global-allocation phloem model with a local, pressure-driven Münch flow system where sugar transport emerges from per-edge osmotic pressure differentials rather than from a central scheduler.

**Architecture:** Each vascular node holds a phloem pressure derived from its local sugar concentration and water availability. Each tick, every parent-child edge computes flow proportional to the pressure difference and pipe conductance. Sugar dissolved in the flowing water moves with the stream. No global demand/supply aggregation is required.

**Scope:** Phloem only. The xylem pass (water + cytokinin) keeps its Phase 1/Phase 2 structure — it is a different physical system (transpiration pull) and is addressed in the world-physics milestone.

---

## 1. The Problem with the Current Allocation Model

The current `vascular_transport` is a two-phase global scheduler:

- **Phase 1 (post-order, leaves → seed):** Every node aggregates its subtree's supply and demand into the seed node. After Phase 1, `flat[0].supply` and `flat[0].demand` represent the *entire plant's* sugar surplus and deficit.
- **Phase 2 (pre-order, seed → leaves):** The seed distributes `min(supply, demand)` outward, proportionally by conductance weight, iterating to water-fill capped children.

This is a central planner model. It has the following structural problems:

### 1.1 Sugar teleports across the plant in one tick

A leaf at the far canopy tip and a meristem at the deep root tip communicate in a single tick because both contribute to the seed's global totals in Phase 1 and both draw from the seed's global pool in Phase 2. The seed mediates all exchange. In real plants, a sugar molecule moves from leaf mesophyll → sieve tube → phloem → adjacent sieve tube, one cell at a time. A large plant has tens of thousands of cells between source and sink. The Münch model moves sugar at most one edge per tick, naturally matching the physical constraint.

### 1.2 Meristems need artificial demand caps

`meristem_sink_fraction` exists because active meristems can declare a deficit up to `sugar_cap_meristem`, and the Phase 2 pass will try to fill it. Without the cap, a single hungry meristem outbids all other sinks at the seed junction and drains available supply. This cap is an engineering fix for a model artifact. In Münch flow, a meristem's pull is naturally bounded: it can only receive what its local pressure gradient and pipe conductance allow in one tick. No explicit cap is needed.

### 1.3 The leaf reserve floor is a static hack

`phloem_reserve_fraction` prevents leaves from contributing all their sugar to Phase 1 supply. This is necessary in the allocation model because leaves have no way to "hold back" what they need — they declare a surplus and the planner takes it. In Münch flow, a leaf holds its sugar as osmotic pressure. It loads into the phloem stream only at the loading rate. If a leaf's concentration is near its reserve floor, the loading rate drops naturally. No explicit floor is needed.

### 1.4 Flow direction is fixed as source → sink

The current model classifies every node as either source, sink, or conduit once per tick. Flow always goes from classified sources to classified sinks. In real phloem, flow direction is determined by pressure — it is fully bidirectional. A stem node with accumulated starch that mobilizes to sugar becomes a local source, pushing outward to adjacent sinks, regardless of whether it was classified a "conduit" during Phase 1. Münch flow handles this automatically: if `pressure_parent < pressure_child`, flow goes child → parent.

### 1.5 Water availability doesn't affect phloem flow

Münch flow is fundamentally osmotic: sieve tubes need xylem water to build turgor pressure at sources. If a leaf is drought-stressed, its phloem pressure should drop even if it has abundant sugar, because there's no water entering the sieve tubes. The current model has no such coupling — sugar transport ignores water status entirely (water only affects photosynthesis via stomatal conductance). The Münch model couples phloem pressure to water availability at each node, completing the hydraulic feedback loop.

---

## 2. How Münch Pressure Flow Works (Biology)

### 2.1 The physical mechanism

Phloem sieve tubes run from leaves to sinks as continuous pipes. At each source leaf, companion cells actively pump sucrose from the mesophyll into the sieve tube using ATP-driven symporters. This raises sucrose concentration inside the sieve tube well above the surrounding apoplast. The high solute concentration draws water in from adjacent xylem by osmosis. Water entry raises turgor pressure inside the sieve tube — typically 1–3 MPa at a productive source leaf, compared to 0.1–0.5 MPa at a consuming sink.

At sinks (meristems, growing roots, storage organs), sucrose is consumed or converted to starch. Sieve tube sucrose concentration drops. Water exits back to xylem via osmosis. Turgor drops. The resulting pressure gradient between source and sink drives bulk flow through the sieve tube — sucrose dissolved in water moves with the stream at rates of 0.3–1.5 m/hour in real trees.

Key properties:

- **Local:** Each sieve element only "knows" its pressure and the pressure of its immediate neighbors. There is no plant-wide sensor or scheduler.
- **Bidirectional:** Flow goes from high pressure to low pressure regardless of orientation in the plant. Phloem can flow upward (root starch → shoot during spring flush) or downward (leaf sugar → root growth) depending on where sources and sinks are.
- **Distance-sensitive:** Pressure drops along the path. A nearby sink gets a steeper gradient than a distant sink with equal demand. Proximity is a competitive advantage.
- **Water-coupled:** Osmotic pressure cannot build without water. A dehydrated source leaf cannot load phloem even if it has abundant sugar.

### 2.2 Loading and unloading

**Loading (active, source leaves):** Companion cells spend ATP to move sucrose against a concentration gradient from mesophyll into sieve tubes. Loading is energy-costly. Loading rate scales with: leaf photosynthesis rate (more production → more to load), sugar concentration above the storage reserve floor (can't load what you don't have), and water availability (no water = no osmotic pressure = no flow).

**Unloading (passive, sinks):** Sugar exits sieve tubes down its own concentration gradient into consuming cells. Fast-growing cells (meristems, expanding leaves) consume sugar rapidly, maintaining a steep unloading gradient. Resting cells have lower demand and unload slower. No energy cost — the plant doesn't pay to unload.

**Conduits (stems, older roots):** Maintain sugar at a low background concentration that keeps the sieve tube functional. They don't strongly load or unload. Their turgor equilibrates between adjacent source and sink pressures, transmitting the pressure gradient along the pipe.

---

## 3. Proposed Algorithm

### 3.1 Overview

Replace the `run_vascular` call for sugar with a new `run_munch_phloem` function. The new function makes two O(N) passes:

1. **Pressure pass:** For each vascular node, compute phloem pressure from sugar concentration, water availability, and (for leaves) active loading boost.
2. **Flow pass:** For each parent-child edge, compute bulk flow from the pressure differential and pipe conductance, determine how much sugar moves with the flow, and write to accumulator buffers.
3. **Apply:** Flush accumulator buffers into node `chemicals` map.

No Phase 1 demand/supply aggregation. No Phase 2 top-down distribution. Each edge is independent.

### 3.2 Pressure computation

```
// For each vascular node n:
sugar_conc  = n.sugar / sugar_cap(n, g)
water_avail = clamp(n.water / water_cap(n, g), 0.0, 1.0)
base_pressure = sugar_conc × g.phloem_osmotic_coefficient × water_avail

// Leaves add active loading boost
if n.type == LEAF and n.sugar > phloem_reserve × sugar_cap:
    surplus_frac = (n.sugar - phloem_reserve × sugar_cap - n.sugar_reserved_for_growth) / sugar_cap
    loading_boost = g.phloem_loading_rate × surplus_frac × water_avail
    phloem_pressure[n] = base_pressure + loading_boost
else:
    phloem_pressure[n] = base_pressure
```

The `water_avail` factor appears twice: once in the base pressure (the node's own osmotic potential depends on how full its phloem sieve tubes are with water) and in the loading boost (the companion cell pump is ineffective without water to draw into the tube). A drought-stressed leaf with full sugar but empty water storage has near-zero phloem pressure — it cannot export.

Non-vascular nodes (leaves, meristems) do not compute phloem pressure. They participate in last-mile delivery via local diffusion, same as today.

### 3.3 Flow computation

```
// For each parent-child edge where BOTH parent and child are vascular:
dp = phloem_pressure[parent] - phloem_pressure[child]
pipe_cap = π × parent.radius² × g.phloem_conductance
flow_volume = dp × pipe_cap   // positive = parent→child, negative = child→parent

// Sugar moves dissolved in bulk flow
if flow_volume > 0:   // parent exports to child
    source_conc = parent.sugar / sugar_cap(parent)
    raw_sugar = flow_volume × source_conc
    // Bound: respect grow reserve and phloem reserve floor
    max_export = parent.sugar - parent.sugar_reserved_for_growth
                             - g.phloem_reserve_fraction × sugar_cap(parent)
    sugar_moved = clamp(raw_sugar, 0.0, max(0.0, max_export))
    // Loading cost (ATP for active loading) — only applies when leaf is source
    if parent.type == LEAF and sugar_moved > 0:
        cost = sugar_moved × g.phloem_loading_cost
        sugar_moved -= cost      // net export is reduced
        delta_parent -= sugar_moved + cost
        delta_child  += sugar_moved
    else:
        delta_parent -= sugar_moved
        delta_child  += sugar_moved
else:  // child→parent
    source_conc = child.sugar / sugar_cap(child)
    raw_sugar = -flow_volume × source_conc   // positive magnitude
    max_export = child.sugar - child.sugar_reserved_for_growth
                            - g.phloem_reserve_fraction × sugar_cap(child)
    sugar_moved = clamp(raw_sugar, 0.0, max(0.0, max_export))
    delta_parent += sugar_moved
    delta_child  -= sugar_moved

// Accumulate — do NOT write yet (all flows use start-of-tick concentrations)
```

All edges are computed from start-of-tick concentrations and pressures. Deltas accumulate per node. After all edges are processed, apply deltas atomically:

```
for each node:
    n.sugar = clamp(n.sugar + delta[n], 0.0, sugar_cap(n))
```

This is functionally equivalent to using `transport_received` buffers — it prevents order-dependent artifacts where a node's sugar changes mid-pass would affect downstream pressure calculations.

### 3.4 What happens to non-vascular last-mile delivery?

The existing `Node::transport_with_children` handles local diffusion for auxin, GA, stress, and last-mile delivery of vascular chemicals to/from leaves and meristems. **This is unchanged.** Leaves and meristems are not in the vascular pressure network — they connect to the nearest vascular conduit via local diffusion (low sugar in a meristem node creates a diffusion gradient pulling sugar from the adjacent stem). The Münch pass moves sugar through the stem pipe network; diffusion handles the final step into the leaf mesophyll or meristem tissue.

### 3.5 Single-pass vs. iterative relaxation

**Single-pass** (recommended): Compute all edge flows once from start-of-tick pressures, apply all changes. This means pressure gradients propagate one edge per tick. A leaf surplus reaches the adjacent stem in one tick, the next stem in the next tick, and so on. This matches the "anti-teleportation" philosophy already in the sim — chemicals move at most one hop per tick per pass. The approximation error is bounded by the distance between source and sink.

**Iterative relaxation** (future option): Run 3–5 pressure-flow iterations per tick, re-computing pressures after each flow application. This allows multi-hop pressure propagation within a single tick, at the cost of O(N × k) instead of O(N). Useful if single-pass produces noticeable sourcing artifacts (e.g., the leaf directly below an apex gets all the sugar while a leaf one node further starves the apex). Implement only if single-pass shows this problem in practice.

The current vascular pass also uses a single-pass structure (Phase 2 is a single pre-order walk). Moving to Münch does not introduce a new performance class.

---

## 4. New Genome Parameters

| Parameter | Default | Biological basis |
|---|---|---|
| `phloem_osmotic_coefficient` | 3.0 | Maps sugar concentration [g/g] to dimensionless pressure. At full saturation (conc=1.0), base pressure = 3.0. This is the primary calibration knob — adjust to match observed flow rates. |
| `phloem_loading_rate` | 1.0 | Pressure boost per unit surplus fraction above phloem_reserve. A leaf at 100% surplus above reserve adds 1.0 to its pressure, amplifying the source gradient above conduit baseline. |
| `phloem_loading_cost` | 0.01 | Fraction of exported sugar consumed as ATP cost (active companion cell pump). 1% loading cost is on the low end of real estimates (1–5%) but won't significantly stress the sugar economy. |
| `phloem_water_coupling` | 0.8 | Fraction of phloem pressure that is suppressed under zero water. At 0.8: a fully dehydrated node runs phloem pressure at 20% of its sugar-only value. At 0: water availability is ignored. |

### Parameters to remove

| Parameter | Why removed |
|---|---|
| `meristem_sink_fraction` | Meristem pull is bounded by pressure gradient × pipe conductance, not by an explicit cap. Removing the cap lets small, nearby meristems compete based on concentration difference. |

### Parameters to simplify

| Parameter | Change |
|---|---|
| `phloem_reserve_fraction` | Kept but semantics shift: it is now the sugar concentration floor below which leaves do *not* load (i.e., the minimum they retain for their own cell osmotic needs). In the allocation model it was an explicit source-exclusion floor; in Münch it becomes the threshold that determines whether `loading_boost > 0`. Numeric value (0.3) is unchanged. |

---

## 5. Integration with Existing Systems

### 5.1 Xylem (water + cytokinin)

Xylem keeps Phase 1/Phase 2 unchanged. Xylem runs first in `vascular_transport`, before Münch phloem. This ordering is correct and critical: the Münch pressure calculation requires `n.water / water_cap(n, g)` at start-of-tick, after xylem delivery. Running xylem first ensures water is distributed before phloem pressure is computed.

The coupling is unidirectional in this design: xylem affects phloem (water availability shapes phloem pressure), but phloem does not affect xylem. In real plants, phloem flow returns water to xylem at sinks, creating a circulation. This coupling can be added as a future refinement by having the Münch pass record water flux alongside sugar flux — wherever sugar enters a sink via phloem, a proportional water volume exits back to xylem.

### 5.2 PIN transport and canalization

The Münch model uses the same `π × radius² × phloem_conductance` pipe capacity as the current Phase 2. The canalization feedback loop is unchanged:

```
PIN auxin flux → auxin_flow_bias → cambium activation → radius grows →
π × radius² increases → more phloem flow reaches this branch →
more sugar → more growth → more auxin → more PIN flux
```

The only change is that pipe capacity now drives a pressure-gradient flow rather than a conductance-weighted allocation. The competitive dynamics are the same: thicker pipes have higher conductance and carry more flow at the same pressure differential. The self-reinforcing hierarchy of main trunks vs. lateral branches is preserved and sharpened — a thick trunk conducts efficiently across a pressure gradient, a thin lateral has high resistance and receives less per unit gradient.

### 5.3 Starch storage (storage/plan.md)

Münch flow interacts with starch mobilization naturally, without special cases:

**Starch synthesis (surplus → storage):** A leaf or stem node with sugar above `storage_threshold × sugar_cap` converts surplus to starch. Local sugar drops. Sugar concentration drops. `phloem_pressure` at this node decreases. Reduced pressure means this node imports less from adjacent high-pressure sources and exports less to adjacent sinks. The node partially "bows out" of the phloem stream while it is banking reserves — exactly correct biology.

**Starch mobilization (storage → surplus):** A node under sugar deficit mobilizes starch to sugar. Local sugar rises. `phloem_pressure` rises. This node now pushes outward to adjacent lower-pressure sinks. A stem node with large starch reserves becomes a local phloem source during periods of photosynthetic deficit, without any new allocation logic. The Münch model automatically routes sugar from wherever pressure is highest to wherever it is lowest — mobilizing stems are just new local pressure peaks.

**Seed bootstrap:** When the seed mobilizes its starch reserve (via the GA-triggered mechanism described in storage/plan.md), the seed's sugar concentration rises, raising its phloem pressure. Adjacent stem and root nodes have lower initial sugar and therefore lower pressure. Sugar flows outward from the seed toward early growing regions. This replaces the current `seed_sugar` classification as a special phloem source in the allocation model. The seed becomes a natural source just by having higher concentration — no explicit source classification needed.

### 5.4 grow-before-transport reserve (`sugar_reserved_for_growth`)

The `pre_transport_growth` mechanism is unchanged. Each node computes `sugar_reserved_for_growth` before `vascular_transport` runs. The Münch pressure and export computations exclude reserved sugar:

```
max_export = n.sugar - n.sugar_reserved_for_growth - g.phloem_reserve_fraction × sugar_cap(n)
```

This ensures that a leaf that has reserved sugar for its own expansion this tick does not export that sugar via the phloem pass. Call order in `tick_tree` stays:

```
pre_transport_growth()      // nodes claim sugar_reserved_for_growth
vascular_transport()        // xylem Phase1/Phase2, then Münch phloem
pin_transport()             // PIN-mediated polar auxin
tick_recursive()            // DFS: grow, diffuse, canalize
```

### 5.5 Anti-teleportation correctness

The current allocation model modifies `n.chemical(ChemicalID::Sugar)` directly during Phase 2, which runs before any DFS tick. The Münch approach computes all edge flows from start-of-tick state, accumulates deltas, then applies them. This is strictly more correct: no edge's flow changes another edge's starting conditions mid-pass. The result is equivalent to applying all flows simultaneously from a snapshot, which is the physically correct interpretation for a discrete-time simulator.

Local diffusion (`transport_with_children`) continues to use `transport_received` buffers for its own hop-by-hop sugar exchange. These are two separate passes — vascular Münch runs first (global simultaneous update), then DFS local diffusion runs per-node (buffered per-hop). The two-buffer design is preserved.

---

## 6. Implementation Steps

### Task 1: Add genome parameters

**Files:**
- Modify: `src/engine/genome.h`
- Modify: `src/evolution/genome_bridge.cpp` (add/remove bridged genes)

- [ ] Add `phloem_osmotic_coefficient`, `phloem_loading_rate`, `phloem_loading_cost`, `phloem_water_coupling` to `Genome` struct in `genome.h`
- [ ] Set defaults in `default_genome()`: osmotic_coefficient=3.0, loading_rate=1.0, loading_cost=0.01, water_coupling=0.8
- [ ] Remove `meristem_sink_fraction` from struct and default (do a grep to find all use sites first: `grep -r meristem_sink_fraction src/`)
- [ ] Run build: `/usr/local/bin/cmake --build build` — fix all compilation errors from removed parameter
- [ ] Run tests: `./build/botany_tests` — all 149 tests should still pass (behavior unchanged until Task 3)
- [ ] Commit: `feat: add Münch phloem genome params, remove meristem_sink_fraction`

### Task 2: Add `compute_phloem_pressure` helper and unit test

**Files:**
- Modify: `src/engine/vascular.cpp`
- Modify: `src/engine/vascular.h`
- Modify: `tests/test_vascularization.cpp`

- [ ] Add static helper to `vascular.cpp`:

```cpp
static float compute_phloem_pressure(const Node& n, const Genome& g, float sugar_cap_val) {
    float sugar_conc = (sugar_cap_val > 0) ? n.chemical(ChemicalID::Sugar) / sugar_cap_val : 0.0f;
    float water_avail = std::clamp(
        n.chemical(ChemicalID::Water) / std::max(1e-6f, water_cap(n, g)),
        0.0f, 1.0f);
    // Base osmotic pressure: concentration × coefficient × water coupling
    float water_factor = 1.0f - g.phloem_water_coupling * (1.0f - water_avail);
    float base = sugar_conc * g.phloem_osmotic_coefficient * water_factor;
    // Leaves add active loading boost when above reserve
    float loading_boost = 0.0f;
    if (n.type == NodeType::LEAF) {
        float reserve = g.phloem_reserve_fraction * sugar_cap_val;
        float avail = n.chemical(ChemicalID::Sugar) - reserve - n.sugar_reserved_for_growth;
        if (avail > 0) {
            float surplus_frac = avail / sugar_cap_val;
            loading_boost = g.phloem_loading_rate * surplus_frac * water_avail;
        }
    }
    return base + loading_boost;
}
```

- [ ] Write unit test `test_phloem_pressure_from_concentration`: create two nodes (leaf with sugar=0.8×cap, stem with sugar=0.1×cap, both at full water), verify that leaf pressure > stem pressure, and that leaf pressure drops toward zero when water is zero
- [ ] Write unit test `test_phloem_pressure_water_coupling`: same two nodes, reduce water to zero on the leaf, verify pressure drops by `phloem_water_coupling` fraction
- [ ] Run tests: `./build/botany_tests` — new tests pass
- [ ] Commit: `feat: add compute_phloem_pressure helper`

### Task 3: Write `run_munch_phloem` function

**Files:**
- Modify: `src/engine/vascular.cpp`

This replaces `run_vascular` for sugar only. Xylem passes are unchanged.

- [ ] Add `run_munch_phloem` after `run_vascular` in `vascular.cpp`:

```cpp
static void run_munch_phloem(std::vector<VascNodeInfo>& flat, const Genome& g) {
    if (flat.empty()) return;

    // Pass 1: compute phloem pressure at each vascular node
    std::vector<float> pressure(flat.size(), 0.0f);
    for (int i = 0; i < (int)flat.size(); ++i) {
        const Node& n = *flat[i].node;
        if (!flat[i].is_conduit) continue;
        float cap = sugar_cap(n, g);
        pressure[i] = compute_phloem_pressure(n, g, cap);
    }

    // Pass 2: compute per-edge flows, accumulate deltas
    std::vector<float> delta(flat.size(), 0.0f);
    for (int i = 1; i < (int)flat.size(); ++i) {
        auto& child_info = flat[i];
        if (!child_info.is_conduit) continue;
        int pi = child_info.parent_idx;
        if (pi < 0 || !flat[pi].is_conduit) continue;

        float dp = pressure[pi] - pressure[i];
        float cap_child = pipe_capacity(*child_info.node, g.phloem_conductance);
        float cap_parent = pipe_capacity(*flat[pi].node, g.phloem_conductance);
        float pipe = std::min(cap_parent, cap_child);
        float flow_volume = dp * pipe;   // positive = parent→child

        if (flow_volume > 0.0f) {
            // Parent exports to child
            const Node& src = *flat[pi].node;
            float src_cap = sugar_cap(src, g);
            float src_conc = (src_cap > 0) ? src.chemical(ChemicalID::Sugar) / src_cap : 0.0f;
            float raw = flow_volume * src_conc;
            float max_export = src.chemical(ChemicalID::Sugar)
                - src.sugar_reserved_for_growth
                - g.phloem_reserve_fraction * src_cap;
            float moved = std::clamp(raw, 0.0f, std::max(0.0f, max_export));
            float cost = (src.type == NodeType::LEAF) ? moved * g.phloem_loading_cost : 0.0f;
            delta[pi] -= (moved + cost);
            delta[i]  += moved;
        } else if (flow_volume < 0.0f) {
            // Child exports to parent
            const Node& src = *child_info.node;
            float src_cap = sugar_cap(src, g);
            float src_conc = (src_cap > 0) ? src.chemical(ChemicalID::Sugar) / src_cap : 0.0f;
            float raw = -flow_volume * src_conc;
            float max_export = src.chemical(ChemicalID::Sugar)
                - src.sugar_reserved_for_growth
                - g.phloem_reserve_fraction * src_cap;
            float moved = std::clamp(raw, 0.0f, std::max(0.0f, max_export));
            float cost = (src.type == NodeType::LEAF) ? moved * g.phloem_loading_cost : 0.0f;
            delta[i]  -= (moved + cost);
            delta[pi] += moved;
        }
    }

    // Apply deltas atomically from snapshot
    for (int i = 0; i < (int)flat.size(); ++i) {
        if (delta[i] == 0.0f) continue;
        Node& n = *flat[i].node;
        float cap = sugar_cap(n, g);
        n.chemical(ChemicalID::Sugar) = std::clamp(
            n.chemical(ChemicalID::Sugar) + delta[i], 0.0f, cap);
    }
}
```

- [ ] Write unit test `test_munch_flow_leaf_to_meristem`: build a 3-node chain (leaf → stem → meristem), give leaf high sugar + high water, give meristem low sugar. Run `run_munch_phloem` once. Verify stem has gained sugar and leaf has lost sugar.
- [ ] Write unit test `test_munch_flow_bidirectional`: same topology, but reverse: give meristem high sugar (simulate starch mobilization), leaf low sugar. Verify flow goes meristem→stem direction.
- [ ] Write unit test `test_munch_no_flow_without_water`: leaf has high sugar but zero water. Verify flow is suppressed by `phloem_water_coupling` fraction.
- [ ] Write unit test `test_munch_reserve_floor_respected`: leaf at exactly `phloem_reserve_fraction × cap`. Verify no export occurs (`max_export` = 0 after reserve + growth reserve).
- [ ] Run tests: `./build/botany_tests`
- [ ] Commit: `feat: add run_munch_phloem function with unit tests`

### Task 4: Wire Münch into `vascular_transport`

**Files:**
- Modify: `src/engine/vascular.cpp`

- [ ] In `vascular_transport`, locate the sugar pass:
  ```cpp
  run_vascular(flat, ChemicalID::Sugar, g.phloem_conductance, g, log, world.current_tick);
  ```
- [ ] Replace it with:
  ```cpp
  // Phloem: Münch pressure-driven flow (replaces global allocation)
  // Each vascular edge carries flow proportional to local osmotic pressure differential.
  run_munch_phloem(flat, g);
  ```
- [ ] Remove the debug log parameter from the sugar call (Münch doesn't use `VascNodeInfo.demand/supply` fields). Update or disable the CSV log for the sugar system — xylem log lines remain.
- [ ] Run build: `/usr/local/bin/cmake --build build`
- [ ] Run all tests: `./build/botany_tests` — existing tests may need pressure-based sugar expectations adjusted
- [ ] Run visual check: `./build/botany_realtime --color sugar` and observe sugar distribution gradient (high at leaves, lower toward roots/meristems, with visible path-length effect)
- [ ] Run sugar economy check: `./build/botany_sugar_test --tree seedling --ticks 100` — production/maintenance ratio should remain in the ~4x range from earlier tuning
- [ ] Commit: `feat: replace phloem allocation with Münch pressure-driven flow`

### Task 5: Update `VascNodeInfo` — remove stale fields

**Files:**
- Modify: `src/engine/vascular.cpp`

After Münch replaces the sugar pass, the `supply` and `demand` fields of `VascNodeInfo` are only used by the xylem pass. The reset loops between passes are still correct. But the struct comment and any sugar-related field initialization should be cleaned up.

- [ ] Update `VascNodeInfo` comment to note that `supply`/`demand` are xylem-only fields now
- [ ] Verify the xylem reset loops (`info.supply = 0; info.demand = 0;`) still run correctly before the water and cytokinin passes
- [ ] Remove any sugar-related seed source logic from Phase 1 / Phase 2 of `run_vascular` that is now dead code (the seed's sugar surplus was handled there; Münch handles it implicitly via concentration)
- [ ] Grep for `meristem_sink_fraction` to confirm it is fully removed: `grep -r meristem_sink_fraction src/ tests/`
- [ ] Run build and tests
- [ ] Commit: `refactor: clean up VascNodeInfo after Münch phloem migration`

### Task 6: Tune osmotic coefficient against observed plant behavior

**Files:**
- Modify: `src/engine/genome.h` (adjust defaults only)

This step is calibration, not code architecture. The goal is that:
1. A mature leaf exporting at peak can supply 1–2 active meristems without starving neighboring conduit stems
2. A deep root meristem receives less per tick than a meristem adjacent to the seed — distance penalty is visible
3. Drought-stressed plants (soil_moisture = 0.1) show reduced sugar export rate from leaves (water_coupling working)

- [ ] Run `botany_realtime --color sugar` at default genome (phloem_osmotic_coefficient=3.0). Observe: do leaves maintain distinctly higher sugar than stems? Do meristems show lower sugar levels that recover over time?
- [ ] If leaves equilibrate too quickly to stem level (too much flow per tick): reduce `phloem_osmotic_coefficient` or `phloem_loading_rate`
- [ ] If meristems are starving despite healthy leaves: increase `phloem_loading_rate`
- [ ] If drought makes no visible difference to sugar distribution: increase `phloem_water_coupling`
- [ ] Update defaults in `default_genome()` to match observed good behavior
- [ ] Commit: `tune: calibrate Münch phloem osmotic coefficient defaults`

### Task 7: Integration tests for Münch-specific behavior

**Files:**
- Modify: `tests/test_vascularization.cpp`

- [ ] Write test `test_munch_distance_penalty`: grow a chain of 5 stem nodes with one leaf at tip and one active meristem at base. After N ticks, verify the intermediate stem nodes have a sugar gradient (high near leaf, low near meristem) — not uniform distribution as allocation would produce.
- [ ] Write test `test_munch_water_drought_reduces_flow`: grow same chain, set water=0 on leaf. Verify meristem receives less sugar than when water is full.
- [ ] Write test `test_munch_no_teleport`: single tick, leaf at far tip with sugar, meristem at far base with no sugar. Verify meristem sugar does not jump by more than one-hop flow rate (adjacent-to-leaf stem gains some, meristem gains nothing yet).
- [ ] Write test `test_munch_bidirectional_starch_source`: simulate starch mobilization in a mid-stem by setting its sugar to 0.9×cap (high), while both its parent (seed) and child (meristem) have low sugar. Verify sugar flows both toward seed AND toward meristem (bidirectional).
- [ ] Run: `./build/botany_tests`
- [ ] Commit: `test: add Münch pressure flow integration tests`

### Task 8: Update documentation and comments

**Files:**
- Modify: `src/engine/vascular.cpp` (block comment)
- Modify: `CLAUDE.md` (Chemical Transport Model section)
- Modify: `src/engine/vascular.h`

- [ ] Update the top comment in `vascular.cpp` to describe the Münch model, explain that phloem is now pressure-driven (not allocated), and note that xylem remains Phase 1/Phase 2
- [ ] Update the "Chemical Transport Model" and phloem/xylem table in `CLAUDE.md` to reflect the new phloem mechanism
- [ ] Update `is_vascular_chemical` comment in `vascular.h` if needed
- [ ] Commit: `docs: update Chemical Transport Model for Münch phloem`

---

## 7. Expected Behavior Changes

### What improves

**Source-proximity advantage.** A meristem adjacent to a productive leaf will be supplied faster than an equally-sized meristem two meters of stem away. The current allocation model ignores path length (sugar reaches the seed junction and distributes globally). Münch flow naturally attenuates with distance, matching the real competitive dynamics of branch architecture.

**Bidirectional phloem.** When a stem node has high sugar from starch mobilization (e.g., spring flush simulations), it can export toward the seed just as easily as toward root tips. Currently, stems are classified as conduits and cannot source sugar in the allocation pass — only leaves and the seed can source. Münch removes this classification constraint.

**Drought coupling.** A leaf with 20% water fills will export at roughly 20% of its hydrated rate even if sugar is abundant. This completes the hydraulic loop: xylem drought reduces phloem flow, reducing growth rate, reducing auxin production, reducing canalization flux — the plant correctly slows down under drought without requiring additional logic.

**No meristem tuning hacks.** With `meristem_sink_fraction` removed, meristems compete based on their actual size (pipe conductance to them) and their sugar deficit (concentration gradient). A microscopic dormant lateral bud gets almost nothing from phloem — it's tiny, its concentration gradient is small (it's not consuming rapidly), and its pipe is thin. An active growing tip gets a lot — it's consuming continuously, maintaining a steep concentration gradient, and the pipe connecting it to the main trunk has been canalized. The biologically correct competition emerges from the model.

**Natural loading saturation.** When a leaf's phloem concentration approaches the stem concentration (fully loaded conduit, low sink demand), the pressure differential collapses and loading slows. Currently, a leaf with surplus sugar always exports at full rate until it hits `phloem_reserve_fraction`. Münch naturally throttles: well-supplied conduits don't need more.

### What may need monitoring

**Early plant (small, 3–10 nodes).** The distance-penalty effect is negligible at small scale since every node is adjacent to every other. The main behavioral difference at seedling scale is the removal of `meristem_sink_fraction`. The active SAM should still receive sugar — it maintains a steep concentration gradient by growing and consuming. If the seedling starves more frequently in early ticks, tune `phloem_osmotic_coefficient` or `phloem_loading_rate` upward.

**Long-range phloem during rapid growth.** A tall plant (50+ nodes from canopy to root tips) may show root meristems receiving sugar more slowly than today because single-pass Münch only propagates one hop per tick. If root growth stalls due to sugar starvation in tall plants, consider adding 2–3 relaxation iterations within `run_munch_phloem`. This is the primary tradeoff of single-pass vs. relaxation.

**Sugar accumulation in stems.** Stems are now true conduits — their sugar concentration equilibrates between adjacent source and sink pressures. If leaves are very productive and sinks are limited, stems will accumulate sugar to high concentrations. This is biologically correct (parenchyma loading) and sets up the starch synthesis pathway naturally: when stems are above `storage_threshold × sugar_cap`, they synthesize starch. Monitor that stems don't permanently sit at cap (which would block phloem flow by eliminating the pressure gradient).

**Tests that check exact sugar transfer amounts.** Any test that asserts `sugar_delivered == X` from the allocation model will need its expected values recalculated. The test structure remains the same; the numeric expectations change because Münch distributes by pressure × conductance, not by demand-proportional allocation.

---

## 8. What Is Not Changing

- Xylem (water + cytokinin) keeps Phase 1 / Phase 2 global allocation. Transpiration pull is a different physical mechanism and is a future milestone.
- Local diffusion (`transport_with_children`) is unchanged. Auxin, GA, stress, and last-mile sugar/water delivery remain per-node.
- PIN transport pass is unchanged.
- Canalization model (`auxin_flow_bias`, radius-driven thickening) is unchanged.
- The `sugar_reserved_for_growth` / `pre_transport_growth` mechanic is unchanged.
- Gibberellin, ethylene, stress models are unchanged.
- The seed node classification as always-vascular is unchanged.

---

## Relationship to Other Milestone 2 Work

**Starch storage (storage/plan.md):** Münch is a prerequisite for starch storage to work correctly. Starch mobilization at a stem node should push sugar outward to sinks — this only works if the stem can be a phloem source, which requires pressure-driven flow. Under the allocation model, stems cannot source; under Münch, a high-sugar stem is automatically a source. Implement Münch first, then starch.

**World physics / soil model:** The gradient-based root absorption from vascularization/plan.md Step 2 ties into `world.soil_moisture`. The Münch `phloem_water_coupling` parameter ties the phloem flow to the xylem-delivered water level. Both systems use the same `water / water_cap` normalized metric. If the soil model adds spatial moisture variation, both xylem and phloem automatically respond — no additional coupling is needed.

**Production throttle:** The vascularization plan notes a future production throttle: when storage nodes are near capacity, leaves reduce photosynthesis. Under Münch, this signal is already partially implicit: when conduit stems are at high sugar concentration, the pressure gradient from leaf to stem collapses, and loading_boost from the leaf drops. Leaf phloem pressure naturally decreases relative to the stem, reducing outflow. The photosynthesis-reduction throttle would be an additional layer on top of this implicit signal, not a replacement for it.

**Evolution:** `meristem_sink_fraction` removal reduces the genome by one parameter. `phloem_osmotic_coefficient`, `phloem_loading_rate`, `phloem_loading_cost`, and `phloem_water_coupling` are added. Net: +3 evolvable parameters. All four new parameters have clear biological interpretations and should be registered in `genome_bridge.cpp` under the `sugar_economy` linkage group.
