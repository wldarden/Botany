# Compartmented Vascular Model

**Status:** Design spec, awaiting review
**Date:** 2026-04-19
**Supersedes:** [2026-04-19-demand-driven-phloem.md](../plans/2026-04-19-demand-driven-phloem.md) (a tactical fix that would have masked the structural issue this spec addresses)

---

## Motivation

The current model keeps a single chemical pool per node (`node.chemical(id)`) that conflates three physically distinct things: the tissue's own metabolic supply, the sieve-tube content carrying sugar longitudinally, and the vessel content carrying water. Every vascular bug we've hit recently has the same shape: a fix that protects one of those uses breaks another.

Recent history of the symptom:

- `1f68afc` disabled local sugar diffusion so phloem could "fully own" transport — reverted, broke other flows.
- `954d85b` removed conduit demand from phloem allocation so meristems could be fed — reverted, stems starved.
- The demand-driven plan at [`plans/2026-04-19-demand-driven-phloem.md`](../plans/2026-04-19-demand-driven-phloem.md) would have replaced pairwise Jacobi with a global supply/demand allocator. Correct in outcome but abandons Münch pressure-flow physics in favor of accounting, which diverges further from biology and forfeits the natural "distal sinks get less" behavior that real plants show.

The structural root cause is that we don't have separate compartments to model. A stem's parenchyma sugar and its sieve-tube sugar are different physical pools with different rates of change, interacting via radial flow. Once compartments are modeled explicitly, Münch pressure-flow with sub-stepped Jacobi resolution becomes the natural algorithm — no global allocation, no ambiguous reserves, no leaks.

## Goals

1. Represent phloem, xylem, and parenchyma (local) as three separate pools, with compile-time guarantees about which chemicals can live in which pool.
2. Model longitudinal transport as Münch pressure-flow, using sub-stepped Jacobi iteration to propagate pressure waves across the full plant within one tick.
3. Keep all routing decisions local: sources inject, sinks extract, pressure gradients do the rest. No global supply/demand matching.
4. Eliminate the classes of bugs where fixing one compartment's behavior breaks another, by making compartments distinct objects.
5. Preserve existing canalization, cambial thickening, gibberellin, ethylene, and hormone signaling unchanged.

## Non-goals

- Water osmotic potential modeling (xylem "pressure" remains approximated as a concentration-like scalar).
- Per-chemical active transport mechanisms beyond the two already identified (leaf phloem loading, sink-side active pull).
- Changes to the node class hierarchy beyond adding pool members to STEM and ROOT subclasses.
- Changes to evolution, fitness functions, genome serialization, renderer, or any app beyond what's needed to keep them compiling.

---

## Design

### 1. Compartment types

Two new value types:

```cpp
struct LocalEnv {
    std::unordered_map<ChemicalID, float> chemicals;
};

struct TransportPool {
    std::unordered_map<ChemicalID, float> chemicals;
};
```

**Per-node-type placement:**

| Node type        | local_env | phloem | xylem |
|------------------|:---------:|:------:|:-----:|
| `StemNode`       | ✓         | ✓      | ✓     |
| `RootNode`       | ✓         | ✓      | ✓     |
| `LeafNode`       | ✓         | —      | —     |
| `ApicalNode`     | ✓         | —      | —     |
| `RootApicalNode` | ✓         | —      | —     |

**Per-chemical residency:**

| Chemical     | local_env | phloem | xylem |
|--------------|:---------:|:------:|:-----:|
| Sugar        | ✓         | ✓      | —     |
| Water        | ✓         | —      | ✓     |
| Cytokinin    | ✓         | —      | ✓     |
| Auxin        | ✓         | —      | —     |
| Gibberellin  | ✓         | —      | —     |
| Stress       | ✓         | —      | —     |
| Ethylene     | —         | —      | —     | *(spatial gas, separate system, unchanged)*

### 2. Node interface

```cpp
class Node {
public:
    LocalEnv& local();
    virtual TransportPool* phloem() { return nullptr; }
    virtual TransportPool* xylem()  { return nullptr; }
    // ... rest unchanged
};

class StemNode : public Node {
    TransportPool phloem_pool_;
    TransportPool xylem_pool_;
public:
    TransportPool* phloem() override { return &phloem_pool_; }
    TransportPool* xylem()  override { return &xylem_pool_;  }
};

class RootNode : public Node {
    // identical pool members and overrides as StemNode
};
```

`LeafNode`, `ApicalNode`, `RootApicalNode` do not override. A call to `leaf.phloem()` returns `nullptr`, which is an immediately visible bug in debug builds.

**All existing `node.chemical(id)` call sites migrate to `node.local().chemical(id)`.** This is mechanical find-and-replace. Explicit compartment access at every call site is self-documenting: future readers see which compartment each piece of code is reading or writing.

### 3. Walk-up resolution

Specialty nodes (leaves, meristems) need to locate the nearest ancestor with a transport pool in order to load or unload. Since a leaf's direct parent may be an `ApicalNode` (no pools), the lookup walks up the parent chain:

```cpp
TransportPool* Node::nearest_phloem_upstream() {
    for (Node* n = parent; n; n = n->parent) {
        if (auto p = n->phloem()) return p;
    }
    return nullptr;
}
```

Same pattern for `nearest_xylem_upstream()`. Recomputed on demand each sub-step (caching considered and rejected until profiling shows a need).

### 4. Per-tick flow

**Ordering B: tick-then-vascular.**

```cpp
void Plant::tick(const WorldParams& world) {
    tick_tree(seed_, world);            // Phase 1: per-node metabolism
    vascular_sub_stepped(*this, world); // Phase 2: sub-stepped transport
}
```

**Phase 1 (unchanged from today's behavior):** DFS walk from seed. Each node's `tick()`:
1. `age++`, sync world position
2. Produce (photosynthesis for leaves, water absorption for roots, hormone production)
3. Pay maintenance from `local_env`; starvation check
4. `update_tissue()` — type-specific growth; consumes from `local_env`; may spawn children
5. Sync position, update physics
6. Local signaling diffusion (auxin, gibberellin, stress) between neighboring `local_env`s
7. Decay, ethylene spatial pass

**Phase 2 (new):** Sub-stepped vascular transport. See section 5.

**Consequence — 1-tick lag:** Sugar photosynthesized in tick N's Phase 1 sits in leaf `local_env` until tick N+1 Phase 2 loads it into phloem. Sinks fed in tick N Phase 2 consume that sugar in tick N+1 Phase 1. At 1-hour tick granularity this is physically defensible — real plants buffer more than an hour of sugar in mesophyll cells.

### 5. Sub-stepped vascular transport

**Part A — Budget snapshot** (runs once, before the sub-step loop)

Walk the plant to compute each node's supply and/or demand for this tick. Budgets are frozen here; the sub-step loop amortizes them without re-evaluating.

**Phloem (sugar):**

| Node              | Role   | Budget |
|-------------------|--------|--------|
| Leaf              | Source | `max(0, local.sugar − leaf_reserve)` |
| Active meristem   | Sink   | `max(0, sink_target − local.sugar)` |
| Dormant meristem  | —      | 0 |
| Stem / root       | (radial sink only, not direct) | see Radial Flow |

**Xylem (water + cytokinin):**

| Node              | Role   | Budget |
|-------------------|--------|--------|
| Root              | (radial source only) | see Radial Flow |
| Root apical       | Source | cytokinin produced this tick, sitting in local_env |
| Leaf              | Sink   | `max(0, turgor_target − local.water)` — local_env was depleted by Phase 1 transpiration and photosynthesis |
| Active meristem   | Sink   | `max(0, turgor_target − local.water)` |

**Part B — Sub-step loop**

Compute `N = max_chain_length(plant) + 2` (plant walks itself once before the loop to find its deepest shoot and root chains).

Each of the N iterations runs four steps, in order:

1. **Inject (sources):** sources that inject into a *different* node's conduit do it here. Two cases:
   - Leaf → walk-up parent's phloem.sugar: transfer `budget/N`. **Active**, meaning transfer happens regardless of gradient (real phloem loading is ATP-driven).
   - Root apical → walk-up parent's xylem.cytokinin: transfer `budget/N`. Passive-gradient, capped by budget slice.

   Root water loading is *not* a cross-node injection — water moves from the root's own `local_env` into the root's own `xylem` via radial flow (step 3). It's a same-node transfer, not a source injection into a different node.

2. **Extract (sinks):** every sink transfers up to `budget/N` from its walk-up parent's phloem/xylem into its own `local_env`. The transfer is capped by what's actually present in the parent conduit: `pulled = min(budget/N, parent_pool.chemical)`. If the local phloem is dry at this sub-step, the sink gets less than requested — that's the correct physics, not a bug.

3. **Radial flow** (stems and roots only): bidirectional passive gradient between a node's own `local_env` and its own `phloem` / `xylem`:
   ```
   Δsugar = radial_permeability_sugar × (phloem.sugar/phloem_cap − local.sugar/local_cap) / N
   ```
   Rate-capped by pool sizes. Same pattern for water between `local_env` and `xylem`. This is how stem parenchyma refills from its sieve tubes and how root-absorbed water moves from `local_env` into the xylem stream.

4. **Longitudinal Jacobi:** one pass over every (parent, child) pair of transport nodes along the stem/root backbone:
   ```
   p_parent = parent.phloem.sugar / phloem_capacity(parent)
   p_child  = child.phloem.sugar  / phloem_capacity(child)
   flow     = phloem_conductance(parent, child) × (p_parent − p_child)
   parent.phloem.sugar −= flow
   child.phloem.sugar  += flow
   ```
   Same formulation for xylem with its own capacity and conductance. **Jacobi has no awareness of sources or sinks.** It is a pure neighbor equalizer; routing emerges from the pressure field created by injection and extraction at step-1 and step-2 locations.

### 6. Capacity and conductance

For a stem or root of radius `r` and length `L`:

```
phloem_capacity   = π · r² · L · phloem_fraction
phloem_conductance(A, B) = phloem_capacity(min(A, B)) × (1 + canalization_weight × auxin_flow_bias[A → B])
```

Same pattern for xylem with its own `xylem_fraction`. `auxin_flow_bias` is the existing single-layer canalization bias on the parent node, keyed by child pointer (already in the code today; see `node.cpp:567` for update logic). Thin newborn stems have tiny capacity — they participate in the network but don't dominate until their radius grows.

### 7. Retirements

- **`has_vasculature()`** — no longer needed. Pool ownership is type-driven; `vascular_radius_threshold` has no job.
- **`vascular_radius_threshold` genome param** — remove.
- **`phloem_iterations` world param** — already noted as deprecated; remove.
- **`transport_received` buffers on every node** — were needed for the anti-teleportation pattern under per-node transport. The new model does all vascular work in a single centralized function with explicit pool reads/writes, so there is no same-tick cascade risk. Remove.
- **Old `phloem_resolve` / `xylem_resolve` functions** — replaced entirely by `vascular_sub_stepped`.

### 8. Mass conservation

Every transfer in the vascular phase is a paired `+= x` / `-= x` between two named pools:
- Inject: source `local_env` → parent conduit
- Extract: parent conduit → sink `local_env`
- Radial flow: node's conduit ↔ same node's `local_env`
- Jacobi: parent conduit ↔ child conduit

No step creates or destroys chemical. A single assertion at the end of Phase 2 — sum the chemical across all pools and compare to the pre-phase-2 sum — catches any arithmetic bug in one line. This assertion is a new test.

---

## Migration plan

Each phase leaves the test suite green. Commit after each.

**A. Additive types.** Create `LocalEnv`, `TransportPool` structs. Add `Node::local()` returning the existing chemical map as a `LocalEnv` alias. Add virtual `Node::phloem()` and `Node::xylem()` returning `nullptr`. Zero behavior change.

**B. Chemical access sweep.** Find-and-replace `node.chemical(id)` → `node.local().chemical(id)` across the codebase. Mechanical; each site is trivial. Zero behavior change.

**C. Pool members on StemNode / RootNode.** Add `TransportPool phloem_pool_`, `TransportPool xylem_pool_` to both classes. Override `phloem()` and `xylem()` to return them. Pools are unused by any transport code yet. Additive.

**D. Implement `vascular_sub_stepped`.** New function in `vascular.cpp`. Reads and writes the new pools. Write new tests (see Test Strategy below). Not wired into `Plant::tick()` yet — the existing `vascular_transport` still runs. Both algorithms compile, both test suites pass.

**E. Cutover commit.** In `Plant::tick()`, swap the call order and replace `vascular_transport()` with `vascular_sub_stepped()`. This is the behavior-changing commit. Expect several existing tests to need updates — they asserted end-states specific to the old pairwise Jacobi or to vascular-before-tick ordering. Update them to express the equivalent invariants under the new model.

**F. Dead-code removal.** Delete `has_vasculature()`, `vascular_radius_threshold`, old `phloem_resolve`, old `xylem_resolve`, `transport_received` buffers, `phloem_iterations`.

**G. Documentation.** Rewrite the "Chemical Transport Model", "Canalization Model", and "Tick Control Flow" sections of `CLAUDE.md` to reflect the new architecture. Remove the obsolete two-layer canalization description.

---

## Test strategy

**Regression guards (kept unchanged, must continue passing):**

- 1000-tick integration smoke test (plant survives, grows, reasonable shoot/root ratio)
- Long-chain phloem delivery test from the superseded plan doc — apex on a 15-stem chain must receive non-zero sugar within one tick
- Auxin sensitivity, gibberellin, ethylene, serializer, evolution tests

**New tests (written during Phase D):**

1. **Compartment invariants.** After a tick, no sugar in any xylem pool; no water in any phloem pool; leaves, apicals, and root apicals return `nullptr` from `phloem()` and `xylem()`.
2. **Vascular-phase mass conservation.** For every chemical, sum across all pools (local + phloem + xylem) before Phase 2 equals the sum after Phase 2, within float epsilon. No production or consumption runs during Phase 2, so the equality is exact modulo rounding.
3. **Full-tick mass conservation.** Total sugar after tick = total sugar before tick + photosynthesis_this_tick − maintenance_this_tick − growth_consumption_this_tick. Same accounting for water.
4. **Pressure propagation.** Inject sugar at one end of a 30-stem chain (seed with high sugar, empty elsewhere). After one tick, the shoot-tip phloem pool has non-zero sugar.
5. **Radial flow equilibration.** A stem with full phloem and empty local_env converges toward equal concentrations over several ticks.
6. **Walk-up correctness.** A leaf attached to an apical whose parent is a stem identifies the stem as the phloem provider.
7. **Dormant meristems are inert.** A dormant meristem has zero phloem-demand budget, pulls nothing, consumes nothing.
8. **Adaptive N.** For a plant with a 100-stem chain, `vascular_sub_stepped` chooses `N ≥ chain_length` and the apex receives sugar in one tick.

**Tests requiring updates at Phase E:**

- `test_vascularization.cpp` — the four integration tests assume specific pairwise-Jacobi behavior (canalization ratchet, conductance-weighted distribution). Restate each in terms of the new algorithm's equivalent invariant. Conductance-weighted distribution, for instance, now manifests as higher-bias edges carrying more Jacobi flow per sub-step — the invariant is preserved, but the assertion needs to measure it differently.
- Any sugar-economy test measuring equilibrium states under old Jacobi dynamics — replace with equivalent invariant under the new model.

---

## Risks and open questions

- **Tuning.** New params introduced: `radial_permeability_sugar`, `radial_permeability_water`, `phloem_fraction`, `xylem_fraction`, `leaf_reserve_fraction`, `sink_target` (per meristem type), `turgor_target`. First-pass values can be derived from the current `phloem_conductance` and `xylem_conductance` constants so behavior roughly matches. Expect a tuning pass using `botany_sugar_test` after Phase E.
- **Performance.** Sub-stepped vascular does `O(N × plant_size)` work per tick where N is the longest chain. On a 1000-node plant with N=50, that's ~50k pair-ops per tick — well within budget for C++, but worth measuring on the realtime viewer after cutover. If the viewer stutters, the fix is almost certainly reducing how often the full plant is rebuilt into a flat array, not changing the algorithm.
- **Interaction with `spawn_child` during tick.** Phase 1 can add new nodes. Phase 2 must walk the plant after those additions are flushed — same as today. No change in this ordering.
- **Edge case: plant has no leaves.** Sources list is empty; injection does nothing; extraction finds empty phloem and delivers nothing. Plant starves, as today. Not a bug, test for clean behavior.

---

## Self-review

- **Placeholder scan:** no TBD, TODO, or unfilled sections. Every numerical formula is concrete.
- **Internal consistency:** pool placement table matches the interface definition. Tick ordering description matches the migration plan. Budget table in section 5A matches the sink/source distinctions in section 4.
- **Scope:** this is one coherent refactor covering compartmentalization + vascular algorithm + tick ordering. All three are entangled; doing them separately would require temporary shims. A single plan is appropriate.
- **Ambiguity:** "parent" in the walk-up definition means immediate tree parent, as in existing code. "Adjacent" in the Jacobi description means every (parent, child) edge in the tree, not spatial adjacency. Both clarified inline.
