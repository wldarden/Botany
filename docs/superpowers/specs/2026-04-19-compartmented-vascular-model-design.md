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

`N = world.vascular_substeps` — a fixed world parameter (starting value ~20–40, tuned empirically). Fixed N means plants whose longest source-to-sink chain exceeds N will have multi-tick pressure propagation delays, and distal apices will be partially under-supplied. **This is intentional and biologically realistic** — real tall plants have hydraulic limitations that cap their practical height, and a fixed N reproduces that constraint naturally. Plants that can't deliver to their apex will self-limit or show drought-stress patterns; evolution can favor genomes whose radial-permeability curve keeps trunks pipe-like enough to deliver despite the fixed-N budget (see section 6 for how the curve gives evolution control over this).

Each of the N iterations runs four steps, in order:

1. **Inject (sources):** sources that inject into a *different* node's conduit do it here. Two cases:
   - Leaf → walk-up parent's phloem.sugar: transfer `budget/N`. **Active**, meaning transfer happens regardless of gradient (real phloem loading is ATP-driven).
   - Root apical → walk-up parent's xylem.cytokinin: transfer `budget/N`. Passive-gradient, capped by budget slice.

   Root water loading is *not* a cross-node injection — water moves from the root's own `local_env` into the root's own `xylem` via radial flow (step 3). It's a same-node transfer, not a source injection into a different node.

2. **Extract (sinks):** every sink transfers up to `budget/N` from its walk-up parent's phloem/xylem into its own `local_env`. The transfer is capped by what's actually present in the parent conduit: `pulled = min(budget/N, parent_pool.chemical)`. If the local phloem is dry at this sub-step, the sink gets less than requested — that's the correct physics, not a bug.

3. **Radial flow** (stems and roots only): bidirectional passive gradient between a node's own `local_env` and its own `phloem` / `xylem`:
   ```
   Δsugar = radial_permeability_sugar(r) × (phloem.sugar/phloem_cap − local.sugar/local_cap) / N
   ```
   Rate-capped by pool sizes. Same pattern for water between `local_env` and `xylem`. This is how stem parenchyma refills from its sieve tubes and how root-absorbed water moves from `local_env` into the xylem stream.

   **Radial permeability is radius-dependent** (see section 6) so young thin stems leak freely into/out of their own tissue while mature thick trunks mostly pass flow through to downstream sinks. This is the hydraulic asymmetry that makes tall plants viable — without it, a mature trunk would siphon all the xylem water into its own parenchyma and starve the apex.

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

**Longitudinal (pipe) capacity and conductance.** For a stem or root of radius `r` and length `L`:

```
phloem_capacity   = π · r² · L · phloem_fraction
phloem_conductance(A, B) = phloem_capacity(min(A, B)) × (1 + canalization_weight × auxin_flow_bias[A → B])
```

Same pattern for xylem with its own `xylem_fraction`. `auxin_flow_bias` is the existing single-layer canalization bias on the parent node, keyed by child pointer (already in the code today; see `node.cpp:567` for update logic). Thin newborn stems have tiny capacity — they participate in the network but don't dominate until their radius grows.

**Radial permeability (asymmetric — decreases with radius).** The radial-flow step (5B.3) moves chemicals between a stem/root's own `local_env` and its own `phloem`/`xylem`. Permeability for this exchange is radius-dependent so young stems leak freely (get plenty of nutrients, grow fast) and mature trunks mostly conduct through (not stealing everything from distal sinks):

```
radial_permeability(r) = base_radial_permeability × ( radial_floor_fraction
                         + (1 − radial_floor_fraction) / (1 + (r / radial_half_radius)²) )
```

Shape:
- `r → 0`: permeability ≈ `base_radial_permeability` (newborn stems fully leaky)
- `r = radial_half_radius`: permeability ≈ `base × (floor + (1 - floor)/2)` (halfway through transition)
- `r → ∞`: permeability ≈ `base_radial_permeability × radial_floor_fraction` (mature trunks retain floor-level leak — enough to maintain themselves, not enough to starve distal sinks)

Independent instances of the formula for phloem-radial and xylem-radial (each with its own `base`, `floor`, and `half_radius` params), so the shape of each pipe's asymmetry can evolve independently.

**Why this matters for fixed-N sub-stepping.** With a fixed world-param N smaller than the plant's longest chain, a pressure wave takes several ticks to propagate from root to apex. Without the radius-dependent permeability, every intermediate stem would bleed significant flow into its own parenchyma as the wave passes through — starving distal sinks. With the curve: mature trunks act as near-pipes, young stems still participate in radial exchange, and the sim naturally produces the kind of hydraulic architecture real tall plants have.

**Evolution knobs.** Species can evolve different curves:
- Low `radial_half_radius` → rapidly seal the trunk at thin widths (aggressive tall-growing architectures)
- High `radial_half_radius` → keep stems leaky at larger widths (shrubby, herbaceous)
- High `radial_floor_fraction` → retain more nutrient access in mature tissue (long-lived trees with living inner wood)
- Low `radial_floor_fraction` → near-sealed mature pipes (apical-delivery specialists, possibly at cost of trunk longevity)

The fixed-N world param and the three radial-permeability genome params together form the "hydraulic budget" of the plant — they determine how tall a given genome can reliably grow before apical supply falls apart.

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

1. **Compartment invariants.** Enforce the residency table from section 1 — chemicals must only appear in the compartments they're allowed in. Specifically: sugar never appears in any `xylem` pool (it belongs in phloem and local_env only); water and cytokinin never appear in any `phloem` pool (they belong in xylem and local_env only); auxin, gibberellin, and stress never appear in any transport pool at all (they're local-signaling chemicals). Leftover sugar in phloem and leftover water in xylem between ticks is normal and expected — that's the in-transit content the Jacobi pass will keep moving next tick. Leaves, apicals, and root apicals return `nullptr` from `phloem()` and `xylem()`.
2. **Vascular-phase mass conservation.** For every chemical, sum across all pools (local + phloem + xylem) before Phase 2 equals the sum after Phase 2, within float epsilon. No production or consumption runs during Phase 2, so the equality is exact modulo rounding.
3. **Full-tick mass conservation.** Total sugar after tick = total sugar before tick + photosynthesis_this_tick − maintenance_this_tick − growth_consumption_this_tick. Same accounting for water.
4. **Pressure propagation within N.** Inject sugar at one end of a chain shorter than `vascular_substeps` (e.g., 10 stems with N=20). After one tick, the shoot-tip phloem pool has received meaningful sugar — confirms Jacobi reaches across chains that fit within N sub-steps.
5. **Radial flow equilibration.** A stem with full phloem and empty local_env converges toward equal concentrations over several ticks.
6. **Walk-up correctness.** A leaf attached to an apical whose parent is a stem identifies the stem as the phloem provider.
7. **Dormant meristems are inert.** A dormant meristem has zero phloem-demand budget, pulls nothing, consumes nothing.
8. **Distance-dependent supply under fixed N.** For a plant whose chain length exceeds `vascular_substeps`, distal apices receive less than proximal ones in the same tick — assert that the gradient of received sugar vs chain position is monotonically decreasing. Confirms hydraulic limitation is active.

**Tests requiring updates at Phase E:**

- `test_vascularization.cpp` — the four integration tests assume specific pairwise-Jacobi behavior (canalization ratchet, conductance-weighted distribution). Restate each in terms of the new algorithm's equivalent invariant. Conductance-weighted distribution, for instance, now manifests as higher-bias edges carrying more Jacobi flow per sub-step — the invariant is preserved, but the assertion needs to measure it differently.
- Any sugar-economy test measuring equilibrium states under old Jacobi dynamics — replace with equivalent invariant under the new model.

---

## Risks and open questions

- **Tuning.** New genome/world params introduced:
  - Radial-flow shape (per chemical class, two sets — one for phloem-radial, one for xylem-radial):
    - `base_radial_permeability_sugar`, `base_radial_permeability_water`
    - `radial_floor_fraction_sugar`, `radial_floor_fraction_water`
    - `radial_half_radius_sugar`, `radial_half_radius_water`
  - Pipe geometry: `phloem_fraction`, `xylem_fraction`
  - Source/sink targets: `leaf_reserve_fraction`, `sink_target` (per meristem type), `turgor_target`
  - Sub-step count: `vascular_substeps` (world param — fixed per tick, e.g., 20–40 as starting value)

  First-pass values for longitudinal conductance can be derived from the current `phloem_conductance` and `xylem_conductance` constants so behavior roughly matches. Radial-flow shape params are new; starting values in the spec body (`base=1.0, floor=0.1, half_radius=0.3`) are the suggested default curve and will need calibration via `botany_sugar_test` and watching realtime viewer after Phase E.
- **Performance.** Each sub-step visits every tree edge once for Jacobi (O(P)) plus every source/sink/conduit once for inject/extract/radial (also O(P)). N sub-steps per tick gives total work of `O(N × P)` where P = node count and N = longest chain depth. For a typical branched plant, N grows much more slowly than P (bushy: ~log P, typical: ~√P), so the per-tick cost is effectively linear-ish in plant size. The pathological case is a single unbranched chain where N = P and work becomes quadratic — but real plants branch, so this shouldn't happen in normal use. For a 1000-node plant with depth ~50, that's ~50k pair-ops per tick, well within C++'s budget. Worth measuring on the realtime viewer after cutover; if it stutters, the fix is likely reducing how often the plant is re-flattened into an array between sub-steps, not changing the algorithm.
- **Interaction with `spawn_child` during tick.** Phase 1 can add new nodes. Phase 2 must walk the plant after those additions are flushed — same as today. No change in this ordering.
- **Edge case: plant has no leaves.** Sources list is empty; injection does nothing; extraction finds empty phloem and delivers nothing. Plant starves, as today. Not a bug, test for clean behavior.

---

## Self-review

- **Placeholder scan:** no TBD, TODO, or unfilled sections. Every numerical formula is concrete.
- **Internal consistency:** pool placement table matches the interface definition. Tick ordering description matches the migration plan. Budget table in section 5A matches the sink/source distinctions in section 4.
- **Scope:** this is one coherent refactor covering compartmentalization + vascular algorithm + tick ordering. All three are entangled; doing them separately would require temporary shims. A single plan is appropriate.
- **Ambiguity:** "parent" in the walk-up definition means immediate tree parent, as in existing code. "Adjacent" in the Jacobi description means every (parent, child) edge in the tree, not spatial adjacency. Both clarified inline.
