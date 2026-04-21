# Node Inspector Panel & Overlay Redesign

**Date:** 2026-04-21
**Status:** design
**Scope:** `src/app_realtime.cpp`, `src/renderer/`, `src/engine/node/`, `src/engine/vascular_sub_stepped.*`, `src/engine/pin_transport.*`, plus new unit tests

## Motivation

Post-compartmented-vascular refactor (2026-04-19) and post-hormone-biology refactor (2026-04-20), the realtime node inspector panel and overlays only read `node.local()`. Stems that carry sugar in their `phloem` and cytokinin/water in their `xylem` show as empty in both the panel and the heatmaps. Transport-capacity cells use the diffusion formula even on conduit edges, overstating capacity for thick stems and understating it for thin ones. Several chemicals have no per-tick production/consumption tracking at all — only auxin, cytokinin, and sugar do, and of those only sugar has a full budget (maintenance / activity / transport).

This project brings the realtime UI into agreement with the current engine model and expands it to show the information the user actually wants when inspecting a plant.

## Goals

1. **Accurate chemical readout.** Show Local, Vasculature, and Total levels per chemical per node. For stems/roots, Vasculature reads their own xylem/phloem; for leaves/meristems/seed, Vasculature reads the nearest upstream stem/root's conduit pools.
2. **Per-tick production and consumption for every chemical at every node.** Replace scattered one-off counters with a generalized per-chemical pair on `Node`.
3. **Correct transport capacity readout on each pathway.** Jacobi (sugar/water/CK), PIN (auxin), and diffusion (GA/Eth/Stress) each have their own capacity formula; the UI must use the right one per chemical.
4. **Overlay modes that make the new information spatial.** Top-level categories with sub-pickers for chemical and scope, instead of a flat 11-item dropdown.
5. **Richer identity info on the panel.** Height in meters, nodes-to-seed, maturity (three senses of the word), and activation/activity satisfaction %.

Out of scope: changing engine mechanics, changing the genome, changing the renderer pipeline beyond what the overlay categories require.

## Architecture

Two-pronged:

### A. Engine instrumentation (backend)

Backend changes stay read-only with respect to engine dynamics. No mass-balance behavior changes; only new counters and helper getters.

### B. UI redesign (frontend)

All changes live in `src/app_realtime.cpp` (panel + overlay selector) and small helper functions. The renderer's `set_color_mode(ChemicalAccessor)` API is unchanged — we just build more accessor closures.

## Section 1 — Engine instrumentation

### 1.1 Generalized per-tick chemical counters

Add to `Node` (in `src/engine/node/node.h`):

```cpp
// Reset to zero at the start of Plant::tick_tree().
std::array<float, static_cast<size_t>(ChemicalID::Count)> tick_chem_produced{};
std::array<float, static_cast<size_t>(ChemicalID::Count)> tick_chem_consumed{};
```

`ChemicalID::Count` must exist as the enum's count sentinel. If it does not today, add it.

**Migration of existing counters:**

| Existing field | Migration |
|---|---|
| `tick_auxin_produced` | Also increment `tick_chem_produced[Auxin]`. Keep the named field for backward compatibility with existing call sites, but panel reads the array. After migration is stable, remove the named field. |
| `tick_cytokinin_produced` | Same — mirror into `tick_chem_produced[Cytokinin]`. |
| `tick_sugar_maintenance` | Keep (signed, carries intent). Also mirror into `tick_chem_consumed[Sugar]`. |
| `tick_sugar_activity` | Keep (signed net for photosynthesis vs growth). When positive → `tick_chem_produced[Sugar] += amount`; when negative → `tick_chem_consumed[Sugar] += -amount`. |
| `tick_sugar_transport` | Keep (signed net). This is NOT mirrored into produced/consumed — it's a separate "transport column" concept, since transport moves mass, not creates/destroys it. |

**New instrumentation sites:**

| Chemical | Produced at | Consumed at |
|---|---|---|
| Sugar | `LeafNode::photosynthesize` | `Node::pay_maintenance`, `*::elongate`/`::thicken`/`::expand`, `meristem::grow` |
| Water | `RootNode::absorb`, `RootApicalNode::absorb` | `LeafNode::transpire`, `LeafNode::photosynthesize` (water cost) |
| Auxin | `ApicalNode::tick`, `LeafNode::tick`, `RootApicalNode::tick` | `Node::decay_chemicals` (decay fraction), PIN outflow residuals |
| Cytokinin | `RootApicalNode::tick` | `Node::decay_chemicals` (decay fraction in local() only; xylem pool has no decay) |
| GA | `compute_gibberellin` (leaf emitter) | Whole-chemical reset each tick — count the pre-reset value as consumed |
| Ethylene | `compute_ethylene` (leaf emitter) | Same reset-each-tick rule |
| Stress | `update_physics` (mechanical load) | `Node::decay_chemicals` (if stress has a decay rate; otherwise reset-each-tick semantics apply) |

Each site increments at the exact point where the chemical value moves.

### 1.2 Vascular-scope read helper

Add to `src/engine/node/node.h` or a new small header:

```cpp
// Returns the conduit pool this chemical lives in for this node.
// Stems/roots -> own pool; leaves/meristems/seed specialty -> walk up to
// nearest ancestor stem/root. Returns nullptr for signaling chemicals
// (Auxin, Gibberellin, Ethylene, Stress) which have no vascular pool, or
// when no upstream conduit exists.
const TransportPool* vascular_scope(const Node& n, ChemicalID chem);
```

Implementation: dispatch by `chem`:
- `Sugar` → `nearest_phloem_upstream(n)`
- `Water`, `Cytokinin` → `nearest_xylem_upstream(n)`
- `Auxin`, `Gibberellin`, `Ethylene`, `Stress` → `nullptr`

`nearest_phloem_upstream()` / `nearest_xylem_upstream()` are already used internally by `vascular_sub_stepped.cpp`. Expose them (or wrap them) so the UI can call them.

### 1.3 Per-edge transport flux and capacity

For the "Transport Capacity Used" overlay. Per parent node, keyed by child pointer:

```cpp
// On Node, parallel in shape to the existing `last_auxin_flux`.
std::array<std::unordered_map<Node*, float>, ChemicalID::Count> tick_edge_flux;
std::array<std::unordered_map<Node*, float>, ChemicalID::Count> tick_edge_cap;
```

`tick_edge_flux[chem][child]` = signed amount that moved across (parent → child) this tick. `tick_edge_cap[chem][child]` = the theoretical cap for that edge and pathway this tick.

**Instrumentation sites:**

| Pathway | Chemicals | flux source | cap source |
|---|---|---|---|
| Jacobi phloem | Sugar | longitudinal phloem Jacobi pass in `vascular_sub_stepped.cpp` | `phloem_conductance × cross_section × dt_per_substep` summed across sub-steps |
| Jacobi xylem | Water, Cytokinin | longitudinal xylem Jacobi pass | `xylem_conductance × cross_section × dt_per_substep` |
| PIN | Auxin | `pin_transport.cpp` actual throughput | `r² × pin_capacity_per_area` (child radius squared, using edge connection radius) |
| Diffusion | GA, Ethylene, Stress | `compute_transport_flow` actual | same formula `base + radius_factor × scale` |

The sum across sub-steps applies to Jacobi: each sub-step carries ≤ `conductance × cross_section × Δp`; we add them all up, so the tick-total cap is the sum across sub-steps of whatever cap applied in each one. Simplest correct cap = `conductance × cross_section × N_substeps` (treats every sub-step as if it could have moved cap).

**Aggregated node utilization** for a chemical (for the overlay):

```
util(node, chem) = sum over edges touching node of tick_edge_flux[chem][edge]
                 / max(sum over edges touching node of tick_edge_cap[chem][edge], epsilon)
```

`epsilon` prevents div-by-zero for chemicals the node cannot transport at all (returns 0, rendered neutral gray).

### 1.4 Reset discipline

At the top of `Plant::tick_tree()` (before Phase 0 PIN transport):

```cpp
for_each_node([](Node& n) {
    n.tick_chem_produced.fill(0.0f);
    n.tick_chem_consumed.fill(0.0f);
    for (auto& m : n.tick_edge_flux) m.clear();
    for (auto& m : n.tick_edge_cap)  m.clear();
    n.last_auxin_flux.clear();
    // existing scalar resets
    n.tick_sugar_maintenance = 0;
    n.tick_sugar_activity    = 0;
    n.tick_sugar_transport   = 0;
    n.tick_auxin_produced    = 0;
    n.tick_cytokinin_produced = 0;
});
```

Single loop. Fold all existing ad-hoc resets into it if they aren't already here.

### 1.5 Testing

New test file `tests/test_tick_counters.cpp`:

- Fresh tree after one tick: assert every `tick_chem_produced` and `tick_chem_consumed` is `≥ 0`, and that `tick_chem_produced[Sugar]` matches the photosynthesis value for a single-leaf plant.
- Mass balance per chemical on a small hand-built tree: `Δlocal(chem) == produced - consumed + net_transport_in` within a float tolerance.
- Per-edge cap > 0 for any conduit edge with a pool pair for that chemical.
- Existing counters (`tick_sugar_maintenance`, `tick_auxin_produced`, etc.) continue to return the same numbers as before the refactor on a smoke-test scenario.

## Section 2 — Panel content

Panel order, top to bottom:

### 2.1 Identity section

- `ID: #<id>   Type: <TYPE>`
- `Age: <ticks> (<days>)`
- `Length: <fmt_dist>   Radius: <fmt_dist>`
- `Cross-section: <π r² in dm²>`
- `Height (y): <y/10.0> m` — world y, converted dm → m
- `Nodes to seed: <walk_up_count>` — count `parent` hops until `parent == nullptr`
- Stems/roots only: three maturity lines
  - `Elongation: mature (age X/Y)` or `growing (age X/Y)`
  - `Hydraulic maturity: NN% closed` — computed from `1 - radial_permeability_sugar(radius, genome) / base_radial_permeability_sugar`
  - `PIN saturation: 0.XX` — `get_parent_auxin_flow_bias()`
- `Starvation: N ticks`
- Leaves only: `Senescence: N / duration ticks` or `healthy`

Remove the standalone `Children: N` line (redundant with the navigation buttons).

### 2.2 Chemicals table

Fixed shape: 7 rows × 10 columns.

Column layout: `[Chem] [Local-lvl] [Local-+] [Local--] [Vasc-lvl] [Vasc-+] [Vasc--] [Total-lvl] [Total-+] [Total--]`

Rows (in order): Sugar, Water, Auxin, Cyt, GA, Eth, Stress.

**Value sources:**

| Column | Source |
|---|---|
| Local-lvl | `node.local().chemical(chem)` |
| Local-+ | `node.tick_chem_produced[chem]` |
| Local-- | `node.tick_chem_consumed[chem]` |
| Vasc-lvl | `vascular_scope(node, chem) ? that->chemical(chem) : "—"` |
| Vasc-+ | For stems/roots that own the conduit: the produced amount into their own pool. For nodes that use an upstream pool: the upstream pool's produced-into attribution. (Simpler v1: always show "—" for Vasc-+/- on non-conduit nodes, and only show Vasc-+/- on the owning stem/root.) |
| Vasc-- | Same treatment as Vasc-+ |
| Total-lvl | `Local-lvl + (Vasc-lvl if not "—" else 0)` |
| Total-+ | `Local-+ + (Vasc-+ if not "—" else 0)` |
| Total-- | `Local-- + (Vasc-- if not "—" else 0)` |

For signaling chems (Auxin, GA, Ethylene, Stress) all Vasc cells render as `—`. Total then equals Local.

Formatting reuses existing `fmt_mass`, `fmt_au`, `fmt_vol` helpers.

### 2.3 Activity section

Per activity, show satisfaction % and any qualifier info. Existing panel code already computes most of these; reuse and expand.

- **Stem**: `Thicken` (bias, sugar-gf, stress-boost multiplier), `Elongate` (sugar-gf, water-gf, rate, GA/Eth/Auxin/Stress multipliers) — same as today.
- **Root**: `Absorb` (area, gradient, fill) + `Thicken` + `Elongate` — same as today.
- **Leaf**: `Photo` (stomatal, light, angle), `Transpire` (rate), `Expand` (sugar-gf, water-gf, auxin-boost), `Carbon deficit`, `Senescence` — same as today.
- **Active meristem**: `Growth` (total %, Sugar/Cyt/Water breakdown) — same as today.
- **Dormant meristem**: `Activation readiness` (overall READY/blocked flag + Auxin/Cyt/Sugar % bars) — same as today.

No formula changes; any inaccuracies in the current panel carry forward. If they need correction, that's a separate follow-up driven by what the new counters reveal.

### 2.4 Navigation section

Unchanged: Parent button if present, then one Child button per child, each showing type. Hover highlights the target node. Click selects it.

## Section 3 — Overlay modes

### 3.1 Selector UI

Replace the flat dropdown in `app_realtime.cpp` (search: the long `if (ImGui::Selectable(...))` chain near line 571) with:

```
[ Radio: Default | Type | Light | Level | Capacity | Growth | Activation | Starvation ]
(sub-picker row, appears only for Level and Capacity)
```

For `Level`: `Chemical: [ Sugar ▾ ]   Scope: [ Local | Vasculature ]`
For `Capacity`: `Chemical: [ Sugar ▾ ]`

### 3.2 Value function per mode

| Mode | Value function per node | Scope |
|---|---|---|
| Default | — (genome colors) | all |
| Type | categorical | all |
| Light | `leaf.light_exposure` | leaves; others gray |
| Level (chem, Local) | `node.local().chemical(chem) / chem_cap(node, chem, genome)` | all |
| Level (chem, Vasc) | `vascular_scope(node, chem) ? pool->chemical(chem) / pool_cap : NaN` | conduit users; others gray |
| Capacity (chem) | `util(node, chem)` as defined in 1.3 | nodes with a cap for that pathway; others gray |
| Growth | see 3.3 | all node types |
| Activation | see 3.4 | stems/roots with a dormant child meristem; others gray |
| Starvation | see 3.5 | all nodes |

`chem_cap(node, chem, genome)` is the existing per-chemical capacity helper or a new one if that chemical lacks one. For signaling chems lacking an explicit cap, use a fixed reference value (e.g. `1.0 AU` for auxin/CK/GA) so the overlay gradient has a meaningful range.

### 3.3 Growth overlay

Per node type, compute `actual_rate_now / genome_max_rate`:

- **Stem**: if `age < internode_maturation_ticks`: `eff_elongation_rate(node) / internode_elongation_rate`. Else: `thicken_rate_now / cambium_responsiveness` (bias-capped at 1.0).
- **Root**: same formula with root genome params.
- **Leaf**: `actual_expansion_rate_this_tick / (leaf_growth_rate × (1 + leaf_auxin_max_boost))`.
- **Active meristem**: `growth_fraction` (already computed — product of sugar-gf × [cyt-gf] × water-gf).
- **Dormant meristem**: 0 (rendered gray, not dark blue — to distinguish "off" from "0%").

### 3.4 Activation overlay

Color *parent stems/roots of dormant meristems*, not the meristems themselves (too small to see).

For each dormant SA or RA, compute its readiness = `min(auxin_%, cyt_%, sugar_%)` using the same formulas the panel shows today. Propagate that value to the parent stem/root (if parent has multiple dormant-meristem children, take the max readiness). Nodes with no dormant-meristem child render gray.

### 3.5 Starvation overlay

Four named color stops with two gradients:

```
red     -- starving gradient --     orange  |  yellow  -- sugar-coverage gradient --  green
 │                                    │     │                                          │
 starvation_ticks / max_starvation    0     │    0%                                  100%
```

Decision logic per node:
- If `starvation_ticks > 0`: interpolate between **red** (`t = starvation_ticks / max_starvation_ticks = 1`) and **orange** (`t = 0`). Color = lerp(orange, red, t).
- Else if `sugar < maintenance_cost_this_tick`:
  - `coverage = sugar / maintenance_cost_this_tick` ∈ `[0, 1)`
  - Color = lerp(yellow, green, coverage). At `coverage = 0`: pure yellow. At `coverage → 1`: green.
- Else: green.

`maintenance_cost_this_tick` is available via `tick_sugar_maintenance` once it's incurred. For the overlay we compute it fresh as `pay_maintenance_preview(node, genome, world)` — a helper that returns what maintenance would cost without actually deducting it. (Alternative: sample the previous tick's value, which is what's in `tick_sugar_maintenance` after Phase 1 has run. For rendering this is fine and cheaper.)

The hard boundary between orange and yellow is intentional: those states (currently-starving vs. about-to-starve) represent qualitatively different plant states.

### 3.6 Implementation plumbing

Replace the `enum class Overlay { ... }` with:

```cpp
enum class OverlayCategory { Default, Type, Light, Level, Capacity, Growth, Activation, Starvation };
static OverlayCategory g_overlay_category = OverlayCategory::Default;
static ChemicalID g_overlay_chem = ChemicalID::Sugar;
enum class OverlayScope { Local, Vasc };
static OverlayScope g_overlay_scope = OverlayScope::Local;
```

A single `rebuild_color_accessor()` function reads these three and produces the appropriate `ChemicalAccessor` lambda for `renderer.set_color_mode()`. Called whenever any selector changes.

Nodes returning `NaN` from the accessor are rendered neutral gray by the renderer; the renderer needs a one-line change to detect `std::isnan` and use a gray color path if it doesn't already.

## Dependencies and risks

- **Depends on `ChemicalID::Count` existing.** If it doesn't, add it (trivial).
- **TransportPool public API.** UI needs to call `pool->chemical(id)`. That accessor already exists internally (`vascular_sub_stepped.cpp` calls it), so exposing it from `TransportPool` is just a header visibility change.
- **`nearest_phloem_upstream()` / `nearest_xylem_upstream()` symbol visibility.** Same — exists internally, move to a public header.
- **Per-edge maps on every node.** 7 chems × 2 maps × ~2 entries ≈ 28 entries per node avg. On a 4000-node tree, ~112k `std::unordered_map` entries. Measure before and after; if hot, swap `unordered_map<Node*, float>` for a small vector of `{child_ptr, flux}` pairs.
- **Maintenance-cost preview.** Adding a no-op preview path to `pay_maintenance` risks divergence from the real maintenance calc. Option: extract the cost computation into a `compute_maintenance_cost()` free function that both `pay_maintenance` and the overlay call, ensuring one source of truth.
- **Panel density.** The 7 × 10 chemical table is dense. If it feels cramped at 320-px panel width, allow the panel to widen to ~480 px or switch to Alternative C (levels-only table with produced/consumed in a side block) — but start with the full grid and confirm visually first.

## Testing plan

1. `tests/test_tick_counters.cpp` (new) — mass balance per chemical on a hand-built tree.
2. Existing tests (`test_sugar_economy`, `test_cytokinin_transport`, `test_meristem`) remain green, and each gets one new assertion confirming the relevant new counter fires.
3. Realtime smoke test: launch `./build/botany_realtime`, spawn tree, click leaf, verify photosynthesis/transpiration values are nonzero and match genome constants × light.
4. Overlay smoke test: switch through all 8 categories, confirm no crashes, confirm Default-mode colors match pre-refactor.

## Rollout

Land in this order:

1. **Engine instrumentation** (Section 1) — add counters, public API helpers, reset discipline, tests. No UI changes yet.
2. **Panel content** (Section 2) — rewrite the panel body to consume new counters. Keep existing overlay selector as-is in this step.
3. **Overlay modes** (Section 3) — new two-tier selector, new accessor builder, NaN-aware renderer gray path, mode value functions.
4. **Polish + density pass** — if the table is cramped, apply the Alternative C fallback or widen the panel.

Each step is independently shippable. Steps 1 and 2 together are already a clear improvement.
