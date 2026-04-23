# Botany Simulator

A plant growth simulator using hormone-driven meristem mechanics. Plants grow from a seed node that spawns one shoot apical meristem (upward) and one root apical meristem (downward). Growth is governed by auxin/cytokinin signaling.

## Build & Run

```bash
# Build (cmake is at /usr/local/bin/cmake)
# IMPORTANT: always rebuild before running tests — LSP diagnostics can be stale
/usr/local/bin/cmake --build build

# Run realtime viewer (must run from project root for shader path)
./build/botany_realtime [--color auxin|cytokinin|sugar|gibberellin|ethylene|type]

# Load a previously saved snapshot and continue simulating it
./build/botany_realtime --load-plant saves/plant_<timestamp>_tick<N>.tree [--genome-override some_genome.txt]

# Save the current plant from inside botany_realtime
Cmd+S (macOS) / Ctrl+S (Linux/Windows) writes saves/plant_<timestamp>_tick<N>.tree

# Run evolution app
./build/botany_evolve

# Run tests
./build/botany_tests

# Run a headless sim with sparse recording + chain-profile dump
# (see Diagnostics section for what these produce)
./build/botany_headless 10000 /tmp/rec.bin \
    --recording-interval 500 \
    --chain-profile /tmp/profile.csv \
    [--debug-log /tmp/debug.csv] \
    [--economy-log /tmp/economy.csv]

# External dependency: evolve library at /Users/wldarden/repos/evolve
# (linked via add_subdirectory in CMakeLists.txt)
```

## Architecture

### Engine (`src/engine/`)
- **genome.h** - `Genome` struct with all tunable parameters + `default_genome()`
- **node/** - Node subfolder:
  - `node.h/cpp` — `Node` base class (position, radius, parent/children, chemicals) with non-virtual `tick()` and virtual `update_tissue()`. `NodeType` enum: `STEM, ROOT, LEAF, APICAL, ROOT_APICAL`. Downcasting via `as_stem()`, `as_root()`, `as_leaf()`, `as_apical()`, `as_root_apical()`. `is_meristem()` returns true for APICAL and ROOT_APICAL.
  - `stem_node.h/cpp` — `StemNode` (secondary thickening via `auxin_flow_bias`, owns `phloem` and `xylem` TransportPools, intercalary elongation)
  - `root_node.h/cpp` — `RootNode` (same with root params, gradient-based water absorption; owns `phloem` and `xylem` TransportPools)
  - `tissues/leaf.h/cpp` — `LeafNode` (owns `leaf_size`, `light_exposure`, `senescence_ticks`; photosynthesis gated by stomatal conductance, phototropism, expansion)
  - `tissues/apical.h/cpp` — `ApicalNode` (chain growth via self-reparenting, phyllotaxis, auxin production)
  - `tissues/root_apical.h/cpp` — `RootApicalNode` (gravitropism, root chain growth, cytokinin production from local metabolic state)
  - `meristems/helpers.h` — shared growth helpers (`turgor_fraction`, `auxin_growth_factor`, `sugar_growth_fraction`, etc.)
- **gibberellin.h/cpp** - `compute_gibberellin()` — local GA production by young leaves (reset each tick)
- **ethylene.h/cpp** - `compute_ethylene()` + `process_abscission()` — spatial gas diffusion, leaf abscission (reset each tick)
- **vascular_sub_stepped.h/cpp** - `vascular_sub_stepped()` — N-sub-stepped vascular loop that runs after `tick_tree()`. Each sub-step: inject (sources push into conduits) → longitudinal Jacobi (pressure equalization) → radial flow (local_env ⇄ phloem/xylem in stems/roots) → extract (sinks pull from conduits). Leaves and meristems access conduits via `nearest_phloem_upstream()` / `nearest_xylem_upstream()` walk-up helpers.
- **pin_transport.h/cpp** - `pin_transport()` for polar auxin transport (shoot post-order → seed_collected → root pre-order), plus `diffuse_auxin_across_seed_junction()` — an explicit shoot↔root auxin bridge that runs BEFORE `pin_transport()`. The bridge does a Fickian diffusion between every (shoot-child, root-child) pair of the seed, bypassing `seed.local`. It's a safety net so auxin from the shoot always reaches the root system even if PIN phase ordering misbehaves.
- **sugar.h/cpp** - `sugar_cap()`, `transport_sugar()` helpers. Diffusion itself lives in `Node::transport_chemicals()`.
- **hormone.h/cpp** - Empty placeholder (auxin/cytokinin transport moved to `Node::transport_chemicals()`)
- **world_params.h** - `WorldParams` struct (light level, construction costs, sugar_production_rate, maintenance rates, soil_moisture 0-1) — non-genetic physical constants
- **plant.h/cpp** - `Plant` class owns all nodes, hard-caps root meristems at 10 000. `Plant::tick_tree()` orchestrates: per-node metabolism DFS walk, then `vascular_sub_stepped()` for vascular transport (runs after metabolism — 1-tick lag). `queue_removal()` / `flush_removals()` for deferred node cleanup.
- **engine.h/cpp** - `Engine` iterates plants and calls `plant.tick(world_params)`

### Node Class Hierarchy
```
Node (base — position, local_env compartment, non-virtual tick, virtual update_tissue)
├── StemNode   (thicken via auxin_flow_bias, elongate; owns phloem + xylem TransportPools)
├── RootNode   (same with root params, gradient water absorption; owns phloem + xylem TransportPools)
├── LeafNode   (leaf_size, light_exposure, photosynthesize, expand, phototropism; local_env only)
├── ApicalNode (chain growth, phyllotaxis, auxin production; local_env only)
└── RootApicalNode (gravitropism, chain growth, cytokinin production from local metabolic state; local_env only)
```

All 5 node types extend `Node` directly — flat hierarchy, no intermediate classes. Meristems are real nodes in the tree graph, participating in chemical transport naturally. No axillary node types.

### Renderer (`src/renderer/`)
- OpenGL 4.1 core profile, GLFW window, orbit camera
- Draws cylinders between parent-child nodes, leaves as quads
- Color modes: default (brown stem/root), chemical heatmap (auxin/cytokinin/sugar/gibberellin/ethylene), type (green=shoot, orange=root, red/blue=meristems)
- **GPU shadow system** — deep shadow maps (opacity shadow maps) with power-curve slice distribution, 5 samples per leaf for soft penumbras, arbitrary sun angle via ImGui sliders, ground shadow plane rendering
- **Dormant-bud render skip** — `Node::ever_active` is a sticky flag (default `true` on non-meristems; `false` on freshly spawned lateral SA/RA buds; flips to `true` on first `activate()` and stays). Renderer `for_each_node` loops early-return for meristems where `is_meristem() && !ever_active`. Mature trees have thousands of permanently-dormant lateral buds — skipping them removes that many draw calls. Sim semantics are unchanged; only the render pass skips.
- **Selected-node info panel** — chemicals table has columns `Self | Transporting | Max Trans.` per chemical. `Transporting` shows last-tick flux where tracked (auxin via `last_auxin_flux`, sugar via `tick_sugar_transport`); CK/water/GA show `-` (no per-tick flux counters yet). `Max Trans.` is `base_transport + radius_factor × transport_scale` summed across edges — **note:** this is the diffusion formula, which is wrong for vascular chems on conduit edges (real Jacobi throughput is `xylem_conductance × cross-section`). Known issue, see Known Issues section.

### Evolution (`src/evolution/`)
- **genome_bridge.h/cpp** — `build_genome_template()`, `to_structured()`, `from_structured()` convert between `botany::Genome` and `evolve::StructuredGenome`. Each genome field is a named gene with per-gene mutation config (strength = % of valid range). 8 linkage groups for crossover (auxin, cytokinin, shoot_growth, root_growth, geometry, sugar_economy, gibberellin, ethylene).
- **fitness.h/cpp** — `PlantStats` (8 objectives: survival, biomass, leaves, sugar, height, crown_ratio, branch_depth, leaf_height_spread). `evaluate_plant()` runs a headless sim. `evaluate_group()` runs multiple competing plants in one shared Engine. `compute_fitness()` does per-generation normalized weighted scoring.
- **evolution_runner.h/cpp** — `EvolutionRunner` manages the GA lifecycle: init population from default genome with mutations, evaluate in parallel (std::thread), normalize + score fitness, tournament selection + elitism + crossover + mutation. Supports competition mode (competitors > 1 = shared sim with light competition). Tracks best genome + its competitor group for display. Autosaves genome on fitness improvement.

### Apps
- **app_realtime.cpp** - Interactive viewer with pause/speed controls, production/maintenance sugar stats
- **app_headless.cpp** - Headless precompute, saves binary recording
- **app_playback.cpp** - Playback viewer with ImGui scrubbing
- **app_evolve.cpp** - Evolution app with ImGui config (population, competitors, max ticks, threads, 8 fitness weight sliders), progress bar, live stats table, fitness history plot, best plant + competitors rendering (best at full color, competitors dimmed), Export Best button. Autosaves best_genome.txt on improvement.
- **app_sugar_test.cpp** - Headless sugar economy tester. Builds 3 hardcoded static trees (seedling/medium/large), freezes growth, runs N ticks of production/maintenance/transport. Reports production/maintenance ratios. Usage: `./build/botany_sugar_test [--ticks N] [--csv] [--tree seedling|medium|large]`

Snapshot format: `src/serialization/plant_snapshot.{h,cpp}` writes `.tree`
binary files (magic "BTNT", single plant, full node + conduit + meristem
state). Playback recordings (`botany_headless --recording-interval N` → `.bin`)
are a separate, playback-only format owned by `serializer.*`.

### Tests (`tests/`)
218 tests / 583 assertions covering: node, plant, sugar, water, hormone, gibberellin, ethylene, meristem, engine, serializer, evolution, auxin sensitivity, vascularization, and cytokinin transport.
Key files:
- `tests/test_vascularization.cpp` — 6 integration tests: no-bias no-thicken, bias-proportional growth rate, canalization ratchet, conductance-weighted phloem distribution.
- `tests/test_cytokinin_transport.cpp` — production, transport through xylem, elongation gate, dormant-meristem lifecycle.
- `tests/test_meristem.cpp` — dormant activation paths (apical dominance, xylem CK sensing, `ever_active` flag lifecycle).

### Diagnostics / validation tooling (`app_headless.cpp`)
- `--chain-profile <path.csv>` — at end of sim, walks from seed along the longest shoot and root chain and dumps a CSV with columns: `chain_side, depth, node_id, type, y, radius, auxin, cytokinin, sugar, water, xyl_cytokinin, xyl_water, phloem_sugar`. Consumer for validation scripts that want trajectory and final-state hormone gradients.
- `--recording-interval N` — save the full tick state every Nth tick instead of every tick. Default 1 (backward-compatible, every tick). `N=0` disables recording entirely. For long validation sims a value like 500 keeps the `.bin` under a few MB and gives ~20 snapshots viewable via `botany_dump` / `botany_playback`.
- `--debug-log <path.csv>` — per-tick per-node full state CSV (large; ~1 GB for 10k ticks on mature plants).
- `--economy-log <path.csv>` — per-tick plant-wide chemical mass balance (production, consumption, storage).

### Planning docs (`docs/long-term-plan/`)
Milestone folders with design + review documents. Current vascularization plan and code-review docs live at `docs/long-term-plan/milestone-2/vascularization/`.

## Chemical Transport Model

The plant uses a **compartmented dual transport system** matching real plant biology.

**Three compartments per node class:**
- `local_env` — every node owns one. Holds the node's own parenchyma chemicals: sugar, water, auxin, cytokinin, gibberellin, stress (anything the node itself uses for metabolism, growth, or signaling).
- `phloem` — `StemNode` and `RootNode` only. Sieve-tube pool, carries sugar longitudinally.
- `xylem` — `StemNode` and `RootNode` only. Vessel pool, carries water and cytokinin longitudinally.

Specialty nodes (`LeafNode`, `ApicalNode`, `RootApicalNode`) have only `local_env`. They interact with their nearest ancestor's phloem/xylem via `nearest_phloem_upstream()` and `nearest_xylem_upstream()` walk-up helpers.

All chemical access goes through explicit compartment accessors:
- `node.local().chemical(id)` for the local compartment
- `stem.phloem()->chemical(id)` and `stem.xylem()->chemical(id)` for transport pools

**Two transport pathways:**

1. **Sub-stepped vascular bulk flow** (`vascular_sub_stepped.cpp`) — N iterations per tick in this order:
   1. **Inject** — sources transfer `budget/N` into their walk-up parent's conduit.  Sources: leaves push sugar into parent phloem; seed pushes its reserve sugar into its own phloem; root nodes and root apicals push water into their own / parent's xylem (root pressure); root apicals also push their produced cytokinin into parent xylem.
   2. **Longitudinal Jacobi** — one pass of neighbor pressure equalization across every (parent, child) edge that has matching conduit pools. Pushes CK through xylem chain from root to shoot; pushes sugar through phloem; pushes water through xylem.
   3. **Radial flow** — stem/root `local_env` ⇄ own `phloem`/`xylem`, rate-limited by radius-dependent `radial_permeability(r)` (see Radial Permeability section). Cytokinin is NOT subject to radial flow (vascular-only).
   4. **Extract** — sinks transfer `budget/N` from their walk-up parent's conduit.  Sinks: apical and root apical meristems pull sugar from phloem; leaves (transpiration) and meristems (turgor) pull water from xylem. Cytokinin rides along proportionally with water drawn from the xylem pool.

   Ordering changed 2026-04-20: previously `inject → radial → extract → Jacobi`, which stranded freshly-injected CK in root xylem because radial immediately absorbed it into root `local_env` before Jacobi could propagate it upward. New order lets pressure waves traverse the xylem chain before radial absorption pulls chemicals back out. Jacobi is a pure neighbor equalizer — no global awareness of sources or sinks. Pressure fields created by inject and extract drive routing automatically.

2. **Local diffusion** (`Node::transport_with_children()`) — per-node during DFS walk, handles only signaling chemicals: auxin, gibberellin, stress. Short-range cell-to-cell transport. Uses `transport_received` buffer for anti-teleportation (chemicals move at most one hop per tick).
   - **Cytokinin is owned end-to-end by the vascular system.** Local diffusion explicitly skips CK on any edge touching a vascular conduit (stem/root at or above `vascular_radius_threshold`). Without this skip, local diffusion on specialty-to-conduit edges (e.g., RA → parent ROOT) would drain the RA's freshly produced CK into the parent's `local` pool before vascular inject could move it into the xylem — where it would be stranded with no routing.
   - **Auxin is PIN-transported, not diffused.** `pin_transport()` (phase-A/B/C custom algorithm) runs separately for auxin. Local diffusion still runs on auxin for the short-range "last mile" within non-conduit neighborhoods, but the main flow is PIN.
   - **`diffuse_auxin_across_seed_junction()` runs before `pin_transport()`** as an explicit shoot↔root auxin bridge. Directly transfers a Fickian diffusion fraction from every shoot-child to every root-child of the seed, bypassing `seed.local`. Belts-and-suspenders path guaranteeing the seed junction always carries auxin across.

| Chemical | Pathway | Notes |
|----------|---------|-------|
| Auxin | PIN transport + local diffusion | `pin_transport.cpp` handles polar shoot→seed→root cascade (Phase A/B/C); `diffuse_auxin_across_seed_junction` bridges shoot↔root at the seed; local diffusion handles last-mile within non-conduit neighborhoods. Basipetal bias (-0.1) on shoot side, inverted to (+0.1) on root side so the shifted equilibrium always favors flowing toward tips. |
| Gibberellin | Local diffusion | Short-range, local to producing leaf |
| Stress | Local diffusion | Local mechanical alarm signal |
| Ethylene | Spatial gas diffusion | Global 3D pass (reset each tick) |
| Sugar | Vascular (phloem): inject + Jacobi + radial + extract | Leaves (source) and seed (reserve source) push into phloem via inject; Jacobi propagates along phloem chain; radial exchanges with stem/root `local_env`; meristems pull from phloem via extract |
| Water | Vascular (xylem): inject + Jacobi + radial + extract | Roots absorb from soil into `local_env`; root nodes and root apicals push water into xylem via inject (root pressure); Jacobi propagates; radial exchanges with stem/root `local_env`; leaves and meristems pull from xylem via extract |
| Cytokinin | Vascular (xylem): inject + Jacobi + extract (no radial) | Root apicals produce CK locally from sugar+water metabolic state; inject pushes CK into nearest upstream xylem; Jacobi propagates along xylem chain; extract moves CK into sink `local_env` proportionally with water drawn. No radial — CK stays in conduits until a sink explicitly pulls it. |

**Vasculature admission:** Inline check in `transport_with_children()` — stems/roots with `radius >= vascular_radius_threshold` (0.01 dm) have conduit pools. The seed is always a conduit junction. Leaves and meristems never have conduits; they access the nearest ancestor's pools via walk-up helpers. Since `initial_radius` (0.015 dm) is above the threshold, newly spawned internodes qualify from birth.

**Budgets are frozen at the start of each vascular pass.** `compute_budget(node, g, world)` classifies each node as source/sink for each chemical using the state at the moment vascular starts. Amortization across N sub-steps is pure division — no dynamic re-evaluation within the loop.

**Fixed N (`world.vascular_substeps`, default 25).** Each Jacobi iteration propagates pressure by ~1 hop. Plants whose longest source-to-sink chain exceeds N will show distance-dependent apical supply — distal apices get less sugar than proximal ones. This is intentional and biologically realistic (real tall plants have hydraulic limitations). Evolution can select for genomes whose radial-permeability curve and tissue geometry fit within a given N budget.

### Vascular Sources and Sinks

**Phloem (sugar):**
- Sources: leaves (surplus over `leaf_reserve_fraction_sugar × sugar_cap`), seed (mobilizing its reserve — the STEM with no parent actively pushes sugar above reserve into its own phloem, analogous to root pressure for water).
- Sinks: apical and root apical meristems (refilling to `meristem_sink_target_fraction × sugar_cap`).

**Xylem (water):**
- Sources: root nodes and root apicals (root pressure — actively push `water_local − reserve` into their own / parent's xylem each tick; models the osmotic pumping real plants use before transpiration ramps up).
- Sinks: leaves (transpiration demand, refilling to `leaf_turgor_target_fraction × water_cap`) and shoot apical meristems.

**Xylem (cytokinin):**
- Source: root apicals only (inject their full `local().Cytokinin` into the nearest upstream xylem each tick; CK production at the RA is driven by `root_cytokinin_production_rate × mf_cyto(sugar, water)` — NOT gated on local auxin post-2026-04-20).
- Sink: any node drawing water from xylem extract — CK rides along proportionally with the water drawn.

### Auxin, Cytokinin, Sugar in local_env

**Auxin** (PIN transport + local diffusion, bias -0.1 shoot / +0.1 root-inverted):
- Persistent across ticks (not reset)
- Produced by `ApicalNode` each tick (modulated by growth rate)
- Also produced by growing `LeafNode` proportional to growth rate (zero when full-size)
- Also produced at a small floor by `RootApicalNode` (`root_tip_auxin_production_rate ≈ 0.002`, self-equilibrium ~0.017 — just enough for PIN recycling and lateral root initiation, below `root_auxin_activation_threshold = 0.01` so deep laterals activate on their own floor without needing shoot-delivered auxin)
- Shoot axillary buds sense **parent stem's** auxin level to judge apical dominance; activate when it drops
- Primary flow path is PIN (`pin_transport.cpp`), not local diffusion — diffusion handles "last mile" within non-conduit neighborhoods only
- Decays by `auxin_decay_rate` per tick (0.12, half-life ~5.4h). Intentionally short-range: long-distance shoot↔root signaling is handled by CK/sugar in vascular, not auxin.

**Cytokinin** (long-distance via xylem; produced at RAs from local metabolic state):
- Persistent across ticks (not reset)
- Produced by `RootApicalNode` from local sugar+water metabolic state — **NOT** gated on local auxin.  See `docs/superpowers/specs/2026-04-20-hormone-biology-tree-scale-design.md` for rationale; the prior auxin coupling broke scaling beyond small seedlings.
- Travels up in xylem via sub-stepped Jacobi pressure flow; scale-free with plant height
- Drives SA elongation (gate Km = `cytokinin_growth_threshold`) and RA elongation (gate Km = `root_ck_growth_floor`)
- Lateral root axillary activation is suppressed by high local CK (`root_cytokinin_inhibition_threshold`) — prevents runaway root branching when the root system is already well-established
- Decays by `cytokinin_decay_rate` per tick (5%) in `local()` compartments only; xylem conduits have no decay

**Sugar** (vascular phloem; local diffusion handles last-hop within non-conduit nodes):
- Persistent across ticks
- Produced by `LeafNode` via photosynthesis (gated by stomatal conductance — see Water Model)
- Consumed by all nodes (maintenance costs in WorldParams)

## Radial Permeability (radius-dependent)

Radial flow between a stem/root's own `local_env` and its own `phloem`/`xylem` is rate-limited by a radius-dependent permeability. Young thin stems are leaky (get plenty of nutrients, grow fast). Mature thick trunks asymptote to a floor (enough to maintain themselves, not enough to siphon flow away from distal apices).

```
radial_permeability(r) = base × (floor + (1 - floor) / (1 + (r / half_radius)²))
```

Independent curves for phloem-radial (sugar) and xylem-radial (water + cytokinin), each with its own genome params: `base_radial_permeability_*`, `radial_floor_fraction_*`, `radial_half_radius_*`.

Combined with fixed N sub-steps, the curve shape is what makes tall plants viable — a mature trunk with low radial permeability acts as a hydraulic highway delivering water/sugar to the canopy rather than bleeding it into intermediate stem tissue. The `radial_half_radius` inflection point determines how quickly a stem "closes off" as it matures.

## Secondary Thickening (Cambium)

Radial growth is driven by PIN saturation history (`auxin_flow_bias`), not node age:

```
delta_radius = cambium_responsiveness × auxin_flow_bias × sugar_gf × stress_boost
```

- `auxin_flow_bias` stored on the parent, keyed by child pointer, range `[0, 1]` (see Canalization Model) — connections where auxin has flowed repeatedly have higher saturation and thicker cambium
- `sugar_gf` = `sugar / max_cost`, capped at 1 — sugar-limited when starved
- `stress_boost = 1 + stress × stress_thickening_boost` — mechanical load accelerates thickening
- Early return if `auxin_flow_bias < 1e-6` — zero-flux connections never thicken (monocot behavior)
- **No age gate** — `cambium_maturation_ticks` is gone; the bias threshold is the gate

The self-reinforcing loop: main axis carries most flux → highest bias → fastest thickening → widest pipe → more flow → more bias. Lateral branches with weak canalization stay thin automatically.

## Canalization Model

Auxin flow through parent-child connections builds transport preference over time. Single-layer PIN saturation model, stored on the parent node, keyed by child pointer.

**`auxin_flow_bias`** — represents PIN-transporter saturation on the parent→child connection, in `[0, 1]`.
- Each tick, `saturation = flux / (child.radius² × pin_capacity_per_area)`, clamped to `[0, 1]`
- Bias chases saturation exponentially: `flow_bias += (saturation - flow_bias) × smoothing_rate`
- Natural decay when flux stops — no separate decay param; bias drifts back toward 0 as saturation drops to 0.
- Persistent across ticks (not reset).

**Effect on transport:** Bias multiplies the conductance weight for each child connection. `weight = pipe_cap × (1 + canalization_weight × auxin_flow_bias)`. Redistributes flow (all chemicals, not just auxin) among siblings proportionally — does **not** amplify total flow.

**Effect on cambium:** Drives secondary thickening (see Secondary Thickening section).

**Lifecycle:** Bias transfers on `replace_child` (chain growth preserves branch history). New children start at 0. Cleaned up on `die()`.

## Water Model

Water is a **persistent** capacity-based resource (like sugar):

- **Absorbed** by `RootNode` and `RootApicalNode` proportional to surface area and `soil_moisture`. Absorption is gradient-based: `desired = (soil_moisture - water_conc) × absorption_rate × surface_area` — self-limiting when the root is full (concentration approaches soil moisture).
- **Lost** by leaves via transpiration (proportional to leaf area and light exposure) and photosynthesis water cost
- **Stomatal conductance** — `stomatal = clamp(water / water_cap, 0.2, 1.0)` scales photosynthesis output. Water deficit partially closes stomata, reducing sugar production. (Minimum 20% conductance — stressed leaves still photosynthesize slowly.)
- **Transported** via vascular xylem bulk flow + local diffusion with upward bias (`water_bias = 0.05`)

Surface area: `2π r L` for root segments, `2π r²` (hemisphere) for root apical tips.

## Gibberellin Model

GA is **reset to zero every tick** (signal model):

- Produced by young LEAF nodes only (`leaf_age < ga_leaf_age_max`)
- Applied locally — `ga_production_rate × leaf_size` added to the leaf's parent and grandparent stem nodes
- **Effects:** boosts internode elongation rate (`× ga_elongation_sensitivity`) and max internode length (`× ga_length_sensitivity`)

## Ethylene Model

Ethylene is **reset to zero every tick** (signal model). Four production triggers per leaf:

1. **Sugar starvation** — leaf sugar below maintenance threshold
2. **Shade** — leaf `light_exposure < ethylene_shade_threshold`
3. **Old age** — leaf age exceeds species maximum
4. **Crowding** — local node density above threshold

**Spatial diffusion** — spreads as a gas through 3D space (NOT via the tree graph). Each emitting node contributes to all nodes within `diffusion_radius`, attenuated by distance.

**Effects:**
- **Leaf abscission** — ethylene > threshold → senescence flag set → leaf yellows → removed after `senescence_duration` ticks
- **Elongation inhibition** — high ethylene suppresses internode elongation in nearby stem nodes

## Tick Control Flow

`Plant::tick_tree()` runs three phases in order:

**Phase 0 — Polar auxin transport** (before per-node metabolism):
1. `diffuse_auxin_across_seed_junction()` — explicit shoot↔root auxin bridge.  Takes a Fickian diffusion fraction from every shoot-child to every root-child of the seed, bypassing `seed.local`.  Runs FIRST so it can move auxin while it's still at the shoot children's local pools.
2. `pin_transport()` — PIN-mediated polar auxin transport. Phase A: shoot post-order, each shoot node pumps its local auxin into the parent's `transport_received` buffer (or `seed_collected` local accumulator for direct shoot children of seed).  Phase B: seed distributes `seed_collected` to root children of seed by radius share.  Phase C: root pre-order, each root node receives its forwarded auxin and passes the remainder deeper, accumulating any overflow in its `local()`.

**Phase 1 — Per-node metabolism** (`tick_recursive()`): DFS walk from seed. Parent ticks before children. Each `Node::tick()`:
1. `age++`
2. `sync_world_position()` — recompute world position from offset chain
3. `pay_maintenance()` — deduct sugar from `local_env`; die and return if starvation ticks exceed max
4. `update_tissue(plant, world)` — **virtual** type-specific growth (see below); produces into `local_env`, consumes from `local_env`
5. `sync_world_position()` — re-sync after tissue may have changed offset
6. `update_physics()` — compute mass, stress; apply droop and break if overstressed
7. `transport_chemicals()` → `transport_with_children()` (local diffusion: auxin, gibberellin, stress between neighboring `local_env`s) + `update_canalization()` + `decay_chemicals()`
8. Flush `transport_received` buffer into `chemical()`

Children are ticked recursively after step 8 using a snapshot of the children list (safe for reparenting).

**Phase 2 — Vascular transport** (`vascular_sub_stepped()`): N-sub-stepped loop of inject → Jacobi → radial → extract — see Chemical Transport Model above.

Ordering is **metabolism first, vascular after**. This gives a 1-tick lag between a leaf producing sugar (into `local_env` during Phase 1 of tick N) and that sugar being loaded into phloem (during Phase 2 of tick N, which reads the just-produced budget). At 1-hour tick granularity this is biologically defensible — real plants buffer an hour of sugar in mesophyll cells easily.

**`update_tissue` convention** — subclass-specific growth happens entirely inside this virtual method. Recommended call order within `update_tissue`: produce → consume → evaluate → grow → orient. Growth verbs: `elongate` (extend offset), `thicken` (increase radius), `expand` (increase leaf_size). Meristem chain growth and axillary activation also happen inside `update_tissue` via `Plant::queue_removal()`.

## Key Design Decisions
- **Meristems are nodes** — 2 meristem types (`ApicalNode`, `RootApicalNode`) extend `Node` directly. Real children in the tree graph, participate in chemical transport naturally. No axillary node types.
- **Compartmented chemical model** — every node has a `local_env`; stems/roots additionally have typed `phloem` and `xylem` TransportPools. Specialty nodes use walk-up helpers (`nearest_phloem_upstream()` / `nearest_xylem_upstream()`) to reach their nearest ancestor's conduit.
- **Sub-stepped vascular with fixed N** — N iterations per tick of inject → Jacobi → radial → extract. Budgets computed once before the loop and amortized evenly. Fixed N produces realistic hydraulic limitation on tall plants — distal apices get less supply when chains exceed N.
- **Radial permeability asymmetry** — young thin stems are leaky (lots of radial exchange), mature trunks asymptote toward a low floor. This is what makes the vascular system act like a highway: old trunk tissue stops bleeding sugar/water so it reaches the canopy.
- **Tick-then-vascular ordering** — per-node metabolism runs first (Phase 1), vascular transport after (Phase 2). The 1-tick lag between photosynthesis and phloem loading is biologically defensible at 1-hour granularity.
- **Vascular admission is radius-gated inline** — stems/roots with `radius >= vascular_radius_threshold` (0.01 dm) have conduit pools. No age gate. Since `initial_radius` (0.015 dm) is above the threshold, new internodes join the network immediately. Inline check in `transport_with_children()`, not a separate function.
- **Local diffusion rates are low by design** — auxin/GA/stress use low diffusion rates (0.05–0.10) so they remain cell-to-cell signals, not whole-plant gradients. Radius bottleneck handles the rest.
- **Thickening is bias-driven** — `cambium_responsiveness × auxin_flow_bias` replaces old fixed `thickening_rate`. Zero-flux branches never thicken.
- **Canalization redistributes, does not amplify** — `canalization_weight × auxin_flow_bias` shifts proportional weight among sibling edges. Does not increase total flow budget.
- **Flat node hierarchy** — `Node` base class with 5 direct subclasses. Each owns type-specific fields and behavior via `virtual update_tissue()`. Downcasting via `as_stem()`/`as_root()`/`as_leaf()`/`as_apical()`/`as_root_apical()` — fast `static_cast` gated on `NodeType` enum, no RTTI.
- **Root elongation is gated by local cytokinin (not auxin)** — RA `elongate()` uses `growth_fraction(sugar, max_cost, local_cytokinin, root_ck_growth_floor)`. CK at the RA is produced locally from sugar+water metabolic state, so the gate is scale-free: deep RAs in a 10m tree self-gate on their own CK, no dependency on distant-shoot signals.
- **Dormant meristems cost zero maintenance** — only active meristems pay `sugar_maintenance_meristem`
- **Dormant meristems sense parent xylem directly** — SA and RA `can_activate()` read `max(own_local_CK, nearest_xylem_upstream_CK)` for the cytokinin threshold. Dormant buds never run vascular extract (no water demand), so their `local()` CK stays ~0; without reading upstream xylem the SA CK-permissive check was unreachable (dormant SAs never woke) and the RA CK-inhibition check was dead code (lateral RAs over-proliferated in CK-rich regions). Flow is unchanged — this is read-only sensing, analogous to how leaves access parent-stem vasculature.
- **Sugar economy: physics in WorldParams, strategy in Genome** — production rate, maintenance costs, and construction costs are physical constants in WorldParams (not evolvable). Evolution optimizes plant shape/structure, not metabolism.
- **Hormone scale separation** — auxin is a short-range signal (decay 12%/tick, one hop per tick via PIN transport, effective range ~5–20 internodes), used for apical dominance, canalization, and local tropism.  Cytokinin is the long-distance root→shoot signal, produced at root tips from local metabolic state (sugar + water, independent of shoot auxin) and carried by xylem bulk flow.  Sugar is the long-distance shoot→root channel via phloem bulk flow.  **No mechanism depends on auxin reaching a distant node** — this is the fundamental scaling invariant that lets the model support plants from 1 dm seedlings to 10 m trees.  See spec/plan docs for the 2026-04-20 hormone-biology-tree-scale refactor.

## Hormone Architecture Validation (2026-04-20)

The hormone-biology-tree-scale plan was validated at two scales:

- **Seedling (1000 ticks):** 688 nodes, primary SA cytokinin = 0.126 (125× the `cytokinin_growth_threshold`), primary RA depth 42 internodes, root xylem CK sum = 54. Visible auxin gradient along shoot chain, visible CK gradient in root xylem.
- **Adolescent/mature (10000 ticks):** 3748 nodes, primary SA at 5.77 m height, cytokinin = 0.105, auxin = 0.016 (active producer), deepest RA at 3.02 m with active CK production (0.002), seed.xylem CK = 7.62 (760× the > 0.01 spec threshold), seed radius 7.97 cm (mature trunk). Lateral break pattern correct: 256 STEM + 1 SA on primary shoot chain (tight apical dominance). Growth self-regulates via sugar/CK feedback — equilibrium plateaus at ~3200 nodes with periodic growth bursts.

Checkpoints are reproducible with: `./build/botany_headless <ticks> /tmp/out.bin --recording-interval 500 --chain-profile /tmp/profile.csv`.

## Tuning Parameters (genome.h)
- `auxin_threshold` (0.15) - lower = fewer shoot branches, higher = more
- `cytokinin_threshold` (0.15) - lower = fewer root branches, higher = more
- `auxin_bias` (-0.1) - shifted equilibrium for basipetal flow (negative = toward root)
- `stem_auxin_max_boost` (0.5) - max elongation promotion from auxin (saturating)
- `stem_auxin_half_saturation` (0.2) - auxin level for half-max stem effect
- `root_auxin_max_boost` (-0.3) - max elongation inhibition from auxin (negative = inhibits)
- `root_auxin_half_saturation` (0.1) - auxin level for half-max root effect (roots very sensitive)
- `leaf_auxin_max_boost` (0.3) - max leaf expansion promotion from auxin
- `leaf_auxin_half_saturation` (0.2) - auxin level for half-max leaf effect
- `apical_auxin_max_boost` (0.2) - max tip extension promotion from auxin
- `apical_auxin_half_saturation` (0.3) - auxin level for half-max apical effect
- `root_apical_auxin_max_boost` (-0.2) - max root tip extension inhibition from auxin
- `root_apical_auxin_half_saturation` (0.1) - auxin level for half-max root apical effect
- `cytokinin_bias` (0.1) - shifted equilibrium for acropetal flow (positive = toward tips); negated for root-type children (see transport model)
- `root_cytokinin_production_rate` (0.5) - baseline cytokinin produced per tick at active RAs, multiplied by a sugar+water metabolic factor (NOT by local auxin, post 2026-04-20 decoupling); drives whole-plant CK supply via xylem bulk flow
- `root_auxin_activation_threshold` (0.01) - minimum local auxin to activate a dormant root meristem; below RA self-eq (~0.017) so deep lateral RAs can activate on their own auxin floor without needing distant shoot auxin
- `root_cytokinin_inhibition_threshold` (0.15) - cytokinin above this inhibits new root meristem activation
- `root_tip_auxin_production_rate` (0.002) - auxin produced per tick at active root tips (PIN recycling, lateral root initiation); self-equilibrium ~0.017 with default decay; below activation threshold so it doesn't spuriously activate laterals, but sufficient for canalization
- `root_ck_growth_floor` (0.001) - Km for cytokinin-gated RA elongation; RA grows at sugar-rate-limited pace when local CK ≥ this floor, stops when CK drops to zero
- `auxin_diffusion_rate` (0.05) - slow polar transport (cell-to-cell only)
- `cytokinin_diffusion_rate` (0.1) - moderate; bulk xylem flow handles the rest
- `hormone_base_transport` (0.5) - throughput floor for hormones (even thin tips can signal)
- `hormone_transport_scale` (1.0) - how much radius amplifies hormone throughput
- `sugar_base_transport` (0.1) - throughput floor for sugar
- `sugar_transport_scale` (5.0) - how much radius amplifies sugar throughput
- `apical_auxin_baseline` (0.15) - base auxin output of shoot apical meristem per tick
- `apical_growth_auxin_multiplier` (2.0) - growth bonus: total = baseline × (1 + multiplier × growth_fraction)
- `leaf_auxin_baseline` (0.15) - scaling constant for leaf auxin production
- `leaf_growth_auxin_multiplier` (0.1) - fraction of leaf_auxin_baseline at max growth (10% of apical baseline)
- `cambium_responsiveness` (0.00002 dm/hr·bias) - thickening rate per unit auxin_flow_bias (PIN saturation in [0,1]). A fully-saturated trunk edge (bias ~1.0) thickens at ~0.00002 dm/hr; laterals with lower saturation thicken proportionally less.
- `branch_angle` (0.785 rad / 45°) - angle of shoot branches from parent stem
- `root_branch_angle` (0.35 rad / 20°) - angle of root branches
- `ga_production_rate` (0.5) - GA per dm leaf_size per tick from young leaves
- `ga_leaf_age_max` (168) - only leaves younger than 7 days produce GA
- `ga_elongation_sensitivity` (2.0) - GA boost to elongation rate
- `ga_length_sensitivity` (1.5) - GA boost to max internode length
- `ethylene_abscission_threshold` (0.5) - triggers leaf senescence
- `ethylene_shade_threshold` (0.3) - light_exposure below which shade-ethylene kicks in
- `senescence_duration` (48) - ticks from senescence to leaf drop
- `smoothing_rate` (0.1) - exponential lerp rate for auxin_flow_bias toward current PIN saturation (~20 tick response; natural decay when flux stops)
- `pin_capacity_per_area` (500.0 AU/(dm²·tick)) - max auxin transport per unit cross-section at full efficiency; also the denominator in `saturation = flux / (r² × pin_capacity_per_area)`
- `pin_base_efficiency` ([0,1]) - cold-start PIN efficiency when auxin_flow_bias = 0 (constitutively active PINs)
- `canalization_weight` (1.0) - global scaling on bias effect in transport weighting. 0 = canalization disabled entirely
- `vascular_radius_threshold` (0.01 dm) - minimum stem/root radius for bulk vascular admission. Below `initial_radius` (0.015 dm) so all new internodes qualify from birth.
- `water_absorption_rate` (0.5) - ml / (dm² root surface · hr) per unit soil moisture
- `transpiration_rate` (0.04) - ml / (dm² leaf area · hr) at full light
- `photosynthesis_water_ratio` (0.5) - ml water consumed per g sugar produced
- `water_storage_density_stem` (800.0) - ml / dm³ stem/root tissue
- `water_storage_density_leaf` (3.0) - ml / dm² leaf area
- `water_cap_meristem` (1.0) - ml fixed cap for meristems
- `water_diffusion_rate` (0.9) - gradient responsiveness (faster than sugar)
- `water_bias` (0.05) - slight upward equilibrium shift (transpiration pull)
- `water_base_transport` (0.2) - throughput floor (higher than sugar — xylem is open pipes)
- `water_transport_scale` (4.0) - radius scaling on throughput
- `xylem_conductance` (100.0) - xylem throughput per dm² cross-section per tick (water + cytokinin); primary long-range water mover
- `phloem_conductance` (8.0) - phloem throughput per dm² cross-section per tick (sugar)
- `phloem_reserve_fraction` (0.3) - legacy name, superseded by `leaf_reserve_fraction_sugar`; fraction of sugar_cap leaves keep for themselves
- `vascular_radius_threshold` (0.01 dm) - minimum stem/root radius for conduit pool admission. Below `initial_radius` (0.015 dm) so all new internodes qualify from birth
- `phloem_fraction` (0.05) - fraction of `π × r²` cross-section that is sieve-tube phloem
- `xylem_fraction` (0.20) - fraction of `π × r²` cross-section that is xylem vessels
- `base_radial_permeability_sugar` (1.0) - peak radial permeability for phloem (at r = 0, i.e., youngest stem)
- `radial_floor_fraction_sugar` (0.1) - asymptotic floor as `r → ∞`, as a fraction of base; mature trunks bleed ~10% of young-stem rate
- `radial_half_radius_sugar` (0.3 dm) - inflection radius for phloem radial curve; stems thicker than this are past half-max permeability
- `base_radial_permeability_water` (1.0) - peak radial permeability for xylem (at r = 0)
- `radial_floor_fraction_water` (0.1) - asymptotic floor for xylem radial curve
- `radial_half_radius_water` (0.3 dm) - inflection radius for xylem radial curve
- `leaf_reserve_fraction_sugar` (0.3) - leaf keeps this fraction of sugar_cap as structural glucose reserve before loading into phloem
- `meristem_sink_target_fraction` (0.05) - meristem refills from phloem up to this fraction of its sugar_cap per vascular pass
- `leaf_turgor_target_fraction` (0.7) - leaf refills from xylem up to this fraction of water_cap per vascular pass
- *(WorldParam)* `vascular_substeps` (25) - N sub-steps in vascular loop; each step propagates pressure ~1 hop; plants with chains longer than N show distance-dependent apical supply (intentional hydraulic limit)

## Known Issues / TODO

These are known gaps flagged during the 2026-04-20 hormone-biology refactor. None block current functionality but each will mislead or be incomplete until addressed.

- **UI heatmap shows `local()` only for vascular chemicals.** After the compartmented-vascular refactor, sugar/water/cytokinin live primarily in the `xylem` and `phloem` pools, not in stems'/roots' `local`. The realtime heatmap colors by `local().chemical(id)`, so stems carrying CK in their xylem show as "empty" blue when they should show the color of their actual CK throughput. Fix would be either a `local + xylem/phloem` sum for conduit nodes or dedicated `xyl_cytokinin` / `phloem_sugar` heatmap modes.
- **Selected-node panel doesn't show xylem/phloem pools.** The chemicals table shows only `local` values. Stems can have empty local CK but meaningful xylem CK, and the user has no way to see the xylem value. Needs extra rows (or columns) for `xylem` and `phloem` chemicals on conduit nodes.
- **`Transporting` / `Max Trans.` columns use diffusion formula.** Wrong for vascular chemicals on conduit edges where real throughput is Jacobi `xylem_conductance × cross-section` (water/CK) or `phloem_conductance × cross-section` (sugar). Currently overstates transport capacity for thick stems and understates for thin ones on the vascular pathway.
- **`cytokinin_production_rate` genome parameter is dead code.** Declared with the comment "cytokinin per g sugar produced by leaves", serialized, evolved — but nothing in the engine uses it (only `root_cytokinin_production_rate` is wired). The realtime selected-node panel at `app_realtime.cpp:748` displays a "predicted leaf CK production" number computed from this dead param — visually misleading. Either wire up leaf CK production (biologically real: young leaves are a secondary CK source) or remove the dead parameter.
- **No per-tick flux counters for cytokinin, water, gibberellin, ethylene.** `tick_auxin_produced`, `tick_cytokinin_produced`, and `tick_sugar_transport` exist but analogous counters for the other chemicals don't. The `Transporting` column shows `-` for these. Low priority; add per chemical as needed.
