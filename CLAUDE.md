# Botany Simulator

A plant growth simulator using hormone-driven meristem mechanics. Plants grow from a seed node that spawns one shoot apical meristem (upward) and one root apical meristem (downward). Growth is governed by auxin/cytokinin signaling.

## Build & Run

```bash
# Build (cmake is at /usr/local/bin/cmake)
# IMPORTANT: always rebuild before running tests — LSP diagnostics can be stale
/usr/local/bin/cmake --build build

# Run realtime viewer (must run from project root for shader path)
./build/botany_realtime [--color auxin|cytokinin|sugar|gibberellin|ethylene|type]

# Run evolution app
./build/botany_evolve

# Run tests
./build/botany_tests

# External dependency: evolve library at /Users/wldarden/repos/evolve
# (linked via add_subdirectory in CMakeLists.txt)
```

## Architecture

### Engine (`src/engine/`)
- **genome.h** - `Genome` struct with all tunable parameters + `default_genome()`
- **node/** - Node subfolder:
  - `node.h/cpp` — `Node` base class (position, radius, parent/children, chemicals) with non-virtual `tick()` and virtual `update_tissue()`. `NodeType` enum: `STEM, ROOT, LEAF, APICAL, ROOT_APICAL`. Downcasting via `as_stem()`, `as_root()`, `as_leaf()`, `as_apical()`, `as_root_apical()`. `is_meristem()` returns true for APICAL and ROOT_APICAL.
  - `stem_node.h/cpp` — `StemNode` (secondary thickening via structural_flow_bias, intercalary elongation)
  - `root_node.h/cpp` — `RootNode` (same with root params, gradient-based water absorption)
  - `tissues/leaf.h/cpp` — `LeafNode` (owns `leaf_size`, `light_exposure`, `senescence_ticks`; photosynthesis gated by stomatal conductance, phototropism, expansion)
  - `tissues/apical.h/cpp` — `ApicalNode` (chain growth via self-reparenting, phyllotaxis, auxin production)
  - `tissues/root_apical.h/cpp` — `RootApicalNode` (gravitropism, root chain growth, auxin-gated cytokinin production)
  - `meristems/helpers.h` — shared growth helpers (`turgor_fraction`, `auxin_growth_factor`, `sugar_growth_fraction`, etc.)
- **gibberellin.h/cpp** - `compute_gibberellin()` — local GA production by young leaves (reset each tick)
- **ethylene.h/cpp** - `compute_ethylene()` + `process_abscission()` — spatial gas diffusion, leaf abscission (reset each tick)
- **vascular.h/cpp** - `phloem_resolve()` + `xylem_resolve()` — two separate bulk-flow passes that run AFTER the DFS tick. `phloem_resolve()`: Münch pressure-driven sugar transport — (3a) leaves load sugar above phloem reserve into parent stem, (3b) BFS on conduit network with time budget and pressure-flow, (3c) meristems unload from parent stem via permeability. `xylem_resolve()`: Phase 1/Phase 2 water+cytokinin bulk flow. `has_vasculature()`: returns true for STEM/ROOT with radius >= `vascular_radius_threshold` (0.01 dm); seed always true.
- **sugar.h/cpp** - `sugar_cap()`, `transport_sugar()` helpers. Diffusion itself lives in `Node::transport_chemicals()`.
- **hormone.h/cpp** - Empty placeholder (auxin/cytokinin transport moved to `Node::transport_chemicals()`)
- **world_params.h** - `WorldParams` struct (light level, construction costs, sugar_production_rate, maintenance rates, soil_moisture 0-1) — non-genetic physical constants
- **plant.h/cpp** - `Plant` class owns all nodes, hard-caps root meristems at 10 000. `Plant::tick()` orchestrates: `pin_transport()` → `tick_tree()` (DFS from seed) → `phloem_resolve()` → `xylem_resolve()` → `flush_removals()`. `queue_removal()` / `flush_removals()` for deferred node cleanup.
- **engine.h/cpp** - `Engine` iterates plants and calls `plant.tick(world_params)`

### Node Class Hierarchy
```
Node (base — position, chemicals, non-virtual tick, virtual update_tissue)
├── StemNode   (thicken via structural_flow_bias, elongate)
├── RootNode   (same with root params, gradient water absorption)
├── LeafNode   (leaf_size, light_exposure, photosynthesize, expand, phototropism)
├── ApicalNode (chain growth, phyllotaxis, auxin production)
└── RootApicalNode (gravitropism, chain growth, auxin-gated cytokinin production)
```

All 5 node types extend `Node` directly — flat hierarchy, no intermediate classes. Meristems are real nodes in the tree graph, participating in chemical transport naturally. No axillary node types.

### Renderer (`src/renderer/`)
- OpenGL 4.1 core profile, GLFW window, orbit camera
- Draws cylinders between parent-child nodes, leaves as quads
- Color modes: default (brown stem/root), chemical heatmap (auxin/cytokinin/sugar/gibberellin/ethylene), type (green=shoot, orange=root, red/blue=meristems)
- **GPU shadow system** — deep shadow maps (opacity shadow maps) with power-curve slice distribution, 5 samples per leaf for soft penumbras, arbitrary sun angle via ImGui sliders, ground shadow plane rendering

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

### Tests (`tests/`)
168 tests / 1285 assertions covering: node, plant, sugar, water, hormone, gibberellin, ethylene, meristem, engine, serializer, evolution, auxin sensitivity, cytokinin transport, and vascularization.  
Key file: `tests/test_vascularization.cpp` — 5 integration tests: no-bias no-thicken, bias-proportional growth rate, canalization ratchet, Münch phloem distribution, overlay accessor.

### Planning docs (`docs/long-term-plan/`)
Milestone folders with design + review documents. Current vascularization plan and code-review docs live at `docs/long-term-plan/milestone-2/vascularization/`.

## Chemical Transport Model

The plant uses a **dual transport system** matching real plant biology:

1. **Vascular bulk flow** (`vascular.h/cpp`) — two passes that run AFTER the DFS tick:
   - `phloem_resolve()` — Münch pressure-driven sugar. (3a) leaves load surplus above phloem reserve into parent stem; (3b) BFS propagates pressure-driven flow through conduit network with time budget; (3c) meristems unload from nearest stem via permeability constant. Sugar does NOT diffuse locally.
   - `xylem_resolve()` — Phase 1/Phase 2 bulk water + cytokinin from roots to shoot. Water does NOT diffuse locally.

2. **Local diffusion** (`Node::transport_with_children()`) — per-node during DFS walk. Handles **auxin, gibberellin, cytokinin, and stress** only (sugar and water are skipped). On vascular-to-vascular edges, cytokinin is also skipped (xylem handles it). Leaves, meristems, and non-vascular nodes still get cytokinin via last-mile local diffusion from the nearest vascular stem.

| Chemical | Transport | Notes |
|----------|-----------|-------|
| Auxin | Local diffusion only | Polar cell-to-cell, basipetal bias (-0.1) |
| Gibberellin | Local diffusion only | Short-range, local to producing leaf |
| Stress | Local diffusion only | Local mechanical alarm signal |
| Ethylene | Spatial gas diffusion | Global 3D pass (reset each tick) |
| Sugar | **Phloem (phloem_resolve)** only | Münch pressure-flow; no local diffusion |
| Water | **Xylem (xylem_resolve)** only | Root → shoot bulk flow; no local diffusion |
| Cytokinin | **Xylem** + local (last-mile) | Skipped on vascular-to-vascular edges; local diffusion handles tips/leaves |

**Vasculature admission:** `has_vasculature(node, genome)` — STEM or ROOT with `radius >= vascular_radius_threshold` (0.01 dm). Seed always true. Leaves, meristems, and tiny nodes stay on local diffusion only.

### Phloem (Münch pressure-flow)

`pressure = (sugar / phloem_vol) × osmotic_coeff × water_frac`

where `phloem_vol` = ring area × length for stems, sphere volume for meristems, leaf area × thickness for leaves. Pressure flows from high (sources = sugar-loaded leaves/stems) to low (sinks = empty tips). BFS uses a time budget (1.0 tick) with `time_cost = length / speed`; speed scales with `r² / r_ref²` (Hagen-Poiseuille).

### Vascular Sources and Sinks

**Phloem (sugar):** Leaves with sugar above `phloem_reserve_fraction × sugar_cap` load sugar into parent stem pre-BFS. After BFS propagates pressure gradients, meristems (APICAL, ROOT_APICAL) unload from their parent stem via `phloem_unloading_meristem` permeability. Pipe capacity = `π × radius² × phloem_conductance`.

**Xylem (water + cytokinin):** Root nodes and root apicals are sources. Leaves (transpiration demand) and shoot apicals (turgor + cytokinin for growth) are sinks. Pipe capacity = `π × radius² × xylem_conductance`.

### Local Diffusion: Concentration + Shifted Equilibrium + Throughput Cap + Equalization

```
concentration = has_capacity ? level / capacity : level
diff = (parent_conc - child_conc) + bias
desired = diff * diffusion_rate [* avg_cap if capped] * radius_factor
desired = clamp(desired, -max_transport, max_transport)
equalize = flow_to_reach_bias_adjusted_equilibrium   // never overshoot
desired = clamp(desired, min(0, equalize), max(0, equalize))
```

Five mechanisms work together:

1. **Concentration gradient** — flow driven by relative fullness (sugar) or raw level (hormones). Prevents large nodes from draining small ones at equal concentration.
2. **Shifted equilibrium (bias)** — offsets where "equal" is. `bias < 0` = chemical accumulates root-ward (auxin). `bias > 0` = chemical accumulates tip-ward (cytokinin). **Root-type inversion:** bias sign flips for `ROOT` and `ROOT_APICAL` children so the seed is a correct transit junction: cytokinin flows root-tip → seed → shoot, auxin flows shoot → seed → root-tip, without cycling.
3. **Throughput cap** — `max_transport = base_transport + radius_factor * transport_scale`. Bottlenecked by the thinner of two connected nodes. `base_transport` floor ensures even thin tips can signal.
4. **Equalization clamp** — flow can never overshoot the bias-adjusted equilibrium. For uncapped chemicals: `equalize = (parent - child + bias) / 2`. For capped: solved from equal-concentration constraint. Prevents oscillation.
5. **Multi-way equalization cap** — limits total Phase 2 outflow when per-child equalization clamps would together drain the parent to zero. Computes multi-way equilibrium (parent retains proportional share). Only activates for branching nodes; chain topologies use per-pair clamps.

**Anti-teleportation:** All Phase 2 outflows go into `transport_received` buffer, not `chemical()`. Flushed after the node's own transport completes. Chemicals move at most one hop per tick.

**Capacity model:**
- **Sugar / Water** (resources): per-node capacity proportional to tissue volume. Concentration-based diffusion, destination headroom clamped.
- **Hormones** (signals): no capacity. Concentration = raw level. Decay limits accumulation.

**Auxin** (bias -0.1, basipetal):
- Persistent across ticks (not reset)
- Produced by `ApicalNode` each tick (modulated by growth rate)
- Also produced by growing `LeafNode` proportional to growth rate (zero when full-size)
- Shoot axillary buds sense **parent's** auxin level; activate when low
- Decays by `auxin_decay_rate` per tick

**Cytokinin** (bias +0.1, acropetal):
- Persistent across ticks (not reset)
- Produced by `RootApicalNode` gated by auxin (`root_cytokinin_production_rate × local_auxin`)
- Feedback loop: shoot auxin → flows to root tips → cytokinin produced → flows back to shoot
- Decays by `cytokinin_decay_rate` per tick
- Root axillary activation gated by `root_cytokinin_inhibition_threshold`

**Sugar** (bias 0, gradient):
- Persistent across ticks
- Produced by `LeafNode` via photosynthesis (gated by stomatal conductance — see Water Model)
- Concentration-based diffusion with radius-bottlenecked throughput cap
- Consumed by all nodes (maintenance costs in WorldParams)

## Secondary Thickening (Cambium)

Radial growth is driven by vascular history, not node age:

```
delta_radius = cambium_responsiveness × structural_flow_bias × sugar_gf × stress_boost
```

- `structural_flow_bias` stored on the parent, keyed by child pointer — connections where auxin has flowed repeatedly have stronger canalization and thicker cambium
- `sugar_gf` = `sugar / max_cost`, capped at 1 — sugar-limited when starved
- `stress_boost = 1 + stress × stress_thickening_boost` — mechanical load accelerates thickening
- Early return if `structural_flow_bias < 1e-6` — zero-flux connections never thicken (monocot behavior)
- **No age gate** — `cambium_maturation_ticks` is gone; the bias threshold is the gate

The self-reinforcing loop: main axis carries most flux → highest bias → fastest thickening → widest pipe → more flow → more bias. Lateral branches with weak canalization stay thin automatically.

## Canalization Model

Auxin flow through parent-child connections builds transport preference over time. Two layers of memory on the parent node, keyed by child pointer:

**Transient bias (auxin_flow_bias)** — fast, reversible. Represents PIN protein polarization.
- `target = auxin_flux × transient_gain`, bias chases target exponentially at `transient_rate`
- Decays toward 0 when flux stops.

**Structural bias (structural_flow_bias)** — slow, permanent. Represents built xylem/phloem.
- Grows by `structural_growth_rate` per tick when auxin flux exceeds `structural_threshold`
- Never decays. Capped at `structural_max`.
- Also gates vasculature admission (`has_vasculature`) and cambial thickening rate.

**Effect on transport:** Both biases multiply the conductance weight for each child connection. `weight = pipe_cap × (1 + canalization_weight × (flow_bias + structural_bias))`. This redistributes flow (all chemicals, not just auxin) among siblings proportionally — does **not** amplify total flow.

**Lifecycle:** Biases transfer on `replace_child` (chain growth preserves branch history). New children start at 0, 0. Cleaned up on `die()`.

## Water Model

Water is a **persistent** capacity-based resource (like sugar):

- **Absorbed** by `RootNode` and `RootApicalNode` proportional to surface area and `soil_moisture`. Absorption is gradient-based: `desired = (soil_moisture - water_conc) × absorption_rate × surface_area` — self-limiting when the root is full (concentration approaches soil moisture).
- **Lost** by leaves via transpiration (proportional to leaf area and light exposure) and photosynthesis water cost
- **Stomatal conductance** — `stomatal = clamp(water / water_cap, 0.2, 1.0)` scales photosynthesis output. Water deficit partially closes stomata, reducing sugar production. (Minimum 20% conductance — stressed leaves still photosynthesize slowly.)
- **Transported** via `xylem_resolve()` bulk flow only — no local diffusion (skipped in `transport_with_children`)

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

`Plant::tick()` runs four phases in order:
1. `pin_transport(plant, genome)` — polar auxin transport (global pre-pass)
2. `tick_recursive(seed)` — DFS walk: parent ticks before children (grow, consume, produce)
3. `phloem_resolve(plant, genome, world)` — Münch pressure-driven sugar transport
4. `xylem_resolve(plant, genome, world)` — Phase 1/Phase 2 water + cytokinin

Each `Node::tick()` (non-virtual, same for all types):
1. `age++`
2. `sync_world_position()` — recompute world position from offset chain
3. `pay_maintenance()` — deduct sugar; die and return if starvation ticks exceed max
4. `update_tissue(plant, world)` — **virtual** type-specific growth (see below)
5. `sync_world_position()` — re-sync after tissue may have changed offset
6. `update_physics()` — compute mass, stress; apply droop and break if overstressed
7. `transport_chemicals()` → `transport_with_children()` + `update_canalization()` + `decay_chemicals()`
8. Flush `transport_received` buffer into `chemical()`

Children are ticked recursively after step 8 using a snapshot of the children list (safe for reparenting).

**`update_tissue` convention** — subclass-specific growth happens entirely inside this virtual method. Recommended call order within `update_tissue`: produce → consume → evaluate → grow → orient. Growth verbs: `elongate` (extend offset), `thicken` (increase radius), `expand` (increase leaf_size). Meristem chain growth and axillary activation also happen inside `update_tissue` via `Plant::queue_removal()`.

## Key Design Decisions
- **Meristems are nodes** — 2 meristem types (`ApicalNode`, `RootApicalNode`) extend `Node` directly. Real children in the tree graph, participate in chemical transport naturally. No axillary node types.
- **Local diffusion rates are low by design** — auxin/GA/stress use low diffusion rates (0.05–0.10) so they remain cell-to-cell signals, not whole-plant gradients. Radius bottleneck handles the rest.
- **Vascular admission is radius-gated** — `has_vasculature()` returns true for STEM/ROOT with `radius >= vascular_radius_threshold` (0.01 dm). No bias-map needed; new internodes qualify immediately at default radius (0.015 dm).
- **Thickening is bias-driven** — `cambium_responsiveness × auxin_flow_bias` drives radial growth. Zero-flux branches never thicken.
- **Phloem is Münch pressure-flow** — sugar flows from high osmotic pressure (sugar-loaded sources) to low pressure (empty sinks) via BFS on conduit network. Flow is NOT conductance-weighted by auxin_flow_bias.
- **Flat node hierarchy** — `Node` base class with 5 direct subclasses. Each owns type-specific fields and behavior via `virtual update_tissue()`. Downcasting via `as_stem()`/`as_root()`/`as_leaf()`/`as_apical()`/`as_root_apical()` — fast `static_cast` gated on `NodeType` enum, no RTTI.
- **Root elongation is sugar-gated only** — root apical `elongate()` uses `sugar_growth_fraction()` with no auxin gate. Shoot-derived auxin affects root *activation*, not elongation.
- **Dormant meristems cost zero maintenance** — only active meristems pay `sugar_maintenance_meristem`
- **Sugar economy: physics in WorldParams, strategy in Genome** — production rate, maintenance costs, and construction costs are physical constants in WorldParams (not evolvable). Evolution optimizes plant shape/structure, not metabolism.

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
- `root_cytokinin_production_rate` (0.15) - baseline cytokinin produced per tick by root apicals, scaled by local auxin
- `root_auxin_activation_threshold` (0.05) - minimum auxin to activate a dormant root meristem
- `root_cytokinin_inhibition_threshold` (0.15) - cytokinin above this inhibits new root meristem activation
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
- `cambium_responsiveness` (0.00002 dm/hr·bias) - thickening rate per unit structural_flow_bias. Main trunk (bias ~2.0) thickens at ~0.00004 dm/hr; laterals with weak canalization thicken proportionally less.
- `branch_angle` (0.785 rad / 45°) - angle of shoot branches from parent stem
- `root_branch_angle` (0.35 rad / 20°) - angle of root branches
- `ga_production_rate` (0.5) - GA per dm leaf_size per tick from young leaves
- `ga_leaf_age_max` (168) - only leaves younger than 7 days produce GA
- `ga_elongation_sensitivity` (2.0) - GA boost to elongation rate
- `ga_length_sensitivity` (1.5) - GA boost to max internode length
- `ethylene_abscission_threshold` (0.5) - triggers leaf senescence
- `ethylene_shade_threshold` (0.3) - light_exposure below which shade-ethylene kicks in
- `senescence_duration` (48) - ticks from senescence to leaf drop
- `transient_gain` (2.0) - target transient bias per unit of auxin flux
- `transient_rate` (0.2) - exponential chase speed for transient bias (~87% in 8 hours)
- `structural_threshold` (0.15) - minimum auxin flux per tick to grow structural bias (filters noise)
- `structural_growth_rate` (0.005) - structural bias increment per tick above threshold (~8 days to reach 1.0)
- `structural_max` (2.0) - cap on structural bias. At max: connection weight = 1 + 2.0 = 3×
- `canalization_weight` (1.0) - global scaling on bias effect. 0 = canalization disabled entirely
- `vascular_conductance_threshold` (0.005) - minimum structural_flow_bias for vascular admission. Just below the initial stamp so new internodes join immediately; zero-flux nodes stay excluded.
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
- `xylem_conductance` (10.0) - vascular throughput per dm² cross-section per tick (water + cytokinin)
- `phloem_conductance` (8.0) - vascular throughput per dm² cross-section per tick (sugar)
- `phloem_reserve_fraction` (0.3) - fraction of sugar_cap leaves keep for themselves (not loaded into phloem)
