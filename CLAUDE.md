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
  - `node.h/cpp` — `Node` base class (position, radius, parent/children, hormones, sugar) with virtual `tick()`, `transport_chemicals()`. `NodeType` enum: `STEM, ROOT, LEAF, SHOOT_APICAL, SHOOT_AXILLARY, ROOT_APICAL, ROOT_AXILLARY`. Downcasting via `as_stem()`, `as_root()`, `as_leaf()`, `as_shoot_apical()`, `as_shoot_axillary()`, `as_root_apical()`, `as_root_axillary()`. `is_meristem()` type-group check.
  - `stem_node.h/cpp` — `StemNode` (thickening, intercalary elongation)
  - `root_node.h/cpp` — `RootNode` (same with root params)
  - `leaf_node.h/cpp` — `LeafNode` (owns `leaf_size`, `light_exposure`, `senescence_ticks`; `produce()` for photosynthesis/GA/abscission, `grow()` for phototropism + size growth)
  - **meristems/** — Meristem node subfolder:
    - `meristem_types.h` — convenience umbrella, includes all 4 type headers
    - `shoot_apical.h/cpp` — `ShootApicalNode` extends `Node` (chain growth via self-reparenting, phyllotaxis, auxin production; owns `ticks_since_last_node`)
    - `shoot_axillary.h/cpp` — `ShootAxillaryNode` extends `Node` (auxin-gated activation, replaces self with apical; owns `active`)
    - `root_apical.h/cpp` — `RootApicalNode` extends `Node` (gravitropism, root chain growth, cytokinin production; owns `ticks_since_last_node`)
    - `root_axillary.h/cpp` — `RootAxillaryNode` extends `Node` (cytokinin-gated activation; owns `active`)
    - `helpers.h` — shared helper functions (growth_direction, branch_direction, perturb, sugar_growth_fraction, etc.)
- **gibberellin.h/cpp** - `compute_gibberellin()` — local GA production by young leaves (reset each tick)
- **ethylene.h/cpp** - `compute_ethylene()` + `process_abscission()` — spatial gas diffusion, leaf abscission (reset each tick)
- **sugar.h/cpp** - `transport_sugar()` — sugar production (leaves), maintenance consumption, starvation pruning. Diffusion moved to `Node::transport_chemicals()`.
- **hormone.h/cpp** - Empty placeholder (auxin/cytokinin transport moved to `Node::transport_chemicals()`)
- **world_params.h** - `WorldParams` struct (light level, construction costs, sugar_production_rate, maintenance rates) — non-genetic physical constants
- **plant.h/cpp** - `Plant` class owns all nodes, has root meristem cap (100). `Plant::tick()` orchestrates per-plant tick order: global chemical passes, then `tick_tree()` (recursive DFS walk from seed, snapshots children for safe reparenting, flushes removals). `queue_removal()` / `flush_removals()` for deferred node cleanup.
- **engine.h/cpp** - `Engine` iterates plants and calls `plant.tick(world_params)`

### Node Class Hierarchy
```
Node (base — position, chemicals, tick, transport_chemicals)
├── StemNode (thickening, intercalary elongation)
├── RootNode (same with root params)
├── LeafNode (leaf_size, light_exposure, phototropism, growth)
├── ShootApicalNode (chain growth, phyllotaxis, auxin production)
├── ShootAxillaryNode (auxin-gated activation, owns `active`)
├── RootApicalNode (gravitropism, chain growth, cytokinin production)
└── RootAxillaryNode (cytokinin-gated activation, owns `active`)
```

All 7 node types extend `Node` directly — flat hierarchy, no intermediate classes. Meristems are real nodes in the tree graph, participating in chemical transport naturally.

### Renderer (`src/renderer/`)
- OpenGL 4.1 core profile, GLFW window, orbit camera
- Draws cylinders between parent-child nodes, leaves as quads
- Color modes: default (brown), chemical heatmap (auxin/cytokinin/sugar/gibberellin/ethylene), type (green=shoot, orange=root, red/blue=meristems)

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

## Chemical Transport Model

All chemicals use a **unified transport function** (`transport_chemical()` in `hormone.h`), called per-node per-chemical during `Node::transport_chemicals()`. The recursive tick walks the tree from seed outward (DFS pre-order), so each node handles its own transport with its parent. Every edge is processed exactly once.

### Unified Transport: Concentration + Shifted Equilibrium + Throughput Cap

```
concentration = has_capacity ? level / capacity : level
effective_diff = (my_conc - parent_conc) - bias
desired_flow = effective_diff * diffusion_rate * (has_cap ? avg_cap : 1) * radius_factor
flow = clamp(desired_flow, -max_transport, max_transport)
flow = clamp(flow, source_available, destination_headroom)
```

Three mechanisms work together:

1. **Concentration gradient** — flow driven by relative fullness (sugar) or raw level (hormones). Prevents large nodes from draining small ones at equal concentration.
2. **Shifted equilibrium (bias)** — offsets where "equal" is. `bias < 0` = chemical accumulates root-ward (auxin). `bias > 0` = chemical accumulates tip-ward (cytokinin). Self-correcting: gradient still governs flow, just from an offset resting point.
3. **Throughput cap** — `max_transport = base_transport + radius_factor * transport_scale`. Bottlenecked by the thinner of two connected nodes. `base_transport` floor ensures even thin tips can transport.

**Capacity model:**
- **Sugar** (resource): per-node capacity from `sugar_cap()` — proportional to tissue volume. Concentration-based diffusion, destination headroom clamped.
- **Hormones** (signals): no capacity (cap=0). Concentration = raw level. No destination clamp — hormones can accumulate freely. Decay limits accumulation instead.

**Auxin** (bias -0.1, basipetal):
- **Persistent** across ticks (not reset)
- Produced by active `ShootApicalNode` during its `tick()`, modulated by growth rate (growth_gf multiplier)
- Also produced by growing `LeafNode` during `grow_size()` — proportional to growth rate, zero when full-size
- Shifted equilibrium pushes auxin root-ward; decay creates a gradient from apex to base
- Decays by `auxin_decay_rate` per tick
- Shoot axillary buds sense **parent's** auxin level; activate when low

**Cytokinin** (bias +0.1, acropetal):
- **Persistent** across ticks (not reset)
- Produced by active `RootApicalNode` during its `tick()`
- Shifted equilibrium pushes cytokinin tip-ward
- Decays by `cytokinin_decay_rate` per tick
- Root axillary buds sense **parent's** cytokinin level; activate when low

**Sugar** (bias 0, gradient):
- **Persistent** across ticks
- Produced by leaf nodes via `LeafNode::produce()` (photosynthesis)
- Concentration-based diffusion with radius-bottlenecked throughput cap
- Consumed by all nodes (maintenance costs in WorldParams)
- Leaves use effective petiole radius (`leaf_size * initial_radius`) for transport — guarantees they can export full production

## Gibberellin Model

GA is **reset to zero every tick** (signal model):

- Produced by young LEAF nodes only (`leaf_age < ga_leaf_age_max`)
- Applied locally — `ga_production_rate * leaf_size` is added to the leaf's parent and grandparent stem nodes
- Not transported through the rest of the tree; effect is purely local to the internode being elongated
- **Effects:** boosts internode elongation rate (`* ga_elongation_sensitivity`) and max internode length (`* ga_length_sensitivity`) on nodes that receive GA

## Ethylene Model

Ethylene is **reset to zero every tick** (signal model). Four production triggers per leaf:

1. **Sugar starvation** — leaf sugar below maintenance threshold
2. **Shade** — leaf `light_exposure < ethylene_shade_threshold`
3. **Old age** — leaf age exceeds species maximum
4. **Crowding** — local node density above threshold

**Spatial diffusion** — ethylene spreads as a gas through 3D space (NOT via the tree graph). Each emitting node contributes to all nodes within `diffusion_radius`, attenuated by distance.

**Effects:**
- **Leaf abscission** — if a leaf's ethylene exceeds `ethylene_abscission_threshold`, senescence begins; the leaf yellows and is removed after `senescence_duration` ticks
- **Elongation inhibition** — high ethylene suppresses internode elongation in nearby stem nodes

**Abscission lifecycle:** ethylene > threshold → senescence flag set → leaf gradually yellows (visual only) → removed from graph after `senescence_duration` ticks have elapsed

## Canalization Model

Auxin flow through parent-child connections builds transport preference over time. Two layers of memory:

**Transient bias (auxin_flow_bias)** — fast, reversible. Represents PIN protein polarization.
- Each tick: `target = auxin_flux * transient_gain`, bias chases target exponentially at `transient_rate`
- Decays toward 0 when flux stops. Responds within hours/days.

**Structural bias (structural_flow_bias)** — slow, permanent. Represents built xylem/phloem.
- Grows by `structural_growth_rate` per tick when auxin flux exceeds `structural_threshold`
- Never decays. Capped at `structural_max`.

**Effect on transport:** Both biases stored on the parent node, keyed by child pointer. During `transport_with_children()`, each child's proportional weight is multiplied by `1 + canalization_weight * (flow_bias + structural_bias)`. This redistributes chemical flow (all chemicals, not just auxin) among siblings — biased connections get a larger share of the same total. Does not amplify total flow.

**Lifecycle:** Biases transfer on `replace_child` (chain growth preserves branch history). New children start at 0, 0. Cleaned up on `die()`.

## Tick Control Flow

The recursive tick walks the tree from seed outward. Each node's `tick()`:
1. Position (recomputed from parent + offset)
2. `age++`
3. Maintenance (sugar deducted based on WorldParams rates)
4. Sugar cap clamp
5. Starvation tracking (resets to 0 if sugar > 0; death after `starvation_ticks_max`)
6. `produce()` — virtual: LeafNode does photosynthesis + GA + abscission tracking; base is no-op
7. `grow()` — virtual: type-specific growth (tip extension, thickening, elongation, leaf size)
8. Mass & stress computation
9. Droop & break (stems only)
10. `transport_chemicals()` — unified diffusion for all chemicals
11. Children ticked recursively (snapshot of children list for safe reparenting)

Meristem chain growth: the meristem node inserts an internode above itself (self-reparenting). Axillary activation: the node replaces itself in the parent's children with a new apical node, queues itself for deferred removal.

## Key Design Decisions
- **Meristems are nodes** — 4 meristem types extend `Node` directly (flat hierarchy). Real children in the tree graph, participate in chemical diffusion naturally. No intermediate `MeristemNode` class.
- **Local chemical transport** — each node handles its own transport during `tick()` via `transport_chemicals()`. No global tree passes for auxin/cytokinin/sugar diffusion.
- **Unified transport with shifted equilibrium** — single `transport_chemical()` function handles all chemicals. Concentration-based diffusion (sugar uses per-node capacity, hormones use raw levels) with shifted equilibrium bias and radius-bottlenecked throughput cap. Same function for basipetal (auxin), acropetal (cytokinin), and bidirectional (sugar).
- **Flat node hierarchy** — `Node` base class with 7 direct subclasses (`StemNode`, `RootNode`, `LeafNode`, `ShootApicalNode`, `ShootAxillaryNode`, `RootApicalNode`, `RootAxillaryNode`). Each subclass owns its type-specific fields and growth behavior via `virtual tick()`. Downcasting via `as_stem()`/`as_root()`/`as_leaf()`/`as_shoot_apical()`/etc. (fast `static_cast` gated on `NodeType` enum, no RTTI).
- Leaves are real `LeafNode` graph nodes — they own `leaf_size`, `light_exposure`, `senescence_ticks`
- Chain growth inserts an internode above the meristem: parent → new_internode → [meristem, axillary, leaf]
- Axillary buds check their **parent's** hormone level, not their own node
- Root meristems are hard-capped at 100 (`Plant::max_root_meristems`) to prevent runaway root growth
- Sugar persists across ticks (resource model); auxin/cytokinin persist too (not reset each tick)
- Gibberellin and ethylene still use reset-each-tick signal model (global passes)
- **Produce vs grow** — `Node::produce()` handles resource/signal generation (photosynthesis, GA), `Node::grow()` handles structural change (tip extension, thickening, leaf expansion). Called as separate steps in tick pipeline — production isn't gated by growth.
- **Sugar economy: physics in WorldParams, strategy in Genome** — production rate, maintenance costs, and construction costs are physical constants in WorldParams (not evolvable). Evolution optimizes plant *shape* (leaf size, branching, growth rates), not metabolism. Storage capacity and transport params remain evolvable.

## Tuning Parameters (genome.h)
- `auxin_threshold` (0.15) - lower = fewer shoot branches, higher = more
- `auxin_thickening_threshold` (0.03) - auxin level for full-speed cambial thickening (much lower than branching threshold — cambium is very sensitive to auxin)
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
- `cytokinin_bias` (0.1) - shifted equilibrium for acropetal flow (positive = toward tips)
- `auxin_diffusion_rate` / `cytokinin_diffusion_rate` (0.3) - gradient responsiveness
- `hormone_base_transport` (0.5) - throughput floor for hormones (even thin tips can signal)
- `hormone_transport_scale` (1.0) - how much radius amplifies hormone throughput
- `sugar_base_transport` (0.01) - throughput floor for sugar (small — sugar is radius-dependent)
- `sugar_transport_scale` (5.0) - how much radius amplifies sugar throughput
- `apical_auxin_baseline` (0.15) - base auxin output of shoot apical meristem per tick
- `apical_growth_auxin_multiplier` (2.0) - growth bonus: total = baseline * (1 + multiplier * growth_fraction). 0 = no bonus, 2 = 3x at max growth
- `leaf_auxin_baseline` (0.15) - scaling constant for leaf auxin production (decoupled from apical rate)
- `leaf_growth_auxin_multiplier` (0.1) - fraction of leaf_auxin_baseline at max growth. Single leaf at max growth = 10% of apical baseline
- `branch_angle` (0.785 rad / 45 deg) - angle of shoot branches from parent stem
- `root_branch_angle` (0.35 rad / 20 deg) - angle of root branches
- `ga_production_rate` (0.5) - GA per dm leaf_size per tick from young leaves
- `ga_leaf_age_max` (168) - only leaves younger than 7 days produce GA
- `ga_elongation_sensitivity` (2.0) - GA boost to elongation rate
- `ga_length_sensitivity` (1.5) - GA boost to max internode length
- `ethylene_abscission_threshold` (0.5) - triggers leaf senescence
- `ethylene_shade_threshold` (0.3) - light_exposure below which shade-ethylene kicks in
- `senescence_duration` (48) - ticks from senescence to leaf drop
- `transient_gain` (2.0) - target transient bias per unit of auxin flux. Higher = stronger short-term flow reinforcement
- `transient_rate` (0.2) - exponential chase speed for transient bias (0.2 = ~87% in 8 hours)
- `structural_threshold` (0.05) - minimum auxin flux to build structural bias (filters noise)
- `structural_growth_rate` (0.005) - structural bias growth per tick when flux exceeds threshold (~8 days to reach 1.0)
- `structural_max` (2.0) - cap on structural bias. At max: connection gets 1 + 2.0 = 3x weight
- `canalization_weight` (1.0) - global scaling on bias effect. 0 = canalization disabled entirely
