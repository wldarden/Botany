# Botany Simulator

A plant growth simulator using hormone-driven meristem mechanics. Plants grow from a seed node that spawns one shoot apical meristem (upward) and one root apical meristem (downward). Growth is governed by auxin/cytokinin signaling.

## Build & Run

```bash
# Build (cmake is at /usr/local/bin/cmake)
# IMPORTANT: always rebuild before running tests — LSP diagnostics can be stale
/usr/local/bin/cmake --build build

# Run realtime viewer (must run from project root for shader path)
./build/botany_realtime [--color auxin|cytokinin|sugar|gibberellin|ethylene|type]

# Run tests
./build/botany_tests
```

## Architecture

### Engine (`src/engine/`)
- **genome.h** - `Genome` struct with all tunable parameters + `default_genome()`
- **node/** - Node subfolder:
  - `node.h/cpp` — `Node` base class (position, radius, parent/children, hormones, sugar) with virtual `tick()`, `transport_chemicals()`. `NodeType` enum: `STEM, ROOT, LEAF, SHOOT_APICAL, SHOOT_AXILLARY, ROOT_APICAL, ROOT_AXILLARY`. Downcasting via `as_stem()`, `as_root()`, `as_leaf()`, `as_meristem()`, `as_shoot_apical()`, etc.
  - `stem_node.h/cpp` — `StemNode` (thickening, intercalary elongation)
  - `root_node.h/cpp` — `RootNode` (same with root params)
  - `leaf_node.h/cpp` — `LeafNode` (owns `leaf_size`, `light_exposure`, `senescence_ticks`; phototropism + size growth)
  - `meristem_node.h/cpp` — `MeristemNode` base class (owns `active`, `ticks_since_last_node`; `is_tip()`)
  - **meristems/** — Meristem node subfolder:
    - `meristem.h/cpp` — `tick_meristems()` recursively walks the tree from seed, calling `node.tick()` on each; snapshots children before iterating (meristems may reparent). Calls `plant.flush_removals()` after.
    - `meristem_types.h` — convenience umbrella, includes all 4 type headers
    - `shoot_apical.h/cpp` — `ShootApicalNode` extends `MeristemNode` (chain growth via self-reparenting, phyllotaxis, auxin production)
    - `shoot_axillary.h/cpp` — `ShootAxillaryNode` extends `MeristemNode` (auxin-gated activation, replaces self with apical)
    - `root_apical.h/cpp` — `RootApicalNode` extends `MeristemNode` (gravitropism, root chain growth, cytokinin production)
    - `root_axillary.h/cpp` — `RootAxillaryNode` extends `MeristemNode` (cytokinin-gated activation)
    - `helpers.h` — shared helper functions (growth_direction, branch_direction, perturb, sugar_growth_fraction, etc.)
- **gibberellin.h/cpp** - `compute_gibberellin()` — local GA production by young leaves (reset each tick)
- **ethylene.h/cpp** - `compute_ethylene()` + `process_abscission()` — spatial gas diffusion, leaf abscission (reset each tick)
- **sugar.h/cpp** - `transport_sugar()` — sugar production (leaves), maintenance consumption, starvation pruning. Diffusion moved to `Node::transport_chemicals()`.
- **hormone.h/cpp** - Empty placeholder (auxin/cytokinin transport moved to `Node::transport_chemicals()`)
- **world_params.h** - `WorldParams` struct (light level, construction costs) — non-genetic sim parameters
- **plant.h/cpp** - `Plant` class owns all nodes, has root meristem cap (100). `Plant::tick()` orchestrates per-plant tick order. `queue_removal()` / `flush_removals()` for deferred node cleanup.
- **engine.h/cpp** - `Engine` iterates plants and calls `plant.tick(world_params)`

### Node Class Hierarchy
```
Node (base — position, chemicals, tick, transport_chemicals)
├── StemNode (thickening, intercalary elongation)
├── RootNode (same with root params)
├── LeafNode (leaf_size, light_exposure, phototropism, growth)
└── MeristemNode (active, ticks_since_last_node, is_tip)
    ├── ShootApicalNode (chain growth, phyllotaxis, auxin production)
    ├── ShootAxillaryNode (auxin-gated activation)
    ├── RootApicalNode (gravitropism, chain growth, cytokinin production)
    └── RootAxillaryNode (cytokinin-gated activation)
```

Meristems are real nodes in the tree graph, not separate objects. They participate in chemical transport naturally.

### Renderer (`src/renderer/`)
- OpenGL 4.1 core profile, GLFW window, orbit camera
- Draws cylinders between parent-child nodes, leaves as quads
- Color modes: default (brown), chemical heatmap (auxin/cytokinin/sugar/gibberellin/ethylene), type (green=shoot, orange=root, red/blue=meristems)

### Apps
- **app_realtime.cpp** - Interactive viewer with pause/speed controls
- **app_headless.cpp** - Headless precompute, saves binary recording
- **app_playback.cpp** - Playback viewer with ImGui scrubbing

## Chemical Transport Model

All chemicals are transported **locally** by each node during `Node::transport_chemicals()`, called from `Node::tick()`. The recursive tick walks the tree from seed outward (DFS pre-order), so each node handles its own transport with its parent.

### Generic Transport with Directional Bias

Each chemical uses a blend of gradient diffusion and directional push, controlled by `directional_bias` (-1 to +1):

```
gradient_weight = 1 - |bias|
directional_weight = |bias|
gradient_flow = (my_val - parent_val) * gradient_weight * rate
directional_flow = bias < 0 ? my_val * dw * rate : -parent_val * dw * rate
flow = gradient_flow + directional_flow  (clamped to available supply)
```

**Auxin** (bias -0.9, basipetal):
- **Persistent** across ticks (not reset)
- Produced by active `ShootApicalNode` during its `tick()`
- 90% directional push toward parent (seed-ward), 10% gradient
- Decays by `auxin_decay_rate` per tick
- Shoot axillary buds sense **parent's** auxin level; activate when low

**Cytokinin** (bias +0.9, acropetal):
- **Persistent** across ticks (not reset)
- Produced by active `RootApicalNode` during its `tick()`
- 90% directional pull from parent (shoot-ward), 10% gradient
- Decays by `cytokinin_decay_rate` per tick
- Root axillary buds sense **parent's** cytokinin level; activate when low

**Sugar** (bias 0, gradient):
- **Persistent** across ticks
- Produced by leaf nodes (global `produce_sugar()` pass — involves shadow casting)
- 100% gradient diffusion with radius-based conductance and cap clamping
- Consumed by all nodes (global `consume_sugar()` pass)

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

## Tick Control Flow

The recursive tick walks the tree from seed outward. Each node:
1. Produces chemicals if applicable (meristem nodes produce auxin/cytokinin)
2. Calls base `Node::tick()` → `age++`, `transport_chemicals()`
3. Does type-specific work (growth, activation, etc.)
4. Children are ticked recursively (snapshot of children list for safety)

Meristem chain growth: the meristem node inserts an internode above itself (self-reparenting). Axillary activation: the node replaces itself in the parent's children with a new apical node, queues itself for deferred removal.

## Key Design Decisions
- **Meristems are nodes** — `MeristemNode` base with 4 subclasses. Real children in the tree graph, participate in chemical diffusion naturally. No separate `Meristem` objects.
- **Local chemical transport** — each node handles its own transport during `tick()` via `transport_chemicals()`. No global tree passes for auxin/cytokinin/sugar diffusion.
- **Directional bias** — generic transport function blends gradient diffusion with directional push. Allows same code for basipetal (auxin), acropetal (cytokinin), and bidirectional (sugar).
- **Node class hierarchy** — `Node` base class with `StemNode`, `RootNode`, `LeafNode`, `MeristemNode` subclasses. Each subclass owns its type-specific fields and growth behavior via `virtual tick()`. Downcasting via `as_stem()`/`as_root()`/`as_leaf()`/`as_meristem()` (fast `static_cast` gated on `NodeType` enum, no RTTI).
- Leaves are real `LeafNode` graph nodes — they own `leaf_size`, `light_exposure`, `senescence_ticks`
- Chain growth inserts an internode above the meristem: parent → new_internode → [meristem, axillary, leaf]
- Axillary buds check their **parent's** hormone level, not their own node
- Root meristems are hard-capped at 100 (`Plant::max_root_meristems`) to prevent runaway root growth
- Sugar persists across ticks (resource model); auxin/cytokinin persist too (not reset each tick)
- Gibberellin and ethylene still use reset-each-tick signal model (global passes)

## Tuning Parameters (genome.h)
- `auxin_threshold` (0.15) - lower = fewer shoot branches, higher = more
- `cytokinin_threshold` (0.15) - lower = fewer root branches, higher = more
- `auxin_directional_bias` (-0.9) - basipetal bias strength
- `cytokinin_directional_bias` (0.9) - acropetal bias strength
- `auxin_transport_rate` / `cytokinin_transport_rate` (0.3) - how fast hormones flow per tick
- `branch_angle` (0.785 rad / 45 deg) - angle of shoot branches from parent stem
- `root_branch_angle` (0.35 rad / 20 deg) - angle of root branches
- `ga_production_rate` (0.5) - GA per dm leaf_size per tick from young leaves
- `ga_leaf_age_max` (168) - only leaves younger than 7 days produce GA
- `ga_elongation_sensitivity` (2.0) - GA boost to elongation rate
- `ga_length_sensitivity` (1.5) - GA boost to max internode length
- `ethylene_abscission_threshold` (0.5) - triggers leaf senescence
- `ethylene_shade_threshold` (0.3) - light_exposure below which shade-ethylene kicks in
- `senescence_duration` (48) - ticks from senescence to leaf drop
