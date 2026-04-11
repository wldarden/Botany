# Botany Simulator

A plant growth simulator using hormone-driven meristem mechanics. Plants grow from a seed node that spawns one shoot apical meristem (upward) and one root apical meristem (downward). Growth is governed by auxin/cytokinin signaling.

## Build & Run

```bash
# Build (cmake is at /usr/local/bin/cmake)
/usr/local/bin/cmake --build build

# Run realtime viewer (must run from project root for shader path)
./build/botany_realtime [--color auxin|cytokinin|sugar|gibberellin|ethylene|type]

# Run tests
./build/botany_tests
```

## Architecture

### Engine (`src/engine/`)
- **genome.h** - `Genome` struct with all tunable parameters + `default_genome()`
- **node.h** - `Node` (position, radius, parent/children, hormone levels, sugar, leaf_size), `Meristem` base class, `MeristemType` enum. `NodeType`: STEM, ROOT, LEAF
- **meristems/** - Meristem subfolder:
  - `meristem.h/cpp` — `tick_meristems()` dispatch + secondary/intercalary growth
  - `meristem_types.h` — convenience umbrella, includes all 4 type headers
  - `shoot_apical.h/cpp` — `ShootApicalMeristem` (chain growth, phyllotaxis)
  - `shoot_axillary.h/cpp` — `ShootAxillaryMeristem` (auxin-gated activation)
  - `root_apical.h/cpp` — `RootApicalMeristem` (gravitropism, root chain growth)
  - `root_axillary.h/cpp` — `RootAxillaryMeristem` (cytokinin-gated activation)
  - `helpers.h` — shared helper functions (growth_direction, branch_direction, perturb, etc.)
- **hormone.cpp** - `transport_auxin()` and `transport_cytokinin()` - per-tick hormone reset + transport
- **gibberellin.h/cpp** - `compute_gibberellin()` — local GA production by young leaves
- **ethylene.h/cpp** - `compute_ethylene()` + `process_abscission()` — spatial gas diffusion, leaf abscission
- **sugar.h/cpp** - `transport_sugar()` — sugar production, gradient-based diffusion, maintenance consumption
- **world_params.h** - `WorldParams` struct (light level, diffusion iterations) — non-genetic sim parameters
- **plant.h/cpp** - `Plant` class owns all nodes/meristems, has root meristem cap (100)
- **engine.h/cpp** - `Engine` orchestrates tick order: auxin → cytokinin → GA → sugar → ethylene → abscission → meristems → positions

### Renderer (`src/renderer/`)
- OpenGL 4.1 core profile, GLFW window, orbit camera
- Draws cylinders between parent-child nodes, leaves as quads
- Color modes: default (brown), chemical heatmap (auxin/cytokinin/sugar/gibberellin/ethylene), type (green=shoot, orange=root)

### Apps
- **app_realtime.cpp** - Interactive viewer with pause/speed controls
- **app_headless.cpp** - Headless precompute, saves binary recording
- **app_playback.cpp** - Playback viewer with ImGui scrubbing

## Hormone Model

Hormones are **reset to zero every tick** then recomputed as a fresh signal snapshot:

**Auxin** (shoot branching control):
- Produced by active shoot apical meristems only
- Flows basipetally (child -> parent, toward seed) via `auxin_collect`
- Small spillback fraction redistributes from junctions back into branches (`auxin_spillback`)
- Shoot axillary buds sense **parent stem node's** auxin level (not their own)
- Activate when `parent.auxin < auxin_threshold` (far from any active shoot tip)

**Cytokinin** (root branching control):
- Produced by active root apical meristems only
- Flows toward seed via `cytokinin_collect`, then distributes to children via `cytokinin_distribute`
- Root axillary buds sense **parent root node's** cytokinin level
- Activate when `parent.cytokinin < cytokinin_threshold` (far from any active root tip)

## Gibberellin Model

GA is **reset to zero every tick** (signal model, same as auxin/cytokinin):

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

## Sugar Model

Sugar **persists across ticks** (NOT reset like hormones). Four phases per tick:

1. **Production** — LEAF nodes produce: `sugar += light_level * leaf_size * sugar_production_rate`
   - Feedback inhibition: production skipped if node sugar >= storage cap
2. **Diffusion** — Gradient-based bidirectional flow through all node connections:
   - Transport capacity = `min_radius^2 * PI * sugar_transport_conductance` (thicker = more)
   - LEAF connections use baseline capacity (leaf radius is 0)
   - Runs `sugar_diffusion_iterations` passes per tick (WorldParams, default 5)
   - Cap-aware: transfers clamped by receiver's available headroom
3. **Leaf growth** — Growing leaves spend sugar to expand toward max_leaf_size
4. **Consumption** — Every node deducts volume-based maintenance cost:
   - LEAF: `sugar_maintenance_leaf * leaf_size²` (scales with leaf area)
   - STEM: `sugar_maintenance_stem * π * r² * internode_length` (scales with tissue volume)
   - ROOT: `sugar_maintenance_root * π * r² * internode_length` (scales with tissue volume)
   - Active meristem tips: `+ sugar_maintenance_meristem` (flat per-tip)
   - Safety clamp: sugar capped to node storage limit after consumption

**Storage caps** — Each node has a maximum sugar capacity proportional to its tissue volume:
- STEM/ROOT: `π * r² * internode_length * sugar_storage_density_wood`
- LEAF: `leaf_size² * sugar_storage_density_leaf`
- Minimum cap of `sugar_cap_minimum` for tiny/new nodes
- Seed node cap is at least `seed_sugar` to hold initial reserves

**WorldParams** (non-genetic, on Engine):
- `light_level` (1.0) — global light intensity, controls sugar production
- `sugar_diffusion_iterations` (5) — simulation quality for diffusion smoothing

## Key Design Decisions
- Leaves are real LEAF graph nodes (NodeType::LEAF), not struct properties — they have position, sugar, leaf_size
- Chain growth creates 3 children on interior STEM nodes: continuation tip, axillary meristem, and LEAF node
- Axillary buds check their **parent's** hormone level, not their own node (hormones flow through the stem, not into side branches)
- Root meristems are hard-capped at 100 (`Plant::max_root_meristems`) to prevent runaway root growth
- The cap only gates creation of new axillary buds during chain growth, not activation of existing ones
- Hormone reset each tick prevents accumulation artifacts (base of trunk having more auxin than tip)
- Sugar persists across ticks (resource model, not signal model like hormones)

## Tuning Parameters (genome.h)
- `auxin_threshold` (0.15) - lower = fewer shoot branches, higher = more
- `cytokinin_threshold` (0.15) - lower = fewer root branches, higher = more
- `auxin_spillback_rate` (0.1) - how much junction auxin spills back into branches
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
