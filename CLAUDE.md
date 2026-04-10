# Botany Simulator

A plant growth simulator using hormone-driven meristem mechanics. Plants grow from a seed node that spawns one shoot apical meristem (upward) and one root apical meristem (downward). Growth is governed by auxin/cytokinin signaling.

## Build & Run

```bash
# Build (cmake is at /usr/local/bin/cmake)
/usr/local/bin/cmake --build build

# Run realtime viewer (must run from project root for shader path)
./build/botany_realtime [--color auxin|cytokinin|sugar|type]

# Run tests
./build/botany_tests
```

## Architecture

### Engine (`src/engine/`)
- **genome.h** - `Genome` struct with all tunable parameters + `default_genome()`
- **node.h** - `Node` (position, radius, parent/children, hormone levels, sugar, leaf_size), `Meristem` base class, `MeristemType` enum. `NodeType`: STEM, ROOT, LEAF
- **meristem_types.h** - Concrete meristem subclasses: `ShootApicalMeristem`, `ShootAxillaryMeristem`, `RootApicalMeristem`, `RootAxillaryMeristem`
- **meristem.cpp** - All meristem tick logic and `tick_meristems()` dispatch
- **hormone.cpp** - `transport_auxin()` and `transport_cytokinin()` - per-tick hormone reset + transport
- **sugar.h/cpp** - `transport_sugar()` — sugar production, gradient-based diffusion, maintenance consumption
- **world_params.h** - `WorldParams` struct (light level, diffusion iterations) — non-genetic sim parameters
- **plant.h/cpp** - `Plant` class owns all nodes/meristems, has root meristem cap (100)
- **engine.h/cpp** - `Engine` orchestrates tick order: hormones, sugar, then meristems

### Renderer (`src/renderer/`)
- OpenGL 4.1 core profile, GLFW window, orbit camera
- Draws cylinders between parent-child nodes, leaves as quads
- Color modes: default (brown), chemical heatmap (auxin/cytokinin/sugar), type (green=shoot, orange=root)

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

## Sugar Model

Sugar **persists across ticks** (NOT reset like hormones). Three phases per tick:

1. **Production** — LEAF nodes produce: `sugar += light_level * leaf_size * sugar_production_rate`
2. **Diffusion** — Gradient-based bidirectional flow through all node connections:
   - Transport capacity = `min_radius^2 * PI * sugar_transport_conductance` (thicker = more)
   - LEAF connections use baseline capacity (leaf radius is 0)
   - Runs `sugar_diffusion_iterations` passes per tick (WorldParams, default 5)
   - Compute-then-apply: all flows calculated first, then applied atomically (order-independent)
3. **Consumption** — Every node deducts maintenance cost. `sugar = max(0, sugar - cost)`:
   - LEAF: `sugar_maintenance_leaf * leaf_size`
   - STEM: `sugar_maintenance_stem * radius`
   - ROOT: `sugar_maintenance_root * radius`
   - Active meristem tips: `+ sugar_maintenance_meristem`

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
