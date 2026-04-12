# Mass & Stress System Design

## Overview

Add gravitational physics to the plant simulator. Each above-ground node tracks its mass and the accumulated mass of its subtree, computes structural stress from gravity and lever arm distance, and responds via a new stress hormone. Branches droop under moderate load and snap under severe load. Plants evolve wood density, flexibility, and stress hormone sensitivity to develop diverse structural strategies.

## Mass Computation

Each node computes its own mass and accumulates subtree totals during the existing DFS tick (seed→tip). Values read from children are one tick stale — acceptable for a slowly changing system.

### Per-Node Fields

| Field | Type | Stored? | Description |
|---|---|---|---|
| `self_mass` | float | No — local variable in tick | This node's mass (volume × density for stems, area × density for leaves, constant for meristems) |
| `total_mass` | float | Yes — read by parent next tick | self_mass + Σ child.total_mass |
| `mass_moment` | vec3 | Yes — read by parent next tick | self_mass × position + Σ child.mass_moment |
| `stress` | float | Yes — used by growth logic and renderer | Computed torque / cross-section |

`total_mass` and `mass_moment` must be node fields because parents read them (one tick stale). `stress` is stored so growth logic and the renderer can access it.

### Mass Formulas

**Stem/Root nodes:**
```
volume = π × radius² × length
self_mass = volume × genome.wood_density
```

**Leaf nodes:**
```
area = leaf_size²
self_mass = area × leaf_mass_density    (WorldParams constant, not evolvable)
```

**Meristem nodes:**
```
self_mass = meristem_mass               (WorldParams constant, small fixed value)
```

### Accumulation (during DFS tick, before growth)

```
node.total_mass  = self_mass + Σ child.total_mass
node.mass_moment = self_mass × node.position + Σ child.mass_moment
```

Each node only reads its direct children's values — O(number of direct children), not O(subtree). Children's values are from last tick.

## Stress Computation

Stress measures how hard gravity is trying to break this node, accounting for lever arm distance via center of mass.

```
child_mass = total_mass - self_mass
child_com  = (mass_moment - self_mass × position) / child_mass
lever_arm  = horizontal_distance(position, child_com)
torque     = child_mass × gravity × lever_arm
cross_section = π × radius²
stress     = torque / cross_section
```

Key property: a heavy node far from the branch point creates proportionally more stress than the same mass nearby, because it pulls the center of mass outward, increasing the lever arm.

**Ground support:** If a node's Y position is below `ground_support_height` (WorldParams, default 0.5 dm), its stress is zeroed out. The ground bears the weight. This enables creeping/vine strategies to emerge naturally.

**Scope:** Stress only applies to above-ground nodes (STEM, LEAF, and above-ground meristems). ROOT nodes skip stress computation entirely.

## Stress Response

### Droop

When stress exceeds the droop threshold, the node's offset rotates toward gravity. The droop rate is proportional to how far past the threshold the stress is.

```
droop_threshold = wood_density × break_strength_factor × wood_flexibility
```

`wood_flexibility` is a genome parameter (0.1–1.0). High flexibility means droop starts at a lower fraction of break stress — the branch bends significantly before snapping. Low flexibility means almost no droop before sudden breakage.

Droop implementation: rotate the offset vector toward (0, -1, 0) by a small angle per tick, proportional to `(stress - droop_threshold) / droop_threshold`. Cap the rotation rate so branches don't teleport downward.

### Break

When stress exceeds the break threshold, the branch snaps. Call `plant.remove_subtree(node)` — the node and all its children are removed.

```
break_stress = wood_density × break_strength_factor
```

`break_strength_factor` is a WorldParams constant (same for all plants — material physics, not genetics). Dense wood is harder to break because `break_stress` scales with `wood_density`.

### Sugar Cost Scaling

Growth sugar costs scale with wood density. In the existing WorldParams sugar costs (shoot growth, thickening, elongation), multiply by `wood_density / reference_density`:

```
effective_cost = base_cost × (genome.wood_density / world.reference_wood_density)
```

`reference_wood_density` is a WorldParams constant representing the density at which the existing sugar costs were calibrated. Denser wood costs proportionally more sugar to produce.

## Stress Hormone

A new chemical signal (`ChemicalID::Stress`) produced by nodes under mechanical stress. Uses the existing hormone transport infrastructure.

### Production

Each above-ground node produces stress hormone proportional to its computed stress:

```
production = stress × genome.stress_hormone_production_rate
```

### Transport

Uses the existing generic transport system (`Node::transport_chemicals`). Registered in the chemical registry with its own transport rate, directional bias, and decay rate — all evolvable genome parameters.

Default transport characteristics: seedward bias (-0.7), moderate transport rate, high decay rate (local signal that doesn't travel far).

### Effects on Growth

Three separate evolvable sensitivities control how the stress hormone modifies growth:

**Thickening boost:** Stress hormone increases thickening rate in `StemNode::thicken()`.
```
thickening_multiplier = 1.0 + stress_hormone_level × genome.stress_thickening_boost
```

**Elongation inhibition:** Stress hormone suppresses internode elongation in `StemNode::elongate()`.
```
stress_inhibit = max(0, 1.0 - stress_hormone_level × genome.stress_elongation_inhibition)
effective_rate *= stress_inhibit
```

**Gravitropism boost:** Stress hormone increases the tendency for new growth to angle upward in meristem growth direction calculations.
```
vertical_pull = base_gravitropism + stress_hormone_level × genome.stress_gravitropism_boost
```

## Genome Parameters

### New Fields

| Parameter | Default | Bounds | Unit | Role |
|---|---|---|---|---|
| `wood_density` | 50.0 | [10, 200] | g/dm³ | Mass per volume. Heavier = stronger but more stress. Scales sugar growth costs. |
| `wood_flexibility` | 0.5 | [0.1, 1.0] | ratio | Droop threshold as fraction of break threshold. High = bendy, low = rigid. |
| `stress_hormone_production_rate` | 0.1 | [0.0, 1.0] | hormone/stress | Signal output per unit stress. Low = stoic, high = reactive. |
| `stress_hormone_transport_rate` | 0.15 | [0.01, 0.5] | fraction/tick | Diffusion speed to neighbors. |
| `stress_hormone_directional_bias` | -0.7 | [-1.0, 1.0] | bias | Seedward bias. Negative = flows toward trunk. |
| `stress_hormone_decay_rate` | 0.2 | [0.01, 0.5] | fraction/tick | Signal fade rate. High = very local. |
| `stress_thickening_boost` | 1.0 | [0.0, 5.0] | multiplier/hormone | How much stress hormone boosts thickening. |
| `stress_elongation_inhibition` | 1.0 | [0.0, 5.0] | multiplier/hormone | How much stress hormone suppresses elongation. |
| `stress_gravitropism_boost` | 0.5 | [0.0, 5.0] | multiplier/hormone | How much stress hormone pulls growth vertical. |

### Linkage Group

All 9 stress parameters form one linkage group (`stress`) for crossover: wood properties and stress hormone response are inherited together.

### Evolution Registration

All parameters registered in `genome_bridge.cpp` with proportional mutation strength (`mutation_pct × range`), same as existing genes.

## WorldParams Constants

| Parameter | Default | Unit | Role |
|---|---|---|---|
| `gravity` | 9.81 | m/s² | Gravitational acceleration for torque computation |
| `break_strength_factor` | 5.0 | stress units / (g/dm³) | How strong wood is per unit density |
| `reference_wood_density` | 50.0 | g/dm³ | Density at which existing sugar costs are calibrated |
| `leaf_mass_density` | 5.0 | g/dm² | Mass per unit leaf area |
| `meristem_mass` | 0.1 | g | Fixed mass for meristem tips |
| `ground_support_height` | 0.5 | dm | Below this Y, stress is zeroed (ground support) |
| `droop_rate` | 0.01 | radians/tick | Max angular droop per tick when overstressed |

## Evolutionary Strategies This Enables

| Strategy | Key genome values | Behavior |
|---|---|---|
| Thick trunk builder | High thickening_boost, moderate production | Grows wide, invests sugar in structural wood |
| Compact vertical | High elongation_inhibition + gravitropism_boost | Stays short and upright |
| Late bloomer | High sensitivity across all three effects | Only spreads wide once trunk can handle stress |
| Ground creeper | Low production + low sensitivity | Branches droop to ground, get supported, spread sideways |
| Light and flexible | Low wood_density, high wood_flexibility | Bends under load but doesn't break, cheap to grow |
| Dense and rigid | High wood_density, low wood_flexibility | Strong but expensive, snaps suddenly if overloaded |

## Files Changed / Created

**Modified files:**
- `src/engine/genome.h` — 9 new genome fields + defaults
- `src/engine/world_params.h` — 7 new WorldParams constants
- `src/engine/node/node.h` — add `total_mass`, `mass_moment`, `stress` fields
- `src/engine/node/node.cpp` — mass/stress computation in `tick()`
- `src/engine/node/stem_node.cpp` — stress hormone effects on thickening/elongation
- `src/engine/node/leaf_node.cpp` — leaf self_mass computation
- `src/engine/node/meristems/shoot_apical.cpp` — gravitropism boost from stress hormone
- `src/engine/node/meristems/helpers.h` — droop logic
- `src/engine/chemical/chemical.h` — add `ChemicalID::Stress`
- `src/engine/chemical/chemical_registry.h` — register stress hormone transport config
- `src/engine/plant.cpp` — droop and break checks after stress computation
- `src/evolution/genome_bridge.cpp` — register 9 new genes + stress linkage group
- `src/renderer/renderer.cpp` — stress heatmap color mode
- `src/app_realtime.cpp` — stress overlay button
- `tests/test_evolution.cpp` — update round-trip test for new genome fields

**No new files** — stress integrates into the existing node tick pipeline and hormone transport system.
