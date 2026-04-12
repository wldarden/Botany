# Evolution System Design

## Overview

Add a genetic algorithm system that evolves plant genomes using the external `evolve` library (`/Users/wldarden/repos/evolve`). A new `botany_evolve` app runs populations of plants headlessly, scores them on multi-objective fitness, and uses `evolve::StructuredGenome` to mutate/crossover the population across generations. The best plant is rendered live in the viewport.

## Architecture

Four new components, all under `src/evolution/`:

```
app_evolve
├── EvolutionRunner (owns population, orchestrates generations)
│   ├── genome_bridge (botany::Genome ↔ evolve::StructuredGenome)
│   └── fitness (runs Engine per-plant, computes score)
├── botany_renderer (renders best plant)
└── imgui (config panel + stats)
```

The engine and renderer are unchanged except one small addition: `Plant` gains a `total_sugar_produced` accumulator incremented during the existing sugar production pass.

## Genome Bridge (`src/evolution/genome_bridge.h/cpp`)

### Conversion Functions

```cpp
evolve::StructuredGenome to_structured(const botany::Genome& g);
botany::Genome from_structured(const evolve::StructuredGenome& sg);
```

A template `StructuredGenome` is built once at startup from `default_genome()`. Each field of `botany::Genome` becomes a named gene with a `GeneMutationConfig` tuned to that field's scale and valid range.

### Per-Gene Mutation Configs

| Gene group | Example fields | Mutation strength | Bounds |
|---|---|---|---|
| Hormone rates | `auxin_production_rate`, `auxin_transport_rate` | 0.05 | [0.01, 2.0] |
| Hormone bias | `auxin_directional_bias` | 0.1 | [-1.0, 1.0] |
| Thresholds | `auxin_threshold`, `ethylene_abscission_threshold` | 0.03 | [0.01, 1.0] |
| Growth rates | `growth_rate`, `leaf_growth_rate` | 0.001 | [0.001, 0.05] |
| Lengths | `max_internode_length`, `max_leaf_size` | 0.05 | [0.05, 3.0] |
| Angles | `branch_angle`, `growth_noise` | 0.1 | [0.05, 1.57] |
| Sugar params | `sugar_production_rate`, maintenance rates | 0.002 | [0.001, 1.0] |
| Tick counts | `internode_maturation_ticks`, `senescence_duration` | 10.0 | [12, 2000] |

### Linkage Groups

During crossover, all genes in a linkage group are inherited from the same parent:

- `auxin` — auxin_production_rate, auxin_transport_rate, auxin_directional_bias, auxin_decay_rate, auxin_threshold
- `cytokinin` — cytokinin_production_rate, cytokinin_transport_rate, cytokinin_directional_bias, cytokinin_decay_rate, cytokinin_threshold
- `shoot_growth` — growth_rate, max_internode_length, min_internode_length, branch_angle, thickening_rate, internode_elongation_rate, internode_maturation_ticks
- `root_growth` — root_growth_rate, root_max_internode_length, root_min_internode_length, root_branch_angle, root_internode_elongation_rate, root_internode_maturation_ticks, root_gravitropism_strength, root_gravitropism_depth
- `sugar_economy` — sugar_production_rate, sugar_transport_conductance, sugar_maintenance_leaf, sugar_maintenance_stem, sugar_maintenance_root, sugar_maintenance_meristem, seed_sugar, sugar_storage_density_wood, sugar_storage_density_leaf, sugar_cap_minimum, sugar_cap_meristem, sugar_save_shoot, sugar_save_root, sugar_save_stem, sugar_activation_shoot, sugar_activation_root
- `gibberellin` — ga_production_rate, ga_leaf_age_max, ga_elongation_sensitivity, ga_length_sensitivity, ga_transport_rate, ga_directional_bias, ga_decay_rate
- `ethylene` — ethylene_starvation_rate, ethylene_shade_rate, ethylene_shade_threshold, ethylene_age_rate, ethylene_age_onset, ethylene_crowding_rate, ethylene_crowding_radius, ethylene_diffusion_radius, ethylene_abscission_threshold, ethylene_elongation_inhibition, senescence_duration
- `geometry` — max_leaf_size, leaf_growth_rate, leaf_bud_size, initial_radius, root_initial_radius, tip_offset, growth_noise, leaf_phototropism_rate

## Fitness Evaluator (`src/evolution/fitness.h/cpp`)

### Simulation

Each plant is evaluated independently: one `Engine` + one `Plant`, run tick-by-tick up to `max_ticks` (default 17,520 = 2 years). Early termination if the plant has zero active meristems (no growth potential).

### Stats Collected

```cpp
struct PlantStats {
    uint32_t survival_ticks;
    uint32_t node_count;
    uint32_t leaf_count;
    float total_sugar_produced;
    float height;
    float crown_ratio;       // canopy_width / height
    uint32_t branch_depth;   // max branching generations
    float leaf_height_spread; // std dev of leaf Y positions
};
```

Computed at evaluation end from one pass over the node tree:
- `height` — max node Y minus seed Y
- `crown_ratio` — max horizontal bounding box extent / height (returns 0.0 if height < 0.01 to avoid div-by-zero)
- `branch_depth` — max depth of branch-from-branch chain (trunk = 0, first branch = 1, sub-branch = 2, etc.)
- `leaf_height_spread` — standard deviation of living leaf Y positions

### Fitness Computation

```cpp
struct FitnessWeights {
    float survival    = 1.0f;
    float biomass     = 1.0f;
    float sugar       = 1.0f;
    float leaves      = 1.0f;
    float height      = 1.0f;
    float crown_ratio = 1.0f;
    float branch_depth = 1.0f;
    float leaf_spread  = 1.0f;
};
```

Each objective is normalized to [0, 1] by dividing by the max value seen in the current generation. Final fitness is the weighted sum:

```
fitness = sum(weight_i * (stat_i / gen_max_i))
```

Normalizing per-generation ensures no single objective dominates by virtue of larger raw numbers.

```cpp
PlantStats evaluate_plant(const Genome& genome, const WorldParams& world, uint32_t max_ticks);
float compute_fitness(const PlantStats& stats, const PlantStats& gen_max, const FitnessWeights& weights);
```

### Shape-Aware Objectives Rationale

Raw width and height reward extremes (single tall stem, flat fan of branches). The revised objectives prevent gaming:

| Objective | Rewards | Prevents |
|---|---|---|
| survival_ticks | Longevity | Fragile genomes |
| node_count | Biomass | Tiny plants |
| leaf_count | Canopy | All-wood, no leaves |
| total_sugar_produced | Photosynthetic efficiency | Parasitic layouts |
| height | Tall growth | Staying short |
| crown_ratio | Proportional shape | Poles and pancakes |
| branch_depth | Recursive branching | Single-stem hacks |
| leaf_height_spread | Distributed canopy | All leaves in one cluster |

## Evolution Runner (`src/evolution/evolution_runner.h/cpp`)

### Configuration

```cpp
struct EvolutionConfig {
    uint32_t population_size = 100;
    uint32_t max_ticks = 17520;
    uint32_t num_threads = 4;
    FitnessWeights weights;
};
```

### Generation Loop

1. **Randomize environment** — roll new `WorldParams` for this generation. Randomize `light_level` (uniform in [0.5, 1.0]) and `light_direction` (tilt the Y-up vector by a random angle up to 30 degrees). All other WorldParams fields stay at defaults. All plants in the generation share the same conditions.
2. **Evaluate in parallel** — partition population across `num_threads` worker threads. Each thread creates its own `Engine` instance. Each worker calls `evaluate_plant()` and stores the `PlantStats`. No shared mutable state between threads.
3. **Normalize and score** — find per-generation max for each objective, call `compute_fitness()` for every individual, set fitness on the `StructuredGenome`.
4. **Record best** — stash the best genome and stats for the renderer.
5. **Evolve** — call the evolve library's `mutate()`/`crossover()` on the StructuredGenome population for the next generation.

### Threading Model

`std::vector<std::thread>` with work partitioned by index range. Each thread owns its own `Engine` — no shared mutable state. The runner joins all threads before proceeding to the scoring step.

### Best-Plant Rendering

After each generation completes, the runner converts the best genome back to `botany::Genome` via `from_structured()`, re-simulates it in a dedicated `Engine` on the main thread, and hands that engine to the renderer. The renderer draws it while the next generation evaluates in background threads.

### State Exposed to UI

```cpp
uint32_t generation;
float best_fitness;
PlantStats best_stats;
std::vector<float> fitness_history;  // best fitness per generation, for line plot
```

## App: botany_evolve (`src/app_evolve.cpp`)

GLFW + OpenGL + ImGui, same stack as the other apps. Single window with 3D viewport and ImGui side panel.

### Startup Flow

App opens with config panel visible, evolution paused. User sets parameters, clicks "Start". Evolution runs in background. User can pause/resume anytime.

### ImGui Side Panel

**Config section** (editable while paused):
- Population size, max ticks per plant, thread count
- 8 fitness weight sliders
- WorldParams variation range (light level min/max)

**Stats section** (live during evolution):
- Generation counter, best fitness, evaluation throughput
- Best plant stats table (all 8 objectives, raw and normalized values)
- Fitness history line plot (`ImGui::PlotLines`)

**Controls:**
- Start / Pause / Reset buttons
- "Export Best" button — writes a `.genome` text file (one field per line, `name=value`). Loading genomes into other apps is out of scope for this spec.

### Viewport

Renders the best plant from the last completed generation using the existing `Renderer`. Orbit camera, same controls as `app_realtime`. Updates when a new generation completes.

## Engine Change

One addition to `Plant`: a `float total_sugar_produced_` accumulator, incremented during the existing sugar production pass in `sugar.cpp`. Exposed via `Plant::total_sugar_produced() const`.

## Build Integration

CMakeLists.txt additions:
- Add `evolve` library via `add_subdirectory(/Users/wldarden/repos/evolve)`
- New `botany_evolution` static library target with the 3 source files
- New `botany_evolve` executable linking `botany_engine`, `botany_evolution`, `botany_renderer`, `evolve`, `imgui`

## Files Changed / Created

**New files:**
- `src/evolution/genome_bridge.h` / `src/evolution/genome_bridge.cpp`
- `src/evolution/fitness.h` / `src/evolution/fitness.cpp`
- `src/evolution/evolution_runner.h` / `src/evolution/evolution_runner.cpp`
- `src/app_evolve.cpp`
- `tests/test_evolution.cpp`

**Modified files:**
- `src/engine/plant.h` / `src/engine/plant.cpp` — add `total_sugar_produced_` accumulator
- `src/engine/sugar.cpp` — increment accumulator during production pass
- `CMakeLists.txt` — add evolve dependency, new library target, new executable, new test file
