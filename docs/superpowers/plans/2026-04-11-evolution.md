# Evolution System Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a genetic algorithm that evolves plant genomes using the `evolve` library's `StructuredGenome`, with a new `botany_evolve` app showing live stats and the best plant.

**Architecture:** New `src/evolution/` folder with genome bridge, fitness evaluator, and evolution runner. One small engine change (sugar production accumulator on Plant). New `app_evolve.cpp` with ImGui config + renderer. Integrates `/Users/wldarden/repos/evolve` via `add_subdirectory`.

**Tech Stack:** C++17, evolve library (StructuredGenome, mutate, crossover), GLFW, OpenGL 4.1, ImGui, Catch2, std::thread

**Build:** `/usr/local/bin/cmake --build build` (always rebuild before running tests)

**Test:** `./build/botany_tests`

---

## File Structure

**New files:**
- `src/evolution/genome_bridge.h` / `.cpp` — `build_genome_template()`, `to_structured()`, `from_structured()`
- `src/evolution/fitness.h` / `.cpp` — `PlantStats`, `FitnessWeights`, `evaluate_plant()`, `compute_fitness()`
- `src/evolution/evolution_runner.h` / `.cpp` — `EvolutionRunner` class with threaded generation loop
- `src/app_evolve.cpp` — ImGui config panel + 3D viewport showing best plant
- `tests/test_evolution.cpp` — tests for genome bridge, fitness, runner

**Modified files:**
- `src/engine/plant.h` / `src/engine/plant.cpp` — add `total_sugar_produced_` accumulator
- `src/engine/node/leaf_node.cpp` — increment accumulator in `grow()`
- `CMakeLists.txt` — add evolve dep, `botany_evolution` lib, `botany_evolve` exe, test file

---

### Task 1: CMake — add evolve library and scaffold new targets

**Files:**
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Add evolve library and new targets to CMakeLists.txt**

Add after the existing `FetchContent_MakeAvailable` block and library definitions. Add `evolve` via `add_subdirectory`, create the `botany_evolution` static library (empty source files for now), the `botany_evolve` executable, and add the test file:

```cmake
# Evolution library (external)
add_subdirectory(/Users/wldarden/repos/evolve ${CMAKE_BINARY_DIR}/evolve)

# Evolution library (project)
add_library(botany_evolution STATIC
    src/evolution/genome_bridge.cpp
    src/evolution/fitness.cpp
    src/evolution/evolution_runner.cpp
)
target_include_directories(botany_evolution PUBLIC src)
target_link_libraries(botany_evolution PUBLIC botany_engine evolve)

# Evolution app
add_executable(botany_evolve src/app_evolve.cpp)
target_link_libraries(botany_evolve PRIVATE botany_engine botany_evolution botany_renderer imgui)
```

Add `tests/test_evolution.cpp` to the existing `botany_tests` target's source list.

- [ ] **Step 2: Create stub source files so the build succeeds**

Create minimal stubs for all new source files:

`src/evolution/genome_bridge.h`:
```cpp
#pragma once

namespace botany {} // namespace botany
```

`src/evolution/genome_bridge.cpp`:
```cpp
#include "evolution/genome_bridge.h"
```

`src/evolution/fitness.h`:
```cpp
#pragma once

namespace botany {} // namespace botany
```

`src/evolution/fitness.cpp`:
```cpp
#include "evolution/fitness.h"
```

`src/evolution/evolution_runner.h`:
```cpp
#pragma once

namespace botany {} // namespace botany
```

`src/evolution/evolution_runner.cpp`:
```cpp
#include "evolution/evolution_runner.h"
```

`src/app_evolve.cpp`:
```cpp
int main() { return 0; }
```

`tests/test_evolution.cpp`:
```cpp
#include <catch2/catch_test_macros.hpp>

TEST_CASE("evolution stub", "[evolution]") {
    REQUIRE(true);
}
```

- [ ] **Step 3: Build and verify**

Run: `/usr/local/bin/cmake --build build`
Expected: Build succeeds with no errors. The evolve library compiles and links.

- [ ] **Step 4: Run tests to verify stub passes**

Run: `./build/botany_tests "[evolution]"`
Expected: 1 test passes.

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt src/evolution/ src/app_evolve.cpp tests/test_evolution.cpp
git commit -m "feat(evolution): scaffold cmake targets and stub files for evolution system"
```

---

### Task 2: Sugar production accumulator on Plant

**Files:**
- Modify: `src/engine/plant.h`
- Modify: `src/engine/node/leaf_node.cpp`
- Test: `tests/test_evolution.cpp`

- [ ] **Step 1: Write the failing test**

In `tests/test_evolution.cpp`, add:

```cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/engine.h"
#include "engine/genome.h"

using Catch::Matchers::WithinAbs;

TEST_CASE("Plant tracks total sugar produced", "[evolution]") {
    botany::Engine engine;
    botany::Genome g = botany::default_genome();
    engine.create_plant(g, glm::vec3(0.0f));

    // Run enough ticks for leaves to grow and photosynthesize
    for (int i = 0; i < 200; i++) {
        engine.tick();
    }

    float produced = engine.get_plant(0).total_sugar_produced();
    REQUIRE(produced > 0.0f);
}
```

- [ ] **Step 2: Build and verify it fails**

Run: `/usr/local/bin/cmake --build build`
Expected: Compile error — `total_sugar_produced()` not declared in `Plant`.

- [ ] **Step 3: Add accumulator to Plant**

In `src/engine/plant.h`, add to the public section after `root_meristems_at_cap()`:

```cpp
float total_sugar_produced() const { return total_sugar_produced_; }
void add_sugar_produced(float amount) { total_sugar_produced_ += amount; }
```

Add to the private section after `pending_removals_`:

```cpp
float total_sugar_produced_ = 0.0f;
```

- [ ] **Step 4: Increment accumulator in LeafNode::grow()**

In `src/engine/node/leaf_node.cpp`, in the `grow()` method, wrap the `photosynthesize` call to track the delta:

Replace:
```cpp
    photosynthesize(g, world);
```

With:
```cpp
    float sugar_before = chemical(ChemicalID::Sugar);
    photosynthesize(g, world);
    float delta = chemical(ChemicalID::Sugar) - sugar_before;
    if (delta > 0.0f) plant.add_sugar_produced(delta);
```

- [ ] **Step 5: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[evolution]"`
Expected: All tests pass. The sugar accumulator tracks production.

- [ ] **Step 6: Commit**

```bash
git add src/engine/plant.h src/engine/node/leaf_node.cpp tests/test_evolution.cpp
git commit -m "feat(plant): add total_sugar_produced accumulator for fitness tracking"
```

---

### Task 3: Genome Bridge — template, to_structured, from_structured

**Files:**
- Modify: `src/evolution/genome_bridge.h`
- Modify: `src/evolution/genome_bridge.cpp`
- Test: `tests/test_evolution.cpp`

- [ ] **Step 1: Write the failing test**

In `tests/test_evolution.cpp`, add:

```cpp
#include "evolution/genome_bridge.h"

TEST_CASE("Genome round-trips through StructuredGenome", "[evolution]") {
    botany::Genome original = botany::default_genome();
    auto sg = botany::to_structured(original);
    botany::Genome restored = botany::from_structured(sg);

    // Spot-check representative fields from different groups
    REQUIRE_THAT(restored.auxin_production_rate,
                 WithinAbs(original.auxin_production_rate, 1e-6));
    REQUIRE_THAT(restored.branch_angle,
                 WithinAbs(original.branch_angle, 1e-6));
    REQUIRE_THAT(restored.sugar_production_rate,
                 WithinAbs(original.sugar_production_rate, 1e-6));
    REQUIRE_THAT(restored.ga_elongation_sensitivity,
                 WithinAbs(original.ga_elongation_sensitivity, 1e-6));
    REQUIRE_THAT(restored.ethylene_abscission_threshold,
                 WithinAbs(original.ethylene_abscission_threshold, 1e-6));
    REQUIRE(restored.internode_maturation_ticks == original.internode_maturation_ticks);
    REQUIRE(restored.senescence_duration == original.senescence_duration);
}

TEST_CASE("Genome template has linkage groups", "[evolution]") {
    auto tmpl = botany::build_genome_template(botany::default_genome());
    auto& groups = tmpl.linkage_groups();
    REQUIRE(groups.size() == 8);

    // Verify auxin group has 5 genes
    bool found_auxin = false;
    for (auto& g : groups) {
        if (g.name == "auxin") {
            found_auxin = true;
            REQUIRE(g.gene_tags.size() == 5);
        }
    }
    REQUIRE(found_auxin);
}
```

- [ ] **Step 2: Build and verify it fails**

Run: `/usr/local/bin/cmake --build build`
Expected: Compile error — `to_structured`, `from_structured`, `build_genome_template` not declared.

- [ ] **Step 3: Write genome_bridge.h**

```cpp
#pragma once

#include "engine/genome.h"
#include <evolve/structured_genome.h>

namespace botany {

// Build a StructuredGenome template with all gene definitions,
// mutation configs, and linkage groups. Values come from the given genome.
evolve::StructuredGenome build_genome_template(const Genome& g);

// Convert a botany Genome to a StructuredGenome (uses build_genome_template internally).
evolve::StructuredGenome to_structured(const Genome& g);

// Convert a StructuredGenome back to a botany Genome.
Genome from_structured(const evolve::StructuredGenome& sg);

} // namespace botany
```

- [ ] **Step 4: Write genome_bridge.cpp**

```cpp
#include "evolution/genome_bridge.h"

namespace botany {

// Helper: register a single-value gene with mutation config
static void reg(evolve::StructuredGenome& sg, const std::string& tag, float val,
                float rate, float strength, float min_val, float max_val) {
    sg.add_gene({tag, {val}, {rate, strength, min_val, max_val}});
}

evolve::StructuredGenome build_genome_template(const Genome& g) {
    evolve::StructuredGenome sg;

    // --- Auxin (linkage group) ---
    reg(sg, "auxin_production_rate",    g.auxin_production_rate,    0.1f, 0.05f, 0.01f, 2.0f);
    reg(sg, "auxin_transport_rate",     g.auxin_transport_rate,     0.1f, 0.05f, 0.01f, 2.0f);
    reg(sg, "auxin_directional_bias",   g.auxin_directional_bias,   0.1f, 0.1f, -1.0f, 1.0f);
    reg(sg, "auxin_decay_rate",         g.auxin_decay_rate,         0.1f, 0.02f, 0.001f, 0.5f);
    reg(sg, "auxin_threshold",          g.auxin_threshold,          0.1f, 0.03f, 0.01f, 1.0f);

    // --- Cytokinin (linkage group) ---
    reg(sg, "cytokinin_production_rate",    g.cytokinin_production_rate,    0.1f, 0.05f, 0.01f, 2.0f);
    reg(sg, "cytokinin_transport_rate",     g.cytokinin_transport_rate,     0.1f, 0.05f, 0.01f, 2.0f);
    reg(sg, "cytokinin_directional_bias",   g.cytokinin_directional_bias,   0.1f, 0.1f, -1.0f, 1.0f);
    reg(sg, "cytokinin_decay_rate",         g.cytokinin_decay_rate,         0.1f, 0.02f, 0.001f, 0.5f);
    reg(sg, "cytokinin_threshold",          g.cytokinin_threshold,          0.1f, 0.03f, 0.01f, 1.0f);

    // --- Shoot growth (linkage group) ---
    reg(sg, "growth_rate",               g.growth_rate,               0.1f, 0.001f, 0.001f, 0.05f);
    reg(sg, "max_internode_length",      g.max_internode_length,      0.1f, 0.05f, 0.05f, 3.0f);
    reg(sg, "min_internode_length",      g.min_internode_length,      0.1f, 0.03f, 0.01f, 1.0f);
    reg(sg, "branch_angle",             g.branch_angle,             0.1f, 0.1f, 0.05f, 1.57f);
    reg(sg, "thickening_rate",          g.thickening_rate,          0.1f, 0.00002f, 0.00001f, 0.001f);
    reg(sg, "internode_elongation_rate", g.internode_elongation_rate, 0.1f, 0.001f, 0.0005f, 0.02f);
    reg(sg, "internode_maturation_ticks", static_cast<float>(g.internode_maturation_ticks), 0.1f, 10.0f, 12.0f, 500.0f);

    // --- Root growth (linkage group) ---
    reg(sg, "root_growth_rate",               g.root_growth_rate,               0.1f, 0.001f, 0.001f, 0.05f);
    reg(sg, "root_max_internode_length",      g.root_max_internode_length,      0.1f, 0.05f, 0.05f, 3.0f);
    reg(sg, "root_min_internode_length",      g.root_min_internode_length,      0.1f, 0.03f, 0.01f, 1.0f);
    reg(sg, "root_branch_angle",             g.root_branch_angle,             0.1f, 0.05f, 0.05f, 1.57f);
    reg(sg, "root_internode_elongation_rate", g.root_internode_elongation_rate, 0.1f, 0.001f, 0.0005f, 0.02f);
    reg(sg, "root_internode_maturation_ticks", static_cast<float>(g.root_internode_maturation_ticks), 0.1f, 10.0f, 12.0f, 500.0f);
    reg(sg, "root_gravitropism_strength",    g.root_gravitropism_strength,    0.1f, 0.2f, 0.1f, 10.0f);
    reg(sg, "root_gravitropism_depth",       g.root_gravitropism_depth,       0.1f, 0.1f, 0.1f, 5.0f);

    // --- Geometry (linkage group) ---
    reg(sg, "max_leaf_size",           g.max_leaf_size,           0.1f, 0.03f, 0.05f, 1.0f);
    reg(sg, "leaf_growth_rate",        g.leaf_growth_rate,        0.1f, 0.0005f, 0.0001f, 0.01f);
    reg(sg, "leaf_bud_size",           g.leaf_bud_size,           0.1f, 0.005f, 0.005f, 0.1f);
    reg(sg, "initial_radius",          g.initial_radius,          0.1f, 0.01f, 0.01f, 0.2f);
    reg(sg, "root_initial_radius",     g.root_initial_radius,     0.1f, 0.005f, 0.005f, 0.1f);
    reg(sg, "tip_offset",              g.tip_offset,              0.1f, 0.005f, 0.001f, 0.1f);
    reg(sg, "growth_noise",            g.growth_noise,            0.1f, 0.05f, 0.01f, 0.8f);
    reg(sg, "leaf_phototropism_rate",  g.leaf_phototropism_rate,  0.1f, 0.005f, 0.001f, 0.1f);

    // --- Sugar economy (linkage group) ---
    reg(sg, "sugar_production_rate",       g.sugar_production_rate,       0.1f, 0.002f, 0.001f, 0.1f);
    reg(sg, "sugar_transport_conductance", g.sugar_transport_conductance, 0.1f, 2.0f, 1.0f, 100.0f);
    reg(sg, "sugar_maintenance_leaf",      g.sugar_maintenance_leaf,      0.1f, 0.002f, 0.001f, 0.1f);
    reg(sg, "sugar_maintenance_stem",      g.sugar_maintenance_stem,      0.1f, 0.005f, 0.001f, 0.2f);
    reg(sg, "sugar_maintenance_root",      g.sugar_maintenance_root,      0.1f, 0.02f, 0.01f, 0.5f);
    reg(sg, "sugar_maintenance_meristem",  g.sugar_maintenance_meristem,  0.1f, 0.001f, 0.0001f, 0.01f);
    reg(sg, "seed_sugar",                  g.seed_sugar,                  0.1f, 5.0f, 5.0f, 200.0f);
    reg(sg, "sugar_storage_density_wood",  g.sugar_storage_density_wood,  0.1f, 50.0f, 50.0f, 2000.0f);
    reg(sg, "sugar_storage_density_leaf",  g.sugar_storage_density_leaf,  0.1f, 0.1f, 0.05f, 5.0f);
    reg(sg, "sugar_cap_minimum",           g.sugar_cap_minimum,           0.1f, 0.01f, 0.005f, 0.5f);
    reg(sg, "sugar_cap_meristem",          g.sugar_cap_meristem,          0.1f, 0.5f, 0.1f, 10.0f);
    reg(sg, "sugar_save_shoot",            g.sugar_save_shoot,            0.1f, 0.005f, 0.001f, 0.2f);
    reg(sg, "sugar_save_root",             g.sugar_save_root,             0.1f, 0.003f, 0.001f, 0.1f);
    reg(sg, "sugar_save_stem",             g.sugar_save_stem,             0.1f, 0.005f, 0.001f, 0.2f);
    reg(sg, "sugar_activation_shoot",      g.sugar_activation_shoot,      0.1f, 0.1f, 0.05f, 5.0f);
    reg(sg, "sugar_activation_root",       g.sugar_activation_root,       0.1f, 0.05f, 0.05f, 3.0f);

    // --- Gibberellin (linkage group) ---
    reg(sg, "ga_production_rate",       g.ga_production_rate,       0.1f, 0.05f, 0.01f, 2.0f);
    reg(sg, "ga_leaf_age_max",          static_cast<float>(g.ga_leaf_age_max), 0.1f, 20.0f, 24.0f, 2000.0f);
    reg(sg, "ga_elongation_sensitivity", g.ga_elongation_sensitivity, 0.1f, 0.2f, 0.1f, 5.0f);
    reg(sg, "ga_length_sensitivity",    g.ga_length_sensitivity,    0.1f, 0.2f, 0.1f, 5.0f);
    reg(sg, "ga_transport_rate",        g.ga_transport_rate,        0.1f, 0.05f, 0.01f, 1.0f);
    reg(sg, "ga_directional_bias",      g.ga_directional_bias,      0.1f, 0.1f, -1.0f, 1.0f);
    reg(sg, "ga_decay_rate",            g.ga_decay_rate,            0.1f, 0.03f, 0.01f, 0.5f);

    // --- Ethylene (linkage group) ---
    reg(sg, "ethylene_starvation_rate",       g.ethylene_starvation_rate,       0.1f, 0.05f, 0.01f, 2.0f);
    reg(sg, "ethylene_shade_rate",            g.ethylene_shade_rate,            0.1f, 0.05f, 0.01f, 2.0f);
    reg(sg, "ethylene_shade_threshold",       g.ethylene_shade_threshold,       0.1f, 0.05f, 0.05f, 1.0f);
    reg(sg, "ethylene_age_rate",              g.ethylene_age_rate,              0.1f, 0.01f, 0.001f, 0.5f);
    reg(sg, "ethylene_age_onset",             static_cast<float>(g.ethylene_age_onset), 0.1f, 50.0f, 100.0f, 5000.0f);
    reg(sg, "ethylene_crowding_rate",         g.ethylene_crowding_rate,         0.1f, 0.02f, 0.01f, 1.0f);
    reg(sg, "ethylene_crowding_radius",       g.ethylene_crowding_radius,       0.1f, 0.1f, 0.1f, 3.0f);
    reg(sg, "ethylene_diffusion_radius",      g.ethylene_diffusion_radius,      0.1f, 0.2f, 0.1f, 5.0f);
    reg(sg, "ethylene_abscission_threshold",  g.ethylene_abscission_threshold,  0.1f, 0.1f, 0.05f, 2.0f);
    reg(sg, "ethylene_elongation_inhibition", g.ethylene_elongation_inhibition, 0.1f, 0.2f, 0.1f, 5.0f);
    reg(sg, "senescence_duration",            static_cast<float>(g.senescence_duration), 0.1f, 10.0f, 12.0f, 500.0f);

    // --- Linkage groups ---
    sg.add_linkage_group({"auxin", {
        "auxin_production_rate", "auxin_transport_rate", "auxin_directional_bias",
        "auxin_decay_rate", "auxin_threshold"
    }});
    sg.add_linkage_group({"cytokinin", {
        "cytokinin_production_rate", "cytokinin_transport_rate", "cytokinin_directional_bias",
        "cytokinin_decay_rate", "cytokinin_threshold"
    }});
    sg.add_linkage_group({"shoot_growth", {
        "growth_rate", "max_internode_length", "min_internode_length", "branch_angle",
        "thickening_rate", "internode_elongation_rate", "internode_maturation_ticks"
    }});
    sg.add_linkage_group({"root_growth", {
        "root_growth_rate", "root_max_internode_length", "root_min_internode_length",
        "root_branch_angle", "root_internode_elongation_rate", "root_internode_maturation_ticks",
        "root_gravitropism_strength", "root_gravitropism_depth"
    }});
    sg.add_linkage_group({"geometry", {
        "max_leaf_size", "leaf_growth_rate", "leaf_bud_size", "initial_radius",
        "root_initial_radius", "tip_offset", "growth_noise", "leaf_phototropism_rate"
    }});
    sg.add_linkage_group({"sugar_economy", {
        "sugar_production_rate", "sugar_transport_conductance",
        "sugar_maintenance_leaf", "sugar_maintenance_stem", "sugar_maintenance_root",
        "sugar_maintenance_meristem", "seed_sugar",
        "sugar_storage_density_wood", "sugar_storage_density_leaf",
        "sugar_cap_minimum", "sugar_cap_meristem",
        "sugar_save_shoot", "sugar_save_root", "sugar_save_stem",
        "sugar_activation_shoot", "sugar_activation_root"
    }});
    sg.add_linkage_group({"gibberellin", {
        "ga_production_rate", "ga_leaf_age_max", "ga_elongation_sensitivity",
        "ga_length_sensitivity", "ga_transport_rate", "ga_directional_bias", "ga_decay_rate"
    }});
    sg.add_linkage_group({"ethylene", {
        "ethylene_starvation_rate", "ethylene_shade_rate", "ethylene_shade_threshold",
        "ethylene_age_rate", "ethylene_age_onset", "ethylene_crowding_rate",
        "ethylene_crowding_radius", "ethylene_diffusion_radius",
        "ethylene_abscission_threshold", "ethylene_elongation_inhibition", "senescence_duration"
    }});

    return sg;
}

evolve::StructuredGenome to_structured(const Genome& g) {
    return build_genome_template(g);
}

Genome from_structured(const evolve::StructuredGenome& sg) {
    Genome g{};
    g.auxin_production_rate    = sg.get("auxin_production_rate");
    g.auxin_transport_rate     = sg.get("auxin_transport_rate");
    g.auxin_directional_bias   = sg.get("auxin_directional_bias");
    g.auxin_decay_rate         = sg.get("auxin_decay_rate");
    g.auxin_threshold          = sg.get("auxin_threshold");

    g.cytokinin_production_rate    = sg.get("cytokinin_production_rate");
    g.cytokinin_transport_rate     = sg.get("cytokinin_transport_rate");
    g.cytokinin_directional_bias   = sg.get("cytokinin_directional_bias");
    g.cytokinin_decay_rate         = sg.get("cytokinin_decay_rate");
    g.cytokinin_threshold          = sg.get("cytokinin_threshold");

    g.growth_rate               = sg.get("growth_rate");
    g.max_internode_length      = sg.get("max_internode_length");
    g.min_internode_length      = sg.get("min_internode_length");
    g.branch_angle              = sg.get("branch_angle");
    g.thickening_rate           = sg.get("thickening_rate");
    g.internode_elongation_rate = sg.get("internode_elongation_rate");
    g.internode_maturation_ticks = static_cast<uint32_t>(sg.get("internode_maturation_ticks"));

    g.root_growth_rate               = sg.get("root_growth_rate");
    g.root_max_internode_length      = sg.get("root_max_internode_length");
    g.root_min_internode_length      = sg.get("root_min_internode_length");
    g.root_branch_angle              = sg.get("root_branch_angle");
    g.root_internode_elongation_rate = sg.get("root_internode_elongation_rate");
    g.root_internode_maturation_ticks = static_cast<uint32_t>(sg.get("root_internode_maturation_ticks"));
    g.root_gravitropism_strength     = sg.get("root_gravitropism_strength");
    g.root_gravitropism_depth        = sg.get("root_gravitropism_depth");

    g.max_leaf_size          = sg.get("max_leaf_size");
    g.leaf_growth_rate       = sg.get("leaf_growth_rate");
    g.leaf_bud_size          = sg.get("leaf_bud_size");
    g.initial_radius         = sg.get("initial_radius");
    g.root_initial_radius    = sg.get("root_initial_radius");
    g.tip_offset             = sg.get("tip_offset");
    g.growth_noise           = sg.get("growth_noise");
    g.leaf_phototropism_rate = sg.get("leaf_phototropism_rate");

    g.sugar_production_rate       = sg.get("sugar_production_rate");
    g.sugar_transport_conductance = sg.get("sugar_transport_conductance");
    g.sugar_maintenance_leaf      = sg.get("sugar_maintenance_leaf");
    g.sugar_maintenance_stem      = sg.get("sugar_maintenance_stem");
    g.sugar_maintenance_root      = sg.get("sugar_maintenance_root");
    g.sugar_maintenance_meristem  = sg.get("sugar_maintenance_meristem");
    g.seed_sugar                  = sg.get("seed_sugar");
    g.sugar_storage_density_wood  = sg.get("sugar_storage_density_wood");
    g.sugar_storage_density_leaf  = sg.get("sugar_storage_density_leaf");
    g.sugar_cap_minimum           = sg.get("sugar_cap_minimum");
    g.sugar_cap_meristem          = sg.get("sugar_cap_meristem");
    g.sugar_save_shoot            = sg.get("sugar_save_shoot");
    g.sugar_save_root             = sg.get("sugar_save_root");
    g.sugar_save_stem             = sg.get("sugar_save_stem");
    g.sugar_activation_shoot      = sg.get("sugar_activation_shoot");
    g.sugar_activation_root       = sg.get("sugar_activation_root");

    g.ga_production_rate       = sg.get("ga_production_rate");
    g.ga_leaf_age_max          = static_cast<uint32_t>(sg.get("ga_leaf_age_max"));
    g.ga_elongation_sensitivity = sg.get("ga_elongation_sensitivity");
    g.ga_length_sensitivity    = sg.get("ga_length_sensitivity");
    g.ga_transport_rate        = sg.get("ga_transport_rate");
    g.ga_directional_bias      = sg.get("ga_directional_bias");
    g.ga_decay_rate            = sg.get("ga_decay_rate");

    g.ethylene_starvation_rate       = sg.get("ethylene_starvation_rate");
    g.ethylene_shade_rate            = sg.get("ethylene_shade_rate");
    g.ethylene_shade_threshold       = sg.get("ethylene_shade_threshold");
    g.ethylene_age_rate              = sg.get("ethylene_age_rate");
    g.ethylene_age_onset             = static_cast<uint32_t>(sg.get("ethylene_age_onset"));
    g.ethylene_crowding_rate         = sg.get("ethylene_crowding_rate");
    g.ethylene_crowding_radius       = sg.get("ethylene_crowding_radius");
    g.ethylene_diffusion_radius      = sg.get("ethylene_diffusion_radius");
    g.ethylene_abscission_threshold  = sg.get("ethylene_abscission_threshold");
    g.ethylene_elongation_inhibition = sg.get("ethylene_elongation_inhibition");
    g.senescence_duration            = static_cast<uint32_t>(sg.get("senescence_duration"));

    return g;
}

} // namespace botany
```

- [ ] **Step 5: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[evolution]"`
Expected: All tests pass. Round-trip preserves all genome values. Template has 8 linkage groups.

- [ ] **Step 6: Commit**

```bash
git add src/evolution/genome_bridge.h src/evolution/genome_bridge.cpp tests/test_evolution.cpp
git commit -m "feat(evolution): genome bridge with StructuredGenome conversion and linkage groups"
```

---

### Task 4: Fitness — PlantStats + evaluate_plant

**Files:**
- Modify: `src/evolution/fitness.h`
- Modify: `src/evolution/fitness.cpp`
- Test: `tests/test_evolution.cpp`

- [ ] **Step 1: Write the failing test**

In `tests/test_evolution.cpp`, add:

```cpp
#include "evolution/fitness.h"

TEST_CASE("evaluate_plant returns populated stats", "[evolution]") {
    botany::Genome g = botany::default_genome();
    botany::WorldParams world = botany::default_world_params();

    auto stats = botany::evaluate_plant(g, world, 500);

    // Plant should survive some ticks and produce sugar
    REQUIRE(stats.survival_ticks > 0);
    REQUIRE(stats.node_count > 3);  // at least seed + 2 meristems
    REQUIRE(stats.total_sugar_produced > 0.0f);
    REQUIRE(stats.height > 0.0f);
}

TEST_CASE("evaluate_plant respects max_ticks", "[evolution]") {
    botany::Genome g = botany::default_genome();
    botany::WorldParams world = botany::default_world_params();

    auto stats = botany::evaluate_plant(g, world, 50);
    REQUIRE(stats.survival_ticks <= 50);
}
```

- [ ] **Step 2: Build and verify it fails**

Run: `/usr/local/bin/cmake --build build`
Expected: Compile error — `evaluate_plant` not declared, `PlantStats` not defined.

- [ ] **Step 3: Write fitness.h**

```cpp
#pragma once

#include <cstdint>
#include "engine/genome.h"
#include "engine/world_params.h"

namespace botany {

struct PlantStats {
    uint32_t survival_ticks = 0;
    uint32_t node_count = 0;
    uint32_t leaf_count = 0;
    float total_sugar_produced = 0.0f;
    float height = 0.0f;
    float crown_ratio = 0.0f;      // canopy_width / height
    uint32_t branch_depth = 0;      // max branching generations
    float leaf_height_spread = 0.0f; // std dev of leaf Y positions
};

struct FitnessWeights {
    float survival     = 1.0f;
    float biomass      = 1.0f;
    float sugar        = 1.0f;
    float leaves       = 1.0f;
    float height       = 1.0f;
    float crown_ratio  = 1.0f;
    float branch_depth = 1.0f;
    float leaf_spread  = 1.0f;
};

// Run a single plant simulation and collect stats.
PlantStats evaluate_plant(const Genome& genome, const WorldParams& world, uint32_t max_ticks);

// Compute weighted fitness score. gen_max contains the maximum value seen
// in the current generation for each stat (used for normalization).
float compute_fitness(const PlantStats& stats, const PlantStats& gen_max, const FitnessWeights& weights);

} // namespace botany
```

- [ ] **Step 4: Write fitness.cpp**

```cpp
#include "evolution/fitness.h"
#include "engine/engine.h"
#include "engine/node/node.h"
#include "engine/node/leaf_node.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace botany {

// Count branch depth: how many branch-from-branch generations.
// Trunk nodes (first child chain from seed) are depth 0.
// First branch off trunk is depth 1, branch off that is depth 2, etc.
static uint32_t compute_branch_depth(const Plant& plant) {
    uint32_t max_depth = 0;

    // Walk the tree, tracking branching depth.
    // A node is a "branch point" if its parent has multiple children.
    struct Entry { const Node* node; uint32_t depth; };
    std::vector<Entry> stack;
    stack.push_back({plant.seed(), 0});

    while (!stack.empty()) {
        auto [node, depth] = stack.back();
        stack.pop_back();
        max_depth = std::max(max_depth, depth);

        if (node->children.empty()) continue;

        // First child continues the same branch; others start new branches
        bool first = true;
        for (const Node* child : node->children) {
            if (child->type == NodeType::LEAF) continue; // leaves aren't branches
            if (child->is_meristem()) continue;           // meristems aren't branches
            if (first) {
                stack.push_back({child, depth});
                first = false;
            } else {
                stack.push_back({child, depth + 1});
            }
        }
    }
    return max_depth;
}

static bool has_active_meristems(const Plant& plant) {
    bool found = false;
    plant.for_each_node([&](const Node& n) {
        if (n.type == NodeType::SHOOT_APICAL || n.type == NodeType::ROOT_APICAL) {
            found = true;
        }
    });
    return found;
}

PlantStats evaluate_plant(const Genome& genome, const WorldParams& world, uint32_t max_ticks) {
    Engine engine;
    engine.world_params_mut() = world;
    PlantID pid = engine.create_plant(genome, glm::vec3(0.0f));

    uint32_t ticks = 0;
    for (; ticks < max_ticks; ticks++) {
        engine.tick();
        if (!has_active_meristems(engine.get_plant(pid))) break;
    }

    const Plant& plant = engine.get_plant(pid);
    PlantStats stats;
    stats.survival_ticks = ticks;
    stats.node_count = plant.node_count();
    stats.total_sugar_produced = plant.total_sugar_produced();

    // Collect per-node stats in one pass
    float min_y = 1e9f, max_y = -1e9f;
    float min_x = 1e9f, max_x = -1e9f;
    float min_z = 1e9f, max_z = -1e9f;
    std::vector<float> leaf_ys;

    plant.for_each_node([&](const Node& n) {
        min_y = std::min(min_y, n.position.y);
        max_y = std::max(max_y, n.position.y);
        min_x = std::min(min_x, n.position.x);
        max_x = std::max(max_x, n.position.x);
        min_z = std::min(min_z, n.position.z);
        max_z = std::max(max_z, n.position.z);

        if (n.type == NodeType::LEAF) {
            auto* leaf = n.as_leaf();
            if (leaf && leaf->senescence_ticks == 0) {
                stats.leaf_count++;
                leaf_ys.push_back(n.position.y);
            }
        }
    });

    // Height: max Y above seed (seed is at y=0)
    stats.height = std::max(0.0f, max_y);

    // Crown ratio: horizontal extent / height
    float width_x = max_x - min_x;
    float width_z = max_z - min_z;
    float canopy_width = std::max(width_x, width_z);
    stats.crown_ratio = stats.height > 0.01f ? canopy_width / stats.height : 0.0f;

    // Branch depth
    stats.branch_depth = compute_branch_depth(plant);

    // Leaf height spread: std dev of leaf Y positions
    if (leaf_ys.size() > 1) {
        float mean = 0.0f;
        for (float y : leaf_ys) mean += y;
        mean /= static_cast<float>(leaf_ys.size());

        float variance = 0.0f;
        for (float y : leaf_ys) {
            float d = y - mean;
            variance += d * d;
        }
        variance /= static_cast<float>(leaf_ys.size());
        stats.leaf_height_spread = std::sqrt(variance);
    }

    return stats;
}

float compute_fitness(const PlantStats& stats, const PlantStats& gen_max, const FitnessWeights& w) {
    auto norm = [](float val, float max_val) -> float {
        return max_val > 1e-9f ? val / max_val : 0.0f;
    };

    float score = 0.0f;
    score += w.survival     * norm(static_cast<float>(stats.survival_ticks), static_cast<float>(gen_max.survival_ticks));
    score += w.biomass      * norm(static_cast<float>(stats.node_count), static_cast<float>(gen_max.node_count));
    score += w.sugar        * norm(stats.total_sugar_produced, gen_max.total_sugar_produced);
    score += w.leaves       * norm(static_cast<float>(stats.leaf_count), static_cast<float>(gen_max.leaf_count));
    score += w.height       * norm(stats.height, gen_max.height);
    score += w.crown_ratio  * norm(stats.crown_ratio, gen_max.crown_ratio);
    score += w.branch_depth * norm(static_cast<float>(stats.branch_depth), static_cast<float>(gen_max.branch_depth));
    score += w.leaf_spread  * norm(stats.leaf_height_spread, gen_max.leaf_height_spread);
    return score;
}

} // namespace botany
```

- [ ] **Step 5: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[evolution]"`
Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/evolution/fitness.h src/evolution/fitness.cpp tests/test_evolution.cpp
git commit -m "feat(evolution): fitness evaluator with PlantStats and multi-objective scoring"
```

---

### Task 5: Fitness — compute_fitness unit test

**Files:**
- Test: `tests/test_evolution.cpp`

- [ ] **Step 1: Write tests for compute_fitness**

In `tests/test_evolution.cpp`, add:

```cpp
TEST_CASE("compute_fitness normalizes and weights correctly", "[evolution]") {
    botany::PlantStats stats;
    stats.survival_ticks = 100;
    stats.node_count = 50;
    stats.leaf_count = 10;
    stats.total_sugar_produced = 5.0f;
    stats.height = 2.0f;
    stats.crown_ratio = 0.5f;
    stats.branch_depth = 3;
    stats.leaf_height_spread = 1.0f;

    // gen_max = same values -> all normalized to 1.0
    botany::PlantStats gen_max = stats;

    botany::FitnessWeights w;  // all weights = 1.0
    float fitness = botany::compute_fitness(stats, gen_max, w);
    // 8 objectives, each normalized to 1.0, weight 1.0 -> sum = 8.0
    REQUIRE_THAT(fitness, WithinAbs(8.0, 1e-4));
}

TEST_CASE("compute_fitness handles zero gen_max gracefully", "[evolution]") {
    botany::PlantStats stats;
    stats.survival_ticks = 100;
    stats.height = 2.0f;

    botany::PlantStats gen_max;  // all zeros

    botany::FitnessWeights w;
    float fitness = botany::compute_fitness(stats, gen_max, w);
    // Zero max means all normalized to 0
    REQUIRE_THAT(fitness, WithinAbs(0.0, 1e-4));
}
```

- [ ] **Step 2: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[evolution]"`
Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
git add tests/test_evolution.cpp
git commit -m "test(evolution): add compute_fitness unit tests"
```

---

### Task 6: Evolution Runner

**Files:**
- Modify: `src/evolution/evolution_runner.h`
- Modify: `src/evolution/evolution_runner.cpp`
- Test: `tests/test_evolution.cpp`

- [ ] **Step 1: Write the failing test**

In `tests/test_evolution.cpp`, add:

```cpp
#include "evolution/evolution_runner.h"

TEST_CASE("EvolutionRunner advances generations", "[evolution]") {
    botany::EvolutionConfig config;
    config.population_size = 10;
    config.max_ticks = 200;
    config.num_threads = 2;

    botany::EvolutionRunner runner(config);
    REQUIRE(runner.generation() == 0);

    runner.run_generation();
    REQUIRE(runner.generation() == 1);
    REQUIRE(runner.best_fitness() > 0.0f);
    REQUIRE(runner.best_stats().survival_ticks > 0);

    runner.run_generation();
    REQUIRE(runner.generation() == 2);
}
```

- [ ] **Step 2: Build and verify it fails**

Run: `/usr/local/bin/cmake --build build`
Expected: Compile error — `EvolutionRunner`, `EvolutionConfig` not declared.

- [ ] **Step 3: Write evolution_runner.h**

```cpp
#pragma once

#include <cstdint>
#include <random>
#include <vector>
#include "evolution/fitness.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include <evolve/structured_genome.h>

namespace botany {

struct EvolutionConfig {
    uint32_t population_size = 100;
    uint32_t max_ticks = 17520;
    uint32_t num_threads = 4;
    uint32_t elitism_count = 2;
    uint32_t tournament_size = 5;
    float light_level_min = 0.5f;
    float light_level_max = 1.0f;
    float light_tilt_max = 0.52f;  // ~30 degrees max tilt
    FitnessWeights weights;
};

struct Individual {
    evolve::StructuredGenome genome;
    PlantStats stats;
    float fitness = 0.0f;
};

class EvolutionRunner {
public:
    explicit EvolutionRunner(const EvolutionConfig& config, uint32_t seed = 42);

    // Run one full generation: evaluate all, score, evolve.
    void run_generation();

    // Reset to generation 0 with fresh random population.
    void reset();

    // Accessors
    uint32_t generation() const { return generation_; }
    float best_fitness() const { return best_fitness_; }
    const PlantStats& best_stats() const { return best_stats_; }
    const evolve::StructuredGenome& best_genome() const { return best_genome_; }
    const std::vector<float>& fitness_history() const { return fitness_history_; }
    const EvolutionConfig& config() const { return config_; }
    EvolutionConfig& config_mut() { return config_; }

    // Convert best genome back to botany::Genome for rendering.
    Genome best_as_botany_genome() const;

private:
    void init_population();
    void evaluate_all();
    void score_all();
    void evolve_population();
    WorldParams randomize_world();
    const Individual& tournament_select();

    EvolutionConfig config_;
    std::mt19937 rng_;
    uint32_t generation_ = 0;

    std::vector<Individual> population_;
    evolve::StructuredGenome genome_template_;

    float best_fitness_ = 0.0f;
    PlantStats best_stats_;
    evolve::StructuredGenome best_genome_;
    std::vector<float> fitness_history_;
};

} // namespace botany
```

- [ ] **Step 4: Write evolution_runner.cpp**

```cpp
#include "evolution/evolution_runner.h"
#include "evolution/genome_bridge.h"
#include <algorithm>
#include <cmath>
#include <thread>
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>

namespace botany {

EvolutionRunner::EvolutionRunner(const EvolutionConfig& config, uint32_t seed)
    : config_(config), rng_(seed)
{
    genome_template_ = build_genome_template(default_genome());
    init_population();
}

void EvolutionRunner::init_population() {
    population_.clear();
    population_.resize(config_.population_size);

    // Start from the default genome with small random perturbations
    for (auto& ind : population_) {
        ind.genome = genome_template_;
        evolve::mutate(ind.genome, rng_);
    }
    // Keep one individual as the unperturbed default
    population_[0].genome = genome_template_;
}

void EvolutionRunner::reset() {
    generation_ = 0;
    best_fitness_ = 0.0f;
    best_stats_ = {};
    fitness_history_.clear();
    init_population();
}

WorldParams EvolutionRunner::randomize_world() {
    WorldParams world = default_world_params();

    // Randomize light level
    std::uniform_real_distribution<float> level_dist(config_.light_level_min, config_.light_level_max);
    world.light_level = level_dist(rng_);

    // Randomize light direction: tilt the up vector by a random angle
    std::uniform_real_distribution<float> angle_dist(0.0f, config_.light_tilt_max);
    std::uniform_real_distribution<float> azimuth_dist(0.0f, 6.2831853f);
    float tilt = angle_dist(rng_);
    float azimuth = azimuth_dist(rng_);
    float st = std::sin(tilt);
    world.light_direction = glm::normalize(glm::vec3(st * std::cos(azimuth), std::cos(tilt), st * std::sin(azimuth)));

    return world;
}

void EvolutionRunner::evaluate_all() {
    WorldParams world = randomize_world();
    uint32_t pop_size = static_cast<uint32_t>(population_.size());
    uint32_t num_threads = std::min(config_.num_threads, pop_size);

    if (num_threads <= 1) {
        for (auto& ind : population_) {
            Genome g = from_structured(ind.genome);
            ind.stats = evaluate_plant(g, world, config_.max_ticks);
        }
        return;
    }

    std::vector<std::thread> threads;
    uint32_t chunk = pop_size / num_threads;
    uint32_t remainder = pop_size % num_threads;
    uint32_t start = 0;

    for (uint32_t t = 0; t < num_threads; t++) {
        uint32_t end = start + chunk + (t < remainder ? 1 : 0);
        threads.emplace_back([this, &world, start, end]() {
            for (uint32_t i = start; i < end; i++) {
                Genome g = from_structured(population_[i].genome);
                population_[i].stats = evaluate_plant(g, world, config_.max_ticks);
            }
        });
        start = end;
    }

    for (auto& t : threads) t.join();
}

void EvolutionRunner::score_all() {
    // Find per-generation max for each stat
    PlantStats gen_max;
    for (const auto& ind : population_) {
        gen_max.survival_ticks = std::max(gen_max.survival_ticks, ind.stats.survival_ticks);
        gen_max.node_count = std::max(gen_max.node_count, ind.stats.node_count);
        gen_max.leaf_count = std::max(gen_max.leaf_count, ind.stats.leaf_count);
        gen_max.total_sugar_produced = std::max(gen_max.total_sugar_produced, ind.stats.total_sugar_produced);
        gen_max.height = std::max(gen_max.height, ind.stats.height);
        gen_max.crown_ratio = std::max(gen_max.crown_ratio, ind.stats.crown_ratio);
        gen_max.branch_depth = std::max(gen_max.branch_depth, ind.stats.branch_depth);
        gen_max.leaf_height_spread = std::max(gen_max.leaf_height_spread, ind.stats.leaf_height_spread);
    }

    // Score all individuals
    for (auto& ind : population_) {
        ind.fitness = compute_fitness(ind.stats, gen_max, config_.weights);
    }

    // Sort by fitness descending
    std::sort(population_.begin(), population_.end(),
              [](const Individual& a, const Individual& b) { return a.fitness > b.fitness; });

    // Record best
    best_fitness_ = population_[0].fitness;
    best_stats_ = population_[0].stats;
    best_genome_ = population_[0].genome;
    fitness_history_.push_back(best_fitness_);
}

const Individual& EvolutionRunner::tournament_select() {
    std::uniform_int_distribution<uint32_t> dist(0, static_cast<uint32_t>(population_.size()) - 1);
    uint32_t best_idx = dist(rng_);

    for (uint32_t i = 1; i < config_.tournament_size; i++) {
        uint32_t idx = dist(rng_);
        if (population_[idx].fitness > population_[best_idx].fitness) {
            best_idx = idx;
        }
    }
    return population_[best_idx];
}

void EvolutionRunner::evolve_population() {
    std::vector<Individual> next_gen;
    next_gen.reserve(config_.population_size);

    // Elitism: keep top N
    uint32_t elite = std::min(config_.elitism_count, static_cast<uint32_t>(population_.size()));
    for (uint32_t i = 0; i < elite; i++) {
        next_gen.push_back(population_[i]);  // already sorted by fitness
    }

    // Fill rest with crossover + mutation
    while (next_gen.size() < config_.population_size) {
        const Individual& parent_a = tournament_select();
        const Individual& parent_b = tournament_select();

        Individual child;
        child.genome = evolve::crossover(parent_a.genome, parent_b.genome, rng_);
        evolve::mutate(child.genome, rng_);
        next_gen.push_back(std::move(child));
    }

    population_ = std::move(next_gen);
}

void EvolutionRunner::run_generation() {
    evaluate_all();
    score_all();
    evolve_population();
    generation_++;
}

Genome EvolutionRunner::best_as_botany_genome() const {
    return from_structured(best_genome_);
}

} // namespace botany
```

- [ ] **Step 5: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[evolution]"`
Expected: All tests pass. Runner advances generations and reports positive fitness.

- [ ] **Step 6: Commit**

```bash
git add src/evolution/evolution_runner.h src/evolution/evolution_runner.cpp tests/test_evolution.cpp
git commit -m "feat(evolution): EvolutionRunner with threaded evaluation and tournament selection"
```

---

### Task 7: App — botany_evolve

**Files:**
- Modify: `src/app_evolve.cpp`

No unit test for this task — it's GUI code. Verification is visual: launch the app, start evolution, confirm stats update and best plant renders.

- [ ] **Step 1: Write app_evolve.cpp**

Follow the same GLFW + ImGui pattern as `src/app_realtime.cpp`. The app has:
- An `EvolutionRunner` for the GA loop
- A separate `Engine` for re-simulating + rendering the best plant
- ImGui panel with config sliders, controls, and stats
- The existing `Renderer` for the viewport

```cpp
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <glm/vec3.hpp>
#include <cfloat>
#include <iostream>
#include <thread>
#include <atomic>
#include "engine/engine.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "evolution/evolution_runner.h"
#include "evolution/genome_bridge.h"
#include "renderer/renderer.h"

using namespace botany;

static Renderer* g_renderer = nullptr;

static void scroll_callback(GLFWwindow*, double, double yoffset) {
    if (g_renderer) g_renderer->camera().zoom(static_cast<float>(yoffset));
}

int main() {
    Renderer renderer;
    if (!renderer.init(1280, 800, "shaders")) {
        std::cerr << "Failed to initialize renderer" << std::endl;
        return 1;
    }
    g_renderer = &renderer;
    GLFWwindow* window = renderer.window();
    glfwSetScrollCallback(window, scroll_callback);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 410");
    ImGui::StyleColorsDark();

    EvolutionConfig evo_config;
    evo_config.population_size = 50;
    evo_config.max_ticks = 5000;
    evo_config.num_threads = std::max(1u, std::thread::hardware_concurrency());
    EvolutionRunner runner(evo_config);

    // Display engine: re-simulates the best plant for rendering
    Engine display_engine;
    PlantID display_plant = display_engine.create_plant(default_genome(), glm::vec3(0.0f));
    bool display_needs_update = false;

    bool running = false;
    std::atomic<bool> gen_in_progress{false};
    std::thread evo_thread;

    // Config state (editable while paused)
    int pop_size = static_cast<int>(evo_config.population_size);
    int max_ticks = static_cast<int>(evo_config.max_ticks);
    int num_threads = static_cast<int>(evo_config.num_threads);

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, true);

        // Camera controls
        const float rotate_speed = 6.0f;
        if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) renderer.camera().rotate(-rotate_speed, 0.0f);
        if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) renderer.camera().rotate(rotate_speed, 0.0f);
        if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) renderer.camera().rotate(0.0f, rotate_speed);
        if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) renderer.camera().rotate(0.0f, -rotate_speed);
        const float pan_speed = 0.05f;
        if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS)
            renderer.camera().set_target(renderer.camera().target() + glm::vec3(0.0f, pan_speed, 0.0f));
        if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS)
            renderer.camera().set_target(renderer.camera().target() - glm::vec3(0.0f, pan_speed, 0.0f));

        // Check if background generation finished
        if (gen_in_progress && !evo_thread.joinable()) {
            gen_in_progress = false;
        }
        if (gen_in_progress && evo_thread.joinable()) {
            // non-blocking: thread still running, skip
        }
        if (!gen_in_progress && evo_thread.joinable()) {
            evo_thread.join();
            display_needs_update = true;
        }

        // Auto-launch next generation if running
        if (running && !gen_in_progress) {
            gen_in_progress = true;
            evo_thread = std::thread([&runner, &gen_in_progress]() {
                runner.run_generation();
                gen_in_progress = false;
            });
        }

        // Re-simulate best plant for display
        if (display_needs_update && runner.generation() > 0) {
            display_engine.reset();
            Genome best_g = runner.best_as_botany_genome();
            display_plant = display_engine.create_plant(best_g, glm::vec3(0.0f));
            uint32_t replay_ticks = runner.best_stats().survival_ticks;
            for (uint32_t i = 0; i < replay_ticks; i++) {
                display_engine.tick();
            }
            display_needs_update = false;
        }

        // ImGui
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_Once);
        ImGui::Begin("Evolution", nullptr, ImGuiWindowFlags_AlwaysAutoResize);

        // --- Config section ---
        if (ImGui::CollapsingHeader("Configuration", ImGuiTreeNodeFlags_DefaultOpen)) {
            bool paused = !running;
            if (paused) {
                ImGui::SliderInt("Population", &pop_size, 10, 500);
                ImGui::SliderInt("Max Ticks", &max_ticks, 100, 20000);
                ImGui::SliderInt("Threads", &num_threads, 1, 16);
            } else {
                ImGui::Text("Population: %d  Max Ticks: %d  Threads: %d", pop_size, max_ticks, num_threads);
            }

            ImGui::SeparatorText("Fitness Weights");
            FitnessWeights& w = runner.config_mut().weights;
            ImGui::SliderFloat("Survival", &w.survival, 0.0f, 5.0f);
            ImGui::SliderFloat("Biomass", &w.biomass, 0.0f, 5.0f);
            ImGui::SliderFloat("Sugar", &w.sugar, 0.0f, 5.0f);
            ImGui::SliderFloat("Leaves", &w.leaves, 0.0f, 5.0f);
            ImGui::SliderFloat("Height", &w.height, 0.0f, 5.0f);
            ImGui::SliderFloat("Crown Ratio", &w.crown_ratio, 0.0f, 5.0f);
            ImGui::SliderFloat("Branch Depth", &w.branch_depth, 0.0f, 5.0f);
            ImGui::SliderFloat("Leaf Spread", &w.leaf_spread, 0.0f, 5.0f);
        }

        // --- Controls ---
        ImGui::SeparatorText("Controls");
        if (!running) {
            if (ImGui::Button("Start")) {
                runner.config_mut().population_size = static_cast<uint32_t>(pop_size);
                runner.config_mut().max_ticks = static_cast<uint32_t>(max_ticks);
                runner.config_mut().num_threads = static_cast<uint32_t>(num_threads);
                if (runner.generation() == 0) runner.reset();
                running = true;
            }
        } else {
            if (ImGui::Button("Pause")) {
                running = false;
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Reset")) {
            running = false;
            if (evo_thread.joinable()) evo_thread.join();
            gen_in_progress = false;
            runner.config_mut().population_size = static_cast<uint32_t>(pop_size);
            runner.config_mut().max_ticks = static_cast<uint32_t>(max_ticks);
            runner.config_mut().num_threads = static_cast<uint32_t>(num_threads);
            runner.reset();
            display_engine.reset();
            display_plant = display_engine.create_plant(default_genome(), glm::vec3(0.0f));
        }

        // --- Stats ---
        if (runner.generation() > 0) {
            ImGui::SeparatorText("Stats");
            ImGui::Text("Generation: %u", runner.generation());
            ImGui::Text("Best Fitness: %.3f", runner.best_fitness());
            if (gen_in_progress) ImGui::Text("(evaluating...)");

            const PlantStats& s = runner.best_stats();
            if (ImGui::BeginTable("stats", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                ImGui::TableSetupColumn("Metric", ImGuiTableColumnFlags_WidthFixed, 120);
                ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthFixed, 100);
                ImGui::TableHeadersRow();

                auto row = [](const char* label, const char* fmt, auto val) {
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("%s", label);
                    ImGui::TableSetColumnIndex(1); ImGui::Text(fmt, val);
                };
                row("Survival", "%u ticks", s.survival_ticks);
                row("Nodes", "%u", s.node_count);
                row("Leaves", "%u", s.leaf_count);
                row("Sugar", "%.1f g", s.total_sugar_produced);
                row("Height", "%.2f dm", s.height);
                row("Crown Ratio", "%.2f", s.crown_ratio);
                row("Branch Depth", "%u", s.branch_depth);
                row("Leaf Spread", "%.2f dm", s.leaf_height_spread);

                ImGui::EndTable();
            }

            // Fitness history plot
            const auto& hist = runner.fitness_history();
            if (hist.size() > 1) {
                ImGui::PlotLines("Fitness", hist.data(), static_cast<int>(hist.size()),
                                 0, nullptr, FLT_MAX, FLT_MAX, ImVec2(0, 80));
            }
        }

        ImGui::End();
        ImGui::Render();

        // Render
        renderer.begin_frame();
        renderer.draw_ground();
        renderer.draw_plant(display_engine.get_plant(display_plant));
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        renderer.end_frame();
    }

    // Cleanup
    running = false;
    if (evo_thread.joinable()) evo_thread.join();

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    renderer.shutdown();
    return 0;
}
```

- [ ] **Step 2: Build and verify**

Run: `/usr/local/bin/cmake --build build`
Expected: Build succeeds. The `botany_evolve` binary is produced.

- [ ] **Step 3: Smoke test**

Run: `./build/botany_evolve`
Expected: Window opens with config panel. Click "Start". Generations advance, stats update, best plant renders in viewport. Click "Pause" to stop.

- [ ] **Step 4: Commit**

```bash
git add src/app_evolve.cpp
git commit -m "feat(evolution): botany_evolve app with ImGui config, live stats, and best-plant rendering"
```

---

### Task 8: Export Best Genome

**Files:**
- Modify: `src/app_evolve.cpp`

- [ ] **Step 1: Add export button to ImGui panel**

In `app_evolve.cpp`, after the stats table inside the `if (runner.generation() > 0)` block, add:

```cpp
if (ImGui::Button("Export Best")) {
    Genome best_g = runner.best_as_botany_genome();
    std::ofstream out("best_genome.txt");
    if (out) {
        out << "auxin_production_rate=" << best_g.auxin_production_rate << "\n";
        out << "auxin_transport_rate=" << best_g.auxin_transport_rate << "\n";
        out << "auxin_directional_bias=" << best_g.auxin_directional_bias << "\n";
        out << "auxin_decay_rate=" << best_g.auxin_decay_rate << "\n";
        out << "auxin_threshold=" << best_g.auxin_threshold << "\n";
        out << "cytokinin_production_rate=" << best_g.cytokinin_production_rate << "\n";
        out << "cytokinin_transport_rate=" << best_g.cytokinin_transport_rate << "\n";
        out << "cytokinin_directional_bias=" << best_g.cytokinin_directional_bias << "\n";
        out << "cytokinin_decay_rate=" << best_g.cytokinin_decay_rate << "\n";
        out << "cytokinin_threshold=" << best_g.cytokinin_threshold << "\n";
        out << "growth_rate=" << best_g.growth_rate << "\n";
        out << "max_internode_length=" << best_g.max_internode_length << "\n";
        out << "min_internode_length=" << best_g.min_internode_length << "\n";
        out << "branch_angle=" << best_g.branch_angle << "\n";
        out << "thickening_rate=" << best_g.thickening_rate << "\n";
        out << "internode_elongation_rate=" << best_g.internode_elongation_rate << "\n";
        out << "internode_maturation_ticks=" << best_g.internode_maturation_ticks << "\n";
        out << "root_growth_rate=" << best_g.root_growth_rate << "\n";
        out << "root_max_internode_length=" << best_g.root_max_internode_length << "\n";
        out << "root_min_internode_length=" << best_g.root_min_internode_length << "\n";
        out << "root_branch_angle=" << best_g.root_branch_angle << "\n";
        out << "root_internode_elongation_rate=" << best_g.root_internode_elongation_rate << "\n";
        out << "root_internode_maturation_ticks=" << best_g.root_internode_maturation_ticks << "\n";
        out << "root_gravitropism_strength=" << best_g.root_gravitropism_strength << "\n";
        out << "root_gravitropism_depth=" << best_g.root_gravitropism_depth << "\n";
        out << "max_leaf_size=" << best_g.max_leaf_size << "\n";
        out << "leaf_growth_rate=" << best_g.leaf_growth_rate << "\n";
        out << "leaf_bud_size=" << best_g.leaf_bud_size << "\n";
        out << "initial_radius=" << best_g.initial_radius << "\n";
        out << "root_initial_radius=" << best_g.root_initial_radius << "\n";
        out << "tip_offset=" << best_g.tip_offset << "\n";
        out << "growth_noise=" << best_g.growth_noise << "\n";
        out << "leaf_phototropism_rate=" << best_g.leaf_phototropism_rate << "\n";
        out << "sugar_production_rate=" << best_g.sugar_production_rate << "\n";
        out << "sugar_transport_conductance=" << best_g.sugar_transport_conductance << "\n";
        out << "sugar_maintenance_leaf=" << best_g.sugar_maintenance_leaf << "\n";
        out << "sugar_maintenance_stem=" << best_g.sugar_maintenance_stem << "\n";
        out << "sugar_maintenance_root=" << best_g.sugar_maintenance_root << "\n";
        out << "sugar_maintenance_meristem=" << best_g.sugar_maintenance_meristem << "\n";
        out << "seed_sugar=" << best_g.seed_sugar << "\n";
        out << "sugar_storage_density_wood=" << best_g.sugar_storage_density_wood << "\n";
        out << "sugar_storage_density_leaf=" << best_g.sugar_storage_density_leaf << "\n";
        out << "sugar_cap_minimum=" << best_g.sugar_cap_minimum << "\n";
        out << "sugar_cap_meristem=" << best_g.sugar_cap_meristem << "\n";
        out << "sugar_save_shoot=" << best_g.sugar_save_shoot << "\n";
        out << "sugar_save_root=" << best_g.sugar_save_root << "\n";
        out << "sugar_save_stem=" << best_g.sugar_save_stem << "\n";
        out << "sugar_activation_shoot=" << best_g.sugar_activation_shoot << "\n";
        out << "sugar_activation_root=" << best_g.sugar_activation_root << "\n";
        out << "ga_production_rate=" << best_g.ga_production_rate << "\n";
        out << "ga_leaf_age_max=" << best_g.ga_leaf_age_max << "\n";
        out << "ga_elongation_sensitivity=" << best_g.ga_elongation_sensitivity << "\n";
        out << "ga_length_sensitivity=" << best_g.ga_length_sensitivity << "\n";
        out << "ga_transport_rate=" << best_g.ga_transport_rate << "\n";
        out << "ga_directional_bias=" << best_g.ga_directional_bias << "\n";
        out << "ga_decay_rate=" << best_g.ga_decay_rate << "\n";
        out << "ethylene_starvation_rate=" << best_g.ethylene_starvation_rate << "\n";
        out << "ethylene_shade_rate=" << best_g.ethylene_shade_rate << "\n";
        out << "ethylene_shade_threshold=" << best_g.ethylene_shade_threshold << "\n";
        out << "ethylene_age_rate=" << best_g.ethylene_age_rate << "\n";
        out << "ethylene_age_onset=" << best_g.ethylene_age_onset << "\n";
        out << "ethylene_crowding_rate=" << best_g.ethylene_crowding_rate << "\n";
        out << "ethylene_crowding_radius=" << best_g.ethylene_crowding_radius << "\n";
        out << "ethylene_diffusion_radius=" << best_g.ethylene_diffusion_radius << "\n";
        out << "ethylene_abscission_threshold=" << best_g.ethylene_abscission_threshold << "\n";
        out << "ethylene_elongation_inhibition=" << best_g.ethylene_elongation_inhibition << "\n";
        out << "senescence_duration=" << best_g.senescence_duration << "\n";
        std::cout << "Exported best genome to best_genome.txt" << std::endl;
    }
}
```

Add `#include <fstream>` at the top of the file.

- [ ] **Step 2: Build and verify**

Run: `/usr/local/bin/cmake --build build`
Expected: Build succeeds.

- [ ] **Step 3: Commit**

```bash
git add src/app_evolve.cpp
git commit -m "feat(evolution): add Export Best button to save genome as text file"
```
