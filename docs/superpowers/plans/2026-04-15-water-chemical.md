# Water Chemical Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add water as a capacity-based resource chemical — roots absorb it, leaves lose it via transpiration and photosynthesis, transport uses existing unified system.

**Architecture:** Water follows the exact same pattern as sugar: capacity-based storage, concentration-driven diffusion via `transport_chemical()`, production in `produce()` methods. New `ChemicalID::Water` entry, new genome params, new `WorldParams::soil_moisture`, new `water_cap()` function mirroring `sugar_cap()`.

**Tech Stack:** C++17, Catch2, existing chemical transport framework.

---

### Task 1: Add Water to ChemicalID and Chemical Registry

**Files:**
- Modify: `src/engine/chemical/chemical.h:8-15`
- Modify: `src/engine/chemical/chemical_registry.h:12-40`
- Modify: `src/engine/node/node.cpp:28-34`

- [ ] **Step 1: Add `Water` to `ChemicalID` enum**

In `src/engine/chemical/chemical.h`, add `Water` after `Stress`:

```cpp
enum class ChemicalID : uint8_t {
    Auxin,
    Cytokinin,
    Gibberellin,
    Sugar,
    Ethylene,
    Stress,
    Water,
};
```

- [ ] **Step 2: Add Water to `all_chemical_ids` array**

In `src/engine/chemical/chemical_registry.h`, update the array size from 6 to 7 and add `Water`:

```cpp
inline constexpr std::array<ChemicalID, 7> all_chemical_ids = {
    ChemicalID::Auxin,
    ChemicalID::Cytokinin,
    ChemicalID::Gibberellin,
    ChemicalID::Sugar,
    ChemicalID::Ethylene,
    ChemicalID::Stress,
    ChemicalID::Water,
};
```

- [ ] **Step 3: Initialize Water to zero in Node constructor**

In `src/engine/node/node.cpp`, add to the chemical initialization block (after line 34):

```cpp
chemicals[ChemicalID::Water] = 0.0f;
```

- [ ] **Step 4: Build and verify**

Run: `/usr/local/bin/cmake --build build`
Expected: Clean build, no errors.

- [ ] **Step 5: Commit**

```bash
git add src/engine/chemical/chemical.h src/engine/chemical/chemical_registry.h src/engine/node/node.cpp
git commit -m "feat: add Water to ChemicalID enum and registry"
```

---

### Task 2: Add Water Genome Parameters

**Files:**
- Modify: `src/engine/genome.h:83-92` (add after sugar economy section)

- [ ] **Step 1: Add water fields to `Genome` struct**

In `src/engine/genome.h`, add after the sugar cap fields (after line 92, before `leaf_phototropism_rate`):

```cpp
    // Water economy (ml)
    float water_absorption_rate;          // ml / (dm² root surface · hr) per unit soil_moisture
    float transpiration_rate;             // ml / (dm² leaf area · hr) at full light
    float photosynthesis_water_ratio;     // ml water consumed per g sugar produced
    float water_storage_density_stem;     // ml / dm³ of stem/root tissue
    float water_storage_density_leaf;     // ml / dm² of leaf area
    float water_cap_meristem;             // ml — fixed cap for meristem nodes
    float water_diffusion_rate;           // fraction diffused per tick
    float water_bias;                     // upward equilibrium shift (positive = toward tips)
    float water_base_transport;           // throughput floor
    float water_transport_scale;          // radius scaling on throughput
```

- [ ] **Step 2: Add default values in `default_genome()`**

In `src/engine/genome.h`, add after `.sugar_cap_meristem = 2.0f,` (line 215):

```cpp
        // Water economy
        .water_absorption_rate = 0.05f,          // ml / (dm² · hr) — moderate absorption
        .transpiration_rate = 0.04f,             // ml / (dm² · hr) — slightly less than absorption
        .photosynthesis_water_ratio = 0.5f,      // 0.5 ml water per g sugar (small cost)
        .water_storage_density_stem = 800.0f,    // ml / dm³ — wood is ~80% water by volume
        .water_storage_density_leaf = 3.0f,      // ml / dm² — leaves hold water in vacuoles
        .water_cap_meristem = 1.0f,              // ml — small active reserve
        .water_diffusion_rate = 0.9f,            // faster than sugar (0.8) — water moves easily
        .water_bias = 0.05f,                     // slight upward bias (transpiration pull)
        .water_base_transport = 0.2f,            // higher floor than sugar — xylem is open pipes
        .water_transport_scale = 4.0f,           // radius matters but less than sugar
```

- [ ] **Step 3: Build and verify**

Run: `/usr/local/bin/cmake --build build`
Expected: Clean build, no errors.

- [ ] **Step 4: Commit**

```bash
git add src/engine/genome.h
git commit -m "feat: add water economy genome parameters"
```

---

### Task 3: Add `soil_moisture` to WorldParams

**Files:**
- Modify: `src/engine/world_params.h`

- [ ] **Step 1: Add `soil_moisture` field**

In `src/engine/world_params.h`, add after `light_level`:

```cpp
    float soil_moisture = 1.0f;           // 0-1 — soil water availability
```

- [ ] **Step 2: Build and verify**

Run: `/usr/local/bin/cmake --build build`
Expected: Clean build, no errors.

- [ ] **Step 3: Commit**

```bash
git add src/engine/world_params.h
git commit -m "feat: add soil_moisture to WorldParams"
```

---

### Task 4: Add `water_cap()` Function

**Files:**
- Modify: `src/engine/sugar.h:10` (add declaration)
- Modify: `src/engine/sugar.cpp` (add implementation)
- Create: `tests/test_water.cpp`

- [ ] **Step 1: Write failing tests for `water_cap()`**

Create `tests/test_water.cpp`:

```cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <glm/geometric.hpp>
#include "engine/plant.h"
#include "engine/node/tissues/leaf.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/sugar.h"

using namespace botany;
using Catch::Matchers::WithinAbs;

// === Water cap tests ===

TEST_CASE("water_cap scales with stem volume", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.1f);
    plant.seed_mut()->add_child(stem);

    float volume = 3.14159f * 0.1f * 0.1f * 1.0f;
    float expected = volume * g.water_storage_density_stem;
    REQUIRE_THAT(water_cap(*stem, g), WithinAbs(expected, 1e-4));
}

TEST_CASE("water_cap scales with leaf area", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    plant.seed_mut()->add_child(leaf);

    float area = 0.5f * 0.5f;
    float expected = area * g.water_storage_density_leaf;
    REQUIRE_THAT(water_cap(*leaf, g), WithinAbs(expected, 1e-6));
}

TEST_CASE("water_cap for root scales with volume", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* root = plant.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.8f, 0.0f), 0.025f);
    plant.seed_mut()->add_child(root);

    float length = 0.8f;
    float volume = 3.14159f * 0.025f * 0.025f * length;
    float expected = volume * g.water_storage_density_stem;
    REQUIRE_THAT(water_cap(*root, g), WithinAbs(expected, 1e-5));
}

TEST_CASE("water_cap for meristem is fixed", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* apical = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::APICAL) apical = &n;
    });
    REQUIRE(apical != nullptr);
    REQUIRE_THAT(water_cap(*apical, g), WithinAbs(g.water_cap_meristem, 1e-6));
}

TEST_CASE("water_cap returns minimum for tiny nodes", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* tiny = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.001f, 0.0f), 0.001f);
    plant.seed_mut()->add_child(tiny);

    REQUIRE_THAT(water_cap(*tiny, g), WithinAbs(g.sugar_cap_minimum, 1e-6));
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[water]" -v`
Expected: FAIL — `water_cap` not declared.

- [ ] **Step 3: Add `water_cap` declaration to `sugar.h`**

In `src/engine/sugar.h`, add after the `sugar_cap` declaration:

```cpp
float water_cap(const Node& node, const Genome& g);
```

- [ ] **Step 4: Implement `water_cap` in `sugar.cpp`**

In `src/engine/sugar.cpp`, add after the `sugar_cap` function:

```cpp
float water_cap(const Node& node, const Genome& g) {
    float length = std::max(glm::length(node.offset), 0.01f);

    switch (node.type) {
        case NodeType::LEAF: {
            float ls = node.as_leaf()->leaf_size;
            float area = ls * ls;
            return std::max(g.sugar_cap_minimum, area * g.water_storage_density_leaf);
        }
        case NodeType::STEM:
        case NodeType::ROOT: {
            float volume = 3.14159f * node.radius * node.radius * length;
            float cap = volume * g.water_storage_density_stem;
            return std::max(g.sugar_cap_minimum, cap);
        }
        case NodeType::APICAL:
        case NodeType::ROOT_APICAL:
            return g.water_cap_meristem;
    }
    return g.sugar_cap_minimum;
}
```

- [ ] **Step 5: Add test file to CMakeLists.txt**

Find the test target in `CMakeLists.txt` and add `tests/test_water.cpp` to its sources.

- [ ] **Step 6: Run tests to verify they pass**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[water]" -v`
Expected: All 5 tests PASS.

- [ ] **Step 7: Commit**

```bash
git add src/engine/sugar.h src/engine/sugar.cpp tests/test_water.cpp CMakeLists.txt
git commit -m "feat: add water_cap() with capacity tests"
```

---

### Task 5: Add Water to Transport System

**Files:**
- Modify: `src/engine/chemical/chemical_registry.h:23-40`
- Modify: `src/engine/node/node.cpp:298` (has_cap check)

- [ ] **Step 1: Write failing test for water transport conservation**

Append to `tests/test_water.cpp`:

```cpp
TEST_CASE("Water diffuses from high to low concentration", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* child = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.1f, 0.0f), 0.02f);
    plant.seed_mut()->add_child(child);

    plant.seed_mut()->chemical(ChemicalID::Water) = 50.0f;
    child->chemical(ChemicalID::Water) = 0.0f;

    plant.seed_mut()->transport_with_children(g);

    REQUIRE(child->chemical(ChemicalID::Water) > 0.0f);
    REQUIRE(plant.seed()->chemical(ChemicalID::Water) < 50.0f);
}

TEST_CASE("Water transport conserves total water", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    for (int i = 0; i < 3; i++) {
        Node* child = plant.create_node(NodeType::STEM,
            glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
        plant.seed_mut()->add_child(child);
    }

    plant.seed_mut()->chemical(ChemicalID::Water) = 20.0f;
    plant.seed_mut()->children[0]->chemical(ChemicalID::Water) = 5.0f;

    float total_before = 0.0f;
    plant.for_each_node([&](const Node& n) { total_before += n.chemical(ChemicalID::Water); });

    for (int i = 0; i < 10; i++) {
        plant.for_each_node_mut([&](Node& n) {
            n.transport_with_children(g);
        });
    }

    float total_after = 0.0f;
    plant.for_each_node([&](const Node& n) { total_after += n.chemical(ChemicalID::Water); });

    REQUIRE_THAT(total_after, WithinAbs(total_before, 1e-4));
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[water]" -v`
Expected: FAIL — water is not in `diffusion_params`, so no transport occurs.

- [ ] **Step 3: Add Water to `diffusion_params`**

In `src/engine/chemical/chemical_registry.h`, update the array size from 5 to 6 and add the Water entry:

```cpp
inline std::array<ChemicalDiffusionParams, 6> diffusion_params(const Genome& g) {
    return {{
        {ChemicalID::Auxin,       g.auxin_diffusion_rate,     g.auxin_decay_rate,
         g.auxin_bias,     g.hormone_base_transport, g.hormone_transport_scale},

        {ChemicalID::Cytokinin,   g.cytokinin_diffusion_rate, g.cytokinin_decay_rate,
         g.cytokinin_bias, g.hormone_base_transport, g.hormone_transport_scale},

        {ChemicalID::Gibberellin, g.ga_diffusion_rate,        g.ga_decay_rate,
         0.0f,             g.hormone_base_transport, g.hormone_transport_scale},

        {ChemicalID::Sugar,       g.sugar_diffusion_rate,     0.0f,
         0.0f,             g.sugar_base_transport,   g.sugar_transport_scale},

        {ChemicalID::Stress,      g.stress_hormone_diffusion_rate, g.stress_hormone_decay_rate,
         0.0f,             g.hormone_base_transport, g.hormone_transport_scale},

        {ChemicalID::Water,       g.water_diffusion_rate,     0.0f,
         g.water_bias,     g.water_base_transport,   g.water_transport_scale},
    }};
}
```

- [ ] **Step 4: Update `has_cap` check to include Water**

In `src/engine/node/node.cpp`, in `transport_with_children()`, change the `has_cap` line (line 298):

From:
```cpp
bool has_cap = (dp.id == ChemicalID::Sugar);
```
To:
```cpp
bool has_cap = (dp.id == ChemicalID::Sugar || dp.id == ChemicalID::Water);
```

- [ ] **Step 5: Update `parent_cap` computation to use water_cap**

In the same function, after line 295 (`float parent_cap_sugar = sugar_cap(*this, g);`), add:

```cpp
float parent_cap_water = water_cap(*this, g);
```

Then in the child loop (around line 321), update the parent_cap selection:

From:
```cpp
float parent_cap = has_cap ? parent_cap_sugar : 0.0f;
```
To:
```cpp
float parent_cap = has_cap
    ? (dp.id == ChemicalID::Water ? parent_cap_water : parent_cap_sugar)
    : 0.0f;
```

And update the child_cap line similarly:

From:
```cpp
float child_cap = has_cap ? sugar_cap(*child, g) : 0.0f;
```
To:
```cpp
float child_cap = has_cap
    ? (dp.id == ChemicalID::Water ? water_cap(*child, g) : sugar_cap(*child, g))
    : 0.0f;
```

Also update the `parent_headroom` line (around line 342):

From:
```cpp
float parent_headroom = has_cap ? std::max(0.0f, parent_cap_sugar - parent_val) : 1e30f;
```
To:
```cpp
float current_parent_cap = (dp.id == ChemicalID::Water) ? parent_cap_water : parent_cap_sugar;
float parent_headroom = has_cap ? std::max(0.0f, current_parent_cap - parent_val) : 1e30f;
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[water]" -v`
Expected: All 7 tests PASS.

- [ ] **Step 7: Run full test suite to check for regressions**

Run: `./build/botany_tests -v`
Expected: All tests PASS (existing sugar tests still pass).

- [ ] **Step 8: Commit**

```bash
git add src/engine/chemical/chemical_registry.h src/engine/node/node.cpp tests/test_water.cpp
git commit -m "feat: add water to transport system with capacity model"
```

---

### Task 6: Root Water Absorption

**Files:**
- Modify: `src/engine/node/root_node.h`
- Modify: `src/engine/node/root_node.cpp`
- Modify: `src/engine/node/tissues/root_apical.cpp`
- Modify: `tests/test_water.cpp`

- [ ] **Step 1: Write failing tests for root absorption**

Append to `tests/test_water.cpp`:

```cpp
#include "engine/node/root_node.h"
#include "engine/node/tissues/root_apical.h"
#include "engine/light.h"

TEST_CASE("Root nodes absorb water proportional to surface area", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* root = plant.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.5f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(root);
    root->chemical(ChemicalID::Water) = 0.0f;
    root->chemical(ChemicalID::Sugar) = 10.0f;  // prevent starvation

    WorldParams wp = default_world_params();
    wp.soil_moisture = 1.0f;

    root->tick(plant, wp);

    REQUIRE(root->chemical(ChemicalID::Water) > 0.0f);
}

TEST_CASE("Root absorption scales with soil_moisture", "[water]") {
    Genome g = default_genome();

    // Plant with high soil moisture
    Plant plant_wet(g, glm::vec3(0.0f));
    Node* root_wet = plant_wet.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.5f, 0.0f), 0.05f);
    plant_wet.seed_mut()->add_child(root_wet);
    root_wet->chemical(ChemicalID::Water) = 0.0f;
    root_wet->chemical(ChemicalID::Sugar) = 10.0f;

    WorldParams wp_wet = default_world_params();
    wp_wet.soil_moisture = 1.0f;
    root_wet->tick(plant_wet, wp_wet);

    // Plant with low soil moisture
    Plant plant_dry(g, glm::vec3(0.0f));
    Node* root_dry = plant_dry.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.5f, 0.0f), 0.05f);
    plant_dry.seed_mut()->add_child(root_dry);
    root_dry->chemical(ChemicalID::Water) = 0.0f;
    root_dry->chemical(ChemicalID::Sugar) = 10.0f;

    WorldParams wp_dry = default_world_params();
    wp_dry.soil_moisture = 0.2f;
    root_dry->tick(plant_dry, wp_dry);

    REQUIRE(root_wet->chemical(ChemicalID::Water) > root_dry->chemical(ChemicalID::Water));
}

TEST_CASE("Zero soil_moisture means zero water absorption", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* root = plant.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.5f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(root);
    root->chemical(ChemicalID::Water) = 0.0f;
    root->chemical(ChemicalID::Sugar) = 10.0f;

    WorldParams wp = default_world_params();
    wp.soil_moisture = 0.0f;

    root->tick(plant, wp);

    REQUIRE(root->chemical(ChemicalID::Water) < 1e-6f);
}

TEST_CASE("Root apical tips also absorb water", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Find a root apical in the default plant
    Node* root_apical = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::ROOT_APICAL) root_apical = &n;
    });
    REQUIRE(root_apical != nullptr);

    root_apical->chemical(ChemicalID::Water) = 0.0f;
    root_apical->chemical(ChemicalID::Sugar) = 10.0f;

    WorldParams wp = default_world_params();
    wp.soil_moisture = 1.0f;

    root_apical->tick(plant, wp);

    REQUIRE(root_apical->chemical(ChemicalID::Water) > 0.0f);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[water]" -v`
Expected: FAIL — root nodes don't absorb water yet.

- [ ] **Step 3: Add water absorption to `RootNode::tissue_tick()`**

In `src/engine/node/root_node.cpp`, add `#include "engine/sugar.h"` at the top. Then add water absorption at the start of `tissue_tick()` (before `thicken`):

```cpp
void RootNode::tissue_tick(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    // Water absorption: proportional to root surface area and soil moisture
    float length = std::max(glm::length(offset), 0.01f);
    float surface_area = 2.0f * 3.14159f * radius * length;
    float absorbed = g.water_absorption_rate * surface_area * world.soil_moisture;
    float cap = water_cap(*this, g);
    chemical(ChemicalID::Water) = std::min(chemical(ChemicalID::Water) + absorbed, cap);

    thicken(g, world);
    elongate(g, world);
}
```

- [ ] **Step 4: Add water absorption to `RootApicalNode::tissue_tick()`**

In `src/engine/node/tissues/root_apical.cpp`, add `#include "engine/sugar.h"` at the top. Then add water absorption at the start of `tissue_tick()` (before the active check):

```cpp
void RootApicalNode::tissue_tick(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    // Water absorption: hemisphere approximation for tip surface area
    float surface_area = 2.0f * 3.14159f * radius * radius;
    float absorbed = g.water_absorption_rate * surface_area * world.soil_moisture;
    float cap = water_cap(*this, g);
    chemical(ChemicalID::Water) = std::min(chemical(ChemicalID::Water) + absorbed, cap);

    if (!active) {
        if (can_activate(g, world)) activate(g, world);
        return;
    }

    // Chain growth
    ticks_since_last_node++;
    grow_tip(g, world);

    // Time-based spawning (plastochron). Don't spawn if starving.
    if (parent && ticks_since_last_node >= g.root_plastochron && starvation_ticks == 0) {
        spawn_internode(plant, g);
    }
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[water]" -v`
Expected: All 11 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/root_node.cpp src/engine/node/tissues/root_apical.cpp tests/test_water.cpp
git commit -m "feat: root nodes absorb water proportional to surface area"
```

---

### Task 7: Leaf Transpiration and Photosynthesis Water Cost

**Files:**
- Modify: `src/engine/node/tissues/leaf.cpp`
- Modify: `tests/test_water.cpp`

- [ ] **Step 1: Write failing tests for transpiration and photosynthesis water cost**

Append to `tests/test_water.cpp`:

```cpp
TEST_CASE("Leaves lose water via transpiration", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    plant.seed_mut()->add_child(leaf);

    float initial_water = 5.0f;
    leaf->chemical(ChemicalID::Water) = initial_water;
    leaf->chemical(ChemicalID::Sugar) = 10.0f;  // prevent starvation

    WorldParams wp = default_world_params();
    wp.light_level = 2.0f;
    compute_light_exposure(plant, wp);
    leaf->tick(plant, wp);

    REQUIRE(leaf->chemical(ChemicalID::Water) < initial_water);
}

TEST_CASE("Transpiration scales with light exposure", "[water]") {
    Genome g = default_genome();

    // Leaf in bright light
    Plant plant_bright(g, glm::vec3(0.0f));
    Node* leaf_bright = plant_bright.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf_bright->as_leaf()->leaf_size = 0.5f;
    plant_bright.seed_mut()->add_child(leaf_bright);
    leaf_bright->chemical(ChemicalID::Water) = 5.0f;
    leaf_bright->chemical(ChemicalID::Sugar) = 10.0f;

    WorldParams wp_bright = default_world_params();
    wp_bright.light_level = 5.0f;
    compute_light_exposure(plant_bright, wp_bright);
    leaf_bright->tick(plant_bright, wp_bright);

    // Leaf in dim light
    Plant plant_dim(g, glm::vec3(0.0f));
    Node* leaf_dim = plant_dim.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf_dim->as_leaf()->leaf_size = 0.5f;
    plant_dim.seed_mut()->add_child(leaf_dim);
    leaf_dim->chemical(ChemicalID::Water) = 5.0f;
    leaf_dim->chemical(ChemicalID::Sugar) = 10.0f;

    WorldParams wp_dim = default_world_params();
    wp_dim.light_level = 0.5f;
    compute_light_exposure(plant_dim, wp_dim);
    leaf_dim->tick(plant_dim, wp_dim);

    // Brighter leaf should lose more water
    REQUIRE(leaf_bright->chemical(ChemicalID::Water) < leaf_dim->chemical(ChemicalID::Water));
}

TEST_CASE("Water does not go below zero from transpiration", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 1.0f;
    plant.seed_mut()->add_child(leaf);

    leaf->chemical(ChemicalID::Water) = 0.0001f;  // almost zero
    leaf->chemical(ChemicalID::Sugar) = 10.0f;

    WorldParams wp = default_world_params();
    wp.light_level = 10.0f;  // very bright — high transpiration
    compute_light_exposure(plant, wp);
    leaf->tick(plant, wp);

    REQUIRE(leaf->chemical(ChemicalID::Water) >= 0.0f);
}

TEST_CASE("Photosynthesis consumes water proportional to sugar produced", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    plant.seed_mut()->add_child(leaf);

    // Give plenty of water so transpiration doesn't zero it out
    leaf->chemical(ChemicalID::Water) = 100.0f;
    leaf->chemical(ChemicalID::Sugar) = 0.0f;

    WorldParams wp = default_world_params();
    wp.light_level = 2.0f;
    compute_light_exposure(plant, wp);

    float water_before = leaf->chemical(ChemicalID::Water);
    leaf->tick(plant, wp);
    float water_after = leaf->chemical(ChemicalID::Water);
    float sugar_produced = leaf->chemical(ChemicalID::Sugar);  // approximate — some consumed by maintenance

    // Water loss should exceed transpiration alone (photosynthesis adds cost)
    float leaf_area = 0.5f * 0.5f;
    float transpiration_only = g.transpiration_rate * leaf_area * leaf->as_leaf()->light_exposure;
    float total_loss = water_before - water_after;

    // Total loss > transpiration alone (photosynthesis cost adds to it)
    REQUIRE(total_loss > transpiration_only - 1e-6f);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[water]" -v`
Expected: FAIL — leaves don't consume water yet.

- [ ] **Step 3: Add transpiration and photosynthesis water cost to `LeafNode`**

In `src/engine/node/tissues/leaf.cpp`, add `#include "engine/sugar.h"` if not already present.

Modify `photosynthesize()` to deduct water after sugar production. Add transpiration to `tissue_tick()` after the photosynthesis call.

In `tissue_tick()`, add transpiration after the photosynthesis block (after line 28):

```cpp
    // Transpiration: water loss proportional to leaf area and light exposure
    float leaf_area = leaf_size * leaf_size;
    float transpired = g.transpiration_rate * leaf_area * light_exposure;
    chemical(ChemicalID::Water) = std::max(0.0f, chemical(ChemicalID::Water) - transpired);
```

In `photosynthesize()`, add water cost after sugar production (after line 74 where sugar is clamped):

```cpp
    // Photosynthesis water cost: small deduction proportional to sugar produced
    float water_cost = sugar_produced * g.photosynthesis_water_ratio;
    chemical(ChemicalID::Water) = std::max(0.0f, chemical(ChemicalID::Water) - water_cost);
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[water]" -v`
Expected: All 15 tests PASS.

- [ ] **Step 5: Run full test suite**

Run: `./build/botany_tests -v`
Expected: All tests PASS.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/tissues/leaf.cpp tests/test_water.cpp
git commit -m "feat: leaf transpiration and photosynthesis water cost"
```

---

### Task 8: Seed Water Initialization

**Files:**
- Modify: `src/engine/plant.h` or `src/engine/plant.cpp` (wherever seed sugar is initialized)
- Modify: `tests/test_water.cpp`

- [ ] **Step 1: Find where seed sugar is initialized**

Search for `seed_sugar` in `plant.cpp` to find where the seed node gets its initial sugar. The seed needs initial water too.

- [ ] **Step 2: Write failing test**

Append to `tests/test_water.cpp`:

```cpp
TEST_CASE("Seed node starts with water", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    REQUIRE(plant.seed()->chemical(ChemicalID::Water) > 0.0f);
}
```

- [ ] **Step 3: Run test to verify it fails**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "Seed node starts with water" -v`
Expected: FAIL — seed water is 0.

- [ ] **Step 4: Initialize seed water alongside seed sugar**

Find the line where `seed_sugar` is set on the seed node and add a corresponding water initialization. A reasonable initial value is the seed's water cap (start full):

```cpp
seed->chemical(ChemicalID::Water) = water_cap(*seed, g);
```

- [ ] **Step 5: Run tests to verify it passes**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[water]" -v`
Expected: All 16 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add src/engine/plant.cpp tests/test_water.cpp
git commit -m "feat: seed node starts with full water capacity"
```

---

### Task 9: Evolution Bridge — Water Linkage Group

**Files:**
- Modify: `src/evolution/genome_bridge.cpp`

- [ ] **Step 1: Add water genes to `build_genome_template()`**

In `src/evolution/genome_bridge.cpp`, add the water economy gene registrations after the sugar economy group (after line 101) and before the gibberellin group:

```cpp
    // --- Water economy group (10 genes) ---
    reg(sg, "water_absorption_rate",       g.water_absorption_rate,       r, 0.001f,  0.5f, p);
    reg(sg, "transpiration_rate",          g.transpiration_rate,          r, 0.001f,  0.5f, p);
    reg(sg, "photosynthesis_water_ratio",  g.photosynthesis_water_ratio,  r, 0.01f,   5.0f, p);
    reg(sg, "water_storage_density_stem",  g.water_storage_density_stem,  r, 100.0f,  2000.0f, p);
    reg(sg, "water_storage_density_leaf",  g.water_storage_density_leaf,  r, 0.1f,    10.0f, p);
    reg(sg, "water_cap_meristem",          g.water_cap_meristem,          r, 0.1f,    10.0f, p);
    reg(sg, "water_diffusion_rate",        g.water_diffusion_rate,        r, 0.1f,    1.0f, p);
    reg(sg, "water_bias",                  g.water_bias,                  r, 0.0f,    0.5f, p);
    reg(sg, "water_base_transport",        g.water_base_transport,        r, 0.01f,   1.0f, p);
    reg(sg, "water_transport_scale",       g.water_transport_scale,       r, 0.5f,    20.0f, p);
```

- [ ] **Step 2: Add water linkage group**

Add after the sugar_economy linkage group (after line 196):

```cpp
    sg.add_linkage_group({"water_economy", {
        "water_absorption_rate", "transpiration_rate", "photosynthesis_water_ratio",
        "water_storage_density_stem", "water_storage_density_leaf", "water_cap_meristem",
        "water_diffusion_rate", "water_bias", "water_base_transport", "water_transport_scale"
    }});
```

- [ ] **Step 3: Add water fields to `from_structured()`**

In `from_structured()`, add after the sugar economy section (after line 314):

```cpp
    // Water economy
    g.water_absorption_rate       = sg.get("water_absorption_rate");
    g.transpiration_rate          = sg.get("transpiration_rate");
    g.photosynthesis_water_ratio  = sg.get("photosynthesis_water_ratio");
    g.water_storage_density_stem  = sg.get("water_storage_density_stem");
    g.water_storage_density_leaf  = sg.get("water_storage_density_leaf");
    g.water_cap_meristem          = sg.get("water_cap_meristem");
    g.water_diffusion_rate        = sg.get("water_diffusion_rate");
    g.water_bias                  = sg.get("water_bias");
    g.water_base_transport        = sg.get("water_base_transport");
    g.water_transport_scale       = sg.get("water_transport_scale");
```

- [ ] **Step 4: Build and run all tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests -v`
Expected: All tests PASS (including evolution tests).

- [ ] **Step 5: Commit**

```bash
git add src/evolution/genome_bridge.cpp
git commit -m "feat: add water params to evolution bridge with linkage group"
```

---

### Task 10: Renderer Water Color Mode

**Files:**
- Modify: `src/app_realtime.cpp:289-310` (--color parsing)
- Modify: `src/app_realtime.cpp:523-575` (ImGui overlay buttons)

- [ ] **Step 1: Add `water` to `--color` CLI argument parsing**

In `src/app_realtime.cpp`, add a case in the color_chemical parsing block (around line 302, after the ethylene case):

```cpp
            } else if (color_chemical == "water") {
                accessor = [](const Node& n) { return n.chemical(ChemicalID::Water); };
```

Update the error message to include `water`:

```cpp
                          << " (available: auxin, cytokinin, sugar, gibberellin, ethylene, water, type)" << std::endl;
```

- [ ] **Step 2: Add `WATER` to the `Overlay` enum**

In `src/app_realtime.cpp`, update the enum (line 323):

```cpp
    enum class Overlay { NONE, NODE_TYPE, AUXIN, CYTOKININ, SUGAR, LIGHT, GIBBERELLIN, ETHYLENE, STRESS, WATER };
```

- [ ] **Step 3: Add Water button to ImGui overlay panel**

In the overlay button section (around line 571, after the Stress button):

```cpp
            ImGui::SameLine();
            if (ImGui::Button("Water")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.chemical(ChemicalID::Water); });
                active_overlay = Overlay::WATER;
            }
```

- [ ] **Step 4: Build and verify**

Run: `/usr/local/bin/cmake --build build`
Expected: Clean build. Test with: `./build/botany_realtime --color water`

- [ ] **Step 5: Commit**

```bash
git add src/app_realtime.cpp
git commit -m "feat: add water heatmap color mode to renderer"
```

---

### Task 11: Update CLAUDE.md Documentation

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Add Water Model section to CLAUDE.md**

Add after the Canalization Model section:

```markdown
## Water Model

Water is a **persistent** capacity-based resource (like sugar):

- **Absorbed** by all root tissue (`RootNode` + `RootApicalNode`) proportional to surface area and `soil_moisture`
- **Lost** by leaves via transpiration (proportional to leaf area and light exposure) and photosynthesis water cost
- **Transported** through the tree graph via the unified transport system with a slight upward bias
- **No step-one effects** — water deficit does not gate growth or cause wilting (future enhancement)

Surface area: `2 * pi * r * length` for root segments, `2 * pi * r * r` (hemisphere) for root apical tips.
```

- [ ] **Step 2: Add water params to tuning section**

Add water parameters to the Tuning Parameters section with their default values and descriptions.

- [ ] **Step 3: Add `soil_moisture` to WorldParams description**

In the World Params section, add `soil_moisture` (default 1.0, range 0-1).

- [ ] **Step 4: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: add water model documentation to CLAUDE.md"
```
