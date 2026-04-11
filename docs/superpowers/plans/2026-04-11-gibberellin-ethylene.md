# Gibberellin & Ethylene Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add gibberellin (local internode elongation hormone) and ethylene (spatial gas stress signal with leaf abscission) to the plant growth engine.

**Architecture:** Two new hormone systems as standalone .h/.cpp pairs, following the existing hormone.h/sugar.h pattern. GA is reset-and-recompute with local-only transport (leaf → parent node). Ethylene is reset-and-recompute with spatial proximity diffusion (O(n^2) brute-force). Both integrate into the engine tick loop and modify existing intercalary elongation. Ethylene adds a new abscission phase that removes senesced leaves.

**Tech Stack:** C++17, GLM, Catch2 (tests), OpenGL/GLFW/ImGui (renderer)

**Spec:** `docs/superpowers/specs/2026-04-10-gibberellin-ethylene-design.md`

---

### Task 1: Add Node Fields

**Files:**
- Modify: `src/engine/node.h:30-53`
- Modify: `src/engine/node.cpp:5-16`

- [ ] **Step 1: Add gibberellin, ethylene, senescence_ticks fields to Node**

In `src/engine/node.h`, add three new fields after `starvation_ticks`:

```cpp
uint32_t starvation_ticks = 0;
float gibberellin = 0.0f;            // GA concentration (reset each tick)
float ethylene = 0.0f;               // ethylene concentration (reset each tick)
uint32_t senescence_ticks = 0;       // 0 = healthy, >0 = senescing (irreversible)
```

- [ ] **Step 2: Initialize new fields in Node constructor**

In `src/engine/node.cpp`, the constructor already value-initializes `sugar`, `leaf_size`, `light_exposure`, `starvation_ticks` via in-class defaults. The new fields use in-class defaults too (= 0.0f / = 0), so no constructor change is needed. Verify this compiles.

Run: `/usr/local/bin/cmake --build build 2>&1 | tail -5`
Expected: Build succeeds

- [ ] **Step 3: Commit**

```bash
git add src/engine/node.h
git commit -m "feat: add gibberellin, ethylene, senescence_ticks fields to Node"
```

---

### Task 2: Add Genome Parameters

**Files:**
- Modify: `src/engine/genome.h:11-78` (struct) and `src/engine/genome.h:80-139` (default_genome)

- [ ] **Step 1: Add GA genome parameters**

In `src/engine/genome.h`, add after the `sugar_activation_root` field (line 77):

```cpp
    // Gibberellin — promotes internode elongation, produced by young leaves
    float ga_production_rate;         // GA produced per dm of leaf_size per tick
    uint32_t ga_leaf_age_max;         // ticks — only leaves younger than this produce GA
    float ga_elongation_sensitivity;  // how strongly GA boosts elongation rate
    float ga_length_sensitivity;      // how strongly GA boosts target internode length
```

- [ ] **Step 2: Add ethylene genome parameters**

In `src/engine/genome.h`, add after the GA parameters:

```cpp
    // Ethylene — stress/crowding gas signal, triggers leaf abscission
    float ethylene_starvation_rate;       // production when sugar = 0
    float ethylene_shade_rate;            // production from low light
    float ethylene_shade_threshold;       // light_exposure below which shade-ethylene kicks in
    float ethylene_age_rate;              // production ramp from old age
    uint32_t ethylene_age_onset;          // tick age when age-ethylene starts
    float ethylene_crowding_rate;         // production per nearby node
    float ethylene_crowding_radius;       // dm — radius for crowding density count
    float ethylene_diffusion_radius;      // dm — gas cloud spread distance
    float ethylene_abscission_threshold;  // ethylene level triggering leaf senescence
    float ethylene_elongation_inhibition; // strength of elongation suppression
    uint32_t senescence_duration;         // ticks from senescence start to leaf drop
```

- [ ] **Step 3: Add defaults to default_genome()**

In `default_genome()`, add after `.sugar_activation_root = 0.3f,`:

```cpp
        // Gibberellin
        .ga_production_rate = 0.5f,
        .ga_leaf_age_max = 168,               // 7 days
        .ga_elongation_sensitivity = 2.0f,
        .ga_length_sensitivity = 1.5f,

        // Ethylene
        .ethylene_starvation_rate = 0.3f,
        .ethylene_shade_rate = 0.2f,
        .ethylene_shade_threshold = 0.3f,
        .ethylene_age_rate = 0.05f,
        .ethylene_age_onset = 720,            // 30 days
        .ethylene_crowding_rate = 0.1f,
        .ethylene_crowding_radius = 0.5f,     // dm
        .ethylene_diffusion_radius = 1.0f,    // dm
        .ethylene_abscission_threshold = 0.5f,
        .ethylene_elongation_inhibition = 1.0f,
        .senescence_duration = 48,            // 2 days
```

- [ ] **Step 4: Build to verify**

Run: `/usr/local/bin/cmake --build build 2>&1 | tail -5`
Expected: Build succeeds

- [ ] **Step 5: Commit**

```bash
git add src/engine/genome.h
git commit -m "feat: add gibberellin and ethylene genome parameters"
```

---

### Task 3: Gibberellin Computation + Tests

**Files:**
- Create: `src/engine/gibberellin.h`
- Create: `src/engine/gibberellin.cpp`
- Create: `tests/test_gibberellin.cpp`
- Modify: `CMakeLists.txt:69-81` (engine lib) and `CMakeLists.txt:116-125` (tests)

- [ ] **Step 1: Add source files to CMakeLists.txt**

In `CMakeLists.txt`, add `src/engine/gibberellin.cpp` to the `botany_engine` library sources (after `src/engine/sugar.cpp`):

```
    src/engine/gibberellin.cpp
```

Add `tests/test_gibberellin.cpp` to the `botany_tests` sources (after `tests/test_sugar.cpp`):

```
    tests/test_gibberellin.cpp
```

- [ ] **Step 2: Create gibberellin.h**

Create `src/engine/gibberellin.h`:

```cpp
#pragma once

namespace botany {

class Plant;

void compute_gibberellin(Plant& plant);

} // namespace botany
```

- [ ] **Step 3: Write failing tests**

Create `tests/test_gibberellin.cpp`:

```cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/gibberellin.h"

using namespace botany;
using Catch::Matchers::WithinAbs;

TEST_CASE("Young leaf produces GA on parent node", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Create: seed -> stem -> leaf
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    leaf->leaf_size = 0.2f;
    leaf->age = 10; // young
    stem->add_child(leaf);

    compute_gibberellin(plant);

    float expected = leaf->leaf_size * g.ga_production_rate;
    REQUIRE_THAT(stem->gibberellin, WithinAbs(expected, 1e-6));
}

TEST_CASE("Old leaf produces no GA", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    leaf->leaf_size = 0.2f;
    leaf->age = g.ga_leaf_age_max + 1; // too old
    stem->add_child(leaf);

    compute_gibberellin(plant);

    REQUIRE(stem->gibberellin == 0.0f);
}

TEST_CASE("GA reaches grandparent at reduced fraction", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // seed -> grandparent -> parent -> leaf
    Node* grandparent = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(grandparent);

    Node* parent = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    grandparent->add_child(parent);

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    leaf->leaf_size = 0.2f;
    leaf->age = 10;
    parent->add_child(leaf);

    compute_gibberellin(plant);

    float parent_expected = leaf->leaf_size * g.ga_production_rate;
    float grandparent_expected = leaf->leaf_size * g.ga_production_rate * 0.3f;

    REQUIRE_THAT(parent->gibberellin, WithinAbs(parent_expected, 1e-6));
    REQUIRE_THAT(grandparent->gibberellin, WithinAbs(grandparent_expected, 1e-6));
}

TEST_CASE("GA does not spread beyond grandparent", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // seed -> great_grandparent -> grandparent -> parent -> leaf
    Node* ggp = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(ggp);

    Node* gp = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    ggp->add_child(gp);

    Node* parent = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    gp->add_child(parent);

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    leaf->leaf_size = 0.2f;
    leaf->age = 10;
    parent->add_child(leaf);

    compute_gibberellin(plant);

    // Great-grandparent and seed should have no GA
    REQUIRE(ggp->gibberellin == 0.0f);
    REQUIRE(plant.seed()->gibberellin == 0.0f);
}

TEST_CASE("GA resets to zero before recomputing", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);

    // Manually set GA to something nonzero
    stem->gibberellin = 999.0f;

    // No leaves — compute should reset to 0
    compute_gibberellin(plant);

    REQUIRE(stem->gibberellin == 0.0f);
}
```

- [ ] **Step 4: Run tests to verify they fail**

Run: `/usr/local/bin/cmake --build build 2>&1 | tail -3`
Expected: Linker error — `compute_gibberellin` not defined

- [ ] **Step 5: Implement compute_gibberellin**

Create `src/engine/gibberellin.cpp`:

```cpp
#include "engine/gibberellin.h"
#include "engine/plant.h"
#include "engine/node.h"

namespace botany {

void compute_gibberellin(Plant& plant) {
    const Genome& g = plant.genome();

    // Phase 1: Reset all GA to zero
    plant.for_each_node_mut([](Node& node) {
        node.gibberellin = 0.0f;
    });

    // Phase 2: Young leaves produce GA on their parent (and grandparent)
    plant.for_each_node_mut([&](Node& node) {
        if (node.type != NodeType::LEAF) return;
        if (node.age >= g.ga_leaf_age_max) return;
        if (node.leaf_size < 1e-6f) return;

        float production = node.leaf_size * g.ga_production_rate;

        // Apply to parent
        if (node.parent) {
            node.parent->gibberellin += production;

            // Apply reduced fraction to grandparent
            if (node.parent->parent) {
                node.parent->parent->gibberellin += production * 0.3f;
            }
        }
    });
}

} // namespace botany
```

- [ ] **Step 6: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[gibberellin]" -v`
Expected: All 5 gibberellin tests pass

- [ ] **Step 7: Commit**

```bash
git add src/engine/gibberellin.h src/engine/gibberellin.cpp tests/test_gibberellin.cpp CMakeLists.txt
git commit -m "feat: gibberellin hormone — local production by young leaves"
```

---

### Task 4: Ethylene Computation + Tests

**Files:**
- Create: `src/engine/ethylene.h`
- Create: `src/engine/ethylene.cpp`
- Create: `tests/test_ethylene.cpp`
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Add source files to CMakeLists.txt**

Add `src/engine/ethylene.cpp` to `botany_engine` (after `gibberellin.cpp`).
Add `tests/test_ethylene.cpp` to `botany_tests` (after `test_gibberellin.cpp`).

- [ ] **Step 2: Create ethylene.h**

Create `src/engine/ethylene.h`:

```cpp
#pragma once

namespace botany {

class Plant;
struct WorldParams;

void compute_ethylene(Plant& plant, const WorldParams& world);
void process_abscission(Plant& plant);

} // namespace botany
```

- [ ] **Step 3: Write failing tests for ethylene production triggers**

Create `tests/test_ethylene.cpp`:

```cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/ethylene.h"

using namespace botany;
using Catch::Matchers::WithinAbs;

TEST_CASE("Starvation produces ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);
    stem->sugar = 0.0f;

    // Place far from other nodes so spatial diffusion doesn't interfere
    stem->position = glm::vec3(100.0f, 100.0f, 100.0f);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    REQUIRE_THAT(stem->ethylene, WithinAbs(g.ethylene_starvation_rate, 0.1));
}

TEST_CASE("Fed node produces no starvation ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);
    stem->sugar = 5.0f;
    stem->position = glm::vec3(100.0f, 100.0f, 100.0f);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    // No starvation trigger, no shade trigger (not a leaf), no age trigger (not a leaf),
    // no crowding (isolated). But may receive spatial diffusion from other starved nodes.
    // This stem is isolated so should be ~0.
    REQUIRE(stem->ethylene < 0.01f);
}

TEST_CASE("Shaded leaf produces ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.2f;
    leaf->sugar = 1.0f; // not starved
    leaf->light_exposure = 0.1f; // heavily shaded
    leaf->position = glm::vec3(100.0f, 100.0f, 100.0f);
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    float expected = g.ethylene_shade_rate * (1.0f - 0.1f);
    REQUIRE(leaf->ethylene > expected * 0.5f); // at least significant shade production
}

TEST_CASE("Well-lit leaf produces no shade ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.2f;
    leaf->sugar = 1.0f;
    leaf->light_exposure = 0.8f; // above shade threshold
    leaf->position = glm::vec3(100.0f, 100.0f, 100.0f);
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    // No shade ethylene (0.8 > 0.3 threshold), no starvation, young age, isolated
    REQUIRE(leaf->ethylene < 0.01f);
}

TEST_CASE("Old leaf produces age ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.2f;
    leaf->sugar = 1.0f;
    leaf->light_exposure = 1.0f;
    leaf->age = g.ethylene_age_onset + 360; // 15 days past onset
    leaf->position = glm::vec3(100.0f, 100.0f, 100.0f);
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    float expected = g.ethylene_age_rate * 360.0f / static_cast<float>(g.ethylene_age_onset);
    REQUIRE(leaf->ethylene > expected * 0.5f);
}

TEST_CASE("Crowded nodes produce ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Create a cluster of nodes very close together
    glm::vec3 center(5.0f, 5.0f, 5.0f);
    Node* target = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
    target->sugar = 1.0f;
    target->position = center;
    plant.seed_mut()->add_child(target);

    for (int i = 0; i < 5; i++) {
        Node* n = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
        n->sugar = 1.0f;
        n->position = center + glm::vec3(0.1f * i, 0.0f, 0.0f);
        plant.seed_mut()->add_child(n);
    }

    plant.recompute_world_positions();

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    // Target should have crowding ethylene since 5 nodes are within 0.5 dm
    REQUIRE(target->ethylene > 0.1f);
}

TEST_CASE("Spatial diffusion spreads ethylene to nearby nodes", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Source: starving node producing ethylene
    Node* source = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
    source->sugar = 0.0f;
    source->position = glm::vec3(10.0f, 10.0f, 10.0f);
    plant.seed_mut()->add_child(source);

    // Nearby: within diffusion radius (1.0 dm)
    Node* nearby = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
    nearby->sugar = 1.0f;
    nearby->position = glm::vec3(10.5f, 10.0f, 10.0f); // 0.5 dm away
    plant.seed_mut()->add_child(nearby);

    // Far: outside diffusion radius
    Node* far = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
    far->sugar = 1.0f;
    far->position = glm::vec3(15.0f, 10.0f, 10.0f); // 5 dm away
    plant.seed_mut()->add_child(far);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    REQUIRE(nearby->ethylene > 0.0f); // received diffused ethylene
    REQUIRE(far->ethylene < 0.01f);   // too far, nothing received
}

TEST_CASE("Ethylene resets to zero before recomputing", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Set ethylene manually, place nodes far away
    plant.seed_mut()->ethylene = 999.0f;
    plant.seed_mut()->sugar = 1.0f; // not starved
    plant.seed_mut()->position = glm::vec3(100.0f, 100.0f, 100.0f);

    // Remove children so no triggers fire
    // (seed has shoot + root children that might trigger starvation)
    plant.for_each_node_mut([](Node& n) {
        n.sugar = 1.0f; // prevent starvation trigger
    });

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    // Ethylene should have been reset, not accumulated
    REQUIRE(plant.seed()->ethylene < 999.0f);
}
```

- [ ] **Step 4: Run tests to verify they fail**

Run: `/usr/local/bin/cmake --build build 2>&1 | tail -3`
Expected: Linker error — `compute_ethylene` not defined

- [ ] **Step 5: Implement compute_ethylene and process_abscission**

Create `src/engine/ethylene.cpp`:

```cpp
#include "engine/ethylene.h"
#include "engine/plant.h"
#include "engine/node.h"
#include "engine/world_params.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <glm/geometric.hpp>

namespace botany {

void compute_ethylene(Plant& plant, const WorldParams& /*world*/) {
    const Genome& g = plant.genome();

    // Collect node pointers and positions for spatial queries
    struct NodeInfo {
        Node* node;
        glm::vec3 pos;
    };
    std::vector<NodeInfo> nodes;
    plant.for_each_node_mut([&](Node& node) {
        nodes.push_back({&node, node.position});
    });

    // Phase 1: Reset and compute local production
    for (auto& info : nodes) {
        Node& node = *info.node;
        node.ethylene = 0.0f;

        // Trigger 1: Sugar starvation
        if (node.sugar <= 0.0f) {
            node.ethylene += g.ethylene_starvation_rate;
        }

        // Trigger 2: Low light (LEAF only)
        if (node.type == NodeType::LEAF &&
            node.light_exposure < g.ethylene_shade_threshold) {
            node.ethylene += g.ethylene_shade_rate * (1.0f - node.light_exposure);
        }

        // Trigger 3: Old age (LEAF only)
        if (node.type == NodeType::LEAF &&
            node.age > g.ethylene_age_onset) {
            float age_past = static_cast<float>(node.age - g.ethylene_age_onset);
            node.ethylene += g.ethylene_age_rate * age_past
                           / static_cast<float>(g.ethylene_age_onset);
        }

        // Trigger 4: Crowding — count nearby nodes within crowding_radius
        float cr2 = g.ethylene_crowding_radius * g.ethylene_crowding_radius;
        int nearby_count = 0;
        for (const auto& other : nodes) {
            if (other.node == &node) continue;
            glm::vec3 diff = other.pos - info.pos;
            if (glm::dot(diff, diff) < cr2) {
                nearby_count++;
            }
        }
        node.ethylene += g.ethylene_crowding_rate * static_cast<float>(nearby_count);
    }

    // Phase 2: Spatial gas diffusion (compute-then-apply)
    std::vector<float> received(nodes.size(), 0.0f);
    float dr = g.ethylene_diffusion_radius;
    float dr2 = dr * dr;

    for (size_t i = 0; i < nodes.size(); i++) {
        if (nodes[i].node->ethylene <= 0.0f) continue;
        for (size_t j = 0; j < nodes.size(); j++) {
            if (i == j) continue;
            glm::vec3 diff = nodes[j].pos - nodes[i].pos;
            float dist2 = glm::dot(diff, diff);
            if (dist2 >= dr2) continue;
            float dist = std::sqrt(dist2);
            float falloff = 1.0f - dist / dr;
            received[j] += nodes[i].node->ethylene * falloff;
        }
    }

    for (size_t i = 0; i < nodes.size(); i++) {
        nodes[i].node->ethylene += received[i];
    }
}

void process_abscission(Plant& plant) {
    const Genome& g = plant.genome();

    // Increment senescence on senescing leaves, collect leaves to remove
    std::vector<Node*> to_remove;
    plant.for_each_node_mut([&](Node& node) {
        if (node.type != NodeType::LEAF) return;

        // Start senescence if ethylene exceeds threshold and not yet senescing
        if (node.senescence_ticks == 0 &&
            node.ethylene > g.ethylene_abscission_threshold) {
            node.senescence_ticks = 1;
        }

        // Advance senescence
        if (node.senescence_ticks > 0) {
            node.senescence_ticks++;
            if (node.senescence_ticks >= g.senescence_duration) {
                to_remove.push_back(&node);
            }
        }
    });

    for (Node* n : to_remove) {
        plant.remove_subtree(n);
    }
}

} // namespace botany
```

- [ ] **Step 6: Build and run ethylene tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[ethylene]" -v`
Expected: All 8 ethylene tests pass

- [ ] **Step 7: Commit**

```bash
git add src/engine/ethylene.h src/engine/ethylene.cpp tests/test_ethylene.cpp CMakeLists.txt
git commit -m "feat: ethylene hormone — spatial gas diffusion with four stress triggers"
```

---

### Task 5: Abscission Tests

**Files:**
- Modify: `tests/test_ethylene.cpp`

- [ ] **Step 1: Add abscission lifecycle tests**

Append to `tests/test_ethylene.cpp`:

```cpp
TEST_CASE("Leaf above ethylene threshold begins senescence", "[ethylene][abscission]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.2f;
    leaf->ethylene = g.ethylene_abscission_threshold + 0.1f;
    plant.seed_mut()->add_child(leaf);

    process_abscission(plant);

    REQUIRE(leaf->senescence_ticks > 0);
}

TEST_CASE("Leaf below ethylene threshold stays healthy", "[ethylene][abscission]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.2f;
    leaf->ethylene = g.ethylene_abscission_threshold * 0.5f;
    plant.seed_mut()->add_child(leaf);

    process_abscission(plant);

    REQUIRE(leaf->senescence_ticks == 0);
}

TEST_CASE("Senescing leaf is removed after senescence_duration", "[ethylene][abscission]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.2f;
    leaf->ethylene = g.ethylene_abscission_threshold + 0.1f;
    leaf->senescence_ticks = g.senescence_duration - 1; // almost done
    plant.seed_mut()->add_child(leaf);

    uint32_t count_before = plant.node_count();
    process_abscission(plant);

    REQUIRE(plant.node_count() < count_before);
}

TEST_CASE("Non-leaf nodes do not senesce", "[ethylene][abscission]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    stem->ethylene = g.ethylene_abscission_threshold + 10.0f; // very high
    plant.seed_mut()->add_child(stem);

    uint32_t count_before = plant.node_count();
    process_abscission(plant);

    // Stem should not be removed and should not have senescence_ticks
    REQUIRE(stem->senescence_ticks == 0);
    REQUIRE(plant.node_count() == count_before);
}
```

- [ ] **Step 2: Run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[abscission]" -v`
Expected: All 4 abscission tests pass

- [ ] **Step 3: Commit**

```bash
git add tests/test_ethylene.cpp
git commit -m "test: add abscission lifecycle tests"
```

---

### Task 6: Engine Tick Integration

**Files:**
- Modify: `src/engine/engine.cpp:1-24`

- [ ] **Step 1: Wire GA, ethylene, and abscission into tick loop**

In `src/engine/engine.cpp`, add includes:

```cpp
#include "engine/gibberellin.h"
#include "engine/ethylene.h"
```

Replace the tick loop body (lines 17-22) with:

```cpp
        transport_auxin(*plant);
        transport_cytokinin(*plant);
        compute_gibberellin(*plant);
        transport_sugar(*plant, world_params_);
        compute_ethylene(*plant, world_params_);
        process_abscission(*plant);
        tick_meristems(*plant, world_params_);
        plant->recompute_world_positions();
```

- [ ] **Step 2: Build and run all tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests -v`
Expected: All tests pass (gibberellin, ethylene, and all existing tests)

- [ ] **Step 3: Commit**

```bash
git add src/engine/engine.cpp
git commit -m "feat: wire gibberellin, ethylene, abscission into engine tick loop"
```

---

### Task 7: GA/Ethylene Elongation Modifiers

**Files:**
- Modify: `src/engine/meristems/meristem.cpp:34-56`
- Modify: `tests/test_meristem.cpp`

- [ ] **Step 1: Write failing test for GA-boosted elongation**

Add to `tests/test_meristem.cpp`:

```cpp
TEST_CASE("GA boosts intercalary elongation rate", "[meristem][gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
    WorldParams wp = default_world_params();

    // Create an interior stem node eligible for elongation
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    stem->age = 1; // young enough
    stem->sugar = 5.0f;
    plant.seed_mut()->add_child(stem);

    // Run without GA
    stem->gibberellin = 0.0f;
    float offset_before = glm::length(stem->offset);
    tick_meristems(plant, wp);
    float growth_no_ga = glm::length(stem->offset) - offset_before;

    // Reset
    stem->offset = glm::vec3(0.0f, 0.5f, 0.0f);
    stem->age = 1;
    stem->sugar = 5.0f;

    // Run with GA
    stem->gibberellin = 1.0f;
    offset_before = glm::length(stem->offset);
    tick_meristems(plant, wp);
    float growth_with_ga = glm::length(stem->offset) - offset_before;

    REQUIRE(growth_with_ga > growth_no_ga);
}

TEST_CASE("Ethylene inhibits elongation", "[meristem][ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
    WorldParams wp = default_world_params();

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    stem->age = 1;
    stem->sugar = 5.0f;
    plant.seed_mut()->add_child(stem);

    // Run without ethylene
    stem->ethylene = 0.0f;
    float offset_before = glm::length(stem->offset);
    tick_meristems(plant, wp);
    float growth_no_eth = glm::length(stem->offset) - offset_before;

    // Reset
    stem->offset = glm::vec3(0.0f, 0.5f, 0.0f);
    stem->age = 1;
    stem->sugar = 5.0f;

    // Run with high ethylene
    stem->ethylene = 2.0f;
    offset_before = glm::length(stem->offset);
    tick_meristems(plant, wp);
    float growth_with_eth = glm::length(stem->offset) - offset_before;

    REQUIRE(growth_with_eth < growth_no_eth);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[meristem][gibberellin]" -v`
Expected: Test fails (GA has no effect yet on elongation)

- [ ] **Step 3: Modify intercalary elongation in meristem.cpp**

In `src/engine/meristems/meristem.cpp`, replace the intercalary growth block (lines 36-56) with:

```cpp
        // Intercalary growth: young interior nodes elongate their internode.
        // Eligible = no meristem (not a tip or dormant bud), not a leaf, has parent, young enough.
        if (!n.meristem && n.type != NodeType::LEAF && n.parent) {
            float elong_rate = (n.type == NodeType::STEM) ? g.internode_elongation_rate
                                                          : g.root_internode_elongation_rate;
            uint32_t mat_ticks = (n.type == NodeType::STEM) ? g.internode_maturation_ticks
                                                            : g.root_internode_maturation_ticks;
            float save = (n.type == NodeType::STEM) ? g.sugar_save_stem : g.sugar_save_root;

            if (n.age < mat_ticks && elong_rate > 1e-8f) {
                // GA boosts elongation rate and max target length
                float ga_boost = 1.0f + n.gibberellin * g.ga_elongation_sensitivity;
                // Ethylene inhibits elongation
                float eth_inhibit = std::max(0.0f, 1.0f - n.ethylene * g.ethylene_elongation_inhibition);
                float effective_rate = elong_rate * ga_boost * eth_inhibit;

                // GA-modulated max internode length
                float max_len = (n.type == NodeType::STEM) ? g.max_internode_length
                                                           : g.root_max_internode_length;
                max_len *= (1.0f + n.gibberellin * g.ga_length_sensitivity);
                float current_len = glm::length(n.offset);
                if (current_len >= max_len) continue;

                float max_cost = effective_rate * world.sugar_cost_elongation;
                float gf = sugar_growth_fraction(n.sugar, save, max_cost);
                if (gf > 1e-6f) {
                    float actual_rate = effective_rate * gf;
                    float actual_cost = actual_rate * world.sugar_cost_elongation;
                    n.sugar -= actual_cost;
                    float len = glm::length(n.offset);
                    if (len > 1e-4f) {
                        n.offset += (n.offset / len) * actual_rate;
                    }
                }
            }
        }
```

Add `#include <algorithm>` to the includes in `meristem.cpp` if not already present (for `std::max`).

- [ ] **Step 4: Run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[meristem]" -v`
Expected: All meristem tests pass, including the two new GA/ethylene tests

- [ ] **Step 5: Commit**

```bash
git add src/engine/meristems/meristem.cpp tests/test_meristem.cpp
git commit -m "feat: GA boosts and ethylene inhibits intercalary elongation"
```

---

### Task 8: Senescing Leaves Stop Producing Sugar

**Files:**
- Modify: `src/engine/sugar.cpp:55-95`
- Modify: `tests/test_sugar.cpp`

- [ ] **Step 1: Write failing test**

Add to `tests/test_sugar.cpp`:

```cpp
TEST_CASE("Senescing leaf produces no sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.5f;
    leaf->senescence_ticks = 1; // senescing
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 2.0f;

    produce_sugar(plant, wp);

    REQUIRE(leaf->sugar == 0.0f);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "Senescing leaf produces no sugar" -v`
Expected: FAIL — senescing leaf still produces sugar

- [ ] **Step 3: Add senescence check to produce_sugar**

In `src/engine/sugar.cpp`, in the `produce_sugar` function, inside the leaf collection block (line 55), change:

```cpp
        if (node.type == NodeType::LEAF && node.leaf_size > 1e-6f) {
```

to:

```cpp
        if (node.type == NodeType::LEAF && node.leaf_size > 1e-6f && node.senescence_ticks == 0) {
```

This excludes senescing leaves from both shadow-casting as producers and from the production loop.

- [ ] **Step 4: Run test**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "Senescing leaf produces no sugar" -v`
Expected: PASS

- [ ] **Step 5: Run all tests to verify no regressions**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests -v`
Expected: All tests pass

- [ ] **Step 6: Commit**

```bash
git add src/engine/sugar.cpp tests/test_sugar.cpp
git commit -m "feat: senescing leaves stop producing sugar"
```

---

### Task 9: Renderer — Senescence Coloring + New Overlays

**Files:**
- Modify: `src/renderer/renderer.cpp:112-125` (leaf coloring)
- Modify: `src/app_realtime.cpp:106-125` (CLI color modes) and `src/app_realtime.cpp:137` (Overlay enum) and `src/app_realtime.cpp:286-346` (overlay buttons)

- [ ] **Step 1: Add senescence leaf coloring in renderer**

In `src/renderer/renderer.cpp`, in `draw_plant()`, replace the leaf coloring block (lines 113-125):

```cpp
        if (node.type == NodeType::LEAF) {
            glm::vec3 dir = glm::normalize(node.position - node.parent->position);
            // Sunlit leaves are bright green, shaded leaves darken
            glm::vec3 sun_color(0.2f, 0.6f, 0.15f);
            glm::vec3 shade_color(0.08f, 0.25f, 0.06f);
            glm::vec3 leaf_color = glm::mix(shade_color, sun_color, node.light_exposure);
            if (node.starvation_ticks > 0) {
                float stress = static_cast<float>(node.starvation_ticks) / 50.0f;
                stress = glm::clamp(stress, 0.0f, 1.0f);
                glm::vec3 dead_color(0.4f, 0.3f, 0.1f);
                leaf_color = glm::mix(leaf_color, dead_color, stress);
            }
            draw_leaf(node.position, dir, node.leaf_size, leaf_color);
            return;
        }
```

with:

```cpp
        if (node.type == NodeType::LEAF) {
            glm::vec3 dir = glm::normalize(node.position - node.parent->position);
            glm::vec3 leaf_color;
            if (chemical_accessor_) {
                float v = chemical_accessor_(node);
                leaf_color = heatmap(v / max_val);
            } else {
                // Sunlit leaves are bright green, shaded leaves darken
                glm::vec3 sun_color(0.2f, 0.6f, 0.15f);
                glm::vec3 shade_color(0.08f, 0.25f, 0.06f);
                leaf_color = glm::mix(shade_color, sun_color, node.light_exposure);
                if (node.starvation_ticks > 0) {
                    float stress = static_cast<float>(node.starvation_ticks) / 50.0f;
                    stress = glm::clamp(stress, 0.0f, 1.0f);
                    glm::vec3 dead_color(0.4f, 0.3f, 0.1f);
                    leaf_color = glm::mix(leaf_color, dead_color, stress);
                }
            }
            // Senescence: green -> yellow -> brown (overrides other coloring)
            if (node.senescence_ticks > 0) {
                float progress = static_cast<float>(node.senescence_ticks) / 48.0f;
                progress = glm::clamp(progress, 0.0f, 1.0f);
                glm::vec3 yellow(0.8f, 0.7f, 0.1f);
                glm::vec3 brown(0.4f, 0.25f, 0.05f);
                glm::vec3 senesce_color = glm::mix(yellow, brown, progress);
                leaf_color = glm::mix(leaf_color, senesce_color, progress);
            }
            draw_leaf(node.position, dir, node.leaf_size, leaf_color);
            return;
        }
```

- [ ] **Step 2: Add gibberellin and ethylene CLI color modes**

In `src/app_realtime.cpp`, add after the `sugar` accessor block (around line 115):

```cpp
            } else if (color_chemical == "gibberellin") {
                accessor = [](const Node& n) { return n.gibberellin; };
            } else if (color_chemical == "ethylene") {
                accessor = [](const Node& n) { return n.ethylene; };
```

Update the error message (line 118-119) to include the new options:

```cpp
                std::cerr << "Unknown chemical: " << color_chemical
                          << " (available: auxin, cytokinin, sugar, gibberellin, ethylene, type)" << std::endl;
```

- [ ] **Step 3: Add overlay buttons in ImGui**

In `src/app_realtime.cpp`, update the `Overlay` enum (line 137):

```cpp
    enum class Overlay { NONE, NODE_TYPE, AUXIN, CYTOKININ, SUGAR, LIGHT, GIBBERELLIN, ETHYLENE };
```

After the "Light" button block (around line 321), add:

```cpp
            if (ImGui::Button("GA")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.gibberellin; });
                active_overlay = Overlay::GIBBERELLIN;
            }
            ImGui::SameLine();
            if (ImGui::Button("Ethylene")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.ethylene; });
                active_overlay = Overlay::ETHYLENE;
            }
```

- [ ] **Step 4: Add GA/ethylene to node inspector**

In `src/app_realtime.cpp`, in the chemicals table (around line 447), add two new rows after the cytokinin row:

```cpp
                    // Gibberellin
                    float parent_ga = sel.parent ? sel.parent->gibberellin : 0.0f;
                    float child_ga = 0.0f;
                    if (!sel.children.empty()) {
                        for (const Node* c : sel.children) child_ga += c->gibberellin;
                        child_ga /= static_cast<float>(sel.children.size());
                    }
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("GA (au)");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%s", ratio_str(parent_ga, sel.gibberellin).c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%.4f", sel.gibberellin);
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%s", ratio_str(child_ga, sel.gibberellin).c_str());

                    // Ethylene
                    float parent_eth = sel.parent ? sel.parent->ethylene : 0.0f;
                    float child_eth = 0.0f;
                    if (!sel.children.empty()) {
                        for (const Node* c : sel.children) child_eth += c->ethylene;
                        child_eth /= static_cast<float>(sel.children.size());
                    }
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("Eth (au)");
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%s", ratio_str(parent_eth, sel.ethylene).c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%.4f", sel.ethylene);
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%s", ratio_str(child_eth, sel.ethylene).c_str());
```

Also add senescence display after "Starvation" (around line 399):

```cpp
                if (sel.senescence_ticks > 0) {
                    ImGui::Text("Senescence: %u / %u ticks", sel.senescence_ticks, g.senescence_duration);
                }
```

(This requires access to the genome — add `const Genome& g = engine.get_plant(plant_id).genome();` at the start of the node inspector block if not already available.)

- [ ] **Step 5: Build**

Run: `/usr/local/bin/cmake --build build 2>&1 | tail -5`
Expected: Build succeeds

- [ ] **Step 6: Commit**

```bash
git add src/renderer/renderer.cpp src/app_realtime.cpp
git commit -m "feat: renderer senescence coloring + gibberellin/ethylene overlays"
```

---

### Task 10: Update CLAUDE.md

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Update architecture documentation**

Add to the Engine section in CLAUDE.md:

- `gibberellin.h/cpp` — `compute_gibberellin()` — local GA production by young leaves
- `ethylene.h/cpp` — `compute_ethylene()` + `process_abscission()` — spatial gas diffusion, leaf abscission

Update the tick order:
```
auxin → cytokinin → GA → sugar → ethylene → abscission → meristems → positions
```

Add a Gibberellin Model section and an Ethylene Model section describing the biology and mechanics.

Add the new genome parameters to the Tuning Parameters section.

Add `gibberellin` and `ethylene` to the `--color` options for the realtime viewer.

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: add gibberellin and ethylene to CLAUDE.md"
```

---

### Task 11: Full Build + Test Verification

- [ ] **Step 1: Clean build**

Run: `/usr/local/bin/cmake --build build --clean-first 2>&1 | tail -5`
Expected: Build succeeds with no warnings

- [ ] **Step 2: Run full test suite**

Run: `./build/botany_tests -v`
Expected: All tests pass — gibberellin, ethylene, abscission, meristem modifiers, sugar senescence, plus all original tests

- [ ] **Step 3: Verify realtime viewer launches**

Run: `./build/botany_realtime --color gibberellin`
Expected: Window opens, plant grows, GA heatmap visible on stem nodes near young leaves

Run: `./build/botany_realtime --color ethylene`
Expected: Window opens, ethylene heatmap visible on stressed/shaded/crowded nodes

---

### Task 12: Self-Thinning Integration Test

**Files:**
- Modify: `tests/test_ethylene.cpp`

- [ ] **Step 1: Add self-thinning cascade test**

Append to `tests/test_ethylene.cpp`:

```cpp
TEST_CASE("Self-thinning cascade prunes shaded interior leaves", "[ethylene][integration]") {
    Genome g = default_genome();
    // Make abscission fast for test
    g.senescence_duration = 5;
    g.ethylene_shade_rate = 1.0f;      // strong shade response
    g.ethylene_shade_threshold = 0.5f;
    g.ethylene_abscission_threshold = 0.3f; // easy to trigger

    Plant plant(g, glm::vec3(0.0f));

    // Build a small branch with shaded interior leaves and sunlit outer leaves.
    // seed -> stem -> inner_leaf (shaded) + outer_leaf (sunlit)
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);
    plant.seed_mut()->sugar = 100.0f; // plenty of sugar

    Node* inner_leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    inner_leaf->leaf_size = 0.2f;
    inner_leaf->light_exposure = 0.1f; // heavily shaded
    inner_leaf->position = glm::vec3(0.1f, 1.1f, 0.0f);
    stem->add_child(inner_leaf);

    Node* outer_leaf = plant.create_node(NodeType::LEAF, glm::vec3(-0.1f, 0.1f, 0.0f), 0.0f);
    outer_leaf->leaf_size = 0.2f;
    outer_leaf->light_exposure = 0.9f; // well-lit
    outer_leaf->position = glm::vec3(-0.1f, 1.1f, 0.0f);
    stem->add_child(outer_leaf);

    uint32_t inner_id = inner_leaf->id;
    uint32_t outer_id = outer_leaf->id;

    WorldParams wp = default_world_params();

    // Run ethylene + abscission for enough ticks to trigger senescence and removal
    for (int i = 0; i < 20; i++) {
        // Keep light_exposure fixed (simulate persistent shade)
        plant.for_each_node_mut([&](Node& n) {
            if (n.id == inner_id) n.light_exposure = 0.1f;
            if (n.id == outer_id) n.light_exposure = 0.9f;
        });
        compute_ethylene(plant, wp);
        process_abscission(plant);
    }

    // Inner (shaded) leaf should have been removed
    bool inner_exists = false;
    bool outer_exists = false;
    plant.for_each_node([&](const Node& n) {
        if (n.id == inner_id) inner_exists = true;
        if (n.id == outer_id) outer_exists = true;
    });

    REQUIRE_FALSE(inner_exists); // shaded leaf pruned
    REQUIRE(outer_exists);       // sunlit leaf survived
}
```

- [ ] **Step 2: Run test**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[ethylene][integration]" -v`
Expected: Self-thinning cascade test passes

- [ ] **Step 3: Commit**

```bash
git add tests/test_ethylene.cpp
git commit -m "test: self-thinning cascade integration test"
```
