# Sugar Milestone Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a sugar/photosynthesis resource system where leaves produce sugar from light, sugar diffuses bidirectionally through the plant, and all nodes consume maintenance sugar proportional to their mass/type.

**Architecture:** Leaves become real LEAF graph nodes (new NodeType). Sugar is a persistent per-node float (NOT reset each tick like auxin/cytokinin). Production happens at LEAF nodes, gradient-based diffusion runs multiple iterations per tick through all node connections, and every node consumes maintenance sugar proportional to its type and size. A WorldParams struct holds non-genetic simulation parameters (light level, diffusion iterations).

**Tech Stack:** C++17, GLM, Catch2 (tests), OpenGL 4.1 + GLFW + ImGui (renderer)

---

## File Structure

**New files:**
- `src/engine/world_params.h` — WorldParams struct (light level, diffusion iterations)
- `src/engine/sugar.h` — Sugar transport function declarations
- `src/engine/sugar.cpp` — Sugar production, diffusion, consumption logic
- `tests/test_sugar.cpp` — Tests for the sugar system

**Modified files:**
- `src/engine/node.h` — Add LEAF to NodeType, add `float sugar` + `float leaf_size`, remove Leaf struct + `Leaf* leaf`
- `src/engine/plant.h` — Remove `create_leaf()` and `leaves_` vector
- `src/engine/plant.cpp` — Remove `create_leaf()` implementation
- `src/engine/genome.h` — Add sugar genome parameters to Genome struct and default_genome()
- `src/engine/meristem.cpp` — Chain growth creates LEAF child nodes, skip thickening for LEAF
- `src/engine/engine.h` — Add `WorldParams` member and accessors
- `src/engine/engine.cpp` — Call `transport_sugar()` in tick loop
- `src/renderer/renderer.cpp` — Handle LEAF NodeType in draw_plant, draw LEAF as quad not cylinder
- `src/serialization/serializer.h` — Add `sugar` and `leaf_size` to NodeSnapshot, remove `has_leaf`
- `src/serialization/serializer.cpp` — Update binary read/write for new format
- `src/app_realtime.cpp` — Add `--color sugar`, add ImGui sugar stats + sliders
- `src/app_dump.cpp` — Add LEAF to node_type_str, print sugar values
- `CMakeLists.txt` — Add sugar.cpp to botany_engine, test_sugar.cpp to botany_tests
- `tests/test_node.cpp` — Update leaf tests for LEAF NodeType
- `tests/test_plant.cpp` — Remove create_leaf tests, add LEAF node tests
- `tests/test_meristem.cpp` — Update chain growth tests for LEAF child nodes
- `tests/test_genome.cpp` — Add sugar parameter validation
- `tests/test_serializer.cpp` — Update round-trip tests for new format

---

## Sprint 1: Sugar Foundation

### Task 1: Add LEAF to NodeType and sugar/leaf_size fields to Node

**Files:**
- Modify: `src/engine/node.h`

- [ ] **Step 1: Add LEAF to NodeType enum and new fields to Node**

In `src/engine/node.h`, add LEAF to the enum and add `sugar` + `leaf_size` to Node:

```cpp
enum class NodeType { STEM, ROOT, LEAF };
```

In the Node struct, after `float cytokinin;` add:

```cpp
float sugar = 0.0f;
float leaf_size = 0.0f;
```

Remove the `Leaf` struct entirely (lines 30-32):
```cpp
// DELETE:
// struct Leaf {
//     float size;
// };
```

Remove `Leaf* leaf = nullptr;` from the Node struct.

- [ ] **Step 2: Build to check for compile errors**

Run: `/usr/local/bin/cmake --build build 2>&1 | grep -i error | head -20`

Expected: Many compile errors from files still referencing `Leaf` and `node.leaf`. This confirms the field was removed and we need to update all references.

- [ ] **Step 3: Commit the node.h changes**

```bash
git add src/engine/node.h
git commit -m "refactor: add LEAF NodeType, sugar/leaf_size fields, remove Leaf struct"
```

---

### Task 2: Remove create_leaf from Plant

**Files:**
- Modify: `src/engine/plant.h`
- Modify: `src/engine/plant.cpp`

- [ ] **Step 1: Remove create_leaf and leaves_ from Plant**

In `src/engine/plant.h`, remove:
```cpp
// DELETE: Leaf* create_leaf(float size);
// DELETE from private: std::vector<std::unique_ptr<Leaf>> leaves_;
```

In `src/engine/plant.cpp`, remove the entire `create_leaf` function:
```cpp
// DELETE:
// Leaf* Plant::create_leaf(float size) {
//     auto l = std::make_unique<Leaf>(Leaf{size});
//     Leaf* ptr = l.get();
//     leaves_.push_back(std::move(l));
//     return ptr;
// }
```

- [ ] **Step 2: Commit**

```bash
git add src/engine/plant.h src/engine/plant.cpp
git commit -m "refactor: remove create_leaf and leaves_ from Plant"
```

---

### Task 3: Update chain growth to create LEAF nodes

**Files:**
- Modify: `src/engine/meristem.cpp`

- [ ] **Step 1: Update ShootApicalMeristem::tick chain growth**

In `ShootApicalMeristem::tick()`, replace the chain growth section (the block inside `if (dist > g.max_internode_length)`). The old code creates an axillary node with a leaf property. The new code creates the axillary node WITHOUT a leaf, and creates a separate LEAF child node on the interior STEM node:

```cpp
if (dist > g.max_internode_length) {
    // This node becomes interior — give it an axillary meristem
    glm::vec3 branch_dir = branch_direction(dir, g.branch_angle, node.id);
    glm::vec3 ax_pos = node.position + branch_dir * g.tip_offset;
    Node* axillary = plant.create_node(NodeType::STEM, ax_pos, g.initial_radius * 0.5f);
    axillary->meristem = plant.create_meristem<ShootAxillaryMeristem>();
    node.add_child(axillary);

    // Create leaf as a separate LEAF node on this interior node
    glm::vec3 leaf_dir = branch_direction(dir, g.branch_angle, node.id + 1000);
    glm::vec3 leaf_pos = node.position + leaf_dir * g.tip_offset;
    Node* leaf_node = plant.create_node(NodeType::LEAF, leaf_pos, 0.0f);
    leaf_node->leaf_size = g.leaf_size;
    node.add_child(leaf_node);

    // Create new tip node slightly ahead along current direction
    Node* new_tip = plant.create_node(NodeType::STEM, node.position + dir * g.tip_offset, node.radius);
    new_tip->meristem = node.meristem;
    node.meristem = nullptr;
    node.add_child(new_tip);
}
```

- [ ] **Step 2: Skip thickening for LEAF nodes in tick_meristems**

In `tick_meristems()`, change the thickening guard from:
```cpp
if (!is_active_tip) {
    n.radius += g.thickening_rate;
}
```
to:
```cpp
if (!is_active_tip && n.type != NodeType::LEAF) {
    n.radius += g.thickening_rate;
}
```

- [ ] **Step 3: Commit**

```bash
git add src/engine/meristem.cpp
git commit -m "refactor: chain growth creates LEAF child nodes instead of Leaf property"
```

---

### Task 4: Update renderer to handle LEAF nodes

**Files:**
- Modify: `src/renderer/renderer.cpp`

- [ ] **Step 1: Update draw_plant to handle LEAF NodeType**

In `draw_plant()`, restructure the `for_each_node` lambda. LEAF nodes should draw as leaf quads, not cylinders. Replace the entire lambda body:

```cpp
plant.for_each_node([&](const Node& node) {
    if (!node.parent) return;

    // LEAF nodes draw as quads, not cylinders
    if (node.type == NodeType::LEAF) {
        glm::vec3 dir = glm::normalize(node.position - node.parent->position);
        draw_leaf(node.position, dir, node.leaf_size);
        return;
    }

    glm::vec3 color;
    if (chemical_accessor_) {
        float v = chemical_accessor_(node);
        color = heatmap(v / max_val);
    } else if (color_by_type_) {
        if (node.type == NodeType::STEM) {
            color = glm::vec3(0.2f, 0.8f, 0.2f);
        } else {
            color = glm::vec3(0.8f, 0.4f, 0.1f);
        }
    } else if (node.type == NodeType::STEM) {
        color = glm::vec3(0.45f, 0.3f, 0.15f);
    } else {
        color = glm::vec3(0.35f, 0.2f, 0.1f);
    }

    draw_cylinder(node.parent->position, node.position,
                  node.parent->radius, node.radius, color);
});
```

- [ ] **Step 2: Update draw_snapshot for LEAF nodes**

In `draw_snapshot()`, add a LEAF check before the cylinder drawing:

```cpp
for (const auto& ns : snapshot.nodes) {
    if (ns.parent_id == UINT32_MAX) continue;

    auto it = id_to_idx.find(ns.parent_id);
    if (it == id_to_idx.end()) continue;
    const auto& parent = snapshot.nodes[it->second];

    if (ns.type == NodeType::LEAF) {
        glm::vec3 dir = glm::normalize(ns.position - parent.position);
        draw_leaf(ns.position, dir, ns.leaf_size);
        continue;
    }

    glm::vec3 color;
    if (ns.type == NodeType::STEM) {
        color = glm::vec3(0.45f, 0.3f, 0.15f);
    } else {
        color = glm::vec3(0.35f, 0.2f, 0.1f);
    }

    draw_cylinder(parent.position, ns.position,
                  parent.radius, ns.radius, color);
}
```

- [ ] **Step 3: Commit**

```bash
git add src/renderer/renderer.cpp
git commit -m "refactor: renderer draws LEAF nodes as quads, skips cylinders for leaves"
```

---

### Task 5: Update serializer for new node format

**Files:**
- Modify: `src/serialization/serializer.h`
- Modify: `src/serialization/serializer.cpp`

- [ ] **Step 1: Update NodeSnapshot struct**

In `serializer.h`, replace `bool has_leaf` with `float leaf_size` and add `float sugar`:

```cpp
struct NodeSnapshot {
    uint32_t id;
    uint32_t parent_id;
    NodeType type;
    glm::vec3 position;
    float radius;
    float auxin;
    float cytokinin;
    float sugar;
    float leaf_size;
};
```

- [ ] **Step 2: Update save_tick in serializer.cpp**

Replace the `has_leaf` / `leaf_size` writes with:

```cpp
write_val(out, node.sugar);
write_val(out, node.leaf_size);
```

Remove the old `has_leaf` bool write and the conditional leaf_size write.

- [ ] **Step 3: Update load_tick in serializer.cpp**

Replace the `has_leaf` / `leaf_size` reads with:

```cpp
ns.sugar = read_val<float>(in);
ns.leaf_size = read_val<float>(in);
```

- [ ] **Step 4: Commit**

```bash
git add src/serialization/serializer.h src/serialization/serializer.cpp
git commit -m "refactor: update serializer for LEAF nodes, add sugar field to snapshots"
```

---

### Task 6: Update app_dump for LEAF nodes

**Files:**
- Modify: `src/app_dump.cpp`

- [ ] **Step 1: Update node_type_str**

```cpp
static const char* node_type_str(NodeType t) {
    switch (t) {
        case NodeType::STEM: return "STEM";
        case NodeType::ROOT: return "ROOT";
        case NodeType::LEAF: return "LEAF";
    }
    return "???";
}
```

- [ ] **Step 2: Update print functions**

In `print_tick_full()`, replace `has_leaf` checks with `type == NodeType::LEAF` checks. Add sugar to the output line alongside auxin and cytokinin.

In `print_tick_stats()`, count LEAF nodes via `type == NodeType::LEAF` instead of `has_leaf`. Add sugar min/max/avg.

- [ ] **Step 3: Commit**

```bash
git add src/app_dump.cpp
git commit -m "refactor: app_dump handles LEAF NodeType and prints sugar values"
```

---

### Task 7: Update all tests for LEAF node conversion

**Files:**
- Modify: `tests/test_node.cpp`
- Modify: `tests/test_plant.cpp`
- Modify: `tests/test_meristem.cpp`
- Modify: `tests/test_serializer.cpp`

- [ ] **Step 1: Update test_node.cpp**

Replace the "Node can have a leaf attached" test with:

```cpp
TEST_CASE("LEAF node stores leaf_size", "[node]") {
    Node node(1, NodeType::LEAF, glm::vec3(0.0f), 0.0f);
    node.leaf_size = 0.3f;
    REQUIRE(node.type == NodeType::LEAF);
    REQUIRE(node.leaf_size == 0.3f);
    REQUIRE(node.sugar == 0.0f);
}
```

- [ ] **Step 2: Update test_plant.cpp**

Remove the `create_leaf` test. Replace with:

```cpp
TEST_CASE("Plant can create LEAF nodes", "[plant]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 1.0f, 0.0f), 0.0f);
    leaf->leaf_size = 0.3f;
    REQUIRE(leaf->type == NodeType::LEAF);
    REQUIRE(leaf->leaf_size == 0.3f);
}
```

- [ ] **Step 3: Update test_meristem.cpp**

In the chain growth test that checks for leaves, replace `n.leaf != nullptr` with checking for a LEAF child node:

```cpp
// Instead of checking n.leaf, check for LEAF child
bool has_leaf_sibling = false;
if (n.parent) {
    for (const Node* sibling : n.parent->children) {
        if (sibling->type == NodeType::LEAF) {
            has_leaf_sibling = true;
        }
    }
}
```

Update the "interior nodes have at most N children" assertion — interior STEM nodes now have up to 3 children (continuation tip, axillary, and LEAF).

- [ ] **Step 4: Update test_serializer.cpp**

Add `sugar` assertions to the round-trip test:

```cpp
REQUIRE(snap.nodes[0].sugar == 0.0f);
REQUIRE(snap.nodes[0].leaf_size == 0.0f);
```

- [ ] **Step 5: Build and run all tests**

Run: `/usr/local/bin/cmake --build build 2>&1 | tail -5`
Run: `cd /Users/wldarden/learning/botany/build && ./botany_tests 2>&1 | tail -5`

Expected: All tests pass. Zero compile errors.

- [ ] **Step 6: Commit**

```bash
git add tests/test_node.cpp tests/test_plant.cpp tests/test_meristem.cpp tests/test_serializer.cpp
git commit -m "test: update all tests for LEAF node conversion"
```

---

### Task 8: Add WorldParams and sugar genome parameters

**Files:**
- Create: `src/engine/world_params.h`
- Modify: `src/engine/genome.h`
- Modify: `src/engine/engine.h`
- Modify: `src/engine/engine.cpp`

- [ ] **Step 1: Create world_params.h**

```cpp
// src/engine/world_params.h
#pragma once

namespace botany {

struct WorldParams {
    float light_level = 1.0f;
    int sugar_diffusion_iterations = 5;
};

inline WorldParams default_world_params() {
    return WorldParams{
        .light_level = 1.0f,
        .sugar_diffusion_iterations = 5,
    };
}

} // namespace botany
```

- [ ] **Step 2: Add sugar parameters to Genome**

In `src/engine/genome.h`, add after the Geometry section:

```cpp
// Sugar / photosynthesis
float sugar_production_rate;
float sugar_transport_conductance;
float sugar_maintenance_leaf;
float sugar_maintenance_stem;
float sugar_maintenance_root;
float sugar_maintenance_meristem;
```

In `default_genome()`, add:

```cpp
.sugar_production_rate = 0.5f,
.sugar_transport_conductance = 0.1f,
.sugar_maintenance_leaf = 0.02f,
.sugar_maintenance_stem = 0.01f,
.sugar_maintenance_root = 0.01f,
.sugar_maintenance_meristem = 0.005f,
```

- [ ] **Step 3: Add WorldParams to Engine**

In `src/engine/engine.h`, add `#include "engine/world_params.h"` and:

```cpp
const WorldParams& world_params() const { return world_params_; }
WorldParams& world_params_mut() { return world_params_; }
```

Add private member:
```cpp
WorldParams world_params_;
```

- [ ] **Step 4: Update test_genome.cpp**

Add a section for sugar parameters:

```cpp
SECTION("sugar parameters are positive") {
    REQUIRE(g.sugar_production_rate > 0.0f);
    REQUIRE(g.sugar_transport_conductance > 0.0f);
    REQUIRE(g.sugar_maintenance_leaf > 0.0f);
    REQUIRE(g.sugar_maintenance_stem > 0.0f);
    REQUIRE(g.sugar_maintenance_root > 0.0f);
    REQUIRE(g.sugar_maintenance_meristem > 0.0f);
}
```

- [ ] **Step 5: Build and test**

Run: `/usr/local/bin/cmake --build build 2>&1 | tail -5`
Run: `cd /Users/wldarden/learning/botany/build && ./botany_tests 2>&1 | tail -5`

Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/engine/world_params.h src/engine/genome.h src/engine/engine.h src/engine/engine.cpp tests/test_genome.cpp
git commit -m "feat: add WorldParams and sugar genome parameters"
```

---

### Task 9: Write failing tests for sugar production

**Files:**
- Create: `tests/test_sugar.cpp`
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Create test_sugar.cpp with production tests**

```cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/sugar.h"

using namespace botany;
using Catch::Matchers::WithinAbs;

TEST_CASE("LEAF nodes produce sugar proportional to light and leaf_size", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Manually create a LEAF node attached to the seed
    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.5f;
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 2.0f;

    produce_sugar(plant, wp);

    float expected = wp.light_level * leaf->leaf_size * g.sugar_production_rate;
    REQUIRE_THAT(leaf->sugar, WithinAbs(expected, 1e-6));
}

TEST_CASE("Zero light produces zero sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.5f;
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 0.0f;

    produce_sugar(plant, wp);

    REQUIRE(leaf->sugar == 0.0f);
}

TEST_CASE("Non-LEAF nodes do not produce sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    WorldParams wp = default_world_params();
    produce_sugar(plant, wp);

    plant.for_each_node([](const Node& n) {
        if (n.type != NodeType::LEAF) {
            REQUIRE(n.sugar == 0.0f);
        }
    });
}
```

- [ ] **Step 2: Add test_sugar.cpp and sugar.cpp to CMakeLists.txt**

In the `botany_engine` sources, add `src/engine/sugar.cpp`.

In the `botany_tests` sources, add `tests/test_sugar.cpp`.

- [ ] **Step 3: Create sugar.h header**

```cpp
// src/engine/sugar.h
#pragma once

namespace botany {

class Plant;
struct WorldParams;

void produce_sugar(Plant& plant, const WorldParams& world);
void consume_sugar(Plant& plant);
void diffuse_sugar(Plant& plant, const WorldParams& world);
void transport_sugar(Plant& plant, const WorldParams& world);

} // namespace botany
```

- [ ] **Step 4: Create sugar.cpp with empty stubs**

```cpp
// src/engine/sugar.cpp
#include "engine/sugar.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include "engine/node.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace botany {

void produce_sugar(Plant& plant, const WorldParams& world) {
    // TODO: implement
}

void consume_sugar(Plant& plant) {
    // TODO: implement
}

void diffuse_sugar(Plant& plant, const WorldParams& world) {
    // TODO: implement
}

void transport_sugar(Plant& plant, const WorldParams& world) {
    produce_sugar(plant, world);
    diffuse_sugar(plant, world);
    consume_sugar(plant);
}

} // namespace botany
```

- [ ] **Step 5: Build and run tests to see production tests fail**

Run: `/usr/local/bin/cmake --build build 2>&1 | tail -5`
Run: `cd /Users/wldarden/learning/botany/build && ./botany_tests "[sugar]" 2>&1`

Expected: Tests FAIL because `produce_sugar` is a stub.

- [ ] **Step 6: Commit failing tests**

```bash
git add src/engine/sugar.h src/engine/sugar.cpp tests/test_sugar.cpp CMakeLists.txt
git commit -m "test: add failing sugar production tests and sugar module stubs"
```

---

### Task 10: Implement sugar production

**Files:**
- Modify: `src/engine/sugar.cpp`

- [ ] **Step 1: Implement produce_sugar**

Replace the stub with:

```cpp
void produce_sugar(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();
    plant.for_each_node_mut([&](Node& node) {
        if (node.type == NodeType::LEAF) {
            node.sugar += world.light_level * node.leaf_size * g.sugar_production_rate;
        }
    });
}
```

- [ ] **Step 2: Run tests to verify they pass**

Run: `cd /Users/wldarden/learning/botany/build && /usr/local/bin/cmake --build . 2>&1 | tail -3 && ./botany_tests "[sugar]" -v 2>&1`

Expected: All 3 sugar tests pass.

- [ ] **Step 3: Commit**

```bash
git add src/engine/sugar.cpp
git commit -m "feat: implement sugar production at LEAF nodes"
```

---

### Task 11: Write failing tests for sugar consumption

**Files:**
- Modify: `tests/test_sugar.cpp`

- [ ] **Step 1: Add consumption tests**

Append to `test_sugar.cpp`:

```cpp
TEST_CASE("consume_sugar deducts maintenance cost by node type", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Set known sugar on seed (STEM)
    Node* seed = plant.seed_mut();
    seed->sugar = 10.0f;

    consume_sugar(plant);

    float expected_cost = g.sugar_maintenance_stem * seed->radius;
    REQUIRE_THAT(seed->sugar, WithinAbs(10.0f - expected_cost, 1e-6));
}

TEST_CASE("LEAF maintenance cost uses leaf_size", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.5f;
    leaf->sugar = 10.0f;
    plant.seed_mut()->add_child(leaf);

    consume_sugar(plant);

    float expected_cost = g.sugar_maintenance_leaf * leaf->leaf_size;
    REQUIRE_THAT(leaf->sugar, WithinAbs(10.0f - expected_cost, 1e-6));
}

TEST_CASE("Sugar cannot go below zero after consumption", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    seed->sugar = 0.0001f; // Less than any maintenance cost

    consume_sugar(plant);

    REQUIRE(seed->sugar >= 0.0f);
}

TEST_CASE("Active meristem tips have additional maintenance cost", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Find the shoot tip (has active apical meristem)
    Node* shoot_tip = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.meristem && n.meristem->is_tip() && n.meristem->active &&
            n.type == NodeType::STEM) {
            shoot_tip = &n;
        }
    });
    REQUIRE(shoot_tip != nullptr);
    shoot_tip->sugar = 10.0f;

    consume_sugar(plant);

    float expected_cost = g.sugar_maintenance_stem * shoot_tip->radius
                        + g.sugar_maintenance_meristem;
    REQUIRE_THAT(shoot_tip->sugar, WithinAbs(10.0f - expected_cost, 1e-6));
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/wldarden/learning/botany/build && /usr/local/bin/cmake --build . 2>&1 | tail -3 && ./botany_tests "[sugar]" -v 2>&1`

Expected: New consumption tests fail (consume_sugar is a stub).

- [ ] **Step 3: Commit failing tests**

```bash
git add tests/test_sugar.cpp
git commit -m "test: add failing sugar consumption tests"
```

---

### Task 12: Implement sugar consumption

**Files:**
- Modify: `src/engine/sugar.cpp`

- [ ] **Step 1: Implement consume_sugar**

Replace the stub with:

```cpp
void consume_sugar(Plant& plant) {
    const Genome& g = plant.genome();
    plant.for_each_node_mut([&](Node& node) {
        float cost = 0.0f;
        switch (node.type) {
            case NodeType::LEAF:
                cost = g.sugar_maintenance_leaf * node.leaf_size;
                break;
            case NodeType::STEM:
                cost = g.sugar_maintenance_stem * node.radius;
                break;
            case NodeType::ROOT:
                cost = g.sugar_maintenance_root * node.radius;
                break;
        }
        if (node.meristem && node.meristem->is_tip() && node.meristem->active) {
            cost += g.sugar_maintenance_meristem;
        }
        node.sugar = std::max(0.0f, node.sugar - cost);
    });
}
```

- [ ] **Step 2: Run tests**

Run: `cd /Users/wldarden/learning/botany/build && /usr/local/bin/cmake --build . 2>&1 | tail -3 && ./botany_tests "[sugar]" -v 2>&1`

Expected: All sugar tests pass.

- [ ] **Step 3: Commit**

```bash
git add src/engine/sugar.cpp
git commit -m "feat: implement sugar consumption with per-type maintenance costs"
```

---

### Task 13: Write failing tests for sugar diffusion

**Files:**
- Modify: `tests/test_sugar.cpp`

- [ ] **Step 1: Add diffusion tests**

Append to `test_sugar.cpp`:

```cpp
TEST_CASE("Sugar diffuses from high to low concentration", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Give seed node a lot of sugar, children start at 0
    Node* seed = plant.seed_mut();
    seed->sugar = 100.0f;

    WorldParams wp = default_world_params();
    wp.sugar_diffusion_iterations = 1;

    float seed_before = seed->sugar;
    diffuse_sugar(plant, wp);

    // Seed should have lost some sugar
    REQUIRE(seed->sugar < seed_before);
    // Children should have gained some sugar
    for (const Node* child : seed->children) {
        REQUIRE(child->sugar > 0.0f);
    }
}

TEST_CASE("Sugar diffusion preserves total sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Set known sugar on seed only
    plant.seed_mut()->sugar = 100.0f;

    float total_before = 0.0f;
    plant.for_each_node([&](const Node& n) { total_before += n.sugar; });

    WorldParams wp = default_world_params();
    wp.sugar_diffusion_iterations = 5;
    diffuse_sugar(plant, wp);

    float total_after = 0.0f;
    plant.for_each_node([&](const Node& n) { total_after += n.sugar; });

    REQUIRE_THAT(total_after, WithinAbs(total_before, 1e-4));
}

TEST_CASE("Multiple diffusion iterations produce smoother distribution", "[sugar]") {
    Genome g = default_genome();

    // Run with 1 iteration
    Plant plant1(g, glm::vec3(0.0f));
    plant1.seed_mut()->sugar = 100.0f;
    WorldParams wp1 = default_world_params();
    wp1.sugar_diffusion_iterations = 1;
    diffuse_sugar(plant1, wp1);

    // Run with 10 iterations
    Plant plant2(g, glm::vec3(0.0f));
    plant2.seed_mut()->sugar = 100.0f;
    WorldParams wp2 = default_world_params();
    wp2.sugar_diffusion_iterations = 10;
    diffuse_sugar(plant2, wp2);

    // With more iterations, children should have more sugar (more spreading)
    float child_sugar_1iter = 0.0f;
    float child_sugar_10iter = 0.0f;
    for (const Node* c : plant1.seed()->children) { child_sugar_1iter += c->sugar; }
    for (const Node* c : plant2.seed()->children) { child_sugar_10iter += c->sugar; }

    REQUIRE(child_sugar_10iter > child_sugar_1iter);
}
```

- [ ] **Step 2: Run to see failures**

Run: `cd /Users/wldarden/learning/botany/build && /usr/local/bin/cmake --build . 2>&1 | tail -3 && ./botany_tests "[sugar]" -v 2>&1`

Expected: Diffusion tests fail (diffuse_sugar is a stub).

- [ ] **Step 3: Commit failing tests**

```bash
git add tests/test_sugar.cpp
git commit -m "test: add failing sugar diffusion tests"
```

---

### Task 14: Implement sugar diffusion

**Files:**
- Modify: `src/engine/sugar.cpp`

- [ ] **Step 1: Implement diffuse_sugar**

Replace the stub with the gradient-based bidirectional diffusion:

```cpp
void diffuse_sugar(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    // Build edge list once — all parent-child connections
    struct Edge {
        Node* a;
        Node* b;
        float capacity;
    };
    std::vector<Edge> edges;
    plant.for_each_node_mut([&](Node& node) {
        if (node.parent) {
            float min_radius = std::min(node.radius, node.parent->radius);
            // LEAF nodes have radius 0, use a small baseline for leaf connections
            if (node.type == NodeType::LEAF || node.parent->type == NodeType::LEAF) {
                min_radius = std::max(min_radius, 0.01f);
            }
            float capacity = min_radius * min_radius * 3.14159f * g.sugar_transport_conductance;
            edges.push_back({&node, node.parent, capacity});
        }
    });

    // Run multiple diffusion iterations
    for (int iter = 0; iter < world.sugar_diffusion_iterations; iter++) {
        // Compute all flows first, then apply (avoids order-dependent artifacts)
        struct Flow {
            Node* from;
            Node* to;
            float amount;
        };
        std::vector<Flow> flows;

        for (const Edge& e : edges) {
            float gradient = e.a->sugar - e.b->sugar;
            if (std::abs(gradient) < 1e-6f) continue;

            float flow = gradient * e.capacity;
            if (flow > 0.0f) {
                flow = std::min(flow, e.a->sugar);
                flows.push_back({e.a, e.b, flow});
            } else {
                flow = std::min(-flow, e.b->sugar);
                flows.push_back({e.b, e.a, flow});
            }
        }

        for (const Flow& f : flows) {
            f.from->sugar -= f.amount;
            f.to->sugar += f.amount;
        }
    }
}
```

- [ ] **Step 2: Run tests**

Run: `cd /Users/wldarden/learning/botany/build && /usr/local/bin/cmake --build . 2>&1 | tail -3 && ./botany_tests "[sugar]" -v 2>&1`

Expected: All sugar tests pass.

- [ ] **Step 3: Commit**

```bash
git add src/engine/sugar.cpp
git commit -m "feat: implement gradient-based bidirectional sugar diffusion"
```

---

### Task 15: Wire sugar transport into Engine tick loop

**Files:**
- Modify: `src/engine/engine.cpp`

- [ ] **Step 1: Add transport_sugar call**

Add `#include "engine/sugar.h"` at the top.

In `Engine::tick()`, add `transport_sugar` after cytokinin and before meristems:

```cpp
void Engine::tick() {
    for (auto& plant : plants_) {
        transport_auxin(*plant);
        transport_cytokinin(*plant);
        transport_sugar(*plant, world_params_);
        tick_meristems(*plant);
    }
    tick_++;
}
```

- [ ] **Step 2: Write integration test**

Add to `test_sugar.cpp`:

```cpp
TEST_CASE("transport_sugar runs full produce-diffuse-consume cycle", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Create a LEAF with known leaf_size
    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 1.0f;
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 1.0f;

    transport_sugar(plant, wp);

    // Leaf should have produced sugar, some diffused away, some consumed
    // Net result: leaf has some sugar > 0
    REQUIRE(leaf->sugar > 0.0f);

    // Adjacent seed node should have received some sugar via diffusion
    REQUIRE(plant.seed()->sugar > 0.0f);
}
```

- [ ] **Step 3: Build and run all tests**

Run: `cd /Users/wldarden/learning/botany/build && /usr/local/bin/cmake --build . 2>&1 | tail -5 && ./botany_tests 2>&1 | tail -5`

Expected: ALL tests pass (sugar tests + existing hormone/meristem/engine tests).

- [ ] **Step 4: Commit**

```bash
git add src/engine/engine.cpp tests/test_sugar.cpp
git commit -m "feat: wire transport_sugar into engine tick loop"
```

---

### Task 16: Add --color sugar visualization and ImGui stats

**Files:**
- Modify: `src/app_realtime.cpp`

- [ ] **Step 1: Add sugar to --color argument**

In the `--color` parsing block, add a sugar case:

```cpp
} else if (color_chemical == "sugar") {
    accessor = [](const Node& n) { return n.sugar; };
```

Update the error message to include sugar:
```cpp
std::cerr << "Unknown chemical: " << color_chemical
          << " (available: auxin, cytokinin, sugar, type)" << std::endl;
```

- [ ] **Step 2: Add sugar stats and controls to ImGui panel**

After the existing `ImGui::Text("Tick: %d  Nodes: %d", ...)` line, add:

```cpp
// Sugar stats
float total_sugar = 0.0f;
float max_sugar = 0.0f;
int leaf_count = 0;
engine.get_plant(plant_id).for_each_node([&](const Node& n) {
    total_sugar += n.sugar;
    if (n.sugar > max_sugar) max_sugar = n.sugar;
    if (n.type == NodeType::LEAF) leaf_count++;
});
ImGui::Text("Leaves: %d", leaf_count);
ImGui::Text("Sugar: total=%.1f max=%.3f", total_sugar, max_sugar);
ImGui::Separator();
ImGui::SliderFloat("Light", &engine.world_params_mut().light_level, 0.0f, 2.0f);
ImGui::SliderInt("Diffusion Iters",
    &engine.world_params_mut().sugar_diffusion_iterations, 1, 20);
```

- [ ] **Step 3: Build**

Run: `/usr/local/bin/cmake --build build 2>&1 | tail -5`

Expected: Compiles cleanly.

- [ ] **Step 4: Commit**

```bash
git add src/app_realtime.cpp
git commit -m "feat: add --color sugar visualization and ImGui sugar controls"
```

---

### Task 17: Final verification and docs

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Run all tests**

Run: `cd /Users/wldarden/learning/botany/build && /usr/local/bin/cmake --build . 2>&1 | tail -5 && ./botany_tests 2>&1 | tail -5`

Expected: All tests pass.

- [ ] **Step 2: Visual verification**

Run: `cd /Users/wldarden/learning/botany && ./build/botany_realtime --color sugar`

Expected:
- Leaves show bright colors (high sugar)
- Sugar fades through stem nodes away from leaves
- Roots show low sugar
- ImGui panel shows leaf count, sugar stats, light slider, diffusion slider
- Setting light to 0 via slider causes sugar to drain over time

Run: `cd /Users/wldarden/learning/botany && ./build/botany_realtime --color type`

Expected: No visual regression. Green stems, orange roots, leaf quads visible.

- [ ] **Step 3: Update CLAUDE.md**

Add a Sugar Model section documenting:
- Sugar persists across ticks (not reset like hormones)
- Produced by LEAF nodes: `light_level * leaf_size * sugar_production_rate`
- Consumed by all nodes: maintenance cost proportional to type and size
- Diffused bidirectionally via gradient, multiple iterations per tick
- Transport capacity proportional to cross-sectional area (thicker = more)
- WorldParams holds light_level and sugar_diffusion_iterations

- [ ] **Step 4: Commit everything**

```bash
git add CLAUDE.md
git commit -m "docs: document sugar model in CLAUDE.md"
```

---

## Sprint 2: Activity-Based Sugar Consumption (Outline)

**Goal:** Growing costs sugar. No sugar = no growth/activation.

**New genome parameters:**
- `sugar_cost_activation` — cost for axillary meristem to activate
- `sugar_cost_growth` — cost per unit of extension growth
- `sugar_cost_thickening` — cost per unit of radius increase

**Changes to meristem.cpp:**
- `ShootAxillaryMeristem::tick()`: check `node.sugar >= sugar_cost_activation` before activating, deduct on activation
- `ShootApicalMeristem::tick()`: check sugar before extending position, deduct proportional to distance moved
- `tick_meristems()` thickening: check sugar before thickening, deduct proportional to radius increase
- Same patterns for root meristems

---

## Sprint 3: Sugar Modifies Growth Rates (Outline)

**Goal:** Growth rate scales continuously with sugar availability.

**New genome parameters (per node type):**
- `sugar_save_threshold` — sugar reserve that shouldn't be used for growth
- `sugar_max_growth_cost` — sugar amount for maximum growth rate

**Growth formula:**
```
available = max(local_sugar - save_threshold, 0)
growth_fraction = min(available / max_growth_cost, 1.0)
actual_growth = max_growth_rate * growth_fraction
sugar_cost = actual_growth * cost_per_unit
```

---

## Sprint 4: Starvation and Node Death (Outline)

**Goal:** Nodes starving for sugar die and fall off.

**New fields:** `uint32_t starvation_ticks` on Node
**New WorldParams:** `uint32_t starvation_ticks_max` (default ~50)

**Behavior:**
- Each tick with `sugar == 0`: increment `starvation_ticks`
- Each tick with `sugar > 0`: reset `starvation_ticks = 0`
- Visual: interpolate node color toward brown/grey based on `starvation_ticks / starvation_ticks_max`
- At max starvation: remove node and all descendants from tree

**Critical complexity:** Node removal requires `Plant::remove_subtree()` — first deletion in the codebase. Needs careful pointer cleanup (parent's children vector, meristem ownership, etc.).
