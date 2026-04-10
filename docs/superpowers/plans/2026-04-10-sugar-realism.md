# Sugar Realism Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make sugar economics realistic — volume-based maintenance (~50% of production) and per-node storage caps so a plant survives 2-4 weeks in darkness, not a year.

**Architecture:** Two changes to the sugar system: (1) maintenance cost formula changes from `rate * radius` to `rate * volume` where volume = π·r²·length, making thick old wood expensive; (2) each node has a sugar storage cap proportional to its tissue volume, enforced in production, diffusion, and a safety clamp.

**Tech Stack:** C++17, Catch2 tests, glm vectors

**Spec:** `docs/superpowers/specs/2026-04-10-sugar-realism-design.md`

---

### Task 1: Add genome params and sugar_cap() helper

**Files:**
- Modify: `src/engine/genome.h:52-118` — update maintenance comments/values, add 3 storage params
- Modify: `src/engine/sugar.h:1-16` — declare `sugar_cap()`
- Modify: `src/engine/sugar.cpp:1-8` — add `#include <glm/geometric.hpp>`, implement `sugar_cap()`
- Test: `tests/test_sugar.cpp`

- [ ] **Step 1: Write failing tests for sugar_cap**

Add these tests at the end of `tests/test_sugar.cpp`:

```cpp
// === Sugar cap tests ===

TEST_CASE("sugar_cap scales with stem volume", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Create a stem with known dimensions: radius 0.1, internode length 1.0
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.1f);
    plant.seed_mut()->add_child(stem);

    float volume = 3.14159f * 0.1f * 0.1f * 1.0f;  // π * r² * length
    float expected = volume * g.sugar_storage_density_wood;
    REQUIRE_THAT(sugar_cap(*stem, g), WithinAbs(expected, 1e-4));
}

TEST_CASE("sugar_cap scales with leaf area", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.3f;
    plant.seed_mut()->add_child(leaf);

    float area = 0.3f * 0.3f;
    float expected = area * g.sugar_storage_density_leaf;
    REQUIRE_THAT(sugar_cap(*leaf, g), WithinAbs(expected, 1e-6));
}

TEST_CASE("sugar_cap returns minimum for tiny nodes", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Tiny node: radius 0.001, offset 0.001 — volume-based cap is below minimum
    Node* tiny = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.001f, 0.0f), 0.001f);
    plant.seed_mut()->add_child(tiny);

    REQUIRE_THAT(sugar_cap(*tiny, g), WithinAbs(g.sugar_cap_minimum, 1e-6));
}

TEST_CASE("sugar_cap for seed node covers seed_sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Seed node (no parent) should have cap >= seed_sugar
    REQUIRE(sugar_cap(*plant.seed(), g) >= g.seed_sugar);
}

TEST_CASE("sugar_cap for root scales with volume", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* root = plant.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.8f, 0.0f), 0.025f);
    plant.seed_mut()->add_child(root);

    float length = 0.8f;  // glm::length of offset
    float volume = 3.14159f * 0.025f * 0.025f * length;
    float expected = volume * g.sugar_storage_density_wood;
    // Volume cap (0.0785) > minimum (0.01), so volume wins
    REQUIRE_THAT(sugar_cap(*root, g), WithinAbs(expected, 1e-5));
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[sugar]" -c "sugar_cap"`

Expected: compile error — `sugar_cap` not declared, `sugar_storage_density_wood` not a member of Genome.

- [ ] **Step 3: Add genome params**

In `src/engine/genome.h`, update the sugar section comments and add new params. Replace the sugar params block (lines 53-66) with:

```cpp
    // Sugar / photosynthesis (g glucose)
    float sugar_production_rate;      // g glucose / (dm leaf_size · hr) at full sun
    float sugar_transport_conductance;// conductance scaling for diffusion between nodes
    float sugar_maintenance_leaf;     // g glucose / (dm² leaf area · hr)
    float sugar_maintenance_stem;     // g glucose / (dm³ volume · hr)
    float sugar_maintenance_root;     // g glucose / (dm³ volume · hr)
    float sugar_maintenance_meristem; // g glucose / hr per active meristem tip

    float seed_sugar;                 // g glucose — initial reserves in the seed

    // Sugar storage caps — maximum sugar a node can hold, proportional to tissue volume.
    // Prevents unbounded accumulation; represents finite starch storage capacity.
    float sugar_storage_density_wood; // g glucose max / dm³ of stem/root tissue
    float sugar_storage_density_leaf; // g glucose max / dm² of leaf area
    float sugar_cap_minimum;          // g glucose — floor for tiny/new nodes

    // Sugar save thresholds — minimum reserve before growth occurs (g glucose)
    float sugar_save_shoot;           // reserve for shoot apical meristems
    float sugar_save_root;            // reserve for root apical meristems
    float sugar_save_stem;            // reserve for stem thickening
```

Update the default values in `default_genome()` (replace lines 106-118):

```cpp
        .sugar_production_rate = 0.012f,    // g glucose / (dm leaf · hr) — ~7 g/m²/day
        .sugar_transport_conductance = 1.0f,  // ~5% gradient transfer per iter at 5mm radius
        .sugar_maintenance_leaf = 0.013f,   // g / (dm² · hr) — leaves dominate maintenance budget
        .sugar_maintenance_stem = 0.028f,   // g / (dm³ · hr) — wood is cheap per volume
        .sugar_maintenance_root = 0.135f,   // g / (dm³ · hr) — fine roots expensive (high turnover)
        .sugar_maintenance_meristem = 0.001f, // high per mass, small organ

        .seed_sugar = 8.0f,                 // ~15 days heterotrophic growth

        .sugar_storage_density_wood = 50.0f,  // g glucose max / dm³ — ~5% of dry mass as starch
        .sugar_storage_density_leaf = 0.5f,   // g glucose max / dm² — thin tissue, less storage
        .sugar_cap_minimum = 0.01f,           // floor for new/tiny nodes

        .sugar_save_shoot = 0.01f,          // buffer before growth
        .sugar_save_root = 0.005f,
        .sugar_save_stem = 0.02f,
```

- [ ] **Step 4: Declare sugar_cap in sugar.h**

In `src/engine/sugar.h`, add the declaration after the existing function declarations (before the closing `}`):

```cpp
float sugar_cap(const Node& node, const Genome& g);
```

Also add the forward declarations needed at the top of the file (after `namespace botany {`):

```cpp
struct Node;
struct Genome;
```

- [ ] **Step 5: Implement sugar_cap in sugar.cpp**

Add `#include <glm/geometric.hpp>` to the includes in `src/engine/sugar.cpp`.

Add the implementation after the namespace opening, before `produce_sugar`:

```cpp
float sugar_cap(const Node& node, const Genome& g) {
    float length = std::max(glm::length(node.offset), 0.01f);

    switch (node.type) {
        case NodeType::LEAF: {
            float area = node.leaf_size * node.leaf_size;
            return std::max(g.sugar_cap_minimum, area * g.sugar_storage_density_leaf);
        }
        case NodeType::STEM:
        case NodeType::ROOT: {
            float volume = 3.14159f * node.radius * node.radius * length;
            float cap = volume * g.sugar_storage_density_wood;
            // Seed node (no parent): cap must hold initial reserves
            if (!node.parent) {
                cap = std::max(cap, g.seed_sugar);
            }
            return std::max(g.sugar_cap_minimum, cap);
        }
    }
    return g.sugar_cap_minimum;
}
```

- [ ] **Step 6: Run tests to verify sugar_cap tests pass**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[sugar]" -c "sugar_cap"`

Expected: All 5 new sugar_cap tests PASS. Some existing consume tests may fail (genome value changes) — that's expected and fixed in Task 2.

- [ ] **Step 7: Commit**

```bash
git add src/engine/genome.h src/engine/sugar.h src/engine/sugar.cpp tests/test_sugar.cpp
git commit -m "feat: add sugar storage cap params and sugar_cap() helper"
```

---

### Task 2: Volume-based maintenance in consume_sugar

**Files:**
- Modify: `src/engine/sugar.cpp:59-86` — update `consume_sugar()` formula
- Modify: `tests/test_sugar.cpp` — fix broken tests, add volume-based test

- [ ] **Step 1: Update existing maintenance tests for volume-based formula**

In `tests/test_sugar.cpp`, add `#include <glm/geometric.hpp>` at the top.

Replace the `consume_sugar deducts maintenance cost by node type` test (the one checking stem maintenance on seed):

```cpp
TEST_CASE("consume_sugar deducts volume-based maintenance cost", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Create a stem with known dimensions for predictable cost
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.1f);
    stem->sugar = 10.0f;
    plant.seed_mut()->add_child(stem);

    // Zero out other nodes' sugar so we only check our stem
    plant.seed_mut()->sugar = 10.0f;

    consume_sugar(plant);

    float length = glm::length(stem->offset);  // 1.0
    float volume = 3.14159f * 0.1f * 0.1f * length;
    float expected_cost = g.sugar_maintenance_stem * volume;
    REQUIRE_THAT(stem->sugar, WithinAbs(10.0f - expected_cost, 1e-5));
}
```

Replace the `LEAF maintenance cost uses leaf_size` test:

```cpp
TEST_CASE("LEAF maintenance cost uses leaf area", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.5f;
    leaf->sugar = 10.0f;
    plant.seed_mut()->add_child(leaf);

    consume_sugar(plant);

    float expected_cost = g.sugar_maintenance_leaf * 0.5f * 0.5f;
    REQUIRE_THAT(leaf->sugar, WithinAbs(10.0f - expected_cost, 1e-6));
}
```

Replace the `Active meristem tips have additional maintenance cost` test:

```cpp
TEST_CASE("Active meristem tips have additional maintenance cost", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

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

    float length = std::max(glm::length(shoot_tip->offset), 0.01f);
    float volume = 3.14159f * shoot_tip->radius * shoot_tip->radius * length;
    float expected_cost = g.sugar_maintenance_stem * volume
                        + g.sugar_maintenance_meristem;
    REQUIRE_THAT(shoot_tip->sugar, WithinAbs(10.0f - expected_cost, 1e-5));
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[sugar]" -c "maintenance"`

Expected: FAIL — tests expect volume-based costs but `consume_sugar` still uses linear formula.

- [ ] **Step 3: Update consume_sugar to volume-based formula**

In `src/engine/sugar.cpp`, add `#include <glm/geometric.hpp>` if not already added in Task 1.

Replace the `consume_sugar` function body (lines 59-86):

```cpp
void consume_sugar(Plant& plant) {
    const Genome& g = plant.genome();
    plant.for_each_node_mut([&](Node& node) {
        float cost = 0.0f;
        float length = std::max(glm::length(node.offset), 0.01f);

        switch (node.type) {
            case NodeType::LEAF:
                cost = g.sugar_maintenance_leaf * node.leaf_size * node.leaf_size;
                break;
            case NodeType::STEM: {
                float volume = 3.14159f * node.radius * node.radius * length;
                cost = g.sugar_maintenance_stem * volume;
                break;
            }
            case NodeType::ROOT: {
                float volume = 3.14159f * node.radius * node.radius * length;
                cost = g.sugar_maintenance_root * volume;
                break;
            }
        }
        if (node.meristem && node.meristem->is_tip() && node.meristem->active) {
            cost += g.sugar_maintenance_meristem;
        }
        node.sugar = std::max(0.0f, node.sugar - cost);

        // Track starvation
        if (node.sugar <= 0.0f) {
            node.starvation_ticks++;
        } else {
            node.starvation_ticks = 0;
        }
    });
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[sugar]"`

Expected: All sugar tests PASS including updated maintenance tests.

- [ ] **Step 5: Commit**

```bash
git add src/engine/sugar.cpp tests/test_sugar.cpp
git commit -m "feat: volume-based maintenance — cost scales with tissue volume not radius"
```

---

### Task 3: Production cap enforcement

**Files:**
- Modify: `src/engine/sugar.cpp:12-57` — add cap check in `produce_sugar()`
- Test: `tests/test_sugar.cpp`

- [ ] **Step 1: Write failing tests for production cap**

Add to `tests/test_sugar.cpp`:

```cpp
TEST_CASE("Production skipped when leaf sugar is at cap", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.3f;
    plant.seed_mut()->add_child(leaf);

    // Fill leaf to its cap
    float cap = sugar_cap(*leaf, g);
    leaf->sugar = cap;

    WorldParams wp = default_world_params();
    produce_sugar(plant, wp);

    // Sugar should not have increased
    REQUIRE_THAT(leaf->sugar, WithinAbs(cap, 1e-6));
}

TEST_CASE("Production works normally when leaf is below cap", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.5f;
    leaf->sugar = 0.0f;
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 2.0f;

    produce_sugar(plant, wp);

    REQUIRE(leaf->sugar > 0.0f);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[sugar]" -c "cap"`

Expected: "Production skipped when leaf sugar is at cap" FAILS — sugar increases past cap.

- [ ] **Step 3: Add cap enforcement in produce_sugar**

In `src/engine/sugar.cpp`, in the `produce_sugar` function, add a cap check inside the leaf loop. Replace the production line (around line 54):

```cpp
        leaf.node->light_exposure = std::exp(-k * shade);

        // Feedback inhibition: full storage stops photosynthesis
        float cap = sugar_cap(*leaf.node, g);
        if (leaf.node->sugar >= cap) continue;

        leaf.node->sugar += leaf.node->light_exposure * world.light_level
                          * leaf.node->leaf_size * g.sugar_production_rate;
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[sugar]"`

Expected: All sugar tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/engine/sugar.cpp tests/test_sugar.cpp
git commit -m "feat: production cap — leaves stop photosynthesizing when sugar storage is full"
```

---

### Task 4: Cap-aware diffusion

**Files:**
- Modify: `src/engine/sugar.cpp:88-138` — cap-aware transfers in `diffuse_sugar()`
- Test: `tests/test_sugar.cpp`

- [ ] **Step 1: Write failing test for cap-aware diffusion**

Add to `tests/test_sugar.cpp`:

```cpp
TEST_CASE("Diffusion does not push sugar past receiver cap", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Create a child node with small cap (tiny radius, short offset)
    Node* child = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.1f, 0.0f), 0.02f);
    plant.seed_mut()->add_child(child);

    float child_cap = sugar_cap(*child, g);

    // Give seed lots of sugar, child already at cap
    plant.seed_mut()->sugar = 100.0f;
    child->sugar = child_cap;

    WorldParams wp = default_world_params();
    wp.sugar_diffusion_iterations = 10;
    diffuse_sugar(plant, wp);

    // Child should not exceed its cap
    REQUIRE(child->sugar <= child_cap + 1e-6f);
}

TEST_CASE("Diffusion still conserves sugar when caps are not hit", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Build a tree with nodes that have large caps (big radii, long offsets)
    for (int i = 0; i < 3; i++) {
        Node* child = plant.create_node(NodeType::STEM,
            glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
        plant.seed_mut()->add_child(child);
    }

    // Give modest sugar amounts (well below caps)
    plant.seed_mut()->sugar = 5.0f;
    plant.seed_mut()->children[0]->sugar = 2.0f;

    float total_before = 0.0f;
    plant.for_each_node([&](const Node& n) { total_before += n.sugar; });

    WorldParams wp = default_world_params();
    wp.sugar_diffusion_iterations = 10;
    diffuse_sugar(plant, wp);

    float total_after = 0.0f;
    plant.for_each_node([&](const Node& n) { total_after += n.sugar; });

    REQUIRE_THAT(total_after, WithinAbs(total_before, 1e-4));
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[sugar]" -c "Diffusion does not push"`

Expected: FAIL — child sugar exceeds cap.

- [ ] **Step 3: Make diffusion cap-aware**

In `src/engine/sugar.cpp`, update the `diffuse_sugar` function. Add cap fields to the Edge struct and cap-aware clamping.

Replace the Edge struct and edge-building loop:

```cpp
    struct Edge {
        Node* a;
        Node* b;
        float conductance;
        float a_cap;
        float b_cap;
    };
    std::vector<Edge> edges;
    plant.for_each_node_mut([&](Node& node) {
        if (node.parent) {
            float min_radius = std::min(node.radius, node.parent->radius);
            if (node.type == NodeType::LEAF || node.parent->type == NodeType::LEAF) {
                min_radius = std::max(min_radius, 0.01f);
            }
            float area = min_radius * min_radius * 3.14159f;
            float conductance = std::min(area * g.sugar_transport_conductance, 0.25f);
            edges.push_back({&node, node.parent, conductance,
                             sugar_cap(node, g), sugar_cap(*node.parent, g)});
        }
    });
```

Replace the inner diffusion loop body:

```cpp
    for (int iter = 0; iter < world.sugar_diffusion_iterations; iter++) {
        for (Edge& e : edges) {
            float diff = e.a->sugar - e.b->sugar;
            float transfer = diff * e.conductance;

            if (transfer > 0.0f) {
                // a -> b: clamp by a's available sugar and b's headroom
                float b_headroom = std::max(0.0f, e.b_cap - e.b->sugar);
                transfer = std::min({transfer, e.a->sugar, b_headroom});
            } else {
                // b -> a: clamp by b's available sugar and a's headroom
                float a_headroom = std::max(0.0f, e.a_cap - e.a->sugar);
                transfer = std::max({transfer, -e.b->sugar, -a_headroom});
            }

            e.a->sugar -= transfer;
            e.b->sugar += transfer;
        }
    }
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[sugar]"`

Expected: All sugar tests PASS, including both new cap-aware tests and the existing conservation tests.

- [ ] **Step 5: Commit**

```bash
git add src/engine/sugar.cpp tests/test_sugar.cpp
git commit -m "feat: cap-aware diffusion — sugar transfer respects receiver storage limits"
```

---

### Task 5: Safety clamp + CLAUDE.md update

**Files:**
- Modify: `src/engine/sugar.cpp` — add safety clamp in `consume_sugar()`
- Modify: `CLAUDE.md` — update Sugar Model documentation
- Test: `tests/test_sugar.cpp`

- [ ] **Step 1: Write failing test for safety clamp**

Add to `tests/test_sugar.cpp`:

```cpp
TEST_CASE("Safety clamp caps sugar after consumption", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.1f);
    plant.seed_mut()->add_child(stem);

    float cap = sugar_cap(*stem, g);
    // Set sugar way above cap (simulating an edge case)
    stem->sugar = cap * 10.0f;

    consume_sugar(plant);

    REQUIRE(stem->sugar <= cap + 1e-6f);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[sugar]" -c "Safety clamp"`

Expected: FAIL — sugar remains above cap (only maintenance was deducted).

- [ ] **Step 3: Add safety clamp in consume_sugar**

In `src/engine/sugar.cpp`, in the `consume_sugar` function, add the clamp after the `max(0, sugar - cost)` line, before starvation tracking:

```cpp
        node.sugar = std::max(0.0f, node.sugar - cost);

        // Safety clamp: cap sugar to node's storage capacity
        float cap = sugar_cap(node, g);
        node.sugar = std::min(node.sugar, cap);

        // Track starvation
```

- [ ] **Step 4: Run all tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`

Expected: ALL tests PASS (sugar tests, meristem tests, genome tests).

- [ ] **Step 5: Update CLAUDE.md sugar model documentation**

In `CLAUDE.md`, replace the `## Sugar Model` section with:

```markdown
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
```

- [ ] **Step 6: Run full build and all tests one final time**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`

Expected: ALL tests PASS. Zero warnings.

- [ ] **Step 7: Commit**

```bash
git add src/engine/sugar.cpp tests/test_sugar.cpp CLAUDE.md
git commit -m "feat: safety clamp + docs — complete sugar realism implementation"
```
