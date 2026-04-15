# Auxin Canalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add per-connection transport biases driven by auxin flow history — transient (fast, reversible) and structural (slow, permanent) — that redistribute chemical flow among siblings at junctions.

**Architecture:** Two `unordered_map<Node*, float>` on each parent node store bias values per child. During transport, biases multiply the proportional weights used to split flow among children (redistribution, not amplification). Auxin flux captured during transport updates biases after each tick.

**Tech Stack:** C++17, Catch2, CMake

---

### Task 1: Genome parameters for canalization

**Files:**
- Modify: `src/engine/genome.h`

- [ ] **Step 1: Add 6 new fields to the Genome struct**

In `src/engine/genome.h`, add after the `elastic_recovery_rate` field (end of the Stress section):

```cpp
    // Canalization — auxin flow history biases transport
    float transient_gain;                 // target bias per unit of auxin flux
    float transient_rate;                 // how fast transient bias chases its target (0-1)
    float structural_threshold;           // minimum auxin flux to grow structural bias
    float structural_growth_rate;         // structural bias increment per tick above threshold
    float structural_max;                 // cap on structural bias
    float canalization_weight;            // global scaling on combined bias effect (0 = disabled)
```

- [ ] **Step 2: Add defaults to default_genome()**

Add after the `.elastic_recovery_rate` initializer:

```cpp
        // Canalization
        .transient_gain = 2.0f,               // target = flux * 2 — responsive
        .transient_rate = 0.2f,               // ~87% in 8 hours
        .structural_threshold = 0.05f,        // min flux for vascular development
        .structural_growth_rate = 0.005f,     // ~8 days to reach 1.0
        .structural_max = 2.0f,               // at max: 1 + 2.0 = 3.0x weight
        .canalization_weight = 1.0f,          // full effect by default
```

- [ ] **Step 3: Build and run tests**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests
```

Expected: all tests pass. No behavioral change — fields exist but aren't used yet.

- [ ] **Step 4: Commit**

```bash
git add src/engine/genome.h
git commit -m "feat: add canalization genome parameters (6 new fields)"
```

---

### Task 2: Bias maps and lifecycle management on Node

**Files:**
- Modify: `src/engine/node/node.h`
- Modify: `src/engine/node/node.cpp`
- Modify: `tests/test_node.cpp` (or `tests/test_meristem.cpp` if test_node doesn't exist)

- [ ] **Step 1: Write tests for bias lifecycle**

Check if `tests/test_node.cpp` exists. If not, create it with proper Catch2 includes. Add these tests:

```cpp
#include <catch2/catch_test_macros.hpp>
#include "engine/plant.h"
#include "engine/world_params.h"

using namespace botany;

TEST_CASE("Canalization: replace_child transfers biases", "[canalization]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* parent_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* old_child = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* new_child = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);

    plant.seed_mut()->add_child(parent_node);
    parent_node->add_child(old_child);

    // Manually set biases on the parent for old_child
    parent_node->auxin_flow_bias[old_child] = 0.5f;
    parent_node->structural_flow_bias[old_child] = 0.3f;

    // Replace old_child with new_child
    parent_node->replace_child(old_child, new_child);

    // Biases should transfer to new_child
    REQUIRE(parent_node->auxin_flow_bias.count(new_child) == 1);
    REQUIRE(parent_node->auxin_flow_bias[new_child] == 0.5f);
    REQUIRE(parent_node->structural_flow_bias[new_child] == 0.3f);

    // Old child entries should be gone
    REQUIRE(parent_node->auxin_flow_bias.count(old_child) == 0);
    REQUIRE(parent_node->structural_flow_bias.count(old_child) == 0);
}

TEST_CASE("Canalization: die cleans up parent bias entries", "[canalization]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* parent_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* child_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);

    plant.seed_mut()->add_child(parent_node);
    parent_node->add_child(child_node);

    // Set biases
    parent_node->auxin_flow_bias[child_node] = 0.7f;
    parent_node->structural_flow_bias[child_node] = 0.4f;

    // Kill the child
    child_node->die(plant);

    // Parent's bias entries for dead child should be gone
    REQUIRE(parent_node->auxin_flow_bias.count(child_node) == 0);
    REQUIRE(parent_node->structural_flow_bias.count(child_node) == 0);
}

TEST_CASE("Canalization: get_bias_multiplier returns 1.0 with no bias", "[canalization]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* parent_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* child_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    plant.seed_mut()->add_child(parent_node);
    parent_node->add_child(child_node);

    // No biases set — multiplier should be 1.0
    REQUIRE(parent_node->get_bias_multiplier(child_node, g) == 1.0f);
}

TEST_CASE("Canalization: get_bias_multiplier scales with canalization_weight", "[canalization]") {
    Genome g = default_genome();
    g.canalization_weight = 2.0f;
    Plant plant(g, glm::vec3(0.0f));

    Node* parent_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* child_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    plant.seed_mut()->add_child(parent_node);
    parent_node->add_child(child_node);

    parent_node->auxin_flow_bias[child_node] = 0.5f;
    parent_node->structural_flow_bias[child_node] = 0.3f;

    // 1.0 + 2.0 * (0.5 + 0.3) = 2.6
    float expected = 1.0f + 2.0f * (0.5f + 0.3f);
    REQUIRE(parent_node->get_bias_multiplier(child_node, g) == expected);
}
```

- [ ] **Step 2: Build and verify tests fail**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests "[canalization]"
```

Expected: compilation fails — `auxin_flow_bias`, `structural_flow_bias`, `get_bias_multiplier` don't exist yet.

- [ ] **Step 3: Add bias maps and helper to node.h**

In `src/engine/node/node.h`, add to the public section of Node, after the `chemicals` map:

```cpp
    // Canalization — per-child transport bias (stored on parent, keyed by child pointer)
    std::unordered_map<Node*, float> auxin_flow_bias;       // transient — fast, decays
    std::unordered_map<Node*, float> structural_flow_bias;  // persistent — slow, permanent
    std::unordered_map<Node*, float> last_auxin_flux;       // transient per-tick: auxin moved per child

    float get_bias_multiplier(Node* child, const Genome& g) const;
```

- [ ] **Step 4: Implement get_bias_multiplier, update replace_child and die in node.cpp**

Add `get_bias_multiplier` implementation:

```cpp
float Node::get_bias_multiplier(Node* child, const Genome& g) const {
    float flow = 0.0f, structural = 0.0f;
    auto it_f = auxin_flow_bias.find(child);
    if (it_f != auxin_flow_bias.end()) flow = it_f->second;
    auto it_s = structural_flow_bias.find(child);
    if (it_s != structural_flow_bias.end()) structural = it_s->second;
    return 1.0f + g.canalization_weight * (flow + structural);
}
```

In `replace_child()`, after `old_child->parent = nullptr;`, add bias transfer:

```cpp
            // Transfer canalization biases to new child entry
            auto it_flow = auxin_flow_bias.find(old_child);
            if (it_flow != auxin_flow_bias.end()) {
                auxin_flow_bias[new_child] = it_flow->second;
                auxin_flow_bias.erase(it_flow);
            }
            auto it_struct = structural_flow_bias.find(old_child);
            if (it_struct != structural_flow_bias.end()) {
                structural_flow_bias[new_child] = it_struct->second;
                structural_flow_bias.erase(it_struct);
            }
```

In `die()`, before the `auto& siblings = parent->children;` line, add:

```cpp
        parent->auxin_flow_bias.erase(this);
        parent->structural_flow_bias.erase(this);
```

- [ ] **Step 5: Build and run tests**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests
```

Expected: ALL tests pass including the 4 new canalization tests.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/node.h src/engine/node/node.cpp tests/test_node.cpp
git commit -m "feat: add canalization bias maps and lifecycle management on Node"
```

---

### Task 3: Bias-weighted transport distribution

**Files:**
- Modify: `tests/test_node.cpp`
- Modify: `src/engine/node/node.cpp`

- [ ] **Step 1: Write failing test for bias-weighted distribution**

Add to the canalization test file:

```cpp
TEST_CASE("Canalization: biased child gets larger sugar share", "[canalization]") {
    Genome g = default_genome();
    g.canalization_weight = 1.0f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world_params();

    // Build: seed -> parent -> [childA, childB]
    Node* parent_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* childA = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* childB = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    plant.seed_mut()->add_child(parent_node);
    parent_node->add_child(childA);
    parent_node->add_child(childB);

    // Give parent lots of sugar, children none — parent will distribute to both
    parent_node->chemical(ChemicalID::Sugar) = 10.0f;
    childA->chemical(ChemicalID::Sugar) = 0.0f;
    childB->chemical(ChemicalID::Sugar) = 0.0f;

    // Bias childA's connection strongly
    parent_node->auxin_flow_bias[childA] = 1.0f;    // multiplier = 1 + 1*(1+0) = 2.0
    parent_node->auxin_flow_bias[childB] = 0.0f;    // multiplier = 1 + 1*(0+0) = 1.0

    // Run one tick of transport only (not full plant tick to avoid side effects)
    parent_node->transport_with_children(g);

    // ChildA should have received more sugar than childB
    REQUIRE(childA->chemical(ChemicalID::Sugar) > childB->chemical(ChemicalID::Sugar) * 1.5f);
}

TEST_CASE("Canalization: zero canalization_weight disables bias", "[canalization]") {
    Genome g = default_genome();
    g.canalization_weight = 0.0f;  // disabled
    Plant plant(g, glm::vec3(0.0f));

    Node* parent_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* childA = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* childB = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    plant.seed_mut()->add_child(parent_node);
    parent_node->add_child(childA);
    parent_node->add_child(childB);

    parent_node->chemical(ChemicalID::Sugar) = 10.0f;
    childA->chemical(ChemicalID::Sugar) = 0.0f;
    childB->chemical(ChemicalID::Sugar) = 0.0f;

    // Set big biases — but weight=0 should ignore them
    parent_node->auxin_flow_bias[childA] = 5.0f;
    parent_node->auxin_flow_bias[childB] = 0.0f;

    parent_node->transport_with_children(g);

    // With weight=0, both children should get roughly equal sugar
    float ratio = childA->chemical(ChemicalID::Sugar) /
                  std::max(childB->chemical(ChemicalID::Sugar), 1e-8f);
    REQUIRE(ratio > 0.8f);
    REQUIRE(ratio < 1.2f);
}

TEST_CASE("Canalization: bias affects Phase 1 (children giving to parent)", "[canalization]") {
    Genome g = default_genome();
    g.canalization_weight = 1.0f;
    Plant plant(g, glm::vec3(0.0f));

    // Build: seed -> parent -> [childA, childB]
    Node* parent_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* childA = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* childB = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    plant.seed_mut()->add_child(parent_node);
    parent_node->add_child(childA);
    parent_node->add_child(childB);

    // Children have lots of auxin, parent has none — both want to give
    // Set parent sugar cap very low to create headroom constraint
    childA->chemical(ChemicalID::Auxin) = 5.0f;
    childB->chemical(ChemicalID::Auxin) = 5.0f;
    parent_node->chemical(ChemicalID::Auxin) = 0.0f;

    // Bias childA strongly
    parent_node->auxin_flow_bias[childA] = 2.0f;  // multiplier = 1+1*(2+0) = 3.0
    parent_node->auxin_flow_bias[childB] = 0.0f;  // multiplier = 1.0

    float a_before = childA->chemical(ChemicalID::Auxin);
    float b_before = childB->chemical(ChemicalID::Auxin);

    parent_node->transport_with_children(g);

    // ChildA should have given more auxin than childB (lost more)
    float a_gave = a_before - childA->chemical(ChemicalID::Auxin);
    float b_gave = b_before - childB->chemical(ChemicalID::Auxin);

    // Both should have given something
    REQUIRE(a_gave > 0.0f);
    REQUIRE(b_gave > 0.0f);
    // ChildA (higher bias) should have given more
    REQUIRE(a_gave > b_gave * 1.3f);
}
```

- [ ] **Step 2: Build and verify tests fail**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests "[canalization]"
```

Expected: new tests FAIL — transport doesn't use biases yet.

- [ ] **Step 3: Implement bias-weighted transport in transport_with_children**

In `src/engine/node/node.cpp`, modify `transport_with_children()`:

**3a. Add `bias_mult` to ChildFlow struct:**

```cpp
        struct ChildFlow {
            Node* child;
            float desired;   // positive = child wants to receive, negative = child wants to give
            float child_cap;
            float bias_mult; // canalization weight multiplier for this connection
        };
```

**3b. Compute bias_mult when creating each ChildFlow:**

After computing `desired`, change the push_back to include bias_mult:

```cpp
            if (std::abs(desired) > 1e-8f) {
                float bm = get_bias_multiplier(child, g);
                flows.push_back({child, desired, child_cap, bm});
            }
```

**3c. Replace uniform Phase 1 scaling with bias-weighted redistribution:**

Replace the Phase 1 scaling block:

```cpp
        // Scale if total exceeds parent's capacity
        if (total_inflow > parent_headroom && total_inflow > 1e-8f) {
            float scale = parent_headroom / total_inflow;
            for (auto& cf : flows) {
                if (cf.desired < 0.0f) cf.desired *= scale;
            }
            total_inflow = parent_headroom;
        }
```

with:

```cpp
        // Scale if total exceeds parent's capacity — bias-weighted redistribution
        if (total_inflow > parent_headroom && total_inflow > 1e-8f) {
            float total_weighted = 0.0f;
            for (auto& cf : flows) {
                if (cf.desired >= 0.0f) continue;
                total_weighted += (-cf.desired) * cf.bias_mult;
            }
            if (total_weighted > 1e-8f) {
                for (auto& cf : flows) {
                    if (cf.desired >= 0.0f) continue;
                    float raw_give = -cf.desired;
                    float share = parent_headroom * (raw_give * cf.bias_mult / total_weighted);
                    cf.desired = -std::min(share, raw_give);
                }
            }
            total_inflow = parent_headroom;
        }
```

**3d. Use bias-weighted weights in Phase 2 receiver setup:**

Change the Receiver struct and population:

```cpp
        struct Receiver {
            Node* child;
            float raw_weight;   // unbiased desired flow
            float weight;       // bias-adjusted for proportional split
            float headroom;
        };
        std::vector<Receiver> receivers;
        for (auto& cf : flows) {
            if (cf.desired <= 0.0f) continue;
            float headroom = has_cap
                ? std::max(0.0f, cf.child_cap - cf.child->chemical(dp.id))
                : 1e30f;
            if (headroom > 1e-8f) {
                receivers.push_back({cf.child, cf.desired, cf.desired * cf.bias_mult, headroom});
            }
        }
```

**3e. Use raw_total for budget, bias_total for split in Phase 2 loop:**

Replace:

```cpp
        while (!receivers.empty() && available > 1e-8f) {
            float total_weight = 0.0f;
            for (const auto& r : receivers) total_weight += r.weight;
            if (total_weight < 1e-8f) break;

            float to_distribute = std::min(available, total_weight);
```

with:

```cpp
        while (!receivers.empty() && available > 1e-8f) {
            float raw_total = 0.0f;
            float bias_total = 0.0f;
            for (const auto& r : receivers) {
                raw_total += r.raw_weight;
                bias_total += r.weight;
            }
            if (bias_total < 1e-8f) break;

            float to_distribute = std::min(available, raw_total);
```

And update the share calculation:

```cpp
                float share = to_distribute * (r.weight / bias_total);
```

(This line is unchanged from current code — `r.weight` is now the bias-adjusted weight.)

- [ ] **Step 4: Build and run tests**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests
```

Expected: ALL tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/engine/node/node.cpp tests/test_node.cpp
git commit -m "feat: bias-weighted transport distribution (canalization Phase 1 & 2)"
```

---

### Task 4: Auxin flux capture and canalization bias update

**Files:**
- Modify: `tests/test_node.cpp`
- Modify: `src/engine/node/node.h`
- Modify: `src/engine/node/node.cpp`

- [ ] **Step 1: Write tests for bias update mechanics**

Add to canalization test file:

```cpp
TEST_CASE("Canalization: transient bias builds from auxin flux", "[canalization]") {
    Genome g = default_genome();
    g.transient_gain = 2.0f;
    g.transient_rate = 0.5f;  // fast for testing
    g.canalization_weight = 1.0f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world_params();

    // Give shoot apical lots of sugar so it produces auxin and grows
    for (int i = 0; i < 10; i++) {
        plant.for_each_node_mut([](Node& n) {
            n.chemical(ChemicalID::Sugar) = 100.0f;
            n.chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant.tick(world);
    }

    // The seed should have non-zero transient bias toward the shoot apical child
    // (auxin flows basipetally from shoot to seed)
    bool found_nonzero_bias = false;
    Node* seed = plant.seed_mut();
    for (auto& [child, bias] : seed->auxin_flow_bias) {
        if (bias > 0.01f) found_nonzero_bias = true;
    }
    REQUIRE(found_nonzero_bias);
}

TEST_CASE("Canalization: transient bias decays without flux", "[canalization]") {
    Genome g = default_genome();
    g.transient_gain = 2.0f;
    g.transient_rate = 0.5f;
    g.apical_auxin_baseline = 0.0f;  // no auxin production
    g.canalization_weight = 1.0f;
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    Node* child = nullptr;
    for (Node* c : seed->children) {
        if (c->type == NodeType::APICAL) { child = c; break; }
    }
    REQUIRE(child != nullptr);

    // Manually set a transient bias
    seed->auxin_flow_bias[child] = 1.0f;

    // Tick with no auxin production — bias should decay toward 0
    WorldParams world = default_world_params();
    for (int i = 0; i < 10; i++) {
        plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; });
        plant.tick(world);
    }

    // Bias should have decayed significantly (target=0, rate=0.5)
    REQUIRE(seed->auxin_flow_bias[child] < 0.1f);
}

TEST_CASE("Canalization: structural bias builds with sustained flux", "[canalization]") {
    Genome g = default_genome();
    g.structural_threshold = 0.01f;  // low for testing
    g.structural_growth_rate = 0.1f; // fast for testing
    g.structural_max = 5.0f;
    g.canalization_weight = 1.0f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world_params();

    // Tick with auxin production active
    for (int i = 0; i < 10; i++) {
        plant.for_each_node_mut([](Node& n) {
            n.chemical(ChemicalID::Sugar) = 100.0f;
            n.chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant.tick(world);
    }

    // Seed should have structural bias toward shoot child
    bool found_structural = false;
    Node* seed = plant.seed_mut();
    for (auto& [child, bias] : seed->structural_flow_bias) {
        if (bias > 0.01f) found_structural = true;
    }
    REQUIRE(found_structural);
}

TEST_CASE("Canalization: structural bias does NOT decay", "[canalization]") {
    Genome g = default_genome();
    g.structural_threshold = 0.01f;
    g.structural_growth_rate = 0.1f;
    g.structural_max = 5.0f;
    g.apical_auxin_baseline = 0.0f;  // will disable after building bias
    g.canalization_weight = 1.0f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world_params();

    Node* seed = plant.seed_mut();
    Node* shoot = nullptr;
    for (Node* c : seed->children) {
        if (c->type == NodeType::APICAL) { shoot = c; break; }
    }
    REQUIRE(shoot != nullptr);

    // Manually set a structural bias (simulating past buildup)
    seed->structural_flow_bias[shoot] = 1.5f;
    float bias_before = seed->structural_flow_bias[shoot];

    // Tick with zero auxin — no flux, but structural should NOT decay
    for (int i = 0; i < 20; i++) {
        plant.for_each_node_mut([](Node& n) {
            n.chemical(ChemicalID::Sugar) = 100.0f;
            n.chemical(ChemicalID::Auxin) = 0.0f;
        });
        plant.tick(world);
    }

    REQUIRE(seed->structural_flow_bias[shoot] >= bias_before);
}

TEST_CASE("Canalization: structural bias capped at structural_max", "[canalization]") {
    Genome g = default_genome();
    g.structural_max = 1.0f;
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    Node* child = nullptr;
    for (Node* c : seed->children) { child = c; break; }
    REQUIRE(child != nullptr);

    // Set bias above max
    seed->structural_flow_bias[child] = 5.0f;

    // update_canalization should clamp to max
    seed->update_canalization(plant.genome());

    REQUIRE(seed->structural_flow_bias[child] <= g.structural_max);
}
```

- [ ] **Step 2: Build and verify tests fail**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests "[canalization]"
```

Expected: compilation fails — `update_canalization` doesn't exist yet. Some tests may fail even after compile because flux isn't captured.

- [ ] **Step 3: Add update_canalization declaration to node.h**

In `src/engine/node/node.h`, add to the public section:

```cpp
    void update_canalization(const Genome& g);
```

- [ ] **Step 4: Implement auxin flux capture in transport_with_children**

In `src/engine/node/node.cpp`, at the start of `transport_with_children()`, after `if (children.empty()) return;`, add:

```cpp
    last_auxin_flux.clear();
```

After the Phase 1 "Apply inflows" block (where `cf.child->chemical(dp.id) -= give;`), add flux capture:

```cpp
            // Capture auxin flux for canalization
            if (dp.id == ChemicalID::Auxin) {
                last_auxin_flux[cf.child] += give;
            }
```

In the Phase 2 distribution loop, after `r.child->chemical(dp.id) += actual;`, add:

```cpp
                if (dp.id == ChemicalID::Auxin) {
                    last_auxin_flux[r.child] += actual;
                }
```

- [ ] **Step 5: Implement update_canalization**

Add to `src/engine/node/node.cpp`:

```cpp
void Node::update_canalization(const Genome& g) {
    for (Node* child : children) {
        // Get auxin flux for this child (0 if none recorded)
        float flux = 0.0f;
        auto it = last_auxin_flux.find(child);
        if (it != last_auxin_flux.end()) flux = it->second;

        // Transient bias: exponential chase toward flux-derived target
        float target = flux * g.transient_gain;
        float& flow_bias = auxin_flow_bias[child];
        flow_bias += (target - flow_bias) * g.transient_rate;

        // Structural bias: slow ratchet, never decays
        float& struct_bias = structural_flow_bias[child];
        if (flux > g.structural_threshold) {
            struct_bias += g.structural_growth_rate;
        }
        struct_bias = std::min(struct_bias, g.structural_max);
    }
}
```

- [ ] **Step 6: Wire update_canalization into update_chemicals**

In `src/engine/node/node.cpp`, change `update_chemicals()`:

```cpp
void Node::update_chemicals(const Genome& g) {
    transport_with_children(g);
    update_canalization(g);
    decay_chemicals(g);
}
```

- [ ] **Step 7: Build and run tests**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests
```

Expected: ALL tests pass.

- [ ] **Step 8: Commit**

```bash
git add src/engine/node/node.h src/engine/node/node.cpp tests/test_node.cpp
git commit -m "feat: auxin flux capture and canalization bias update"
```

---

### Task 5: Evolution bridge update

**Files:**
- Modify: `src/evolution/genome_bridge.cpp`

- [ ] **Step 1: Register 6 new genes in build_genome_template()**

Add after the stress group `reg()` calls (before the linkage groups):

```cpp
    // --- Canalization group (6 genes) ---
    reg(sg, "transient_gain",           g.transient_gain,           r, 0.0f, 20.0f, p);
    reg(sg, "transient_rate",           g.transient_rate,           r, 0.0f, 1.0f, p);
    reg(sg, "structural_threshold",     g.structural_threshold,     r, 0.0f, 1.0f, p);
    reg(sg, "structural_growth_rate",   g.structural_growth_rate,   r, 0.0f, 0.1f, p);
    reg(sg, "structural_max",           g.structural_max,           r, 0.0f, 10.0f, p);
    reg(sg, "canalization_weight",      g.canalization_weight,      r, 0.0f, 5.0f, p);
```

- [ ] **Step 2: Add from_structured() mappings**

Add after the stress section in `from_structured()`:

```cpp
    // Canalization
    g.transient_gain           = sg.get("transient_gain");
    g.transient_rate           = sg.get("transient_rate");
    g.structural_threshold     = sg.get("structural_threshold");
    g.structural_growth_rate   = sg.get("structural_growth_rate");
    g.structural_max           = sg.get("structural_max");
    g.canalization_weight      = sg.get("canalization_weight");
```

- [ ] **Step 3: Add canalization linkage group**

Add after the stress linkage group:

```cpp
    sg.add_linkage_group({"canalization", {
        "transient_gain", "transient_rate",
        "structural_threshold", "structural_growth_rate", "structural_max",
        "canalization_weight"
    }});
```

- [ ] **Step 4: Build and run tests**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests
```

Expected: ALL tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/evolution/genome_bridge.cpp
git commit -m "feat: add canalization params to evolution bridge"
```

---

### Task 6: Update CLAUDE.md

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Add canalization section**

In `CLAUDE.md`, add a new section after the "Ethylene Model" section:

```
## Canalization Model

Auxin flow through parent-child connections builds transport preference over time. Two layers of memory:

**Transient bias (auxin_flow_bias)** — fast, reversible. Represents PIN protein polarization.
- Each tick: `target = auxin_flux * transient_gain`, bias chases target exponentially at `transient_rate`
- Decays toward 0 when flux stops. Responds within hours/days.

**Structural bias (structural_flow_bias)** — slow, permanent. Represents built xylem/phloem.
- Grows by `structural_growth_rate` per tick when auxin flux exceeds `structural_threshold`
- Never decays. Capped at `structural_max`.

**Effect on transport:** Both biases stored on the parent node, keyed by child pointer. During `transport_with_children()`, each child's proportional weight is multiplied by `1 + canalization_weight * (flow_bias + structural_bias)`. This redistributes chemical flow (all chemicals, not just auxin) among siblings — biased connections get a larger share of the same total. Does not amplify total flow.

**Lifecycle:** Biases transfer on `replace_child` (chain growth preserves branch history). New children start at 0, 0. Cleaned up on `die()`.
```

- [ ] **Step 2: Add genome parameters to tuning section**

In the "Tuning Parameters" section, add after the ethylene entries:

```
- `transient_gain` (2.0) - target transient bias per unit of auxin flux. Higher = stronger short-term flow reinforcement
- `transient_rate` (0.2) - exponential chase speed for transient bias (0.2 = ~87% in 8 hours)
- `structural_threshold` (0.05) - minimum auxin flux to build structural bias (filters noise)
- `structural_growth_rate` (0.005) - structural bias growth per tick when flux exceeds threshold (~8 days to reach 1.0)
- `structural_max` (2.0) - cap on structural bias. At max: connection gets 1 + 2.0 = 3x weight
- `canalization_weight` (1.0) - global scaling on bias effect. 0 = canalization disabled entirely
```

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: add canalization model to CLAUDE.md"
```
