# Plant Compression (Node Fusion) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Periodically merge adjacent stem/root nodes that pass structural similarity checks so mature plants shed redundant bookkeeping without losing architectural meaning.

**Architecture:** New `src/engine/compression.{h,cpp}` module owns all merge logic as a pure function `compress_plant(Plant&, params)` that runs *between* ticks. `Engine` gets autocompress plumbing (enable/interval/params/trigger/last-result) and hooks the call into `tick()` end-of-cycle. Realtime viewer gets an ImGui CollapsingHeader with checkbox+interval or button+stats.

**Tech Stack:** C++17, CMake with FetchContent, Catch2 v3, ImGui, GLFW.

---

## File structure

**Create:**
- `src/engine/compression.h` — `CompressionParams`, `CompressionResult`, `can_merge()`, `compress_plant()`.
- `src/engine/compression.cpp` — merge predicate, merge execution, multi-pass scan.
- `tests/test_compression.cpp` — unit tests for every merge rule + end-to-end compress_plant tests + post-compression save/load smoke.

**Modify:**
- `src/engine/engine.h` / `engine.cpp` — new compression controls on Engine + tick() hook.
- `src/app_realtime.cpp` — ImGui CollapsingHeader "Compression" inside Controls.
- `CMakeLists.txt` — add `src/engine/compression.cpp` to `botany_engine`, add `tests/test_compression.cpp` to `botany_tests`.
- `CLAUDE.md` — document compression feature.

---

## Task 1: Module skeleton + stub `compress_plant`

**Files:**
- Create: `src/engine/compression.h`
- Create: `src/engine/compression.cpp`
- Create: `tests/test_compression.cpp`
- Modify: `CMakeLists.txt`

- [ ] **Step 1.1: Write the failing test**

Create `tests/test_compression.cpp`:
```cpp
#include <catch2/catch_test_macros.hpp>
#include <glm/vec3.hpp>

#include "engine/compression.h"
#include "engine/plant.h"
#include "engine/genome.h"

using namespace botany;

TEST_CASE("compress_plant on a fresh young plant reports no merges", "[compression][stub]") {
    Plant plant(default_genome(), glm::vec3(0.0f));
    CompressionParams params;
    CompressionResult r = compress_plant(plant, params);
    REQUIRE(r.merges_performed == 0u);
    REQUIRE(r.passes_run == 0u);
}
```

- [ ] **Step 1.2: Update CMake**

Edit `CMakeLists.txt`:
1. In the `botany_engine` source list (the `add_library(botany_engine STATIC ...)` block), add `src/engine/compression.cpp` in alphabetical order (after `compartments.*` entries if any, otherwise after `vascular.cpp`).
2. In the `botany_tests` source list (the `add_executable(botany_tests ...)` block), add `tests/test_compression.cpp`.

- [ ] **Step 1.3: Run test to verify it fails**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5
```
Expected: compile error — `engine/compression.h` not found.

- [ ] **Step 1.4: Create the header**

Create `src/engine/compression.h`:
```cpp
#pragma once

#include <cstdint>

namespace botany {

class Plant;
class Node;
struct Genome;

// Tunable thresholds for compression.  See design doc
// docs/superpowers/specs/2026-04-23-plant-compression-design.md.
struct CompressionParams {
    float    max_angle_rad       = 0.175f;  // ≈ 10°
    float    max_radius_ratio    = 0.20f;   // |r1-r2|/max(r1,r2)
    float    max_combined_length = 0.0f;    // 0 → auto = 2 × g.max_internode_length
    uint32_t max_passes          = 5;       // multi-pass convergence cap
};

// Summary of a compress_plant run.  delta_* are negative when cap clamping
// discards chemicals; zero if the merge was mass-conservative.
struct CompressionResult {
    uint32_t merges_performed = 0;
    uint32_t passes_run       = 0;
    float    delta_sugar      = 0.0f;
    float    delta_water      = 0.0f;
    float    delta_auxin      = 0.0f;
    float    delta_cytokinin  = 0.0f;
};

// Predicate: can parent P absorb its single structural child C under params?
// Used by compress_plant's scan and exposed here for fine-grained testing.
bool can_merge(const Node& parent, const Node& child,
               const Genome& g, const CompressionParams& params);

// Main entry.  Scans plant and merges every accepted adjacent (parent, child)
// stem/root pair in up to params.max_passes iterations, flushing deferred
// removals between passes.  Runs mutations directly on the plant; must only
// be called between full ticks (never from inside Plant::tick / Engine::tick).
CompressionResult compress_plant(Plant& plant, const CompressionParams& params);

} // namespace botany
```

- [ ] **Step 1.5: Create the stub implementation**

Create `src/engine/compression.cpp`:
```cpp
#include "engine/compression.h"
#include "engine/plant.h"
#include "engine/node/node.h"
#include "engine/genome.h"

namespace botany {

bool can_merge(const Node& /*parent*/, const Node& /*child*/,
               const Genome& /*g*/, const CompressionParams& /*params*/) {
    return false; // Filled in by Task 2.
}

CompressionResult compress_plant(Plant& /*plant*/, const CompressionParams& /*params*/) {
    CompressionResult r;
    return r; // Filled in by Tasks 3-4.
}

} // namespace botany
```

- [ ] **Step 1.6: Run test to verify it passes**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[compression]" 2>&1 | tail -5
```
Expected: 1 test case, 2 assertions pass.

- [ ] **Step 1.7: Full-suite regression**

```
./build/botany_tests 2>&1 | tail -3
```
Expected: all tests pass (count up by 1).

- [ ] **Step 1.8: Commit**

```bash
git add src/engine/compression.h src/engine/compression.cpp tests/test_compression.cpp CMakeLists.txt
git commit -m "$(cat <<'EOF'
compression: module skeleton with stub compress_plant

Introduces src/engine/compression.{h,cpp} with CompressionParams,
CompressionResult, and stub compress_plant() that returns a zero
result.  Adds tests/test_compression.cpp under CMake.  Later tasks
fill in the predicate and scan.

EOF
)"
```

---

## Task 2: `can_merge` predicate with one test per rule

**Files:**
- Modify: `src/engine/compression.cpp`
- Test: `tests/test_compression.cpp`

Goal: implement the 7-rule predicate, each gate tested in isolation. Hand-build small plants via `Plant::from_empty()` + `install_node()` so every rule's trigger can be set deterministically.

- [ ] **Step 2.1: Write failing tests — one per rule**

Append to `tests/test_compression.cpp`:
```cpp
#include <memory>
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"

// Helpers: build a 3-node chain seed → A → B so we can test predicates on the
// (A, B) pair.  Returns raw pointers to the three nodes; plant owns them.
struct Chain3 {
    std::unique_ptr<Plant> plant;
    Node* seed;
    Node* A;
    Node* B;
};

static Chain3 build_stem_chain(const Genome& g,
                               float radiusA, float radiusB,
                               const glm::vec3& offsetA, const glm::vec3& offsetB,
                               uint32_t ageA, uint32_t ageB) {
    Chain3 out;
    out.plant = Plant::from_empty(g);

    auto seed = std::make_unique<StemNode>(0, glm::vec3(0.0f), g.initial_radius);
    out.seed = seed.get();
    seed->offset = glm::vec3(0.0f);

    auto a = std::make_unique<StemNode>(1, glm::vec3(0.0f), radiusA);
    out.A = a.get();
    a->parent = out.seed;
    a->offset = offsetA;
    a->age = ageA;
    out.seed->children.push_back(out.A);

    auto b = std::make_unique<StemNode>(2, glm::vec3(0.0f), radiusB);
    out.B = b.get();
    b->parent = out.A;
    b->offset = offsetB;
    b->age = ageB;
    out.A->children.push_back(out.B);

    out.plant->install_node(std::move(seed));
    out.plant->install_node(std::move(a));
    out.plant->install_node(std::move(b));
    out.plant->set_next_id(3);
    return out;
}

TEST_CASE("can_merge: same-type + mature + similar radius + straight passes", "[compression][predicate]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g,
        /*rA=*/0.030f, /*rB=*/0.030f,
        /*offA=*/glm::vec3(0.0f, 0.2f, 0.0f),
        /*offB=*/glm::vec3(0.0f, 0.2f, 0.0f),
        /*ageA=*/g.internode_maturation_ticks + 10,
        /*ageB=*/g.internode_maturation_ticks + 10);
    CompressionParams params;
    REQUIRE(can_merge(*c.A, *c.B, g, params));
}

TEST_CASE("can_merge: seed (no grandparent) is rejected", "[compression][predicate]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.03f, 0.03f,
        glm::vec3(0, 0.2f, 0), glm::vec3(0, 0.2f, 0),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    CompressionParams params;
    // Try to merge seed + A (seed has no parent → must be refused)
    REQUIRE_FALSE(can_merge(*c.seed, *c.A, g, params));
}

TEST_CASE("can_merge: bent chain (angle too large) is rejected", "[compression][predicate]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.03f, 0.03f,
        glm::vec3(0.0f, 0.2f, 0.0f),
        glm::vec3(0.2f, 0.0f, 0.0f), // 90° bend
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    CompressionParams params;
    REQUIRE_FALSE(can_merge(*c.A, *c.B, g, params));
}

TEST_CASE("can_merge: radius ratio above threshold is rejected", "[compression][predicate]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g,
        /*rA=*/0.030f, /*rB=*/0.060f, // 2× mismatch
        glm::vec3(0, 0.2f, 0), glm::vec3(0, 0.2f, 0),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    CompressionParams params;
    REQUIRE_FALSE(can_merge(*c.A, *c.B, g, params));
}

TEST_CASE("can_merge: still-elongating (young) nodes are rejected", "[compression][predicate]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.03f, 0.03f,
        glm::vec3(0, 0.2f, 0), glm::vec3(0, 0.2f, 0),
        /*ageA=*/10, /*ageB=*/10);
    CompressionParams params;
    REQUIRE_FALSE(can_merge(*c.A, *c.B, g, params));
}

TEST_CASE("can_merge: combined length above threshold is rejected", "[compression][predicate]") {
    Genome g = default_genome();
    // Two segments each at max_internode_length → combined = 2× max; threshold auto is 2×.
    // Make each slightly over max/2 so the sum exceeds the cap.
    float L = g.max_internode_length * 1.1f;
    auto c = build_stem_chain(g, 0.03f, 0.03f,
        glm::vec3(0, L, 0), glm::vec3(0, L, 0),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    CompressionParams params;
    // Default resolves to 2 * max_internode_length; 2 * 1.1 * max > 2 * max.
    REQUIRE_FALSE(can_merge(*c.A, *c.B, g, params));
}

TEST_CASE("can_merge: parent with 2 structural children is rejected", "[compression][predicate]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.03f, 0.03f,
        glm::vec3(0, 0.2f, 0), glm::vec3(0, 0.2f, 0),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    // Add a SECOND structural child to A — another stem branch.
    auto branch = std::make_unique<StemNode>(99, glm::vec3(0.0f), 0.03f);
    branch->parent = c.A;
    branch->offset = glm::vec3(0.1f, 0.1f, 0.0f);
    branch->age = g.internode_maturation_ticks + 10;
    c.A->children.push_back(branch.get());
    c.plant->install_node(std::move(branch));

    CompressionParams params;
    REQUIRE_FALSE(can_merge(*c.A, *c.B, g, params));
}

TEST_CASE("can_merge: parent with leaf + bud + one structural child passes", "[compression][predicate]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.03f, 0.03f,
        glm::vec3(0, 0.2f, 0), glm::vec3(0, 0.2f, 0),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    // Attach a leaf and a dormant axillary bud to A — these must NOT count as
    // structural children (leaves + meristems don't block compression).
    auto leaf = std::make_unique<LeafNode>(50, glm::vec3(0.0f), 0.0f);
    leaf->parent = c.A;
    leaf->offset = glm::vec3(0.1f, 0.0f, 0.0f);
    c.A->children.push_back(leaf.get());
    c.plant->install_node(std::move(leaf));

    auto bud = std::make_unique<ApicalNode>(51, glm::vec3(0.0f), 0.01f);
    bud->parent = c.A;
    bud->offset = glm::vec3(-0.1f, 0.0f, 0.0f);
    if (auto* a = bud->as_apical()) a->active = false;
    c.A->children.push_back(bud.get());
    c.plant->install_node(std::move(bud));

    CompressionParams params;
    REQUIRE(can_merge(*c.A, *c.B, g, params));
}

TEST_CASE("can_merge: stem cannot absorb root-type child", "[compression][predicate]") {
    Genome g = default_genome();
    auto plant = Plant::from_empty(g);
    auto seed = std::make_unique<StemNode>(0, glm::vec3(0.0f), g.initial_radius);
    auto stem = std::make_unique<StemNode>(1, glm::vec3(0.0f), 0.03f);
    stem->parent = seed.get();
    stem->offset = glm::vec3(0.0f, 0.2f, 0.0f);
    stem->age = g.internode_maturation_ticks + 10;
    seed->children.push_back(stem.get());
    auto root = std::make_unique<RootNode>(2, glm::vec3(0.0f), 0.03f);
    root->parent = stem.get();
    root->offset = glm::vec3(0.0f, 0.2f, 0.0f);
    root->age = g.internode_maturation_ticks + 10;
    stem->children.push_back(root.get());

    Node* stem_ptr = stem.get();
    Node* root_ptr = root.get();
    plant->install_node(std::move(seed));
    plant->install_node(std::move(stem));
    plant->install_node(std::move(root));
    plant->set_next_id(3);

    CompressionParams params;
    REQUIRE_FALSE(can_merge(*stem_ptr, *root_ptr, g, params));
}
```

- [ ] **Step 2.2: Run tests to verify they fail**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[compression][predicate]" 2>&1 | tail -5
```
Expected: most fail because `can_merge` currently returns false (the one that expects true fails immediately).

- [ ] **Step 2.3: Implement the predicate**

Replace the `can_merge` stub in `src/engine/compression.cpp` with:
```cpp
#include "engine/node/node.h"
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include <cmath>
#include <glm/geometric.hpp>

namespace botany {
namespace {

float effective_max_combined_length(const Genome& g, const CompressionParams& p) {
    return p.max_combined_length > 0.0f ? p.max_combined_length
                                        : 2.0f * g.max_internode_length;
}

// Count children of N that participate in the structural chain.  Leaves,
// apicals, and root_apicals are not structural — they decorate an internode
// but do not branch the trunk/root chain.
uint32_t count_structural_children(const Node& n) {
    uint32_t k = 0;
    for (const Node* c : n.children) {
        if (c->type == NodeType::STEM || c->type == NodeType::ROOT) ++k;
    }
    return k;
}

bool is_conduit_type(NodeType t) {
    return t == NodeType::STEM || t == NodeType::ROOT;
}

} // namespace

bool can_merge(const Node& parent, const Node& child,
               const Genome& g, const CompressionParams& params) {
    // Rule 1: same type, and must be a conduit type.
    if (parent.type != child.type) return false;
    if (!is_conduit_type(parent.type)) return false;

    // Rule 2: parent must have a grandparent (parent != seed).
    if (!parent.parent) return false;

    // Rule 3: parent has exactly ONE structural child and it's this child.
    if (count_structural_children(parent) != 1) return false;
    if (child.parent != &parent) return false;

    // Rule 4: angle between parent.offset and child.offset is small.
    float lenP = glm::length(parent.offset);
    float lenC = glm::length(child.offset);
    if (lenP < 1e-6f || lenC < 1e-6f) return false;
    float cos_gate = std::cos(params.max_angle_rad);
    float cos_pair = glm::dot(parent.offset, child.offset) / (lenP * lenC);
    if (cos_pair < cos_gate) return false;

    // Rule 5: radius ratio within tolerance.
    float rmax = std::max(parent.radius, child.radius);
    if (rmax < 1e-6f) return false;
    if (std::fabs(parent.radius - child.radius) / rmax > params.max_radius_ratio) return false;

    // Rule 6: both past elongation maturation.
    if (parent.age < g.internode_maturation_ticks) return false;
    if (child.age  < g.internode_maturation_ticks) return false;

    // Rule 7: combined length below threshold.
    if (lenP + lenC >= effective_max_combined_length(g, params)) return false;

    return true;
}

} // namespace botany
```

- [ ] **Step 2.4: Run tests to verify they pass**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[compression][predicate]" 2>&1 | tail -5
```
Expected: 9 test cases pass.

- [ ] **Step 2.5: Full-suite regression**

```
./build/botany_tests 2>&1 | tail -3
```
Expected: all tests pass.

- [ ] **Step 2.6: Commit**

```bash
git add src/engine/compression.cpp tests/test_compression.cpp
git commit -m "$(cat <<'EOF'
compression: can_merge predicate with per-rule tests

Seven-rule gate: same type, not seed, single structural child,
angle, radius ratio, maturation, combined length.  Leaves and
dormant meristems don't count as structural children so they
can ride along with a merge.

EOF
)"
```

---

## Task 3: `merge_pair()` execution — geometry, chemicals, reparent, bias, clamp

**Files:**
- Modify: `src/engine/compression.cpp`
- Test: `tests/test_compression.cpp`

Implement the per-pair merge. Expose a test-visible entry point `merge_pair()` for unit tests. The real `compress_plant` will call it from its scan in Task 4.

- [ ] **Step 3.1: Write failing tests**

Append to `tests/test_compression.cpp`:
```cpp
#include "engine/compression.h"

// Forward-declaration of the test-visible helper.  compress_plant calls this
// internally but the tests exercise it directly for precise assertions.
namespace botany {
void merge_pair(Plant& plant, Node& parent, Node& child, CompressionResult& result);
}

TEST_CASE("merge_pair absorbs C's geometry into P", "[compression][merge]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.04f, 0.04f,
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);

    CompressionResult result;
    merge_pair(*c.plant, *c.A, *c.B, result);
    c.plant->flush_removals();

    // Post-merge A's offset is the grandparent→C vector.
    REQUIRE(c.A->offset == glm::vec3(0.0f, 2.0f, 0.0f));
    // Volume-preserving radius of two equal-radius, equal-length segments is the same radius.
    REQUIRE(c.A->radius == Approx(0.04f).epsilon(1e-4));
    REQUIRE(c.plant->node_count() == 2u); // seed + A; C removed
    REQUIRE(result.merges_performed == 1u);
}

TEST_CASE("merge_pair uses volume-preserving radius on unequal r/L", "[compression][merge]") {
    Genome g = default_genome();
    // r1=0.05 L1=1.0, r2=0.03 L2=1.0 → r_merged = sqrt((0.0025 + 0.0009) / 2) = sqrt(0.0017)
    auto c = build_stem_chain(g, 0.05f, 0.03f,
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);

    // Temporarily loosen the radius ratio gate so the test can proceed.
    // (can_merge defaults reject |0.05-0.03|/0.05 = 0.4 > 0.20.)
    // merge_pair doesn't re-check; it's the caller's responsibility.
    CompressionResult result;
    merge_pair(*c.plant, *c.A, *c.B, result);
    c.plant->flush_removals();

    float expected = std::sqrt((0.05f*0.05f + 0.03f*0.03f) / 2.0f);
    REQUIRE(c.A->radius == Approx(expected).epsilon(1e-4));
}

TEST_CASE("merge_pair sums local chemicals", "[compression][merge]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.04f, 0.04f,
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    c.A->local().chemical(ChemicalID::Sugar) = 0.003f;
    c.B->local().chemical(ChemicalID::Sugar) = 0.005f;
    c.A->local().chemical(ChemicalID::Water) = 2.0f;
    c.B->local().chemical(ChemicalID::Water) = 3.0f;

    CompressionResult result;
    merge_pair(*c.plant, *c.A, *c.B, result);
    c.plant->flush_removals();

    // Caps on a 0.04 dm radius, L=2 dm stem are large enough to not clamp these values.
    REQUIRE(c.A->local().chemical(ChemicalID::Sugar) == Approx(0.008f).epsilon(1e-5));
    REQUIRE(c.A->local().chemical(ChemicalID::Water) == Approx(5.0f).epsilon(1e-5));
    REQUIRE(result.delta_sugar == 0.0f); // no clamping
    REQUIRE(result.delta_water == 0.0f);
}

TEST_CASE("merge_pair sums conduit pool chemicals", "[compression][merge]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.04f, 0.04f,
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    c.A->phloem()->chemical(ChemicalID::Sugar) = 0.001f;
    c.B->phloem()->chemical(ChemicalID::Sugar) = 0.002f;
    c.A->xylem()->chemical(ChemicalID::Water) = 0.5f;
    c.B->xylem()->chemical(ChemicalID::Water) = 0.7f;

    CompressionResult result;
    merge_pair(*c.plant, *c.A, *c.B, result);
    c.plant->flush_removals();

    REQUIRE(c.A->phloem()->chemical(ChemicalID::Sugar) == Approx(0.003f).epsilon(1e-5));
    REQUIRE(c.A->xylem()->chemical(ChemicalID::Water)  == Approx(1.2f).epsilon(1e-5));
}

TEST_CASE("merge_pair reparents children of C to P with same offsets", "[compression][merge]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.04f, 0.04f,
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    // Attach a leaf and a dormant bud to C — both should migrate to A.
    auto leaf = std::make_unique<LeafNode>(100, glm::vec3(0.0f), 0.0f);
    leaf->parent = c.B;
    leaf->offset = glm::vec3(0.1f, 0.0f, 0.0f);
    Node* leaf_ptr = leaf.get();
    c.B->children.push_back(leaf_ptr);
    c.plant->install_node(std::move(leaf));

    auto bud = std::make_unique<ApicalNode>(101, glm::vec3(0.0f), 0.01f);
    bud->parent = c.B;
    bud->offset = glm::vec3(-0.1f, 0.0f, 0.0f);
    Node* bud_ptr = bud.get();
    if (auto* a = bud_ptr->as_apical()) a->active = false;
    c.B->children.push_back(bud_ptr);
    c.plant->install_node(std::move(bud));

    CompressionResult result;
    merge_pair(*c.plant, *c.A, *c.B, result);
    c.plant->flush_removals();

    // Both leaf and bud now have A as parent; their offsets unchanged.
    REQUIRE(leaf_ptr->parent == c.A);
    REQUIRE(leaf_ptr->offset == glm::vec3(0.1f, 0.0f, 0.0f));
    REQUIRE(bud_ptr->parent == c.A);
    REQUIRE(bud_ptr->offset == glm::vec3(-0.1f, 0.0f, 0.0f));

    // A's children now include leaf + bud (and no B).
    bool found_leaf = false, found_bud = false;
    for (Node* cc : c.A->children) {
        if (cc == leaf_ptr) found_leaf = true;
        if (cc == bud_ptr)  found_bud = true;
    }
    REQUIRE(found_leaf);
    REQUIRE(found_bud);
    REQUIRE(c.A->children.size() == 2u); // leaf + bud; B removed
}

TEST_CASE("merge_pair carries auxin_flow_bias from B→child into A→child", "[compression][merge]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.04f, 0.04f,
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    // Attach a grandchild to B so bias transfer has somewhere to land.
    auto gc = std::make_unique<StemNode>(200, glm::vec3(0.0f), 0.04f);
    gc->parent = c.B;
    gc->offset = glm::vec3(0.0f, 1.0f, 0.0f);
    gc->age    = g.internode_maturation_ticks + 10;
    Node* gc_ptr = gc.get();
    c.B->children.push_back(gc_ptr);
    c.plant->install_node(std::move(gc));

    c.B->auxin_flow_bias[gc_ptr] = 0.73f;
    c.A->auxin_flow_bias[c.B]    = 0.44f; // will be destroyed

    CompressionResult result;
    merge_pair(*c.plant, *c.A, *c.B, result);
    c.plant->flush_removals();

    auto it = c.A->auxin_flow_bias.find(gc_ptr);
    REQUIRE(it != c.A->auxin_flow_bias.end());
    REQUIRE(it->second == 0.73f);
    // A→B entry is gone.
    REQUIRE(c.A->auxin_flow_bias.find(c.B) == c.A->auxin_flow_bias.end());
}
```

- [ ] **Step 3.2: Run tests to verify they fail**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5
```
Expected: link errors — `merge_pair` undeclared.

- [ ] **Step 3.3: Implement `merge_pair`**

In `src/engine/compression.cpp` add includes:
```cpp
#include "engine/sugar.h"              // sugar_cap, water_cap
#include "engine/vascular_sub_stepped.h" // phloem_capacity, xylem_capacity
#include <algorithm>
```

Define `merge_pair` (outside the anonymous namespace — the tests need its symbol):

```cpp
void merge_pair(Plant& plant, Node& parent, Node& child, CompressionResult& result) {
    const Genome& g = plant.genome();

    const float lenP = glm::length(parent.offset);
    const float lenC = glm::length(child.offset);
    const float total_len = lenP + lenC;

    // Volume-preserving radius: sqrt((rP²·LP + rC²·LC) / (LP + LC))
    const float r_merged = (total_len > 1e-8f)
        ? std::sqrt((parent.radius * parent.radius * lenP +
                     child.radius  * child.radius  * lenC) / total_len)
        : parent.radius;

    // Extend P's geometry so it now occupies C's old world position.
    const glm::vec3 new_offset      = parent.offset      + child.offset;
    const glm::vec3 new_rest_offset = parent.rest_offset + child.rest_offset;

    // Sum local chemicals.
    for (const auto& kv : child.local().chemicals) {
        parent.local().chemical(kv.first) += kv.second;
    }

    // Sum conduit pool chemicals.  Both are STEM or both ROOT (can_merge gate).
    if (auto* pp = parent.phloem()) {
        if (auto* cp = child.phloem()) {
            for (const auto& kv : cp->chemicals) {
                pp->chemical(kv.first) += kv.second;
            }
        }
    }
    if (auto* px = parent.xylem()) {
        if (auto* cx = child.xylem()) {
            for (const auto& kv : cx->chemicals) {
                px->chemical(kv.first) += kv.second;
            }
        }
    }

    // Reparent C's children onto P.  Their world-space positions stay put
    // (P now sits where C did), so their offsets (= child_position - parent_position)
    // are unchanged.
    for (Node* gc : child.children) {
        gc->parent = &parent;
        parent.children.push_back(gc);
    }
    child.children.clear(); // prevent die() cascade

    // Canalization bias transfer:
    // 1. Every C→grandchild entry migrates into parent.auxin_flow_bias.
    for (const auto& kv : child.auxin_flow_bias) {
        parent.auxin_flow_bias[kv.first] = kv.second;
    }
    // 2. The parent→child edge no longer exists — drop its entry.
    parent.auxin_flow_bias.erase(&child);

    // Apply new geometry BEFORE clamping (caps depend on radius + length).
    parent.offset      = new_offset;
    parent.rest_offset = new_rest_offset;
    parent.radius      = r_merged;

    // Clamp local + pool chemicals to new caps; log clamp losses.
    auto clamp_local = [&](ChemicalID id, float cap, float& running_delta) {
        float& v = parent.local().chemical(id);
        if (v > cap) {
            running_delta -= (v - cap);
            v = cap;
        }
    };
    const float sc = sugar_cap(parent, g);
    const float wc = water_cap(parent, g);
    clamp_local(ChemicalID::Sugar, sc, result.delta_sugar);
    clamp_local(ChemicalID::Water, wc, result.delta_water);
    // Auxin + cytokinin: no explicit per-node cap in the current model;
    // record nothing but keep the field wired for future clamps.
    (void)result.delta_auxin;
    (void)result.delta_cytokinin;

    // Pool caps.
    if (auto* pp = parent.phloem()) {
        const float pc = phloem_capacity(parent, g);
        float& v = pp->chemical(ChemicalID::Sugar);
        if (v > pc) { result.delta_sugar -= (v - pc); v = pc; }
    }
    if (auto* px = parent.xylem()) {
        const float xc = xylem_capacity(parent, g);
        float& vw = px->chemical(ChemicalID::Water);
        if (vw > xc) { result.delta_water -= (vw - xc); vw = xc; }
        // Cytokinin in xylem: no independent cap; Jacobi drives concentration.
    }

    // Queue child for deferred removal.
    plant.queue_removal(&child);
    result.merges_performed++;
}
```

`Approx` lives in `catch2/matchers/catch_matchers_floating_point.hpp` in Catch2 v3, but older tests use `Catch::Approx` from `catch_approx.hpp`. Check the existing test files for which include + namespace the codebase uses and match that — add the appropriate include and `using` line at the top of `test_compression.cpp`.

- [ ] **Step 3.4: Ensure `Approx` is visible in the test file**

Look at another test file (e.g. `tests/test_plant_snapshot.cpp`) for how it imports `Approx`. Match that style. If the codebase doesn't already include Approx, add:
```cpp
#include <catch2/catch_approx.hpp>
using Catch::Approx;
```
at the top of `tests/test_compression.cpp`.

- [ ] **Step 3.5: Run tests to verify they pass**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[compression][merge]" 2>&1 | tail -5
```
Expected: 6 merge tests pass.

- [ ] **Step 3.6: Full-suite regression**

```
./build/botany_tests 2>&1 | tail -3
```
Expected: all tests pass.

- [ ] **Step 3.7: Commit**

```bash
git add src/engine/compression.cpp tests/test_compression.cpp
git commit -m "$(cat <<'EOF'
compression: merge_pair execution — geometry, chemicals, bias, clamp

Implements the per-pair merge: volume-preserving radius, chemical
summation on local + conduit pools, child reparenting (offsets
unchanged since merged-P occupies C's old position), canalization
bias transfer (C→grandchild entries migrate to parent, parent→C
entry discarded), cap clamping with loss logged into
CompressionResult.delta_* fields.

EOF
)"
```

---

## Task 4: `compress_plant` scan + multi-pass loop + idempotence test

**Files:**
- Modify: `src/engine/compression.cpp`
- Test: `tests/test_compression.cpp`

- [ ] **Step 4.1: Write failing tests**

Append:
```cpp
TEST_CASE("compress_plant merges a straight mature chain", "[compression][scan]") {
    // Build seed → A → B → C → D, all mature same radius, all straight.
    Genome g = default_genome();
    auto plant = Plant::from_empty(g);
    auto seed = std::make_unique<StemNode>(0, glm::vec3(0.0f), g.initial_radius);
    Node* seed_ptr = seed.get();
    plant->install_node(std::move(seed));

    auto add_child = [&](Node* parent, uint32_t id, const glm::vec3& offset) -> Node* {
        auto n = std::make_unique<StemNode>(id, glm::vec3(0.0f), 0.04f);
        n->parent = parent;
        n->offset = offset;
        n->age    = g.internode_maturation_ticks + 50;
        Node* raw = n.get();
        parent->children.push_back(raw);
        plant->install_node(std::move(n));
        return raw;
    };
    Node* A = add_child(seed_ptr, 1, glm::vec3(0.0f, 0.1f, 0.0f));
    Node* B = add_child(A,        2, glm::vec3(0.0f, 0.1f, 0.0f));
    Node* C = add_child(B,        3, glm::vec3(0.0f, 0.1f, 0.0f));
    Node* D = add_child(C,        4, glm::vec3(0.0f, 0.1f, 0.0f));
    (void)B; (void)D;
    plant->set_next_id(5);
    REQUIRE(plant->node_count() == 5u); // seed + 4

    CompressionParams params;
    CompressionResult r = compress_plant(*plant, params);
    REQUIRE(r.merges_performed > 0u);
    REQUIRE(plant->node_count() < 5u);
}

TEST_CASE("compress_plant is idempotent — second run merges nothing", "[compression][scan]") {
    Genome g = default_genome();
    auto plant = Plant::from_empty(g);
    auto seed = std::make_unique<StemNode>(0, glm::vec3(0.0f), g.initial_radius);
    Node* seed_ptr = seed.get();
    plant->install_node(std::move(seed));
    auto add = [&](Node* p, uint32_t id) {
        auto n = std::make_unique<StemNode>(id, glm::vec3(0.0f), 0.04f);
        n->parent = p; n->offset = glm::vec3(0.0f, 0.1f, 0.0f);
        n->age = g.internode_maturation_ticks + 50;
        Node* raw = n.get();
        p->children.push_back(raw);
        plant->install_node(std::move(n));
        return raw;
    };
    Node* A = add(seed_ptr, 1);
    Node* B = add(A, 2);
    (void)B;
    plant->set_next_id(3);

    CompressionParams params;
    CompressionResult r1 = compress_plant(*plant, params);
    CompressionResult r2 = compress_plant(*plant, params);
    REQUIRE(r1.merges_performed >= 1u);
    REQUIRE(r2.merges_performed == 0u);
}

TEST_CASE("compress_plant respects max_passes cap", "[compression][scan]") {
    Genome g = default_genome();
    auto plant = Plant::from_empty(g);
    auto seed = std::make_unique<StemNode>(0, glm::vec3(0.0f), g.initial_radius);
    Node* seed_ptr = seed.get();
    plant->install_node(std::move(seed));
    Node* last = seed_ptr;
    for (uint32_t id = 1; id < 20; ++id) {
        auto n = std::make_unique<StemNode>(id, glm::vec3(0.0f), 0.04f);
        n->parent = last; n->offset = glm::vec3(0.0f, 0.05f, 0.0f);
        n->age = g.internode_maturation_ticks + 50;
        Node* raw = n.get();
        last->children.push_back(raw);
        plant->install_node(std::move(n));
        last = raw;
    }
    plant->set_next_id(20);

    CompressionParams params;
    params.max_passes = 1;
    CompressionResult r = compress_plant(*plant, params);
    REQUIRE(r.passes_run == 1u);
}
```

- [ ] **Step 4.2: Run tests to verify they fail**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[compression][scan]" 2>&1 | tail -5
```
Expected: `merges_performed == 0` because `compress_plant` is still a stub.

- [ ] **Step 4.3: Implement `compress_plant`**

Replace the stub in `src/engine/compression.cpp`:
```cpp
namespace {

// Single pass.  For each candidate (parent, child) pair where parent has
// exactly one structural child, invoke merge_pair if can_merge accepts it.
// Collects pointers first to avoid mutating the tree mid-traversal.
uint32_t compress_plant_single_pass(Plant& plant,
                                    const Genome& g,
                                    const CompressionParams& params,
                                    CompressionResult& result) {
    // Snapshot candidate pairs before mutation.
    std::vector<std::pair<Node*, Node*>> candidates;
    plant.for_each_node_mut([&](Node& parent) {
        if (!is_conduit_type(parent.type)) return;
        if (count_structural_children(parent) != 1) return;
        for (Node* child : parent.children) {
            if (child->type == parent.type && is_conduit_type(child->type)) {
                candidates.emplace_back(&parent, child);
                break;
            }
        }
    });

    uint32_t merges_this_pass = 0;
    // Track nodes that are already consumed this pass so we don't try to
    // merge a node that just got absorbed.
    std::unordered_set<Node*> consumed;
    for (auto& pair : candidates) {
        Node* p = pair.first;
        Node* c = pair.second;
        if (consumed.count(p) || consumed.count(c)) continue;
        if (!can_merge(*p, *c, g, params)) continue;
        merge_pair(plant, *p, *c, result);
        consumed.insert(c);
        merges_this_pass++;
    }
    return merges_this_pass;
}

} // namespace

CompressionResult compress_plant(Plant& plant, const CompressionParams& params) {
    CompressionResult r;
    const Genome& g = plant.genome();
    for (uint32_t pass = 0; pass < params.max_passes; ++pass) {
        uint32_t merged = compress_plant_single_pass(plant, g, params, r);
        r.passes_run++;
        plant.flush_removals();
        if (merged == 0) break;
    }
    return r;
}
```

Add `#include <unordered_set>` and `#include <utility>` to the file's includes.

- [ ] **Step 4.4: Run tests to verify they pass**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[compression]" 2>&1 | tail -5
```
Expected: all compression tests pass.

- [ ] **Step 4.5: Full-suite regression**

```
./build/botany_tests 2>&1 | tail -3
```
Expected: all tests pass.

- [ ] **Step 4.6: Commit**

```bash
git add src/engine/compression.cpp tests/test_compression.cpp
git commit -m "$(cat <<'EOF'
compression: compress_plant scan + multi-pass loop

Scans the tree for mergeable (parent, child) pairs, invokes
merge_pair on accepted candidates, and loops up to max_passes to
catch newly-adjacent pairs created by earlier merges.  Calls
flush_removals between passes so pointer state stays clean.

EOF
)"
```

---

## Task 5: Engine integration — autocompress controls + tick() hook

**Files:**
- Modify: `src/engine/engine.h`
- Modify: `src/engine/engine.cpp`
- Test: `tests/test_compression.cpp`

- [ ] **Step 5.1: Write failing tests**

Append:
```cpp
#include "engine/engine.h"

TEST_CASE("Engine::trigger_compression runs immediately", "[compression][engine]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));
    for (int i = 0; i < 200; ++i) engine.tick(); // grow enough internodes
    uint32_t before = engine.get_plant(0).node_count();

    engine.trigger_compression();
    uint32_t after = engine.get_plant(0).node_count();
    // Either the plant has mergeable pairs (after < before) or it doesn't
    // (after == before).  Both are acceptable; we're testing the mechanism.
    REQUIRE(after <= before);
    REQUIRE(engine.last_compression().passes_run >= 1u);
}

TEST_CASE("Engine autocompress fires on interval", "[compression][engine]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));
    engine.enable_autocompress(true);
    engine.set_compression_interval(5);

    // Tick 20 times — interval 5 should trigger 4 compressions (ticks 5,10,15,20).
    for (int i = 0; i < 20; ++i) engine.tick();
    REQUIRE(engine.last_compression().passes_run >= 1u);
}

TEST_CASE("Engine autocompress disabled does not run compression", "[compression][engine]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));
    engine.enable_autocompress(false);
    engine.set_compression_interval(5);
    for (int i = 0; i < 20; ++i) engine.tick();
    REQUIRE(engine.last_compression().passes_run == 0u);
}
```

- [ ] **Step 5.2: Add declarations to `engine.h`**

At the top of `engine.h` add:
```cpp
#include "engine/compression.h"
```

In the `public:` section of `Engine`:
```cpp
    // --- Compression controls ---
    void enable_autocompress(bool enabled);
    void set_compression_interval(uint32_t ticks);
    void set_compression_params(const CompressionParams& params);
    CompressionResult trigger_compression();
    const CompressionResult& last_compression() const;
```

In the `private:` section:
```cpp
    bool              compression_enabled_  = false;
    uint32_t          compression_interval_ = 1000;
    CompressionParams compression_params_;
    CompressionResult last_compression_;
```

- [ ] **Step 5.3: Implement in `engine.cpp`**

At the bottom of `engine.cpp` (still inside `namespace botany { ... }`):
```cpp
void Engine::enable_autocompress(bool enabled) {
    compression_enabled_ = enabled;
}

void Engine::set_compression_interval(uint32_t ticks) {
    compression_interval_ = ticks == 0 ? 1 : ticks;
}

void Engine::set_compression_params(const CompressionParams& params) {
    compression_params_ = params;
}

CompressionResult Engine::trigger_compression() {
    if (plants_.empty()) {
        last_compression_ = CompressionResult{};
        return last_compression_;
    }
    last_compression_ = compress_plant(*plants_[0], compression_params_);
    return last_compression_;
}

const CompressionResult& Engine::last_compression() const {
    return last_compression_;
}
```

Find the existing `Engine::tick()` body. After the per-plant tick loop, BEFORE `tick_++`, add:
```cpp
    // Autocompress: runs at tick-interval boundaries, between full ticks.
    if (compression_enabled_
        && tick_ > 0
        && (tick_ % compression_interval_) == 0
        && !plants_.empty()) {
        last_compression_ = compress_plant(*plants_[0], compression_params_);
    }
```

If you aren't sure where `tick_++` lives in the current `Engine::tick()`, search for `tick_++` or `++tick_` in `engine.cpp` — it's right before the function returns.

- [ ] **Step 5.4: Run tests to verify they pass**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[compression][engine]" 2>&1 | tail -5
```
Expected: 3 engine tests pass.

- [ ] **Step 5.5: Full-suite regression**

```
./build/botany_tests 2>&1 | tail -3
```
Expected: all tests pass.

- [ ] **Step 5.6: Commit**

```bash
git add src/engine/engine.h src/engine/engine.cpp tests/test_compression.cpp
git commit -m "$(cat <<'EOF'
engine: autocompress hook + manual trigger + last-result accessor

Engine gains enable_autocompress(bool), set_compression_interval,
set_compression_params, trigger_compression, last_compression.
Engine::tick() checks at end: if autocompress is enabled and the
elapsed tick count hits a multiple of the interval, run
compress_plant on the primary plant.  The compression always runs
strictly between full ticks (after per-plant tick completes,
before tick counter advances).

EOF
)"
```

---

## Task 6: Integration tests — leaf survival + save/load after compression

**Files:**
- Test: `tests/test_compression.cpp`

End-to-end tests with a real grown plant.

- [ ] **Step 6.1: Write failing tests**

Append:
```cpp
#include "serialization/plant_snapshot.h"
#include <filesystem>
#include <optional>

TEST_CASE("compression preserves leaves world-position after next tick", "[compression][integration]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));
    for (int i = 0; i < 300; ++i) engine.tick();
    Plant& plant = engine.get_plant_mut(0);

    // Snapshot every LEAF's world position pre-compression.
    std::unordered_map<uint32_t, glm::vec3> leaf_pos_before;
    plant.for_each_node([&](const Node& n) {
        if (n.type == NodeType::LEAF) leaf_pos_before[n.id] = n.position;
    });
    REQUIRE(!leaf_pos_before.empty());

    CompressionResult r = engine.trigger_compression();
    // One more tick so sync_world_position runs on the modified chain.
    engine.tick();

    uint32_t checked = 0;
    plant.for_each_node([&](const Node& n) {
        if (n.type != NodeType::LEAF) return;
        auto it = leaf_pos_before.find(n.id);
        if (it == leaf_pos_before.end()) return; // leaf that was created post-compression
        float d = glm::distance(n.position, it->second);
        // Allow one tick of growth movement (leaves drift as stems elongate).
        REQUIRE(d < 0.01f);
        checked++;
    });
    // Only require that we actually checked something if compression ran.
    if (r.merges_performed > 0u) REQUIRE(checked > 0u);
}

TEST_CASE("compress → save → load → continue ticks cleanly", "[compression][integration]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));
    for (int i = 0; i < 300; ++i) engine.tick();
    engine.trigger_compression();
    uint32_t nodes_after_compress = engine.get_plant(0).node_count();

    auto tmp = std::filesystem::temp_directory_path() / "botany_compress_snapshot_test";
    std::filesystem::remove_all(tmp);
    SaveResult sr = save_plant_snapshot(engine.get_plant(0), engine.get_tick(), tmp.string());
    REQUIRE(sr.ok);

    LoadedPlant lp = load_plant_snapshot(sr.path, std::nullopt);
    REQUIRE(lp.plant->node_count() == nodes_after_compress);

    Engine engineB;
    engineB.adopt_plant(std::move(lp.plant));
    engineB.set_tick(static_cast<uint32_t>(lp.engine_tick));
    for (int i = 0; i < 20; ++i) engineB.tick();
    REQUIRE(engineB.get_plant(0).node_count() > 0u);

    std::filesystem::remove_all(tmp);
}
```

- [ ] **Step 6.2: Run tests to verify they pass**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[compression][integration]" 2>&1 | tail -5
```
Expected: both integration tests pass.

- [ ] **Step 6.3: Full-suite regression**

```
./build/botany_tests 2>&1 | tail -3
```
Expected: all tests pass.

- [ ] **Step 6.4: Commit**

```bash
git add tests/test_compression.cpp
git commit -m "compression: integration tests for leaf survival + save/load round-trip"
```

---

## Task 7: Realtime viewer — Compression panel in Controls

**Files:**
- Modify: `src/app_realtime.cpp`

No unit tests — the viewer panel is exercised manually. Keep the changes tight.

- [ ] **Step 7.1: Locate the Controls window**

In `src/app_realtime.cpp`, find the `ImGui::Begin("Controls", ...)` line (there's only one). The window already contains a few `CollapsingHeader` sections (e.g. "Overlays", "Info"). Add the new Compression section alongside them.

Near the top of `main()` (among other viewer state like `save_toast_msg`), add:
```cpp
bool compression_autocompress = false;
int  compression_interval     = 1000;
```

Right after `engine` is constructed, seed the engine state to match:
```cpp
engine.enable_autocompress(compression_autocompress);
engine.set_compression_interval(static_cast<uint32_t>(compression_interval));
```

- [ ] **Step 7.2: Add the Compression header inside the Controls window**

Inside the `ImGui::Begin("Controls", ...)` block, after the existing headers (find a natural position — e.g., after "Overlays" and before "Info"), add:
```cpp
if (ImGui::CollapsingHeader("Compression")) {
    if (ImGui::Checkbox("Auto-compress", &compression_autocompress)) {
        engine.enable_autocompress(compression_autocompress);
    }
    if (compression_autocompress) {
        if (ImGui::InputInt("Every N ticks", &compression_interval)) {
            if (compression_interval < 1) compression_interval = 1;
            engine.set_compression_interval(static_cast<uint32_t>(compression_interval));
        }
    } else {
        if (ImGui::Button("Compress Now")) {
            engine.trigger_compression();
        }
    }
    const CompressionResult& r = engine.last_compression();
    ImGui::Separator();
    ImGui::Text("Last run: %u merges / %u passes",
                r.merges_performed, r.passes_run);
    if (r.delta_sugar != 0.0f || r.delta_water != 0.0f) {
        ImGui::Text("Loss: sugar=%.4g g, water=%.4g ml",
                    r.delta_sugar, r.delta_water);
    } else if (r.passes_run > 0u) {
        ImGui::TextUnformatted("No cap-clamp loss.");
    }
}
```

- [ ] **Step 7.3: Build**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5
```
Expected: clean build. The realtime viewer is GUI-gated so you can't automate a deeper check here; launch it manually to sanity-check the new panel appears.

- [ ] **Step 7.4: Full-suite regression**

```
./build/botany_tests 2>&1 | tail -3
```
Expected: all tests pass.

- [ ] **Step 7.5: Commit**

```bash
git add src/app_realtime.cpp
git commit -m "app_realtime: Compression panel — checkbox, interval, button, stats"
```

---

## Task 8: Manual smoke + CLAUDE.md update + push

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 8.1: Manual smoke**

1. `./build/botany_realtime` — grow a plant for a couple thousand ticks.
2. Expand the Controls window → Compression header.
3. Click "Compress Now". Node count drops; last-run stats populate.
4. Toggle "Auto-compress" on; the interval input replaces the button.
5. Keep ticking; watch the stats update every N ticks.
6. Press Cmd/Ctrl+S to save. Quit. Relaunch with `--load-plant <path>`. Verify plant loads and continues.

If anything misbehaves (e.g., popping leaves, crash), root-cause before moving on.

- [ ] **Step 8.2: Update CLAUDE.md**

In `CLAUDE.md`, in the `### Apps` section near where the Snapshot format is described, add a new paragraph:
```
Compression: `src/engine/compression.{h,cpp}` collapses adjacent stem/root
nodes that share type, alignment, and size.  Runs between ticks — manually
from the Controls → Compression panel in `botany_realtime`, or automatically
every N ticks when auto-compress is enabled.  Lossy by design: cap-clamping
can trim trace sugar/water per merge (logged in `CompressionResult`).
```

- [ ] **Step 8.3: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: plant compression feature in CLAUDE.md"
```

- [ ] **Step 8.4: Push**

```
git push origin main
```
