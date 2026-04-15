# Auxin Growth Sensitivity Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add per-tissue-type auxin growth sensitivity using a saturating Michaelis-Menten multiplier that chains into existing growth rate computations.

**Architecture:** One shared helper function (`auxin_growth_factor`) in `helpers.h`, 10 new genome fields (5 tissue types x 2 params), one-line insertions in 5 existing growth functions, evolution bridge registration. Auxin shapes desire, sugar gates ability, physical caps unchanged.

**Tech Stack:** C++17, Catch2 (testing), CMake (build)

**Spec:** `docs/superpowers/specs/2026-04-15-auxin-growth-sensitivity-design.md`

---

### Task 1: Helper Function + Unit Tests

**Files:**
- Modify: `src/engine/node/meristems/helpers.h` (add function at bottom of `meristem_helpers` namespace)
- Create: `tests/test_auxin_sensitivity.cpp`
- Modify: `CMakeLists.txt:139-151` (add test file to `botany_tests`)

- [ ] **Step 1: Write the test file**

Create `tests/test_auxin_sensitivity.cpp`:

```cpp
#include <catch2/catch_test_macros.hpp>
#include "engine/node/meristems/helpers.h"
#include <cmath>

using namespace botany::meristem_helpers;

TEST_CASE("auxin_growth_factor: zero auxin returns 1.0", "[auxin_sensitivity]") {
    REQUIRE(auxin_growth_factor(0.0f, 0.5f, 0.2f) == 1.0f);
    REQUIRE(auxin_growth_factor(0.0f, -0.3f, 0.1f) == 1.0f);
}

TEST_CASE("auxin_growth_factor: positive max_boost promotes growth", "[auxin_sensitivity]") {
    float factor = auxin_growth_factor(0.2f, 0.5f, 0.2f);
    // At half-saturation, factor should be 1.0 + 0.5 * 0.5 = 1.25
    REQUIRE(factor > 1.0f);
    REQUIRE(factor < 1.5f);
    REQUIRE(std::abs(factor - 1.25f) < 0.01f);
}

TEST_CASE("auxin_growth_factor: negative max_boost inhibits growth", "[auxin_sensitivity]") {
    float factor = auxin_growth_factor(0.1f, -0.3f, 0.1f);
    // At half-saturation, factor should be 1.0 + (-0.3) * 0.5 = 0.85
    REQUIRE(factor < 1.0f);
    REQUIRE(factor > 0.7f);
    REQUIRE(std::abs(factor - 0.85f) < 0.01f);
}

TEST_CASE("auxin_growth_factor: high auxin asymptotes at 1.0 + max_boost", "[auxin_sensitivity]") {
    float factor = auxin_growth_factor(1000.0f, 0.5f, 0.2f);
    REQUIRE(std::abs(factor - 1.5f) < 0.01f);

    float factor_neg = auxin_growth_factor(1000.0f, -0.3f, 0.1f);
    REQUIRE(std::abs(factor_neg - 0.7f) < 0.01f);
}

TEST_CASE("auxin_growth_factor: zero max_boost always returns 1.0", "[auxin_sensitivity]") {
    REQUIRE(auxin_growth_factor(0.5f, 0.0f, 0.2f) == 1.0f);
    REQUIRE(auxin_growth_factor(100.0f, 0.0f, 0.1f) == 1.0f);
}

TEST_CASE("auxin_growth_factor: monotonically increases with auxin for positive boost", "[auxin_sensitivity]") {
    float f1 = auxin_growth_factor(0.05f, 0.5f, 0.2f);
    float f2 = auxin_growth_factor(0.1f, 0.5f, 0.2f);
    float f3 = auxin_growth_factor(0.5f, 0.5f, 0.2f);
    REQUIRE(f1 < f2);
    REQUIRE(f2 < f3);
}

TEST_CASE("auxin_growth_factor: monotonically decreases with auxin for negative boost", "[auxin_sensitivity]") {
    float f1 = auxin_growth_factor(0.05f, -0.3f, 0.1f);
    float f2 = auxin_growth_factor(0.1f, -0.3f, 0.1f);
    float f3 = auxin_growth_factor(0.5f, -0.3f, 0.1f);
    REQUIRE(f1 > f2);
    REQUIRE(f2 > f3);
}
```

- [ ] **Step 2: Add test file to CMakeLists.txt**

In `CMakeLists.txt`, add `tests/test_auxin_sensitivity.cpp` to the `botany_tests` target. Insert after line 150 (`tests/test_evolution.cpp`):

```
    tests/test_auxin_sensitivity.cpp
```

- [ ] **Step 3: Build and verify tests fail**

Run: `/usr/local/bin/cmake --build build`

Expected: Compile error — `auxin_growth_factor` not defined.

- [ ] **Step 4: Implement the helper function**

In `src/engine/node/meristems/helpers.h`, add at the bottom of the `meristem_helpers` namespace (before the closing `}`):

```cpp
// Saturating auxin growth multiplier (Michaelis-Menten).
// Returns 1.0 at zero auxin, asymptotes to 1.0 + max_boost.
// Positive max_boost = promotion, negative = inhibition.
inline float auxin_growth_factor(float auxin, float max_boost, float half_sat) {
    if (std::abs(max_boost) < 1e-8f) return 1.0f;
    float saturation = auxin / (auxin + std::max(half_sat, 1e-6f));
    return 1.0f + max_boost * saturation;
}
```

- [ ] **Step 5: Build and verify tests pass**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[auxin_sensitivity]"`

Expected: All 7 tests pass.

- [ ] **Step 6: Commit**

```bash
git add tests/test_auxin_sensitivity.cpp src/engine/node/meristems/helpers.h CMakeLists.txt
git commit -m "feat: add auxin_growth_factor helper with unit tests"
```

---

### Task 2: Genome Fields

**Files:**
- Modify: `src/engine/genome.h` (add 10 fields to struct + defaults)

- [ ] **Step 1: Add fields to Genome struct**

In `src/engine/genome.h`, add these 10 fields after the `auxin_bias` field (after line 21):

```cpp
    // Auxin growth sensitivity — saturating Michaelis-Menten per tissue type.
    // max_boost: signed ceiling (positive = promotes, negative = inhibits growth).
    // half_saturation: auxin level for half-max effect.
    float stem_auxin_max_boost;
    float stem_auxin_half_saturation;
    float root_auxin_max_boost;
    float root_auxin_half_saturation;
    float leaf_auxin_max_boost;
    float leaf_auxin_half_saturation;
    float apical_auxin_max_boost;
    float apical_auxin_half_saturation;
    float root_apical_auxin_max_boost;
    float root_apical_auxin_half_saturation;
```

- [ ] **Step 2: Add defaults in default_genome()**

In `default_genome()`, add after the `.auxin_bias = -0.1f,` line:

```cpp
        .stem_auxin_max_boost = 0.5f,            // auxin promotes stem elongation by up to 50%
        .stem_auxin_half_saturation = 0.2f,
        .root_auxin_max_boost = -0.3f,            // auxin inhibits root elongation by up to 30%
        .root_auxin_half_saturation = 0.1f,       // roots are very sensitive
        .leaf_auxin_max_boost = 0.3f,             // auxin promotes leaf expansion by up to 30%
        .leaf_auxin_half_saturation = 0.2f,
        .apical_auxin_max_boost = 0.2f,           // mild promotion of tip extension
        .apical_auxin_half_saturation = 0.3f,     // high half-sat — apicals sit in high auxin
        .root_apical_auxin_max_boost = -0.2f,     // mild inhibition of root tip extension
        .root_apical_auxin_half_saturation = 0.1f,
```

- [ ] **Step 3: Build to verify struct compiles**

Run: `/usr/local/bin/cmake --build build`

Expected: Clean build. All existing tests still pass.

- [ ] **Step 4: Commit**

```bash
git add src/engine/genome.h
git commit -m "feat: add per-tissue auxin growth sensitivity genome fields"
```

---

### Task 3: Wire Up Stem Elongation

**Files:**
- Modify: `src/engine/node/stem_node.cpp:39-66` (`elongate()` method)

- [ ] **Step 1: Add auxin multiplier to StemNode::elongate()**

In `src/engine/node/stem_node.cpp`, in the `elongate()` method, add the auxin boost line after the stress_inhibit line (after line 46) and chain it into effective_rate (modify line 47):

Current line 47:
```cpp
    float effective_rate = g.internode_elongation_rate * ga_boost * eth_inhibit * stress_inhibit;
```

Replace with:
```cpp
    float auxin_boost = meristem_helpers::auxin_growth_factor(
        chemical(ChemicalID::Auxin), g.stem_auxin_max_boost, g.stem_auxin_half_saturation);
    float effective_rate = g.internode_elongation_rate * ga_boost * eth_inhibit * stress_inhibit * auxin_boost;
```

Also add the include at the top of the file if not already present:
```cpp
#include "engine/node/meristems/helpers.h"
```

- [ ] **Step 2: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`

Expected: All tests pass. No existing behavior changes since default `stem_auxin_max_boost = 0.5` is a mild boost, not a breaking change.

- [ ] **Step 3: Commit**

```bash
git add src/engine/node/stem_node.cpp
git commit -m "feat: wire auxin growth sensitivity into stem elongation"
```

---

### Task 4: Wire Up Root Elongation

**Files:**
- Modify: `src/engine/node/root_node.cpp:37-60` (`elongate()` method)

- [ ] **Step 1: Add auxin multiplier to RootNode::elongate()**

In `src/engine/node/root_node.cpp`, in the `elongate()` method, add the auxin boost after the eth_inhibit line (after line 43) and chain it into effective_rate (modify line 44):

Current line 44:
```cpp
    float effective_rate = g.root_internode_elongation_rate * ga_boost * eth_inhibit;
```

Replace with:
```cpp
    float auxin_boost = meristem_helpers::auxin_growth_factor(
        chemical(ChemicalID::Auxin), g.root_auxin_max_boost, g.root_auxin_half_saturation);
    float effective_rate = g.root_internode_elongation_rate * ga_boost * eth_inhibit * auxin_boost;
```

Also add the include at the top of the file:
```cpp
#include "engine/node/meristems/helpers.h"
```

- [ ] **Step 2: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`

Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
git add src/engine/node/root_node.cpp
git commit -m "feat: wire auxin growth sensitivity into root elongation"
```

---

### Task 5: Wire Up Leaf Expansion

**Files:**
- Modify: `src/engine/node/tissues/leaf.cpp:111-135` (`grow_size()` method)

- [ ] **Step 1: Add auxin multiplier to LeafNode::grow_size()**

In `src/engine/node/tissues/leaf.cpp`, in the `grow_size()` method, add the auxin boost after the first line and modify `max_growth`:

Current lines 114-115:
```cpp
    float max_growth = g.leaf_growth_rate;
    float remaining = g.max_leaf_size - leaf_size;
```

Replace with:
```cpp
    float auxin_boost = meristem_helpers::auxin_growth_factor(
        chemical(ChemicalID::Auxin), g.leaf_auxin_max_boost, g.leaf_auxin_half_saturation);
    float max_growth = g.leaf_growth_rate * auxin_boost;
    float remaining = g.max_leaf_size - leaf_size;
```

The `helpers.h` include is already present via `meristem_types.h` or other meristem includes. If not, add:
```cpp
#include "engine/node/meristems/helpers.h"
```

- [ ] **Step 2: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`

Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
git add src/engine/node/tissues/leaf.cpp
git commit -m "feat: wire auxin growth sensitivity into leaf expansion"
```

---

### Task 6: Wire Up Shoot Apical Tip Extension

**Files:**
- Modify: `src/engine/node/tissues/apical.cpp:117-169` (`grow_tip()` method)

- [ ] **Step 1: Add auxin multiplier to ApicalNode::grow_tip()**

In `src/engine/node/tissues/apical.cpp`, in the `grow_tip()` method, modify the actual_rate computation near the end of the function.

Current line 166:
```cpp
    float actual_rate = g.growth_rate * gf;
```

Replace with:
```cpp
    float auxin_boost = meristem_helpers::auxin_growth_factor(
        chemical(ChemicalID::Auxin), g.apical_auxin_max_boost, g.apical_auxin_half_saturation);
    float actual_rate = g.growth_rate * gf * auxin_boost;
```

The `helpers.h` include is already present (line 2).

- [ ] **Step 2: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`

Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
git add src/engine/node/tissues/apical.cpp
git commit -m "feat: wire auxin growth sensitivity into shoot apical tip extension"
```

---

### Task 7: Wire Up Root Apical Tip Extension

**Files:**
- Modify: `src/engine/node/tissues/root_apical.cpp:49-60` (`grow_tip()` method)

- [ ] **Step 1: Add auxin multiplier to RootApicalNode::grow_tip()**

In `src/engine/node/tissues/root_apical.cpp`, in the `grow_tip()` method, modify the actual_rate computation:

Current line 57:
```cpp
    float actual_rate = g.root_growth_rate * gf;
```

Replace with:
```cpp
    float auxin_boost = meristem_helpers::auxin_growth_factor(
        chemical(ChemicalID::Auxin), g.root_apical_auxin_max_boost, g.root_apical_auxin_half_saturation);
    float actual_rate = g.root_growth_rate * gf * auxin_boost;
```

The `helpers.h` include is already present (line 2).

- [ ] **Step 2: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`

Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
git add src/engine/node/tissues/root_apical.cpp
git commit -m "feat: wire auxin growth sensitivity into root apical tip extension"
```

---

### Task 8: Evolution Bridge Registration

**Files:**
- Modify: `src/evolution/genome_bridge.cpp` (add 10 gene registrations + linkage group entries)

- [ ] **Step 1: Add gene registrations**

In `src/evolution/genome_bridge.cpp`, in `build_genome_template()`, add after the `auxin_bias` registration (after line 27):

```cpp
    reg(sg, "stem_auxin_max_boost",              g.stem_auxin_max_boost,              r, -1.0f, 2.0f, p);
    reg(sg, "stem_auxin_half_saturation",        g.stem_auxin_half_saturation,        r, 0.01f, 1.0f, p);
    reg(sg, "root_auxin_max_boost",              g.root_auxin_max_boost,              r, -1.0f, 2.0f, p);
    reg(sg, "root_auxin_half_saturation",        g.root_auxin_half_saturation,        r, 0.01f, 1.0f, p);
    reg(sg, "leaf_auxin_max_boost",              g.leaf_auxin_max_boost,              r, -1.0f, 2.0f, p);
    reg(sg, "leaf_auxin_half_saturation",        g.leaf_auxin_half_saturation,        r, 0.01f, 1.0f, p);
    reg(sg, "apical_auxin_max_boost",            g.apical_auxin_max_boost,            r, -1.0f, 2.0f, p);
    reg(sg, "apical_auxin_half_saturation",      g.apical_auxin_half_saturation,      r, 0.01f, 1.0f, p);
    reg(sg, "root_apical_auxin_max_boost",       g.root_apical_auxin_max_boost,       r, -1.0f, 2.0f, p);
    reg(sg, "root_apical_auxin_half_saturation", g.root_apical_auxin_half_saturation, r, 0.01f, 1.0f, p);
```

- [ ] **Step 2: Add to auxin linkage group**

In the same file, add the 10 new gene names to the `"auxin"` linkage group (around line 126-131). Add after `"auxin_bias"`:

```cpp
        "stem_auxin_max_boost", "stem_auxin_half_saturation",
        "root_auxin_max_boost", "root_auxin_half_saturation",
        "leaf_auxin_max_boost", "leaf_auxin_half_saturation",
        "apical_auxin_max_boost", "apical_auxin_half_saturation",
        "root_apical_auxin_max_boost", "root_apical_auxin_half_saturation"
```

- [ ] **Step 3: Add to from_structured()**

In `from_structured()`, add after the `g.auxin_bias` line (after line 210):

```cpp
    g.stem_auxin_max_boost              = sg.get("stem_auxin_max_boost");
    g.stem_auxin_half_saturation        = sg.get("stem_auxin_half_saturation");
    g.root_auxin_max_boost              = sg.get("root_auxin_max_boost");
    g.root_auxin_half_saturation        = sg.get("root_auxin_half_saturation");
    g.leaf_auxin_max_boost              = sg.get("leaf_auxin_max_boost");
    g.leaf_auxin_half_saturation        = sg.get("leaf_auxin_half_saturation");
    g.apical_auxin_max_boost            = sg.get("apical_auxin_max_boost");
    g.apical_auxin_half_saturation      = sg.get("apical_auxin_half_saturation");
    g.root_apical_auxin_max_boost       = sg.get("root_apical_auxin_max_boost");
    g.root_apical_auxin_half_saturation = sg.get("root_apical_auxin_half_saturation");
```

- [ ] **Step 4: Build and run all tests (including evolution tests)**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`

Expected: All tests pass, including `[evolution]` tests.

- [ ] **Step 5: Commit**

```bash
git add src/evolution/genome_bridge.cpp
git commit -m "feat: register auxin growth sensitivity genes in evolution bridge"
```

---

### Task 9: Update CLAUDE.md

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Update tuning parameters section**

In `CLAUDE.md`, in the "Tuning Parameters" section, add after the existing auxin parameters:

```
- `stem_auxin_max_boost` (0.5) - max elongation promotion from auxin (saturating)
- `stem_auxin_half_saturation` (0.2) - auxin level for half-max stem effect
- `root_auxin_max_boost` (-0.3) - max elongation inhibition from auxin (negative = inhibits)
- `root_auxin_half_saturation` (0.1) - auxin level for half-max root effect (roots very sensitive)
- `leaf_auxin_max_boost` (0.3) - max leaf expansion promotion from auxin
- `leaf_auxin_half_saturation` (0.2) - auxin level for half-max leaf effect
- `apical_auxin_max_boost` (0.2) - max tip extension promotion from auxin
- `apical_auxin_half_saturation` (0.3) - auxin level for half-max apical effect
- `root_apical_auxin_max_boost` (-0.2) - max root tip extension inhibition from auxin
- `root_apical_auxin_half_saturation` (0.1) - auxin level for half-max root apical effect
```

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: add auxin growth sensitivity params to CLAUDE.md"
```
