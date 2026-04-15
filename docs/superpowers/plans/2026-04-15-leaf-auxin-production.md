# Leaf Auxin Production Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add growth-coupled auxin production from young leaves and growth-modulated apical auxin production.

**Architecture:** Rename `auxin_production_rate` to `apical_auxin_baseline`, add 3 new genome params, modify `ApicalNode::produce_auxin()` to accept a growth fraction multiplier, add auxin production to `LeafNode::grow_size()`.

**Tech Stack:** C++17, Catch2, CMake

---

### Task 1: Genome rename and new parameters

**Files:**
- Modify: `src/engine/genome.h`
- Modify: `src/engine/node/tissues/apical.cpp`
- Modify: `src/evolution/genome_bridge.cpp`
- Modify: `tests/test_hormone.cpp`
- Modify: `tests/test_meristem.cpp`

- [ ] **Step 1: Rename `auxin_production_rate` to `apical_auxin_baseline` and add new fields in genome.h**

In `src/engine/genome.h`, replace the first field in the Genome struct:

```cpp
    // was: auxin_production_rate
    float apical_auxin_baseline;
    float apical_growth_auxin_multiplier; // growth-scaled bonus: total = baseline * (1 + multiplier * gf)
```

Add after the `auxin_bias` field:

```cpp
    float leaf_auxin_baseline;            // scaling constant for leaf auxin production (decoupled from apical)
    float leaf_growth_auxin_multiplier;   // fraction of leaf_auxin_baseline produced at max leaf growth
```

In `default_genome()`, replace:

```cpp
        .apical_auxin_baseline = 0.15f,
        .apical_growth_auxin_multiplier = 2.0f,  // total = baseline * 3 at max growth
```

And add after the `.auxin_bias` line:

```cpp
        .leaf_auxin_baseline = 0.15f,             // same scale as apical, but multiplier keeps it at 10%
        .leaf_growth_auxin_multiplier = 0.1f,     // single leaf at max growth = 10% of apical baseline
```

- [ ] **Step 2: Fix apical.cpp reference**

In `src/engine/node/tissues/apical.cpp`, in `produce_auxin()`, change:

```cpp
    float base = g.auxin_production_rate;
```

to:

```cpp
    float base = g.apical_auxin_baseline;
```

- [ ] **Step 3: Fix genome_bridge.cpp references**

In `src/evolution/genome_bridge.cpp`, in `build_genome_template()`, change:

```cpp
    reg(sg, "auxin_production_rate",      g.auxin_production_rate,      r, 0.01f, 2.0f, p);
```

to:

```cpp
    reg(sg, "apical_auxin_baseline",      g.apical_auxin_baseline,      r, 0.01f, 2.0f, p);
```

In `from_structured()`, change:

```cpp
    g.auxin_production_rate   = sg.get("auxin_production_rate");
```

to:

```cpp
    g.apical_auxin_baseline   = sg.get("apical_auxin_baseline");
```

In the auxin linkage group, change:

```cpp
        "auxin_production_rate", "auxin_diffusion_rate",
```

to:

```cpp
        "apical_auxin_baseline", "auxin_diffusion_rate",
```

- [ ] **Step 4: Fix test_hormone.cpp references**

In `tests/test_hormone.cpp`, replace all 4 occurrences of `g.auxin_production_rate` with `g.apical_auxin_baseline`:

Line 22: `g.apical_auxin_baseline = 1.0f;`
Line 39: `g.apical_auxin_baseline = 1.0f;`
Line 57: `g.apical_auxin_baseline = 1.0f;`
Line 74: `g.apical_auxin_baseline = 1.0f;`

- [ ] **Step 5: Fix test_meristem.cpp reference**

In `tests/test_meristem.cpp`, line 179, change:

```cpp
    g.auxin_production_rate = 0.0f; // disable auxin production so manual values stick
```

to:

```cpp
    g.apical_auxin_baseline = 0.0f; // disable auxin production so manual values stick
```

- [ ] **Step 6: Build and run tests**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests
```

Expected: all existing tests pass. No behavioral change yet — just a rename + unused new fields.

- [ ] **Step 7: Commit**

```bash
git add src/engine/genome.h src/engine/node/tissues/apical.cpp src/evolution/genome_bridge.cpp tests/test_hormone.cpp tests/test_meristem.cpp
git commit -m "refactor: rename auxin_production_rate to apical_auxin_baseline, add leaf auxin genome params"
```

---

### Task 2: Apical growth-modulated auxin production (TDD)

**Files:**
- Modify: `tests/test_hormone.cpp`
- Modify: `src/engine/node/tissues/apical.h`
- Modify: `src/engine/node/tissues/apical.cpp`

- [ ] **Step 1: Write failing test**

Add to `tests/test_hormone.cpp`:

```cpp
TEST_CASE("Auxin: apical growth multiplier boosts production with sugar", "[hormone]") {
    WorldParams world = default_world_params();

    // Plant 1: growth multiplier disabled
    Genome g1 = default_genome();
    g1.apical_auxin_baseline = 1.0f;
    g1.apical_growth_auxin_multiplier = 0.0f;  // no growth bonus
    Plant plant1(g1, glm::vec3(0.0f));

    // Plant 2: growth multiplier enabled
    Genome g2 = default_genome();
    g2.apical_auxin_baseline = 1.0f;
    g2.apical_growth_auxin_multiplier = 5.0f;  // large for clear signal
    Plant plant2(g2, glm::vec3(0.0f));

    // Saturate sugar + cytokinin on both → max growth → max multiplier
    for (int i = 0; i < 3; i++) {
        plant1.for_each_node_mut([](Node& n) {
            n.chemical(ChemicalID::Sugar) = 100.0f;
            n.chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant2.for_each_node_mut([](Node& n) {
            n.chemical(ChemicalID::Sugar) = 100.0f;
            n.chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant1.tick(world);
        plant2.tick(world);
    }

    // Sum total auxin in each plant
    float total1 = 0.0f, total2 = 0.0f;
    plant1.for_each_node([&](const Node& n) { total1 += n.chemical(ChemicalID::Auxin); });
    plant2.for_each_node([&](const Node& n) { total2 += n.chemical(ChemicalID::Auxin); });

    // Growth-boosted plant should have meaningfully more total auxin
    REQUIRE(total2 > total1 * 1.5f);
}

TEST_CASE("Auxin: apical growth multiplier has no effect without sugar", "[hormone]") {
    WorldParams world = default_world_params();

    Genome g1 = default_genome();
    g1.apical_auxin_baseline = 1.0f;
    g1.apical_growth_auxin_multiplier = 0.0f;
    Plant plant1(g1, glm::vec3(0.0f));

    Genome g2 = default_genome();
    g2.apical_auxin_baseline = 1.0f;
    g2.apical_growth_auxin_multiplier = 5.0f;
    Plant plant2(g2, glm::vec3(0.0f));

    // Zero sugar → zero growth → multiplier should have no effect
    for (int i = 0; i < 3; i++) {
        plant1.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 0.0f; });
        plant2.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 0.0f; });
        plant1.tick(world);
        plant2.tick(world);
    }

    float total1 = 0.0f, total2 = 0.0f;
    plant1.for_each_node([&](const Node& n) { total1 += n.chemical(ChemicalID::Auxin); });
    plant2.for_each_node([&](const Node& n) { total2 += n.chemical(ChemicalID::Auxin); });

    // Both should produce similar auxin (baseline only, growth_gf ≈ 0)
    // Allow 20% tolerance for sugar_factor floor (0.1 minimum) interacting with multiplier
    float ratio = (total1 > 1e-8f) ? total2 / total1 : 1.0f;
    REQUIRE(ratio < 1.2f);
}
```

- [ ] **Step 2: Build and verify tests fail**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests "[hormone]"
```

Expected: both new tests FAIL — `apical_growth_auxin_multiplier` is not yet used in production code, so both plants produce the same auxin regardless of the multiplier value.

- [ ] **Step 3: Update apical.h signature**

In `src/engine/node/tissues/apical.h`, change:

```cpp
    float produce_auxin(const Plant& plant) const;
```

to:

```cpp
    float produce_auxin(const Plant& plant, float growth_gf) const;
```

- [ ] **Step 4: Implement growth-modulated production in apical.cpp**

In `src/engine/node/tissues/apical.cpp`, in `tissue_tick()`, replace:

```cpp
    // Auxin production
    chemical(ChemicalID::Auxin) += produce_auxin(plant);
```

with:

```cpp
    // Compute growth fraction before auxin production — same inputs grow_tip() will use
    float max_cost = g.growth_rate * world.sugar_cost_meristem_growth;
    float growth_gf = meristem_helpers::growth_fraction(
        chemical(ChemicalID::Sugar), max_cost,
        chemical(ChemicalID::Cytokinin), g.cytokinin_growth_threshold);

    // Auxin production (growth-modulated)
    chemical(ChemicalID::Auxin) += produce_auxin(plant, growth_gf);
```

In `produce_auxin()`, change the signature and add the multiplier. Replace:

```cpp
float ApicalNode::produce_auxin(const Plant& plant) const {
    const Genome& g = plant.genome();
    float base = g.apical_auxin_baseline;

    // Light: shade boosts production (shade-avoidance / TAA-YUCCA upregulation)
    float local_light = estimate_local_light();
    float light_factor = 1.0f + g.auxin_shade_boost * (1.0f - local_light);

    // Sugar: well-fed meristems produce more (TOR kinase pathway)
    // Floor of 10% — even starving meristems have stored metabolic precursors
    float sugar = chemical(ChemicalID::Sugar);
    float sugar_factor = 0.1f + 0.9f * sugar / (sugar + g.auxin_sugar_half_saturation);

    // Age: young meristems are more biosynthetically active
    float age_factor = 1.0f / (1.0f + static_cast<float>(age) / g.auxin_age_half_life);

    return base * light_factor * sugar_factor * age_factor;
}
```

with:

```cpp
float ApicalNode::produce_auxin(const Plant& plant, float growth_gf) const {
    const Genome& g = plant.genome();
    float base = g.apical_auxin_baseline;

    // Light: shade boosts production (shade-avoidance / TAA-YUCCA upregulation)
    float local_light = estimate_local_light();
    float light_factor = 1.0f + g.auxin_shade_boost * (1.0f - local_light);

    // Sugar: well-fed meristems produce more (TOR kinase pathway)
    // Floor of 10% — even starving meristems have stored metabolic precursors
    float sugar = chemical(ChemicalID::Sugar);
    float sugar_factor = 0.1f + 0.9f * sugar / (sugar + g.auxin_sugar_half_saturation);

    // Age: young meristems are more biosynthetically active
    float age_factor = 1.0f / (1.0f + static_cast<float>(age) / g.auxin_age_half_life);

    float modulated_baseline = base * light_factor * sugar_factor * age_factor;

    // Growth modulation: actively growing meristems produce more auxin
    return modulated_baseline * (1.0f + g.apical_growth_auxin_multiplier * growth_gf);
}
```

- [ ] **Step 5: Build and run tests**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests
```

Expected: ALL tests pass, including the two new ones.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/tissues/apical.h src/engine/node/tissues/apical.cpp tests/test_hormone.cpp
git commit -m "feat: growth-modulated apical auxin production"
```

---

### Task 3: Leaf auxin production (TDD)

**Files:**
- Modify: `tests/test_hormone.cpp`
- Modify: `src/engine/node/tissues/leaf.cpp`

- [ ] **Step 1: Write failing tests**

Add to `tests/test_hormone.cpp` (include `leaf.h` at the top if not already present):

```cpp
#include "engine/node/tissues/leaf.h"
```

Then add:

```cpp
TEST_CASE("Auxin: growing leaf produces auxin", "[hormone]") {
    Genome g = default_genome();
    g.apical_auxin_baseline = 0.0f;          // disable meristem auxin
    g.leaf_auxin_baseline = 1.0f;             // high for easy detection
    g.leaf_growth_auxin_multiplier = 0.5f;
    g.growth_rate = 0.5f;
    g.shoot_plastochron = 1;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world_params();

    // Grow plant until leaves exist
    for (int i = 0; i < 5; i++) {
        plant.for_each_node_mut([](Node& n) {
            n.chemical(ChemicalID::Sugar) = 100.0f;
            n.chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant.tick(world);
    }

    // Find a growing leaf
    LeafNode* growing_leaf = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (auto* leaf = n.as_leaf()) {
            if (leaf->leaf_size < g.max_leaf_size) growing_leaf = leaf;
        }
    });
    REQUIRE(growing_leaf != nullptr);

    // Zero all auxin, give sugar for leaf growth
    plant.for_each_node_mut([](Node& n) {
        n.chemical(ChemicalID::Auxin) = 0.0f;
        n.chemical(ChemicalID::Sugar) = 100.0f;
    });

    plant.tick(world);

    // Growing leaf should have produced auxin (some may have transported, but should retain some)
    REQUIRE(growing_leaf->chemical(ChemicalID::Auxin) > 0.0f);
}

TEST_CASE("Auxin: full-size leaf produces zero auxin", "[hormone]") {
    Genome g = default_genome();
    g.apical_auxin_baseline = 0.0f;          // disable meristem auxin
    g.leaf_auxin_baseline = 1.0f;
    g.leaf_growth_auxin_multiplier = 0.5f;
    g.growth_rate = 0.5f;
    g.shoot_plastochron = 1;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world_params();

    // Grow plant until leaves exist
    for (int i = 0; i < 5; i++) {
        plant.for_each_node_mut([](Node& n) {
            n.chemical(ChemicalID::Sugar) = 100.0f;
            n.chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant.tick(world);
    }

    // Force ALL leaves to max size
    plant.for_each_node_mut([&](Node& n) {
        if (auto* leaf = n.as_leaf()) {
            leaf->leaf_size = g.max_leaf_size;
        }
    });

    // Zero all auxin everywhere
    plant.for_each_node_mut([](Node& n) {
        n.chemical(ChemicalID::Auxin) = 0.0f;
        n.chemical(ChemicalID::Sugar) = 100.0f;
    });

    plant.tick(world);

    // No auxin source anywhere → all auxin should be zero
    float total_auxin = 0.0f;
    plant.for_each_node([&](const Node& n) { total_auxin += n.chemical(ChemicalID::Auxin); });
    REQUIRE(total_auxin < 1e-6f);
}

TEST_CASE("Auxin: leaf auxin scales with growth amount", "[hormone]") {
    WorldParams world = default_world_params();

    // Plant 1: low leaf multiplier
    Genome g1 = default_genome();
    g1.apical_auxin_baseline = 0.0f;
    g1.leaf_auxin_baseline = 1.0f;
    g1.leaf_growth_auxin_multiplier = 0.1f;
    g1.growth_rate = 0.5f;
    g1.shoot_plastochron = 1;
    Plant plant1(g1, glm::vec3(0.0f));

    // Plant 2: high leaf multiplier
    Genome g2 = default_genome();
    g2.apical_auxin_baseline = 0.0f;
    g2.leaf_auxin_baseline = 1.0f;
    g2.leaf_growth_auxin_multiplier = 0.9f;
    g2.growth_rate = 0.5f;
    g2.shoot_plastochron = 1;
    Plant plant2(g2, glm::vec3(0.0f));

    // Grow both until leaves exist, then measure auxin
    for (int i = 0; i < 5; i++) {
        plant1.for_each_node_mut([](Node& n) {
            n.chemical(ChemicalID::Sugar) = 100.0f;
            n.chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant2.for_each_node_mut([](Node& n) {
            n.chemical(ChemicalID::Sugar) = 100.0f;
            n.chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant1.tick(world);
        plant2.tick(world);
    }

    // Zero auxin, give sugar, tick once more
    plant1.for_each_node_mut([](Node& n) {
        n.chemical(ChemicalID::Auxin) = 0.0f;
        n.chemical(ChemicalID::Sugar) = 100.0f;
    });
    plant2.for_each_node_mut([](Node& n) {
        n.chemical(ChemicalID::Auxin) = 0.0f;
        n.chemical(ChemicalID::Sugar) = 100.0f;
    });
    plant1.tick(world);
    plant2.tick(world);

    float total1 = 0.0f, total2 = 0.0f;
    plant1.for_each_node([&](const Node& n) { total1 += n.chemical(ChemicalID::Auxin); });
    plant2.for_each_node([&](const Node& n) { total2 += n.chemical(ChemicalID::Auxin); });

    // Higher multiplier → more leaf auxin in the system
    REQUIRE(total2 > total1 * 2.0f);
}
```

- [ ] **Step 2: Build and verify tests fail**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests "[hormone]"
```

Expected: the "growing leaf produces auxin" and "leaf auxin scales" tests FAIL (leaf doesn't produce auxin yet). The "full-size leaf produces zero" test may pass trivially (no auxin source).

- [ ] **Step 3: Implement leaf auxin production in leaf.cpp**

In `src/engine/node/tissues/leaf.cpp`, in `grow_size()`, add auxin production after the `leaf_size += growth;` line and before the petiole extension block. Replace:

```cpp
    leaf_size += growth;
    chemical(ChemicalID::Sugar) -= cost;

    // Extend petiole proportionally as leaf grows
```

with:

```cpp
    leaf_size += growth;
    chemical(ChemicalID::Sugar) -= cost;

    // Auxin production: growing leaves produce auxin proportional to growth rate.
    // No growth (full size, stressed, starved) → zero auxin.
    float growth_fraction = growth / g.leaf_growth_rate;
    chemical(ChemicalID::Auxin) += growth_fraction * g.leaf_growth_auxin_multiplier * g.leaf_auxin_baseline;

    // Extend petiole proportionally as leaf grows
```

- [ ] **Step 4: Build and run tests**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests
```

Expected: ALL tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/engine/node/tissues/leaf.cpp tests/test_hormone.cpp
git commit -m "feat: growth-coupled auxin production from young leaves"
```

---

### Task 4: Evolution bridge update

**Files:**
- Modify: `src/evolution/genome_bridge.cpp`

- [ ] **Step 1: Add new genes to build_genome_template()**

In `src/evolution/genome_bridge.cpp`, in `build_genome_template()`, add after the existing `"apical_auxin_baseline"` line (already renamed in Task 1):

```cpp
    reg(sg, "apical_growth_auxin_multiplier", g.apical_growth_auxin_multiplier, r, 0.0f, 10.0f, p);
```

Add after the `"auxin_bias"` line:

```cpp
    reg(sg, "leaf_auxin_baseline",            g.leaf_auxin_baseline,            r, 0.01f, 2.0f, p);
    reg(sg, "leaf_growth_auxin_multiplier",   g.leaf_growth_auxin_multiplier,   r, 0.0f, 1.0f, p);
```

- [ ] **Step 2: Add new genes to from_structured()**

In `from_structured()`, add after the `g.apical_auxin_baseline` line:

```cpp
    g.apical_growth_auxin_multiplier = sg.get("apical_growth_auxin_multiplier");
```

Add after the `g.auxin_bias` line:

```cpp
    g.leaf_auxin_baseline              = sg.get("leaf_auxin_baseline");
    g.leaf_growth_auxin_multiplier     = sg.get("leaf_growth_auxin_multiplier");
```

- [ ] **Step 3: Update auxin linkage group**

In the auxin linkage group, replace:

```cpp
    sg.add_linkage_group({"auxin", {
        "apical_auxin_baseline", "auxin_diffusion_rate",
        "auxin_decay_rate", "auxin_threshold",
        "auxin_shade_boost", "auxin_sugar_half_saturation", "auxin_age_half_life",
        "auxin_bias"
    }});
```

with:

```cpp
    sg.add_linkage_group({"auxin", {
        "apical_auxin_baseline", "apical_growth_auxin_multiplier",
        "auxin_diffusion_rate",
        "auxin_decay_rate", "auxin_threshold",
        "auxin_shade_boost", "auxin_sugar_half_saturation", "auxin_age_half_life",
        "auxin_bias",
        "leaf_auxin_baseline", "leaf_growth_auxin_multiplier"
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
git commit -m "feat: add leaf auxin genome params to evolution bridge"
```

---

### Task 5: Update CLAUDE.md

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Update auxin production description**

In `CLAUDE.md`, in the "Chemical Transport Model" section under **Auxin**, replace:

```
- Produced by active `ShootApicalNode` during its `tick()`
```

with:

```
- Produced by active `ShootApicalNode` during its `tick()`, modulated by growth rate (growth_gf multiplier)
- Also produced by growing `LeafNode` during `grow_size()` — proportional to growth rate, zero when full-size
```

- [ ] **Step 2: Update tuning parameters**

In the "Tuning Parameters" section, replace:

```
- `auxin_threshold` (0.15) - lower = fewer shoot branches, higher = more
```

with:

```
- `apical_auxin_baseline` (0.15) - base auxin output of shoot apical meristem per tick
- `apical_growth_auxin_multiplier` (2.0) - growth bonus: total = baseline * (1 + multiplier * growth_fraction). 0 = no bonus, 2 = 3x at max growth
- `leaf_auxin_baseline` (0.15) - scaling constant for leaf auxin production (decoupled from apical rate)
- `leaf_growth_auxin_multiplier` (0.1) - fraction of leaf_auxin_baseline at max growth. Single leaf at max growth = 10% of apical baseline
- `auxin_threshold` (0.15) - lower = fewer shoot branches, higher = more
```

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: update CLAUDE.md with leaf auxin production details"
```
