# Root↔Shoot Balance via Metabolic Hormone Feedback — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add sugar-and-water metabolic gating to meristem hormone production, plus a quiescence transition (active→dormant under sustained stress), so shoot and root growth self-balance via local sensing at each meristem.

**Architecture:** Two-mechanism design. Mechanism 1 is a per-tick multiplicative factor `metabolic_factor(sugar, water)` applied at every hormone-production site — sugar and water each contribute a Michaelis-Menten term with a hormone-specific floor (auxin 0.1, cytokinin 0.05). Mechanism 2 reuses the existing `starvation_ticks` counter: at a new `quiescence_threshold` (before the existing death threshold), active meristems flip to dormant instead of dying. When conditions recover, normal `can_activate()` can reactivate them.

**Tech Stack:** C++17, CMake, Catch2 tests. Build via `/usr/local/bin/cmake --build build`. Tests via `./build/botany_tests`.

**Design context:** See [docs/long-term-plan/milestone-2/root-shoot-balance.md](../../long-term-plan/milestone-2/root-shoot-balance.md) for the biological rationale and design reasoning. This plan implements that design.

---

## File Structure

**Create:**
- `docs/superpowers/plans/2026-04-19-root-shoot-balance-feedback.md` — this plan

**Modify:**
- `src/engine/genome.h` — add 4 new genome parameters (Task 1)
- `src/engine/node/meristems/helpers.h` — add `metabolic_factor` helper function (Task 2)
- `src/engine/node/tissues/apical.cpp` — upgrade SA auxin gating, add quiescence (Tasks 3, 7)
- `src/engine/node/tissues/root_apical.cpp` — add RA auxin/cyto gating, add quiescence (Tasks 4, 5, 6)
- `tests/test_meristem.cpp` — add tests for the new behavior (Tasks 2-7)

Each file changes either the genome, the hormone-production formula at one meristem type, or adds quiescence. Changes are additive and don't cross concerns — should be reviewable independently.

---

## Task 1: Add new genome parameters

**Files:**
- Modify: `src/engine/genome.h`

**Rationale:** Before we can use new parameters in code, they need to exist in the `Genome` struct with sensible defaults. Four new floats: two half-saturation constants for the water-MM terms (one per hormone), one new sugar-half-saturation for cytokinin, and one quiescence threshold. Floors (0.1 for auxin, 0.05 for cytokinin) are hard-coded constants in the formulas; they're not per-species.

- [ ] **Step 1: Read existing genome structure**

Run: `grep -n "auxin_sugar_half_saturation\|starvation" /Users/wldarden/learning/botany/src/engine/genome.h | head -10`

Expected output: shows `auxin_sugar_half_saturation` field around line 20 and its default around line 217. Note the formatting pattern (field declarations, then defaults block).

- [ ] **Step 2: Add new fields to the Genome struct**

In `src/engine/genome.h`, find the line declaring `auxin_sugar_half_saturation` (around line 20). Add these fields immediately after it, preserving the existing comment style:

```cpp
float auxin_water_half_saturation;      // ml — water level at which auxin production hits half-max. Floor 0.1 is encoded in the formula, not evolvable.
float cytokinin_sugar_half_saturation;  // g glucose — sugar level at which CK production hits half-max.
float cytokinin_water_half_saturation;  // ml — water level at which CK production hits half-max.
float quiescence_threshold;             // ticks — active meristem reverts to dormant after this many consecutive low-sugar ticks (before starvation_ticks_max death).
```

- [ ] **Step 3: Add defaults in default_genome()**

In the same file, find `.auxin_sugar_half_saturation = 0.3f,` (around line 217). Add four new default lines immediately after it, in the same order as the struct fields:

```cpp
        .auxin_water_half_saturation = 0.1f,     // matches SA stomatal threshold magnitude — half-max at moderate turgor
        .cytokinin_sugar_half_saturation = 0.05f, // CK more sensitive to sugar than auxin; sharper response
        .cytokinin_water_half_saturation = 0.1f,  // matches auxin water sensitivity
        .quiescence_threshold = 150.0f,           // ~6 days — meristem goes dormant after sustained starvation, well before death (starvation_ticks_max = 2200)
```

- [ ] **Step 4: Build to verify no compile errors**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -5`
Expected: `[100%] Built target botany_tests` with no errors.

- [ ] **Step 5: Run existing tests to verify nothing broke**

Run: `./build/botany_tests 2>&1 | tail -3`
Expected: `All tests passed` message. No new tests yet, just verify we didn't break anything.

- [ ] **Step 6: Commit**

```bash
git add src/engine/genome.h
git commit -m "feat: add metabolic_factor and quiescence genome parameters

Four new parameters for the root-shoot balance feedback design:
- auxin_water_half_saturation (0.1 ml)
- cytokinin_sugar_half_saturation (0.05 g)
- cytokinin_water_half_saturation (0.1 ml)
- quiescence_threshold (150 ticks)

No behavior change yet; parameters unused until subsequent tasks wire them
into hormone production formulas and quiescence transitions.

Plan: docs/superpowers/plans/2026-04-19-root-shoot-balance-feedback.md"
```

---

## Task 2: Add `metabolic_factor` helper function

**Files:**
- Modify: `src/engine/node/meristems/helpers.h`
- Modify: `tests/test_meristem.cpp`

**Rationale:** The helper is pure — no state, no side effects. Reused at three production sites (Tasks 3, 4, 5). Writing it first with a focused test establishes the contract and makes the later tasks trivial.

Formula:
```
metabolic_factor(sugar, K_s, floor_s, water, K_w, floor_w)
    = [floor_s + (1 - floor_s) × sugar/(sugar + K_s)]
    × [floor_w + (1 - floor_w) × water/(water + K_w)]
```

Invariants:
- At sugar=0, water=0: result is `floor_s × floor_w` (e.g. 0.1 × 0.1 = 0.01).
- At sugar → ∞, water → ∞: result is 1.0 (saturates).
- At sugar = K_s, water = K_w: each term is `floor + (1-floor) × 0.5`, product is (floor + (1-floor)/2)².
- Monotonically non-decreasing in both sugar and water.

- [ ] **Step 1: Write failing test for the helper**

Add to `tests/test_meristem.cpp` at the end of the file (before any closing braces — all TEST_CASE blocks are at namespace scope):

```cpp
TEST_CASE("metabolic_factor: saturates at high inputs", "[meristem][metabolic]") {
    using meristem_helpers::metabolic_factor;
    float mf = metabolic_factor(/*sugar*/ 1000.0f, /*K_s*/ 0.1f, /*floor_s*/ 0.1f,
                                 /*water*/ 1000.0f, /*K_w*/ 0.1f, /*floor_w*/ 0.1f);
    REQUIRE_THAT(mf, WithinAbs(1.0f, 1e-3f));
}

TEST_CASE("metabolic_factor: floor at zero inputs", "[meristem][metabolic]") {
    using meristem_helpers::metabolic_factor;
    float mf = metabolic_factor(0.0f, 0.1f, 0.1f, 0.0f, 0.1f, 0.05f);
    // floor_s * floor_w = 0.1 * 0.05 = 0.005
    REQUIRE_THAT(mf, WithinAbs(0.005f, 1e-6f));
}

TEST_CASE("metabolic_factor: half-saturation point", "[meristem][metabolic]") {
    using meristem_helpers::metabolic_factor;
    // At sugar = K, sugar term = floor + (1 - floor) * 0.5 = 0.1 + 0.45 = 0.55
    // Same for water. Product = 0.55 * 0.55 = 0.3025.
    float mf = metabolic_factor(0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f);
    REQUIRE_THAT(mf, WithinAbs(0.3025f, 1e-3f));
}

TEST_CASE("metabolic_factor: stresses compound multiplicatively", "[meristem][metabolic]") {
    using meristem_helpers::metabolic_factor;
    // Full sugar, zero water → should reach only water floor (0.1)
    float mf_dry = metabolic_factor(1000.0f, 0.1f, 0.1f, 0.0f, 0.1f, 0.1f);
    REQUIRE_THAT(mf_dry, WithinAbs(0.1f, 1e-3f));
    // Zero sugar, full water → should reach only sugar floor (0.1)
    float mf_starved = metabolic_factor(0.0f, 0.1f, 0.1f, 1000.0f, 0.1f, 0.1f);
    REQUIRE_THAT(mf_starved, WithinAbs(0.1f, 1e-3f));
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | grep -E "error:" | head -5`

Expected: compile error indicating `metabolic_factor` is not declared in `meristem_helpers` namespace. This is the expected failure — proves the test is looking for a function we haven't written yet.

- [ ] **Step 3: Add the `metabolic_factor` helper**

In `src/engine/node/meristems/helpers.h`, add this function inside the `meristem_helpers` namespace, immediately after the existing `auxin_growth_factor` function (around line 122, before the closing `} // namespace meristem_helpers`):

```cpp
// Multiplicative metabolic gating — captures short-timescale biochemistry +
// transcription response to sugar and water at a meristem.  Both stresses
// compound.  Floors encode hormone-specific buffering (conjugate pools for
// auxin, etc.); caller supplies floor per hormone type.
//
// At sugar=0 and water=0, returns floor_s * floor_w (minimum production).
// At saturation (sugar >> K_s, water >> K_w), returns 1.0 (full production).
//
// Used by SA auxin, RA auxin, and RA cytokinin production.
inline float metabolic_factor(float sugar, float K_sugar, float floor_sugar,
                              float water, float K_water, float floor_water) {
    float sf = floor_sugar + (1.0f - floor_sugar) * sugar / (sugar + K_sugar);
    float wf = floor_water + (1.0f - floor_water) * water / (water + K_water);
    return sf * wf;
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[metabolic]" 2>&1 | tail -5`
Expected: all 4 `[metabolic]` test cases pass.

- [ ] **Step 5: Verify all existing tests still pass**

Run: `./build/botany_tests 2>&1 | tail -3`
Expected: `All tests passed` (no regression).

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/meristems/helpers.h tests/test_meristem.cpp
git commit -m "feat: add metabolic_factor helper for hormone production gating

Pure multiplicative MM gating on sugar + water with per-call floor.
Saturates at 1.0 when both inputs are far above K; reaches floor_s*floor_w
at full starvation.  Used at all meristem hormone production sites in
subsequent tasks.

Plan: docs/superpowers/plans/2026-04-19-root-shoot-balance-feedback.md"
```

---

## Task 3: Upgrade SA auxin production to sugar+water gating

**Files:**
- Modify: `src/engine/node/tissues/apical.cpp` (the `produce_auxin` function, around line 32)
- Modify: `tests/test_meristem.cpp`

**Rationale:** Shoot apical already has a sugar-only gate (`0.1 + 0.9 * sugar / (sugar + K)`). Upgrade it to use `metabolic_factor(sugar, ..., water, ...)` which also gates on water. Biologically correct — drought suppresses YUC/TAA1 auxin biosynthesis in shoot apicals via ABA signaling.

The existing formula:
```cpp
float sugar_factor = 0.1f + 0.9f * sugar / (sugar + g.auxin_sugar_half_saturation);
float modulated_baseline = base * light_factor * sugar_factor;
```

becomes:

```cpp
float mf = metabolic_factor(
    chemical(ChemicalID::Sugar), g.auxin_sugar_half_saturation, 0.1f,
    chemical(ChemicalID::Water), g.auxin_water_half_saturation, 0.1f);
float modulated_baseline = base * light_factor * mf;
```

- [ ] **Step 1: Write failing test — SA auxin drops when water is low**

Add to `tests/test_meristem.cpp`:

```cpp
TEST_CASE("SA auxin production drops when water is low", "[meristem][metabolic]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world;

    // Find the primary SA (id 1) and give it abundant sugar but no water
    ApicalNode* sa = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::APICAL && !sa) sa = n.as_apical();
    });
    REQUIRE(sa != nullptr);
    REQUIRE(sa->active);

    // Abundant sugar, zero water: auxin production should be close to sugar-gated
    // max × water floor (0.1), not full rate.
    sa->chemical(ChemicalID::Sugar) = 10.0f;   // far above K_sugar
    sa->chemical(ChemicalID::Water) = 0.0f;
    sa->tick_auxin_produced = 0.0f;

    plant.tick(world);

    float produced_dry = sa->tick_auxin_produced;

    // Reset and test with full water
    sa->chemical(ChemicalID::Sugar) = 10.0f;
    sa->chemical(ChemicalID::Water) = water_cap(*sa, g);
    sa->tick_auxin_produced = 0.0f;

    plant.tick(world);
    float produced_wet = sa->tick_auxin_produced;

    // Dry SA should produce significantly less than wet SA (water gate)
    REQUIRE(produced_dry < produced_wet * 0.5f);
    // But dry SA should still produce something (floor 0.1)
    REQUIRE(produced_dry > 0.0f);
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "SA auxin production drops when water is low" 2>&1 | tail -5`

Expected: test fails — current code doesn't gate on water, so `produced_dry ≈ produced_wet` and the `<` assertion fails. This is the expected failure.

- [ ] **Step 3: Apply the fix to SA auxin production**

In `src/engine/node/tissues/apical.cpp`, find the `produce_auxin` function (around line 32). Replace the `sugar_factor` computation and `modulated_baseline` line. Specifically, find:

```cpp
    float base = g.apical_auxin_baseline;
    float local_light = estimate_local_light();
    float light_factor = 1.0f + g.auxin_shade_boost * (1.0f - local_light);
    float sugar = chemical(ChemicalID::Sugar);
    float sugar_factor = 0.1f + 0.9f * sugar / (sugar + g.auxin_sugar_half_saturation);
    float modulated_baseline = base * light_factor * sugar_factor;
```

Replace with:

```cpp
    float base = g.apical_auxin_baseline;
    float local_light = estimate_local_light();
    float light_factor = 1.0f + g.auxin_shade_boost * (1.0f - local_light);
    // Metabolic gating: sugar and water both required for full auxin synthesis.
    // Floor 0.1 matches conjugate-pool buffering — zero metabolism still yields ~1%.
    float mf = metabolic_factor(
        chemical(ChemicalID::Sugar), g.auxin_sugar_half_saturation, 0.1f,
        chemical(ChemicalID::Water), g.auxin_water_half_saturation, 0.1f);
    float modulated_baseline = base * light_factor * mf;
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "SA auxin production drops when water is low" 2>&1 | tail -5`
Expected: test passes.

- [ ] **Step 5: Run all tests to verify no regression**

Run: `./build/botany_tests 2>&1 | tail -3`
Expected: `All tests passed`. Note: some tests may run with different assertion counts because the MM curve shape changed slightly, but all should pass.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/tissues/apical.cpp tests/test_meristem.cpp
git commit -m "feat: SA auxin production gates on water as well as sugar

Shoot apical auxin production used to be sugar-only gated.  Now uses
metabolic_factor(sugar, water) — sugar + water stresses compound
multiplicatively, matching biology (drought represses YUC/TAA1 sharply
via ABA signaling).

Plan: docs/superpowers/plans/2026-04-19-root-shoot-balance-feedback.md"
```

---

## Task 4: Add metabolic gating to RA auxin production

**Files:**
- Modify: `src/engine/node/tissues/root_apical.cpp` (the `update_tissue` function, around line 35)
- Modify: `tests/test_meristem.cpp`

**Rationale:** Root apical currently produces auxin at a fixed rate (`chemical(Auxin) += g.root_tip_auxin_production_rate`) with no sugar/water gating — a biological and symmetry gap. Real root tip auxin has a shoot-delivered component (PIN, handled elsewhere) and a local synthesis component; our parameter represents the local synthesis, which is energy-dependent. Same gate as SA, same floor 0.1.

- [ ] **Step 1: Write failing test — RA auxin drops when sugar is low**

Add to `tests/test_meristem.cpp`:

```cpp
TEST_CASE("RA auxin production drops when sugar is low", "[meristem][metabolic]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world;

    RootApicalNode* ra = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::ROOT_APICAL && !ra) ra = n.as_root_apical();
    });
    REQUIRE(ra != nullptr);
    REQUIRE(ra->active);  // primary RA is born active

    // Abundant water, zero sugar: auxin production should be water-gate max × sugar floor
    ra->chemical(ChemicalID::Sugar) = 0.0f;
    ra->chemical(ChemicalID::Water) = water_cap(*ra, g);
    ra->tick_auxin_produced = 0.0f;

    plant.tick(world);
    float produced_starved = ra->tick_auxin_produced;

    // Reset with abundant sugar
    ra->chemical(ChemicalID::Sugar) = 1.0f;  // >> K_sugar (0.3)
    ra->chemical(ChemicalID::Water) = water_cap(*ra, g);
    ra->tick_auxin_produced = 0.0f;

    plant.tick(world);
    float produced_fed = ra->tick_auxin_produced;

    REQUIRE(produced_starved < produced_fed * 0.3f);  // sugar-starved meristem hits ~0.1 × water_factor
    REQUIRE(produced_starved > 0.0f);  // floor keeps it non-zero
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "RA auxin production drops when sugar is low" 2>&1 | tail -5`

Expected: test fails — current RA auxin production is a fixed rate, so `produced_starved ≈ produced_fed`.

- [ ] **Step 3: Apply the fix to RA auxin production**

In `src/engine/node/tissues/root_apical.cpp`, find the `update_tissue` function (around line 16). Inside the active branch (after `if (!active) { ... return; }`), find the auxin production lines:

```cpp
    // Active root tips: auxin self-maintenance (local PIN-recycling maximum) and
    // cytokinin production gated by local auxin ("more auxin → stronger cyto
    // signal").  Both productions match the shoot apical pattern.
    chemical(ChemicalID::Auxin) += g.root_tip_auxin_production_rate;
    tick_auxin_produced += g.root_tip_auxin_production_rate;
```

Replace with:

```cpp
    // Active root tips: auxin self-maintenance (local PIN-recycling maximum) and
    // cytokinin production gated by local auxin ("more auxin → stronger cyto
    // signal").  Both productions match the shoot apical pattern.
    // Metabolic gating: sugar + water each contribute an MM factor; floor 0.1
    // matches SA convention (auxin conjugate pools buffer short-term stress).
    float mf_auxin = metabolic_factor(
        chemical(ChemicalID::Sugar), g.auxin_sugar_half_saturation, 0.1f,
        chemical(ChemicalID::Water), g.auxin_water_half_saturation, 0.1f);
    float produced_auxin = g.root_tip_auxin_production_rate * mf_auxin;
    chemical(ChemicalID::Auxin) += produced_auxin;
    tick_auxin_produced += produced_auxin;
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "RA auxin production drops when sugar is low" 2>&1 | tail -5`
Expected: test passes.

- [ ] **Step 5: Run all tests**

Run: `./build/botany_tests 2>&1 | tail -3`
Expected: `All tests passed`.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/tissues/root_apical.cpp tests/test_meristem.cpp
git commit -m "feat: RA auxin production gates on sugar + water

Root apical auxin production was previously a fixed rate whenever the
meristem was active.  Now uses metabolic_factor(sugar, water) matching
the SA formula — closes the symmetry gap between shoot and root apicals.

Plan: docs/superpowers/plans/2026-04-19-root-shoot-balance-feedback.md"
```

---

## Task 5: Add metabolic gating to RA cytokinin production

**Files:**
- Modify: `src/engine/node/tissues/root_apical.cpp` (the `update_tissue` function, the cytokinin lines)
- Modify: `tests/test_meristem.cpp`

**Rationale:** This is the *main* shoot-runaway brake. Roots are the primary CK source; shoot branching gates on cytokinin crossing a threshold. Sugar-gating CK at the source creates the feedback: shoot out-grows root → phloem delivery to roots drops → RA sugar drops → CK production drops → less CK reaches shoot → fewer SA activations → shoot branching slows.

CK floor is **smaller** than auxin's (0.05 vs 0.1) because CK has less biological buffering (smaller conjugate pools, faster turnover). This yields a sharper response curve, intentional.

- [ ] **Step 1: Write failing test — RA cytokinin drops sharply with sugar**

Add to `tests/test_meristem.cpp`:

```cpp
TEST_CASE("RA cytokinin production drops sharply when sugar is low", "[meristem][metabolic]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world;

    RootApicalNode* ra = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::ROOT_APICAL && !ra) ra = n.as_root_apical();
    });
    REQUIRE(ra != nullptr);

    // Zero sugar, full water, abundant auxin — CK should drop to ~floor (0.05)
    ra->chemical(ChemicalID::Sugar) = 0.0f;
    ra->chemical(ChemicalID::Water) = water_cap(*ra, g);
    ra->chemical(ChemicalID::Auxin) = 1.0f;  // high auxin to isolate sugar effect
    ra->tick_cytokinin_produced = 0.0f;

    plant.tick(world);
    float produced_starved = ra->tick_cytokinin_produced;

    // Full sugar, full water, same auxin — CK should be near max
    ra->chemical(ChemicalID::Sugar) = 1.0f;  // >> K_sugar (0.05 for CK)
    ra->chemical(ChemicalID::Water) = water_cap(*ra, g);
    ra->chemical(ChemicalID::Auxin) = 1.0f;
    ra->tick_cytokinin_produced = 0.0f;

    plant.tick(world);
    float produced_fed = ra->tick_cytokinin_produced;

    // CK floor is smaller (0.05) so the gap should be wider than auxin's gap
    REQUIRE(produced_starved < produced_fed * 0.1f);
    REQUIRE(produced_starved > 0.0f);  // floor keeps it non-zero
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "RA cytokinin production drops sharply when sugar is low" 2>&1 | tail -5`

Expected: test fails. Current CK formula is `rate × auxin` — no sugar term at all, so starved produces same as fed.

- [ ] **Step 3: Apply the fix to RA cytokinin production**

In `src/engine/node/tissues/root_apical.cpp`, immediately after the auxin production block you edited in Task 4, find:

```cpp
    float cyto_produced = g.root_cytokinin_production_rate * chemical(ChemicalID::Auxin);
    chemical(ChemicalID::Cytokinin) += cyto_produced;
    tick_cytokinin_produced += cyto_produced;
```

Replace with:

```cpp
    // Cytokinin: floor 0.05 (smaller than auxin's 0.1) — CK has less
    // conjugate-pool buffering than auxin, so its response to sugar is sharper.
    // This is the primary root-to-shoot feedback brake: low root sugar → low
    // CK → less CK delivered to shoot via xylem → fewer SA activations.
    float mf_cyto = metabolic_factor(
        chemical(ChemicalID::Sugar), g.cytokinin_sugar_half_saturation, 0.05f,
        chemical(ChemicalID::Water), g.cytokinin_water_half_saturation, 0.05f);
    float cyto_produced = g.root_cytokinin_production_rate * chemical(ChemicalID::Auxin) * mf_cyto;
    chemical(ChemicalID::Cytokinin) += cyto_produced;
    tick_cytokinin_produced += cyto_produced;
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "RA cytokinin production drops sharply when sugar is low" 2>&1 | tail -5`
Expected: test passes.

- [ ] **Step 5: Run all tests**

Run: `./build/botany_tests 2>&1 | tail -3`
Expected: `All tests passed`.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/tissues/root_apical.cpp tests/test_meristem.cpp
git commit -m "feat: RA cytokinin production gates on sugar + water

Root apical CK production previously was rate × auxin with no metabolic
gating.  Now also multiplies by metabolic_factor(sugar, water) with a
smaller floor (0.05) than auxin — CK responds more sharply to sugar
because it has less conjugate-pool buffering.

This is the primary shoot-runaway brake in the root-shoot feedback design:
low root sugar → low CK → less CK delivered to shoot → fewer SA activations.

Plan: docs/superpowers/plans/2026-04-19-root-shoot-balance-feedback.md"
```

---

## Task 6: Add quiescence transition to RootApicalNode

**Files:**
- Modify: `src/engine/node/tissues/root_apical.cpp` (the `update_tissue` function, near the top)
- Modify: `tests/test_meristem.cpp`

**Rationale:** Mechanism 2 — sustained-stress dormancy. An active meristem that's been starved for `quiescence_threshold` ticks reverts to dormant instead of continuing to pay maintenance and eventually dying at `starvation_ticks_max`. Real root tips do this: they survive weeks of sugar starvation by going quiescent, then reactivate when conditions improve.

Implementation note: `starvation_ticks` is incremented in `Node::check_starvation` (runs before `update_tissue`). By the time `update_tissue` runs, the tick's increment has already happened. So we check `>= quiescence_threshold` at the top of `update_tissue`, flip to dormant, and reset `starvation_ticks` so it doesn't carry over.

- [ ] **Step 1: Write failing test — RA reverts to dormant after sustained starvation**

Add to `tests/test_meristem.cpp`:

```cpp
TEST_CASE("Active RA reverts to dormant after quiescence_threshold ticks of starvation", "[meristem][quiescence]") {
    Genome g = default_genome();
    // Speed up the test: lower quiescence threshold
    g.quiescence_threshold = 10.0f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world;

    RootApicalNode* ra = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::ROOT_APICAL && !ra) ra = n.as_root_apical();
    });
    REQUIRE(ra != nullptr);
    REQUIRE(ra->active);

    // Starve the RA — zero sugar, keep water so we isolate the sugar-based
    // starvation_ticks accumulation.
    for (int i = 0; i < 15; ++i) {
        ra->chemical(ChemicalID::Sugar) = 0.0f;
        plant.tick(world);
    }

    // After 10+ ticks of starvation it should have gone dormant
    REQUIRE_FALSE(ra->active);
    // And starvation_ticks should have reset (so dormant bud can eventually reactivate)
    REQUIRE(ra->starvation_ticks == 0u);
}

TEST_CASE("Dormant RA does not die at starvation_ticks_max", "[meristem][quiescence]") {
    Genome g = default_genome();
    g.quiescence_threshold = 5.0f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world;
    world.starvation_ticks_max = 20;  // also speed up the death check

    RootApicalNode* ra = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::ROOT_APICAL && !ra) ra = n.as_root_apical();
    });
    REQUIRE(ra != nullptr);

    uint32_t ra_id = ra->id;

    // Starve long enough that without quiescence, it would die
    for (int i = 0; i < 40; ++i) {
        plant.for_each_node_mut([&](Node& n) {
            n.chemical(ChemicalID::Sugar) = 0.0f;
        });
        plant.tick(world);
    }

    // The RA should still exist (went quiescent, didn't die)
    bool ra_still_alive = false;
    plant.for_each_node([&](const Node& n) {
        if (n.id == ra_id) ra_still_alive = true;
    });
    REQUIRE(ra_still_alive);
}
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[quiescence]" 2>&1 | tail -10`

Expected: both tests fail. Without the quiescence mechanism, active RAs stay active through starvation and eventually die.

- [ ] **Step 3: Add the quiescence check to RA update_tissue**

In `src/engine/node/tissues/root_apical.cpp`, find the `update_tissue` function. At the very top (after the `absorb_water(g, world);` line, before the `if (!active)` check), add:

```cpp
    // Quiescence: active meristem reverts to dormant after sustained sugar
    // starvation rather than dying.  Matches real root-tip biology — tips
    // survive weeks of stress by going quiescent, reactivate later via
    // normal can_activate() when sugar returns.  Runs before the active
    // check so the newly-dormant RA is handled by the dormant branch below.
    if (active && starvation_ticks >= static_cast<uint32_t>(g.quiescence_threshold)) {
        active = false;
        starvation_ticks = 0;  // dormant meristem pays no maintenance; clean slate for reactivation
    }
```

So the final function top looks like:

```cpp
void RootApicalNode::update_tissue(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    absorb_water(g, world);

    // Quiescence: active meristem reverts to dormant after sustained sugar
    // starvation rather than dying.  Matches real root-tip biology — tips
    // survive weeks of stress by going quiescent, reactivate later via
    // normal can_activate() when sugar returns.  Runs before the active
    // check so the newly-dormant RA is handled by the dormant branch below.
    if (active && starvation_ticks >= static_cast<uint32_t>(g.quiescence_threshold)) {
        active = false;
        starvation_ticks = 0;
    }

    if (!active) {
        // Dormant RAs do not synthesize auxin or cytokinin — matches shoot apical
        // behavior and real plant biology (hormone biosynthesis is an active process
        // gated on metabolism).  Dormant buds depend on PIN-transported auxin from
        // upstream to eventually cross the activation threshold, which is the
        // canalization-driven branching pattern real plants exhibit: laterals
        // activate along well-canalized paths, not uniformly everywhere.
        if (can_activate(g, world)) activate(g, world);
        return;
    }

    // ... (rest of the function unchanged)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[quiescence]" 2>&1 | tail -5`
Expected: both quiescence tests pass.

- [ ] **Step 5: Run all tests**

Run: `./build/botany_tests 2>&1 | tail -3`
Expected: `All tests passed`.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/tissues/root_apical.cpp tests/test_meristem.cpp
git commit -m "feat: RA quiescence — active meristem reverts to dormant under sustained stress

After quiescence_threshold ticks of consecutive sugar starvation, an active
root apical flips to dormant instead of continuing toward death.  Dormant
buds pay no maintenance; they can later reactivate via the normal
can_activate() path when conditions improve.  Matches real plant biology:
root tips survive weeks of stress in quiescent state and resume when conditions
recover, rather than dying and forcing the plant to rebuild from scratch.

Plan: docs/superpowers/plans/2026-04-19-root-shoot-balance-feedback.md"
```

---

## Task 7: Add quiescence transition to ApicalNode

**Files:**
- Modify: `src/engine/node/tissues/apical.cpp` (the `update_tissue` function, near the top)
- Modify: `tests/test_meristem.cpp`

**Rationale:** Same mechanism as Task 6, applied symmetrically to shoot apicals. Real shoot buds also go quiescent under stress (this is how dormant axillary buds survive).

- [ ] **Step 1: Write failing test — SA reverts to dormant after sustained starvation**

Add to `tests/test_meristem.cpp`:

```cpp
TEST_CASE("Active SA reverts to dormant after quiescence_threshold ticks of starvation", "[meristem][quiescence]") {
    Genome g = default_genome();
    g.quiescence_threshold = 10.0f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world;

    ApicalNode* sa = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::APICAL && !sa) sa = n.as_apical();
    });
    REQUIRE(sa != nullptr);
    REQUIRE(sa->active);

    for (int i = 0; i < 15; ++i) {
        sa->chemical(ChemicalID::Sugar) = 0.0f;
        plant.tick(world);
    }

    REQUIRE_FALSE(sa->active);
    REQUIRE(sa->starvation_ticks == 0u);
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "Active SA reverts to dormant" 2>&1 | tail -5`

Expected: test fails — SAs don't quiesce yet.

- [ ] **Step 3: Add the quiescence check to SA update_tissue**

In `src/engine/node/tissues/apical.cpp`, find the `update_tissue` function (around line 17). Currently looks like:

```cpp
void ApicalNode::update_tissue(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    if (!active) {
        if (can_activate(g, world)) activate(g, world);
        return;
    }

    produce_auxin(plant, g, world);
    // ...
```

Insert the quiescence check at the very top (after `const Genome& g = plant.genome();`, before the `if (!active)` check):

```cpp
void ApicalNode::update_tissue(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    // Quiescence: active meristem reverts to dormant after sustained sugar
    // starvation rather than dying.  Symmetric with RA quiescence.
    if (active && starvation_ticks >= static_cast<uint32_t>(g.quiescence_threshold)) {
        active = false;
        starvation_ticks = 0;
    }

    if (!active) {
        if (can_activate(g, world)) activate(g, world);
        return;
    }

    produce_auxin(plant, g, world);
    // ...
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "Active SA reverts to dormant" 2>&1 | tail -5`
Expected: test passes.

- [ ] **Step 5: Run all tests**

Run: `./build/botany_tests 2>&1 | tail -3`
Expected: `All tests passed`.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/tissues/apical.cpp tests/test_meristem.cpp
git commit -m "feat: SA quiescence — symmetric to RA quiescence

Shoot apicals now also revert to dormant under sustained sugar starvation,
matching the RA quiescence behavior added in the previous commit.
Real shoot axillary buds go quiescent under stress and reactivate when
resources return; main-axis SAs rarely hit this but distal SAs on
over-branched shoots will now survive rather than die.

Plan: docs/superpowers/plans/2026-04-19-root-shoot-balance-feedback.md"
```

---

## Task 8: Integration smoke test — long-running plant stays alive and balanced

**Files:**
- Modify: `tests/test_meristem.cpp`

**Rationale:** Unit tests for each piece of the feedback system pass. Now verify the whole thing behaves sensibly end-to-end: a plant run for 1000+ ticks should not crash, should develop both root and shoot, and neither side should have wildly run away. This is a smoke test, not a precise calibration assertion — exact node counts depend on many parameters and tuning.

- [ ] **Step 1: Write the integration smoke test**

Add to `tests/test_meristem.cpp`:

```cpp
TEST_CASE("Plant grows and stays alive for 1000 ticks with metabolic feedback", "[meristem][integration]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world;

    // Run 1000 ticks
    for (int i = 0; i < 1000; ++i) {
        plant.tick(world);
    }

    // Count node types at the end
    int shoot_count = 0;  // STEM + SA + LEAF
    int root_count = 0;   // ROOT + RA
    plant.for_each_node([&](const Node& n) {
        switch (n.type) {
            case NodeType::STEM:   if (n.parent) shoot_count++; break;  // exclude seed
            case NodeType::APICAL:      shoot_count++; break;
            case NodeType::LEAF:        shoot_count++; break;
            case NodeType::ROOT:        root_count++; break;
            case NodeType::ROOT_APICAL: root_count++; break;
        }
    });

    // Plant grew something on both sides — precise counts depend heavily on
    // tuning but both should have at least multiple nodes at 1000 ticks.
    REQUIRE(shoot_count > 5);
    REQUIRE(root_count > 5);

    // Neither side should have run away absurdly.  Before the metabolic
    // feedback a shoot-runaway plant would hit 500+ SAs at this tick count,
    // while roots stayed under 100.  With feedback the ratio should stay
    // in a biologically plausible range (broadly 0.1:1 to 20:1).
    float ratio = static_cast<float>(root_count) / static_cast<float>(shoot_count);
    REQUIRE(ratio > 0.05f);
    REQUIRE(ratio < 50.0f);
}
```

- [ ] **Step 2: Run the test**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[integration]" 2>&1 | tail -5`
Expected: test passes. The assertions are loose — this is a smoke test. If it fails, the plant has a crashing bug or one side is completely dead; investigate the log output and re-tune parameters as noted in the plan doc.

- [ ] **Step 3: Run the full test suite**

Run: `./build/botany_tests 2>&1 | tail -3`
Expected: `All tests passed`. No regressions.

- [ ] **Step 4: Commit**

```bash
git add tests/test_meristem.cpp
git commit -m "test: integration smoke test for 1000-tick plant with metabolic feedback

Verifies that a plant with the new feedback system runs stably for 1000
ticks, develops both root and shoot (>5 nodes each), and stays in a
plausible root:shoot ratio range.  Loose assertions since exact counts
depend on parameter tuning; this catches regressions like crashes,
dead-plant scenarios, and extreme runaway behaviors.

Plan: docs/superpowers/plans/2026-04-19-root-shoot-balance-feedback.md"
```

---

## Post-Implementation: Manual Verification

After all tasks complete, run a real simulation with the default genome and observe:

- [ ] Run a ~1000-tick sim via the realtime app:

```bash
./build/botany_realtime
```

Observe the plant in the viewer. Expected:
- Plant grows both root and shoot.
- Shoot branching does not explode to 500+ SAs within 500 ticks (the previous runaway pattern).
- Water distribution (from the old analysis scripts) should show leaves with non-zero wfrac and shoot with substantially more than 3% of plant water.
- No conservation warnings in the phloem/xylem log SUMMARY rows.

If numbers look off, parameter tuning is expected — start with `cytokinin_sugar_half_saturation` and `quiescence_threshold` in the genome defaults. Refer to [docs/long-term-plan/milestone-2/root-shoot-balance.md](../../long-term-plan/milestone-2/root-shoot-balance.md) "Reference values" section for baseline numbers from prior runs.

---

## Self-Review

**Spec coverage:** Every point from the design doc's "Converged design" section maps to a task:
- `metabolic_factor` helper → Task 2
- SA auxin upgrade → Task 3
- RA auxin gating → Task 4
- RA cytokinin gating → Task 5
- Quiescence RA → Task 6
- Quiescence SA → Task 7
- New genome params → Task 1 (prerequisite)
- Integration smoke test → Task 8

**Types:** `metabolic_factor` signature matches across all call sites. `g.quiescence_threshold` is `float`, cast to `uint32_t` at comparison with `starvation_ticks` (which is `uint32_t`). Consistent everywhere.

**Placeholders:** None — every step has the actual code to write, actual commands to run, actual expected output.
