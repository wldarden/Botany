# Hormone Biology at Tree Scale Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restructure auxin and cytokinin so growth and long-distance signaling work biologically correctly from seedling (1 dm) to mature tree (10 m), ending the scaling failure where shoot-derived auxin must reach the root apex.

**Architecture:** (1) Decouple cytokinin production at root apicals from local auxin — make it a sugar+water metabolic signal. (2) Swap the root apical elongation gate from auxin to local cytokinin. (3) Restore a small root-tip auxin self-production floor for PIN recycling and lateral root initiation. (4) Retune activation and production constants. Auxin remains a short-range local signal (apical dominance, canalization). Cytokinin becomes the root→shoot "healthy" signal carried by xylem bulk flow. Sugar remains the shoot→root energy channel carried by phloem bulk flow.

**Tech Stack:** C++17, Catch2 tests, cmake. Engine library `botany_engine` in `src/engine/`. Reference spec: [docs/superpowers/specs/2026-04-20-hormone-biology-tree-scale-design.md](../specs/2026-04-20-hormone-biology-tree-scale-design.md).

**Build command:** `/usr/local/bin/cmake --build build`
**Test command:** `./build/botany_tests`
**Scaling sim:** `./build/botany_headless <ticks> /tmp/rec.bin --chain-profile /tmp/profile.csv`

---

## Task 1: Decouple cytokinin production at RA from local auxin

Currently the root apical produces cytokinin as `rate × local_auxin × mf_cyto`. This couples CK supply to a signal (auxin) that cannot scale. We change it to `rate × mf_cyto` so CK production depends only on local sugar+water metabolic state.

**Files:**
- Test: `tests/test_hormone.cpp` (add new test case)
- Modify: `src/engine/node/tissues/root_apical.cpp:73`

- [ ] **Step 1: Add failing test for decoupled CK production**

Open `tests/test_hormone.cpp` and add this test case at the end of the file (before the closing brace of the namespace or after the last TEST_CASE):

```cpp
TEST_CASE("Cytokinin: RA produces CK without any local auxin", "[hormone]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world();

    RootApicalNode* ra = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::ROOT_APICAL && !ra) ra = n.as_root_apical();
    });
    REQUIRE(ra != nullptr);
    REQUIRE(ra->active);

    // Prime sugar + water, zero local auxin.
    ra->local().chemical(ChemicalID::Sugar) = 1.0f;
    ra->local().chemical(ChemicalID::Water) = 1.0f;
    ra->local().chemical(ChemicalID::Auxin) = 0.0f;
    ra->tick_cytokinin_produced = 0.0f;

    // Tick the RA directly (bypasses vascular so we observe update_tissue only).
    ra->tick(plant, world);

    // With decoupled production, CK should be produced regardless of auxin.
    REQUIRE(ra->tick_cytokinin_produced > 0.0f);
}
```

- [ ] **Step 2: Build and run the test to verify it FAILS**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests "Cytokinin: RA produces CK without any local auxin"
```

Expected: FAIL (`tick_cytokinin_produced > 0.0f` with actual value 0, because production is gated by local_auxin which is 0).

- [ ] **Step 3: Remove the auxin factor from CK production**

In `src/engine/node/tissues/root_apical.cpp`, find the cytokinin production block around line 70-75:

```cpp
    float mf_cyto = metabolic_factor(
        local().chemical(ChemicalID::Sugar), g.cytokinin_sugar_half_saturation, 0.05f,
        local().chemical(ChemicalID::Water), g.cytokinin_water_half_saturation, 0.05f);
    float cyto_produced = g.root_cytokinin_production_rate * local().chemical(ChemicalID::Auxin) * mf_cyto;
    local().chemical(ChemicalID::Cytokinin) += cyto_produced;
    tick_cytokinin_produced += cyto_produced;
```

Change the `cyto_produced` line (drop the `local().chemical(ChemicalID::Auxin) *` factor) and update the comment on the block above. Replace the entire block with:

```cpp
    // Cytokinin: floor 0.05 (smaller than auxin's 0.1) — CK has less
    // conjugate-pool buffering than auxin, so its response to sugar is sharper.
    // CK production is driven purely by local metabolic state (sugar + water).
    // It is NOT gated on local auxin — in real root tips, CK synthesis reflects
    // "this tip is metabolically healthy," independent of the distant shoot.
    // (Prior versions gated on local_auxin, which broke for tall plants where
    // shoot-derived auxin never reaches the RA; see hormone-biology spec.)
    float mf_cyto = metabolic_factor(
        local().chemical(ChemicalID::Sugar), g.cytokinin_sugar_half_saturation, 0.05f,
        local().chemical(ChemicalID::Water), g.cytokinin_water_half_saturation, 0.05f);
    float cyto_produced = g.root_cytokinin_production_rate * mf_cyto;
    local().chemical(ChemicalID::Cytokinin) += cyto_produced;
    tick_cytokinin_produced += cyto_produced;
```

- [ ] **Step 4: Build and run the new test to verify it PASSES**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests "Cytokinin: RA produces CK without any local auxin"
```

Expected: PASS.

- [ ] **Step 5: Run the full test suite to verify nothing else regressed**

```bash
./build/botany_tests
```

Expected: All tests pass. If the existing test `"Cytokinin: root apical produces cytokinin"` fails because it was using `tick_cytokinin_produced > 0.0f` and now gets an unexpectedly large value, that's not a problem — it's still `> 0`. No changes needed there.

- [ ] **Step 6: Commit**

```bash
git add tests/test_hormone.cpp src/engine/node/tissues/root_apical.cpp
git commit -m "hormone: decouple RA cytokinin production from local auxin

CK now scales only with sugar+water metabolic state. Matches real root-tip
biology and removes the scaling failure where shoot-derived auxin must
reach the RA for CK to be produced. See spec 2026-04-20-hormone-biology-tree-scale.

Co-Authored-By: Claude Opus 4 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Add `root_ck_growth_floor` genome parameter and swap RA elongation gate to cytokinin

Previously RA elongation was gated on local auxin via `growth_fraction(sugar, max_cost, auxin, root_auxin_growth_threshold)`. Since auxin doesn't reach root tips in tall plants, this breaks at scale. We swap the modulator to local cytokinin — which is now produced at the RA from metabolic state (Task 1) and is also delivered to the RA's xylem-adjacent neighborhood. New Km is `root_ck_growth_floor = 0.001` (low — trace CK permits full-rate growth; zero CK stops growth).

**Files:**
- Modify: `src/engine/genome.h` (add field + default)
- Modify: `src/engine/node/tissues/root_apical.cpp:129-138`
- Test: `tests/test_cytokinin_transport.cpp` (new test case)

- [ ] **Step 1: Add failing test for CK-gated RA elongation**

Open `tests/test_cytokinin_transport.cpp` and add this test case at the end of the file:

```cpp
// -----------------------------------------------------------------------
// Test 4: RA elongation halts when local cytokinin drops to zero.
//
// The new gate ties elongation rate to local CK (not auxin). With sugar
// and water plentiful but CK zeroed each tick, the RA should not advance.
// -----------------------------------------------------------------------
TEST_CASE("cytokinin: RA elongation stops with zero local cytokinin", "[cytokinin]") {
    Genome g = cytokinin_test_genome();
    // Re-enable root growth (the test-genome zeroes it).
    g.root_growth_rate = 0.004f;
    // Disable self-production of CK so we can keep local CK at 0 artificially.
    g.root_cytokinin_production_rate = 0.0f;

    Plant plant(g, glm::vec3(0.0f));
    Node* seed = plant.seed_mut();
    RootApicalNode* ra = nullptr;
    for (Node* c : seed->children) {
        if (c->type == NodeType::ROOT_APICAL) { ra = c->as_root_apical(); break; }
    }
    REQUIRE(ra != nullptr);

    // Prime sugar/water so sugar/water gates are satisfied.
    seed->local().chemical(ChemicalID::Sugar) = 100.0f;
    ra->local().chemical(ChemicalID::Sugar)   = 50.0f;
    ra->local().chemical(ChemicalID::Water)   = 1.0f;

    // Capture starting offset length.
    float start_offset_len = glm::length(ra->offset);

    WorldParams world = cytokinin_test_world();
    // Tick N times, zeroing CK at every tick after update_tissue runs.
    for (int i = 0; i < 20; ++i) {
        plant.tick(world);
        ra->local().chemical(ChemicalID::Cytokinin) = 0.0f;
        ra->local().chemical(ChemicalID::Sugar)   = 50.0f;  // keep primed
        ra->local().chemical(ChemicalID::Water)   = 1.0f;
    }

    // With CK pinned to 0, offset should not have grown meaningfully.
    float end_offset_len = glm::length(ra->offset);
    REQUIRE(end_offset_len - start_offset_len < 1e-4f);
}
```

- [ ] **Step 2: Build and run to verify it FAILS**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests "cytokinin: RA elongation stops with zero local cytokinin"
```

Expected: FAIL — the RA still elongates because the current gate uses auxin, which accumulates from RA self-production or shoot delivery.

- [ ] **Step 3: Add `root_ck_growth_floor` to the genome struct**

In `src/engine/genome.h`, find the root-related parameters (around line 80-86, after `root_auxin_growth_threshold`). Add the new field:

```cpp
    float root_ck_growth_floor;            // Km for CK-gated root elongation — low so trace CK permits growth; zero CK stops it
```

- [ ] **Step 4: Add default value in `default_genome()`**

In `src/engine/genome.h`, find the `default_genome()` block where `root_auxin_growth_threshold` is set (around line 287). Add the new default next to it:

```cpp
        .root_ck_growth_floor = 0.001f,             // low Km — trace CK permits growth, zero stops it
```

- [ ] **Step 5: Swap the RA elongation gate from auxin to cytokinin**

In `src/engine/node/tissues/root_apical.cpp`, find `RootApicalNode::elongate` (around line 129). The current body starts:

```cpp
void RootApicalNode::elongate(const Genome& g, const WorldParams& world) {
    float max_cost = g.root_growth_rate * world.sugar_cost_root_growth;
    // Root elongation is gated by auxin just as SAM elongation is gated by
    // cytokinin — Michaelis-Menten rate modulation with Km = root_auxin_growth_threshold.
    // Local PIN-recycling production (root_tip_auxin_production_rate) gives a
    // tiny floor (self-equilibrium ≪ threshold), so shoot-delivered auxin is
    // what actually drives elongation.  Until the stem vasculature thickens
    // enough to transport meaningful auxin downward, root growth is naturally
    // capped — preventing the runaway root/starved-SAM regime.
    float gf = growth_fraction(local().chemical(ChemicalID::Sugar), max_cost,
                               local().chemical(ChemicalID::Auxin), g.root_auxin_growth_threshold);
```

Replace the comment block and the `gf` line with:

```cpp
void RootApicalNode::elongate(const Genome& g, const WorldParams& world) {
    float max_cost = g.root_growth_rate * world.sugar_cost_root_growth;
    // Root elongation is gated by local cytokinin (not auxin).  CK at the RA
    // represents "this tip is metabolically healthy" — it's produced from the
    // RA's own sugar+water (see update_tissue).  Using CK as the modulator
    // keeps the gate scale-free: a 10m tree's deep RAs still produce their
    // own CK and gate their own growth, with no dependency on distant-shoot
    // signals that don't scale (see hormone-biology spec).  Sugar remains the
    // actual rate limiter — CK only permits growth if sugar is also available.
    float gf = growth_fraction(local().chemical(ChemicalID::Sugar), max_cost,
                               local().chemical(ChemicalID::Cytokinin), g.root_ck_growth_floor);
```

(Leave the rest of `elongate` — the water_gf check and everything below — unchanged.)

- [ ] **Step 6: Build and run the new test to verify it PASSES**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests "cytokinin: RA elongation stops with zero local cytokinin"
```

Expected: PASS.

- [ ] **Step 7: Run the full test suite**

```bash
./build/botany_tests
```

Expected: All tests pass.

- [ ] **Step 8: Commit**

```bash
git add src/engine/genome.h src/engine/node/tissues/root_apical.cpp tests/test_cytokinin_transport.cpp
git commit -m "hormone: swap RA elongation gate from auxin to local CK

Adds root_ck_growth_floor genome parameter (default 0.001). Previously
the RA elongation growth_fraction used local_auxin as the modulator,
which broke at tree scale because shoot-derived auxin cannot reach a
distant RA. The new gate uses local cytokinin, which is produced at the
RA itself from its own sugar+water metabolic state.

Co-Authored-By: Claude Opus 4 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Restore small RA auxin self-production floor

While testing the auxin gate we set `root_tip_auxin_production_rate` to 0. With the gate removed, we restore a small floor (0.002) matching the biological PIN-recycling baseline at real root tips. This provides the local auxin a RA needs for (a) its own PIN canalization activity and (b) lateral root initiation. Self-equilibrium with decay 0.12: 0.002 / 0.12 ≈ 0.017 AU, below the activation threshold we'll lower in Task 4.

**Files:**
- Modify: `src/engine/genome.h` (retune default)
- Modify: `tests/test_meristem.cpp` (remove the skip guard on the disabled test)

- [ ] **Step 1: Restore production rate to 0.002**

In `src/engine/genome.h`, find the line with `.root_tip_auxin_production_rate = 0.0f,`. Replace with:

```cpp
        .root_tip_auxin_production_rate = 0.002f,  // small PIN-recycling floor — self-equilibrium ~0.017 supports lateral root initiation without needing distant shoot auxin
```

- [ ] **Step 2: Re-enable the `"RA auxin production drops when sugar is low"` test**

In `tests/test_meristem.cpp`, find the test case `"RA auxin production drops when sugar is low"`. It currently begins with:

```cpp
TEST_CASE("RA auxin production drops when sugar is low", "[meristem][metabolic]") {
    Genome g = default_genome();
    // This test exercises the metabolic gating of RA self-production.  When
    // root_tip_auxin_production_rate is temporarily set to 0 (testing pure
    // external-auxin gating of root elongation), there is no self-production
    // to measure — skip the assertions rather than fail.
    if (g.root_tip_auxin_production_rate <= 0.0f) return;

    Plant plant(g, glm::vec3(0.0f));
```

Remove the 5-line comment and the `if (g.root_tip_auxin_production_rate <= 0.0f) return;` guard so the test body starts directly with `Plant plant(...);`:

```cpp
TEST_CASE("RA auxin production drops when sugar is low", "[meristem][metabolic]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
```

- [ ] **Step 3: Build and run the affected tests**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests "[meristem][metabolic]"
```

Expected: All pass. The `produced_starved < produced_fed * 0.3f` and `produced_starved > 0.0f` assertions should hold because production at sugar=0 uses the `0.1f` floor in `metabolic_factor` and at sugar=1 uses close-to-1 mf.

- [ ] **Step 4: Run the full suite**

```bash
./build/botany_tests
```

Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/engine/genome.h tests/test_meristem.cpp
git commit -m "hormone: restore RA auxin self-production floor (0.002)

Auxin at the root tip is biologically real (PIN recycling, lateral root
initiation). We disabled production to 0 while testing the external-auxin
gate; now that the gate is gone we restore a small floor. Self-equilibrium
with decay 0.12 is ~0.017 AU — below activation threshold so it doesn't
spuriously activate laterals, but enough to support the RA's own canalization.

Co-Authored-By: Claude Opus 4 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Retune activation + production thresholds for the new architecture

Three parameters now need fresh defaults:
- `root_auxin_activation_threshold`: currently 0.05 (above RA self-equilibrium of 0.017). Deep lateral RAs in a tall tree can't activate without shoot-delivered auxin, which fails at scale. Lower to 0.01 so laterals activate on their own auxin floor.
- `root_cytokinin_production_rate`: currently 0.15. Before Task 1, this was multiplied by `local_auxin ~0.05`, so effective rate was ~0.0075/tick/RA. Post-Task 1 we just use the raw rate, so 0.15 is too high and 0.0075 is too low (need to keep CK production meaningful for xylem transport). Set to 0.5.
- `root_cytokinin_inhibition_threshold`: currently 0.15, dead code before because CK was 0 everywhere. Now that CK flows, this threshold suppresses runaway lateral root branching in CK-rich regions. Keep at 0.15 for now but document the validation step.

**Files:**
- Modify: `src/engine/genome.h`

- [ ] **Step 1: Retune three defaults**

In `src/engine/genome.h`, inside `default_genome()`:

Find and replace `.root_auxin_activation_threshold = 0.05f,` → `.root_auxin_activation_threshold = 0.01f,  // below RA self-eq (~0.017) so deep laterals can activate without shoot-delivered auxin`

Find and replace `.root_cytokinin_production_rate = 0.15f,` → `.root_cytokinin_production_rate = 0.5f,  // raw rate (CK no longer gated on local auxin post-Task 1); set so ~20 active RAs feed xylem enough CK to drive SA growth`

(Leave `.root_cytokinin_inhibition_threshold = 0.15f,` unchanged. Just confirm it's present.)

- [ ] **Step 2: Build and run the full test suite**

```bash
/usr/local/bin/cmake --build build && ./build/botany_tests
```

Expected: All tests pass. (Thresholds are not directly tested.)

- [ ] **Step 3: Commit**

```bash
git add src/engine/genome.h
git commit -m "hormone: retune activation and CK production rates for new architecture

- root_auxin_activation_threshold: 0.05 -> 0.01 (below RA self-eq ~0.017
  so deep lateral RAs can activate on their own auxin floor)
- root_cytokinin_production_rate: 0.15 -> 0.5 (was multiplied by local_auxin
  before decoupling; raw rate needed now for meaningful xylem CK flux)
- root_cytokinin_inhibition_threshold: unchanged at 0.15; now engages for
  the first time because CK actually flows post-decoupling

Co-Authored-By: Claude Opus 4 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Seedling-scale validation checkpoint

Run a 1000-tick headless sim with chain-profile output and verify the seedling-scale success criteria from the spec §5.1.

**Files:** No code changes. Retune inline only if criteria fail.

- [ ] **Step 1: Run the seedling sim**

```bash
/usr/local/bin/cmake --build build \
  && ./build/botany_headless 1000 /tmp/seedling.bin --chain-profile /tmp/seedling_profile.csv
```

Expected stdout: `Done. Final node count: <some number between 100 and 400>`.

- [ ] **Step 2: Verify primary SA has meaningful local cytokinin**

```bash
./build/botany_dump /tmp/seedling.bin 999 2>&1 | grep "id=1 " | head -1
```

Expected: The output line contains `cyto=<value>` where `<value>` is ≥ 0.001 (pass). If `cyto` is near zero, tuning must go back to Task 4 — raise `root_cytokinin_production_rate` until SA CK is ≥ `cytokinin_growth_threshold` (0.005) × 2.

- [ ] **Step 3: Verify primary RA is actively elongating**

```bash
./build/botany_dump /tmp/seedling.bin 999 2>&1 | grep "type=ROOT_APICAL" | grep "is_primary=1" | head -1
```

(If the dump format doesn't include `is_primary`, use the root chain profile instead:)

```bash
awk -F, '$1=="root"' /tmp/seedling_profile.csv | tail -1
```

Expected: Last row is type `RA` at depth ≥ 15 (primary chain is growing). If depth ≤ 3, the root grew only a few internodes — raise `root_growth_rate` or investigate sugar delivery.

- [ ] **Step 4: Verify ≥ 50% of RAs are active**

```bash
./build/botany_dump /tmp/seedling.bin 999 2>&1 | grep "type=ROOT_APICAL" | wc -l
```

Record total N. Then:

```bash
./build/botany_dump /tmp/seedling.bin 999 2>&1 | grep "type=ROOT_APICAL" | grep -c "active=1"
```

Expected: active count ≥ 0.5 × N. If below, lowering `root_auxin_activation_threshold` further (e.g., to 0.005) may help.

- [ ] **Step 5: Verify visible auxin gradient along shoot chain**

```bash
awk -F, '$1=="shoot"' /tmp/seedling_profile.csv | awk -F, '{print $2, $7}'
```

Expected: Auxin at depth 1 (near seed) is ≤ 0.8 × auxin at last-depth (near SAM) — shows basipetal gradient. If the gradient is reversed or flat, transport or decay is misconfigured (out of scope for this plan; file separately).

- [ ] **Step 6: Verify visible CK gradient along root chain xylem**

```bash
awk -F, '$1=="root"' /tmp/seedling_profile.csv | awk -F, '{print $2, $11}'
```

Expected: `xyl_cytokinin` is highest at the top of the root chain (near seed) and decreases with depth, OR the total sum across root xylem is > 1e-5. If all values are < 1e-10, CK production+inject is failing; revisit Task 1 implementation.

- [ ] **Step 7: Commit only if validation passed; otherwise retune and repeat**

If all six verifications passed:

```bash
# nothing to commit — checkpoint was validation only
echo "seedling checkpoint passed"
```

If any verification failed, retune the relevant constant in `src/engine/genome.h`, then return to Step 1. Commit any retuning changes separately:

```bash
git add src/engine/genome.h
git commit -m "hormone: retune <param> to pass seedling checkpoint

Co-Authored-By: Claude Opus 4 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: Adolescent-scale validation checkpoint

Run a 10000-tick sim and verify spec §5.2 criteria. At this scale we expect ~500 nodes and ~3m height.

**Files:** No code changes unless criteria fail.

- [ ] **Step 1: Run the adolescent sim**

```bash
/usr/local/bin/cmake --build build \
  && ./build/botany_headless 10000 /tmp/adolescent.bin --chain-profile /tmp/adolescent_profile.csv
```

Expected stdout: `Done. Final node count: <number between 300 and 800>`. Takes up to a couple minutes.

- [ ] **Step 2: Verify primary SA is still alive**

```bash
./build/botany_dump /tmp/adolescent.bin 9999 2>&1 | grep "id=1 " | head -1
```

Expected: Line exists and contains `active=1`. If missing or `active=0`, the primary SAM died — investigate CK delivery or sugar starvation.

- [ ] **Step 3: Verify root growth rate is sane (not runaway, not frozen)**

```bash
awk -F, '$1=="root"' /tmp/adolescent_profile.csv | wc -l
```

Expected: 50–300 nodes in root chain (continuous growth without runaway). If > 500, roots are outpacing shoot — raise `root_cytokinin_inhibition_threshold` to slow lateral initiation. If < 20, root growth is starved — investigate sugar phloem delivery (out of scope for this plan if not fixable via CK tuning).

- [ ] **Step 4: Verify CK arrives at the SA (xylem or local)**

```bash
./build/botany_dump /tmp/adolescent.bin 9999 2>&1 | grep "id=1 " | head -1 | grep -oE "cyto=[0-9.e+-]+"
```

Expected: `cyto=<value>` ≥ 0.001. If below, the xylem CK pipeline isn't moving mass to the SAM — raise `root_cytokinin_production_rate` or consider scaling `vascular_substeps` (deferred per spec §7.1).

- [ ] **Step 5: Verify lateral-break pattern is dense near tips, sparse on trunk**

```bash
awk -F, '$1=="shoot"' /tmp/adolescent_profile.csv | awk -F, 'NR>1 {print $4}'
```

Expected: Mostly `STEM` with `SA` only at the tip. If many `SA` scattered on intermediate stems, apical dominance broke down near the trunk. This is a tuning question on `auxin_threshold` for SA activation; investigate only if the density is wildly off.

- [ ] **Step 6: Commit any tuning changes**

If the checkpoint passed with no changes, skip. Otherwise:

```bash
git add src/engine/genome.h
git commit -m "hormone: retune for adolescent checkpoint

Co-Authored-By: Claude Opus 4 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: Mature tree validation checkpoint

Run a 40000-tick sim and verify spec §5.3 criteria. Expect ~2000 nodes and ~10m. This sim is slow (potentially 10-20 minutes). Run it in the background if possible.

**Files:** No code changes unless criteria fail.

- [ ] **Step 1: Run the mature sim**

```bash
/usr/local/bin/cmake --build build \
  && ./build/botany_headless 40000 /tmp/mature.bin --chain-profile /tmp/mature_profile.csv
```

Expected stdout at completion: `Done. Final node count: <number between 1500 and 5000>`.

- [ ] **Step 2: Verify primary SA is still alive**

```bash
./build/botany_dump /tmp/mature.bin 39999 2>&1 | grep "id=1 " | head -1 | grep -oE "active=[01]"
```

Expected: `active=1`. If the primary SAM died, the plant collapsed — revisit tuning; the architecture failed the spec.

- [ ] **Step 3: Verify the tallest SAM is ≥ 5m above seed**

Positions in the sim are in dm (decimeters). 5m = 50 dm. 10m target = 100 dm.

```bash
awk -F, '$1=="shoot"' /tmp/mature_profile.csv | tail -1 | awk -F, '{print "y_dm="$5}'
```

Expected: `y_dm` ≥ 50 (5 m). Ideal target per spec is ≥ 100 (10 m). If below 30 dm (3 m), tree growth stalled — retune sugar economy or investigate (out of scope for this plan if not CK-related).

- [ ] **Step 4: Verify root tips at depth > 3m are still active**

```bash
./build/botany_dump /tmp/mature.bin 39999 2>&1 | grep "type=ROOT_APICAL" | awk -F'pos=' '{print $2}' | awk -F',' '$2 < -30.0 {print}' | wc -l
```

Expected: at least one RA with `y < -30.0` (i.e., deeper than 3m) that's active. If none, deep roots died off — investigate sugar phloem throughput at scale.

- [ ] **Step 5: Verify CK pressure at seed.xylem is non-trivial**

```bash
./build/botany_dump /tmp/mature.bin 39999 2>&1 | grep "id=0 " | head -1
```

Note: The `dump` output may not separately show xylem pool contents. If needed, read `/tmp/mature_profile.csv`:

```bash
awk -F, '$1=="shoot" && $2==0' /tmp/mature_profile.csv
```

Expected: `xyl_cytokinin` ≥ 0.01 AU in the seed's xylem row. If below, the CK xylem transport couldn't keep pace with decay at scale — consider raising `vascular_substeps` or increasing `root_cytokinin_production_rate` further (capping at 2.0 before treating it as an architecture failure).

- [ ] **Step 6: Commit final tuning**

If validation passed:

```bash
# nothing to commit — checkpoint was validation only
echo "mature checkpoint passed"
```

If tuning was needed and now passes:

```bash
git add src/engine/genome.h
git commit -m "hormone: final tuning for mature-tree checkpoint

Co-Authored-By: Claude Opus 4 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: Update `CLAUDE.md` with the new hormone architecture

The project CLAUDE.md describes the hormone transport model. Update the relevant sections so future readers understand (a) auxin is short-range by design, (b) CK is decoupled from local auxin at the RA, (c) the RA elongation gate uses CK not auxin.

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Find and update the CK description in CLAUDE.md**

Open `CLAUDE.md` and search for the block under `### Auxin, Cytokinin, Sugar in local_env` (or similar — the heading may have changed). Inside the `**Cytokinin**` bullet, replace the body with:

```
**Cytokinin** (long-distance via xylem; produced from local metabolic state):
- Persistent across ticks (not reset)
- Produced by `RootApicalNode` from its own sugar+water metabolic state — **not** gated on local auxin (decoupling was required to scale beyond small plants; see `docs/superpowers/specs/2026-04-20-hormone-biology-tree-scale-design.md`)
- Travels up in xylem via sub-stepped Jacobi pressure flow; scale-free with respect to plant height
- Drives SA elongation (via `cytokinin_growth_threshold`) and RA elongation (via `root_ck_growth_floor`)
- Root axillary activation gated by `root_cytokinin_inhibition_threshold` — high local CK suppresses lateral root initiation
- Decays by `cytokinin_decay_rate` per tick (5%) in `local()` only; no decay in xylem conduits
```

- [ ] **Step 2: Update the key design decisions section**

Find `## Key Design Decisions` and add a new bullet to the end of the list:

```
- **Hormones are scale-separated by design** — auxin is a short-range signal (decay 12%/tick, transport 1 hop/tick, effective range ~5–20 internodes), used for apical dominance, canalization, and local tropism. Cytokinin is the long-distance root→shoot signal, produced at root tips from local metabolic state (sugar + water) and carried by xylem bulk flow. Sugar is the long-distance shoot→root channel via phloem. No mechanism depends on auxin reaching a distant node — this is the fundamental scaling invariant that lets the model support plants from 1 dm seedlings to 10 m trees.
```

- [ ] **Step 3: Update tunable parameters list**

Find the `## Tuning Parameters (genome.h)` section. Search for `root_auxin_growth_threshold` and add a new line immediately after it:

```
- `root_ck_growth_floor` (0.001) - Km for cytokinin-gated RA elongation; RA grows at sugar-rate-limited pace when local CK ≥ this floor, stops when CK drops to zero (spec 2026-04-20 hormone-biology)
```

Also find `root_tip_auxin_production_rate` and update the existing line to reflect the new default:

```
- `root_tip_auxin_production_rate` (0.002) - auxin produced by root tips per tick (PIN recycling, lateral root initiation); self-equilibrium ~0.017 with default decay; too low to gate growth but sufficient for canalization
```

And update `root_auxin_activation_threshold`:

```
- `root_auxin_activation_threshold` (0.01) - minimum auxin to activate a dormant root meristem; lower than RA self-equilibrium so deep laterals can activate without shoot-delivered auxin
```

And update `root_cytokinin_production_rate`:

```
- `root_cytokinin_production_rate` (0.5) - baseline cytokinin produced per tick at active RAs, scaled by local sugar+water metabolic factor (NOT scaled by local auxin since the 2026-04-20 decoupling); drives whole-plant CK supply via xylem bulk flow
```

- [ ] **Step 4: Verify no existing guidance contradicts the new architecture**

```bash
grep -n "auxin.*root tip\|RA.*auxin\|cytokinin.*auxin" CLAUDE.md | head -10
```

Review each line and update any that claim CK production depends on local auxin, or that RA elongation depends on auxin. (Exact edits depend on what exists; keep the intent consistent with the new spec.)

- [ ] **Step 5: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: update hormone model to reflect tree-scale architecture

CK production at RA is no longer gated on local auxin. RA elongation
is gated on local CK (not auxin). Auxin is scoped to short-range
signaling. See spec docs/superpowers/specs/2026-04-20-hormone-biology-tree-scale-design.md.

Co-Authored-By: Claude Opus 4 (1M context) <noreply@anthropic.com>"
```

---

## Completion checklist

- [ ] Task 1 committed: CK production decoupled from auxin
- [ ] Task 2 committed: RA gate swapped to CK, `root_ck_growth_floor` added
- [ ] Task 3 committed: RA auxin self-production restored
- [ ] Task 4 committed: Activation + production constants retuned
- [ ] Task 5 passed: Seedling checkpoint (1000 ticks)
- [ ] Task 6 passed: Adolescent checkpoint (10000 ticks)
- [ ] Task 7 passed: Mature tree checkpoint (40000 ticks)
- [ ] Task 8 committed: CLAUDE.md updated
- [ ] Full test suite passes after final commit: `./build/botany_tests`
