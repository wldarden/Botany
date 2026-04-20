# Demand-Driven Phloem Implementation Plan

> **SUPERSEDED 2026-04-19.** This plan would have replaced the phloem Jacobi with a global demand-driven allocator — a tactical fix that masked the underlying structural problem (no separation between a node's parenchyma, phloem, and xylem pools). Do NOT implement this plan.
>
> The structural rewrite is specified in [`../specs/2026-04-19-compartmented-vascular-model-design.md`](../specs/2026-04-19-compartmented-vascular-model-design.md), which introduces explicit LocalEnv + TransportPool compartments with sub-stepped Münch pressure-flow as the vascular algorithm. Kept here for history only.
>
> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the pairwise-Jacobi pressure-flow in `phloem_resolve` with a demand-driven source-to-sink allocation, mirroring the algorithm already used by `xylem_resolve`. Fixes sugar delivery to apex meristems on long (30+ internode) shoot chains.

**Architecture:** Single-pass proportional allocation. Classify every node as source (has surplus sugar above reserve) or sink (needs sugar for growth/maintenance). Aggregate plant-wide totals. `delivered = min(total_supply, total_demand)`. Sources lose and sinks gain proportional shares. Mass-conservative by construction. Same pattern as `xylem_resolve` so the algorithm is already validated in this codebase.

**Tech Stack:** C++17, Catch2 tests, CMake build (`/usr/local/bin/cmake --build build`, tests via `./build/botany_tests`).

**Design context:** See the bottom of [docs/long-term-plan/milestone-2/root-shoot-balance.md](../../long-term-plan/milestone-2/root-shoot-balance.md) — this solves the phloem attenuation problem observed in the tick-2596 run where the primary SA at y=0.764 had sugar=0 for 1900+ ticks despite the seed having 0.23 g and leaves producing 0.02 g/tick.

---

## File Structure

**Modify:**
- `src/engine/vascular.cpp` — replace `phloem_resolve`'s Phase 2 (Jacobi) and Phase 3 (meristem unloading) with a single demand-driven allocation step. Keep Phase 1 (leaf loading) as-is since it's a biologically meaningful local source-to-pipe transfer.
- `src/engine/world_params.h` — `phloem_iterations` is no longer used by the new algorithm; leave the field in place but update its comment (treat as deprecated/unused).
- `tests/test_vascularization.cpp` — add test for long-chain delivery.
- `tests/test_meristem.cpp` — update the 1000-tick integration smoke test to also assert the primary SA has non-zero sugar at plant height > 0.5 dm.

No new files. All changes are localized to `phloem_resolve` and tests.

---

## Task 1: Add long-chain phloem delivery test (TDD)

**Files:**
- Modify: `/Users/wldarden/learning/botany/tests/test_vascularization.cpp`

**Rationale:** Before touching the implementation, write the test that captures the failure mode we observed. A plant with a long shoot chain and an active apex at the tip must deliver sugar to the apex within a few ticks, regardless of chain length. The current Jacobi Münch attenuates over ~30 hops; the demand-driven version should not.

The test constructs a synthetic 15-stem chain (seed → stem → stem → ... → SA) with sugar stored at the seed and zero elsewhere. After one tick, the apex SA should have measurably non-zero sugar.

- [ ] **Step 1: Add the new test to test_vascularization.cpp**

Add at the end of the file, before the final closing namespace/brace:

```cpp
TEST_CASE("Demand-driven phloem delivers sugar to apex across long shoot chain", "[vascular][phloem][demand]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world;

    // Build a synthetic 15-stem chain above the seed.  Give each stem enough
    // radius to count as a vascular conduit.  The primary SA is already the
    // seed's shoot child from the Plant constructor; we'll reparent it to the
    // tip of the chain so it acts as the apex.
    Node* tip_stem = plant.seed_mut();
    for (int i = 0; i < 15; ++i) {
        Node* stem = plant.create_node(NodeType::STEM,
            glm::vec3(0.0f, 0.05f * (i + 1), 0.0f), 0.015f);
        tip_stem->add_child(stem);
        tip_stem = stem;
    }

    // Find the primary SA and reparent it to the tip of the chain.
    ApicalNode* apex = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (auto sa = n.as_apical(); sa && sa->is_primary) apex = sa;
    });
    REQUIRE(apex != nullptr);
    // Remove apex from its current parent and attach to the chain tip.
    if (apex->parent) {
        auto& siblings = apex->parent->children;
        siblings.erase(std::remove(siblings.begin(), siblings.end(),
                                   static_cast<Node*>(apex)),
                       siblings.end());
    }
    apex->parent = nullptr;
    tip_stem->add_child(apex);

    // Zero all sugar except the seed.  Put a healthy pile at the seed so we
    // can detect any flow reaching the tip.
    plant.for_each_node_mut([&](Node& n) {
        n.chemical(ChemicalID::Sugar) = 0.0f;
    });
    plant.seed_mut()->chemical(ChemicalID::Sugar) = 10.0f;

    float apex_sugar_before = apex->chemical(ChemicalID::Sugar);
    REQUIRE(apex_sugar_before == 0.0f);

    // Single tick — phloem_resolve runs once.  Apex should receive a
    // measurable amount (> 0.001 g, well above float noise).  The precise
    // amount depends on proportional allocation, but anything non-trivial
    // demonstrates that the 15-hop distance did not attenuate delivery.
    plant.tick(world);

    float apex_sugar_after = apex->chemical(ChemicalID::Sugar);
    REQUIRE(apex_sugar_after > 0.001f);
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "Demand-driven phloem delivers sugar to apex across long shoot chain" 2>&1 | tail -10`

Expected: the test FAILS because the current Jacobi phloem with 3 iterations can only propagate sugar 3 hops per tick. After one tick, the apex at the end of a 15-stem chain has essentially no sugar. The `REQUIRE(apex_sugar_after > 0.001f)` assertion fails.

- [ ] **Step 3: Commit the failing test**

```bash
cd /Users/wldarden/learning/botany
git add tests/test_vascularization.cpp
git commit -m "test: failing case for long-chain phloem delivery

Captures the observed regression where the primary SA on a tall plant
stops receiving sugar once the shoot chain grows past ~10 internodes.
Current Jacobi Münch phloem attenuates over each pairwise-equalized
edge; demand-driven rewrite (next commit) is expected to make this pass.

Plan: docs/superpowers/plans/2026-04-19-demand-driven-phloem.md"
```

---

## Task 2: Replace Phase 2 + Phase 3 of `phloem_resolve` with demand-driven allocation

**Files:**
- Modify: `/Users/wldarden/learning/botany/src/engine/vascular.cpp`

**Rationale:** The new algorithm mirrors `xylem_resolve`'s demand-driven pattern — already validated in this codebase, same mass-conservation guarantee. Phase 1 (leaf loading) stays because leaves loading sugar into their parent's sieve tubes is a biologically meaningful local step. Phase 2's Jacobi body and Phase 3's meristem unloading both get replaced by one plant-wide source-to-sink allocation.

**Classification rules:**

| node type & state | source (supply) | sink (demand) |
|---|---|---|
| LEAF | already loaded into parent in Phase 1 — no direct role here | — |
| STEM / ROOT with sugar > reserve | `sugar − reserve` | — |
| STEM / ROOT with sugar < reserve | — | `reserve − sugar` (refill to maintenance buffer) |
| Seed (STEM, no parent) with sugar > reserve | `sugar − reserve` | — |
| Active APICAL / ROOT_APICAL | — | `sink_fraction × cap − sugar` (fill toward per-tick max) |
| Dormant APICAL / ROOT_APICAL | — | — |

`reserve = g.phloem_reserve_fraction × sugar_cap(node, g)`
`sink_fraction = g.meristem_sink_fraction` (existing genome param, default 0.05)
`cap = sugar_cap(node, g)`

**Mass conservation proof:** For every node we either add `supply[i] × (delivered/total_supply)` or subtract `demand[j] × (delivered/total_demand)`. By construction `delivered ≤ total_supply` and `delivered ≤ total_demand`, so no supply can exceed `supply[i]` and no sink can exceed `demand[j]`. Total deducted = `delivered = ` total added. Identical conservation invariant as `xylem_resolve`.

- [ ] **Step 1: Read current phloem_resolve**

Run: `grep -n "Phase 1\|Phase 2\|Phase 3\|^void phloem_resolve" /Users/wldarden/learning/botany/src/engine/vascular.cpp`

Expected output shows:
- `void phloem_resolve(Plant& plant, const Genome& g, const WorldParams& world) {` — the function declaration
- `// ── Phase 1: Leaf loading pass ──` around line 229 — KEEP THIS ENTIRE BLOCK
- `// ── Phase 2: Jacobi simultaneous pipe-network resolution ──` around line 281 — REPLACE
- `// ── Phase 3: Meristem unloading pass ──` around line 424 — REMOVE (absorbed into Phase 2's replacement)

- [ ] **Step 2: Replace Phase 2 + Phase 3 with the demand-driven allocation**

Open `src/engine/vascular.cpp`. Find the block that starts with the comment `// ── Phase 2: Jacobi simultaneous pipe-network resolution ─────────────────` (around line 281). That comment marks the beginning of the Jacobi loop. Find the closing `}` of the `for (uint32_t iter = ...)` loop — it's around line 422.

Then find Phase 3's block starting with `// ── Phase 3: Meristem unloading pass ─────────────────────────────────────` (around line 424) and ending at the blank line before `// ── Phloem debug log ──` (around line 472).

Delete the entire range from the Phase 2 comment through the end of Phase 3 (so `// ── Phase 2:` down to and including the closing brace of the Phase 3 for-loop). Replace it with:

```cpp
    // ── Phase 2: Demand-driven source-to-sink allocation ──────────────────────
    // Replaces the previous pairwise-Jacobi Münch resolution.  Pairwise
    // equilibration attenuates sugar signals exponentially along long shoot
    // chains — observed in a 2596-tick run where the primary SA at the tip
    // of a 30-stem chain had sugar=0 for 1900+ ticks despite the seed holding
    // 0.23 g.  The demand-driven pattern (mirroring xylem_resolve) eliminates
    // hop-by-hop attenuation: sinks pull sugar directly from the plant-wide
    // supply pool in proportion to their demand.
    //
    // Classification (computed once per tick):
    //   SOURCE (supply):
    //     - STEM / ROOT / seed with sugar > phloem_reserve — surplus above the
    //       parenchyma reserve is available to flow to sinks.
    //   SINK (demand):
    //     - Active APICAL / ROOT_APICAL wanting to fill to meristem_sink_fraction
    //       × cap per tick.  Dormant meristems have zero demand.
    //     - STEM / ROOT / seed with sugar < reserve — wants to refill the
    //       parenchyma buffer for maintenance.
    //
    // LEAF loading was already handled in Phase 1 above; by this point leaf
    // surplus is in the parent conduit and participates here as STEM supply.
    //
    // Mass conservation: Σ source_deduct = delivered = Σ sink_add, by
    // construction of proportional scales.  Verified every tick by the SUMMARY
    // row's conservation_error column in phloem_log.csv.
    std::vector<float> supply(N, 0.0f);
    std::vector<float> demand(N, 0.0f);

    for (int i = 0; i < N; ++i) {
        Node& n = *flat[i].node;
        float cap = sugar_cap(n, g);
        if (cap <= 1e-8f) continue;
        float reserve = cap * g.phloem_reserve_fraction;
        float sug = n.chemical(ChemicalID::Sugar);

        if (n.type == NodeType::APICAL) {
            // Active shoot apicals are sinks; dormant ones are neither.
            if (n.as_apical()->active) {
                float target = cap * g.meristem_sink_fraction;
                demand[i] = std::max(0.0f, target - sug);
            }
        } else if (n.type == NodeType::ROOT_APICAL) {
            if (n.as_root_apical()->active) {
                float target = cap * g.meristem_sink_fraction;
                demand[i] = std::max(0.0f, target - sug);
            }
        } else if (n.type == NodeType::STEM || n.type == NodeType::ROOT) {
            // Conduit / storage tissue: surplus above reserve flows, deficit
            // below reserve pulls.  Seed (STEM with no parent) uses the same
            // rule — its reserve is sugar_cap × phloem_reserve_fraction which
            // scales with seed_sugar at plant creation.
            if (sug > reserve) supply[i] = sug - reserve;
            else               demand[i] = reserve - sug;
        }
        // LEAF nodes: participated in Phase 1.  They keep their remaining
        // sugar for their own use (photosynthesis-triggered growth, respiration).
    }

    float total_supply = 0.0f, total_demand = 0.0f;
    for (int i = 0; i < N; ++i) {
        total_supply += supply[i];
        total_demand += demand[i];
    }
    float delivered = std::min(total_supply, total_demand);

    if (delivered > 1e-12f) {
        float supply_scale = (total_supply > 1e-12f) ? delivered / total_supply : 0.0f;
        float demand_scale = (total_demand > 1e-12f) ? delivered / total_demand : 0.0f;

        for (int i = 0; i < N; ++i) {
            Node& n = *flat[i].node;
            if (supply[i] > 0.0f) {
                float deducted = supply[i] * supply_scale;
                n.chemical(ChemicalID::Sugar) -= deducted;
                if (logging) log_flow_out[i] += deducted;
            }
            if (demand[i] > 0.0f) {
                float added = demand[i] * demand_scale;
                n.chemical(ChemicalID::Sugar) += added;
                if (logging) log_flow_in[i] += added;
                // Meristems track "unloaded" separately so the SUMMARY log
                // still distinguishes meristem-sink delivery from conduit
                // refill.  Store under the meristem's parent id, matching the
                // old Phase 3 log convention.
                if ((n.type == NodeType::APICAL || n.type == NodeType::ROOT_APICAL)
                    && flat[i].parent_idx >= 0) {
                    log_unloaded[flat[i].parent_idx] += added;
                }
            }
        }
    }
```

- [ ] **Step 3: Build, verify compile**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | grep -E "error:" | head -5`

Expected: no compile errors. If you see any, re-read the diff in context — the most likely culprit is a missed `}` from the old Phase 2/3 or a stale reference to removed variables (`iter`, `iter_delta`, `pressure`, `edge_final_pressure` outside the logging block — that's OK if the logging block still compiles; just be careful).

The `edge_final_pressure` variable is still allocated in the logging setup block above; leave that alone — the edge-log still wants to write a final pressure column, but under the new algorithm pressures are not meaningful. The field can safely be written as 0 when `logging` is true; no code change needed because `edge_final_pressure[i]` is already initialized to 0.0f and never assigned elsewhere in the replacement code.

- [ ] **Step 4: Run the failing test from Task 1 to verify it now passes**

Run: `./build/botany_tests "Demand-driven phloem delivers sugar to apex across long shoot chain" 2>&1 | tail -10`

Expected: PASS. The apex on a 15-stem chain now receives non-trivial sugar in a single tick.

- [ ] **Step 5: Run the full test suite to verify no regressions**

Run: `./build/botany_tests 2>&1 | tail -3`

Expected: `All tests passed`. A few existing tests may need to be inspected if they assumed specific per-edge Jacobi behavior, but most should pass — they tested the mass-conservation guarantee and end-state sugar distribution, both of which are preserved.

- [ ] **Step 6: Commit**

```bash
cd /Users/wldarden/learning/botany
git add src/engine/vascular.cpp
git commit -m "feat: replace phloem Jacobi with demand-driven source-to-sink allocation

Phase 2 of phloem_resolve was pairwise-Jacobi Münch pressure-flow,
which attenuates sugar signals exponentially along long shoot chains.
A tick-2596 sim showed the primary SA at y=0.764 with sugar=0 for
1900+ ticks despite the seed holding 0.23 g and leaves producing
0.02 g/tick — sugar simply couldn't propagate through 30 internodes.

The rewrite mirrors xylem_resolve's demand-driven pattern:
  - Sources: seed + stems + roots with sugar above phloem_reserve
  - Sinks: active meristems (target = sink_fraction × cap) and
    under-reserve conduit nodes (target = refill to reserve)
  - delivered = min(total_supply, total_demand)
  - Sources lose and sinks gain proportional shares.  Mass-conservative
    by construction (same invariant as xylem_resolve).

Phase 1 (leaf loading into parent stem) stays — it's a biologically
meaningful local source-to-pipe transfer with its own gradient/
permeability model.  Phase 3 (meristem unloading) is removed; meristems
are direct sinks in the new allocation.

world.phloem_iterations is no longer used by the algorithm; field is
kept for backward compatibility and potential future use.

Plan: docs/superpowers/plans/2026-04-19-demand-driven-phloem.md"
```

---

## Task 3: Mark `phloem_iterations` as deprecated/unused in comment

**Files:**
- Modify: `/Users/wldarden/learning/botany/src/engine/world_params.h`

**Rationale:** Leave the field in place (external callers and tests may reference it) but update its comment so the next engineer doesn't wonder why it has no effect. Trivial change.

- [ ] **Step 1: Update the comment**

In `src/engine/world_params.h`, find the `phloem_iterations` field (grep for `phloem_iterations`). Replace its existing `// ...` comment with:

```cpp
    uint32_t phloem_iterations = 3;        // UNUSED since demand-driven phloem rewrite (2026-04-19).
                                           //   Kept for backward-compatible serialization; the new
                                           //   phloem_resolve is single-pass and does not iterate.
```

(The default value 3 is preserved.)

- [ ] **Step 2: Build and run full test suite**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests 2>&1 | tail -3`

Expected: build clean, all tests pass.

- [ ] **Step 3: Commit**

```bash
cd /Users/wldarden/learning/botany
git add src/engine/world_params.h
git commit -m "docs: mark phloem_iterations unused after demand-driven rewrite

Field retained for backward compatibility; the new phloem_resolve
is single-pass and has no inner iteration loop.

Plan: docs/superpowers/plans/2026-04-19-demand-driven-phloem.md"
```

---

## Task 4: Strengthen the 1000-tick integration test

**Files:**
- Modify: `/Users/wldarden/learning/botany/tests/test_meristem.cpp`

**Rationale:** The existing integration smoke test checks plant growth broadly. Add two focused assertions: (1) the primary SA still has non-zero sugar at tick 1000, (2) the primary SA's height is meaningfully above the seed. These specifically regress-test the phloem delivery bug we just fixed.

- [ ] **Step 1: Find the existing integration test**

Run: `grep -n "Plant grows and stays alive for 1000 ticks" /Users/wldarden/learning/botany/tests/test_meristem.cpp`

Expected: shows the line containing the TEST_CASE declaration (around line 996 as of commit e409fc0).

- [ ] **Step 2: Extend the existing test with the new assertions**

Open `tests/test_meristem.cpp` at the test. The current test ends with:

```cpp
    float ratio = static_cast<float>(root_count) / static_cast<float>(shoot_count);
    REQUIRE(ratio > 0.05f);
    REQUIRE(ratio < 50.0f);
}
```

Replace that closing block with:

```cpp
    float ratio = static_cast<float>(root_count) / static_cast<float>(shoot_count);
    REQUIRE(ratio > 0.05f);
    REQUIRE(ratio < 50.0f);

    // Phloem delivery regression check: the primary SA must still have sugar
    // at tick 1000 and have grown above the seed.  This specifically catches
    // the long-chain phloem attenuation bug fixed by the demand-driven
    // rewrite (2026-04-19).
    ApicalNode* primary_sa = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (auto sa = n.as_apical(); sa && sa->is_primary) primary_sa = sa;
    });
    REQUIRE(primary_sa != nullptr);
    REQUIRE(primary_sa->position.y > 0.1f);      // grew above the seed
    REQUIRE(primary_sa->chemical(ChemicalID::Sugar) > 0.0f); // still getting sugar
}
```

- [ ] **Step 3: Run the test to verify it passes**

Run: `cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[integration]" 2>&1 | tail -5`

Expected: PASS.

- [ ] **Step 4: Run the full suite to confirm no regressions**

Run: `./build/botany_tests 2>&1 | tail -3`

Expected: `All tests passed`.

- [ ] **Step 5: Commit**

```bash
cd /Users/wldarden/learning/botany
git add tests/test_meristem.cpp
git commit -m "test: integration smoke test now asserts primary SA gets sugar

Extends the 1000-tick integration test with two phloem-delivery
regression checks: primary SA height > 0.1 dm and sugar > 0.
Directly catches the long-chain attenuation bug that motivated the
demand-driven phloem rewrite.

Plan: docs/superpowers/plans/2026-04-19-demand-driven-phloem.md"
```

---

## Post-Implementation: Manual Verification

After all tasks complete, build once more and run the realtime sim for 2000+ ticks to verify real-world behavior:

- [ ] Build and launch:

```bash
cd /Users/wldarden/learning/botany
/usr/local/bin/cmake --build build
./build/botany_realtime
```

Observe:
- Plant grows continuously past tick 2500 (no death spiral).
- Conservation in `debug/phloem_log.csv` SUMMARY rows stays ≈ 0 every tick.
- The primary SA at the shoot tip is continuously receiving sugar and making new internodes.
- Distal leaves on the shoot chain are healthy (not drought-stressed due to shoot drought).

If any of these fail, the most likely culprits are:
- A missed node type in the sink/source classification — add it to the switch.
- Demand-scale or supply-scale computation edge case at empty plant or all-dormant plant — already guarded with `delivered > 1e-12f`, but verify.
- `phloem_reserve_fraction` or `meristem_sink_fraction` needs tuning for the new algorithm's dynamics. Tune in genome defaults; not an algorithmic fix.

---

## Self-Review

**Spec coverage:**
- Demand-driven replacement of Phase 2 + Phase 3 → Task 2.
- Phase 1 (leaf loading) preserved → explicitly noted in Task 2.
- Mass-conservation guarantee → by construction in Task 2, verified by existing SUMMARY conservation check.
- Long-chain delivery test → Task 1.
- Integration test strengthened → Task 4.
- Deprecated param cleanup → Task 3.

**Placeholder scan:** none found. Every code step has complete code.

**Type consistency:** `supply` and `demand` vectors are `std::vector<float>` with size `N`, indexed by node position in the flat array (same convention as the existing xylem implementation). `delivered`, `supply_scale`, `demand_scale` are `float`. Matches `xylem_resolve`'s style exactly.

**Remaining risk:** existing tests that exercised the Jacobi-specific pressure dynamics (if any) may need to be adjusted if they asserted behavior that's no longer meaningful. If a test fails after Task 2, read its assertions carefully — most likely it was testing the OLD mechanism's edge case rather than a correctness property, and the right fix is to update the assertion to test the new algorithm's equivalent behavior.
