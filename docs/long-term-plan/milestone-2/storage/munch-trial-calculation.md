# Münch Pressure Flow: Trial Calculation

**Purpose:** Verify that the parameter values in `munch-pressure-flow-design.md` produce physically
reasonable results on a representative small plant. Read-only analysis — no code changes.

**Plant state:** tick ~500, moderately grown seedling.

---

## Parameters Extracted from Design Doc

### WorldParams
| Parameter | Value | Unit |
|-----------|-------|------|
| `phloem_ring_thickness` (t) | 0.002 | dm |
| `max_sugar_concentration` | 300.0 | g/dm³ |
| `leaf_thickness` | 0.003 | dm |
| `base_phloem_speed` | 5.0 | dm/tick |
| `phloem_reference_radius` | 0.05 | dm |

### Genome
| Parameter | Value |
|-----------|-------|
| `phloem_osmotic_coefficient` (κ) | 1.0 |
| `phloem_conductance` | 8.0 |
| `phloem_unloading_meristem` | 0.08 |
| `phloem_unloading_leaf` | 0.10 |
| `phloem_unloading_root` | 0.008 |
| `phloem_unloading_stem` | 0.002 |

### Water capacity params (genome, for water_frac calculation)
| Parameter | Value | Unit |
|-----------|-------|------|
| `water_storage_density_stem` | 800.0 | ml/dm³ |
| `water_storage_density_leaf` | 3.0 | ml/dm² |
| `water_cap_meristem` | 1.0 | ml (fixed) |

---

## Test Plant Topology

```
graph (parent → children):

seed (r=0.03, L=0.02, sugar=2.0g)
├── stem1 (r=0.02, L=0.05, sugar=0.5g)    ← near seed, stored sugar
│   └── stem2 (r=0.015, L=0.05, sugar=0.001g)
│       ├── stem3 (r=0.015, L=0.05, sugar=0.001g)
│       │   ├── leaf1 (leaf_size=0.5, sugar=0.03g)
│       │   └── shoot_apical (r=0.005, sugar=0.002g)
│       └── leaf2 (leaf_size=0.3, sugar=0.01g)
└── root1 (r=0.015, L=0.1, sugar=0.001g)
    └── root_apical (r=0.008, sugar=0.005g)
```

---

## Step 1: Phloem Volumes

**Formula:**
- STEM/ROOT: `phloem_ring_area(r, t) × L` where `ring_area = π(2rt − t²)`
- LEAF: `leaf_size² × leaf_thickness`
- APICAL/ROOT_APICAL: `(4/3)π r³`

**`ring_area(r)` computed values:**

| r (dm) | `π(2r×0.002 − 0.000004)` | ring_area (dm²) |
|--------|--------------------------|-----------------|
| 0.005  | π × 0.000016             | 5.027×10⁻⁵     |
| 0.008  | π × 0.000028             | 8.796×10⁻⁵     |
| 0.015  | π × 0.000056             | 1.759×10⁻⁴     |
| 0.020  | π × 0.000076             | 2.388×10⁻⁴     |
| 0.030  | π × 0.000116             | 3.644×10⁻⁴     |

**Node phloem volumes and max sugar capacities:**

| Node | Formula | phloem_vol (dm³) | max_cap (g) |
|------|---------|-----------------|-------------|
| leaf1 | 0.5² × 0.003 | 7.500×10⁻⁴ | 0.2250 |
| stem3 | 1.759×10⁻⁴ × 0.05 | 8.797×10⁻⁶ | 0.00264 |
| stem2 | 1.759×10⁻⁴ × 0.05 | 8.797×10⁻⁶ | 0.00264 |
| stem1 | 2.388×10⁻⁴ × 0.05 | 1.194×10⁻⁵ | 0.00358 |
| seed | 3.644×10⁻⁴ × 0.02 | 7.289×10⁻⁶ | 0.00219 |
| root1 | 1.759×10⁻⁴ × 0.10 | 1.759×10⁻⁵ | 0.00528 |
| root_apical | (4/3)π×0.008³ | 2.145×10⁻⁶ | 0.000643 |
| leaf2 | 0.3² × 0.003 | 2.700×10⁻⁴ | 0.0810 |
| shoot_apical | (4/3)π×0.005³ | 5.236×10⁻⁷ | 0.000157 |

---

## Step 2: Sugar Concentrations

`sugar_conc = sugar / phloem_volume` (g/dm³)

| Node | sugar (g) | phloem_vol (dm³) | conc (g/dm³) | % of 300 cap | Status |
|------|-----------|-----------------|--------------|--------------|--------|
| leaf1 | 0.030 | 7.500×10⁻⁴ | 40.00 | 13.3% | ✓ |
| stem3 | 0.001 | 8.797×10⁻⁶ | 113.68 | 37.9% | ✓ |
| stem2 | 0.001 | 8.797×10⁻⁶ | 113.68 | 37.9% | ✓ |
| stem1 | 0.500 | 1.194×10⁻⁵ | **41,883** | **13,961%** | ⚠️ OVERFLOW |
| seed | 2.000 | 7.289×10⁻⁶ | **274,405** | **91,468%** | ⚠️ OVERFLOW |
| root1 | 0.001 | 1.759×10⁻⁵ | 56.84 | 18.9% | ✓ |
| root_apical | 0.005 | 2.145×10⁻⁶ | **2,331** | **777%** | ⚠️ OVERFLOW |
| leaf2 | 0.010 | 2.700×10⁻⁴ | 37.04 | 12.3% | ✓ |
| shoot_apical | 0.002 | 5.236×10⁻⁷ | **3,820** | **1,273%** | ⚠️ OVERFLOW |

**4 of 9 nodes overflow max_sugar_concentration.** The overflowing nodes have initial sugar that is 8–91,000× more than their phloem ring can physically hold at 300 g/dm³.

---

## Step 3: Water Fractions

`water_frac = clamp(water / water_cap, 0, 1)`

Water cap formulas:
- STEM/ROOT: `π × r² × L × 800`
- LEAF: `leaf_size² × 3.0`
- APICAL/ROOT_APICAL: `1.0 ml`

| Node | water (ml) | water_cap (ml) | water_frac | Status |
|------|-----------|----------------|-----------|--------|
| leaf1 | 3.00 | 0.750 | 1.000 | ⚠️ 4.0× over cap |
| stem3 | 0.50 | 0.0283 | 1.000 | ⚠️ 17.7× over cap |
| stem2 | 0.50 | 0.0283 | 1.000 | ⚠️ 17.7× over cap |
| stem1 | 1.00 | 0.0503 | 1.000 | ⚠️ 19.9× over cap |
| seed | 1.00 | 0.0452 | 1.000 | ⚠️ 22.1× over cap |
| root1 | 2.00 | 0.0566 | 1.000 | ⚠️ 35.4× over cap |
| root_apical | 0.50 | 1.0 | **0.500** | ✓ |
| leaf2 | 1.00 | 0.270 | 1.000 | ⚠️ 3.7× over cap |
| shoot_apical | 0.20 | 1.0 | **0.200** | ✓ |

**7 of 9 nodes have water exceeding water_cap** — all clamp to 1.0. Only the two meristems are within normal range. The test values are unrepresentative of actual node water capacity. Practically: all non-meristem nodes behave as "fully hydrated" (water_frac = 1.0). This is not critical — it simply means drought is not in play here.

---

## Step 4: Phloem Pressures

`phloem_pressure = concentration × κ × water_frac` (κ = 1.0)

| Node | conc (g/dm³) | × water_frac | pressure |
|------|-------------|-------------|---------|
| seed | 274,405 | × 1.000 | **274,405** |
| stem1 | 41,883 | × 1.000 | **41,883** |
| root_apical | 2,331 | × 0.500 | **1,166** |
| shoot_apical | 3,820 | × 0.200 | **764** |
| stem3 | 113.68 | × 1.000 | 113.68 |
| stem2 | 113.68 | × 1.000 | 113.68 |
| root1 | 56.84 | × 1.000 | 56.84 |
| leaf1 | 40.00 | × 1.000 | 40.00 |
| leaf2 | 37.04 | × 1.000 | 37.04 |

**Pressure landscape is biologically inverted on every edge involving leaves or meristems:**
- Seed: 274,000 (reservoir, should be a gradual source but dominates catastrophically)
- Root apical: 1,166 (should be a sink, appears as a massive source)
- Shoot apical: 764 (should be a sink, appears as a source over adjacent stem)
- Leaves: 37–40 (should be sources, are the lowest-pressure nodes in the graph)

---

## Step 5: Pipe Capacities and Phloem Speeds

`pipe_cap = ring_area(r_eff) × conductance` where `r_eff = min(r_parent, r_child)`

`speed = base_phloem_speed × (r_eff / r_ref)²` = `5.0 × (r_eff / 0.05)²`

`time_cost = edge_length / speed`

| Edge | r_eff | pipe_cap (dm³/tick) | speed (dm/tick) | time_cost (tick) |
|------|-------|--------------------|-----------------|--------------------|
| seed → stem1 | 0.020 | 1.910×10⁻³ | 0.800 | 0.063 |
| seed → root1 | 0.015 | 1.407×10⁻³ | 0.450 | 0.222 |
| stem1 → stem2 | 0.015 | 1.407×10⁻³ | 0.450 | 0.111 |
| stem2 → stem3 | 0.015 | 1.407×10⁻³ | 0.450 | 0.111 |
| stem2 → leaf2 | 0.015 | 1.407×10⁻³ | 0.450 | 0.044 |
| stem3 → leaf1 | 0.015 | 1.407×10⁻³ | 0.450 | 0.044 |
| stem3 → shoot_apical | 0.005 | 4.021×10⁻⁴ | 0.050 | 0.100 |
| root1 → root_apical | 0.008 | 7.037×10⁻⁴ | 0.128 | 0.078 |

**Distance budget:** Full path seed→stem1→stem2→stem3→leaf1 costs 0.063+0.111+0.111+0.044 = 0.329 ticks. All edges reachable within the 1.0-tick budget. Sugar can travel the entire shoot in a single tick.

---

## Step 6: Pressure Gradients and Flow Directions

Flow direction: from higher pressure to lower pressure. Expected biological directions are annotated.

| Edge | P_parent | P_child | Actual flow direction | Expected | Match? |
|------|---------|---------|----------------------|----------|--------|
| seed → stem1 | 274,405 | 41,883 | seed → stem1 | seed → stem1 | ✓ |
| seed → root1 | 274,405 | 56.84 | seed → root1 | seed → root1 | ✓ |
| stem1 → stem2 | 41,883 | 113.68 | stem1 → stem2 | stem should relay, not source | ~ |
| stem2 → stem3 | 113.68 | 113.68 | no flow (equal) | should relay | ~ |
| stem2 → leaf2 | 113.68 | 37.04 | **stem2 → leaf2** | leaf2 → stem2 | ✗ |
| stem3 → leaf1 | 113.68 | 40.00 | **stem3 → leaf1** | leaf1 → stem3 | ✗ |
| stem3 → shoot_apical | 113.68 | 763.94 | **shoot_apical → stem3** | stem3 → shoot_apical | ✗ |
| root1 → root_apical | 56.84 | 1,165.69 | **root_apical → root1** | root1 → root_apical | ✗ |

**5 of 8 edges flow in the wrong biological direction.** Only the seed-drain edges are directionally correct (because the seed is so catastrophically overloaded it overwhelms everything).

---

## Step 7: Flow Magnitude — seed → root1 (Full Worked Example)

This is the dominant BFS walk. Source: seed (P=274,405).

**Edge seed → root1:**
```
stream_conc = seed_sugar / seed_phloem_vol = 2.0 / 7.289e-6 = 274,405 g/dm³
pipe_cap    = ring_area(0.015) × 8.0 = 1.7593e-4 × 8.0 = 1.407e-3 dm³/tick
time_cost   = 0.1 dm / 0.45 dm/tick = 0.222 ticks (within 1-tick budget)
time_frac   = 1.0

flow_vol = 274,405 × 1.407e-3 × 1.0 = 386.2 g  ← 193× total plant sugar
```

**Unloading at root1:**
```
root1_phloem_vol = 1.759e-5 dm³
local_conc  = 0.001 / 1.759e-5 = 56.84 g/dm³
gradient    = 274,405 - 56.84 = 274,348 g/dm³
perm        = 0.008
unload_raw  = 274,348 × 0.008 × 386.2 = 847,641 g  →  clamped to 386.2 g

delta[seed]  -= 386.2 g   → seed_sugar = clamp(2.0 - 386.2, 0, 0.00219) = 0 g
delta[root1] += 386.2 g   → root1_sugar = clamp(0.001 + 386.2, 0, 0.00528) = 0.00528 g
```

**Result:** Seed drains to 0g. Root1 fills from 0.001g to 0.00528g (+0.00428g).
**Sugar destroyed by clamping: 386.2 − 0.00428 = 386.2 g.** Total sugar goes from ~2.55g to ~0.01g — 99.6% destroyed in one tick.

This is total conservation failure.

---

## Step 8: Root Cause Analysis

### Issue 1: Reservoir nodes (seed, stem1) overflow phloem ring capacity by 100–91,000×

The seed stores 2.0g of carbohydrate reserve. The seed phloem ring can physically hold at most `7.289×10⁻⁶ dm³ × 300 g/dm³ = 0.00219 g`. The excess (1.998g) is biologically stored as starch in parenchyma cells — but the model tracks it as "sugar" in the same field used for phloem concentration.

**Root cause:** The model conflates total carbohydrate reserves with phloem-loaded sugar. In reality, phloem loading is a rate-limited process — companion cells actively pump sucrose from parenchyma into sieve tubes at ~0.02 g/(dm²·hr). The model has no such rate limit.

### Issue 2: Wrong flow direction — stems always higher concentration than connected leaves

With phloem ring volumes, leaves use `leaf_size² × leaf_thickness` (total mesophyll volume) while stems use ring area × length. This creates a systematic volume asymmetry:

```
leaf1 phloem_vol  = 0.00075 dm³   (85× larger than adjacent stem)
stem3 phloem_vol  = 8.80×10⁻⁶ dm³

Ratio: 85×

For leaf to be higher concentration than stem:
  sugar_leaf / vol_leaf > sugar_stem / vol_stem
  sugar_leaf / sugar_stem > vol_leaf / vol_stem = 85

Test plant: 0.030 / 0.001 = 30 < 85  →  stem wins, flow inverted.
```

A leaf with 0.030g and an adjacent stem with 0.001g gives: leaf concentration 40 g/dm³, stem concentration 114 g/dm³. Sugar flows FROM the stem INTO the leaf. This is always the case unless the leaf has 85× more sugar than the stem.

### Issue 3: Meristems appear as pressure sources, not sinks

Shoot apical (r=0.005): phloem_vol = `(4/3)π×0.005³ = 5.24×10⁻⁷ dm³`. Even 2mg of residual sugar gives concentration 3,820 g/dm³. With water_frac = 0.2, pressure = 764 — far above adjacent stem (113.68). Sugar flows OUT of the meristem into the stem.

The tiny sphere volume makes any nonzero residual sugar produce enormous concentration. Meristems would only become sinks if residual sugar after growth consumption is < ~15 µg (for the test parameters). This requires essentially complete depletion during DFS.

### Issue 4: phloem_conductance = 8.0 is far too high for the ring model

With a young stem ring area of 1.759×10⁻⁴ dm², the pipe_cap = 1.759×10⁻⁴ × 8.0 = 1.41×10⁻³ dm³/tick. At leaf concentration 40 g/dm³:

```
flow_vol = 40 × 1.41e-3 × 1.0 = 0.0563 g/tick from one edge
```

leaf1 only has 0.030g total. One tick drains 188% of the leaf's entire sugar supply through a single edge. The entire leaf would drain to zero in less than one tick.

The conductance was calibrated for the old full-circle model where `π×r² = π×0.015² = 7.07×10⁻⁴ dm²` — 4× larger than the ring area. But the concentration is now much higher with ring volumes, so the net effect is ~8× more flow than intended.

**Target:** leaf should export ~10–20% of surplus sugar per tick (0.003–0.006g).
```
Required pipe_cap = 0.004 / 40.0 = 1.0×10⁻⁴ dm³/tick
Required conductance = 1.0×10⁻⁴ / 1.759×10⁻⁴ = 0.57 dm/tick
```

`phloem_conductance` should be reduced from 8.0 to approximately **0.5–1.0** for stable operation with the ring model.

---

## Step 9: Corrected Scenario — Post-DFS State

The test plant values represent a mix of pre-DFS and steady-state values. The algorithm runs AFTER the DFS tick. Recalculating with realistic post-DFS values:

**Corrected sugar amounts** (stems nearly depleted after growth; meristems near zero; seed capped at phloem capacity; leaves retain photosynthesis surplus):

```
leaf1:         0.030 g  (photosynthesis surplus retained)
stem3:         0.00005 g  (grew this tick, consumed almost all)
stem2:         0.00005 g
stem1:         0.00005 g
seed:          0.00200 g  (capped at phloem capacity: 7.3e-6 × 300 = 0.00219g)
root1:         0.00100 g  (maintenance only, not growing)
root_apical:   0.00003 g  (grew this tick, near zero)
leaf2:         0.010 g
shoot_apical:  0.00003 g  (grew this tick, near zero)

Total: 0.043 g
```

**Corrected concentrations:**

| Node | conc (g/dm³) | % of cap |
|------|-------------|---------|
| seed | 274.4 | 91.5% |
| root1 | 56.84 | 18.9% |
| shoot_apical | 57.3 | 19.1% |
| leaf1 | 40.00 | 13.3% |
| leaf2 | 37.04 | 12.3% |
| root_apical | 14.0 | 4.7% |
| stem3 | 5.68 | 1.9% |
| stem2 | 5.68 | 1.9% |
| stem1 | 4.19 | 1.4% |

**Corrected pressures** (water_frac unchanged from Step 3):

| Node | pressure |
|------|---------|
| seed | 274.4 (seed near phloem saturation — strong source) |
| root1 | 56.84 |
| shoot_apical | 11.46 |
| leaf1 | 40.00 |
| leaf2 | 37.04 |
| root_apical | 6.99 |
| stem3 | 5.68 |
| stem2 | 5.68 |
| stem1 | 4.19 |

**Corrected flow directions:**

| Edge | P_parent | P_child | Flow direction | Bio-correct? |
|------|---------|---------|---------------|-------------|
| seed → stem1 | 274.4 | 4.19 | seed → stem1 | ✓ |
| seed → root1 | 274.4 | 56.84 | seed → root1 | ✓ |
| stem1 → stem2 | 4.19 | 5.68 | **stem2 → stem1** | (depends on topology) |
| stem2 → stem3 | 5.68 | 5.68 | no flow (equal) | — |
| stem2 → leaf2 | 5.68 | 37.04 | **leaf2 → stem2** | ✓ |
| stem3 → leaf1 | 5.68 | 40.00 | **leaf1 → stem3** | ✓ |
| stem3 → shoot_apical | 5.68 | 11.46 | **shoot_apical → stem3** | ✗ |
| root1 → root_apical | 56.84 | 6.99 | root1 → root_apical | ✓ |

**6 of 8 edges now correct** in the corrected scenario. The persistent problem is shoot_apical, which retains enough pressure (P=11.46) to appear as a source over depleted stems (P=5.68).

### Corrected leaf1 → stem3 flow (detailed)

Leaf1 is now correctly a source over stem3:
```
leaf1 conc = 40.00 g/dm3,  stem3 conc = 5.68 g/dm3
gradient   = 34.32 g/dm3
pipe_cap   = 1.759e-4 × 8.0 = 1.407e-3 dm3/tick
flow_vol   = 40.00 × 1.407e-3 × 1.0 = 0.0563 g

unload = min(34.32 × 0.10 × 0.0563, 0.0563) = min(0.193, 0.0563) = 0.0563 g  (full flow unloads)

leaf1 pays:  0.0563 g  →  new_leaf1 = clamp(0.030 - 0.0563, 0, 0.225) = 0 g
stem3 gains: 0.0563 g  →  new_stem3 = clamp(0.00005 + 0.0563, 0, 0.00264) = 0.00264 g

Headroom absorbed: 0.00264 - 0.00005 = 0.00259 g
Sugar destroyed:   0.0563 - 0.00259 = 0.054 g  (95% destroyed)
```

Even in the corrected scenario with proper flow direction, 95% of the moved sugar is destroyed by clamping. Leaf1 drains to zero in one tick but stem3 only absorbs 0.00259g — limited by the tiny phloem ring cap (0.00264g max).

**For stable operation:** phloem_conductance needs to be reduced to ~0.5. At conductance=0.5, flow_vol = 40 × 1.759×10⁻⁴ × 0.5 = 0.00352 g/tick — within stem3's headroom (0.00259g). Leaf1 would drain in ~8 ticks instead of 1 tick. Sugar conservation would hold.

---

## Step 10: Sanity Checks

| Check | Result |
|-------|--------|
| Sugar conserved per tick? | **NO** — 386g destroyed in worst case (seed scenario). 95% destroyed even in corrected scenario due to clamp mismatch. |
| Flow rates reasonable? | **NO** — leaf drains to zero in 1 tick; flow_vol >> node capacity at most edges. |
| Pressure gradient correct direction? | **NO** — 5/8 edges inverted (stems have higher pressure than adjacent leaves). |
| Meristems behave as sinks? | **NO** — both apicals appear as pressure sources due to tiny phloem sphere volume. |
| Would plant survive? | **NO** — seed drains to zero in 1 tick, leaf drains to zero in 1 tick, 99%+ of sugar is destroyed each tick. |
| Any conservation-valid flows? | Corrected scenario has 3 correct flows (leaf→stem2, leaf→stem3, root1→root_apical) but flow rates still too high. |

---

## Summary: Flags and Recommendations

### Flag 1 (CRITICAL): phloem_conductance = 8.0 is ~10× too high for the ring model

**Evidence:** Flow from one leaf to one stem = 0.0563g/tick; leaf has 0.030g total. System drains any node in < 1 tick.

**Fix:** Reduce `phloem_conductance` to **0.5–1.0 dm/tick**. Derivation:
```
target_flow ≈ 10% of leaf sugar per tick = 0.003 g/tick
pipe_cap_needed = 0.003 / 40.0 = 7.5e-5 dm³/tick
conductance = 7.5e-5 / 1.759e-4 = 0.43 dm/tick  →  use 0.5
```

### Flag 2 (CRITICAL): Sugar storage incompatible with phloem ring capacity

**Evidence:** Seed has 2.0g; seed phloem ring holds max 0.00219g. Ratio 914:1. Same for any stem with stored reserves.

**Root cause:** The model conflates total carbohydrate reserves (starch + soluble sugar in vacuoles + apoplast) with phloem-loaded sucrose. In real plants these are distinct compartments.

**Fix options (choose one):**
- **A. Starch reservoir** (from `storage/plan.md`): Store reserves as Starch, not Sugar. After the starch system is implemented, Sugar = active phloem fraction, capped at phloem_volume × max_conc.
- **B. Loading fraction**: Add a WorldParam `phloem_loading_fraction` (e.g., 0.005–0.01). Effective phloem pressure uses `sugar × loading_fraction / phloem_volume`. Starch is implicitly the reservoir.
- **C. Enforce cap at initialization**: Clamp Sugar to phloem_volume × max_sugar_concentration everywhere, and represent the excess differently. Simplest but loses reserve tracking.

**Recommendation:** Wait for starch system (storage/plan.md). Implement Münch first with the understanding that large reserves must be added as Starch, not Sugar.

### Flag 3 (SIGNIFICANT): Volume asymmetry inverts leaf-stem flow direction

**Evidence:** Leaves use total volume (0.00075 dm³) but stems use ring volume (8.8×10⁻⁶ dm³), a 85× ratio. Stem 0.001g / ring_vol → 114 g/dm³ vs. leaf 0.030g / total_vol → 40 g/dm³. Flow goes stem→leaf.

**Condition for leaf to export correctly:**
```
sugar_leaf / sugar_stem > phloem_vol_leaf / phloem_vol_stem = 85
```
At typical leaf/stem sugar ratios (~30:1 here), stems always win.

**Fix options:**
- **A. Use total node volume for stems/roots too** (revert to earlier design): Concentration for stems = sugar / (π×r²×L). Symmetric treatment — no asymmetry. With 0.001g stem sugar and total_vol = 3.53×10⁻⁵ dm³: concentration = 28.3 g/dm³ < leaf's 40 g/dm³. Correct direction.
- **B. Keep ring volume but use higher leaf phloem concentration**: Leaves must have much higher sugar to overcome the asymmetry. This means leaves need to accumulate ~10× more sugar before loading begins — which creates a "loading threshold" that is biologically plausible (leaves do accumulate sugar before exporting) but must be tuned carefully.

**Recommendation:** Option A (use full node volume for concentration) is simpler, preserves correct flow direction, and still uses phloem ring area for pipe_capacity (which is the correct biology for flow throughput). The concentration denominator should reflect the tissue compartment from which sugar is drawn (entire cell content), while the pipe area reflects the actual sieve tube cross-section.

### Flag 4 (MODERATE): Meristems appear as pressure sources

**Evidence:** Shoot_apical with 30µg residual sugar has phloem concentration 57 g/dm³ (water_frac = 0.2, pressure = 11.46), higher than depleted stems at 5.68. Sugar flows meristem → stem.

**Threshold for meristem to be a sink** (relative to a stem at 5.68 pressure):
```
shoot_apical: needs sugar < 14.9 µg (essentially complete depletion during DFS)
root_apical:  needs sugar < 59 µg (larger sphere, more lenient)
```

**Fix:** Combined effect of Flag 3 fix (higher stem volume → lower stem concentration → easier for meristems to be above stems) and realistic post-DFS depletion. If stems have 5.68 pressure, meristems need to be below that. With option A above, depleted stems have ~28 g/dm³ (not 5.68), and meristems need < 7.4µg — even more stringent.

The real fix: DFS growth consumption must nearly completely deplete meristem sugar each tick. Any residual > ~10µg turns the meristem into a pressure source. This is physically realistic (meristems consume aggressively) but numerically fragile.

**Alternative:** Clip meristem sugar at phloem_volume × max_conc immediately after DFS, before Münch runs. This prevents any residual above the saturation point from distorting pressure.

### Flag 5 (MINOR): Water test values exceed water_cap by 4–35×

The test plant water values are realistic for the whole-plant scale but not for individual node volumes. This doesn't affect the core analysis (all non-meristem nodes clamp to water_frac = 1.0), but test values should be regenerated using actual water_cap formulas when running implementation tests.

---

## What Would Make This Work

**Immediate changes to parameters:**

1. `phloem_conductance`: 8.0 → **0.5** (10× reduction for ring model)

2. Stem/root concentration formula: use **full node volume** (π×r²×L) not ring volume for the concentration denominator in pressure. Keep ring area for pipe_capacity.

3. `max_sugar_concentration` would then apply to the full node volume. For a stem (r=0.015, L=0.05): total_vol = 3.53×10⁻⁵ dm³, max_sugar = 3.53×10⁻⁵ × 300 = 0.0106g. Seed: π×0.03²×0.02 = 5.65×10⁻⁵ dm³, max = 0.0170g. Still well below 2.0g — Flag 2 must be resolved by starch system.

4. Enforce `node.sugar = clamp(node.sugar, 0, max_cap)` at the start of phloem_resolve (before BFS), using the corrected volume formula. This ensures initial values don't produce unbounded pressures.

**Corrected flow with recommended parameters (leaf1 → stem3 with conductance = 0.5, full node volume):**
```
stem3 full_vol  = π×0.015²×0.05 = 3.534×10⁻⁵ dm³
stem3 conc      = 0.00005 / 3.534×10⁻⁵ = 1.41 g/dm³
leaf1 conc      = 0.030 / 7.5×10⁻⁴ = 40.0 g/dm³
gradient        = 38.59 g/dm³
pipe_cap        = ring_area(0.015) × 0.5 = 1.759×10⁻⁴ × 0.5 = 8.795×10⁻⁵ dm³/tick
flow_vol        = 40.0 × 8.795×10⁻⁵ = 0.00352 g/tick
unload          = 38.59 × 0.10 × 0.00352 = 0.000136 g/tick → clamp check: 0.00136 < 0.00352 ✓
stem3 max_cap   = 3.534×10⁻⁵ × 300 = 0.0106 g  (well above current 0.00005g)
stem3 gains     = 0.000136 g  (absorbed fully, no clamping)
leaf1 loses     = 0.00352 g  →  new_leaf1 = 0.030 - 0.00352 = 0.0265 g  (12% exported per tick)
```

This is sensible: 12% export per tick, no clamping, conservation holds. Leaf1 would stabilize over ~8 ticks of transport.

---

## Conclusion

The current parameter set (conductance=8.0, phloem ring volume for concentration) **does not produce a working simulation** on this test plant. Three independent parameter/design issues combine to produce catastrophic flow failure:

1. phloem_conductance = 8.0 produces flow rates 10× too high for ring pipe areas
2. Phloem ring volume for concentration denominator inverts leaf-stem flow direction (stems appear as sources over leaves at realistic sugar amounts)
3. Sugar in reservoir nodes (seed, stem1) and meristems exceeds phloem ring capacity by 2–91,000× — requires starch system to separate reserves from active phloem loading

**Priority of fixes before implementation:**
1. Decide on concentration denominator: phloem ring volume (Option 3) vs. full node volume. The trial calculation shows that phloem ring volume creates systematic flow inversion unless complemented by very high leaf sugar (>85× stem sugar). Full node volume gives correct leaf→stem direction at realistic sugar ratios and should be used for concentration.
2. Reduce phloem_conductance to ~0.5 (calibrated for new pipe areas and concentration scale).
3. Implement starch system before or alongside Münch, so seed/stem reserves don't corrupt phloem pressure.
4. Enforce sugar cap in phloem_resolve before BFS runs.
