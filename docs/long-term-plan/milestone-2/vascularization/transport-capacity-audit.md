# Vascular Transport Capacity Audit

**Date:** 2026-04-17  
**Source files:** `vascular.cpp`, `genome.h`, `world_params.h`, `node/leaf.cpp`, `node/root_node.cpp`, `node/root_apical.cpp`, `node/stem_node.cpp`, `node/node.cpp`

---

## Formulas and Parameters

### Pipe Capacity

```cpp
// vascular.cpp — pipe_capacity()
pipe_capacity = π × radius² × conductance
```

| Parameter | Value | Unit |
|-----------|-------|------|
| `phloem_conductance` | 8.0 | g glucose / (dm² · tick) |
| `xylem_conductance` | 10.0 | ml / (dm² · tick) |
| `initial_radius` | 0.015 dm | (1.5 mm) |
| `root_initial_radius` | 0.008 dm | (0.8 mm) |

1 tick = 1 hour.

### Sugar Production

```
sugar_produced = light_exposure × angle_efficiency × world.light_level
               × leaf_size² × sugar_production_rate × stomatal
```

`sugar_production_rate = 0.02 g/(dm²·hr)` at full sun.  
At full conditions (all factors = 1.0) with full-size leaf (leaf_size = 1.5 dm):
```
leaf_area = 1.5² = 2.25 dm²
sugar_per_leaf = 2.25 × 0.02 = 0.045 g/hr
```

### Maintenance Costs

```
leaf_maintenance   = sugar_maintenance_leaf × leaf_size²     = 0.002 × 2.25 = 0.0045 g/hr per leaf
stem_maintenance   = sugar_maintenance_stem × π × r² × L     = 0.010 × π × 0.015² × 1.0 = 0.0000707 g/hr per stem (r=0.015, L=1dm)
root_maintenance   = sugar_maintenance_root × π × r² × L     = 0.004 × π × 0.008² × 1.0 = 0.0000080 g/hr per root segment
apical_maintenance = sugar_maintenance_meristem (if active)  = 0.0005 g/hr per active tip
```

### Water Absorption

```
root segment:  absorbed = water_absorption_rate × 2π × r × L × gradient  (gradient = soil_moisture − fill_fraction)
root apical:   absorbed = water_absorption_rate × 2π × r² × root_hair_multiplier × gradient
```

`water_absorption_rate = 0.5 ml/(dm²·hr)`, `root_hair_multiplier = 20` (hardcoded in `root_apical.cpp`).

At maximum gradient (empty root, soil_moisture = 1.0):
```
root segment (r=0.008, L=1 dm):  surface = 2π × 0.008 × 1.0 = 0.0503 dm²   → 0.5 × 0.0503 = 0.0251 ml/hr
root apical  (r=0.008):           surface = 2π × 0.008² × 20 = 0.00804 dm²  → 0.5 × 0.00804 = 0.0040 ml/hr
```

### Transpiration

```
transpired = transpiration_rate × leaf_size² × light_exposure
```

`transpiration_rate = 0.04 ml/(dm²·hr)`.  
Full-size leaf at full light: `2.25 × 0.04 = 0.09 ml/hr`.

---

## Pipe Capacity at Each Radius

| Radius | r / r₀ | Phloem (g/hr) | Xylem (ml/hr) |
|--------|--------|---------------|---------------|
| 0.015 dm (initial) | 1× | **0.00565** | **0.00707** |
| 0.021 dm | 1.4× | 0.01112 | 0.01385 |
| 0.030 dm | 2× | **0.02262** | **0.02827** |
| 0.050 dm | 3.3× | 0.06283 | 0.07854 |
| 0.075 dm | 5× | **0.14137** | **0.17671** |

---

## Sugar Analysis

### Scenario: Young Plant — 5 stem nodes + 3 leaves + 2 roots + 1 active apical

**Leaf sizes and production:**

| Leaf state | leaf_size | leaf_area | Production (g/hr) | Maintenance (g/hr) | Net export (g/hr) |
|------------|-----------|-----------|-------------------|-------------------|------------------|
| Bud | 0.02 dm | 0.0004 dm² | 0.000008 | 0.0000008 | ~0 |
| Quarter size | 0.375 dm | 0.141 dm² | 0.00281 | 0.000281 | 0.00253 |
| Half size | 0.75 dm | 0.5625 dm² | 0.01125 | 0.00113 | 0.01012 |
| Full size | 1.5 dm | 2.25 dm² | 0.04500 | 0.00450 | 0.04050 |

**Total sugar from 3 full leaves: 3 × 0.04050 = 0.1215 g/hr net export**

**Total maintenance for entire plant:**
```
3 leaves:   3 × 0.0045  = 0.0135 g/hr
5 stems:    5 × 0.0000707 = 0.000354 g/hr
2 roots:    2 × 0.0000080 = 0.000016 g/hr
1 apical:               0.0005 g/hr
Total:                  0.0143 g/hr
```

**Non-leaf maintenance demand (what must reach roots/meristems via phloem):**
```
stems + roots + apical = 0.000354 + 0.000016 + 0.0005 = 0.000870 g/hr
```

**Phloem pipe capacity at initial_radius: 0.00565 g/hr**

| Demand | Production | Pipe capacity | Pipe covers demand? | Pipe covers production? |
|--------|------------|---------------|---------------------|------------------------|
| 0.00087 g/hr (non-leaf nodes) | 0.1215 g/hr (3 full leaves) | 0.00565 g/hr | **Yes — 6.5×** | **No — 4.6%** |

**Finding:** Non-leaf nodes are adequately supplied — the pipe carries 6.5× what stems, roots, and meristems actually need. The bottleneck is not starvation of distant nodes. The bottleneck is **leaf export efficiency**: the pipe can drain only 4.6% of what 3 full leaves produce. Leaves fill to `sugar_cap = 2.0 × leaf_area = 4.5 g` and stop photosynthesizing.

### Effective Photosynthesis Rate Under Pipe Constraint

Once a leaf reaches `sugar_cap`, production halts. In steady state, production = export through pipe + own maintenance:

```
effective_production = pipe_share + own_maintenance
                     = (0.00565 / 3) + 0.0045       (one-third of pipe per leaf)
                     = 0.00188 + 0.00450
                     = 0.00638 g/hr per leaf at cap
```

vs. maximum: 0.045 g/hr per leaf. **Effective utilization: 14%.**

The plant wastes 86% of photosynthetic capacity because the phloem pipe cannot drain leaves fast enough. Plants running the fitness evaluation (200 ticks) hit this limit well before leaf development completes.

### At What Plant Size Does This Become a Problem?

Leaves expand over ~300 hours (0.005 dm/hr × 1.48 dm range). Meanwhile the trunk thickens
at max `cambium_responsiveness × auxin_flow_bias = 0.00002 dm/hr` (at full saturation):

```
r after 300 hrs = 0.015 + 300 × 0.00002 = 0.015 + 0.006 = 0.021 dm
pipe_capacity_phloem at 0.021 = π × 0.021² × 8.0 = 0.01108 g/hr
```

When leaves are full size (production 0.1215 g/hr), the trunk is at 0.021 dm and provides 0.011 g/hr — **9% coverage**. The bottleneck is present throughout leaf development and never resolves within a normal sim run.

**Break-even radius** (pipe capacity = 3-leaf net production of 0.1215 g/hr):
```
0.1215 = π × r² × 8.0
r = √(0.1215 / 25.13) = √0.00483 = 0.0695 dm
```

Reaching 0.0695 dm from initial 0.015 dm at max thickening rate: `(0.0695 − 0.015) / 0.00002 = 2725 hours` (113 days at maximum continuous saturated thickening). Not achievable in typical sim runs.

### Scenario: Mature Plant — 20 full leaves

```
Production:     20 × 0.045    = 0.90 g/hr
Maintenance:    20 × 0.0045   = 0.09 g/hr (leaves)
                + stem/root/apical ≈ 0.005 g/hr
Net to export:  ~0.80 g/hr
```

Pipe capacities:
```
At 1× (0.015 dm): 0.00565 g/hr  →  0.7% of net
At 5× (0.075 dm): 0.14137 g/hr  →  17.7% of net
At 10× (0.150 dm): 0.56549 g/hr →  70.7% of net
```

The main trunk needs to be **10× initial_radius** before phloem throughput approaches the production of a mature 20-leaf canopy. At cambium_responsiveness = 0.00002, that is over 250 days of continuous thickening. Clearly not calibrated to the sim's timescales.

---

## Water Analysis

### Scenario: Young Plant — 5 stems + 3 leaves + 2 root segments + 1 root apical

**Root absorption at maximum gradient (empty roots, soil_moisture = 1.0):**
```
2 root segments (r=0.008, L=1dm):  2 × 0.0251 = 0.0502 ml/hr
1 root apical (r=0.008):                        0.0040 ml/hr
Total absorption (max):                         0.0542 ml/hr
```

**Leaf water demand:**

| Leaf state | leaf_size | Transpiration (per leaf) | Photo water cost (per leaf) | Total per leaf |
|------------|-----------|--------------------------|----------------------------|---------------|
| Quarter size | 0.375 dm | 0.04 × 0.141 = 0.0056 ml/hr | 0.00281 × 0.5 = 0.0014 | 0.0070 |
| Half size | 0.75 dm | 0.04 × 0.5625 = 0.0225 | 0.01125 × 0.5 = 0.0056 | 0.0281 |
| Full size | 1.5 dm | 0.04 × 2.25 = 0.0900 | 0.045 × 0.5 = 0.0225 | 0.1125 |

**3-leaf demand at full size: 3 × 0.1125 = 0.3375 ml/hr**

**Balance:**

| Leaf state | 3-leaf demand (ml/hr) | Root absorption (max) | Status |
|------------|-----------------------|-----------------------|--------|
| Quarter size | 0.0211 | 0.0542 | ✓ Surplus (2.6×) |
| Half size | 0.0843 | 0.0542 | **Deficit 1.6×** |
| Full size | 0.3375 | 0.0542 | **Deficit 6.2×** |

**Root absorption crossover** (2 root segments + 1 apical becomes limiting):
```
0.0542 = 0.12 × leaf_size²   →   leaf_size = √(0.0542 / 0.12) = 0.67 dm
```

At leaf_size > **0.67 dm** (~45% of full size), the 2-root-segment system can't supply water
fast enough. The plant enters water stress before leaves are half-grown.

### Does the Xylem Pipe Matter for Water?

**Xylem pipe at initial_radius: 0.00707 ml/hr**

This is smaller than root absorption (0.0542 ml/hr) and far smaller than transpiration demand.
However, the xylem vascular pass is NOT the primary water mover. Local diffusion handles
bulk water movement:

```
Local diffusion (water_diffusion_rate = 0.9, water_base_transport = 0.2, water_transport_scale = 4.0):
  throughput cap = water_base_transport + radius_factor × water_transport_scale
                 = 0.2 + 1.0 × 4.0 = 4.2 ml/tick at initial_radius

  At maximum concentration gradient, actual desired transfer per hop:
    diff = (1.0 − 0.0) + water_bias = 1.05 (concentration units)
    avg_cap = water_storage_density_stem × volume ≈ 800 × 0.000707 = 0.566 ml
    desired = 1.05 × 0.9 × 0.566 × 1.0 = 0.534 ml/tick
```

Local diffusion can move **0.534 ml/tick** through a single stem hop — **75× more than the xylem pipe**. The xylem pipe capacity is negligible next to local diffusion. The `xylem_conductance = 10.0` parameter has minimal practical effect on water movement.

**The real water bottleneck is root absorption rate vs. transpiration demand**, not pipe capacity.

### Transpiration Has No Stomatal Regulation

```cpp
// leaf.cpp — transpire()
float transpired = g.transpiration_rate * leaf_area * light_exposure;
```

Stomata do NOT gate transpiration. Stomatal conductance only affects photosynthesis:
```cpp
float stomatal = clamp(water / water_cap, 0.2, 1.0);
sugar_produced = ... × stomatal;   // photosynthesis reduced
transpired     = ... ;              // unchanged — no stomatal factor
```

**Effect:** Even when the plant is water-stressed, it continues transpiring at full rate. Only photosynthesis slows. This means the plant continues losing water at maximum rate during stress, accelerating desiccation. Real plants close stomata under water stress, which simultaneously reduces photosynthesis AND transpiration. The current model breaks this coupling.

---

## Summary Table

| Metric | Value | Notes |
|--------|-------|-------|
| Phloem pipe at initial_radius | 0.00565 g/hr | Physical maximum per trunk segment |
| Non-leaf maintenance demand | 0.00087 g/hr | Roots + stems + 1 apical adequately supplied |
| 3-leaf net sugar export | 0.1215 g/hr | **21.5× beyond pipe capacity** |
| Effective photosynthesis at cap | 14% | Leaves fill and stop producing |
| Xylem pipe at initial_radius | 0.00707 ml/hr | Negligible vs local diffusion (0.534 ml/tick) |
| Root absorption (2 segments + apical) | 0.0542 ml/hr | Crossover at leaf_size = 0.67 dm |
| 3-leaf transpiration demand | 0.3375 ml/hr | **6.2× beyond root absorption** |
| Water stress onset | leaf_size ~0.67 dm | Before half-expansion |

---

## Findings and Recommendations

### Finding 1: Phloem Is a Major Bottleneck — But Doesn't Cause Starvation

**The pipe cannot drain leaves.** At initial_radius, phloem pipe = 0.00565 g/hr; 3 full leaves
produce 0.1215 g/hr. Leaves fill to sugar_cap and photosynthesis stalls at 14% of potential.

**However: non-leaf nodes don't starve.** Their total demand (0.00087 g/hr) is well below
the pipe capacity (0.00565 g/hr). Roots and meristems receive adequate sugar. The bottleneck
manifests as wasted photosynthetic potential, not starvation of distant tissue.

**This is probably observable in the sim:** leaves run near their sugar_cap much of the time,
and plant growth is sugar-limited even in full light — not because the plant is starving, but
because leaves can't export fast enough.

**Options if this is a problem:**
- Increase `phloem_conductance` from 8.0 → 80–200. This would allow thin stems to drain leaves more effectively without requiring thick trunks.
- Or: accept this as a design feature — in the biological analogy, this IS why trunk thickening matters (wider phloem supports larger canopies). The feedback loop (leaves fill → auxin production drops → canalization slows) may be working as intended.

**Recommendation:** Check empirically. Run `botany_sugar_test` and look at whether leaves are
spending significant ticks at `sugar_cap`. If photosynthesis is being capped frequently, either
the conductance or the production_rate needs adjustment.

### Finding 2: Xylem Conductance Is Effectively Irrelevant

**Local diffusion dominates water transport.** The xylem pipe at initial_radius moves 0.00707
ml/hr; local diffusion per hop moves up to 0.534 ml/tick. `xylem_conductance = 10.0` has
negligible impact on actual water distribution.

This is not necessarily wrong — xylem may be intended as a long-distance top-up on top of
diffusion, or as a mechanism that becomes meaningful at large scale. But any tuning of
`xylem_conductance` has essentially no effect on current plant behavior.

**No immediate action required.** Note for future: if the design intent is for xylem to be
the PRIMARY water mover (not just supplementary), local diffusion `water_base_transport`
and `water_transport_scale` need to be reduced significantly first.

### Finding 3: Root Absorption Becomes Insufficient Around Half-Leaf-Size

**For a plant with 2 root segments**, water deficit begins when leaves reach 0.67 dm. A
balanced plant needs approximately **1 root segment per leaf** to maintain water equilibrium
at half-leaf-size, or more aggressive root growth to keep up with canopy expansion.

This is a parameter balance question, not a pipe capacity question. `water_absorption_rate`
and `transpiration_rate` are the relevant levers. Current values (0.5 vs 0.04 per dm²) give
a ratio of 12.5:1 — a root segment has 12.5× the water throughput per dm² that a leaf
transpires per dm². But root surface area scales linearly with length while leaf area scales
quadratically, so fast-expanding leaves quickly outpace linear root growth.

**No immediate action required** unless testing reveals the plant is chronically water-stressed
before achieving useful leaf coverage.

### Finding 4: Transpiration Has No Stomatal Gate

**Transpiration ignores stomatal conductance.** When the plant is water-stressed, it continues
transpiring at full rate while photosynthesis slows — the opposite of biological reality. In
real plants, stomata close under water stress, reducing both transpiration and photosynthesis.

This is a known simplification. It means the model's water depletion under stress is faster
than it should be (no emergency conservation response). It does not affect steady-state
(adequately-watered) behavior, only stress trajectories.

**Optional fix** (not blocking): multiply transpiration by stomatal conductance:
```cpp
float stomatal = clamp(chemical(Water) / wcap, 0.2f, 1.0f);
float transpired = g.transpiration_rate * leaf_area * light_exposure * stomatal;
```

---

## What the Numbers Say About Current Configuration

The current default configuration is calibrated for **demonstration correctness, not
biophysical efficiency.** The plant does not starve (non-leaf demand is tiny), but it wastes
most of its photosynthetic potential due to the phloem bottleneck. The xylem conductance
parameter has no measurable effect. Water balance works for seedlings but becomes strained
as leaves develop.

The PIN transport work (adding long-range auxin signal) will drive canalization and thickening,
which directly addresses the phloem bottleneck — a thicker trunk has more pipe capacity. The
current bottleneck may be intentional: it creates selection pressure for trunk thickening,
which is exactly what the PIN + canalization system is designed to produce.

**Before PIN is implemented:** the phloem bottleneck exists but plants cope by running leaves
near their sugar_cap. Growth is sugar-limited but non-zero.

**After PIN is implemented:** auxin reaches roots → more cytokinin → more shoot growth →
more leaves → more auxin → more canalization → more thickening → larger phloem pipe → less
bottleneck. The system should be self-correcting given enough sim time.
