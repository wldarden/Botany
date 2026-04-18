# Münch Pressure Flow — Trial Calculation

**Purpose:** Verify that the parameters and formulas in `munch-pressure-flow-design.md` produce
physically coherent flows before implementation begins. Identifies any design issues early.

**Script:** run `python3 docs/long-term-plan/milestone-2/storage/munch_trial.py` (embedded at the
bottom of this document) to reproduce.

---

## Test Plant Topology

```
leaf1 ──┐
         stem3 ──┐
shoot_apical ──┘  stem2 ──┐
         leaf2 ──┘         stem1 ── seed ── root1 ── root_apical
```

---

## Parameters Used

| Parameter | Value | Source |
|---|---|---|
| `phloem_ring_thickness` t | 0.002 dm | design doc §6 |
| `phloem_osmotic_coefficient` | 1.0 | design doc §6 |
| `max_sugar_concentration` | 300 g/dm³ | design doc §6 |
| `leaf_thickness` | 0.003 dm | design doc §6 |
| `phloem_conductance` | 10.0 | trial script |
| `base_phloem_speed` | 5.0 dm/tick | design doc §6 |
| `phloem_reference_radius` | 0.05 dm | design doc §6 |
| Permeabilities | leaf=0.10, stem=0.002, root=0.008, sa/ra=0.08 | design doc §4.3 |

---

## Node State (pre-flow)

| Name | r (dm) | offset | phloem_vol (dm³) | sugar (g) | conc (g/dm³) | w_frac | pressure | Münch cap (g) | Over cap? |
|---|---|---|---|---|---|---|---|---|---|
| leaf1 | 0.000 | 0.00 | 7.50×10⁻⁴ | 0.0300 | 40.0 | 1.000 | **40.0** | 0.225 | no |
| stem3 | 0.015 | 0.05 | 8.80×10⁻⁶ | 0.0010 | 113.7 | 1.000 | **113.7** | 0.00264 | no |
| stem2 | 0.015 | 0.05 | 8.80×10⁻⁶ | 0.0010 | 113.7 | 1.000 | **113.7** | 0.00264 | no |
| stem1 | 0.020 | 0.05 | 1.19×10⁻⁵ | 0.500 | 41,883 | 1.000 | **41,883** | 0.00358 | **140×** |
| seed | 0.030 | 0.02 | 7.29×10⁻⁶ | 2.000 | 274,405 | 1.000 | **274,405** | 0.00219 | **915×** |
| root1 | 0.015 | 0.10 | 1.76×10⁻⁵ | 0.0010 | 56.8 | 1.000 | **56.8** | 0.00528 | no |
| root_apical | 0.008 | 0.01 | 2.14×10⁻⁶ | 0.0050 | 2,331 | 0.500 | **1,166** | 0.00064 | 8× |
| leaf2 | 0.000 | 0.00 | 2.70×10⁻⁴ | 0.0100 | 37.0 | 1.000 | **37.0** | 0.081 | no |
| shoot_apical | 0.005 | 0.005 | 5.24×10⁻⁷ | 0.0020 | 3,820 | 0.200 | **764** | 0.00016 | 13× |

---

## Results: What Broke

### 1. Seed and stem1 have supraphysical concentrations

The seed has 2 g of sugar in a phloem ring of 7.29×10⁻⁶ dm³, giving
274,405 g/dm³ — nearly **1000× the 300 g/dm³ cap**. stem1 (0.5 g in a thin
conduit) is 140× over cap. These numbers come from the current plant model
calibrating sugar to whole-node volume; the Münch model computes concentration
against the phloem ring only (5–18× smaller volume for young stems, ~1% of
full cross-section for thick ones).

**Root cause:** The test-plant initial conditions were designed for the old
whole-cross-section sugar model. They are inconsistent with Münch phloem ring
volumes. Under the new model, `seed.sugar` should be at most
`phloem_ring_area(0.03, 0.002) × 0.02 × 300 ≈ 0.0022 g` in the phloem ring.
The seed's bulk starch reserve will be in a separate starch field (milestone-2
starch work), not in `node.sugar`.

### 2. BFS source capacity not enforced → conservation failure

The BFS computes `flow_vol = stream_conc × pipe_cap × time_fraction` from the
seed's pressure (274,405 × pipe_cap) and deducts it from `delta[seed]` without
checking whether the source actually has that much sugar. Result: the seed
"sends" −1,137 g in one tick from a 2 g reserve — physically impossible, and
conservation fails by −186 mg (non-zero sum of all deltas).

**Fix required in the implementation:** Before deducting from `delta[cur]`,
clamp `flow_vol` to at most `sugar_0[cur] + delta[cur]` (the running available
balance). This is a missing guard in the algorithm in §4.3.

### 3. Meristems and tiny nodes are sources, not sinks

The shoot_apical has 0.002 g in a 5.24×10⁻⁷ dm³ sphere → 3,820 g/dm³ → pressure 764.
The root_apical has 0.005 g in a 2.14×10⁻⁶ dm³ sphere → 1,166 pressure.
Both have higher pressure than their parent stems (113), so the BFS treats them
as *sources*. They push sugar into the stem instead of receiving it.

**Root cause:** The design spec says "fast-growing nodes have depleted their
sugar during the DFS tick → steep gradient → fast unloading." For the algorithm
to produce the intended sink behavior, meristems must arrive at the phloem
resolve pass with sugar **near zero** — specifically below
`phloem_volume × stream_concentration` of their parent. In a real post-DFS
tick, the meristem would have consumed almost all its sugar for growth. The
test setup using 0.002 g in a volume of 5×10⁻⁷ dm³ is still 13× over the
Münch cap.

**Implication for initial conditions and calibration:**
- Meristem sugar after DFS should be < ~0.0001 g for correct sink behavior.
- The `phloem_unloading_meristem` permeability (0.08) is appropriate once the
  meristem is actually empty; with high sugar it becomes a source instead.

### 4. Leaf radius = 0 produces negative ring area

`phloem_ring_area(0, t) = π × (0 − t²) = −π × t²`, yielding a negative pipe
capacity (−0.000126). Leaves have zero radius in the sim and will always
trigger this path. This is not a problem in the production code because
§4.5 specifies that **leaves are not in the vascular BFS** — they connect to
the nearest vascular stem via local diffusion. The BFS must skip leaf nodes
(and any other node with `has_vasculature == false`).

---

## Results: What Works

### Pressure directions are correct when initial conditions are valid

For nodes whose sugar is within the Münch cap range (leaf1, leaf2, stem3,
stem2, root1), the pressure landscape is coherent:

- leaf1 (40.0) > stem3 (113.7)? Actually stem3 > leaf1 here, which is because
  even 0.001 g in the tiny ring gives high concentration. Once stem conduits
  hold realistic Münch-scale sugar (< 0.003 g), their pressure drops below
  leaf level, and leaves correctly become sources.

- root1 (56.8) is a conduit — lower pressure than any accumulating node above.
  Flow from seed toward root1 direction is physically correct.

### Pipe capacity and speed are in the right ballpark

| Node | Speed (dm/tick) | Speed (m/hr) | Biological range |
|---|---|---|---|
| Young stem r=0.015 | 0.45 | 0.045 | **below** 0.3–1.5 m/hr |
| Reference r=0.05 | 5.00 | 0.50 | middle of range ✓ |
| Seed r=0.03 | 1.80 | 0.18 | slightly low |

Young stems run at ~9% of reference speed (0.045 m/hr vs. 0.3 m/hr minimum).
The design doc anticipates this: "Early seedling transport may require multiple
ticks." If this becomes a problem during calibration (Task 7), lower
`phloem_reference_radius` from 0.05 to 0.015 to make young-stem speed match
reference, or increase `base_phloem_speed`. This is a calibration concern, not
a design defect.

### Unloading permeabilities are well-scaled

With a 3 g/dm³ gradient, a mature stem leaks 0.006 g/tick and a meristem
(when empty) receives 0.24 g/tick — a 40× ratio. The design target is 20–100×
stem leakage. This falls squarely in range and should preserve most flow for
distant sinks rather than being absorbed by conduit walls.

### Leaf → conduit flow works correctly

leaf1 has conc = 40 g/dm³, pressure = 40. Its parent stem3 would have pressure
< 40 once conduit sugar is within Münch range (< 0.003 g → conc < 341). At
equilibrium sugar loads (stem3 ≈ 0.001 g or less), flow correctly runs
leaf→stem and onward.

---

## Summary of Findings

| # | Finding | Severity | Fix |
|---|---|---|---|
| 1 | Stem/seed sugar amounts massively exceed Münch phloem-ring cap | **Critical** | Rethink initial conditions: phloem sugar ≪ old whole-volume sugar; starch goes in separate starch field |
| 2 | BFS does not enforce source capacity → conservation failure | **Critical** | Clamp `flow_vol` to `available[src]` before deducting delta; add test |
| 3 | Meristems with any sugar become sources, not sinks | **Critical** | Post-DFS meristem sugar must be near zero; calibrate DFS tick consumption to fully empty them |
| 4 | Leaf nodes (radius=0) produce negative ring area in BFS | Moderate | Already handled by `has_vasculature` filter in production code; exclude non-vascular nodes from BFS |
| 5 | Young stem phloem velocity below biological range | Minor | Calibration knob: reduce `phloem_reference_radius` from 0.05→0.015 or increase `base_phloem_speed`; addressed in Task 7 |

### The core design is sound

The Münch model's conceptual architecture is correct: local concentrations drive
pressure, BFS propagates flow outward, unloading permeability filters the
stream, and distance-cost limits per-tick reach. The failures are all in
**initial conditions and one missing guard** (source capacity), not in the
algorithm itself.

The key insight surfaced by this calculation:

> **The Münch model requires node.sugar to represent only the phloem sieve-tube
> compartment.** The current model stores sugar as a whole-node resource calibrated
> to total tissue volume. These are different quantities. Implementation must either:
> (a) migrate sugar to represent only the phloem compartment (consistent with
>     `phloem_volume × max_conc` caps), OR
> (b) keep sugar as whole-node storage and derive phloem concentration as
>     `sugar × (phloem_vol / total_vol)` — a partitioning fraction.
>
> Option (a) is what the design doc describes and is the cleaner path. It means
> that at the start of implementation, all existing sugar values on nodes will
> be much smaller than current levels, and any "bulk storage" (seed reserves,
> starch) moves to separate fields introduced in the storage milestone.

---

## Corrected Trial: What the Numbers Look Like With Valid Initial Conditions

To validate the algorithm free of the initialization problem, here is a
manually corrected scenario where conduit sugar is within Münch caps:

| Node | sugar (g) | phloem_vol (dm³) | conc (g/dm³) | pressure |
|---|---|---|---|---|
| leaf1 | 0.030 | 7.5×10⁻⁴ | 40 | 40 (source) |
| stem3 | 0.0003 | 8.8×10⁻⁶ | 34 | 34 (conduit, just below leaf) |
| stem2 | 0.0002 | 8.8×10⁻⁶ | 23 | 23 (conduit) |
| stem1 | 0.001 | 1.2×10⁻⁵ | 84 | 84 (**needs reduction**) |
| seed | 0.002 | 7.3×10⁻⁶ | 274 | 274 (source — near-saturated) |
| root1 | 0.00010 | 1.8×10⁻⁵ | 5.7 | 5.7 (sink) |
| root_apical | 0.00001 | 2.1×10⁻⁶ | 4.7 | 4.7 (sink — nearly empty after DFS) |
| shoot_apical | 0.00001 | 5.2×10⁻⁷ | 19 | **19 (still a source!** → see below) |

Even in the corrected scenario, the shoot apical has a problem: it is a
**sphere of 0.005 dm radius** (vol = 5.24×10⁻⁷ dm³). If it has *any* sugar
above ~0.00001 g, its phloem concentration rapidly exceeds adjacent stems. The
meristem must be truly empty (< 5×10⁻⁵ g) after DFS to behave as a sink.

For root_apical (r=0.008, vol=2.14×10⁻⁶ dm³) the same applies, but there is
more volume so the threshold is higher (~0.0001 g).

**Conclusion:** The DFS tick must completely deplete active meristems. If the
meristem's maintenance + growth costs per tick exceed phloem delivery, it starves
and stops growing. If the DFS tick doesn't fully consume the meristem's sugar,
the Münch system will see it as a source and push sugar back toward the stem.
This means meristem growth and maintenance costs should be calibrated so that
the meristem genuinely empties each tick and relies on phloem_resolve to refill
it.

---

## Reproducible Script

```python
import math
from collections import defaultdict, deque

t            = 0.002   # phloem_ring_thickness (dm)
osmotic      = 1.0
max_conc     = 300.0
leaf_thick   = 0.003
conductance  = 10.0
base_speed   = 5.0
r_ref        = 0.05
perm = {'leaf': 0.10, 'stem': 0.002, 'root': 0.008, 'sa': 0.08, 'ra': 0.08}
W_DENS_STEM  = 800.0
W_DENS_LEAF  = 3.0
W_CAP_MERIS  = 1.0

nodes_raw = [
    ('leaf1',        'leaf', 0.0,   0.0,   0.5,  0.03,  3.0),
    ('stem3',        'stem', 0.015, 0.05,  0.0,  0.001, 0.5),
    ('stem2',        'stem', 0.015, 0.05,  0.0,  0.001, 0.5),
    ('stem1',        'stem', 0.02,  0.05,  0.0,  0.5,   1.0),
    ('seed',         'stem', 0.03,  0.02,  0.0,  2.0,   1.0),
    ('root1',        'root', 0.015, 0.1,   0.0,  0.001, 2.0),
    ('root_apical',  'ra',   0.008, 0.01,  0.0,  0.005, 0.5),
    ('leaf2',        'leaf', 0.0,   0.0,   0.3,  0.01,  1.0),
    ('shoot_apical', 'sa',   0.005, 0.005, 0.0,  0.002, 0.2),
]

edges = [
    ('stem3', 'leaf1'), ('stem3', 'shoot_apical'), ('stem2', 'stem3'),
    ('stem2', 'leaf2'), ('stem1', 'stem2'),
    ('seed',  'stem1'), ('seed',  'root1'), ('root1', 'root_apical'),
]
```

Run `python3 /path/to/munch_trial.py` to reproduce the full output with tables.
