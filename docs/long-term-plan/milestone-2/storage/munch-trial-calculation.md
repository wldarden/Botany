# Münch Pressure Flow — Trial Calculation v2

**Purpose:** Verify that the Münch parameters and three-step algorithm produce physically coherent
flows before any C++ is written. This document supersedes the v1 trial (which exposed critical
issues that were fixed here).

---

## 1. Ring Thickness Verdict

The original question: "what ring thickness lets a young stem carry ~0.035 g/tick without exceeding
300 g/dm³?"

### The key insight: throughput ≠ storage

The ring does **not** need to hold the full leaf surplus each tick. Sugar flows **through** the ring
as a pressurised stream; the ring holds only the in-transit fraction. The throughput formula is:

```
throughput = stream_conc × ring_area × conductance
```

For the young stem (r=0.015 dm, t=0.002 dm, conductance=10):

| t (mm) | ring_area (dm²) | stream_conc needed for 35 mg/tick | % of max_conc |
|--------|-----------------|-----------------------------------|---------------|
| 0.2    | 1.76×10⁻⁴       | **19.9 g/dm³**                    | **6.6%** ✓    |
| 0.5    | 3.93×10⁻⁴       | 8.9 g/dm³                         | 3.0% ✓        |
| 1.0    | 6.28×10⁻⁴       | 5.6 g/dm³                         | 1.9% ✓        |

**All thicknesses satisfy the ≤80% constraint.** The constraint is met at only 6.6% of max concentration
with the original t=0.002 dm (0.2 mm). No need to change ring thickness.

### Why t=0.002 dm is kept

`ring_volume` determines **pressure** (via concentration = sugar/ring_vol), not throughput capacity.
A thinner ring gives higher pressure per gram of sugar stored — matching the real phloem anatomy
where a thin active phloem layer (~0.1–0.5 mm, containing sieve tubes + companion cells + phloem
parenchyma) surrounds the dead heartwood. The conductance multiplier (10) scales throughput
independently. **t=0.002 dm (0.2 mm) retained.**

---

## 2. Revised Architecture: Endpoint Nodes

The v1 trial exposed a fundamental mismatch: leaves and meristems are not conduit pipe sections.
They need to be treated as **boundary conditions** on the phloem network, not internal nodes.

### Three architecture changes

**A. Seed uses total cylinder volume (not ring)**

The seed is a storage organ, not a thin-walled conduit. Using ring volume gives a
phloem_vol of ~7×10⁻⁶ dm³ for a 0.03 dm radius seed, which creates astronomical concentrations
for any biologically plausible sugar reserve (1.5 g → 200,000 g/dm³). Using total volume:

```
phloem_vol_seed = π × r² × L    (full cylinder, not just the outer ring)
```

With r=0.3 dm (enlarged to represent a storage organ), L=0.02 dm:
`phloem_vol = 5.65×10⁻³ dm³` → at 1.5 g sugar and full water: concentration = 265 g/dm³ ✓

**B. Leaves are endpoint sources — pre-BFS loading**

Leaves are not in the phloem BFS graph. They load sugar into their parent stem's phloem ring
**before** the BFS runs:

```
gradient  = max(0, leaf_conc − parent_ring_conc)    [g/dm³]
load_raw  = gradient × phloem_unloading_leaf × ring_area(parent) × phloem_conductance
load      = min(load_raw, room_in_parent_ring, leaf.sugar)
parent.sugar += load;  leaf.sugar -= load
```

The `room_in_parent_ring = phloem_vol(parent) × max_conc − parent.sugar` cap is essential:
a leaf at 160 g/dm³ trying to load into a stem at 15 g/dm³ would otherwise overflow the ring
(gradient=145, raw_load=0.026 g >> ring cap=0.003 g). The cap limits loading to one ring-full
per tick, which is physically correct — the stem pipe is already full before the BFS drains it.

**C. Meristems are terminal sinks — post-BFS unloading**

Meristems are not in the BFS graph. After BFS completes, each meristem unloads from its parent
stem's phloem. The **bottleneck radius** (min of parent and meristem) governs the pipe capacity,
and the unload is capped at the meristem's available room:

```
gradient   = max(0, parent_ring_conc − meristem_conc)
r_eff      = min(parent.radius, meristem.radius)
unload_raw = gradient × phloem_unloading_meristem × ring_area(r_eff) × phloem_conductance
unload     = min(unload_raw, parent.sugar, meristem_room)     ← meristem room cap required
```

Without the meristem room cap, a parent at 300 g/dm³ would transfer its entire sugar into the
meristem's tiny sphere (e.g. 0.003 g into a 0.00016 g cap = 1686% over-cap).

---

## 3. Algorithm Fixes Required

Two fixes to the BFS algorithm in the design doc (§4.3):

### Fix 1: Shared BFS delta for multi-source destination tracking

When multiple sources each run independent BFS walks (design doc §4.4), the second walk uses
the **pre-BFS sugar state** to compute destination room. If source A already filled stem1's ring,
source B doesn't know and tries to fill it again → double over-cap.

**Fix:** Pass a shared running `delta` array into each BFS walk. The destination room calculation
uses `current_sugar + shared_delta[dest]` so later BFS walks see what earlier ones already filled.

### Fix 2: Transit dead-end credit for conservation

The original formula `delta[source] -= flow_vol; delta[dest] += unload` is non-conserving when
BFS dead-ends (no downhill neighbours from dest). The transit fraction `flow_vol − unload` passes
through dest's conduit but has nowhere to go, disappearing from the sugar balance.

**Fix:** When BFS cannot continue from a node (no downhill neighbours or budget exhausted),
credit the transit back to that node's conduit:

```
if not will_continue:
    delta[dest] += (flow_vol − unload)   # transit stays in dest's sieve tube
```

This ensures `delta[source] + delta[dest] = 0` for each hop, giving exact conservation.

---

## 4. Test Plant and Parameters

```
leaf1 ──┐
         stem3 ──┐
shoot_apical ──┘  stem2 ──┐
         leaf2 ──┘         stem1 ── seed ── root1 ── root_apical
```

### Parameters

| Parameter | Value |
|---|---|
| `phloem_ring_thickness` | 0.002 dm (0.2 mm) — confirmed correct |
| `phloem_osmotic_coefficient` | 1.0 |
| `max_sugar_concentration` | 300 g/dm³ |
| `leaf_thickness` | 0.003 dm |
| `phloem_conductance` | 10.0 |
| `base_phloem_speed` | 5.0 dm/tick |
| `phloem_reference_radius` | 0.05 dm |
| Unloading permeabilities | leaf=0.10, stem=0.002, root=0.008, sa/ra=0.08 |

### Initial conditions (Münch-compatible)

Stem ring sugar is set to represent **in-transit concentration** (15–35 g/dm³), not whole-node
sugar as calibrated for the old model:

| Node | Type | sugar (g) | conc (g/dm³) | pressure | Notes |
|---|---|---|---|---|---|
| leaf1 | endpoint | 0.1200 | 160 | — | post-photosynthesis surplus |
| stem3 | conduit | 0.000132 | 15 | 15 | thin distal conduit |
| stem2 | conduit | 0.000220 | 25 | 25 | mid conduit |
| stem1 | conduit | 0.000418 | 35 | 35 | proximal conduit |
| seed | storage | 1.5000 | 265 | 265 | total vol, full water |
| root1 | conduit | 0.000088 | 5 | 5 | low-pressure sink |
| root_apical | endpoint | 0.000010 | ~0 | — | post-DFS, nearly empty |
| leaf2 | endpoint | 0.0600 | 222 | — | post-photosynthesis |
| shoot_apical | endpoint | 0.000010 | ~0 | — | post-DFS, nearly empty |

Pressure gradient (BFS network): **seed (265) > stem1 (35) > stem2 (25) > stem3 (15) > root1 (5)**

After leaf loading (Step 1): stem3 and stem2 both fill to 300 g/dm³ (at ring cap).

---

## 5. Results

### Final sugar state

| Name | Type | before (g) | after (g) | Δ (g) | conc after | % cap | Status |
|---|---|---|---|---|---|---|---|
| leaf1 | leaf | 0.120000 | 0.117493 | −0.002507 | 156.7 | 52% | ok |
| stem3 | stem | 0.000132 | 0.002492 | +0.002360 | 283 | 94% | ok |
| stem2 | stem | 0.000220 | 0.000000 | −0.000220 | 0 | 0% | ok |
| stem1 | stem | 0.000418 | 0.003581 | +0.003163 | 300 | 100% | ok |
| seed | seed_store | 1.500000 | 1.494286 | −0.005714 | 264 | 88% | ok |
| root1 | root | 0.000088 | 0.004644 | +0.004556 | 264 | 88% | ok |
| root_apical | ra | 0.000010 | 0.000643 | +0.000633 | 300 | 100% | ok |
| leaf2 | leaf | 0.060000 | 0.057581 | −0.002419 | 213 | 71% | ok |
| shoot_apical | sa | 0.000010 | 0.000157 | +0.000147 | 300 | 100% | ok |

**Conservation: Σ(Δ) = 0.00e+00 ✓  |  No negatives, no over-cap ✓**

---

## 6. Key Questions Answered

**Does sugar flow from leaf/seed toward hungry meristems?**
Yes. Both meristems received sugar and filled to 100% of their sphere cap:
- shoot_apical: +0.147 mg (filled to 0.157 mg cap)
- root_apical: +0.633 mg (filled to 0.643 mg cap)

The actual maintenance need at these volumes is nanograms/tick (maintenance = ~0.01 g/dm³/tick × 5×10⁻⁷ dm³ ≈ 5×10⁻⁹ g/tick). Inflow exceeds maintenance by ~30,000×. Meristems are abundantly supplied.

**How much sugar reaches each meristem per tick?**
- shoot_apical: 0.147 mg/tick (100% of cap per tick — refills completely each tick)
- root_apical: 0.633 mg/tick (100% of cap per tick — refills completely each tick)

**Is the seed's 1.5 g draining toward sinks?**
Yes. Seed Δ = −0.006 g/tick. Reserve lifetime: ~262 ticks (~11 days at 1 hr/tick). Without
the starch mobilisation model (future milestone), the seed provides steady supply for early
seedling establishment.

**Are any flows unreasonable?**
No. Conservation is exact. No node goes negative. No node exceeds its Münch cap. Flow
directions match the pressure gradient (seed → conduits → sinks).

**Would this plant survive?**
Yes. Both meristems fill completely each tick at 100% capacity. Leaves export 2–2.5 mg/tick
to parent conduits. Seed provides a trickle of ~5.7 mg/tick distributed across the whole
stem network. At these rates, the plant is stable.

**Stem2 drains to zero — is that a problem?**
No. Stem2 loaded from leaf2, then immediately drained all of it into stem1 via the BFS
(stem2 pressure=300 >> stem1 pressure=35). It acts as a pure pass-through conduit in this
tick. Next tick, leaf2 will refill it again. This is normal flowing-pipe behaviour.

---

## 7. Remaining Open Issues for Implementation

| # | Issue | Action |
|---|---|---|
| 1 | Seed total-volume vs ring: requires special-casing seed node type | Add `'seed_store'` branch to `phloem_volume()` |
| 2 | Leaf loading pre-BFS pass not in current design doc | Add to §3 tick order and §4 algorithm |
| 3 | Meristem post-BFS unloading not in current design doc | Add to §4 algorithm |
| 4 | Shared BFS delta required for multi-source correctness | Add note to §4.4 multi-source interaction |
| 5 | Transit dead-end credit required for conservation | Add to §4.3 apply_flow formula |
| 6 | Young stem speed (0.045 m/hr) below biological range | Task 7 calibration: lower `phloem_reference_radius` |

All issues are implementation details — the core Münch architecture (pressure-driven BFS, unloading
permeability table, ring-volume concentration model) is validated and working.

---

## 8. Reproducible Python Script

The calculation can be reproduced with:

```bash
python3 docs/long-term-plan/milestone-2/storage/munch_trial_v2.py
```

Key parameters are at the top of the file. The three steps are clearly labelled: leaf loading,
BFS, meristem unloading.
