# Sugar Economy Audit

Read-only analysis of production, consumption, vascular flow, and caps.
All math uses default genome + world_params values. "Tick" = 1 hour.

---

## 1. Production vs. Consumption Budget

### Production

**Leaf** (`leaf.cpp:70`, `world_params.h:25`)

```
production = light_exposure × angle_efficiency × light_level × leaf_size² × sugar_production_rate × stomatal
```

At full sun, full water, horizontal leaf:

| leaf_size (dm) | leaf area (dm²) | g/tick |
|---|---|---|
| 0.2 (seedling) | 0.04 | 0.0008 |
| 0.5 (growing) | 0.25 | 0.005 |
| 1.0 (mid) | 1.0 | 0.020 |
| 1.5 (max) | 2.25 | 0.045 |

`sugar_production_rate = 0.02 g/(dm²·hr)` calibration check:
real C3 peak ≈ 15 μmol CO2/m²/s = 0.0162 g glucose/(dm²·hr). **The value is correct.**

For a small plant with 5 leaves at avg size 0.5 dm → **~0.025 g/tick total**.

**Stem corticular** (`stem_node.cpp:23`, `genome.h:220`)

```
production = surface_area × light × stem_photosynthesis_rate
surface_area = 2π × r × L
```

Only stems thinner than `stem_green_radius_threshold = 0.04 dm` (all initial-radius stems
qualify: `initial_radius = 0.015 dm < 0.04 dm`). Rate = 0.005 g/(dm²·hr) — 25% of leaf rate
(chlorenchyma, correct biologically).

For 10 stems at r=0.015 dm, L=0.5 dm: 2π × 0.015 × 0.5 = 0.047 dm² → 0.00024 g/tick per stem.
10 stems total: **~0.0024 g/tick — 10% of leaf production, rounding error.**

**Apical** (`apical.cpp:47`): `sugar_meristem_photosynthesis = 0.0` → **zero**.

**Total production, typical small plant: ~0.027 g/tick**

---

### Consumption (Maintenance Only)

| Node | Cost formula | Typical value |
|---|---|---|
| Leaf (size 0.5 dm) | `0.002 × 0.5²` | 0.0005 g/tick |
| Stem (r=0.015, L=0.5) | `0.01 × π × 0.015² × 0.5` | 0.0000035 g/tick |
| Root (r=0.008, L=0.5) | `0.004 × π × 0.008² × 0.5` | 0.0000004 g/tick |
| Active meristem | `0.0005` (flat) | 0.0005 g/tick |

For the typical small plant (5 leaves avg 0.5 dm, 10 stems, 5 roots, 1 SAM + 1 RAM active):

```
Leaf maintenance:  5 × 0.0005   = 0.0025 g/tick
Stem maintenance:  10 × 0.0000035 = 0.000035 g/tick  (negligible)
Root maintenance:  5 × 0.0000004  = 0.000002 g/tick  (negligible)
Meristem maint:    2 × 0.0005    = 0.001 g/tick
──────────────────────────────────────────────────────
Total maintenance: ~0.0035 g/tick
```

**Leaf maintenance-to-production ratio: 0.002 / 0.02 = 10%.** Matches the comment in
`world_params.h:34` ("Leaf: ~10% of gross photosynthesis"). **Correct.**

Stem and root maintenance are essentially zero at juvenile scale. They become meaningful only
on mature trunks: a trunk at r=0.1 dm, L=1.0 dm costs `0.01 × π × 0.01 × 1.0 = 0.000314 g/tick`.
Still small (1.6% of a single mature leaf's output).

---

### Net Surplus

```
Production:   0.027 g/tick
Maintenance:  0.004 g/tick  (including 2 active meristems)
──────────────────────────────────────────
Net surplus:  ~0.023 g/tick   (6.8× safety margin)
```

**The plant is comfortably profitable on maintenance alone.** The actual drain is growth costs,
not respiration. A growing leaf at 0.005 dm/tick × 1.5 g/dm = 0.0075 g/tick per leaf. Five
growing leaves simultaneously cost 0.0375 g/tick — 1.4× total production. The plant must
queue growth behind supply. This is correct and biologically meaningful.

---

## 2. Vascular Sugar Flow

### Phloem Source Classification (`vascular.cpp:103–120`)

Sources (nodes that can export sugar):
1. **Seed** — exports above `phloem_reserve_fraction × cap`
2. **Leaves** — export above `phloem_reserve_fraction × cap + sugar_reserved_for_growth`

That's it. **Stem nodes, root nodes, and meristems are never phloem sources.**

### Phloem Sink Classification (`vascular.cpp:121–146`)

Sinks (nodes that demand sugar):
1. **Active meristems** (SAM, RAM) — demand up to `sugar_cap_meristem × meristem_sink_fraction`
   = `0.1 × 0.05 = 0.005 g/tick` each (after yesterday's cap fix)
2. **Starving nodes** (`starvation_ticks > 0`) — demand up to `cap × 0.5 - current`

No other nodes participate as active phloem sinks. Stems and roots receive sugar only via
local diffusion after vascular delivers to the nearest vascular conduit.

### The Three Deductions on a Leaf

For a growing leaf at size 1.0 dm (cap = 2.0g, current sugar = 1.2g):

```
1. Phloem reserve:  0.3 × 2.0g = 0.6g  protected
2. Growth reserve:  0.005 dm/tick × 1.5 g/dm = 0.0075g reserved
3. Available supply: 1.2 - 0.6 - 0.0075 = 0.5925g
```

At 60% of cap, a leaf still exports >30% of its cap per tick when sinks demand it.
**The three deductions do not eat all the supply.** A leaf is a meaningful source once it
reaches 50% of its cap, which happens within the first 10–20 ticks after the seed starts
funding it.

For a tiny seedling leaf (size 0.2 dm, cap = 0.1g minimum, sugar = 0.08g):

```
1. Phloem reserve:  0.3 × 0.1g = 0.03g
2. Growth reserve:  0.005 × 1.5 × water_gf ≈ 0.005g  (approximate)
3. Available:       0.08 - 0.03 - 0.005 = 0.045g
```

Even tiny leaves contribute a meaningful surplus when they have sugar. The problem is that
tiny leaves hold very little (cap floor = 0.1g) and produce very little (0.0008 g/tick at
size 0.2 dm). They depend almost entirely on seed funding during early growth.

---

## 3. Sugar Caps

### Meristem (`genome.h:262`)

`sugar_cap_meristem = 0.1g` (after fix from 2.0g)

With `meristem_sink_fraction = 0.05`:
- Demand per active meristem: `0.1 × 0.05 = 0.005 g/tick`
- As fraction of 5-leaf production: `0.005 / 0.025 = 20%` — healthy

Multiple simultaneous active meristems could still compete: 5 lateral buds activating at once
= 0.025 g/tick = 100% of small-plant production. This is biologically appropriate (bud break
is expensive) but watch for bursts during early canopy expansion.

### Stem/Root Nodes

Formula: `max(sugar_cap_minimum, π × r² × L × sugar_storage_density_wood)`

| Node | r (dm) | L (dm) | volume (dm³) | cap (g) |
|---|---|---|---|---|
| New internode | 0.015 | 0.01 (before elongation) | 0.0000071 | **0.1 (floor)** |
| Short stem | 0.015 | 0.5 | 0.000354 | **0.177** |
| Full-length stem | 0.015 | 1.0 | 0.000707 | **0.354** |
| Mature trunk | 0.05 | 1.0 | 0.00785 | **3.93** |
| Large trunk | 0.1 | 1.0 | 0.0314 | **15.7** |
| Thin root | 0.008 | 1.0 | 0.000201 | **0.1 (floor)** |

Observation: **All thin stems and roots (initial radius) hit the `sugar_cap_minimum = 0.1g`
floor.** The `sugar_storage_density_wood = 500 g/dm³` only matters for stems thicker than
~0.035 dm (3.5mm radius). This is consistent with the debug log: all juvenile nodes show
`sugar_cap = 0.100000`.

### Leaf

Formula: `max(sugar_cap_minimum, leaf_size² × sugar_storage_density_leaf)`

`sugar_storage_density_leaf = 2.0 g/dm²`

| leaf_size | area | cap |
|---|---|---|
| 0.2 dm (seedling) | 0.04 | **0.1 (floor)** |
| 0.5 dm (growing) | 0.25 | **0.5** |
| 1.0 dm | 1.0 | **2.0** |
| 1.5 dm (max) | 2.25 | **4.5** |

A max-size leaf holds 4.5g — 225 ticks of its own production. Good buffer.

### Seed

Formula: `max(seed_resource_cap(seed, sugar_cap, g), seed_sugar)`

`seed_sugar = 48g`. The seed cap is forced to at least 48g regardless of the tree size.
This makes the seed a permanent 48g reservoir. Phloem reserve on the seed:
`0.3 × 48 = 14.4g` always protected — never exported.

**The seed's exportable budget: 48 - 14.4 = 33.6g.** But exportable doesn't mean demanded.

---

## 4. Parameter Imbalances

### Problem 1: Seed Sugar Is a Frozen Reservoir

The debug log (tick 18168, 9541 nodes, large plant) shows the seed still holding ~44g at
the start of a run. The phloem pass can extract it in theory:
`surplus = sugar - phloem_reserve - growth_reserve = 44 - 14.4 - 0 = 29.6g exportable`.

But sinks only demand 0.005 g/tick (SAM) + starving nodes. A healthy plant with no starving
nodes drains the seed at:

```
44g ÷ 0.005 g/tick = 8800 ticks = 367 days
```

This is the "44g frozen in stems" problem from the storage plan. The seed's sugar cannot
reach the plant fast enough to matter because nobody is hungry enough to demand it.

The phloem_reserve_fraction also permanently locks 14.4g of the seed. That's sugar that
can never leave the seed. For a germinating seedling whose ONLY sugar source is the seed,
this is a significant inefficiency.

**Root cause: `seed_sugar` funds a sugar surplus nobody is hungry enough to demand. The
starch plan (seed_starch → GA-triggered mobilization) addresses this exactly.**

### Problem 2: `sugar_storage_density_wood = 500 g/dm³` Is 20× Unrealistic

Real wood stores starch at ~5% dry mass, with wood density ~500 g/dm³ → ~25 g glucose/dm³.
The sim uses 500 g/dm³ — a 20× overestimate.

In practice this doesn't cause visible bugs at juvenile scale (all small nodes hit the
0.1g floor anyway). It becomes problematic at mature trunk scale: a 1cm-radius trunk holds
15.7g with realistic maintenance of 0.000314 g/tick — 50,000 ticks to drain if totally cut
off. It functions as infinite storage, not a buffer.

The comment says "high cap so stems can pass sugar through." This is a transit-buffer
justification, not a storage-accuracy one. Fine as-is pre-starch, but post-starch the
sugar cap should reflect only the dissolved sugar fraction (small), with starch holding
the actual reserves. Separate starch_storage_density (~25 g/dm³) and sugar_cap density
(~25 g/dm³ dissolved) would be more accurate.

**This is not causing active bugs, but it will need revisiting when starch is implemented.**

### Problem 3: `starvation_ticks_max = 1200` Is Very Long

1200 hours = 50 days. A dead branch can coast on reserves for 50 days before dying and
being removed. During that time it pays maintenance (tiny for stems, but non-zero) and
acts as dead weight in the transport graph.

Real plants shed dead branches in days to weeks, not 50 days. A value of 168–336 (7–14
days) would match real abscission timescales better and keep the transport graph clean.

`leaf_abscission_ticks = 500` hours (20 days) is also on the long side relative to real
deciduous dynamics (7–10 days for stress response), but acceptable for a non-seasonal sim.

### Problem 4: Post-Fix Meristem Demand Is Correct, but Burst Risk Remains

Single active SAM after fix: `0.005 g/tick / 0.025 g/tick = 20%` of small-plant output.
Multiple simultaneous activations (lateral bud burst) could saturate the phloem, but that
is physically correct behavior (bud break IS expensive). Not a bug.

### Problem 5: Phloem Reserve Applied to Seed Is Wrong In Spirit

`phloem_reserve_fraction = 0.3` protects 30% of a node's cap as "structural glucose" —
appropriate for a leaf that needs reserves to stay alive. But the seed is an energy store
whose job is to be depleted. Protecting 14.4g permanently is wasteful.

The seed reserve fraction should be either 0.0 (export everything until empty) or governed
by GA-triggered mobilization (the starch plan). Currently the 14.4g is simply dead weight.

---

## 5. Specific Recommendations

### 5A. The Seed Reserve Problem (High Priority — Blocks Starch)

**Current:** `seed_sugar = 48g`, `phloem_reserve_fraction = 0.3` → 14.4g forever locked.

**Recommendation:** Don't tune this — implement the starch plan's `seed_starch` → GA
mobilization pathway. The seed reserve problem is the primary motivation for the starch epic.
Once seed_starch exists, `seed_sugar` drops to 0 and the 30% reserve lock disappears.

**Short-term workaround (if needed before starch):** Reduce `phloem_reserve_fraction` from
0.3 to 0.1 for seeds only (add a `seed_phloem_reserve` param). Math: 48g × 0.1 = 4.8g
protected (still safe), 43.2g exportable. Drain time: 43.2 / 0.005 = 8640 ticks — still
slow, but 3× faster.

### 5B. `starvation_ticks_max` (Low Priority — Polish)

**Current:** `1200 ticks = 50 days`

**Recommended:** `336 ticks = 14 days`

Math: 14 days is a realistic outer bound for maintaining living cells without sugar. Halving
this would accelerate dead-branch cleanup by 4× without affecting plants that aren't dying.

### 5C. `sugar_storage_density_wood` Post-Starch (Medium Priority — Deferred)

**Defer until starch is implemented.** At that point:
- `sugar_storage_density_wood` → 25–50 g/dm³ (dissolved sugar fraction only)
- `starch_storage_density` → 25–50 g/dm³ (starch fraction)
- Combined cap stays similar to current for large trunks (total = 50–100 g/dm³ vs
  current 500 g/dm³, but nodes small enough to hit the floor anyway at juvenile scale)

The net effect: large trunks store less mobile sugar but more starch. This makes them
naturally dependent on GA-triggered mobilization, which is the intended behavior.

### 5D. Stem Production Is Fine, Leave It

`stem_photosynthesis_rate = 0.005 g/(dm²·hr)` contributes ~8% of plant production at
juvenile scale. Biologically accurate, not worth adjusting.

---

## 6. Summary Table

| Parameter | Current | Status | Priority |
|---|---|---|---|
| `sugar_production_rate` | 0.02 g/(dm²·hr) | Correct (matches real C3) | — |
| `sugar_maintenance_leaf` | 0.002 g/(dm²·hr) | Correct (10% of production) | — |
| `sugar_maintenance_stem` | 0.01 g/(dm³·hr) | Correct | — |
| `sugar_cap_meristem` | 0.1g (fixed) | Good | — |
| `meristem_sink_fraction` | 0.05 | Good (0.005g/tick demand) | — |
| `phloem_reserve_fraction` | 0.3 | OK for leaves, wrong for seed | Post-starch |
| `seed_sugar` | 48g | Frozen reservoir — solve with starch | High |
| `sugar_storage_density_wood` | 500 g/dm³ | 20× unrealistic, benign now | Post-starch |
| `sugar_cap_minimum` | 0.1g | Correct floor | — |
| `starvation_ticks_max` | 1200 ticks | Too long, dead branches persist | Low |
| `stem_photosynthesis_rate` | 0.005 | Correct | — |

The plant's sugar economy is fundamentally sound at the rate level. The structural problems
are all about buffering and transit: too much sugar frozen in the seed reservoir, storage
density too high for large stems, and starvation tolerance too long. All three resolve
naturally once starch is implemented.
