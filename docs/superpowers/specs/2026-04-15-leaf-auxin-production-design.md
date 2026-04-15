# Leaf Auxin Production (Auxin Enhancement #1)

Growth-coupled auxin production from young leaves and growth-modulated apical production.

## Context

Currently only `ApicalNode::produce_auxin()` generates auxin, modulated by light, sugar, and age. In real plants, young expanding leaves are a significant secondary auxin source. Auxin production is tied to active cell division and expansion -- growing tissue produces auxin, quiescent tissue does not.

This is the first of several planned auxin enhancements (see `docs/auxin_enhancement.txt`).

## Design

### New genome parameters

**Rename:** `auxin_production_rate` -> `apical_auxin_baseline`

Same default value (0.15), same meaning -- the base auxin output of a shoot apical meristem per tick before environmental modulation.

**New parameters:**

| Field | Default | Evolvable range | Description |
|---|---|---|---|
| `apical_growth_auxin_multiplier` | 2.0 | [0, inf) | Multiplier on apical baseline from active growth. At max growth: total = baseline * (1 + multiplier). |
| `leaf_auxin_baseline` | 0.15 | (0, inf) | Scaling constant for leaf auxin production (decoupled from apical rate). |
| `leaf_growth_auxin_multiplier` | 0.1 | [0, 1] | Fraction of leaf_auxin_baseline produced at max leaf growth per tick. |

### Apical auxin production formula

Existing environmental modulations (light, sugar, age) remain unchanged. Growth adds a multiplier on top:

```
modulated_baseline = apical_auxin_baseline * light_factor * sugar_factor * age_factor

growth_fraction = growth_gf(sugar, cytokinin)   // 0-1, same inputs as grow_tip()

total = modulated_baseline * (1 + apical_growth_auxin_multiplier * growth_fraction)
```

`growth_fraction` is computed from current sugar/cytokinin levels using the existing `growth_fraction()` helper from `meristems/helpers.h`, called before `produce_auxin()` in `tissue_tick()`. This avoids reordering the tick pipeline -- we just read the same state that `grow_tip()` will read later.

**Magnitudes at default genome:**
- No growth: `baseline * 1.0` (unchanged from today)
- Max growth: `baseline * 3.0` (with default multiplier of 2.0)

### Leaf auxin production formula

Computed inside `LeafNode::grow_size()` after the actual growth amount is determined:

```
growth_fraction = actual_growth / leaf_growth_rate   // 0-1
leaf_auxin = growth_fraction * leaf_growth_auxin_multiplier * leaf_auxin_baseline
```

No growth (full size, stressed, starved) -> zero auxin. Proportional to how fast the leaf is expanding. A leaf that stops growing due to stress stops producing auxin; if growth resumes, auxin production resumes.

**Magnitudes at default genome:**
- Single leaf at max growth: `1.0 * 0.1 * 0.15 = 0.015` auxin/tick (10% of apical baseline)
- Branch with 5 growing leaves: `0.075` auxin/tick (50% of apical baseline)
- Full-grown leaf: `0.0` auxin/tick

### Emergent behavior

- Branches with many actively growing leaves export more auxin at junctions, strengthening their dominance signal relative to branches with fewer/stalled leaves.
- Seasonal growth flushes produce temporary auxin surges from leaves, reinforcing apical dominance during active growth periods.
- A branch under sugar stress (leaves stop growing) loses its leaf-auxin contribution, weakening its competitive position.

## Files changed

### `src/engine/genome.h`
- Rename `auxin_production_rate` -> `apical_auxin_baseline`
- Add `apical_growth_auxin_multiplier` (default 2.0)
- Add `leaf_auxin_baseline` (default 0.15)
- Add `leaf_growth_auxin_multiplier` (default 0.1)
- Update `default_genome()` with new fields and defaults

### `src/engine/node/tissues/apical.cpp`
- `tissue_tick()`: compute `growth_gf` (using `growth_fraction()` helper with current sugar/cytokinin) before calling `produce_auxin()`, pass it as a parameter
- `produce_auxin()`: add `float growth_gf` parameter, multiply result by `(1 + apical_growth_auxin_multiplier * growth_gf)`
- Rename internal reference from `g.auxin_production_rate` to `g.apical_auxin_baseline`

### `src/engine/node/tissues/apical.h`
- Update `produce_auxin()` signature to accept `float growth_gf`

### `src/engine/node/tissues/leaf.cpp`
- `grow_size()`: after computing actual growth amount, add auxin production: `chemical(Auxin) += (growth / g.leaf_growth_rate) * g.leaf_growth_auxin_multiplier * g.leaf_auxin_baseline`

### `src/evolution/genome_bridge.cpp`
- Rename `"auxin_production_rate"` gene to `"apical_auxin_baseline"` (same range: 0.01-2.0)
- Add `"apical_growth_auxin_multiplier"` gene (range: 0.0-10.0) in auxin linkage group
- Add `"leaf_auxin_baseline"` gene (range: 0.01-2.0) in auxin linkage group
- Add `"leaf_growth_auxin_multiplier"` gene (range: 0.0-1.0) in auxin linkage group
- Update `from_structured()` to read new gene names
- Update auxin linkage group to include new genes

### `docs/CLAUDE.md`
- Update auxin production description to mention leaf production and growth modulation
- Update tuning parameters section with new genome fields

### Tests
- Verify leaf auxin production scales with growth (growing leaf > 0, full-size leaf = 0)
- Verify apical growth multiplier: production at max growth > production at zero growth
- Verify full-size leaf produces zero auxin
- Verify leaf auxin resumes if growth resumes after stress
