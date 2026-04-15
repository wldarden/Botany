# Auxin Growth Sensitivity Design

**Date:** 2026-04-15
**Status:** Approved

## Summary

Auxin directly modulates growth rates in all five tissue types via a saturating (Michaelis-Menten) multiplier. Each tissue type gets its own signed sensitivity parameters, allowing auxin to promote growth in shoots and leaves while inhibiting it in roots. The multiplier chains into existing effective-rate computations alongside GA, ethylene, and stress modifiers — auxin shapes desire, sugar gates ability, physical caps remain unchanged.

## Motivation

Real auxin promotes cell elongation in shoots (the classic acid-growth mechanism) while inhibiting root elongation at the same concentrations — roots are far more sensitive. Currently the sim models auxin's role in apical dominance (threshold-based bud inhibition) and cambial thickening, but not its direct effect on elongation and expansion rates. Adding this creates more realistic differential growth between tissue types and gives evolution another axis to optimize.

## Design

### Helper Function

One free function in `src/engine/node/meristems/helpers.h`:

```cpp
inline float auxin_growth_factor(float auxin, float max_boost, float half_sat) {
    if (std::abs(max_boost) < 1e-8f) return 1.0f;
    float saturation = auxin / (auxin + std::max(half_sat, 1e-6f));
    return 1.0f + max_boost * saturation;
}
```

- At zero auxin: returns 1.0 (no effect).
- At infinite auxin: asymptotes to `1.0 + max_boost`.
- Signed `max_boost`: positive = promotion, negative = inhibition.
- `half_sat` controls how much auxin is needed for half-max effect.
- Saturation prevents runaway growth — receptor saturation is biologically real.

### Genome Fields

10 new fields in `Genome` (5 tissue types x 2 params each):

| Field | Default | Description |
|-------|---------|-------------|
| `stem_auxin_max_boost` | +0.5 | Stems: auxin promotes elongation by up to 50% |
| `stem_auxin_half_saturation` | 0.2 | Auxin level for half-max stem effect |
| `root_auxin_max_boost` | -0.3 | Roots: auxin inhibits elongation by up to 30% |
| `root_auxin_half_saturation` | 0.1 | Roots are very sensitive — low half-sat |
| `leaf_auxin_max_boost` | +0.3 | Leaves: auxin promotes expansion by up to 30% |
| `leaf_auxin_half_saturation` | 0.2 | Auxin level for half-max leaf effect |
| `apical_auxin_max_boost` | +0.2 | Shoot tips: mild promotion of tip extension |
| `apical_auxin_half_saturation` | 0.3 | High half-sat — apicals sit in high auxin, need more for effect |
| `root_apical_auxin_max_boost` | -0.2 | Root tips: mild inhibition of tip extension |
| `root_apical_auxin_half_saturation` | 0.1 | Root tips are sensitive |

All values are signed. Evolution can flip any tissue's response (e.g., discover root promotion). Ranges for evolution: max_boost `[-1.0, 2.0]`, half_saturation `[0.01, 1.0]`.

### Application Sites

Each tissue type calls `auxin_growth_factor()` and multiplies it into its existing growth rate computation:

**StemNode::elongate()** — multiplier joins GA boost, ethylene inhibit, stress inhibit:
```
effective_rate = base_rate * ga_boost * eth_inhibit * stress_inhibit * auxin_boost
```

**RootNode::elongate()** — same chain, default params are negative (inhibition):
```
effective_rate = base_rate * ga_boost * eth_inhibit * auxin_boost
```

**LeafNode::grow_size()** — multiplier scales the max_growth rate before sugar gating:
```
max_growth = leaf_growth_rate * auxin_boost
```

**ApicalNode::grow_tip()** — multiplier scales growth rate inside sugar-gated section:
```
actual_rate = growth_rate * gf * auxin_boost
```

**RootApicalNode::grow_tip()** — same pattern, default negative:
```
actual_rate = root_growth_rate * gf * auxin_boost
```

### What Is NOT Touched

- `StemNode::thicken()` and `RootNode::thicken()` — already have auxin-driven cambial growth via `auxin_thickening_threshold`. No double-dipping.
- `ApicalNode::can_activate()` — threshold-based inhibition is a separate mechanism from growth rate modulation.
- Sugar costs, physical caps (max_len, max_leaf_size), and the transport system are unchanged.

### Evolution Integration

All 10 fields registered in `genome_bridge.cpp` in the "auxin" linkage group. Each gets a gene entry with valid range and mutation strength (% of range).

### Testing

Unit test for the helper function:
- `auxin_growth_factor(0, ...)` returns 1.0
- Positive `max_boost` with moderate auxin returns value in `(1.0, 1.0 + max_boost)`
- Negative `max_boost` returns value in `(1.0 + max_boost, 1.0)`
- Very high auxin asymptotes near `1.0 + max_boost`

Existing sim and sugar_test app validate that the economy doesn't break.

## Files Changed

- `src/engine/genome.h` — 10 new fields + defaults
- `src/engine/node/meristems/helpers.h` — `auxin_growth_factor()` helper
- `src/engine/node/stem_node.cpp` — 1 line in `elongate()`
- `src/engine/node/root_node.cpp` — 1 line in `elongate()`
- `src/engine/node/tissues/leaf.cpp` — 1 line in `grow_size()`
- `src/engine/node/tissues/apical.cpp` — 1 line in `grow_tip()`
- `src/engine/node/tissues/root_apical.cpp` — 1 line in `grow_tip()`
- `src/evolution/genome_bridge.cpp` — 10 gene registrations
- `tests/` — new unit test for `auxin_growth_factor()`
