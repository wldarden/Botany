# Sugar Realism: Volume-Based Maintenance & Storage Caps

**Date:** 2026-04-10
**Problem:** With light=0, a mature plant survives 370+ days on accumulated sugar. Real plants survive 2-4 weeks. Two root causes: (1) maintenance scales with radius, not tissue volume, so costs are ~21% of production instead of 40-60%; (2) sugar accumulates in nodes without bound.

## Change 1: Volume-Based Maintenance

### Formula Changes

| Type | Old | New |
|------|-----|-----|
| STEM | `rate * radius` | `rate * pi * r^2 * internode_length` |
| ROOT | `rate * radius` | `rate * pi * r^2 * internode_length` |
| LEAF | `rate * leaf_size` | `rate * leaf_size^2` (area) |
| Meristem | `rate` (flat) | unchanged |

`internode_length` = `glm::length(node.offset)`, floored at 0.01 dm for newly created nodes.

### New Default Values (Genome)

The units change — old values are meaningless under the new formula.

```
sugar_maintenance_stem     = 0.028    // g glucose / (dm^3 * hr) — wood is cheap to maintain
sugar_maintenance_root     = 0.135    // g glucose / (dm^3 * hr) — fine roots are expensive (high turnover)
sugar_maintenance_leaf     = 0.013    // g glucose / (dm^2 * hr) — leaves dominate maintenance budget
sugar_maintenance_meristem = 0.001    // g glucose / hr per active tip (unchanged)
```

### Expected Behavior

For a reference plant (200 stems, 100 leaves, 100 roots):
- Total maintenance: ~0.19 g/hr = **53% of gross production** (0.36 g/hr)
- As trunk thickens, volume grows quadratically with radius, so maintenance rises nonlinearly
- Mature trees naturally shift toward 60-70% maintenance ratio

### Files Changed

- `genome.h` — update default values, update comments to reflect new units
- `sugar.cpp` (`consume_sugar`) — new volume-based cost formula

## Change 2: Per-Node Sugar Storage Cap

### New Genome Parameters

```
sugar_storage_density_wood = 50.0   // g glucose max / dm^3 of stem/root tissue
sugar_storage_density_leaf = 0.5    // g glucose max / dm^2 of leaf area
sugar_cap_minimum          = 0.01   // g glucose — floor for tiny/new nodes
```

### Cap Formula

| Type | Formula |
|------|---------|
| STEM/ROOT | `max(minimum, pi * r^2 * internode_length * density_wood)` |
| LEAF | `max(minimum, leaf_size^2 * density_leaf)` |

### Helper Function

`float sugar_cap(const Node& node, const Genome& g)` — declared in `sugar.h`, used everywhere sugar is modified.

### Enforcement (three places)

1. **Production** (`produce_sugar`): skip production if `node.sugar >= cap`. Binary cutoff — diffusion handles smoothing by draining leaf sugar into stems between ticks.

2. **Diffusion** (`diffuse_sugar`): cap-aware transfers. When computing transfer amount, also clamp by receiver's headroom: `min(transfer, receiver_cap - receiver.sugar)`. Sugar that can't fit stays in sender.

3. **Safety clamp** (`consume_sugar`, end): `node.sugar = min(node.sugar, cap)`. Catches edge cases.

### Expected Behavior

Mature plant total storage capacity: ~85g.
At 0.19 g/hr maintenance drain with zero light: **~450 hours = 19 days**.
Thicker/older plants with more wood volume push toward 3-4 weeks.
Target: 2-4 weeks. Hit.

### Files Changed

- `genome.h` — three new params with defaults
- `sugar.h` — declare `sugar_cap()` helper
- `sugar.cpp` — implement `sugar_cap()`, enforce cap in `produce_sugar`, `diffuse_sugar`, `consume_sugar`

## Test Updates

- `test_sugar.cpp` — update existing maintenance tests for volume-based formula, add cap enforcement tests (production stops at cap, diffusion respects cap, clamp works)

## CLAUDE.md Update

- Update Sugar Model section to document volume-based maintenance and storage caps
