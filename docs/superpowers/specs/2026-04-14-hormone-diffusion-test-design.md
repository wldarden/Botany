# Hormone Diffusion Test App

**Date:** 2026-04-14
**Purpose:** Standalone app to measure hormone transport behavior across different tree sizes, enabling quick parameter tuning by visualizing steady-state gradients.

## Overview

A CLI app that builds 3 hand-crafted static trees (young, medium, large), runs `Plant::tick()` with growth frozen but hormone production active, and reports per-node chemical levels at steady state. Tests one chemical at a time: Auxin, Cytokinin, Gibberellin, or Stress.

## Mock Trees

All trees are hand-built with deterministic topology. Radii taper toward tips. Nodes get spatial positions for consistency but transport is tree-graph only.

### Young Tree (~8 nodes)

Linear chain, one leaf. Tests: does the hormone reach the other end?

```
seed(STEM) ── stem ── stem ── shoot_apical
                        └── leaf
  └── root ── root ── root_apical
```

### Medium Tree (~20 nodes)

One shoot branch, one root branch. Tests: does the hormone split correctly at junctions?

```
seed(STEM) ── stem_1 ── stem_2 ── stem_3 ── shoot_apical_1
                │          │         └── leaf_1
                │          └── stem_4 ── stem_5 ── shoot_apical_2 (branch)
                │                         └── leaf_2
                └── leaf_0 (basal)
  └── root_1 ── root_2 ── root_apical_1
        └── root_3 ── root_apical_2 (branch)
```

### Large Tree (~40-50 nodes)

Multiple branching levels on both sides. Tests: can the hormone penetrate a deep tree? Where does it stall?

```
seed(STEM)
  ├── [main shoot: 5-6 stem segments, 3 branches at different heights, each with 2-3 stems + leaf + apical]
  │   branches get progressively thinner
  │   leaves at each branch tip and along main axis
  └── [root system: 3-4 root segments, 2 branches, each 2-3 roots deep + apical]
```

Exact topology defined in code. ~40-50 nodes total.

## Tree Setup

For each run:

1. **Build tree** from the chosen topology
2. **Freeze growth** — set genome overrides:
   - `growth_rate = 0`, `thickening_rate = 0`, `leaf_growth_rate = 0`
   - Plastochrons set very high (prevent new node creation)
3. **Seed sugar** on all nodes (enough for normal meristem production; avoids Michaelis-Menten starvation floor)
4. **Set light_exposure** on all leaves (e.g., 0.8) — persists since Plant::tick() doesn't recompute it
5. **For stress tests** — manually set stress values on target stem nodes
6. **Zero the target chemical** on all nodes before starting

## Per-Chemical Setup

### Auxin
- **Source:** Shoot apical meristems (produced in `ShootApicalNode::tick()`)
- **Direction:** Basipetal (bias = -0.1, toward root)
- **Key measurement:** Gradient from apex to root; level at axillary bud parents vs `auxin_threshold`
- **Requires:** Sugar on meristems (Michaelis-Menten gating; 10% floor without sugar)

### Cytokinin
- **Source:** Leaf nodes (proportional to `sugar_produced * cytokinin_production_rate`)
- **Direction:** Acropetal (bias = +0.1, toward tips)
- **Key measurement:** Level at growing tips; whether it crosses `cytokinin_growth_threshold`
- **Requires:** Working photosynthesis — leaves need `light_exposure` set and sugar capacity available

### Gibberellin
- **Source:** Young leaf nodes (age < `ga_leaf_age_max`, proportional to `leaf_size`)
- **Direction:** Bidirectional (bias = 0, local effect)
- **Key measurement:** Level on parent/grandparent stem nodes; whether it meaningfully boosts elongation
- **Requires:** Leaves with `leaf_size > 0` and `age < ga_leaf_age_max`

### Stress Hormone
- **Source:** Stem nodes under mechanical load (above `stress_hormone_threshold`)
- **Direction:** Transported via tree graph with genome-defined params
- **Key measurement:** How far the stress signal propagates; level at nearby meristems
- **Requires:** Stress values manually set on target stem nodes (since we're not computing real physics)

## Tick Loop

- Call `Plant::tick(world_params)` each tick — this runs production, transport, and decay naturally
- After each tick, record the target chemical's level on every node
- **Steady-state detection:** When the max per-node absolute delta between consecutive ticks falls below epsilon (e.g., 1e-4), declare steady state reached
- Cap at `--max-ticks` (default 200) if steady state isn't reached

## Output

### Console Table (default)

One table per (tree, chemical) combination. Sorted by tree depth, shoot side first, then root side.

```
=== AUXIN | Medium Tree | default genome | steady state @ tick 47 ===
Node            Type            Depth  Level    Gradient  vs Threshold
shoot_tip_1     SHOOT_APICAL    4      0.482    ---       ---
stem_3          STEM            3      0.341    -0.141    ---
leaf_1          LEAF            4      0.289    ---       ---
stem_2          STEM            2      0.218    -0.123    ---
axillary_1      SHOOT_AXILLARY  3      0.195    ---       ABOVE (0.15)
stem_1          STEM            1      0.134    -0.084    ---
seed            STEM            0      0.078    -0.056    ---
root_1          ROOT            1      0.031    -0.047    ---
root_tip        ROOT_APICAL     2      0.009    -0.022    ---

Summary: apex=0.482  base=0.009  ratio=53.6x  axillary buds above threshold: 1/1
         Time to 90% steady state: 31 ticks
```

**Gradient column:** Difference from parent node's level. Only shown for nodes on the main axis (not branches/leaves).

**vs Threshold column:** For axillary buds only — shows ABOVE/BELOW relative to the relevant threshold (`auxin_threshold` for shoot axillary, `cytokinin_threshold` for root axillary).

**Summary line:** Apex level, base level, apex:base ratio, how many axillary buds are above/below activation threshold, and ticks to reach 90% of final steady-state values.

### CSV (--csv flag)

Columns: `chemical,tree,tick,node_name,node_type,depth,level`

One row per node per final tick (steady state). Easy to diff across param sets or plot in a spreadsheet.

## CLI

```
./build/botany_hormone_test [options]
  --chemical auxin|cytokinin|gibberellin|stress  (default: auxin)
  --tree young|medium|large|all                   (default: all)
  --max-ticks N                                   (default: 200)
  --csv                                           (output CSV instead of table)
```

Genome and WorldParams are hardcoded to defaults initially. Edit source to experiment with different values.

## File

`src/app_hormone_test.cpp` — single file, added to CMakeLists.txt as a new executable target (`botany_hormone_test`).

## Notes on DFS Walk Order

Plant::tick() walks the tree seed → children (DFS pre-order). Nodes closer to seed process transport first. This means:
- Auxin produced at the shoot apex takes multiple ticks to propagate down to the base
- Changes at the base affect tips faster (processed earlier in walk)

This is real behavior, not an artifact — the test output naturally reflects it.

## What This Does NOT Do

- No growth or structural changes during the test
- No sugar transport testing (separate `app_sugar_test` handles that)
- One chemical per run (no cross-chemical interactions)
- No Engine::tick() (avoids shadow casting overhead)
- No ethylene (uses spatial diffusion, not tree transport)
