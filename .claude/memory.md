# Session State - Hormone-Driven Branching (April 2026)

## What Was Done This Session

### Root meristem cap
- Added `root_meristem_count_` to Plant, incremented in `create_meristem` for ROOT_APICAL/ROOT_AXILLARY
- Hard cap at 100 (`Plant::max_root_meristems`), checked via `root_meristems_at_cap()`
- Cap gates new axillary bud creation during chain growth but NOT activation of existing dormant buds

### Auxin transport rewrite
- **Before**: single-pass postorder, auxin accumulated across ticks at the base (wrong - base had more than tip)
- **After**: reset to 0 each tick, then two-phase transport:
  1. `auxin_collect` (postorder): produce at shoot tips, flow basipetally to parent
  2. `auxin_spillback` (preorder): small fraction (10%) redistributes from junctions back into children
- This gives correct gradient: tip highest, falls off with distance

### Cytokinin transport
- Also reset to 0 each tick now
- Two-phase: collect from root tips toward seed, then distribute downward to children

### Shoot axillary activation fix
- Buds now check **parent's** auxin (the stem node they sit on), not their own node
- Removed cytokinin requirement - shoot buds activate on low auxin only
- Root axillaries check parent's cytokinin only
- Old bug: cytokinin never reached shoot nodes in meaningful amounts (diluted across branches), so dual requirement meant buds never activated

### `--color type` mode
- Green = STEM (shoot), orange = ROOT
- Added `color_by_type_` flag to Renderer

## Current Genome Defaults
- auxin_production=1.0, transport=0.3, spillback=0.1, decay=0.05, threshold=0.15
- cytokinin_production=1.0, transport=0.3, decay=0.05, threshold=0.15
- root_initial_radius=0.025 (separate from shoot initial_radius=0.05)
- root_branch_angle=0.35 rad (~20 deg, tighter than shoot's 45 deg)

## Open Questions
- User wants to potentially add cytokinin check back to shoot axillary activation ("activate when auxin low AND cytokinin low") but needs cytokinin delivery to shoot nodes fixed first
- Current shoot activation is auxin-only which works well visually
