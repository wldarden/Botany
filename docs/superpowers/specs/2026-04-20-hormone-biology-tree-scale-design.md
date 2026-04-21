# Hormone Biology at Tree Scale — Design

**Date:** 2026-04-20
**Goal:** Restructure the auxin and cytokinin systems so that plants behave biologically correctly across the full range from seedling (~1 dm) to mature tree (~10 m). Current system works for small plants but has fundamental scaling failures that no amount of parameter tuning can fix.

---

## 1. Motivation

The current model was tuned for seedling-scale plants and has three architectural problems that prevent it from scaling:

1. **Auxin is used as a long-distance signal.** Auxin decays at 12%/tick and moves one hop per tick via PIN. Over 30 hops, signal drops to 2% of source; over 100 hops (a 10m tree), to ~0%. Any mechanism that requires shoot-produced auxin to reach a root tip cannot work at scale.

2. **Cytokinin production is gated on local auxin at the RA.** Primary RA has near-zero auxin (because of #1), so it produces near-zero CK. The rest of the plant's CK-dependent behaviors (SAM growth, lateral activation) starve for lack of source.

3. **Root elongation is gated on external auxin reaching the RA.** Same scaling failure as #1 — works for a seedling, fails for a tree.

In real plants, auxin is genuinely short-range (localized apical dominance, PIN canalization within a few internodes). The long-distance shoot↔root communication channel is **sugar via phloem** for "down" signals and **cytokinin via xylem** for "up" signals. Both are bulk-flow pressure systems that scale to 100m+ trees in nature.

This spec repositions the model to match that biology.

---

## 2. Biological Background

Key findings from tree physiology literature (Domagalska & Leyser 2011; Barbier et al. 2019; Mason et al. 2014; cited here as summaries, not verbatim):

- **Auxin polar transport is slow and short-range.** PIN-mediated flow moves at roughly 1 cm/hr. Auxin conjugate pools buffer short-term flux. Effective signaling range from a single SAM: 5–20 cm (tens of internodes at most). Does not reach root tips of tall trees.
- **Cytokinin is synthesized primarily at root tips**, driven by local metabolic state (sugar, water, nutrients, nitrogen). It is not auxin-gated in real plants.
- **Cytokinin travels up in the xylem stream** via transpiration pull and root pressure. A tall tree's xylem can move CK from root tip to canopy in hours.
- **Lateral bud dormancy is a multi-factor system:**
  1. **Auxin canalization (topological):** a dormant bud's own auxin must canalize a flow path to the main trunk. If the main path is strongly canalized, the bud's auxin has no sink and the bud's own PINs never saturate — it stays dormant. If main flow weakens, laterals escape.
  2. **Cytokinin (permissive):** CK above threshold is required but not sufficient. "The plant has root capacity to support more growth."
  3. **Sugar (energy):** buds cannot grow without an energy supply; recent work has shown sugar depletion is the immediate trigger after decapitation.
- **Root lateral branching is inhibited by local cytokinin.** Too much CK = the root system is already signaling "enough roots." Low CK = space to branch.
- **Apical dominance weakens with distance** naturally in real trees. A lateral bud 30 cm from the apex escapes dominance because auxin no longer reaches it in meaningful concentration. Trees with long internodes (or damaged apices) branch prolifically; trees with compressed apices (many internodes packed into a short stem) suppress all laterals.

---

## 3. Architecture Principles

Three invariants the redesigned system must preserve:

### Principle 1 — Local signals stay local, long-distance signals travel in vascular.

| Signal | Role | Transport | Range |
|---|---|---|---|
| Auxin | Apical dominance, canalization, tropism | PIN (polar, active) | ~5–20 internodes |
| Cytokinin | "Root system healthy" permissive signal | Xylem bulk flow (Jacobi) | Whole plant |
| Sugar | Energy & growth throttle | Phloem bulk flow (Jacobi) | Whole plant |
| Water | Turgor, transpiration | Xylem bulk flow (Jacobi) | Whole plant |

Any behavior that needs to scale to 10m must use sugar or a vascular-transported chemical, never auxin.

### Principle 2 — Growth is sugar-rate-limited, not hormone-gated.

Root growth, shoot growth, and leaf expansion all proceed as fast as their local sugar supply allows, modulated by hormone signals that don't need to be present in any minimum quantity — they set the *rate*, not a *go/no-go* gate.

Sugar supply is set by the phloem network. Root-shoot balance self-regulates: a plant with too many roots runs low on sugar per tip → roots slow → shoot grows more → more leaves → sugar supply increases → roots resume. No external controller needed.

### Principle 3 — Activation (dormancy break) is multi-factor.

No single chemical threshold activates a dormant bud. The conjunction of three conditions is required:

- **Canalization precondition:** the parent's main outgoing flow is weak enough that this bud's own PIN output could establish a path.
- **Permissive signal:** local CK above a threshold indicates the root system is healthy.
- **Energy signal:** local sugar above the activation cost.

For root laterals, the CK condition is inverted (high CK suppresses branching), and a small amount of local auxin is required for meristem initiation.

---

## 4. Component-Level Changes

### 4.1 Cytokinin production at RA

**Current:**
```cpp
cyto_produced = root_cytokinin_production_rate * local_auxin * mf_cyto(sugar, water);
```

**Proposed:**
```cpp
cyto_produced = root_cytokinin_production_rate * mf_cyto(sugar, water);
```

Remove the `local_auxin` factor. Root tips produce CK based on their own metabolic state (sugar and water availability). This matches biology and eliminates the dependency on long-distance auxin.

Dormant RAs still produce nothing (the `if (!active) return;` guard upstream is unchanged).

### 4.2 RA auxin self-production

**Current:** `root_tip_auxin_production_rate = 0.0` (temporarily disabled while testing the gate).

**Proposed:** `root_tip_auxin_production_rate = 0.002f` (~1/75 of SAM rate). This is the small floor real root tips produce for PIN recycling and lateral root initiation. It's not meant to signal anything long-distance; it just maintains the RA's local auxin compartment for canalization and for the tiny sub-gradient that drives lateral root initiation.

### 4.3 RA elongation gate

**Current:**
```cpp
gf = growth_fraction(sugar, max_cost, local_auxin, root_auxin_growth_threshold);
```

**Proposed:**
```cpp
gf = growth_fraction(sugar, max_cost, local_cytokinin, root_ck_growth_floor);
```

Swap the auxin modulator for a cytokinin modulator. Biological reasoning: the root tip needs to sense that the broader root system is metabolically healthy (represented by its own CK level) before committing to further elongation. If CK drops to zero (metabolic crisis), elongation stops. If CK is above the floor, the root grows at sugar-rate-limited speed.

`root_ck_growth_floor` is a low Km (~0.001) — enough to be detectable but not a hard gate. This ensures small amounts of CK permit growth, avoiding the cold-start bootstrap problem the auxin gate had.

### 4.4 SA elongation gate

Unchanged. Current code already uses `growth_fraction(sugar, max_cost, local_cytokinin, cytokinin_growth_threshold)` — this is correct. The problem was CK never reaching the SA, not the gate itself. Fix 4.1 restores CK supply at the source.

### 4.5 Lateral SA activation

**Current:**
```cpp
if (parent_auxin >= auxin_threshold) return false;   // apical dominance
if (local_cytokinin < cytokinin_threshold) return false; // permissive
if (local_sugar < sugar_cost_activation) return false; // energy
```

Keep unchanged. This already encodes the multi-factor model. Validate that the three thresholds are tuned sensibly once CK supply is restored.

### 4.6 Lateral RA activation

**Current:**
```cpp
if (local_auxin < root_auxin_activation_threshold) return false; // needs initiation auxin
if (local_cytokinin > root_cytokinin_inhibition_threshold) return false; // too much CK suppresses
if (local_sugar < sugar_cost_activation) return false; // energy
```

Keep unchanged. Validate and retune thresholds once CK now flows properly (previously CK inhibition threshold was irrelevant because CK was always 0).

### 4.7 Apical dominance scope

Unchanged mechanism, but state explicitly in docs that apical dominance is designed to be **short-range by construction**. The decay rate and transport rate of auxin dictate its range; tuning these changes the spatial extent of dominance but not its nature.

At tree scale, expect apical dominance to suppress laterals within ~5–10 internodes of an active SAM and for distant laterals to escape and become secondary branches naturally. This is the biologically correct behavior.

### 4.8 Diffuse-auxin-across-seed-junction

Keep the junction bridge (`diffuse_auxin_across_seed_junction`). It ensures shoot-produced auxin reaches the top of the root system for the initial few internodes of root auxin flow, which drives early lateral root initiation. The amount that crosses is small and decays within a few root-hops — consistent with the principle that auxin is a local signal.

### 4.9 Local diffusion skip on vascular chemicals

Keep the CK skip on edges touching vascular conduits (already implemented). This prevents local diffusion from draining RA-produced CK into the parent ROOT's local pool where it would be stuck.

---

## 5. Scaling Validation Plan

Three checkpoints, each a headless sim with the chain-profile dump. Run after each implementation step to validate that scaling invariants hold.

### 5.1 Seedling checkpoint (current baseline)

- Plant state: ~150 nodes, ~1 dm height, ~1 dm root depth.
- Ticks: 1000.
- Pass criteria:
  - Primary SA has `local_cytokinin > 0.001` and is growing.
  - Primary RA is actively elongating at a sugar-limited rate.
  - At least 50% of the plant's RAs are active (not dormant).
  - Visible auxin gradient along the shoot chain (SAM end > seed end).
  - Visible CK gradient along the xylem (root side higher than shoot side in xylem pools).

### 5.2 Adolescent checkpoint

- Plant state: ~500 nodes, ~3 m height, ~2 m root depth.
- Ticks: 10000.
- Pass criteria:
  - Plant does not collapse: primary SA still alive and growing.
  - CK delivery to SA scales sub-linearly with tree size (more RAs → more CK, but some attenuation by xylem length is acceptable).
  - Root growth rate stays biologically reasonable (not runaway, not frozen).
  - Lateral break pattern: dense near tips, sparse on main trunk (apical dominance functioning locally).

### 5.3 Mature tree checkpoint

- Plant state: ~2000+ nodes, ~10 m height, ~5 m root depth.
- Ticks: 40000+.
- Pass criteria:
  - Primary SA still alive.
  - Sugar economy balanced: total sugar produced ≈ total sugar consumed over any 24h window.
  - Root tips (RAs) at depth >3m still active and producing CK.
  - CK xylem pressure at seed junction is non-trivial (`> 0.01` AU per pool).
  - Tree does not have pathological pattern like "all growth at trunk base" or "nothing lateral."

If the adolescent or mature checkpoint fails, revisit tuning rather than architecture. If architecture needs revision, reopen this spec.

---

## 6. Tuning Strategy

New tuning parameter added:

| Parameter | Default | Role |
|---|---|---|
| `root_ck_growth_floor` | 0.001 | Km for CK-modulated RA elongation. Low enough that cold-start doesn't stall growth. |

Parameters that will need retuning after the changes:

| Parameter | Consideration |
|---|---|
| `root_cytokinin_production_rate` | Drives total CK supply. Must be high enough that ~20 active RAs produce CK at a rate that keeps SAM's xylem-delivered CK above `cytokinin_growth_threshold`. Expect to raise from 0.15 to ~0.5–1.0. |
| `root_tip_auxin_production_rate` | Reset to 0.002 (small PIN-recycling floor). Not expected to need further tuning. |
| `cytokinin_decay_rate` | Currently 0.05. Only decays in `local()`; xylem has no decay. Keep unchanged. |
| `cytokinin_growth_threshold` | Km for SA growth modulation by CK. Validate it's consistent with expected xylem CK delivery rate. |
| `root_auxin_activation_threshold` | Lower to 0.01 (below the new RA auxin self-equilibrium of ~0.017). Lets deep lateral RAs activate on their own auxin floor without needing shoot-delivered auxin. See §7.2. |
| `root_cytokinin_inhibition_threshold` | Validate after CK supply is restored. With CK now actually flowing, the inhibition mechanism can engage; previously it was dead code. |
| `auxin_diffusion_rate` | Keep low (0.05) by design — auxin is short-range. |
| `auxin_decay_rate` | Keep at 0.12 — biologically correct half-life. |

Retuning approach: run the seedling checkpoint, measure equilibrium CK at the primary SA, compare to `cytokinin_growth_threshold`. Adjust `root_cytokinin_production_rate` until SA CK sits at ~2–5× `cytokinin_growth_threshold` at equilibrium. Then run the adolescent checkpoint and confirm CK supply still meets threshold at scale.

---

## 7. Risks and Open Questions

### 7.1 CK xylem throughput at 10m

The xylem Jacobi moves `N=25` hops per tick. A 10m tree may have 100–200 hops from deepest RA to tallest SA. Propagation will take multiple ticks to fully equilibrate, but per-tick throughput is conserved (Jacobi doesn't lose mass, just smears it).

**Open question:** does `world.vascular_substeps = 25` need to scale with plant height? Possible follow-up: auto-scale N based on longest chain.

### 7.2 Root tip auxin production for lateral initiation

Lateral RA initiation requires local auxin ≥ `root_auxin_activation_threshold` (0.05). With `root_tip_auxin_production_rate = 0.002` and decay 0.12, self-equilibrium = 0.017 — *below* the activation threshold.

This means lateral RAs will only activate when they receive a boost of shoot-derived auxin. In small plants this works (shoot auxin reaches top few root levels). In large trees, deep lateral RAs won't ever activate.

**Mitigation:** lower `root_auxin_activation_threshold` to 0.01 (below self-equilibrium) so deep laterals can activate on their own auxin floor. This matches real biology where root branching is observed at all depths, not just near the stem.

### 7.3 Sugar distribution at scale

Not in scope for this spec, but worth flagging: sugar partitioning between competing sinks becomes non-trivial at tree scale. If phloem Jacobi distributes poorly to distant RAs, the root system will preferentially grow near the trunk and starve deep tips. Validate at the adolescent checkpoint; if a problem, that's a separate phloem tuning spec.

### 7.4 Cambium thickening at tree scale

Cambium responsiveness is currently calibrated so that a trunk with `auxin_flow_bias ~1.0` thickens at ~0.0002 dm/hr. Over 40000 ticks that's ~8 dm radius = ~16 cm diameter. Reasonable for a mature tree. Validate at mature checkpoint.

---

## 8. Out of Scope / Deferred

- **Juvenile/mature phase transitions.** All behavior remains emergent from physics. If we discover a specific realism gap (e.g., tree branching pattern changes dramatically in nature at maturity), we'll add it in a future spec.
- **Strigolactones.** Real lateral inhibition has a strigolactone component. Redundant with the auxin + CK + sugar model we're keeping.
- **Seasonality and dormancy cycles.** No winter behavior yet.
- **Sugar partitioning optimization.** May need its own spec if the adolescent checkpoint reveals problems.
- **Auto-scaling `vascular_substeps`** with plant height.
- **Strigolactone, ABA, jasmonate, brassinosteroid, ethylene's cross-talk with auxin/CK.** Ethylene is modeled as a stress signal, not as hormone crosstalk. Adding more hormones is a future direction.

---

## 9. Summary of Code Touchpoints

| File | Change |
|---|---|
| `src/engine/node/tissues/root_apical.cpp` | Remove `local_auxin` factor from CK production; swap auxin→CK in elongate gate |
| `src/engine/node/meristems/helpers.h` | Confirm `growth_fraction` generic enough to modulate on any chemical |
| `src/engine/genome.h` | Add `root_ck_growth_floor`; retune `root_cytokinin_production_rate`, `root_tip_auxin_production_rate`, `root_auxin_activation_threshold` |
| `src/engine/node/tissues/root_apical.h` | No API change expected |
| `tests/test_meristem.cpp` | Update "RA auxin production drops when sugar is low" to reflect 0.002 production, or leave the early-skip in place |
| `tests/test_cytokinin_transport.cpp` | Validate CK reaches SA at seedling scale after changes |
| `CLAUDE.md` | Update the hormone section reflecting the new architecture |

All changes are small and localized. No new modules needed.

---

## 10. Success Criteria

This spec is successful when:

1. All existing tests pass.
2. Seedling checkpoint passes with visible auxin and CK gradients.
3. Adolescent checkpoint passes — plant grows to several meters without collapse.
4. Mature checkpoint passes — plant reaches 10m with active CK delivery to the SAM.
5. No long-distance behavior (>20 internodes) depends on auxin reaching any node.
6. The sugar economy self-regulates root-shoot balance without explicit controllers.
