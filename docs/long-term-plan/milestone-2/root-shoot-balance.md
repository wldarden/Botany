# Root ↔ Shoot Balance via Hormonal Feedback

Working document. We just switched xylem from a Jacobi pressure-gradient solve to a demand-driven source-to-sink network, which unclogged water flow to the shoot. The plant immediately went from shoot-stalled (15 SAs frozen at tick 200) to shoot-runaway (521 SAs at tick 420, root:shoot ratio inverted to 0.11:1). Water now reaches shoot meristems correctly, but the shoot has no natural throttle when the root system can't keep up — so it keeps branching.

Real plants have layered negative-feedback loops that prevent this. This doc is the scratchpad for designing an equivalent.

---

## What the current sim has (and doesn't have)

**Already in the sim:**
- Active vs dormant meristems. Only active RAs produce cytokinin (after the option D fix from 2026-04-18). Dormant RAs sit silent.
- Cytokinin transport via xylem (from roots → shoots) is now conservative and works correctly.
- Shoot axillary bud activation gated on `cytokinin > cytokinin_threshold`.
- Root axillary bud activation gated on `auxin > root_auxin_activation_threshold` AND `cytokinin < root_cytokinin_inhibition_threshold`.
- PIN auxin transport carries shoot-derived auxin down to root tips.
- Phloem delivers sugar from leaves → source → sinks (including meristems).

**Missing — the feedback this doc is about:**
- **Sugar-gated cytokinin production at roots.** Currently `cyto_produced = rate × local_auxin` at active RAs, regardless of sugar status. A root tip with zero sugar still produces cytokinin as long as it has auxin. This removes the natural brake real plants have.
- **Any "shoot is bigger than root supports" feedback.** When the shoot out-grows root capacity, nothing in the sim tells the shoot to slow down.
- **Other feedbacks real plants have** (strigolactones, ABA, nitrogen signaling, C:N ratio) — we don't model nitrogen at all, and have no strigolactone or ABA chemical.

---

## The biology — short version

In real plants, the chain that prevents shoot runaway is approximately:

```
shoot too big relative to root
  ↓
shoot sinks consume most leaf-produced sugar locally
  ↓
less sugar reaches root tips via phloem
  ↓
root-tip metabolism slows (ATP ↓, TOR kinase ↓, SnRK1 kinase ↑)
  ↓
CK biosynthesis enzymes (IPT, LOG, CYP735A) less active
  ↓
less CK delivered via xylem to shoot buds
  ↓
fewer shoot axillary activations
  ↓
shoot branching slows, plant re-balances
```

This is the strongest feedback and the one most available to us given our current chemical model.

Real plants also have parallel signals we're not modeling yet:
- **Nitrogen sensing** — roots upregulate CK synthesis when nitrate is abundant; that signal goes up to shoot via xylem CK. No nitrate = low CK = less branching.
- **Strigolactones** — roots make more under phosphate stress; they travel to shoots and directly inhibit axillary bud outgrowth (different gate than CK).
- **ABA** — drought-stressed roots make ABA; it closes stomata in shoots and antagonizes CK at the receptor level.
- **C:N stoichiometry** — whole-plant sensing of carbon inflow vs nitrogen inflow; many hormonal adjustments kick in.

For this doc's scope, we focus on the sugar-CK mechanism because it's the one we can build on top of what's already in the sim. Everything else is a later doc.

---

## Timing — immediate vs averaged

Sugar sensing in real plants happens on layered timescales:

| timescale | what's happening | how it modulates CK |
|---|---|---|
| seconds–minutes | ATP level tracks sugar directly via glycolysis | enzyme substrate availability — instantaneous throttle on CK synthesis |
| minutes–hours | Sugar-sensing kinases (TOR, SnRK1, HXK1) rewire gene expression | IPT enzyme abundance adjusts — changes the *max* rate |
| hours–days | Meristem cellular state (active vs quiescent) | CK production capacity rebuilds or tears down over a day+ |

Our 1 tick = 1 hour means each tick has plenty of time for biochemical equilibration. **Instantaneous sugar level at the tick is a reasonable proxy** for the cell's metabolic state — no explicit history buffer needed.

The biologically natural kinetics are **Michaelis-Menten**, not a step threshold:

```
CK_produced = base_rate × auxin × (sugar / (sugar + K_sugar))
```

- `sugar = 0` → factor = 0 → no CK (substrate-starved).
- `sugar ≈ K_sugar` → factor = 0.5 → half-max production.
- `sugar >> K_sugar` → factor → 1 → saturated (more sugar doesn't help).

`K_sugar` becomes a new genome parameter. Tuning it sets where the feedback kicks in.

### Auxin vs CK — same timescale structure, different buffering

Auxin production has the same three-layer timescale dependence on sugar (biochemistry → transcription → meristem state), so the same instantaneous-sugar approximation applies. **But auxin has extra dynamics that make its effective response to sugar smoother than CK's:**

1. **Active degradation is also sugar-dependent.** Auxin is broken down by GH3 (conjugation) and DAO (oxidation) enzymes, all ATP-powered. When sugar drops, both synthesis and degradation slow — the two partially cancel out in terms of pool-level change. CK has mostly active synthesis and relatively passive degradation (CKX), so its response is steeper.

2. **Conjugate storage pools buffer auxin.** Plants maintain significant IAA-amino-acid and IAA-glucose conjugate pools that hydrolyze back to free auxin under stress. This buffers auxin levels for hours to days after synthesis stops. CK has analogous pools (ribosides, glucosides) but smaller and faster-turning — less buffering.

3. **Active transport (PIN) can mask local changes.** Root-tip auxin is dominated by PIN-delivered inflow from shoot, so local synthesis changes have less effect on local auxin level than they do on CK (where transport is passive-via-xylem and responds to concentration, not actively-pumped).

**This is exactly why the existing SA auxin formula has a 0.1 floor:**

```cpp
sugar_factor = 0.1 + 0.9 * sugar / (sugar + K)
```

The floor encodes auxin's inherent buffering — a meristem at zero sugar still makes ~10% of normal auxin because of stored conjugate release, residual enzyme activity, and PIN delivery from upstream. Real plants don't catastrophically drop auxin from short sugar stress.

**CK production (option A) should probably use a smaller floor or none**, because CK is less buffered:

```cpp
sugar_factor_CK = 0.02 + 0.98 * sugar / (sugar + K)   // small floor
// or
sugar_factor_CK = sugar / (sugar + K)                  // no floor
```

A 2–5% floor captures residual CK released from stored ribosides/glucosides while letting the shoot-branching feedback engage sharply when roots are truly sugar-starved.

**Recommended floor values:**
- RA auxin production: 0.1 (match existing SA formula convention)
- RA cytokinin production: 0.05 or 0.0 (CK responds sharper than auxin to sugar)

Meristem active/dormant transitions are the slower hysteresis layer. Right now our meristems only transition between states via death (under starvation) and bud activation (once). No way to *down-shift* an active meristem back to dormant when conditions deteriorate. That might be worth adding — real root tips can survive weeks of low sugar in a quiescent state.

---

## Architecture: two mechanisms capture three biological timescales

Compressing the three biological timescales (biochem seconds-minutes, transcription minutes-hours, meristem state hours-days) into two sim mechanisms:

### Mechanism 1: Per-tick metabolic gating (short + medium timescales)

A single multiplicative factor applied at every hormone production site. Captures both the biochemical (ATP-substrate) and transcriptional (gene expression) timescales collapsed to the tick resolution.

```cpp
float metabolic_factor(float sugar, float K_sugar, float floor_sugar,
                       float water, float K_water, float floor_water) {
    float sf = floor_sugar + (1 - floor_sugar) * sugar / (sugar + K_sugar);
    float wf = floor_water + (1 - floor_water) * water / (water + K_water);
    return sf * wf;  // multiplicative: both stresses compound
}
```

Floor values encode hormone-specific buffering:
- **Auxin**: floor 0.1 (conjugate pool buffer, PIN inflow from neighbors)
- **Cytokinin**: floor 0.02–0.05 (less buffered, sharper response)

No state variables. Sugar/water snapshot per tick = cell's metabolic report.

### Mechanism 2: Sustained-stress dormancy (long timescale)

Reuse `starvation_ticks` — already tracks consecutive low-sugar ticks on active meristems. Add a quiescence threshold before the death threshold:

```
starvation_ticks < quiescence_threshold    → active, mechanism 1 applies
starvation_ticks ≥ quiescence_threshold    → revert to dormant (reset starvation_ticks, no maintenance)
                                              can re-activate via normal can_activate() if conditions improve
starvation_ticks ≥ starvation_ticks_max    → die (existing)
```

Concrete numbers: `quiescence_threshold = 200` (~8 days), `starvation_ticks_max = 2200` (~90 days) for roots. Shorter for shoots (maybe ~96 ticks = 4 days) since shoot meristems tolerate less stress.

This matches real biology: quiescent root tips survive weeks of stress, reactivate when conditions improve. Our current code kills them prematurely.

### Why only two mechanisms (not three)

The middle timescale (transcription, minutes-hours) doesn't get its own explicit system because:

- Our tick = 1 hour, so per-tick snapshot already aggregates over the transcription response time.
- The *effect* of transcription changes (lower max enzyme capacity under sustained low sugar) is captured in the MM curve itself — the curve shape represents enzyme kinetics at the steady state the cell has reached.
- A third explicit system would need an EMA or similar averaging buffer, adding state and tuning complexity for marginal realism benefit.

If we later find oscillation issues that suggest we need EMA-style smoothing, we can add it — but start simpler.

### Different effective timescales fall out of floor values

The user noted auxin has longer effective averaging than CK due to conjugate pools. In this architecture, that difference is captured entirely by the **floor value in mechanism 1**:
- Zero sugar + active meristem → auxin at 10% (floor), CK at ~2% (smaller floor)
- Sustained starvation → both stop (dormancy threshold, mechanism 2)

Auxin's longer effective buffering shows up as "auxin keeps flowing at floor rate during short stress" vs "CK drops to near zero on short stress". Same mechanism, different parameters. Clean.

An alternative with differentiated dormancy thresholds (CK stops producing sooner than auxin as a meristem dies) is possible but adds complexity for marginal realism benefit. Worth considering only if observed behavior demands it.

## Water modulation of meristem hormone production

A separate gap we're closing alongside the sugar feedback: SA and RA hormone production are currently water-independent. This is biologically wrong — drought represses YUC/TAA1 auxin biosynthesis sharply in all meristems via ABA signaling. The `metabolic_factor` helper handles this as a second MM term with matching floors, included in the converged design below.

At full drought + full sugar starvation: auxin production drops to 1% (0.1²), CK production drops to 0.04% (0.02²). Effectively silent. Matches real severely-stressed meristems.

Leaves already have implicit water-gating on auxin via `growth_fraction` in `expand()`, so they're covered without a code change (though we may want to harmonize later).

---

## The emergent self-balancing property

The most elegant consequence of this architecture: **one core feedback loop, operating at the root tip, handles both runaway directions.** We don't need to separately engineer "shoot is too big" detection and "root is too big" detection. Both converge on the same single signal — root-tip sugar status — and get throttled by the same mechanism.

### Scenario 1: Shoot outgrows roots

```
shoot too large → shoot sinks consume most local sugar
                → less sugar in phloem to roots
                → RA sugar drops
                → metabolic_factor drops
                → both auxin AND CK production drop at RAs
                → less CK delivered to shoot via xylem
                → fewer SA activations (cytokinin_threshold gate)
                → shoot branching slows → plant re-balances
```

### Scenario 2: Roots outgrow shoot

```
root too large → shoot too small to photosynthesize enough
               → total plant sugar production drops
               → roots don't get enough sugar even though THEY want it
               → same RA metabolic_factor drops
               → same hormone shutdown
               → root branching slows (less auxin → fewer root axillary activations)
               → BOTH slow until leaves catch up
```

Both scenarios apply the brake at the same node (root tip sugar status). It doesn't matter which side of the plant caused the imbalance — the plant senses "my root tips don't have enough sugar right now" and scales back hormone signaling. Whichever side over-grew ends up starving the feedback loop at the same point.

### Water as an independent dampener

Drought response is a separate axis in the architecture, running orthogonally through the same metabolic_factor:

```
external drought → water drops across plant
                → SA/RA water drops
                → water term of metabolic_factor drops
                → auxin and CK production both down-regulated everywhere
                → whole plant goes quiet, conserving resources
```

This isn't a runaway-correction mechanism — it's a universal-stress dampener. Plant doesn't try to heroically keep growing during drought; it backs off until conditions improve. Matches real plant behavior.

### Quiescence makes it recoverable

Mechanism 2 (quiescence threshold before death) is what lets stressed plants survive and resume later. Without it, any sustained stress kills meristems and the plant has to rebuild from scratch when conditions recover. With it, the plant down-shifts, waits, and picks up where it left off. This is how real perennials survive winter and drought cycles.

The combined system is:
- **Fast response** (per-tick MM): smooth, continuous dampening as stress rises.
- **Slow response** (quiescence): hard shutdown with graceful recovery when stress persists.
- **Bidirectional** (one feedback loop, either direction of imbalance): no special-casing per scenario.
- **Emergent** (balance is the equilibrium, not an engineered target): the plant finds its sustainable configuration from physics + chemistry.

That's what makes it beautiful — the balance isn't programmed. It's the settling point of a system that follows physical and biological rules.

---

## Converged design — what we're actually going to build

### Code changes

**1. New helper function** (`src/engine/node/meristems/helpers.h` or inlined):

```cpp
// Multiplicative metabolic gating — captures short-timescale biochem +
// transcriptional response to sugar and water availability.  Both stresses
// compound.  Floors encode hormone-specific buffering (conjugate pools, etc.).
inline float metabolic_factor(
    float sugar, float K_sugar, float floor_sugar,
    float water, float K_water, float floor_water)
{
    float sf = floor_sugar + (1.0f - floor_sugar) * sugar / (sugar + K_sugar);
    float wf = floor_water + (1.0f - floor_water) * water / (water + K_water);
    return sf * wf;
}
```

**2. Shoot apical — upgrade existing sugar-only gate to sugar+water:**

In `apical.cpp::produce_auxin`, replace the current `sugar_factor` with:

```cpp
float mf = metabolic_factor(
    chemical(ChemicalID::Sugar), g.auxin_sugar_half_saturation, 0.1f,
    chemical(ChemicalID::Water), g.auxin_water_half_saturation, 0.1f);
float modulated_baseline = base * light_factor * mf;  // replaces sugar_factor
```

**3. Root apical — add auxin and cytokinin gating (currently none):**

In `root_apical.cpp::update_tissue`, within the `active` branch:

```cpp
float mf_auxin = metabolic_factor(
    chemical(ChemicalID::Sugar), g.auxin_sugar_half_saturation, 0.1f,
    chemical(ChemicalID::Water), g.auxin_water_half_saturation, 0.1f);
float produced_auxin = g.root_tip_auxin_production_rate * mf_auxin;
chemical(ChemicalID::Auxin) += produced_auxin;
tick_auxin_produced += produced_auxin;

float mf_cyto = metabolic_factor(
    chemical(ChemicalID::Sugar), g.cytokinin_sugar_half_saturation, 0.05f,
    chemical(ChemicalID::Water), g.cytokinin_water_half_saturation, 0.05f);
float produced_cyto = g.root_cytokinin_production_rate * chemical(ChemicalID::Auxin) * mf_cyto;
chemical(ChemicalID::Cytokinin) += produced_cyto;
tick_cytokinin_produced += produced_cyto;
```

**4. Quiescence transition** — add to every meristem's tick logic:

```cpp
// If an active meristem has been starved too long, go dormant instead
// of dying.  Preserves the bud for future conditions.  Applied to both
// ApicalNode and RootApicalNode.
if (active && starvation_ticks >= g.quiescence_threshold) {
    active = false;
    starvation_ticks = 0;
    // No maintenance cost while dormant (existing behavior).
    // Can re-activate via normal can_activate() if sugar/auxin/cyto recover.
}
```

### New genome parameters

| parameter | purpose | suggested default |
|---|---|---|
| `auxin_water_half_saturation` | K for water term in auxin MM | ~0.1 ml (calibrate) |
| `cytokinin_sugar_half_saturation` | K for sugar term in CK MM | ~0.05 g (calibrate) |
| `cytokinin_water_half_saturation` | K for water term in CK MM | ~0.1 ml (calibrate) |
| `quiescence_threshold` | starvation_ticks at which meristem reverts to dormant | ~200 ticks for roots, ~96 for shoots |

All go in `genome.h`; defaults in `default_genome()`; exposed to evolution if we want per-species variation (likely yes for the K values, maybe no for the thresholds).

### Floors (not evolvable)

These represent fundamental biochemistry (buffering via conjugate pools) and don't vary much between species. Set as constants in code:
- Auxin floor: 0.1
- Cytokinin floor: 0.05

### Not in this round

- **Strigolactones** — separate shoot-branching suppressor. Would need new chemical, transport, receptors. Deferred.
- **ABA** — drought signal distinct from metabolic_factor. Needed only if we add meaningful soil_moisture < 1 simulation. Deferred.
- **Nitrogen signaling** — we don't model nitrogen at all. Large addition. Deferred.
- **Differentiated dormancy thresholds** (CK stops before auxin during slow death) — complexity for marginal gain. Only add if observed.
- **EMA-based smoothing of sugar signal** — add only if oscillation is observed.
- **`wcap_meristem` reduction** — an orthogonal tuning issue that also affects water distribution. Consider as a separate change if SAs still over-compete for water after this.

---

## Open questions / things to watch when we run this

### Do existing active SAs need to throttle too, not just new activations?

The CK feedback gates *new* shoot activations via `cytokinin > cytokinin_threshold`. But once active, an SA keeps spawning internodes regardless of current CK level. During a transition period after CK drops, existing active SAs continue running.

This is probably fine — the quiescence mechanism covers it. If the plant is genuinely resource-starved, existing SAs will accumulate `starvation_ticks` and eventually flip to dormant. So the long-timescale feedback handles it even though the short-term feedback doesn't directly.

Watch for: does observed behavior match? Does shoot branching smoothly taper off, or does it crash abruptly when existing SAs all starve at once?

### Oscillation risk

The MM gating is instantaneous. In a stable plant it should settle smoothly. Possible edge case:

- Shoot grows → consumes more sugar → phloem delivery to roots drops → CK drops → shoot growth slows → roots recover → CK rises → shoot grows again...

That oscillation cycle would take hundreds of ticks, probably realistic (real plants do have multi-day growth rhythms). If it goes unstable we can add EMA smoothing on the sugar signal before feeding it to MM. Decide empirically.

### `wcap_meristem` — is 1.0 ml too big?

Orthogonal to the feedback design but relevant because it shapes how many SAs can compete for water. Real meristems are tiny (0.1–1 mm diameter → sub-microliter water volume). Our 1 ml cap may be 100× too generous — causing 521 SAs to demand 365 ml of water per tick against a ~15 ml root supply.

Even with CK feedback working, if `wcap_meristem` stays at 1.0, water competition may still squeeze leaves to wfrac 0.001. Consider as an independent tuning knob.

### Does the root side activation gate need explicit throttling?

Primary RA is always active (starts so by `can_activate`) and produces auxin/CK (now metabolism-gated). Lateral RAs activate when PIN auxin reaches them and cytokinin is low. With the metabolic feedback working, this gate should naturally throttle root branching during stress. But watch for edge cases where PIN delivery overwhelms the cytokinin inhibition.

### When do we need strigolactones, ABA, or nitrogen?

Our converged design handles the core shoot↔root balance cleanly via sugar+water gating. The other real-plant hormones are worth considering later for:
- **Strigolactones**: if we want plants to respond differently to phosphate stress vs general nutrient stress.
- **ABA**: if we enable soil_moisture < 1 scenarios and want correct drought response (stomatal closure signal).
- **Nitrogen**: if we want plants to respond to soil nitrogen variation. Requires adding nitrogen as a tracked nutrient — major addition.

None are needed for the current balance problem. Flagged here so we don't forget when expanding the sim's environmental richness.

### Genome param vs WorldParam for the new K values

Our lean is genome, matching the existing `auxin_sugar_half_saturation` pattern. Real species differ in meristem sensitivity thresholds. Confirmed during implementation.

---

## Implementation sequence

1. **Add the `metabolic_factor` helper** in meristem helpers. Pure function, no state.
2. **Apply to SA auxin production** — upgrade existing sugar-only gate to sugar+water using the new helper. Should produce a small visible change in drought scenarios; no change at default soil moisture.
3. **Apply to RA auxin production** — new gate, inside the `active` branch. Previously unmodulated.
4. **Apply to RA cytokinin production** — new gate, multiplies with the existing `× chemical(Auxin)` factor. This is the main shoot-runaway brake.
5. **Add quiescence transition** in both `ApicalNode` and `RootApicalNode` tick logic — revert to dormant instead of dying when `starvation_ticks ≥ quiescence_threshold`.
6. **Add new genome parameters** with defaults.
7. **Run 1000+ ticks.** Check:
   - Conservation still clean (should be — no algorithm change).
   - Debug log: RA auxin_produced and cytokinin_produced respond to local sugar/water as the MM curve predicts.
   - Trajectory: does root:shoot ratio stay in a reasonable range (0.5–5.0)? Or does it still run away?
   - Over long runs: do we see the plant enter/exit quiescence under stress and recover cleanly?
8. **Tune K values** against observed behavior. Likely 1–2 iterations.

Each step is small and measurable against the existing debug logs (per-node `auxin_produced`, `cytokinin_produced`, `sugar`, `starvation_ticks` are all already logged).

---

## Reference — values observed in the recent runs for calibration

From the 420-tick run after the demand-driven xylem switch:

- Plant nodes: 1732 at tick 420 (3 → 1732 = exponential early phase)
- 521 SAs (120 at tick 200, 521 at tick 420 — doubling every ~100 ticks)
- 86 RAs, only 8 active
- Seed sugar: 48 → 0.24 g over 420 ticks (essentially depleted)
- Production: 0.079 g/tick rising, meristem maintenance total ~0.05 g/tick
- Plant CK total: ~? (need to extract from log)
- Plant auxin total: ~? (need to extract from log)
- RA sugar: ~? (the key missing number — does phloem reach root tips?)
