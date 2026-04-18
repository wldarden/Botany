# Storage System: Starch / Sugar Reserve

## Section 1: The Goal

Plants convert surplus sugar to starch for storage. Starch is insoluble — it doesn't dissolve,
doesn't transport, doesn't affect concentration gradients. It just sits where it was made. When
the plant needs energy later, starch is mobilized back to sugar. This is the battery system that
enables seasonal survival, seedling bootstrap, and buffered growth.

Two chemicals: Sugar (mobile, phloem-transported) and Starch (immobile, stored in place).
Conversion is local, bidirectional, enzyme-driven:

- High local sugar → starch (banking surplus)
- Low local sugar + GA signal → sugar (spending reserves)

## Section 2: How It Works in Real Plants

### Storage triggers (sugar → starch)

- Local sugar concentration exceeds usage — direct enzymatic response, no hormone needed
- ABA promotes storage (future — not modeled yet)

### Mobilization triggers (starch → sugar)

- Low local sugar (demand exceeds supply) — basic feedback
- Gibberellin (GA) — the primary hormonal mobilization trigger. Seeds germinate this way: the
  embryo produces GA, GA activates amylase in the endosperm, starch converts to sugar. We already
  have GA in the sim.
- Cytokinin — indirectly, by promoting growth which creates demand

### The seasonal cycle

Spring mobilization (starch → sugar funds bud break) → summer storage (surplus → starch) →
autumn maximum reserves → winter slow mobilization for maintenance.

## Section 3: What This Fixes In Our Sim

- **The "44g frozen in stems" problem.** That sugar would be starch, mobilizable via GA when
  leaves need it. Currently it just sits at cap and blocks transport.
- **The seedling bootstrap.** The seed IS a starch reserve. GA from the embryo (shoot apical)
  triggers mobilization during germination. We can replace the current `seed_sugar` hack with a
  proper starch reserve that depletes naturally.
- **The stem self-photosynthesis bootstrap.** Young green stems produce sugar; excess converts to
  starch in adjacent mature stems and is available later when leaves take over.
- **Buffered growth.** Nodes don't starve from single-tick supply fluctuations because they carry
  local reserves. The vascular system only needs to cover the delta, not the full maintenance
  cost every tick.

## Section 4: Implementation Sketch

### New chemical

Add `Starch` as a `ChemicalID`. Properties:
- No transport (never enters `transport_chemicals`)
- No diffusion
- No decay
- Purely local — accumulates and depletes in place

### Per-node conversion logic (inside `update_tissue`)

```
// Storage: bank surplus
if sugar > storage_threshold × sugar_cap:
    convert = min(surplus × starch_synthesis_rate, starch_headroom)
    sugar -= convert
    starch += convert

// Mobilization: spend reserves
if sugar < mobilization_threshold × sugar_cap:
    ga_boost = 1.0 + ga_level × starch_mobilization_ga_boost
    convert = min(starch × starch_mobilization_rate × ga_boost, starch_available)
    starch -= convert
    sugar += convert
```

### New genome parameters

| Parameter | Suggested default | Notes |
|---|---|---|
| `starch_synthesis_rate` | 0.1 / tick | Fraction of surplus converted each tick |
| `starch_mobilization_rate` | 0.05 / tick | Fraction of starch mobilized per tick when low |
| `starch_mobilization_ga_boost` | 3.0 | GA multiplier on mobilization rate |
| `storage_threshold` | 0.7 | Sugar fraction of cap that triggers synthesis |
| `mobilization_threshold` | 0.3 | Sugar fraction of cap that triggers mobilization |
| `starch_cap_density` | same as `sugar_storage_density_wood` | g starch / dm³ tissue |

### Starch cap

Use the same volume formula as sugar cap. Big stems store proportionally more — this is correct
(wood parenchyma is the primary starch store in real trees).

### Seed initialization

Replace `seed_sugar` with `seed_starch`. On tick 1, the seed's shoot apical produces GA; GA
triggers mobilization; starch flows to sugar at `starch_mobilization_rate × ga_boost` per tick.
This makes germination a process, not an instantaneous grant.

Migration: the existing `seed_sugar` default genome value can stay for backward compatibility
until this is implemented — it just won't be "realistic."

### GA role expansion

GA currently affects elongation rate and max internode length. Add a third effect: starch
mobilization multiplier (already parameterized above as `starch_mobilization_ga_boost`). No
changes to GA production or transport — just an additional read site in the mobilization check.

## Section 5: Open Questions

1. **Starch cap formula.** Same density as sugar cap (same parenchyma cells store both), or a
   separate `starch_storage_density` genome param? Separate param gives evolution more freedom
   but adds complexity.

2. **Enzymatic vs. rate-limited conversion.** The sketch above uses a simple rate limit. A
   more realistic model would require an "amylase activity" intermediate that GA upregulates over
   several ticks. Is the added fidelity worth the complexity at this stage?

3. **ABA integration.** ABA promotes storage and suppresses growth. When we add ABA (future
   milestone), it should multiply `starch_synthesis_rate` and suppress the mobilization trigger.
   The conversion logic above should be written to accept an ABA term as a zero-default
   parameter so the hook is obvious.

4. **Mass and structural integrity.** Does starch contribute to node mass (affecting stress/droop
   calculations)? Starch granules are dense — a fully loaded parenchyma cell is heavier than an
   empty one. Probably yes, but it's a second-order effect.

5. **Meristem starch.** Should meristems carry starch? Real shoot apical meristems are largely
   starchless — they rely on phloem import, not local reserves. We could exclude meristems from
   synthesis (only allow mobilization) to match this biology.
