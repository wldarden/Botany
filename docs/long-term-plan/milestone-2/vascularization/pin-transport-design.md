# PIN Auxin Transport — Design Document

## Problem

Auxin in the sim is transported by local diffusion: `auxin_diffusion_rate = 0.05` (5% of the
concentration gradient per tick) with `auxin_decay_rate = 0.12` (12% destroyed per tick).

At steady state in a single chain, each node holds a fraction of the node above it:

```
ratio_per_hop ≈ diffusion_rate / (diffusion_rate + decay_rate)
             = 0.05 / (0.05 + 0.12) ≈ 0.29
```

At node 5 below the apex: `0.29^5 ≈ 0.002` of the apex level — essentially zero. The practical
signal reach is **2–3 nodes**. In a mature plant with 30–50 nodes from shoot tip to root tip,
auxin never crosses the seed junction and root tips receive no shoot-derived auxin signal.

This breaks several dependent systems:

- **Apical dominance** requires a gradient across the full shoot, not just the top 2 nodes.
  Lateral buds 4+ nodes below the apex activate immediately, indistinguishable from buds
  directly below the apex.
- **Root cytokinin production** is gated on `local_auxin × root_cytokinin_production_rate`. With
  zero shoot auxin reaching roots, cytokinin production decouples from shoot activity entirely.
- **Canalization** records which connections carry auxin flux. If auxin can't move past the first
  few nodes, only those connections ever build `auxin_flow_bias` and thicken. The main trunk
  never canalizes — and never earns the larger radius — because auxin never travels through it.

Real plants transport auxin across the entire plant axis in hours via PIN efflux carriers — a
fundamentally different mechanism from diffusion, operating at a different scale and requiring
explicit modeling.

---

## Biology Background

**Polar auxin transport (PAT)** is an active, energy-consuming, direction-specific process
mediated by PIN-FORMED (PIN) efflux carrier proteins embedded in the plasma membrane. Key
properties that distinguish it from diffusion:

1. **Asymmetric localization.** PIN proteins are sorted to specific cell faces. PIN1 in shoot
   tissue is basal-face dominant — auxin exits cells preferentially toward the root. PIN2 in
   root epidermis is apical-face dominant — auxin moves toward the root tip. Each cell acts
   as a one-way pump in a direction set by its PIN distribution, not by concentration gradients.

2. **Active, ATP-dependent.** Auxin is pumped against its concentration gradient where needed.
   PAT is blocked by metabolic inhibitors (TIBA, NPA) and by ATPase inhibitors. This means
   the transport requires sugar — starved tissue loses PIN polarity and auxin flux collapses.

3. **Canalization feedback.** Where auxin flux is high, PIN proteins redistribute within cells
   to reinforce that flux direction (Sachs 1969). Where flux is low, PIN polarity randomizes
   and disappears. This is the cellular basis of `auxin_flow_bias` in the sim — the transient
   canalization memory that records sustained directional PIN-mediated flux on each connection.

4. **Speed.** PAT moves auxin at ~5–10 mm/hour in real stems. At the sim's dm/hour scale that
   is ~0.05–0.1 dm/hr — fast enough to traverse a whole plant in 10–20 ticks, far faster than
   diffusion at a 5% gradient per tick.

5. **Tissue-specific direction.** Stems are basipetal (shoot tip → root). Root vasculature is
   acropetal (root base → root tip). The sim needs both.

---

## Design: Dedicated PIN Transport Pass

PIN transport is a new global pass, separate from vascular transport and from local diffusion.
It runs once per tick and moves auxin directionally along the plant axis, with per-connection
conductance scaling with `radius` — thicker stems carry more PIN-mediated auxin flow.

### System Architecture

```
Three separate transport systems, operating at different scales:

1. vascular_transport()  — sugar, water, cytokinin bulk flow (phloem + xylem)
2. pin_transport()       — auxin directed long-range transport (PIN proteins)  ← NEW
3. DFS tree walk         — local diffusion for all chemicals (cell-to-cell)
```

PIN and diffusion **coexist** for auxin. They are complementary:

- **PIN transport** handles long-distance directional flow along the main axis — the auxin
  highway from shoot tip to root. This is the signal that spans the whole plant.
- **Local diffusion** handles short-range lateral spreading — gradients within a few nodes,
  lateral inhibition between competing buds, last-mile delivery from the vascular axis to
  surrounding parenchyma. Diffusion parameters for auxin **do not change**.

### Tick Order

```
1. vascular_transport(plant, genome, world)    — existing (sugar, water, cytokinin)
2. pin_transport(plant, genome)               — NEW (auxin only)
3. tick_recursive(seed)                       — existing DFS walk
     └── per node: update_tissue → diffusion → decay
```

PIN runs before the DFS walk so that auxin is in the correct position before `update_tissue`
reads it this tick. PIN runs after vascular for ordering consistency, though it has no
dependency on vascular output (PIN carries no sugar cost — see Genome Parameters).

### Direction Model

| Node type | PIN direction | What happens |
|-----------|---------------|--------------|
| APICAL, STEM, LEAF | Basipetal | node pumps auxin toward its parent |
| Seed | Junction | collects from shoot children, distributes to root children |
| ROOT, ROOT_APICAL | Acropetal | parent distributes auxin to root children |

The direction flip at the seed mirrors the existing `auxin_bias` sign inversion for ROOT and
ROOT_APICAL children in `transport_with_children`. The same biological logic applies: the seed
is the transit junction. Auxin flows down the shoot and then continues down the root — it does
not reverse direction at the root tip and flow back up.

### Algorithm

**Post-order collection (shoot side):**

Walk shoot nodes from leaves toward seed. For each shoot node (STEM, APICAL, LEAF):

```cpp
// PIN capacity scales with radius — thicker stems have more developed PIN machinery.
// radius_factor = 1.0 at initial_radius, grows with thickening.
float radius_factor = radius / std::max(g.initial_radius, 1e-6f);
float rate = std::clamp(
    g.pin_base_rate + radius_factor * (g.pin_max_rate - g.pin_base_rate),
    0.0f, g.pin_max_rate
);

float moved = chemical(ChemicalID::Auxin) * rate;
chemical(ChemicalID::Auxin) -= moved;
parent->transport_received[ChemicalID::Auxin] += moved;

// Record flux on parent for canalization (feeds auxin_flow_bias / update_canalization)
parent->last_auxin_flux[this] += moved;
```

The existing `transport_received` buffer is used for PIN exactly as in `transport_with_children`
— anti-teleportation is already in place, and the flushing logic in `Node::tick()` applies to
both. PIN records flux into `last_auxin_flux` on the parent, the same map that `update_canalization`
already reads to update `auxin_flow_bias` each tick.

**Seed junction:**

The seed is always at full PIN conductance (oldest, most-canalized node). After flushing its
`transport_received` buffer from shoot-side children:

```cpp
// Weight distribution to root children by radius — thicker roots have more PIN capacity.
float total_root_radius = 0.0f;
for (Node* rc : root_children) total_root_radius += rc->radius;
for (Node* root_child : root_children) {
    float share = (total_root_radius > 1e-8f)
        ? root_child->radius / total_root_radius
        : 1.0f / root_children.size();
    float to_send = chemical(ChemicalID::Auxin) * share * g.pin_base_rate;
    root_child->transport_received[ChemicalID::Auxin] += to_send;
    chemical(ChemicalID::Auxin) -= to_send;
}
```

Root children with larger radius receive a larger share of the auxin transiting through the
seed — the most-developed root branch gets the strongest shoot-derived auxin signal. Radius is
the record of past canalization, so this naturally favors established root paths.

**Pre-order distribution (root side):**

Walk root nodes from seed toward root tips. Each root node flushes its `transport_received`
buffer (already filled by its parent's pass) and then distributes to its own root children using
the same conductance-weighted split. This mirrors how the xylem vascular pass distributes water
pre-order, but with PIN-rate conductance instead of pipe cross-section.

### Canalization Coupling, Flux Recording, and the Structural Memory

PIN is the **primary driver of `auxin_flow_bias`** — and `auxin_flow_bias` IS tracked across
ticks. Each tick it smoothly lerps toward the current PIN saturation: the fraction of available
transport capacity currently being used on each connection.

`pin_transport()` records per-connection auxin flux into the parent node's `last_auxin_flux`
map. After the PIN pass completes, `update_canalization` reads that map and computes:

```
current_saturation = auxin_flux / (radius² × pin_capacity_per_area)
auxin_flow_bias    = lerp(previous_auxin_flow_bias, current_saturation, smoothing_rate)
```

Where:
- `auxin_flux` — auxin moved through this connection this tick (from `last_auxin_flux`)
- `radius²` — cross-sectional area of the connection (proportional to available PIN protein slots)
- `pin_capacity_per_area` — genome param: max auxin flow per unit cross-section at full saturation
- `smoothing_rate` — genome param (~0.1): how fast `auxin_flow_bias` chases the current saturation
- `lerp(a, b, t)` — standard linear interpolation: `a + t × (b − a)`

`current_saturation` is a **0-to-1 value clamped at 1.0** representing instantaneous PIN saturation.
`auxin_flow_bias` is the smoothed, tick-persistent version of that saturation — the value
the cambium and diffusion weighting actually read.

#### Biological Basis: Auxin Flux Density

This models **auxin flux density** — auxin per unit cross-sectional area per unit time. Real
cambium activation responds to flux density, not absolute flux. A thin strand carrying 1 unit
of auxin is more highly activated than a thick trunk carrying the same 1 unit — the thin strand
is working at full capacity while the trunk is barely using its transport area. High flux
density = high PIN polarization = active cambium. Low flux density = dormant cambium.

`auxin_flow_bias` has smoothed memory — it tracks the running saturation via lerp — but no
ratchet. When PIN flux stops, `current_saturation` drops to zero and `auxin_flow_bias` lerps
toward zero over the next ~10 ticks (`smoothing_rate = 0.1`). It decays naturally without a
separate decay parameter. The **permanent** memory is the radius — the stem thickened from past
high-saturation ticks, and that thickness never reverses.

#### Scenario Analysis

With `pin_capacity_per_area = 100` AU/(dm²·tick), the formula produces the following behavior
across typical plant stages. The r² column is the connection's cross-sectional area in dm²:

| Scenario | Flux (AU/tick) | r² (dm²) | Saturation | Cambium response |
|----------|---------------|----------|------------|-----------------|
| Young thin stem, modest auxin | 1.0 | 0.01 | ~1.0 | Fully saturated → fast thickening |
| Old thick trunk, same auxin | 1.0 | 1.0 | ~0.01 | Barely using capacity → cambium barely active |
| Old thick trunk, whole canopy feeding | 50.0 | 1.0 | ~0.5 | Half saturated → moderate thickening |
| Thin side branch, weak auxin | 0.1 | 0.01 | ~0.1 | Low saturation → slow thickening |

**Calibration target:** with `pin_capacity_per_area = 100`, a thin active stem (r² ≈ 0.005 dm²)
carrying ~0.3 AU/tick from a single apical hits saturation ≈ 0.6 — comfortably in the 0.5–0.8
working range. Adjust `pin_capacity_per_area` upward to reduce overall cambium activity, or
downward to increase it. It is a global sensitivity knob, not a per-branch parameter.

#### Self-Limiting Property

A trunk cannot accumulate infinite flow bias because larger cross-section dilutes saturation:

```
saturation = flux / (r² × capacity)
```

As the trunk thickens (r² grows), the same flux produces lower saturation → less cambium
activity → slower thickening. Growth decelerates as the stem becomes adequate for its auxin
load. This is why real trunks reach a steady-state diameter under constant photosynthetic load.

#### Junction Competition

Two sibling branches share their parent's auxin flow. If sibling A is thinner than sibling B
and receives the same absolute flux, sibling A has higher saturation → higher `auxin_flow_bias`
→ thickens faster. No special junction logic is needed. The flux-density formula handles
sibling balancing automatically: the thinner sibling is more PIN-saturated by the same flow
and responds more aggressively.

**The full feedback chain:**

```
PIN moves auxin through a connection
  → PIN records flux in parent.last_auxin_flux[child]
  → auxin_flow_bias = flux / (radius² × pin_capacity_per_area)   [PIN saturation, 0→1]
  → auxin_flow_bias drives cambium (thicken()):
      delta_radius = cambium_responsiveness × auxin_flow_bias × sugar_gf
  → radius grows → r² increases → same flux produces lower future saturation (self-limiting)
  → larger radius → more vascular conductance (π × r²)
  → more sugar and water delivered → tissue grows faster → more auxin produced
  → more PIN flux → saturation stays high → cambium stays active
```

**The structural memory is radius, not a software variable.**

When PIN flux stops, `auxin_flow_bias = 0` → cambium halts → no more thickening. But the
radius gained is permanent — wood does not un-grow. A thick stem retains its pipe capacity
(`π × r²`) even during dormancy or after the apex is shaded. The "memory" of past canalization
is the wood itself, readable from `radius` alone.

`structural_flow_bias` is deleted. Its role as "permanent structural memory" is served by
radius. Its role as "cambium driver" is served by `auxin_flow_bias`. No separate ratchet
variable is needed.

#### Temporal Smoothing

The definitive update rule uses exponential smoothing (lerp) rather than raw per-tick saturation:

```
current_saturation = auxin_flux / (radius² × pin_capacity_per_area)   // [0, 1]
auxin_flow_bias    = lerp(previous_auxin_flow_bias, current_saturation, smoothing_rate)
```

With `smoothing_rate = 0.1`, `auxin_flow_bias` moves 10% toward `current_saturation` each
tick — fast enough to respond in ~20–30 ticks but immune to single-tick flux spikes or
transient gaps (branch sway, momentary shadow, a tick with no growing leaves). Without
smoothing, a single zero-flux tick would immediately drop `auxin_flow_bias` to zero and
halt the cambium — biologically unrealistic and numerically noisy.

`auxin_flow_bias` naturally decays toward zero when PIN flux stops: with no incoming auxin,
`current_saturation = 0`, and `lerp(bias, 0, 0.1)` reduces the bias by 10% per tick.
After ~20 ticks of dormancy it is effectively zero. No separate decay parameter is needed.

### Canalization Behavior After Removing structural_flow_bias

**We are currently double-biasing.** Both `radius` and `structural_flow_bias` weight flow
simultaneously across all three transport systems. The result is that canalization is amplified
beyond what either factor alone would produce — established paths get a bonus from their radius
AND another bonus from their accumulated structural bias. Removing `structural_flow_bias` brings
weighting to single-strength via radius alone. This is a simplification, not a weakening: the
feedback loop is intact, just not double-counted.

| System | Current weighting | After removal | Canalization preserved? |
|--------|------------------|---------------|------------------------|
| Vascular distribution | `r² × (1 + canalization_weight × structural_flow_bias)` | `r²` (pipe_capacity) | **Yes** — thicker stems still get more flow |
| Diffusion sibling weighting | `radius_factor × (1 + canalization_weight × (auxin_flow_bias + structural_flow_bias))` | `radius_factor × (1 + canalization_weight × auxin_flow_bias)` | **Yes** — radius + transient bias both still weight |
| Thickening | `cambium_responsiveness × structural_flow_bias × sugar_gf` | `cambium_responsiveness × auxin_flow_bias × sugar_gf` | **Intentionally changed** — current signal drives cambium, not frozen history |
| PIN transport | `pin_base + structural_flow_bias × range` | `pin_base + radius_factor × range` | **Yes** — thicker stems carry more auxin |

**The full feedback loop is intact.** Removing `structural_flow_bias` does not break the
self-reinforcing hierarchy — it just removes the extra amplifier:

```
auxin flow (PIN)
  → auxin_flow_bias builds (transient)
  → cambium activates → radius grows
  → larger radius → more vascular conductance (π × r²)
  → more sugar + water → faster growth → more auxin
  → more PIN flux → auxin_flow_bias stays high → more thickening
```

The loop closes cleanly through radius. Every step still favors established paths over new
ones. The difference is that thickening now stops when auxin stops — which is correct. Old
wood keeps its conductance advantage permanently; new wood only forms when auxin is currently
flowing.

**If radius-only weighting feels too weak:** increase `cambium_responsiveness`. This makes
thickening respond more aggressively to the transient `auxin_flow_bias` signal, building
radius faster from the same flux history. The double-weighting was masking what should have
been a calibration exercise on `cambium_responsiveness`. The right fix for "canalization feels
too weak" is a larger `cambium_responsiveness`, not a secondary permanent tracking variable.

**PIN capacity scales with radius:**

```cpp
float radius_factor = radius / std::max(g.initial_radius, 1e-6f);
float rate = clamp(pin_base_rate + radius_factor * (pin_max_rate - pin_base_rate), 0, pin_max_rate);
```

Thicker stems have more developed PIN efflux machinery — more cells in the vascular cylinder,
more PIN1 protein, higher auxin throughput. A newly-spawned thin internode starts near
`pin_base_rate`; a mature trunk at 5× initial radius approaches `pin_max_rate`. The
bootstrapping is self-consistent: new internodes carry little auxin → weak `auxin_flow_bias` →
slow cambium → slow thickening → slow PIN rate increase. As the meristem lays down longer,
healthier internodes (larger initial radius), they enter with a higher starting `radius_factor`
and bootstrap faster. No separate "initial bias stamp" parameter is needed — initial radius
at creation is already set by the meristem and already determines initial PIN capacity.

### New Genome Parameters

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `pin_base_rate` | 0.30 | fraction/tick | Fraction of auxin pumped per tick at initial radius |
| `pin_max_rate` | 0.80 | fraction/tick | Upper clamp after radius scaling |
| `pin_capacity_per_area` | 100.0 | AU/(dm²·tick) | Max auxin throughput per unit cross-section; divides flux to get saturation |
| `smoothing_rate` | 0.1 | dimensionless | Lerp rate for `auxin_flow_bias` toward current saturation (~20–30 tick response time) |

PIN transport carries no sugar cost. The sugar economy is still being tuned; adding a PIN
drain would make debugging harder. Sugar cost can be added later as a cost-of-coordination
tradeoff once the economy is stable.

With `pin_base_rate = 0.30` and `auxin_decay_rate = 0.12`, the effective reach in a chain:

```
ratio_per_hop ≈ pin_rate / (pin_rate + decay_rate)
             = 0.30 / (0.30 + 0.12) ≈ 0.71

At node 10:  0.71^10 ≈ 0.035   — meaningful signal
At node 20:  0.71^20 ≈ 0.001   — very weak but present
```

On an established main axis at 5× initial radius (`radius_factor = 5.0`):

```
pin_rate = 0.30 + 5.0 × (0.80 - 0.30) = 0.30 + 2.5 → clamped to 0.80
ratio_per_hop ≈ 0.80 / 0.92 ≈ 0.87

At node 10:  0.87^10 ≈ 0.25    — strong signal even 10 nodes down
At node 20:  0.87^20 ≈ 0.063   — still meaningful at 20 nodes
```

Main axis thickening extends auxin's effective reach from 2–3 nodes (current diffusion alone)
to 15–20+ nodes, while thin lateral branches near `pin_base_rate` stay at shorter range. This
creates the gradient needed for correct apical dominance.

### No Changes to Existing Diffusion Parameters

Auxin diffusion (`diffusion_rate = 0.05`, `decay_rate = 0.12`) is unchanged. Diffusion handles
the lateral and local roles it always has — spreading auxin from the main PIN axis into
surrounding parenchyma, creating local gradients between competing buds, and delivering to leaf
nodes and meristems that sit off the main canalized path. PIN and diffusion are additive
contributions at each node.

---

## Cytokinin: No PIN Needed

Cytokinin transport is already correctly modeled. No changes are required.

In real plants, the dominant long-range cytokinin transport mechanism is passive mass flow in
the xylem. *trans*-Zeatin riboside — the form produced in root tips — is loaded into xylem
vessels and carried upward with the transpiration stream. There is no PIN-like active directional
carrier for cytokinin. Some cytokinin transporters exist (AZA/PUP family) but they operate
locally and are not the primary long-range mechanism.

The sim's model matches this biology exactly:

- Root apicals produce cytokinin gated by local auxin
- `vascular_transport()` carries cytokinin in the xylem stream (`is_vascular_chemical(Cytokinin)
  == true`)
- Last-mile local diffusion delivers it from vascularized internodes to nearby meristems and
  leaves

The cytokinin system requires no architectural change. If the signal feels weak after PIN is
implemented and the whole-plant auxin loop is confirmed working, the right lever is tuning
`cytokinin_decay_rate` downward — not adding a dedicated transport mechanism.

---

## Expected Behavior After Implementation

**`botany_realtime --color auxin`:**
- Bright auxin at shoot tips
- Smooth gradient declining over 15–20 nodes down the main axis
- Non-zero (dim) auxin visible in the upper root internodes
- Lateral branches showing lower auxin than main axis at the same height
- Near-zero auxin in newly-spawned thin internodes until PIN flux raises `auxin_flow_bias` and thickening increases `radius_factor`

**Apical dominance:**
- Buds immediately below the apex remain suppressed (high auxin)
- Buds 5+ nodes below activate more freely (gradient has dropped)
- After the apex is removed (or goes dormant), auxin level drops along the full axis within
  ~10–20 ticks, and lateral buds activate progressively from the top down

**Root cytokinin feedback:**
- Shoot auxin → PIN carries it basipetally → arrives at root internodes → reaches root apicals
  via local diffusion → `root_cytokinin_production_rate × local_auxin` → cytokinin produced →
  xylem carries it to shoot. The full feedback loop is now connected across the whole plant.

**Canalization divergence:**
- Main shoot axis: continuous auxin flux since germination → `auxin_flow_bias` sustained →
  cambium active → radius grows → larger radius → higher PIN rate → more auxin. Self-reinforcing.
- Lateral branches that activated late: start thin → low PIN rate → less auxin → slow thickening.
  Competitive disadvantage that matches biological reality. Their radius is the permanent record
  of that disadvantage — even if they later receive more auxin, they start from a smaller pipe.

---

## Implementation Plan

**Step 1: Files**

Create `src/engine/pin_transport.h` and `src/engine/pin_transport.cpp`. Declaration:
```cpp
void pin_transport(Plant& plant, const Genome& g);
```

**Step 2: Genome parameters**

In `genome.h`: add `pin_base_rate`, `pin_max_rate`, `pin_capacity_per_area`, `smoothing_rate`.
Remove `structural_growth_rate`, `structural_threshold`, `structural_max`,
`vascular_conductance_threshold` (all `structural_flow_bias`-related).
Add `vascular_radius_threshold` (replaces the bias-based admission gate).
In `default_genome()`: set new params to table values above.
In `genome_bridge.cpp`: update `build_genome_template()` accordingly.
In `CLAUDE.md`: update Tuning Parameters section.

**Step 3: Wire into plant.cpp**

In `Plant::tick()`, add `pin_transport(*this, genome_)` immediately after `vascular_transport()`
and before `tick_tree()`. Add `#include "engine/pin_transport.h"`.

**Step 4: Regression test**

Add `tests/test_pin_transport.cpp`:
- Grow a straight-chain plant for 200 ticks
- Assert auxin at seed junction is non-zero (crossed the junction)
- Assert auxin at the root apical is non-zero (reached the root tip)
- Assert auxin gradient direction: shoot apex > seed > root apical
- Assert `auxin_flow_bias` on the main axis connections is non-zero (PIN is recording flux)
- Assert main axis radius > lateral branch radius at comparable age (thickening is happening)
- Confirm `Node` has no `structural_flow_bias` member (grep or static assert)

**Step 5: Update CLAUDE.md**

Add PIN transport to the Chemical Transport Model table:

```
| Auxin | PIN transport (long-range) + local diffusion | Basipetal in stems, acropetal in roots |
```

Add to Tick Control Flow section. Add genome parameter entries for the two new params.
Note that PIN records auxin flux into `last_auxin_flux` and that `update_canalization`
runs after diffusion (same tick), reading that map to update `auxin_flow_bias`.
Note that `structural_flow_bias` is deleted; radius is the structural memory.

---

## Open Questions (Not Blocking)

**Root tip auxin maximum.** Real root tips have a local auxin maximum created by PIN2
(epidermis) recycling auxin from the columella toward the elongation zone — the gradient that
drives gravitropism. The sim uses direct angle adjustment via `meristem_gravitropism_rate`.
Whether the PIN pass needs to model this recirculation, or whether the current gravitropism
is sufficient, is an empirical question. Evaluate after base PIN is live.

**Phototropic PIN redistribution.** Light causes lateral PIN redistribution, bending the shoot
toward light via a lateral auxin gradient. The sim handles phototropism geometrically. Adding
phototropic PIN redistribution is a potential future milestone.

**Wound response.** Stem removal causes auxin accumulation above the cut and vascular
regeneration — a PIN-mediated response. The PIN pass infrastructure would support this later.
