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
  few nodes, only those connections ever develop `structural_flow_bias`. The main trunk never
  fully canalizes because auxin never travels through it.

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
   and disappears. This is the cellular basis of `structural_flow_bias` — the canalization data
   already in the sim represents exactly this history of sustained directional PIN-mediated flux.

4. **Speed.** PAT moves auxin at ~5–10 mm/hour in real stems. At the sim's dm/hour scale that
   is ~0.05–0.1 dm/hr — fast enough to traverse a whole plant in 10–20 ticks, far faster than
   diffusion at a 5% gradient per tick.

5. **Tissue-specific direction.** Stems are basipetal (shoot tip → root). Root vasculature is
   acropetal (root base → root tip). The sim needs both.

---

## Design: Dedicated PIN Transport Pass

PIN transport is a new global pass, separate from vascular transport and from local diffusion.
It runs once per tick and moves auxin directionally along the plant axis using `structural_flow_bias`
as the conductance weight for each connection.

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

PIN runs after vascular so that the sugar needed for PIN's metabolic cost is already in place.
PIN runs before the DFS walk so that auxin is in the correct position before `update_tissue`
reads it this tick.

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
float bias = parent ? parent->structural_flow_bias[this] : 0.0f;
float rate = std::clamp(
    g.pin_base_rate * (1.0f + g.canalization_weight * bias),
    0.0f, g.pin_max_rate
);

// Sugar cost gates the pump — no sugar = no PIN activity
float auxin_available = chemical(ChemicalID::Auxin);
float would_move = auxin_available * rate;
float cost = would_move * g.pin_sugar_cost;

if (chemical(ChemicalID::Sugar) >= cost) {
    chemical(ChemicalID::Sugar) -= cost;
    chemical(ChemicalID::Auxin) -= would_move;
    parent->transport_received[ChemicalID::Auxin] += would_move;
} else {
    // Partial pump proportional to available sugar
    float fraction = chemical(ChemicalID::Sugar) / std::max(cost, 1e-8f);
    float moved = would_move * fraction;
    chemical(ChemicalID::Sugar) = 0.0f;
    chemical(ChemicalID::Auxin) -= moved;
    parent->transport_received[ChemicalID::Auxin] += moved;
}
```

The existing `transport_received` buffer is used for PIN exactly as in `transport_with_children`
— anti-teleportation is already in place, and the flushing logic in `Node::tick()` applies to
both.

**Seed junction:**

The seed is always at full PIN conductance (oldest, most-canalized node). After flushing its
`transport_received` buffer from shoot-side children:

```cpp
float total_root_bias = sum of structural_flow_bias[root_child] for all root children;
for (Node* root_child : root_children) {
    float share = (total_root_bias > 1e-8f)
        ? parent->structural_flow_bias[root_child] / total_root_bias
        : 1.0f / root_children.size();
    float to_send = chemical(ChemicalID::Auxin) * share * g.pin_base_rate;
    root_child->transport_received[ChemicalID::Auxin] += to_send;
    chemical(ChemicalID::Auxin) -= to_send;
}
```

Root children with higher canalization history receive a larger share of the auxin transiting
through the seed — the most-developed root branch gets the strongest shoot-derived auxin signal.

**Pre-order distribution (root side):**

Walk root nodes from seed toward root tips. Each root node flushes its `transport_received`
buffer (already filled by its parent's pass) and then distributes to its own root children using
the same conductance-weighted split. This mirrors how the xylem vascular pass distributes water
pre-order, but with PIN-rate conductance instead of pipe cross-section.

### Canalization Coupling

`structural_flow_bias` on each parent-child edge encodes the history of auxin flux through
that connection. Using it directly as the PIN conductance multiplier closes the Sachs feedback
loop:

```
More auxin flux through a connection
  → structural_flow_bias grows (update_canalization, existing)
  → higher PIN rate on that connection
  → more auxin transported through it
  → more bias accumulation
```

New internodes start with a small initial bias stamp (already implemented: ~0.01 at creation).
Until bias accumulates, PIN transport through them is slow — auxin pools near the apex and
existing canalized paths carry most of the flow. As bias grows, PIN conductance increases,
and auxin reaches progressively further down the plant. This matches real plant development
where procambium differentiation precedes mature PAT capacity.

### New Genome Parameters

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `pin_base_rate` | 0.30 | fraction/tick | Fraction of auxin pumped per tick at bias = 0 |
| `pin_max_rate` | 0.80 | fraction/tick | Upper clamp after canalization scaling |
| `pin_sugar_cost` | 0.001 | g glucose / AU | Sugar cost per unit auxin transported |

With `pin_base_rate = 0.30` and `auxin_decay_rate = 0.12`, the effective reach in a chain:

```
ratio_per_hop ≈ pin_rate / (pin_rate + decay_rate)
             = 0.30 / (0.30 + 0.12) ≈ 0.71

At node 10:  0.71^10 ≈ 0.035   — meaningful signal
At node 20:  0.71^20 ≈ 0.001   — very weak but present
```

With canalization scaling at `structural_flow_bias = 1.0` (established main axis):

```
pin_rate = 0.30 × (1 + 1.0 × 1.0) = 0.60
ratio_per_hop ≈ 0.60 / 0.72 ≈ 0.83

At node 10:  0.83^10 ≈ 0.16    — strong signal even 10 nodes down
At node 20:  0.83^20 ≈ 0.026   — still readable at 20 nodes
```

Main axis canalization extends auxin's effective reach from 2–3 nodes (current diffusion alone)
to 15–20 nodes, while lateral branches with low bias stay at the shorter range. This creates
the gradient needed for correct apical dominance.

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
- Near-zero auxin in uncanalized new internodes until a few ticks of bias accumulation

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
- Main shoot axis: continuous auxin flux since germination → high bias → fast PIN conductance →
  more auxin → more bias. Self-reinforcing dominant pathway.
- Lateral branches that activated late: start with low bias → slow PIN → less auxin → bias
  builds slowly. Competitive disadvantage that matches biological reality.

---

## Implementation Plan

**Step 1: Files**

Create `src/engine/pin_transport.h` and `src/engine/pin_transport.cpp`. Declaration:
```cpp
void pin_transport(Plant& plant, const Genome& g);
```

**Step 2: Genome parameters**

In `genome.h`: add `pin_base_rate`, `pin_max_rate`, `pin_sugar_cost`.
In `default_genome()`: set to table values above.
In `genome_bridge.cpp`: add all three to `build_genome_template()` with mutation config.
In `CLAUDE.md`: add to Tuning Parameters section.

**Step 3: Wire into plant.cpp**

In `Plant::tick()`, add `pin_transport(*this, genome_)` immediately after `vascular_transport()`
and before `tick_tree()`. Add `#include "engine/pin_transport.h"`.

**Step 4: Regression test**

Add `tests/test_pin_transport.cpp`:
- Grow a straight-chain plant for 200 ticks
- Assert auxin at seed junction is non-zero (crossed the junction)
- Assert auxin at the root apical is non-zero (reached the root tip)
- Assert auxin gradient direction: shoot apex > seed > root apical
- Assert `structural_flow_bias` has accumulated on main axis connections (canalization is driven
  by the auxin flux that PIN is now delivering)

**Step 5: Update CLAUDE.md**

Add PIN transport to the Chemical Transport Model table:

```
| Auxin | PIN transport (long-range) + local diffusion | Basipetal in stems, acropetal in roots |
```

Add to Tick Control Flow section. Add genome parameter entries for the three new params.

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
