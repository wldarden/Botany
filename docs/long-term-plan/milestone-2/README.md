# Milestone 2: Physics Sandbox + Complete Tissue Library

**Goal:** Build a world physics sim robust enough to be a real sandbox, and a tissue library general enough that most real plants could be modeled by tuning genome parameters.

These are two parallel tracks that partially gate each other. The tissue library can be built without world physics, but some tissue types (water-storage parenchyma, root hair uptake) only mean something once the water model exists. The world physics sim can be scaffolded without the full tissue library, but the selective pressure it creates won't be interesting until plants can respond to it in structurally distinct ways.

---

## Two Tracks

### Track A: World Physics

What the sim needs to be a real ecological sandbox:

- **Soil water model** — a 2D grid at ground level; rain events add water, roots deplete it locally, soil diffuses it slowly, evaporation removes it. Creates the drought/flood gradient that cactus vs tree strategies are answers to.
- **Evapotranspiration** — photosynthesis costs water. Currently free. Without this cost there's no tradeoff between large leaves (high sugar) and small leaves (low water loss).
- **Light occlusion between plants** — tall plants should shade shorter ones. Currently each plant's light exposure is independent. Without this, height has no competitive advantage.
- **Spatial environmental variation** — a moisture gradient (wet zone vs dry zone) across the world creates distinct niches. Without variation, one strategy dominates everywhere.

[→ Full world physics doc](world-physics.md)

### Track B: Tissue Library

Expanding the current 5 node types to cover the major plant body plans:

- **Intercalary meristem** — basal growth zone; inverted insertion logic vs apical; unlocks grasses
- **Water-storage parenchyma** — turgor/volume model; separate water resource; unlocks succulents and cacti
- **Rhizome** — horizontal underground stem; periodic shoot initiation; unlocks ferns and spreading perennials
- **Storage organ** — high sugar/water capacity, near-zero maintenance, dormancy cycle; unlocks geophytes
- **Photosynthetic stem** — a parameter addition, not a new node type; unlocks leafless cacti
- **Tendril** — contact sensing + differential growth; unlocks vines

[→ Full tissue library doc](tissue-library.md)

---

## Keystone Addition: Water

Water is the working fluid of every plant. Adding it as a resource at every node — not just at leaves — is the single change that makes the most other things possible. It enables turgor-based structural mechanics, the evapotranspiration cost, the water-storage parenchyma tissue type, the soil model, and partial coupling to hormone transport.

[→ Full water model doc](water-model.md)

---

## Visual Layer

A secondary but important goal for this milestone: genome parameters that affect behavior should also affect visual appearance. The stem cuticle parameter is the strongest example — it controls water loss from stem surfaces AND determines whether a plant looks waxy/succulent, dull/barky, or wet/tropical. Several other parameters have the same dual-use potential.

[→ Stem cuticle and visuals](stem-cuticle-and-visuals.md)

---

## Validation Exercise

The test for whether milestone 2 is complete: can you produce a cactus-like and a tree-like plant from the same engine, in a world where their respective environments create genuine selection pressure for their strategies?

[→ Cactus vs tree case study](cactus-vs-tree.md)

---

## Success Criteria

Milestone 2 is done when:

1. Water is a tracked resource at every node, with root uptake from soil and leaf transpiration loss.
2. Light occlusion between plants is computed each tick (or binned efficiently enough to approximate per-tick).
3. The tissue library includes at minimum: intercalary meristem, water-storage parenchyma, rhizome, photosynthetic stem parameter.
4. A cactus-like genome produces a distinct, recognizable morphology in a dry-zone environment.
5. A tree-like genome produces a distinct, recognizable morphology in a well-watered, lit environment.
6. These are not hand-tuned special cases — they emerge from the same engine and physics, just different genome files.
