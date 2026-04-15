# Tissue Library

The current 5 node types model an indeterminate dicot tree reasonably well. Getting to a complete library means covering the major growth strategies that distinguish real plant taxa from one another. This doc catalogs what exists, what's missing, and — critically — which gaps require new coded logic vs which are just tunable parameters on existing types.

---

## Universal Tissue Catalog

Every vascular plant runs on the same four functional systems. Here's how they map to the current sim:

| Functional system | What it does | Current sim coverage | Gaps |
|---|---|---|---|
| Meristematic | Undifferentiated dividing cells; source of all new growth | APICAL, ROOT_APICAL, SHOOT_AXILLARY, ROOT_AXILLARY | Missing: intercalary meristem, cambium |
| Vascular | Xylem (water up) + phloem (sugar down) | Phloem implicitly modeled via radius-bottlenecked sugar transport. No xylem. | No water resource, no pressure-driven upward transport |
| Photosynthetic | Chlorenchyma; converts light + CO2 + water to sugar | LEAF nodes only | No green stem photosynthesis |
| Ground/parenchyma | Bulk tissue: structural support, storage, cortex | STEM and ROOT implicitly serve this role | No distinction between storage parenchyma, structural wood, dead xylem |
| Protective/epidermal | Epidermis, cuticle, bark, periderm | Not modeled | Water loss through surfaces, mechanical protection |
| Absorptive | Root epidermis, root hairs | Not modeled | Roots grow but don't acquire water/nutrients from soil |

The current sim handles meristematic and photosynthetic tissue well, is a rough stand-in for parenchyma, has an implicit phloem model but no xylem, and is missing protective and absorptive tissue entirely.

---

## Currently Covered

| Tissue | Node type | Notes / gaps |
|---|---|---|
| Shoot internode | STEM | Secondary thickening present (auxin-driven). No bark/xylem distinction. |
| Root internode | ROOT | Similar to STEM minus stress physics. No soil water uptake. |
| Leaf blade | LEAF | Photosynthesis, abscission, light response. No sheath, no needle morphology, no compound structure. |
| Shoot apical meristem | SHOOT_APICAL | Chain growth, phyllotaxis. Phyllotaxis hardcoded to golden angle. |
| Root apical meristem | ROOT_APICAL | Gravitropism, chain growth. OK. |
| Lateral shoot bud | SHOOT_AXILLARY | Auxin-gated activation. OK. |
| Lateral root initiation | ROOT_AXILLARY | Cytokinin-gated activation. OK. |

---

## Missing Tissues: What Needs New Code vs What's Tunable

### Parameter-only additions (no new node type)

**Photosynthetic stem.** Some plants — cacti, euphorbias, many succulents — photosynthesize primarily through green stems. Leaves are absent or vestigial. Currently only LEAF nodes call `produce()` for sugar. Adding a `stem_photosynthesis_rate` genome parameter and calling it in `StemNode::produce()` covers this without a new type. Rate scales with stem surface area and light exposure. When `stem_photosynthesis_rate = 0` (default), behavior is unchanged. When > 0, STEM nodes contribute sugar proportionally. This unlocks leafless cacti and allows gradual evolution from leaf-dominant to stem-dominant photosynthesis.

**Root surface area multiplier.** Root hair density dramatically increases water/nutrient uptake surface. Modelable as a genome parameter `root_hair_density` that scales per-tick uptake rate in ROOT nodes. Not a new tissue type; just a ROOT parameter. Relevant once soil resources exist.

**Stem cuticle thickness.** Controls water loss through stem surfaces. See [stem-cuticle-and-visuals.md](stem-cuticle-and-visuals.md). Also a genome parameter on STEM nodes, not a new type.

---

### New coded logic required

These tissues can't be expressed by tuning existing node parameters. Each needs new behavior in the engine.

---

**Intercalary meristem** — unlocks grasses and monocots

The key architectural difference from apical growth: new cells are inserted at the *base* of an internode, not the tip.

Apical model (current): the meristem node sits at the tip. It divides, deposits a new internode below itself, and moves upward. The oldest tissue is at the base; the newest is just below the advancing tip. Nodes' positions are set at creation and don't change (except via stress/droop mechanics).

Intercalary model: a band of dividing cells sits at the *base* of a leaf or internode — embedded between mature tissues. It divides and pushes the segment above it upward. The meristem doesn't move. Multiple intercalary zones can be active simultaneously — bamboo can grow nearly a meter per day because every internode is elongating at once from its own basal zone.

**What needs new code:** In the current model, once a node is placed its position is essentially fixed. Intercalary growth requires a node that continues to grow in *length* after nodes above it have already matured. As the intercalary zone extends, everything above it on that branch gets pushed upward — a position cascade. That's inverted from the current insertion logic. Implementation would need:
- An `INTERCALARY` node type with a `growth_zone_active` flag and a `zone_length` accumulator.
- Each tick the intercalary zone adds length from its base end, then adjusts world-space positions of all nodes above it on the branch.
- That upward position adjustment is the new bookkeeping that doesn't exist today.

Grasses also suppress apical dominance (the main shoot tip differentiates into a flower early; subsequent growth is entirely from intercalary zones and tillering). This is a genome-tunable behavior — very high auxin decay rate near the apex — not a new tissue type.

---

**Water-storage parenchyma** — unlocks succulents and cacti

This is the tissue that makes a cactus a cactus. It cannot be faked with `sugar_storage_density_wood` set high. The differences:

- Water is stored as a *separate resource* from sugar — a distinct chemical with its own pool per node.
- Storage is volumetric: the node physically swells when hydrated (larger radius, visible turgor) and shrinks in drought.
- The stored water supports turgor pressure, which is what keeps the stem rigid and the plant functional during drought.
- Water is consumed during metabolism, transpiration, and cell expansion — it has a usage economy independent of sugar.

**What needs new code:**
- Water as a `ChemicalID` with per-node capacity proportional to node volume (much higher capacity in succulent parenchyma than in woody stems).
- Turgor model: effective node radius or a rigidity modifier scales with `water_level / water_capacity`.
- Storage parenchyma specifically has: high `water_capacity`, low `water_loss_rate` (dense cuticle), near-zero active growth (it's bulk storage, not growing tissue).
- This tissue type is part of the STEM node in the sense that it's a stem structural role, but its parameters are so different from woody tissue that it probably warrants a distinct node type (`SUCCULENT_PARENCHYMA` or just making it a STEM variant with a `parenchyma_mode` flag).

See [water-model.md](water-model.md) for the full water resource design.

---

**Rhizome** — unlocks ferns, grasses, spreading perennials

A horizontal underground stem that spreads laterally and periodically initiates new aerial shoots. Structurally similar to ROOT but behaviorally different:
- Grows horizontally (gravitropism disabled or lateral).
- Not seeking water/nutrients in the same way — it's spreading to colonize space.
- At intervals along its length, initiates new vertical shoot meristems pointing upward.
- Stores sugar as an energy reserve for re-sprouting.

**What needs new code:**
- A `RHIZOME` node type (or STEM variant) with horizontal-preferring direction and periodic shoot initiation logic.
- The shoot initiation interval is genome-controlled (ticks between new shoot meristems, analogous to plastochron but along a horizontal axis).
- Rhizome nodes live underground (y < 0) but aren't root tissue — they need different physical parameters (no root gravitropism, different radius/density).

---

**Storage organ (bulb, corm, tuber)** — unlocks geophytes

A node or cluster of nodes with very high sugar and water storage capacity, near-zero maintenance cost during dormancy, and the ability to re-sprout after a dormancy period. Tulips, dahlias, potatoes.

**What needs new code:**
- A `STORAGE_ORGAN` node type with `dormancy_active` flag.
- During dormancy: maintenance cost drops to near zero, no growth, no chemical production.
- Dormancy trigger: day length shortens or temperature drops below threshold (requires seasonal signal in world).
- Re-sprouting: when dormancy ends, node initiates a new shoot meristem from stored energy.
- The dormancy/regrowth cycle doesn't exist in the current model at all — axillary buds persist but there's no "winter" for them to wait out, and no mechanism for stored-energy regrowth.

---

**Tendril** — unlocks vines and climbers

A modified stem or leaf tip that coils around supports by sensing contact and growing differentially.

**What needs new code:**
- Contact detection: the tendril node queries nearby nodes (other plants' nodes, or world objects) within a small radius each tick.
- On contact: growth rate on the contact side slows while growth on the away side continues — creating a differential that generates the coil.
- Coiled tendrils are rigid once set; they behave like cables under tension, transferring mechanical load from the vine to its support.

This is mechanically the most novel of the missing tissue types — it's the only one that requires sensing and responding to external objects. Probably lowest priority unless vines are a specific near-term target.

---

## Priority Order

1. **Intercalary meristem** — highest structural impact; unlocks grasses (largest biomass on Earth); the position-cascade implementation is hard but finite.
2. **Water-storage parenchyma** — can't build cactus strategy without it; also the keystone of the water model.
3. **Rhizome** — high variety coverage; ferns, grasses, most perennials; fairly contained new code.
4. **Photosynthetic stem** — parameter addition, not a new type; unlocks succulents cheaply.
5. **Storage organ / dormancy** — requires seasonal signal from world; blocked until that exists.
6. **Tendril** — most novel code; lowest immediate coverage gain.
