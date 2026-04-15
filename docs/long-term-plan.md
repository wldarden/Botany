# Long-Term Plan (moved)

This file has been superseded. See [docs/long-term-plan/README.md](long-term-plan/README.md) for the current plan, organized by milestone with per-concept files.

---

<!-- archived content below — kept for reference during transition -->

# Long-Term Plan: From Tree Simulator to Coevolving Ecosystem

This document captures the goal hierarchy for this project and is meant to be reread often. Goals are ordered from farthest away to most immediate. The current work should always be traceable to something higher up this chain.

---

## Far Future: Coevolving Plant Ecosystem

**Goal:** Initialize an environment with many different kinds of plants and watch them coevolve. Success looks like seeing them differentiate into recognizable real-world plant types — without those types being hand-designed.

This is the experiment the whole project is building toward. You populate a world with plants whose genomes start similar (or random), apply ecological pressure — light competition, soil resource limits, wind stress, space — and let evolution run. The hope is that distinct strategies emerge: something that grows fast and tall to shade competitors, something that spreads wide along the ground, something that invests heavily in roots and weathers droughts the tall plants can't. If those strategies correspond to recognizable plant archetypes (tree, grass, vine, succulent), that's genuine emergent morphological diversity.

**What this requires:**

- A physics/world sim with real ecological pressure — competition for light, water/soil resources, space. Without scarcity, there's nothing to differentiate against.
- A genome format expressive enough that evolution can discover structurally different body plans, not just tune parameters within one body plan.
- A fitness signal that rewards ecological success (resource capture relative to neighbors) rather than anything hand-crafted per plant type.
- A speciation pressure mechanism, or at minimum enough genome variability that the population doesn't collapse to one optimum.

**Open questions:**

- What's the fitness signal? Pure sugar accumulation favors fast-growing trees in most conditions. Survival in variable environments (drought cycles, seasonal light) might produce more interesting differentiation.
- How do you prevent one strategy from dominating before others can get established? Real ecosystems have founder effects, microhabitat variation, disturbance regimes. Simulation needs some analog.
- At what population size does interesting coevolution happen? Probably much larger than what's currently practical to simulate in real time.

---

## Future: Genome-File-Driven Plant Types

**Goal:** Craft genome files whose parameters cause the same engine and tissue nodes to generate visually and behaviorally distinct plant types — vines, bushes, trees, grasses, ferns — without adding new code for each type.

The bet here is that a complete-enough tissue library, combined with a genome format that specifies body plan and per-tissue parameters, gives you enough degrees of freedom to approximate most real plant strategies. A vine genome would dial up tendril sensitivity and lean the auxin gradient to favor lateral extension. A grass genome would suppress apical dominance, favor basal meristem activity, and produce many thin tillers. A tree genome would do what the current system already does reasonably well.

**What a genome file needs (sketch, not a spec):**

*Body plan / initial topology.* Today `Plant::Plant()` hardcodes: one STEM seed, one APICAL child, one ROOT_APICAL child. To get grass-like tillering or a bulb that sends up multiple shoots, the genome needs to specify germination topology. Something like: "seed type, initial meristem list with directions and types." This doesn't have to be arbitrary graph structure — a short declarative list of (type, direction, initial-hormone-load) tuples covers most cases.

*Tissue-type selection.* The plant only instantiates tissues it uses. A grass genome references BASAL_MERISTEM; a vine genome references TENDRIL. Tissue types the genome doesn't reference stay out of the simulation for that plant. This is how the engine stays general while individual plants stay lean.

*Per-tissue parameters.* The current 64-field genome applies uniformly. With a tissue library, you need sections: `[stem]`, `[leaf]`, `[rhizome]`, etc., each with its own elongation rate, hormone sensitivities, sugar costs. Parameters shared across all tissues (hormone kinetics, transport rates) stay at the top level.

**What's currently missing:**

- No body plan field in the genome — seed structure is hardcoded.
- No tissue-type selection — all 5 node types are always active.
- No per-tissue parameter sections — the genome is flat.
- Phyllotaxis is hardcoded to golden angle (2.399 rad in `helpers.h`); a genome file can't specify opposite leaves, whorled branching, or distichous (grass-like 2-ranked) arrangement.
- No tendril, rhizome, stolon, basal meristem, or storage tissue types exist yet (see "Soon" section below).

This goal is blocked on the "Soon" goal (tissue library) being substantially complete.

---

## Soon: Physics Sandbox + Complete Tissue Library

This is the next real engineering phase. Two parallel tracks, neither strictly blocking the other.

### Track 1: World / Physics Sandbox

The current simulation has gravity (drooping, stress), light exposure (per-leaf, based on position), and basic sugar economy. That's a real foundation. What it's missing to make plant behaviors genuinely meaningful:

**Light competition.** Right now each leaf gets `light_exposure` based on its position relative to some world light level, but there's no shading of one plant by another. A tall tree should reduce light for shorter plants beneath it. This probably means a simple vertical occlusion model — project each node downward, accumulate shadow at lower nodes. Doesn't need to be ray-traced; even a coarse grid per tick would create real selection pressure.

**Soil / water.** Roots currently grow downward and produce cytokinin, but they don't actually acquire any resource from the soil — there's no soil model. A minimal version: a water/nutrient field on a 2D grid at ground level, roots deplete it locally, depletion slows growth. This creates root competition and rewards root spread.

**Wind / mechanical load.** The stress/droop/break physics exist but are driven by gravity only. Wind would add a directional mechanical load, vary by height (taller = more exposed), and create selection pressure for stem thickness and flexible wood. Even a simple sinusoidal force applied to tall nodes would produce meaningful effects.

**Seasonal / variable environment.** A day/night cycle (light_level oscillates) is probably already easy to add. Seasons (light_level and maybe temperature proxy affecting growth rates) would create dormancy pressure. Without variability, plants just optimize for the steady-state environment they're born into.

**Space competition.** Crowding detection exists via ethylene (ethylene_crowding_radius). Verify whether this currently affects growth of *neighboring plants*, or just within one plant's own nodes. Cross-plant crowding is what creates real space competition.

### Track 2: Complete Tissue Library

The current 5 node types (STEM, ROOT, LEAF, APICAL, ROOT_APICAL) model an indeterminate dicot tree reasonably well. To cover most of the plant kingdom, here's a rough taxonomy of what's needed and what's missing:

**Currently covered (with caveats):**

| Tissue | Current node type | Gaps |
|--------|-------------------|------|
| Shoot internode | STEM | OK |
| Root internode | ROOT | OK |
| Leaf | LEAF | No sheath, no compound leaf, no needle morphology |
| Shoot apical meristem | APICAL | Hardcoded golden-angle phyllotaxis |
| Root apical meristem | ROOT_APICAL | OK |
| Lateral shoot bud | SHOOT_AXILLARY | OK |
| Lateral root initiation | ROOT_AXILLARY | OK |

**Missing tissues for broader plant coverage:**

*Basal / intercalary meristem.* Grasses and many monocots grow from the base of internodes, not the tip. The current model is purely tip-driven. An INTERCALARY_MERISTEM type that inserts itself at the base of an internode rather than the apex would enable grass-like growth. This is architecturally different from APICAL — it changes *where* new cells are added, not just what direction the tip points.

*Rhizome.* A horizontal underground stem that produces shoots at nodes. Structurally similar to ROOT but: grows horizontally (not gravitropically), produces aerial shoots at intervals, stores sugar. Could be modeled with a RHIZOME node type that behaves like STEM but with a horizontal gravitropism override and periodic shoot-meristem spawning. Enables grasses, ferns, many perennials.

*Stolon.* Like a rhizome but above-ground and typically thinner (strawberry runner). Could share a type with RHIZOME, differentiated by genome parameters (above vs below soil surface detection, node spacing, shoot vs root initiation at nodes).

*Tendril.* A modified stem or leaf that coils around supports. Mechanically: senses contact, shortens on the contact side to produce curl. Biologically complex but the key capability is *contact sensing* + *differential growth response*. Enables vines and climbers.

*Storage tissue (bulb, corm, tuber).* A node or cluster of nodes with very high sugar storage capacity, low maintenance cost, and the ability to persist through dormancy and re-sprout. The current sugar_storage_density_wood parameter is global; dedicated storage tissue would have its own much higher storage density and near-zero maintenance rate.

*Cambium / vascular tissue (secondary growth).* Secondary thickening exists in STEM nodes (thickening_rate, auxin-driven). This is already a rough model of cambium activity. The main gap is that real secondary growth also produces bark (protective tissue with different properties) and changes the hydraulic capacity of older stems. May not be necessary for first-pass plant diversity.

*Spine / thorn.* A hardened, non-growing modified stem or leaf tip. Primarily defensive (herbivory modeling) — probably low priority until there are herbivores in the simulation.

*Photosynthetic stem.* Some plants (cacti, many succulents) photosynthesize primarily through green stems, not leaves. Currently only LEAF nodes produce sugar. A flag or tissue variant that allows STEM nodes to photosynthesize would cover succulents and leafless cacti. Probably achievable as a genome parameter (`stem_photosynthesis_rate`) rather than a new node type.

*Root hair zone / mycorrhizae.* Fine roots and fungal associations dramatically increase water/nutrient uptake. Could be modeled as a ROOT node parameter (surface_area_multiplier) rather than a distinct type. Relevant once soil resources are simulated.

**Priority order for tissue addition (suggestion):**

1. Intercalary/basal meristem — unlocks grasses, the largest plant group on earth
2. Rhizome — unlocks most perennial spreading plants, ferns
3. Storage tissue — unlocks geophytes (bulbs, corms, tubers), dormancy cycles
4. Photosynthetic stem parameter — unlocks succulents without new node type
5. Tendril — unlocks vines, adds interesting mechanical complexity
6. Stolon — likely falls out naturally after rhizome

---

## Current Goal: Consolidate and Standardize

**Goal:** Enforce a naming convention and modular organization that gives the project a strong, reliable base before any of the above work begins.

The codebase has grown organically and shows it. The architecture is sound — flat node hierarchy, unified transport, separate produce/grow steps — but the surface-level organization has accumulated inconsistencies that will slow down every future change if not addressed now.

### Naming inconsistencies to resolve

**Node type names in the enum vs class names.** The `NodeType` enum has `APICAL` and `ROOT_APICAL`, but the classes are named `ShootApicalNode` and `RootApicalNode`. Older code sometimes uses `APICAL` to mean "any apical meristem" which creates ambiguity. Convention should be: enum values match class names unambiguously. Proposal: `SHOOT_APICAL`, `ROOT_APICAL`, `SHOOT_AXILLARY`, `ROOT_AXILLARY` across the board. (Verify current enum values in `node.h` — some of these may already be corrected.)

**"Axillary" node split.** The CLAUDE.md lists `SHOOT_AXILLARY` and `ROOT_AXILLARY` as enum values, but the Explore agent found only 5 types. Verify whether the axillary types are in the enum or only in code as class names. If they're missing from the enum, that's a gap to close.

**`helpers.h` in meristems/.** The shared helper functions (`growth_direction`, `branch_direction`, `perturb`, `sugar_growth_fraction`) live in `src/engine/node/meristems/helpers.h`. These aren't meristem-specific — they're used by leaves and stems too (verify). If so, move to `src/engine/node/helpers.h` or `src/engine/node/growth_helpers.h`.

**`hormone.h/cpp` is a placeholder.** The file exists but is described as empty — transport logic moved to `Node::transport_chemicals()`. Either delete it and remove includes, or give it a real role (e.g., the `transport_chemical()` free function lives there as a utility). A file that exists but does nothing is a maintenance hazard.

**`chemical.h` location.** Chemicals are referenced everywhere but the header location isn't canonical in CLAUDE.md (it mentions `src/engine/chemical/chemical.h` in passing). Verify this exists and that the include path is used consistently, not duplicated or re-declared.

**App naming.** The apps are: `app_realtime.cpp`, `app_headless.cpp`, `app_playback.cpp`, `app_evolve.cpp`, `app_sugar_test.cpp`. The binaries are: `botany_realtime`, `botany_evolve`, `botany_tests`, `botany_sugar_test`. The headless/playback apps don't appear to have binary names documented. Clarify and document in CLAUDE.md.

### Organizational rough edges

**`src/engine/node/` vs `src/engine/node/meristems/`.** The split is logical but meristems-only. Leaf, stem, root nodes are directly in `node/`, meristems are in `node/meristems/`. As the tissue library grows, this structure will need to scale — probably to `node/tissues/` with subdirectories per tissue family. Think about that directory structure now, before there are 12 tissue types to reorganize.

**`src/evolution/` tightly coupled to engine internals.** `genome_bridge.h` converts between `Genome` and `evolve::StructuredGenome`. As the Genome struct grows to include per-tissue sections and body plan fields, the bridge will become complex. Consider whether `genome_bridge` should live in `engine/` (it's really a Genome serialization concern) or stay in `evolution/` (it's evolution-specific). Either is fine, but the boundary should be intentional.

**`world_params.json` at repo root.** This file configures the non-genetic physical constants (light level, construction costs, maintenance rates). Having it at the repo root is convenient but informal. Once there are multiple world configurations (different environments for coevolution experiments), these should probably live in `data/` or `worlds/`. Not urgent, but worth noting.

**Test coverage.** `botany_tests` exists — verify what it actually tests. The sugar economy tester (`botany_sugar_test`) is a separate app. Before adding new tissue types, make sure there's a test path for each one. Adding a tissue type without a test for its basic behavior (does it grow? does it transport hormones correctly? does it produce/consume sugar?) means bugs will hide until the visual sim makes something look wrong.

### Recommended first steps

1. Audit the `NodeType` enum against the actual class names — resolve any mismatch, document the canonical list in CLAUDE.md.
2. Decide the fate of `hormone.h/cpp` — delete or give it a clear role.
3. Move `helpers.h` to the right location based on what actually uses it.
4. Sketch the `node/` directory structure that will scale to a full tissue library (even if you don't move files yet — just decide what the target looks like).
5. Read through `botany_tests` and document what's covered and what's not.

---

## Appendix: Tissue Type Reference

Quick reference for the tissue types discussed above, for use when designing the library.

| Tissue | Real botany name | Present in sim? | Priority |
|--------|-----------------|-----------------|----------|
| Shoot apical meristem | SAM | Yes (APICAL) | — |
| Root apical meristem | RAM | Yes (ROOT_APICAL) | — |
| Lateral shoot bud | Axillary bud | Yes (SHOOT_AXILLARY) | — |
| Lateral root initiation | Lateral root primordium | Yes (ROOT_AXILLARY) | — |
| Shoot internode | Internode/stem | Yes (STEM) | — |
| Root internode | Root segment | Yes (ROOT) | — |
| Leaf blade | Lamina | Yes (LEAF) | — |
| Intercalary meristem | Intercalary meristem | No | High |
| Horizontal underground stem | Rhizome | No | High |
| Storage organ | Bulb/corm/tuber | No | Medium |
| Photosynthetic stem | Cladode / chlorenchyma | No (param only) | Medium |
| Above-ground runner | Stolon | No | Medium |
| Coiling contact organ | Tendril | No | Low |
| Protective hardened tip | Spine/thorn | No | Low |

*"Priority" reflects usefulness for morphological diversity, not biological importance.*
