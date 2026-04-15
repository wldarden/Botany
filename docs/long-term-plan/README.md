# Long-Term Plan: From Tree Simulator to Coevolving Ecosystem

Goals are ordered from farthest away to most immediate. Everything being built now should be traceable to something higher up this chain.

---

## The Four Tiers

### Milestone 4 — Far Future: Coevolving Ecosystem

Initialize an environment with many different kinds of plants and watch them coevolve. Success = seeing them differentiate into recognizable real-world plant types without those types being hand-designed.

[→ milestone-4/README.md](milestone-4/README.md)

---

### Milestone 3 — Future: Genome-File-Driven Plant Types

Craft plant genome files whose parameters make the same plant engine tissue nodes generate distinct plant types — vines, bushes, trees, grass, ferns, etc. No new code needed per plant type; just a different genome file.

[→ milestone-3/README.md](milestone-3/README.md)

---

### Milestone 2 — Soon: Physics Sandbox + Complete Tissue Library

Two parallel tracks:

1. A world physics sim robust enough to be a real sandbox — light competition, soil water, evapotranspiration, spatial resource gradients.
2. A tissue library general enough that most real plants could be accurately modeled by tuning genome parameters to match them.

These two tracks gate each other in places: some tissues (water-storage parenchyma) require the water model to mean anything; the water model needs root/soil integration to have pressure.

[→ milestone-2/README.md](milestone-2/README.md)

Concept docs:
- [Water model](milestone-2/water-model.md) — water as a working fluid at every node; turgor, growth, transport, stress
- [Tissue library](milestone-2/tissue-library.md) — universal catalog, gap analysis, what needs new code
- [World physics](milestone-2/world-physics.md) — soil water grid, light occlusion, evapotranspiration, minimum viable set
- [Stem cuticle and visuals](milestone-2/stem-cuticle-and-visuals.md) — cuticle as physics + visual expression; other visual-expressive genome params
- [Cactus vs tree](milestone-2/cactus-vs-tree.md) — divergence case study; validation exercise for whether milestone 2 is complete

---

### Milestone 1 — Current Goal: Consolidate and Standardize

Enforce naming conventions and modular organization before the larger build begins. Resolve existing inconsistencies so they don't compound as the codebase grows.

[→ milestone-1/README.md](milestone-1/README.md)

---

---

## Guiding Principle: Hew to Real Plant Biology

The single most reliable design heuristic this project has found: when in doubt, do what real plants do.

This has come up repeatedly — every time the sim has been brought closer to actual plant hormone signaling (correct auxin directionality, correct apical dominance gating, correct cytokinin-gated root activation), the plants have ended up looking and growing more realistically, without extra tuning. The biology just works. Departures from it tend to work fine in the short term and become design debt when new systems couple in.

The practical rule: when designing a new hormone, resource, producer, consumer, directional transport bias, or activation conditional, look up what real plants do first and default to that. If there's a reason to deviate, note it explicitly so it can be revisited. Shortcuts taken for simplicity now may force a refactor when the next layer of biology is added.

This applies to everything in milestone 2 and beyond: water transport direction (xylem is pressure-driven, acropetal), cytokinin source (roots, not leaves — see `water-model.md`), ABA as a drought-inversion signal, stomatal aperture as a water-gated valve, and intercalary meristem insertion mechanics.

---

## How to Use This Folder

Read the milestone you're actively working on. The concept docs under milestone-2 are the most developed — they capture design decisions and open questions worth revisiting before implementation, not just after. When milestone-2 work begins in earnest, move decisions from those docs into code comments and CLAUDE.md as they're settled.
