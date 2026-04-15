# Milestone 4: Coevolving Ecosystem

**Goal:** Initialize an environment with many different kinds of plants and watch them coevolve. Success = seeing them differentiate into recognizable real-world plant types without those types being hand-designed.

This is the experiment the whole project is building toward. A population of plants with varying genomes, placed in a world with spatial environmental variation and real resource competition, evolved under a fitness signal that rewards ecological success relative to neighbors. The hope: distinct strategies emerge — something fast and tall that shades competitors, something low and spreading, something that weathers droughts the tall plants can't survive. If those strategies correspond to recognizable plant archetypes, that's genuine emergent morphological diversity.

---

## What This Requires

Everything from milestones 2 and 3:
- A complete tissue library (milestone 2)
- Real world physics with resource competition (milestone 2)
- A genome format expressive enough that evolution can discover structurally different body plans, not just tune parameters within one (milestone 3)

And additionally:
- **A fitness signal that rewards ecological success.** Pure sugar accumulation favors fast-growing trees in most conditions. Better: fitness proportional to resource capture *relative to neighbors*, or survival across a variable-environment simulation run with drought cycles and seasonal light.
- **Speciation pressure or population structure.** Without it, one strategy dominates and the rest die out before they can establish. Real ecosystems have founder effects, microhabitat variation, and disturbance. The sim needs an analog — probably spatial isolation with limited seed dispersal, or periodic disturbance events that reset parts of the world.
- **Population scale.** Interesting coevolution probably requires more plants than are currently practical in real time. This milestone likely needs background batch simulation capabilities — not just the interactive realtime viewer.

---

## Open Questions

- At what population size does meaningful divergence happen? Probably large enough that background threading or headless batch runs are required.
- How do you prevent premature convergence to one dominant strategy? This is the speciation problem in miniature.
- What's the right fitness signal? Survival-based (live N ticks without dying) vs resource-based (sugar accumulation over a lifetime) vs competitive-based (relative success vs neighbors in same environment)? These produce very different evolutionary pressures.

---

## Contents TBD

Design for the coevolution system, population dynamics, fitness signal architecture, and world event system will be written once milestone 3 is complete. The far-future goal is stable; the implementation path isn't yet clear enough to pre-specify.
