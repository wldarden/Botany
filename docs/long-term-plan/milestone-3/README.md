# Milestone 3: Genome-File-Driven Plant Types

**Goal:** Craft plant genome files whose parameters make the same plant engine tissue nodes generate distinct plant types — vines, bushes, trees, grass, ferns, etc. No new code needed per plant type; just a different genome file.

This milestone is blocked on milestone 2 (tissue library) being substantially complete. Once the tissue library exists, the work here is about making it *configurable* — moving the "which tissues does this plant use and how" from hardcoded code into the genome file format.

---

## What Needs to Change

**Body plan / initial topology.** Today `Plant::Plant()` hardcodes: one STEM seed, one SHOOT_APICAL child, one ROOT_APICAL child. To support grass (multiple basal meristems), bulb plants (storage organ at seed), or vine (shoot apical immediately seeking support), the genome needs a body plan field. Probably a short declarative list of (tissue type, direction, initial hormone load) tuples for the germination topology.

**Tissue-type selection.** The plant only instantiates tissues referenced in its genome. A grass genome references INTERCALARY_MERISTEM; a vine genome references TENDRIL. Types not referenced stay out of that plant's simulation. This keeps the engine general while individual genomes stay lean.

**Per-tissue parameter sections.** The current 64-field genome is flat — one value for `internode_elongation_rate`, applied to all stems. With a library of distinct tissue types, you need per-tissue sections: `[stem]`, `[rhizome]`, `[leaf]`, each with its own elongation rate, hormone sensitivities, construction costs. Global parameters (hormone kinetics, transport rates) stay at the top level.

**Phyllotaxis options.** Currently hardcoded to golden angle (2.399 rad). Genome should be able to specify: golden angle (default), opposite, whorled, distichous (2-ranked, grass-style). These are real phyllotaxis patterns and they produce very different architectures.

---

## Contents TBD

Detailed design for the genome format extension, body plan schema, and per-tissue parameter structure will be written once milestone 2 is complete and we know exactly what tissue types exist and what their parameters look like. Designing the format before the tissues exist invites premature optimization.
