# Milestone 1: Consolidate and Standardize

**Goal:** Enforce a naming convention and modular organization that gives the project a strong, reliable base before any of the milestone-2 work begins.

The codebase has grown organically and shows it. The architecture is sound — flat node hierarchy, unified transport, separate produce/grow steps — but the surface-level organization has accumulated inconsistencies that will compound as new tissue types and systems are added. Better to resolve them now than to inherit them at scale.

---

## Naming Inconsistencies to Resolve

**Node type enum vs class names.** The `NodeType` enum has `APICAL` and `ROOT_APICAL`, but the classes are `ShootApicalNode` and `RootApicalNode`. Code that uses `APICAL` to mean "any apical meristem" is ambiguous. Convention should be: enum values match class names unambiguously. Proposal: `SHOOT_APICAL`, `ROOT_APICAL`, `SHOOT_AXILLARY`, `ROOT_AXILLARY` consistently across the board. Verify current enum values in `node.h` before changing anything.

**Axillary enum gap.** CLAUDE.md lists `SHOOT_AXILLARY` and `ROOT_AXILLARY` as enum values, but a codebase exploration found only 5 types in the enum. Verify whether axillary types exist in the enum or only as class names. If missing from the enum, that's a gap — the downcast helpers `as_shoot_axillary()` / `as_root_axillary()` would be unreachable via type switch.

**`helpers.h` location.** Shared helpers (`growth_direction`, `branch_direction`, `perturb`, `sugar_growth_fraction`, etc.) live in `src/engine/node/meristems/helpers.h`. These are used by leaves and stems too (verify). If so, they belong at `src/engine/node/helpers.h` — not inside the meristems subfolder. A utility file scoped to the wrong directory is quietly misleading.

**`hormone.h/cpp` is a placeholder.** The transport logic moved to `Node::transport_chemicals()`. The file either needs a real role (home for the `transport_chemical()` free function, chemical enum, or BiasedTransportParams struct) or should be deleted and its includes cleaned up. A file that exists but does nothing is a maintenance hazard — it adds noise to grep results and suggests something is there when it isn't.

**`chemical.h` path.** Chemicals are referenced everywhere; confirm the header is at `src/engine/chemical/chemical.h` and that this path is used consistently, not duplicated or re-declared in multiple places.

**App and binary naming.** Apps: `app_realtime.cpp`, `app_headless.cpp`, `app_playback.cpp`, `app_evolve.cpp`, `app_sugar_test.cpp`. Binaries: `botany_realtime`, `botany_evolve`, `botany_tests`, `botany_sugar_test`. The headless and playback binaries aren't documented. Clarify and add to CLAUDE.md.

---

## Organizational Rough Edges

**`node/` vs `node/meristems/`.** The split is logical but it doesn't scale. Leaf, stem, and root nodes sit directly in `node/`; meristems are in `node/meristems/`. As the tissue library grows to 10–15 types, this structure breaks. Decide now what the target directory looks like — probably `node/tissues/` with subdirectories by functional group (meristems, structural, specialized) — even if you don't move files yet. Having the target in mind prevents a second reorganization later.

**`src/evolution/` coupling.** `genome_bridge.h` converts between `Genome` and `evolve::StructuredGenome`. As Genome grows to include body plan fields and per-tissue parameter sections (milestone 3), the bridge will become the most complex file in the project. Decide now whether it lives in `engine/` (it's a Genome serialization concern) or stays in `evolution/` (it's evolution-specific plumbing). The boundary should be intentional.

**`world_params.json` at repo root.** Convenient but informal. Once there are multiple world configurations (desert, forest, gradient environments for coevolution), these belong in `data/worlds/` or similar. Not urgent, but worth noting now so the path doesn't get baked into scripts.

**Test coverage.** `botany_tests` exists — audit what it actually tests. The sugar economy tester (`botany_sugar_test`) is a standalone app. Before adding any new tissue types in milestone 2, there should be a test path for each. A tissue type without a test means bugs hide until something looks wrong in the visual sim, which is much harder to debug than a unit failure.

---

## Recommended First Steps

1. Read `node.h` and enumerate the actual `NodeType` enum values. Compare against CLAUDE.md. Resolve the mismatch and update CLAUDE.md to match code (not the other way around).
2. Decide the fate of `hormone.h/cpp` — delete or give it a specific role. Pick one and do it.
3. Grep for `helpers.h` includes across the codebase. If stem or leaf nodes include it, move the file up to `node/helpers.h`.
4. Sketch the target `node/` directory structure on paper (or in a doc) that scales to 15 tissue types. Don't move files yet — just agree on the destination.
5. Run `botany_tests` and document (briefly, inline or in a test README) what each test covers and what's missing.
