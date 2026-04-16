# Vascularization

Plants have two fundamentally different transport systems, and the sim should honor that distinction. Getting this right is what makes canalization meaningful — not just an interesting internal detail, but the actual mechanism by which auxin's slow cell-to-cell journey shapes the plant's global logistics infrastructure.

---

## The Goal

The transport architecture should be two distinct systems that operate at different scales and through different mechanisms:

**Vascular transport** — long-distance bulk logistics. Sugar moves from canopy to root tips. Water moves from roots to transpiring leaves. Cytokinin moves from root apices to shoot apices. These flows happen at scale, through pipes, under pressure. The pipes are shaped by auxin canalization history: connections where auxin has flowed repeatedly over time develop structural reinforcement — thicker, more conductive tissue — and those connections carry proportionally more of the global flow budget. Vascular transport is fast relative to diffusion, long-range by design, and structurally committed. Once a pipe is built it persists.

**Local diffusion** — cell-to-cell signaling. Auxin, gibberellin, and stress signals move node-by-node through the tree graph, slowly. The rate is low enough that these chemicals are genuinely local: a signal produced at the shoot tip takes many ticks to reach even a few nodes away, and never reaches the opposite end of the plant under normal conditions. This is what creates the hormone gradients that drive growth decisions. Apical dominance works because auxin diffuses downward slowly and decays before it reaches lateral buds more than a few nodes below — not because it teleports to every bud at once.

**The bridge between them is canalization.** Auxin doesn't travel through the vascular system — it travels by polar active transport, cell-to-cell, via PIN proteins on specific cell faces. This is the local diffusion channel. But wherever auxin has flowed repeatedly along that slow local path, structural reinforcement accumulates. Over time those connections become the plant's vascular highways. Canalization is how a slow local signal builds a fast global infrastructure — it is the developmental history that determines which branches become main axes and which remain lateral twigs. The vascular pass uses that history when distributing flow at junctions.

Biologically, this distinction is precise. Polar auxin transport (PAT) is an active, direction-specific, cell-to-cell process driven by asymmetric PIN efflux carrier proteins on specific faces of each cell. It is slow because it requires a protein-mediated step at every cell boundary. Vascular bulk flow is mechanically driven: phloem by osmotic pressure gradients (sugar loading at source, unloading at sink), xylem by transpiration pull from leaf surfaces. These are physically different phenomena — diffusion versus pressure flow — and they operate at different timescales and distances.

---

## What Changes From Current State

Three specific changes bring the sim into alignment with this architecture.

**1. Vascular Phase 2 weights by structural_flow_bias at junctions.**

The vascular two-phase pass currently distributes flow among children in proportion to their subtree demand. This is correct for a naive pipe model, but misses the structural history. A connection where auxin has flowed heavily over hundreds of ticks has developed more conductive tissue — it should carry more of the bulk flow budget, all else equal. At each junction in Phase 2, when distributing available flow among children, the demand-proportional share should be multiplied by `(1 + canalization_weight * structural_flow_bias)` for that child connection. This is the same weight used in local diffusion's `transport_with_children` — the same bias map, the same formula — just applied now to the vascular pass as well. A branch that has been a consistent auxin highway carries more sugar and water than a branch of equal size that has been dormant. This is what makes branch hierarchy emerge structurally rather than just topologically.

**2. Reduce local diffusion rates for auxin, GA, and stress.**

These chemicals should be genuinely local. In the current parameter defaults, diffusion rates are set broadly. The target behavior: a signal produced at the shoot apex should take on the order of 50+ ticks to traverse even a modest plant, and it should decay substantially along the way. A leaf's GA signal should affect only its own internode and the one above it — not travel to the canopy top or the seed. This is what makes the internode-local elongation model biologically correct: GA produced by a young leaf acts on the internode being elongated, not on the whole plant's elongation budget.

Auxin's slow local journey is what makes apical dominance work at all. If auxin diffuses fast enough to reach lateral buds three or four nodes below the apex in one or two ticks, apical dominance is uniform across the plant and the branching pattern loses spatial structure. The gradient should be steep: high auxin immediately below the apex, falling off within a few internodes, effectively zero at the branch points where lateral buds activate.

Concrete target rates: auxin diffusion rate reduced from current default toward 0.1–0.15; GA diffusion rate reduced toward 0.05 (GA is almost entirely local); stress diffusion rate reduced toward 0.1–0.15 (mechanical stress should propagate a few nodes but not long distances). These numbers need tuning against observed behavior, but the direction is clear: slower, not faster.

**3. Make the architectural separation explicit in code comments.**

The two-system architecture should be legible in the code. `vascular.cpp` already has the right shape, but the comments should make explicit: (a) which chemicals go through vascular, (b) which go through local diffusion, (c) why, and (d) what the relationship between the two is via canalization. `node.cpp`'s `transport_with_children` should note that it handles only non-vascular chemicals plus last-mile delivery. `plant.cpp`'s `tick_tree` should annotate the call order and why vascular runs before the DFS walk. Comments are load-bearing documentation here because the two-system split is non-obvious to anyone reading the code cold.

---

## Biological Reference

**PIN proteins and polar auxin transport.** PIN (PIN-FORMED) proteins are efflux carriers that move auxin out of cells. They are asymmetrically distributed on specific faces of each cell — auxin exits preferentially from the basal face (toward roots) in most shoot tissues. This polarity is what makes auxin transport directional without requiring any active "decision" by each cell. The polarity itself is dynamic: PIN distribution is regulated by auxin concentration, creating a feedback loop. Where auxin flux is high, PIN proteins redistribute to reinforce that flux direction. This is the cellular basis of canalization.

**Transpiration pull.** Xylem water movement is driven by negative pressure at the top. Leaves transpire water vapor through stomata; this creates a tension that propagates down the continuous water column from leaf mesophyll through xylem vessels to root hairs. The whole column is under tension — xylem water is in a metastable state. This is why xylem is pressure-driven toward the transpiration sink (leaves) rather than toward a concentration equilibrium. Water moves upward not because roots push it but because leaves pull it. In the sim, the slight positive bias on water transport (`water_bias = 0.05`) approximates this directional pull without modeling the full tension gradient.

**Phloem osmotic pressure.** Phloem flow direction is reversed from xylem. Sugar is actively loaded into phloem at source tissues (leaves) and actively unloaded at sink tissues (roots, growing tips). Loading at source raises osmotic pressure there, which draws water in and creates turgor pressure. Unloading at sink lowers pressure. The pressure gradient drives bulk flow from source to sink. This is why phloem can flow bidirectionally depending on which ends are loading and unloading — there is no fixed up/down to phloem flow, unlike xylem.

**Sachs 1969 canalization hypothesis.** T. Sachs proposed that vascular strands form by a positive feedback: auxin moving through a cell induces that cell to become better at conducting auxin, which directs more auxin through it, which further increases conductance. Over time, flux concentrates into discrete channels — the future vascular strands. Undirected cells become parenchyma; cells in the early flux path develop into xylem and phloem. This explains why vascular tissue forms in continuous strands connecting auxin sources to auxin sinks, why leaf veins radiate from the midrib, and why wound-healing can redirect vascular strands around damage. The sim's `structural_flow_bias` (slow, never-decaying ratchet that grows with sustained auxin flux) is a direct implementation of the Sachs mechanism.

**The positive feedback loop in full:** auxin flows cell-to-cell through a tissue → flux above threshold at a connection → structural reinforcement accumulates at that connection → reinforced connection carries proportionally more of all bulk transport (sugar, water, cytokinin) → well-supplied tissue is more metabolically active → more active tissue produces more auxin and grows more → more auxin flows through the already-reinforced connection → further reinforcement. The loop makes early flux advantages self-amplifying, which is why vascular patterns are robust and reproducible despite growing in a noisy environment.

---

## Implementation Steps

These are ordered by dependency. Each step has a clear verification criterion so it can be confirmed working before the next begins.

**Step 1: Apply structural_flow_bias in vascular Phase 2.**

In `vascular.cpp`, in the Phase 2 loop, when computing each child's share of available flow, multiply the demand-proportional share by `(1 + canalization_weight * structural_flow_bias[child])`. The `structural_flow_bias` map is already maintained on each parent node (keyed by child pointer) by `update_canalization()` in `node.cpp`. Access it through the node pointer in `VascNodeInfo`. The `canalization_weight` genome parameter already exists.

Verification: run the realtime viewer with `--color sugar` and grow a plant for several hundred ticks. After branching has occurred and canalization has had time to build up, the primary axis (main stem to main root, or whichever branch has been the consistent auxin highway) should visibly carry more sugar than lateral branches of the same size. Without this change, siblings of equal size carry equal sugar. With it, the main axis should be measurably privileged.

**Step 2: Tune local diffusion rates downward.**

In `chemical_registry.h` (or wherever `diffusion_params` is initialized), reduce the diffusion rates for Auxin, Gibberellin, and Stress to the target range (Auxin ~0.1–0.15, GA ~0.05, Stress ~0.1–0.15). Rebuild and run the realtime viewer with `--color auxin`. The auxin gradient from apex to base should be steep and localized: high near shoot tips, rapidly dropping to low within 3–5 nodes, not spreading evenly across the plant. Test with a plant that has developed multiple branches: lateral buds several nodes below an active apex should sit in low-auxin territory and activate normally; buds immediately below an apex should remain suppressed.

Also verify that GA effects on elongation are now clearly internode-local: young leaves a few nodes from the apex should produce elongating internodes, but that elongation should not travel to internodes elsewhere in the canopy.

Verification: the branching pattern should remain qualitatively similar (same hormones, same thresholds), but branches should initiate at positions that reflect the local gradient rather than a uniform signal. Plants should not become unbranched (auxin too high everywhere) or excessively branched (auxin too low everywhere).

**Step 3: Annotate the two-system architecture in code.**

Add a block comment at the top of `vascular.cpp` explaining the vascular/diffusion split, which chemicals use each path, and the canalization bridge. Add a corresponding note in `node.cpp` near `transport_with_children` explaining that this function handles only non-vascular chemicals (Auxin, GA, Stress) plus last-mile delivery. Update `plant.cpp`'s `tick_tree` comment to explain why vascular runs before the DFS walk (vascular fills vascular-capable nodes; the DFS walk handles meristems and non-vascular tissue that depend on those levels). Add a note in `CLAUDE.md`'s Chemical Transport Model section summarizing the two-system split in the same terms used here.

Verification: someone reading the code cold should understand from the comments alone why `vascular_transport` is a separate pre-DFS pass rather than being folded into the node tick.

**Step 4: Integration test — canalization builds observable vascular hierarchy.**

Write a test (or extend an existing one) that: (1) grows a plant for enough ticks that canalization has built up structural bias on the main axis, (2) reads the structural_flow_bias values on the seed node for its two children (shoot-side and root-side), and (3) verifies that the bias values are non-zero and differ between children in a way consistent with a plant that has been actively growing upward (shoot-side connections should have higher auxin flux history than lateral connections that may not have grown). This is a regression test: it ensures that structural bias actually accumulates over time and is not inadvertently zeroed out by a future change.

Verification: test passes. The structural bias values are observable from outside the transport machinery, so the test can be written as a black-box integration test on the plant state after N ticks.

---

## Relationship to Other Milestone 2 Work

Vascularization interacts with water (milestone 2's keystone chemical) in one specific way: the structural_flow_bias weights should eventually apply to water transport in the vascular pass, just as they will for sugar. Water moves through the same xylem pipes that auxin canalized — the same structural history governs how well-watered each branch is. This is already handled once Step 1 is complete, because the vascular pass runs all three chemicals (sugar, water, cytokinin) through the same `run_vascular` function with the same Phase 2 distribution logic.

The tissue library expansion (intercalary meristems, water-storage parenchyma) does not depend on vascularization. New tissue types will automatically participate in vascular transport once they have mature vascular connections, because `has_vasculature()` gates only on node type (STEM/ROOT) and age — new tissue types that extend those base classes will be handled correctly.
