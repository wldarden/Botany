# Vascularization

Plants have two fundamentally different transport systems, and the sim should honor that distinction. Getting this right is what makes canalization meaningful — not just an interesting internal detail, but the actual mechanism by which auxin's slow cell-to-cell journey shapes the plant's global logistics infrastructure. And it unlocks the correct model of thickening: stem width is not a separate growth decision. It is a direct consequence of vascular development.

---

## The Goal

The transport architecture should be two distinct systems that operate at different scales and through different mechanisms:

**Vascular transport** — long-distance bulk logistics. Sugar moves from canopy to root tips. Water moves from roots to transpiring leaves. Cytokinin moves from root apices to shoot apices. These flows happen at scale, through pipes, under pressure. The pipes are shaped by auxin canalization history: connections where auxin has flowed repeatedly over time develop structural reinforcement — thicker, more conductive tissue — and those connections carry proportionally more of the global flow budget. Vascular transport is fast relative to diffusion, long-range by design, and structurally committed. Once a pipe is built it persists.

**Local diffusion** — cell-to-cell signaling. Auxin, gibberellin, and stress signals move node-by-node through the tree graph, slowly. The rate is low enough that these chemicals are genuinely local: a signal produced at the shoot tip takes many ticks to reach even a few nodes away, and never reaches the opposite end of the plant under normal conditions. This is what creates the hormone gradients that drive growth decisions. Apical dominance works because auxin diffuses downward slowly and decays before it reaches lateral buds more than a few nodes below — not because it teleports to every bud at once.

**The bridge between them is canalization.** Auxin doesn't travel through the vascular system — it travels by polar active transport, cell-to-cell, via PIN proteins on specific cell faces. This is the local diffusion channel. But wherever auxin has flowed repeatedly along that slow local path, structural reinforcement accumulates. Over time those connections become the plant's vascular highways. Canalization is how a slow local signal builds a fast global infrastructure — it is the developmental history that determines which branches become main axes and which remain lateral twigs. The vascular pass uses that history when distributing flow at junctions.

Biologically, this distinction is precise. Polar auxin transport (PAT) is an active, direction-specific, cell-to-cell process driven by asymmetric PIN efflux carrier proteins on specific faces of each cell. It is slow because it requires a protein-mediated step at every cell boundary. Vascular bulk flow is mechanically driven: phloem by osmotic pressure gradients (sugar loading at source, unloading at sink), xylem by transpiration pull from leaf surfaces. These are physically different phenomena — diffusion versus pressure flow — and they operate at different timescales and distances.

---

## Thickening as Vascular Consequence

The most important architectural implication of the vascularization model is one that reaches into a part of the sim that has nothing to do with transport: the thickening code in `StemNode`.

**In real plants, secondary thickening is vascular development.** Wood is old xylem. Every annual ring is a year's worth of new xylem tissue laid down by the vascular cambium — the lateral meristem that sits between the xylem and phloem of every vascular stem. The cambium divides when it is activated by auxin flowing through it. Each division adds a new xylem cell inward and a new phloem cell outward, making the stem wider. A stem does not "decide to thicken." It thickens because auxin is flowing through it, activating the cambium. The amount of thickening is directly proportional to the amount of sustained auxin flux — which is exactly what `structural_flow_bias` represents.

**The current model gets this backwards.** `StemNode::thicken()` is an independent growth step that spends sugar on radial growth, controlled by a `thickening_rate` genome parameter. Any stem with adequate sugar and auxin above `auxin_thickening_threshold` can thicken. The auxin threshold is a gate, not a driver — once the gate is open, the `thickening_rate` parameter determines how fast the stem widens, independent of how much auxin is actually flowing. A lateral branch with barely-threshold auxin and a main trunk with heavy auxin flux thicken at the same rate if they have the same sugar supply. This is structurally wrong and erases the hierarchy that canalization should be building.

**The new model:** thickening rate is driven by `structural_flow_bias`. The formula becomes:

```
delta_radius = cambium_responsiveness × structural_flow_bias × sugar_available_fraction
```

Where `cambium_responsiveness` is the new genome parameter (replacing `thickening_rate`), `structural_flow_bias` is the slow-ratchet canalization history on the parent-to-this-node connection (already maintained, never decays), and `sugar_available_fraction` is the fraction of sugar capacity currently filled — thickening still costs sugar, but the rate is set by canalization history, not by a fixed parameter.

The old `auxin_thickening_threshold` gate disappears. It is replaced by the natural minimum: `structural_flow_bias` starts at zero and only accumulates above `structural_threshold`. Zero vascular history = zero cambium activity = no thickening, even with abundant sugar. There is no separate gate to tune.

**This creates the self-reinforcing loop that builds main trunks.** A dominant branch carries heavy auxin flux → `structural_flow_bias` grows rapidly → cambium is highly active → stem widens → wider stem has larger pipe cross-section (π × r²) → carries even more auxin per tick → further bias accumulation → more thickening. A shaded side branch with weak auxin flux barely accumulates bias → cambium is barely active → stays thin → thin pipes → lower auxin throughput → low bias. The feedback separates the main axis from lateral branches without any explicit "this is the main axis" logic. It emerges from auxin flux history.

**Biological variation falls out of `cambium_responsiveness` naturally:**

- High `cambium_responsiveness` — fast thickening response per unit of vascular development. Models willow-like or oak-like rapid wood formation. A plant that builds trunk quickly.
- Low `cambium_responsiveness` — slow thickening even on high-flux paths. Models plants where secondary thickening is minimal but still present.
- Zero `cambium_responsiveness` — no secondary thickening ever. Models monocots (palms, grasses, bamboo) where the vascular bundles are scattered and the cambium is absent or vestigial. Stem width is set at creation and never changes. This is the genome-parameterized model of monocot vs dicot body plan.

The same genome parameter that controls how fast a plant builds wood also determines whether it builds wood at all.

**The maturation gate should also move.** Currently `has_vasculature()` uses an age-based gate (`cambium_maturation_ticks`) to decide whether a node participates in the vascular bulk transport pass. This is wrong for the same reason thickening-by-rate is wrong: a young internode on the main axis of an active meristem has been carrying auxin since day one and should be building vascular tissue fast. A young internode on a dormant lateral that activated two ticks ago should have weak vasculature regardless of age. The data we already track — `structural_flow_bias` — IS the vascular development state. Replace the age gate with a bias threshold: a node participates in the vascular pass when its `structural_flow_bias` (on the connection from its parent) exceeds a minimum conductance threshold. Until then it relies on local diffusion for everything, same as a leaf or meristem. Once it crosses the threshold, it is admitted to the vascular network, and its pipe capacity scales with its bias rather than just with radius. This makes vascular development path-dependent rather than time-dependent, which is biologically accurate.

---

## What Changes From Current State

Five specific changes bring the sim into alignment with this architecture. Changes 3 and 4 are the most structurally significant; changes 1 and 2 feed into them by producing the canalization data that drives the new thickening and maturation logic.

**1. Vascular Phase 2 weights by structural_flow_bias at junctions.**

The vascular two-phase pass currently distributes flow among children in proportion to their subtree demand. This is correct for a naive pipe model, but misses the structural history. A connection where auxin has flowed heavily over hundreds of ticks has developed more conductive tissue — it should carry more of the bulk flow budget, all else equal. At each junction in Phase 2, when distributing available flow among children, the demand-proportional share should be multiplied by `(1 + canalization_weight * structural_flow_bias)` for that child connection. This is the same weight used in local diffusion's `transport_with_children` — the same bias map, the same formula — just applied now to the vascular pass as well. A branch that has been a consistent auxin highway carries more sugar and water than a branch of equal size that has been dormant. This is what makes branch hierarchy emerge structurally rather than just topologically.

**2. Reduce local diffusion rates for auxin, GA, and stress.**

These chemicals should be genuinely local. In the current parameter defaults, diffusion rates are set broadly. The target behavior: a signal produced at the shoot apex should take on the order of 50+ ticks to traverse even a modest plant, and it should decay substantially along the way. A leaf's GA signal should affect only its own internode and the one above it — not travel to the canopy top or the seed. This is what makes the internode-local elongation model biologically correct: GA produced by a young leaf acts on the internode being elongated, not on the whole plant's elongation budget.

Auxin's slow local journey is what makes apical dominance work at all. If auxin diffuses fast enough to reach lateral buds three or four nodes below the apex in one or two ticks, apical dominance is uniform across the plant and the branching pattern loses spatial structure. The gradient should be steep: high auxin immediately below the apex, falling off within a few internodes, effectively zero at the branch points where lateral buds activate.

Concrete target rates: auxin diffusion rate reduced from current default toward 0.1–0.15; GA diffusion rate reduced toward 0.05 (GA is almost entirely local); stress diffusion rate reduced toward 0.1–0.15 (mechanical stress should propagate a few nodes but not long distances). These numbers need tuning against observed behavior, but the direction is clear: slower, not faster.

Slowing diffusion also means `structural_flow_bias` builds up more selectively — only connections on genuine auxin highways accumulate significant bias. This makes the canalization data a better signal for changes 3 and 4.

**3. Refactor thicken() to read structural_flow_bias.**

Remove the `thickening_rate` genome parameter and the `auxin_thickening_threshold` gate. Replace with `cambium_responsiveness`. Change `StemNode::thicken()` (and the equivalent in `RootNode`) to compute:

```
delta_radius = cambium_responsiveness × structural_flow_bias × sugar_available_fraction
```

The `structural_flow_bias` for a node is stored on its parent, keyed by the node's pointer — access it via `node->parent->structural_flow_bias[node]` (the same map that `update_canalization()` maintains). If the node has no parent (seed) or no bias entry yet, `structural_flow_bias` is zero and thickening is zero. The sugar cost of thickening remains proportional to `delta_radius` as before.

**4. Replace age-based vascular maturation gate with bias-based conductance.**

In `has_vasculature()` (currently `node.age >= cambium_maturation_ticks`), replace the age check with a bias threshold: the node participates in the vascular network when the `structural_flow_bias` on its parent-to-node connection exceeds a genome parameter `vascular_conductance_threshold`. Remove `cambium_maturation_ticks` and `root_cambium_maturation_ticks` from the genome.

Additionally, `pipe_capacity()` currently uses only `radius` to compute cross-sectional area. Augment it: `capacity = π × r² × conductance × (1 + bias_conductance_scale × structural_flow_bias)`. Well-developed vascular connections have more open lumen, more plasmodesmata density, more functional xylem vessels — the bias history should increase conductance beyond what radius alone predicts. Add `bias_conductance_scale` as a genome parameter (or start it as a fixed constant until tuning shows it matters).

**5. Make the architectural separation explicit in code comments.**

The two-system architecture should be legible in the code. `vascular.cpp` already has the right shape, but the comments should make explicit: (a) which chemicals go through vascular, (b) which go through local diffusion, (c) why, and (d) what the relationship between the two is via canalization. `node.cpp`'s `transport_with_children` should note that it handles only non-vascular chemicals plus last-mile delivery. `plant.cpp`'s `tick_tree` should annotate the call order and why vascular runs before the DFS walk. `stem_node.cpp`'s thicken function should explain that thickening is cambium activity driven by vascular history, not a separate rate. `CLAUDE.md`'s Chemical Transport Model section should summarize the two-system split.

---

## Biological Reference

**PIN proteins and polar auxin transport.** PIN (PIN-FORMED) proteins are efflux carriers that move auxin out of cells. They are asymmetrically distributed on specific faces of each cell — auxin exits preferentially from the basal face (toward roots) in most shoot tissues. This polarity is what makes auxin transport directional without requiring any active "decision" by each cell. The polarity itself is dynamic: PIN distribution is regulated by auxin concentration, creating a feedback loop. Where auxin flux is high, PIN proteins redistribute to reinforce that flux direction. This is the cellular basis of canalization.

**Vascular cambium and secondary thickening.** The vascular cambium is a thin cylinder of meristematic cells that sits between the primary xylem and primary phloem of every dicot stem. It divides periclinally (parallel to the stem surface) to produce new xylem cells inward and new phloem cells outward. Each division round widens the stem slightly. The signal that drives cambium division is auxin: auxin produced at the shoot apex moves basipetally, passes through the cambium, and activates it. More auxin flux = faster cambium division = faster wood production. This is why a dominant main stem thickens much faster than a lateral branch, even if both are the same age — the main stem has been carrying the full weight of apical auxin production, and the lateral branch carries only a fraction. The cambium doesn't know which branch is "main" — it only responds to how much auxin has been flowing through it.

**Transpiration pull.** Xylem water movement is driven by negative pressure at the top. Leaves transpire water vapor through stomata; this creates a tension that propagates down the continuous water column from leaf mesophyll through xylem vessels to root hairs. The whole column is under tension — xylem water is in a metastable state. This is why xylem is pressure-driven toward the transpiration sink (leaves) rather than toward a concentration equilibrium. Water moves upward not because roots push it but because leaves pull it. In the sim, the slight positive bias on water transport (`water_bias = 0.05`) approximates this directional pull without modeling the full tension gradient.

**Phloem osmotic pressure.** Phloem flow direction is reversed from xylem. Sugar is actively loaded into phloem at source tissues (leaves) and actively unloaded at sink tissues (roots, growing tips). Loading at source raises osmotic pressure there, which draws water in and creates turgor pressure. Unloading at sink lowers pressure. The pressure gradient drives bulk flow from source to sink. This is why phloem can flow bidirectionally depending on which ends are loading and unloading — there is no fixed up/down to phloem flow, unlike xylem.

**Sachs 1969 canalization hypothesis.** T. Sachs proposed that vascular strands form by a positive feedback: auxin moving through a cell induces that cell to become better at conducting auxin, which directs more auxin through it, which further increases conductance. Over time, flux concentrates into discrete channels — the future vascular strands. Undirected cells become parenchyma; cells in the early flux path develop into xylem and phloem. This explains why vascular tissue forms in continuous strands connecting auxin sources to auxin sinks, why leaf veins radiate from the midrib, and why wound-healing can redirect vascular strands around damage. The sim's `structural_flow_bias` (slow, never-decaying ratchet that grows with sustained auxin flux) is a direct implementation of the Sachs mechanism.

**The positive feedback loop in full:** auxin flows cell-to-cell through a tissue → flux above threshold at a connection → `structural_flow_bias` accumulates at that connection → cambium is activated proportionally → stem widens → wider stem has larger pipe cross-section → reinforced connection carries proportionally more of all bulk transport (sugar, water, cytokinin) → well-supplied tissue is more metabolically active → more active tissue produces more auxin → more auxin flows through the already-reinforced connection → more bias accumulation → more thickening. The loop makes early flux advantages self-amplifying, which is why main trunks become main trunks and why vascular patterns are robust and reproducible despite growing in a noisy environment.

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

**Step 3: Refactor thicken() and add cambium_responsiveness.**

In `genome.h`, add `cambium_responsiveness` (default TBD, probably in the range of 0.001–0.01 per tick) and `vascular_conductance_threshold` (minimum `structural_flow_bias` for vascular participation). Remove `thickening_rate`, `auxin_thickening_threshold`, `cambium_maturation_ticks`, and `root_cambium_maturation_ticks`.

In `StemNode::thicken()` (and `RootNode::thicken()`), replace the existing logic with:
```
float bias = parent ? parent->structural_flow_bias[this] : 0.0f;
float sugar_fraction = chemical(ChemicalID::Sugar) / sugar_cap(*this, g);
float delta = g.cambium_responsiveness * bias * sugar_fraction;
radius += delta;
chemical(ChemicalID::Sugar) -= delta * construction_cost_per_radius;
```

In `has_vasculature()` in `vascular.cpp`, replace the age check: a node is vascular when its parent's `structural_flow_bias` entry for it exceeds `vascular_conductance_threshold`. Seed node remains always vascular (it has no parent and is the root of the network).

Verification: grow a plant for 500+ ticks. The main stem should be visibly thicker than lateral branches of similar length. Lateral branches should show thickening proportional to their activity level, not a fixed rate. A branch that has never had a leaf and thus never carried much auxin should remain thin. Setting `cambium_responsiveness = 0.0` in the genome should produce a plant with no secondary thickening at all (stems stay at their initial radius forever) — this is the monocot/palm case.

**Step 4: Verify the self-reinforcing loop with a shading experiment.**

Grow a two-branch plant to the point where both branches have similar thickness. Then modify the lighting to shade one branch heavily (either via the `--color` visualization to identify a shaded configuration, or by temporarily zeroing `light_exposure` on one branch's leaves in the debugger). Observe over the next 100+ ticks:

- The lit branch: continues producing auxin → continues flowing through stem → `structural_flow_bias` keeps accumulating → cambium stays active → stem continues thickening.
- The shaded branch: leaves produce less auxin → auxin flux drops → `structural_flow_bias` stops growing (it never decays, but growth stops) → cambium activity drops to near zero → thickening halts.

This asymmetry is what distinguishes a tree with a dominant trunk from a shrub with many equal stems. The test verifies that the feedback loop works end-to-end: light → auxin production → auxin flux → canalization → thickening, and that interrupting the light input correctly halts the downstream thickening.

Verification: after the shading period, measure (or visually confirm) that the shaded branch's radius has not increased while the lit branch's radius has.

**Step 5: Annotate the two-system architecture in code.**

Add a block comment at the top of `vascular.cpp` explaining the vascular/diffusion split, which chemicals use each path, and the canalization bridge. Add a corresponding note in `node.cpp` near `transport_with_children` explaining that this function handles only non-vascular chemicals (Auxin, GA, Stress) plus last-mile delivery. Update `plant.cpp`'s `tick_tree` comment to explain why vascular runs before the DFS walk. Add a comment in `stem_node.cpp` (and `root_node.cpp`) near `thicken()` explaining that this is cambium activity driven by vascular history. Add a note in `CLAUDE.md`.

Verification: someone reading the code cold should understand from the comments alone why `vascular_transport` is a separate pre-DFS pass, why `thicken()` reads `structural_flow_bias` rather than using a rate parameter, and why the age-based vascular gate was replaced.

**Step 6: Integration test — canalization drives observable vascular hierarchy.**

Write a test (or extend an existing one) that: (1) grows a plant for enough ticks that canalization has built up structural bias on the main axis, (2) reads the `structural_flow_bias` values on the seed node for its two children (shoot-side and root-side), and (3) verifies that the bias values are non-zero and that the thicker of the two stems has higher bias. This is a regression test: it ensures that structural bias actually accumulates over time, that it correlates with radius (which it should, since both are driven by auxin flux), and that neither is inadvertently zeroed out by a future change.

Verification: test passes. Structural bias values are observable from outside the transport machinery, so this can be written as a black-box integration test on plant state after N ticks.

---

## Relationship to Other Milestone 2 Work

Vascularization interacts with water (milestone 2's keystone chemical) in two ways. First, the structural_flow_bias weights apply to water transport in the vascular pass, just as they will for sugar — water moves through the same xylem pipes that auxin canalized, so the same structural history governs how well-watered each branch is. This is already handled once Step 1 is complete, because the vascular pass runs all three chemicals (sugar, water, cytokinin) through the same `run_vascular` function. Second, and more subtle: turgor (water pressure) affects cambium activity in real plants. Well-hydrated tissue has active cambium; drought-stressed tissue slows down. Once water is a tracked resource, `cambium_responsiveness × structural_flow_bias × sugar_available_fraction` could gain a fourth term: `× water_available_fraction`. Leave this for after the base water model is stable.

The tissue library expansion (intercalary meristems, water-storage parenchyma) does not depend on vascularization. New tissue types that extend `StemNode` or `RootNode` will automatically participate in vascular transport and bias-driven thickening once they have accumulated enough `structural_flow_bias` on their parent connections. Types that do not extend those classes (leaf-like or meristem-like tissues) will continue to use local diffusion only, which is the correct behavior — the same rules, no special cases.

The monocot-body-plan case (`cambium_responsiveness = 0.0`) intersects with the tissue library because grasses and palms require not just zero secondary thickening but a different node width model: a grass internode is fixed-width from creation. This is already correct in the current model at `cambium_responsiveness = 0` since `delta_radius = 0` at every tick. The intercalary meristem tissue type can be added to the library independently of this parameter — it inherits the same thickening logic and simply produces zero secondary growth if the genome says so.
