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

## Flow Physics: Deficit, Surplus, and Competition

The vascular distribution model should be physics-driven, not allocation-driven. This distinction matters because it determines whether competitive dynamics between branches emerge naturally from the transport math or need to be engineered separately.

**The wrong framing: demand-proportional allocation.** The naive model — and one step of improvement over it — is to split available flow among children in proportion to their demand, perhaps weighted by canalization history. But this is still a central-allocation model: a scheduler tallies up who wants what and distributes accordingly. Real vascular flow doesn't work this way. There is no scheduler. There is only physics.

**The right framing: flow resistance determines distribution.** At each junction, resources flow toward children based on pipe conductance and resistance. Wider pipes (higher `structural_flow_bias`, larger radius) have lower resistance and receive more flow. A distant thin branch with high demand still gets less than a nearby thick branch with moderate demand, because the pipe geometry limits what can physically reach it. Demand is the pull; resistance is the bottleneck. Flow follows Hagen-Poiseuille: proportional to the pressure gradient and inversely proportional to resistance.

The practical formula for Phase 2 distribution at each junction:

```
child_conductance = pipe_capacity(child)   // π × r² × conductance × (1 + bias_scale × structural_flow_bias)
child_share = available × child_conductance / sum(conductance of all children)
child_share = min(child_share, child_demand)   // don't deliver more than demanded
// redistribute unclaimed capacity to remaining children
```

The weight is pipe conductance, not demand. Demand acts as a ceiling — a child that needs less than its conductance-weighted share takes only what it needs, and the remainder flows to the next most conductive connection. This is the physics of parallel resistors sharing a common pressure source.

The critical difference from the current plan's `demand × (1 + canalization_weight × structural_flow_bias)` formulation: demand is no longer the primary weight. A large-demand branch with thin pipes still gets less than a moderate-demand branch with thick pipes. Canalization history is already encoded in `pipe_capacity` via the `structural_flow_bias` augmentation — it doesn't need to be layered on top of a demand weight separately.

**Self-limiting root absorption under surplus.** The Phase 1 supply calculation should reflect hydraulic reality. Root water absorption is driven by the gradient between soil water potential and root cell water potential. When root cells are near capacity, the gradient collapses and absorption slows. This isn't a design choice — it's how osmotic uptake works. The supply calculation should be:

```
supply = absorption_rate × max(0, soil_water_potential - root_water_level / root_cap)
```

When `root_water_level` approaches `root_cap`, `root_water_level / root_cap → 1`, and if that exceeds the normalized soil potential, supply drops to zero. The plant stops absorbing not because of a guard condition but because the driving force disappears. Under drought (low `soil_water_potential`), absorption starts dropping before cells reach capacity because the soil gradient is weak. Under saturation (high `soil_water_potential`), absorption runs at near-maximum rate until cells are full. The same formula handles both cases.

**Transport bottlenecks create branch competition without a competition system.** The positive and negative spirals both follow from pipe-capacity-first distribution:

*Positive spiral (dominant branch):* active growing tip → high auxin production → sustained auxin flux through stem → `structural_flow_bias` accumulates → pipe capacity grows → more sugar and water arrive → tissue grows faster → more auxin produced → further flux → further bias → further capacity. Each cycle of this loop deepens the branch's advantage.

*Negative spiral (suppressed branch):* tip shaded or dormant → low auxin production → low flux → `structural_flow_bias` stagnates → pipe stays thin → thin pipe competes poorly for flow at parent junction → less sugar and water arrive → tip grows slower, may go dormant → even less auxin → bias stops growing → branch stays marginalized.

Both spirals emerge from the same pipe-conductance distribution math. There is no explicit "this branch wins" decision anywhere. The competition is just physics: thicker pipes carry more fluid. Building a separate branch-competition mechanism would be double-counting something the transport already does. The feedback loop between vascular distribution, canalization, and thickening IS the competition model.

**Local feedback replaces central allocation for leaves.** When a leaf can't get enough water through a thin or long path, stomata partially close, reducing transpiration demand. This is already modeled: `photosynthesis_rate *= clamp(water / water_cap, 0.2, 1.0)`. This scalar is the leaf's feedback signal — it reduces its own demand when supply is short. In the vascular distribution model, reduced demand means the leaf's ceiling is lower, which means its conductance-weighted share from the junction covers its full need and surplus stays available for better-connected competitors. The leaf gracefully "bows out" of competition without any coordination signal from the plant center.

This feedback loop is already in the code. It doesn't need a new mechanism. What it needs is to be connected to a distribution model that actually starves poorly-connected leaves rather than giving them a demand-proportional share regardless of pipe capacity. The stomatal model produces correct behavior only when paired with physics-correct distribution upstream.

**Surplus handling and the future storage feedback.** When supply exceeds demand — a well-watered plant in full sun — excess sugar accumulates. Eventually it should flow to storage: parenchyma cells in stems and roots act as expandable sugar sinks. Storage absorbs surplus, buffering both excess and deficiency. When storage is full, the parenchyma's demand signal drops to zero. Phloem sugar concentration rises. High phloem sugar is a signal to leaves that the whole-plant sugar economy is saturated, and production should throttle back. In real plants this is mediated by hexose signaling pathways; in the sim it can be approximated by: when storage nodes are above some fraction of capacity, leaves reduce photosynthesis proportionally.

This production throttle closes the surplus loop: soil/light → production → distribution → storage → saturation → production downregulation. The plant regulates its own output to match its capacity to use and store resources. This is a future step — implement the storage tissue type first, observe whether surplus accumulation is actually a problem in practice, and add the production throttle only if needed. But the loop should be noted in the plan so it isn't designed around.

---

## What Changes From Current State

Five specific changes bring the sim into alignment with this architecture. Changes 3 and 4 are the most structurally significant; changes 1 and 2 feed into them by producing the canalization data that drives the new thickening and maturation logic.

**1. Vascular Phase 2 uses resistance-based (conductance-first) distribution.**

The current pass distributes flow proportional to subtree demand, which is a central-allocation model. Replace with conductance-first distribution: each child's share is weighted by its `pipe_capacity` (already augmented by `structural_flow_bias` per change 4), with child demand as a ceiling rather than a weight. After each child takes up to its demand, unclaimed capacity is redistributed to remaining children by the same conductance weighting. This mirrors Hagen-Poiseuille flow at a junction: more goes through lower-resistance pipes, regardless of who "needs" it more.

This replaces the earlier plan of multiplying demand by `(1 + canalization_weight * structural_flow_bias)`. That was still demand-first. The new model is conductance-first: pipe geometry drives distribution, and demand only caps how much a node actually absorbs.

**2. Phase 1 supply (root absorption) uses gradient-based self-limiting.**

Root supply in the vascular Phase 1 pass should not be a fixed fraction of current water level. It should reflect the driving gradient: `supply = absorption_rate × max(0, soil_water_potential - root_water_level / root_cap)`. As roots fill, the gradient weakens and supply tapers off naturally. This makes surplus conditions self-correcting (full roots absorb less) and drought conditions physically accurate (weak soil gradient = low absorption rate regardless of root capacity).

**3. Reduce local diffusion rates for auxin, GA, and stress.**

These chemicals should be genuinely local. The target behavior: a signal produced at the shoot apex should take on the order of 50+ ticks to traverse even a modest plant, and it should decay substantially along the way. A leaf's GA signal should affect only its own internode and the one above it — not travel to the canopy top or the seed.

Concrete target rates: auxin diffusion rate reduced toward 0.1–0.15; GA diffusion rate reduced toward 0.05 (GA is almost entirely local); stress diffusion rate reduced toward 0.1–0.15. Slowing diffusion makes `structural_flow_bias` more selective — only genuine auxin highways accumulate significant bias, which makes the canalization data a cleaner signal for changes 4 and 5.

**4. Refactor thicken() to read structural_flow_bias.**

Remove the `thickening_rate` genome parameter and the `auxin_thickening_threshold` gate. Replace with `cambium_responsiveness`. Change `StemNode::thicken()` (and the equivalent in `RootNode`) to compute:

```
delta_radius = cambium_responsiveness × structural_flow_bias × sugar_available_fraction
```

The `structural_flow_bias` for a node is stored on its parent, keyed by the node's pointer. If the node has no parent (seed) or no bias entry yet, `structural_flow_bias` is zero and thickening is zero. The sugar cost of thickening remains proportional to `delta_radius` as before.

**5. Replace age-based vascular maturation gate with bias-based conductance.**

In `has_vasculature()` (currently `node.age >= cambium_maturation_ticks`), replace the age check with a bias threshold: the node participates in the vascular network when the `structural_flow_bias` on its parent-to-node connection exceeds `vascular_conductance_threshold`. Remove `cambium_maturation_ticks` and `root_cambium_maturation_ticks` from the genome.

Augment `pipe_capacity()` to scale with vascular development history: `capacity = π × r² × conductance × (1 + bias_conductance_scale × structural_flow_bias)`. Well-developed connections have more open lumen and more functional xylem vessels — the bias history should increase conductance beyond what radius alone predicts. This augmented capacity is what the resistance-based Phase 2 distribution (change 1) uses as its primary weight, completing the loop: canalization history → pipe capacity → distribution weighting → resource delivery → growth → more auxin → more canalization.

**6. Make the architectural separation explicit in code comments.**

`vascular.cpp` should explain the vascular/diffusion split, which chemicals use each path, and the canalization bridge. `node.cpp`'s `transport_with_children` should note that it handles only non-vascular chemicals plus last-mile delivery. `plant.cpp`'s `tick_tree` should annotate the call order. `stem_node.cpp`'s thicken function should explain that it is cambium activity driven by vascular history. `CLAUDE.md`'s Chemical Transport Model section should summarize the two-system split and the physics-driven distribution model.

---

## Biological Reference

**PIN proteins and polar auxin transport.** PIN (PIN-FORMED) proteins are efflux carriers that move auxin out of cells. They are asymmetrically distributed on specific faces of each cell — auxin exits preferentially from the basal face (toward roots) in most shoot tissues. This polarity is what makes auxin transport directional without requiring any active "decision" by each cell. The polarity itself is dynamic: PIN distribution is regulated by auxin concentration, creating a feedback loop. Where auxin flux is high, PIN proteins redistribute to reinforce that flux direction. This is the cellular basis of canalization.

**Vascular cambium and secondary thickening.** The vascular cambium is a thin cylinder of meristematic cells that sits between the primary xylem and primary phloem of every dicot stem. It divides periclinally (parallel to the stem surface) to produce new xylem cells inward and new phloem cells outward. Each division round widens the stem slightly. The signal that drives cambium division is auxin: auxin produced at the shoot apex moves basipetally, passes through the cambium, and activates it. More auxin flux = faster cambium division = faster wood production. This is why a dominant main stem thickens much faster than a lateral branch, even if both are the same age — the main stem has been carrying the full weight of apical auxin production, and the lateral branch carries only a fraction. The cambium doesn't know which branch is "main" — it only responds to how much auxin has been flowing through it.

**Transpiration pull and xylem hydraulics.** Xylem water movement is driven by negative pressure at the top. Leaves transpire water vapor through stomata; this creates a tension that propagates down the continuous water column from leaf mesophyll through xylem vessels to root hairs. The whole column is under tension — xylem water is in a metastable state. Water moves upward not because roots push it but because leaves pull it. Flow rate through a xylem vessel follows Hagen-Poiseuille: proportional to the fourth power of vessel radius and inversely proportional to length. This is why small increases in vessel radius produce large increases in conductance — doubling radius gives 16× the flow capacity. In the sim, `pipe_capacity = π × r² × conductance` approximates this (dropping the r⁴ sensitivity for stability), and `structural_flow_bias` augments conductance to capture the developmental component.

**Phloem osmotic pressure.** Phloem flow is driven by osmotic pressure gradients (Münch pressure flow). Sugar is actively loaded into phloem at source tissues (leaves), raising osmotic pressure there and drawing water in. Unloading at sinks lowers pressure. The gradient drives bulk flow from source to sink. Phloem can flow bidirectionally depending on where sources and sinks are — there is no fixed up/down to phloem flow, unlike xylem.

**Sachs 1969 canalization hypothesis.** T. Sachs proposed that vascular strands form by positive feedback: auxin moving through a cell induces that cell to become better at conducting auxin, which directs more auxin through it, which further increases conductance. Over time, flux concentrates into discrete channels — the future vascular strands. Undirected cells become parenchyma; cells in the early flux path develop into xylem and phloem. The sim's `structural_flow_bias` (slow, never-decaying ratchet that grows with sustained auxin flux) is a direct implementation of the Sachs mechanism.

**Stomatal regulation and demand-side feedback.** Stomata open to admit CO₂ for photosynthesis and simultaneously lose water vapor. Guard cells control aperture based on turgor: when water is scarce, guard cell turgor drops and stomata close. This reduces transpiration (protecting water reserves) at the cost of reduced photosynthesis. It also reduces the transpiration pull that drives xylem flow, which reduces cytokinin delivery to shoot tips, which can suppress growth. The stomatal response is local — each leaf responds to its own water status independently. There is no signal from a central "plant" deciding how much water each leaf gets. This is the model the sim should implement: conductance-weighted distribution determines how much water each leaf actually receives; stomatal aperture (and thus photosynthesis rate) is a local response to that delivery.

**The positive feedback loop in full:** auxin flows cell-to-cell through a tissue → flux above threshold at a connection → `structural_flow_bias` accumulates → cambium is activated proportionally → stem widens → wider stem has larger pipe cross-section → augmented pipe capacity → conductance-weighted distribution routes more resources to this connection → well-supplied tissue grows faster → produces more auxin → more flux through the same connection → further bias accumulation → further thickening. The loop makes early flux advantages self-amplifying. This is why main trunks become main trunks, why vascular patterns are reproducible despite developmental noise, and why shading a branch collapses its competitive position quickly — the feedback runs in reverse just as readily.

---

## Implementation Steps

These are ordered by dependency. Each step has a clear verification criterion so it can be confirmed working before the next begins.

**Step 1: Conductance-first Phase 2 distribution.**

In `vascular.cpp`, replace the demand-proportional distribution in Phase 2. At each junction, compute each child's conductance weight (the augmented `pipe_capacity` from change 5, even if that augmentation hasn't been implemented yet — use the plain `pipe_capacity` as a first pass). Distribute available flow proportionally to conductance weight. After allocating to each child, cap at that child's demand; redistribute unclaimed capacity to remaining children by the same conductance proportions. Iterate until all excess is absorbed or conductance-weighted capacity is exhausted.

Verification: run the realtime viewer with `--color sugar` and grow a plant for several hundred ticks. The primary axis should carry visibly more sugar than lateral branches of the same size and demand. Grow a plant where one branch is physically shorter and thus has lower path resistance to the seed junction — it should receive a larger share than an equally-demanding but more distant branch.

**Step 2: Gradient-based root absorption (Phase 1 supply).**

Update the Phase 1 supply computation for root nodes in `vascular.cpp`. Replace the current `n.chemical(chem_id) * 0.5f` supply fraction with a gradient-driven calculation:

```
float root_fill = n.chemical(ChemicalID::Water) / water_cap(n, g);
float supply = g.water_absorption_rate * std::max(0.0f, g.soil_water_potential - root_fill);
```

(Where `soil_water_potential` is a world parameter or normalized `WorldParams::soil_moisture`.) Full roots supply nothing regardless of soil conditions. Empty roots supply at full absorption rate. This makes surplus self-correcting without any explicit guard.

Verification: run a plant with high soil moisture. Root water levels should reach a steady state near capacity and then stop climbing — not overflow indefinitely. Reduce soil moisture to near zero; absorption should drop to near zero. Restore soil moisture; uptake should resume at full rate.

**Step 3: Tune local diffusion rates downward.**

In `chemical_registry.h`, reduce diffusion rates for Auxin (~0.1–0.15), GA (~0.05), and Stress (~0.1–0.15). Run the realtime viewer with `--color auxin`. The gradient from apex to base should be steep: high near shoot tips, dropping to low within 3–5 nodes. Lateral buds several nodes below an active apex should be in low-auxin territory and activate; buds immediately below the apex should remain suppressed.

Verification: branching pattern remains qualitatively similar (same genome). GA effects on elongation are internode-local — the elongating zone near a young leaf does not spread to internodes elsewhere in the canopy.

**Step 4: Refactor thicken() and add cambium_responsiveness.**

In `genome.h`, add `cambium_responsiveness` and `vascular_conductance_threshold`. Remove `thickening_rate`, `auxin_thickening_threshold`, `cambium_maturation_ticks`, and `root_cambium_maturation_ticks`.

In `StemNode::thicken()` (and `RootNode::thicken()`):
```cpp
float bias = parent ? parent->structural_flow_bias[this] : 0.0f;
float sugar_fraction = chemical(ChemicalID::Sugar) / sugar_cap(*this, g);
float delta = g.cambium_responsiveness * bias * sugar_fraction;
radius += delta;
chemical(ChemicalID::Sugar) -= delta * construction_cost_per_radius;
```

In `has_vasculature()`, replace the age check with: a node is vascular when its parent's `structural_flow_bias` entry exceeds `vascular_conductance_threshold`. Seed node remains always vascular.

Verification: grow for 500+ ticks. Main stem visibly thicker than equally-old lateral branches. Setting `cambium_responsiveness = 0.0` produces a plant that never thickens (monocot case).

**Step 5: Verify the competitive loop with a shading experiment.**

Grow a two-branch plant to similar thickness on both sides. Shade one branch (zero `light_exposure` on its leaves). Observe over 100+ ticks: shaded branch produces less auxin → flux drops → `structural_flow_bias` stagnation → cambium activity near zero → thickening halts → thinner pipe → less competitive in conductance-weighted distribution → less sugar and water → tip may go dormant. Lit branch continues all positive-feedback loops. The competitive divergence should be visible without any explicit "branch death" mechanism.

This test also validates that the stomatal feedback is wired correctly: shaded leaves should reduce photosynthesis and water demand (`water / water_cap` scalar), which frees resources for the lit branch. The reallocation is automatic — the conductance-weighted distribution redirects freed capacity to whoever has the conductance to use it.

Verification: after shading period, shaded branch radius unchanged; lit branch radius measurably increased. Lit branch water and sugar levels stable or rising; shaded branch levels dropping if demand doesn't fully reduce.

**Step 6: Integration test — canalization drives observable vascular hierarchy.**

Write a test that: (1) grows a plant for enough ticks that `structural_flow_bias` has accumulated on the main axis, (2) reads bias values on the seed node for its shoot-side and root-side children, (3) verifies that values are non-zero and correlate with the corresponding stem radius. This is a regression test ensuring that structural bias accumulates over time, correlates with thickening, and isn't inadvertently zeroed by a future change.

Verification: test passes. The stronger the bias → the thicker the stem. Both should grow together since they are driven by the same flux history.

**Step 7: Code annotation.**

Block comment at the top of `vascular.cpp`: vascular/diffusion split, which chemicals, why, canalization bridge, and the conductance-first distribution model. `node.cpp` near `transport_with_children`: non-vascular chemicals only, last-mile delivery. `plant.cpp` `tick_tree`: call order and why. `stem_node.cpp` `thicken()`: cambium activity driven by vascular history. `CLAUDE.md`: two-system split and physics-driven distribution.

---

## Relationship to Other Milestone 2 Work

**Water and the full hydraulic loop.** The conductance-first distribution model applies to water transport in the vascular pass, just as to sugar — the same `run_vascular` function handles all three chemicals. The gradient-based root absorption (Step 2) ties into the soil water model in [world-physics.md](../world-physics.md): `soil_water_potential` maps directly to the soil grid's local moisture level at the root's position. Once the soil model exists, root uptake becomes spatially heterogeneous: shallow roots in a dry zone take up less than deep roots reaching the water table, automatically, from the same gradient formula.

**Turgor and thickening.** Once water is a tracked resource, `cambium_responsiveness × structural_flow_bias × sugar_available_fraction` should gain a fourth term: `× water_available_fraction`. Well-hydrated tissue has active cambium; drought-stressed tissue slows. Leave this for after the base water model is stable.

**The surplus/storage/production throttle.** When storage nodes (parenchyma, future tissue type) are near capacity, their demand signal drops to zero. Phloem sugar concentration rises. High phloem sugar should downregulate photosynthesis — the production throttle that closes the surplus loop. Add this when the storage tissue type exists. Until then, surplus sugar accumulates in existing nodes, which is bounded by the capacity model already in place.

**Tissue library.** New tissue types that extend `StemNode` or `RootNode` automatically participate in bias-driven thickening and conductance-weighted vascular transport once they accumulate `structural_flow_bias` on their parent connections. Types that don't extend those classes (leaf-like, meristem-like) use local diffusion only. The intercalary meristem can be added independently — it inherits the same thickening logic and produces zero secondary growth if `cambium_responsiveness = 0.0` in the genome.

**The monocot case.** `cambium_responsiveness = 0.0` produces a plant with fixed-width stems from creation. Grasses and palms don't need a separate body plan — they're the same engine with this one parameter at zero. Combine with an intercalary meristem tissue type for grasses; with a fat initial radius for palms. The same distribution model and the same transport physics, just no wood formation.
