# Water Model

Water is the working fluid of every plant. It isn't just a thing leaves use when they photosynthesize — it is present at every living node, drives every growth event, carries every dissolved chemical, and is the primary resource that cactus-like strategies exist to hoard. Getting water right is the keystone of milestone 2.

---

## What Water Does in Plants

**Turgor pressure.** Every living plant cell is essentially a water balloon. Internal water pressure (turgor) keeps cells rigid and tissues stiff. When water is low, turgor drops — cells go flaccid, the plant wilts. This is mechanical collapse before metabolic failure. A thirsty plant droops before it starves. In sim terms: water level at a node should feed into structural rigidity, creating an alternative stress-resistance mechanism that doesn't require sugar-driven mass growth.

**Cell expansion and growth.** Cells don't grow by synthesizing wall material and filling it in. They grow by absorbing water and expanding against a loosened cell wall under turgor pressure. No water = no cell expansion = no elongation growth, regardless of sugar availability. The reason cacti grow slowly is not primarily metabolic — it's that every millimeter of new stem length requires water uptake to drive cell expansion. In sim terms: water level should gate elongation growth, independent of the sugar gate that already exists.

**Transport medium.** Xylem carries water from roots to leaves, but that water is also carrying dissolved minerals (nitrogen, phosphorus, etc.) upward. Phloem carries sugar dissolved in water downward. In the current sim, transport is driven by chemical concentration gradients and radius-bottlenecked throughput — that's a reasonable phloem analog. Xylem doesn't exist as a distinct transport. Adding water as a resource implicitly creates a basis for xylem: the pressure-driven flow of water upward from root to leaf, which would carry cytokinin with it (see hormone coupling below).

**Photosynthesis chemistry.** Water is split during the light reactions to donate electrons and release O2. This is a direct chemical requirement, not a side effect. In practice this rarely limits photosynthesis except under severe dehydration, but it means extremely water-stressed nodes should see photosynthesis drop even if stomata were somehow open.

**Cooling.** Transpiration off leaf surfaces is evaporative cooling. A leaf in full sun without transpiration overheats quickly. This matters if temperature modeling is added later.

**Solvent for everything.** Hormones, sugars, ions — all dissolved in water. This is background plumbing rather than something to model explicitly.

---

## Water as a Stress Resistance Lever

Currently the only way plants counter mechanical stress (droop, breakage) is to grow more mass: thicker stems cost sugar. Water gives a second lever. A well-hydrated stem is stiffer under the same mass load than a dehydrated one — turgor adds rigidity without adding dry weight. This means:

- Plants can resist moderate stress by being well-watered rather than by being thick.
- During drought, the same stem becomes more vulnerable to wind/gravity even if it hasn't lost any dry mass.
- Cactus stems that are fat with water are actually quite strong structurally; when drought-stressed they become softer and may droop.

In the sim, this probably translates to: a `water_stiffness_contribution` parameter that scales how much node water level adds to the effective rigidity used in the stress/droop calculation. At full hydration, it supplements structural mass; at low hydration, the node relies on its dry mass alone.

---

## Hormone Transport Coupling

Different hormones move by different mechanisms, and water matters differently to each:

**Cytokinin — directly water-dependent.** Cytokinin is synthesized in root tips and moves upward through the xylem stream. That stream is driven by transpiration: leaves lose water vapor, creating a negative pressure that pulls water (and cytokinin dissolved in it) up from the roots. Low water = low transpiration = low xylem flow = cytokinin accumulates in roots and doesn't reach shoots. This is biologically meaningful — cytokinin delivery to shoots signals "roots are functional and watered." A drought-stressed plant loses cytokinin signal to its shoot tips, which can cause growth suppression and branch dieback.

In sim terms: cytokinin transport rate scales with plant water status. When water is plentiful, cytokinin reaches the canopy efficiently. Under drought, it stalls near the roots.

**Auxin — indirectly water-dependent.** Auxin moves by polar active transport — cell-to-cell via PIN efflux proteins on specific cell faces. This is metabolically driven and directional, not bulk flow. Water status affects auxin transport indirectly: low turgor slows cell metabolism, which slows the active pumping that drives auxin polarity. The effect is real but not 1:1 with water flux. In sim terms: mild coupling, not the direct scaling that cytokinin gets.

**Abscisic acid (ABA) — the drought hormone, inverse relationship.** Low water triggers ABA synthesis in roots and leaves. ABA moves from roots to leaves via xylem (carried with water, somewhat paradoxically — there's still some residual xylem flow even under drought). ABA causes stomata to close, reducing further water loss. So ABA transport rate increases when water drops — opposite direction from the others. In sim terms: ABA isn't currently modeled as a distinct chemical. When it is added, its transport should invert the water-coupling rule: drought increases ABA signal rather than suppressing it.

**Ethylene — water-independent.** Ethylene is a gas; it diffuses through air spaces in tissue, not through the vascular stream. Water level doesn't meaningfully affect its diffusion. Keep ethylene decoupled from the water model.

**Gibberellins — water-dependent.** GAs move in both xylem and phloem, largely passively with mass flow. Similar to cytokinin but less directional. Water-dependent transport coupling applies.

**Practical implication for sim:** The cleanest rule is "most hormone transport rates scale mildly with plant water status; cytokinin scales directly; ethylene is exempt." This gives water a gameplay-meaningful feedback loop: drought doesn't just slow growth mechanically (via turgor), it disrupts the signaling network that coordinates growth. Under drought, apical dominance may break down as cytokinin fails to reach shoot tips; branch coordination degrades. This is biologically real and creates interesting structural failure modes.

**Planned refactor — cytokinin production source:** The current sim produces cytokinin in leaf nodes, proportional to photosynthesis output (`leaf.cpp:77`). This was a pragmatic choice and works for the current tree model. It becomes wrong once water is modeled. In real plants, cytokinin is produced by root tips in proportion to their water and nutrient access, then transported acropetally to tell shoots "go ahead and grow, the roots have what they need." That signal has actual information content only if its source is at the roots. When water lands, migrate cytokinin production from `leaf.cpp` to `root_apical.cpp`, proportional to root water uptake rate. The acropetal bias already in place (`cytokinin_bias = 0.1`) will then carry the right signal in the right direction.

**Recommended sequencing:** Don't couple hormone transport to water on day one. Get water as a resource working (root uptake, xylem transport, leaf loss, turgor effects on growth and stress) before layering the hormone coupling. Each coupling adds a tuning surface, and you want the base to be stable first.

---

## Implementation Sketch

Water should be a new `ChemicalID` tracked at every node alongside sugar, auxin, cytokinin, and the rest. Unlike sugar (which is a plant-internal resource) water has an external source: the soil.

**Root uptake.** Root nodes absorb water from the soil grid each tick. Uptake rate proportional to root surface area — approximated as `radius × length` (proportional to the lateral surface of a cylinder; the 2π constant is absorbed into the genome rate parameter) — and local soil water level. This is the inflow side of the water economy. Using surface area rather than radius alone is biologically important: a long thin root and a short thick root of equal dry mass have different water uptake. The thin-long root has more absorptive surface per unit mass, making it more water-efficient. This is what allows thin fibrous root architecture to be an evolutionarily viable strategy — without surface area scaling, evolution would always push toward uniformly thick roots.

**Xylem transport.** Water moves upward from root to leaf via pressure-driven flow. Unlike phloem (concentration-gradient-driven), xylem is driven by the negative pressure created by transpiration at the leaf surface. In sim terms: the simplest model is to treat it as a separate upward-biased transport pass, distinct from the unified hormone/sugar transport. Alternatively: add water to the unified transport with a strongly positive bias (acropetal) and high transport scale. The key difference from other chemicals is that water flows *toward* the transpiration sink, not toward a concentration equilibrium.

**Leaf transpiration loss.** Each tick a leaf photosynthesize, it loses water proportional to stomatal aperture × leaf area × vapor pressure deficit (or just a constant rate as a first pass). This is the outflow side. It's also what creates the xylem pull that drives upward flow.

**Turgor feedback.** Node water level (relative to some per-node capacity) determines turgor. Turgor affects: elongation growth rate (gates cell expansion), structural rigidity contribution (counteracts droop), and eventually stomatal aperture (stomata close when turgor drops — a self-protective response).

**Soil grid.** See [world-physics.md](world-physics.md) for the soil water model structure. From the plant's perspective: roots query the soil grid for local water availability and deplete it proportionally to uptake.

---

---

## Step 1: Minimum Viable Water

The incremental first step: add water as a tracked chemical with producers, consumers, and diffusion — and nothing else. No turgor coupling, no hormone coupling, no growth-rate coupling. Just observe water moving through the plant and tune until it feels right before adding any dependent systems.

**What it touches:**

`chemical.h` — add `Water` to the `ChemicalID` enum. The existing registry in `chemical_registry.h` has an `all_chemical_ids` array that needs to grow from 6 to 7, and a `diffusion_params()` function that needs a Water entry (diffusion rate, decay rate, bias, transport base and scale).

`genome.h` — add per-chemical genome fields: `water_diffusion_rate`, `water_decay_rate`, `water_bias`, `water_base_transport`, `water_transport_scale`. Same pattern as the existing cytokinin or auxin fields. Add defaults in `default_genome()`.

`root_node.cpp` — inside `RootNode::tissue_tick()`, add water production proportional to root surface area and a `water_uptake_rate` genome parameter. Both `radius` and `length` are available on the node. Pattern: `chemical(ChemicalID::Water) += radius * length * g.water_uptake_rate`. Surface area (not radius alone) is what allows thin fibrous roots to be evolutionarily competitive with thick roots — a long thin root and a short thick root of the same mass have different uptake. When a soil model exists this will scale with local soil water level; for step 1 it's a constant-rate baseline.

`leaf_node.cpp` — inside `LeafNode::tissue_tick()`, deduct water proportional to photosynthesis output and a `transpiration_rate` genome parameter. The cost is proportional to leaf area × light exposure × the rate. Pattern: deduct from `chemical(ChemicalID::Water)` at the same place sugar is produced, proportional to sugar produced.

`stem_node.cpp` — optionally, add a `lose_water()` step inside `StemNode::tissue_tick()` that deducts water proportional to stem surface area × `(1 - stem_cuticle_thickness)`. Surface area here is also `radius * length` (same approximation as root uptake). Start with cuticle at 1.0 (no loss) and tune down from there so it doesn't immediately drain stems dry. The same surface-area logic applies: long thin stems lose more water per unit mass than short thick ones, so cuticle becomes especially important for tall, thin-stemmed plants.

`leaf_node.cpp` — water deducted proportional to sugar produced (photosynthesis output) × `transpiration_rate`. Sugar production is already proportional to `leaf_area × light_exposure`, so this is implicitly surface-area-based — no correction needed here.

**Flow direction note:** the existing transport model is concentration-gradient-based with an optional directional bias. If root nodes produce water (high level) and leaf nodes consume it (low level), the gradient naturally drives flow from roots outward through the tree graph toward leaves — which is spatially upward in the shoot. No bias is needed to get the direction right; the gradient does it. A mild positive bias (0.1–0.2) can be added later if water doesn't distribute fast enough to leaf tips, but start without it.

**What "tuned" looks like:**

- Under normal conditions: water reaches steady state throughout the plant. Root nodes have highest water, leaves have lowest, stems gradient in between. No node permanently drains to zero while roots are active.
- Under simulated drought (set root production rate to 0): plant water levels drop across all nodes over ~50–100 ticks. Leaves drop first (they're consuming). The gradient flattens as leaves and stems equilibrate at low levels.
- Under re-watering (restore production): plant recovers over ~30–50 ticks, roots fill first, then stems, then leaves.

**Debug visualization:** add `water` as a new `--color` mode in the realtime viewer, alongside the existing auxin/cytokinin/sugar modes. Per-node color from blue (low water) to cyan (high water) gives an immediate spatial readout of the water distribution. This is the primary tuning tool.

**Step 1 boundary — defer these until step 2+:**

- No turgor coupling (water level does not affect stem rigidity or droop yet)
- No growth-rate coupling (water level does not gate elongation growth yet)
- No hormone-transport coupling (cytokinin transport does not scale with water yet)
- No soil grid (roots produce at a constant rate, not from a spatial resource)
- No ABA / drought signaling

**Call step 1 done when:** water distributes visibly and sensibly through the plant in the color mode visualization, responds correctly to production shutdown (simulated drought) and restoration, and doesn't break the existing sugar economy or hormone transport.

---

## Future Layers

Once the base model is stable:

- **Stomatal aperture as a state variable.** Currently photosynthesis runs at full rate whenever light is available. A stomatal model adds: aperture (0–1) controlled by water status, light, and time-of-day. CAM photosynthesis = stomata open at night only, storing CO2 as malic acid. C3 = stomata open proportionally to light. This makes the photosynthesis/water tradeoff explicit.
- **Vapor pressure deficit.** Transpiration rate scales with how dry the air is. Adds a temperature/humidity dimension to the water economy.
- **Hydraulic failure.** Under severe drought, xylem cavitates — air bubbles block water columns. A threshold below which conductance drops catastrophically, accelerating drought death in a realistic way.
