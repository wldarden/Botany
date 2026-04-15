# World Physics

The current sim has gravity (droop, stress, breakage), light exposure per leaf, and a sugar economy. That's a real foundation. What it lacks is any form of resource scarcity or competition *between* plants — which means there's no selective pressure for any particular strategy. A cactus and a tree planted in the same sim today would just grow differently with no ecological reason why one should outcompete the other in any given environment.

This doc describes what needs to be added to make the world a real sandbox, with a clear minimum viable set for the cactus↔tree divergence.

---

## What Exists

- **Gravity** — stems droop under mass load, break if stress exceeds threshold. Wood density and flexibility are genome parameters.
- **Light exposure per leaf** — each LEAF node has a `light_exposure` value. Currently this is based on the node's position in the world and the world's `light_level`, but there's no shading of one plant by another.
- **Sugar economy** — leaves produce sugar, all nodes consume it for maintenance and growth. Starvation leads to node death. This is a real internal economy; what's missing is that the *inputs* to this economy (light, water) aren't yet subject to competition.
- **Crowding via ethylene** — the `ethylene_crowding_radius` parameter creates local density pressure. Verify whether this currently affects neighboring plants' nodes or only nodes within the same plant. Cross-plant crowding is what creates real space competition.

---

## What Needs to Be Added

### Soil water grid

A 2D grid at ground level (y = 0). Each cell holds a water level. The grid evolves each tick:

- **Rain events** — periodic or stochastic additions to grid cells. Uniform rain adds to all cells; localized rain creates moisture gradients. Drought periods are simply rain events that don't happen.
- **Root depletion** — root nodes that overlap a grid cell deplete water proportional to their uptake rate and local soil water. This is the competition mechanism: two plants with overlapping root systems compete for the same pool.
- **Lateral diffusion** — water diffuses slowly between adjacent cells. Prevents roots from creating perfectly sharp depletion zones; creates a more realistic spreading gradient.
- **Evaporation** — bare soil loses water each tick at a constant rate. Shaded soil (under a canopy) evaporates more slowly — a feedback loop where trees reduce evaporation in their own footprint.

Grid resolution: probably 0.5–1.0 dm per cell. Coarser is faster; finer is more accurate root competition. Can start coarse and refine.

### Evapotranspiration

Currently photosynthesis is free — a leaf produces sugar proportional to light × area, no water cost. Adding evapotranspiration closes the key tradeoff: large leaf area = high sugar production but high water cost. Plants in dry environments pay a steep price for big leaves.

Implementation: each tick a LEAF node photosynthesize, it draws water from its node pool (or the nearest root-connected pool) proportional to stomatal aperture × leaf area × a `transpiration_rate` genome parameter. If the node has insufficient water, photosynthesis is throttled.

This single addition creates the cactus strategy's raison d'être: minimize transpiring surface to survive drought, accept lower sugar production, live on stored water between rains.

### Light occlusion between plants

Currently each plant's leaf light exposure is independent. Adding occlusion: for each LEAF node, cast a ray toward the light source (or simply upward in the vertical-sun approximation), accumulate shadow contribution from all nodes the ray passes through, reduce `light_exposure` proportionally.

This doesn't need to be ray-traced per tick. A coarse approximation: project all nodes onto a 2D top-down grid at the start of each tick, assign shadow density to each cell, then sample the cell above each leaf node to compute its light reduction. A cell with many nodes above it is deeply shaded.

With occlusion, tall trees cast shade on shorter competitors. Height becomes a real competitive strategy. Without it, there's no selective pressure to grow tall vs grow wide.

### Spatial environmental variation

For cactus and tree strategies to diverge, the world needs to have distinct ecological niches. The simplest version: a moisture gradient across the x-axis. The left side of the world receives frequent rain; the right side rarely does. Plants near the center experience intermediate conditions.

This means a plant's position in the world now matters to its strategy. A tree genome placed in the dry zone starves or fails to compete. A cactus genome placed in the wet zone gets outcompeted by faster-growing trees. Neither is just "wrong" — they're the wrong tool for that environment.

Gradient steepness and world size determine how sharply the strategies separate. A very steep gradient (desert vs rainforest next to each other) exaggerates divergence; a shallow gradient allows more intermediate strategies.

---

## Future Additions

These aren't needed for the cactus↔tree divergence but become important for broader coevolution:

**Temperature.** Affects evapotranspiration rate (hot air = more water loss), freeze/thaw stress (cold = dormancy pressure, ice crystal damage), and metabolic rates generally. Important for modeling boreal vs tropical vs alpine vs desert. For cactus specifically: desert temperature extremes (hot days, cold nights) matter for CAM photosynthesis timing.

**Seasonality.** Light level and temperature oscillating over a period creates dormancy pressure and the deciduous leaf cycle. Currently plants optimize for steady-state conditions. Seasonality rewards plants that can shed and regrow leaves efficiently, or that store energy and survive winter.

**Wind.** A directional mechanical load varying with height. Creates selection pressure for stem thickness (especially at the base), low center of gravity, and flexible wood vs stiff wood. The stress/droop mechanics already exist — wind is just an additional force input.

**Nutrient availability.** A second soil resource (nitrogen, phosphorus) separate from water. Low-nutrient soils favor plants that grow slowly and invest in long-lived tissue (many desert and tundra plants). High-nutrient soils favor fast-growing strategies. Creates a second axis of variation beyond moisture.

---

## Minimum Viable Set for Cactus↔Tree Divergence

Four things, in rough dependency order:

1. **Soil water grid** — makes water a real soil resource with spatial variation.
2. **Root water uptake** — connects root nodes to the soil grid (also requires water as a chemical in plants; see [water-model.md](water-model.md)).
3. **Evapotranspiration** — closes the leaf-area/water-loss tradeoff; makes cactus leaf reduction meaningful.
4. **Light occlusion** — makes height competitively valuable; makes trees worth being.
5. **Moisture gradient** — gives each strategy an ecological niche where it wins.

With these five things in place, a mixed population of plant genomes in a world with a moisture gradient should, over evolutionary time, develop distinct water-use strategies in the dry zone vs the wet zone. Whether those strategies look recognizably cactus-like and tree-like depends on whether the tissue library supports water-storage parenchyma and intercalary meristems — which is why these two tracks are parallel and mutually reinforcing.
