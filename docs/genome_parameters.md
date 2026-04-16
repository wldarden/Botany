# Genome Parameters Reference

All tunable plant parameters live in [`src/engine/genome.h`](../src/engine/genome.h). Units: distance = **dm** (10 cm), time = **ticks** (1 tick = 1 hour), mass = **grams glucose** (sugar) or **grams** (wood). See [`world_params.h`](../src/engine/world_params.h) for non-genetic simulation constants.

Evolution ranges come from [`genome_bridge.cpp`](../src/evolution/genome_bridge.cpp) — the min/max a gene can reach through mutation.

---

## Auxin

Auxin is produced by shoot apical meristems and flows basipetally (downward toward roots). It inhibits lateral bud activation — high auxin in a parent node keeps axillary buds dormant (apical dominance).

| Parameter | Default | Min | Max | Units | Description |
|-----------|---------|-----|-----|-------|-------------|
| [`auxin_production_rate`](../src/engine/genome.h) | 0.15 | 0.01 | 2.0 | signal/tick | Auxin added per tick by each active shoot apical meristem. Higher = stronger apical dominance, fewer side branches. |
| [`auxin_diffusion_rate`](../src/engine/genome.h) | 0.3 | 0.01 | 1.0 | fraction/tick | Fraction of auxin exchanged with parent per tick. Combined with `auxin_directional_bias` (-0.9) for 90% basipetal push, 10% gradient. |
| [`auxin_decay_rate`](../src/engine/genome.h) | 0.15 | 0.001 | 0.5 | fraction/tick | Fraction of auxin lost per tick. Higher = faster signal clearance, more responsive branching. |
| [`auxin_threshold`](../src/engine/genome.h) | 0.15 | 0.01 | 1.0 | signal | Auxin level in a **parent** node below which a shoot axillary bud activates. Lower = fewer branches (bud needs very low auxin to wake). Higher = more branches (buds activate even with moderate auxin). |

**Key interactions:** Shoot axillary buds check `parent->auxin` against this threshold. When the apical meristem is removed or distant, auxin decays below threshold and lateral buds break.

---

## Cytokinin

Cytokinin is produced by root apical meristems and flows acropetally (upward toward shoots). It gates all growth — nodes need cytokinin above a threshold to grow at full speed.

| Parameter | Default | Min | Max | Units | Description |
|-----------|---------|-----|-----|-------|-------------|
| [`cytokinin_production_rate`](../src/engine/genome.h) | 5.0 | 0.01 | 2.0 | signal/g sugar | Cytokinin produced per gram of sugar produced by leaves, generated at root apical meristems. Links root growth to shoot productivity. |
| [`cytokinin_diffusion_rate`](../src/engine/genome.h) | 0.3 | 0.01 | 1.0 | fraction/tick | Fraction exchanged with parent per tick. Combined with `cytokinin_directional_bias` (+0.9) for 90% acropetal pull, 10% gradient. |
| [`cytokinin_decay_rate`](../src/engine/genome.h) | 0.05 | 0.001 | 0.5 | fraction/tick | Fraction lost per tick. Lower than auxin decay — cytokinin is a slower, more persistent signal. |
| [`cytokinin_threshold`](../src/engine/genome.h) | 0.15 | 0.01 | 1.0 | signal | Cytokinin in a **parent** node below which a root axillary bud activates. Mirrors auxin threshold for root branching. |
| [`cytokinin_growth_threshold`](../src/engine/genome.h) | 0.1 | 0.01 | 1.0 | signal | Cytokinin level needed for full-speed growth at any node. Growth fraction = `min(1, cytokinin / threshold)`. Gates shoot extension, leaf growth, and elongation. |

**Key interactions:** Root-to-shoot signal. Plants with poor root systems have low cytokinin, causing all above-ground growth to slow. This creates a natural shoot:root balance.

---

## Shoot Growth

Controls how the shoot apical meristem extends the stem and creates new nodes.

| Parameter | Default | Min | Max | Units | Description |
|-----------|---------|-----|-----|-------|-------------|
| [`growth_rate`](../src/engine/genome.h) | 0.008 | 0.001 | 0.05 | dm/tick | Shoot tip extension speed. Default ~2 cm/day. The meristem advances this far each tick (scaled by sugar availability and cytokinin). |
| [`shoot_plastochron`](../src/engine/genome.h) | 24 | 6 | 168 | ticks | Ticks between new node creation at the shoot apex. At default, one new internode + leaf + axillary bud per day. Lower = denser nodes, more leaves, higher sugar cost. |
| [`branch_angle`](../src/engine/genome.h) | 0.785 | 0.05 | 1.57 | radians | Angle of lateral branches from the parent stem direction. Default ~45 deg. Affects crown shape — narrow angles give columnar form, wide angles give spreading canopy. |
| [`thickening_rate`](../src/engine/genome.h) | 0.00004 | 0.00001 | 0.001 | dm/tick | Radial growth speed of stem internodes. Default ~3.5 mm radius/year. Thicker stems support more mass but cost more sugar for maintenance. |
| [`internode_elongation_rate`](../src/engine/genome.h) | 0.004 | 0.0005 | 0.02 | dm/tick | Intercalary (mid-stem) stretch rate for young internodes. Distinct from tip growth — this is the internode lengthening after creation. Boosted by gibberellin. |
| [`max_internode_length`](../src/engine/genome.h) | 1.0 | 0.05 | 3.0 | dm | Maximum internode length (elongation target). Default 10 cm. Boosted by gibberellin via `ga_length_sensitivity`. |
| [`internode_maturation_ticks`](../src/engine/genome.h) | 72 | 12 | 500 | ticks | Hours until an internode stops elongating. Default 3 days. After this, the internode is "mature" and only thickens, no longer stretches. |

**Key interactions:** Shoot growth consumes sugar at rate `growth_rate * sugar_cost_shoot_growth` (world param). Growth fraction is gated by both sugar availability and cytokinin level.

---

## Root Growth

Controls root apical meristem behavior. Roots grow downward via gravitropism and create the underground network that produces cytokinin.

| Parameter | Default | Min | Max | Units | Description |
|-----------|---------|-----|-----|-------|-------------|
| [`root_growth_rate`](../src/engine/genome.h) | 0.004 | 0.001 | 0.05 | dm/tick | Root tip extension speed. Default ~1 cm/day (half of shoot rate). |
| [`root_plastochron`](../src/engine/genome.h) | 24 | 6 | 168 | ticks | Ticks between root node creation. Controls root branching density. |
| [`root_branch_angle`](../src/engine/genome.h) | 0.35 | 0.05 | 1.57 | radians | Angle of lateral roots from parent root. Default ~20 deg (narrower than shoot branches). |
| [`root_internode_elongation_rate`](../src/engine/genome.h) | 0.002 | 0.0005 | 0.02 | dm/tick | Intercalary stretch for young root internodes. |
| [`root_internode_maturation_ticks`](../src/engine/genome.h) | 48 | 12 | 500 | ticks | Hours until root internode stops elongating. Default 2 days. |
| [`root_gravitropism_strength`](../src/engine/genome.h) | 0.20 | 0.1 | 10.0 | dimensionless | How strongly root tips turn downward near the soil surface. Applied as a blend toward the down vector, scaled by depth exposure. |
| [`root_gravitropism_depth`](../src/engine/genome.h) | 0.5 | 0.1 | 5.0 | dm | Depth at which gravitropism correction ceases. Above this depth, roots are pulled downward; below it, they grow freely. Default 5 cm. |

**Key interactions:** Root meristems are hard-capped at 100 per plant (`Plant::max_root_meristems`). Root apicals produce cytokinin, linking root mass to shoot growth capacity.

---

## Geometry

Physical dimensions of plant organs — leaf size, stem radius, and growth perturbation.

| Parameter | Default | Min | Max | Units | Description |
|-----------|---------|-----|-----|-------|-------------|
| [`max_leaf_size`](../src/engine/genome.h) | 1.5 | 0.05 | 3.0 | dm | Maximum leaf side-length at maturity. Default 15 cm. Larger leaves capture more light but cost more sugar and are heavier. |
| [`leaf_growth_rate`](../src/engine/genome.h) | 0.005 | 0.0001 | 0.01 | dm/tick | Speed at which leaves grow from bud to max size. Default reaches full size in ~12 days. |
| [`leaf_bud_size`](../src/engine/genome.h) | 0.02 | 0.005 | 0.1 | dm | Initial leaf size at creation. Default 2 mm. |
| [`leaf_petiole_length`](../src/engine/genome.h) | 0.5 | 0.1 | 2.0 | dm | Stalk length holding leaf away from stem. Default 5 cm. Longer petioles spread leaves further, reducing self-shading. |
| [`initial_radius`](../src/engine/genome.h) | 0.05 | 0.01 | 0.2 | dm | Stem radius at creation. Default 5 mm. Determines initial structural capacity and sugar transport conductance. |
| [`root_initial_radius`](../src/engine/genome.h) | 0.025 | 0.005 | 0.1 | dm | Root radius at creation. Default 2.5 mm. |
| [`tip_offset`](../src/engine/genome.h) | 0.01 | 0.001 | 0.1 | dm | Small forward offset when chaining new nodes onto meristems. Prevents zero-length edges. |
| [`growth_noise`](../src/engine/genome.h) | 0.26 | 0.01 | 0.8 | radians | Maximum random angular perturbation per new segment. Default ~15 deg. Creates natural-looking irregular branching. |

---

## Sugar Economy

Photosynthesis, transport, maintenance costs, and storage. Sugar is the universal currency — all growth and survival depends on it.

| Parameter | Default | Min | Max | Units | Description |
|-----------|---------|-----|-----|-------|-------------|
| [`sugar_production_rate`](../src/engine/genome.h) | 0.02 | 0.001 | 0.1 | g/(dm^2 leaf * tick) | Glucose produced per unit leaf area per tick at full sunlight. Actual production is `rate * leaf_area * light_exposure * world.light_level`. |
| [`sugar_diffusion_rate`](../src/engine/genome.h) | 0.8 | 0.01 | 1.0 | fraction/tick | Base diffusion rate. Actual conductance is scaled by `min(radius_a, radius_b)` between nodes, so thin stems throttle sugar flow. |
| [`sugar_maintenance_leaf`](../src/engine/genome.h) | 0.001 | 0.001 | 0.1 | g/(dm^2 * tick) | Glucose consumed per unit leaf area per tick. Default ~5% of max production — leaves are net-positive in good light. |
| [`sugar_maintenance_stem`](../src/engine/genome.h) | 0.028 | 0.001 | 0.2 | g/(dm^3 * tick) | Glucose consumed per unit stem volume per tick. Volume = pi * r^2 * length. Thick, long stems are expensive. |
| [`sugar_maintenance_root`](../src/engine/genome.h) | 0.005 | 0.01 | 0.5 | g/(dm^3 * tick) | Glucose consumed per unit root volume per tick. |
| [`sugar_maintenance_meristem`](../src/engine/genome.h) | 0.0001 | 0.0001 | 0.01 | g/tick | Glucose consumed per active meristem tip per tick. Cheap individually, but many active tips add up. |
| [`seed_sugar`](../src/engine/genome.h) | 48.0 | 5.0 | 200.0 | g glucose | Initial sugar reserves in the seed. Default supports ~15 days of heterotrophic growth before leaves must take over. |
| [`sugar_storage_density_wood`](../src/engine/genome.h) | 500.0 | 50.0 | 2000.0 | g/dm^3 | Max sugar a stem/root node can hold per unit volume. High cap so stems can act as transport conduits. |
| [`sugar_storage_density_leaf`](../src/engine/genome.h) | 2.0 | 0.05 | 5.0 | g/dm^2 | Max sugar a leaf can hold per unit area. Provides export buffer. |
| [`sugar_cap_minimum`](../src/engine/genome.h) | 0.1 | 0.005 | 0.5 | g glucose | Floor cap for tiny/new nodes. Prevents zero-capacity nodes that can't receive sugar. |
| [`sugar_cap_meristem`](../src/engine/genome.h) | 2.0 | 0.1 | 10.0 | g glucose | Cap for meristem tips. Must hold enough sugar for active growth episodes. |
| [`leaf_phototropism_rate`](../src/engine/genome.h) | 0.02 | 0.001 | 0.1 | rad/tick | How fast leaves reorient toward light. Default ~1.1 deg/hr, full correction in ~3 days. |

**Key interactions:** Sugar flows by gradient diffusion (no directional bias). Conductance is throttled by the thinner of two adjacent nodes' radii. Nodes that run out of sugar for `starvation_ticks_max` (world param) die.

---

## Leaf Abscission (Carbon-Balance)

Leaves that are a net sugar drain for too long enter senescence and drop. This is the primary abscission pathway (ethylene is the secondary stress-driven pathway).

| Parameter | Default | Min | Max | Units | Description |
|-----------|---------|-----|-----|-------|-------------|
| [`leaf_abscission_ticks`](../src/engine/genome.h) | 500 | 48 | 2000 | ticks | Consecutive ticks of net sugar deficit before a leaf begins senescence. Default ~21 days. Generous — gives leaves time to recover from temporary shade. |
| [`min_leaf_age_before_abscission`](../src/engine/genome.h) | 240 | 24 | 1000 | ticks | Young leaves below this age are immune to carbon-balance abscission. Default ~10 days. Allows new leaves time to grow to productive size. |

---

## Gibberellin

Gibberellic acid (GA) promotes internode elongation. Produced locally by young leaves — acts as a "this part of the stem is actively growing" signal.

| Parameter | Default | Min | Max | Units | Description |
|-----------|---------|-----|-----|-------|-------------|
| [`ga_production_rate`](../src/engine/genome.h) | 0.5 | 0.01 | 2.0 | GA/dm leaf/tick | GA produced per dm of `leaf_size` per tick by young leaves. Applied directly to the leaf's parent and grandparent stem nodes (local effect only). |
| [`ga_leaf_age_max`](../src/engine/genome.h) | 168 | 24 | 2000 | ticks | Only leaves younger than this produce GA. Default 7 days. Mature leaves stop signaling for elongation. |
| [`ga_elongation_sensitivity`](../src/engine/genome.h) | 2.0 | 0.1 | 5.0 | multiplier | How strongly GA boosts internode elongation rate. Applied as `elongation_rate *= (1 + ga * sensitivity)`. Default doubles elongation at GA=0.5. |
| [`ga_length_sensitivity`](../src/engine/genome.h) | 1.5 | 0.1 | 5.0 | multiplier | How strongly GA boosts max internode length target. Applied as `max_length *= (1 + ga * sensitivity)`. |
| [`ga_diffusion_rate`](../src/engine/genome.h) | 0.2 | 0.01 | 1.0 | fraction/tick | GA diffusion via tree transport. GA is primarily local (applied to parent/grandparent), but does diffuse slightly. |
| [`ga_decay_rate`](../src/engine/genome.h) | 0.15 | 0.01 | 0.5 | fraction/tick | GA decay per tick. GA is **reset to zero every tick** (signal model) then re-produced, so decay applies to the freshly-produced amount during transport. |

**Key interactions:** GA is a signal, not a resource — reset each tick, recalculated. Only affects stem internodes near young leaves. Interacts with stress hormone which can inhibit elongation (see `stress_elongation_inhibition`).

---

## Ethylene

Ethylene is a stress/crowding gas signal that triggers leaf abscission. Unlike auxin/cytokinin, it diffuses spatially through 3D space (not through the tree graph).

| Parameter | Default | Min | Max | Units | Description |
|-----------|---------|-----|-----|-------|-------------|
| [`ethylene_starvation_rate`](../src/engine/genome.h) | 0.3 | 0.01 | 2.0 | signal/tick | Ethylene produced when a leaf has zero sugar. Starving leaves signal distress. |
| [`ethylene_shade_rate`](../src/engine/genome.h) | 0.2 | 0.01 | 2.0 | signal/tick | Ethylene produced when leaf `light_exposure < ethylene_shade_threshold`. Shaded leaves signal for self-thinning. |
| [`ethylene_shade_threshold`](../src/engine/genome.h) | 0.3 | 0.05 | 1.0 | fraction (0-1) | Light exposure below which a leaf produces shade-ethylene. At 0.3, leaves getting <30% of max light start signaling. |
| [`ethylene_age_rate`](../src/engine/genome.h) | 0.05 | 0.001 | 0.5 | signal/tick | Ethylene production ramp from old age. Kicks in after `ethylene_age_onset`. |
| [`ethylene_age_onset`](../src/engine/genome.h) | 720 | 100 | 5000 | ticks | Leaf age (ticks) when age-driven ethylene production begins. Default 30 days. |
| [`ethylene_crowding_rate`](../src/engine/genome.h) | 0.1 | 0.01 | 1.0 | signal/neighbor/tick | Ethylene produced per nearby node within `ethylene_crowding_radius`. Crowded regions produce more ethylene, triggering thinning. |
| [`ethylene_crowding_radius`](../src/engine/genome.h) | 0.5 | 0.1 | 3.0 | dm | Radius for counting nearby nodes in the crowding check. |
| [`ethylene_diffusion_radius`](../src/engine/genome.h) | 1.0 | 0.1 | 5.0 | dm | How far the ethylene gas cloud spreads through 3D space. Each emitting node affects all nodes within this radius, attenuated by distance. |
| [`ethylene_abscission_threshold`](../src/engine/genome.h) | 0.5 | 0.05 | 2.0 | signal | Ethylene level that triggers leaf senescence. Once exceeded, the leaf begins yellowing and is removed after `senescence_duration` ticks. |
| [`ethylene_elongation_inhibition`](../src/engine/genome.h) | 1.0 | 0.1 | 5.0 | multiplier | How strongly ethylene suppresses internode elongation in nearby stems. Applied as `elongation *= max(0, 1 - ethylene * inhibition)`. |
| [`senescence_duration`](../src/engine/genome.h) | 48 | 12 | 500 | ticks | Ticks from senescence start to leaf removal. Default 2 days. During this time the leaf yellows (visual) but is still in the graph. |

**Key interactions:** Ethylene is **reset to zero every tick** (signal model). Four production triggers: starvation, shade, age, crowding. Spatial diffusion means one stressed leaf can trigger abscission in neighbors — enabling cascade self-thinning of dense canopy interiors.

---

## Stress (Mechanical Load Response)

Stress hormone is produced when a branch is under significant mechanical load relative to its breaking capacity. Drives adaptive thickening, elongation inhibition, and negative gravitropism (growing away from the load direction).

| Parameter | Default | Min | Max | Units | Description |
|-----------|---------|-----|-----|-------|-------------|
| [`wood_density`](../src/engine/genome.h) | 50.0 | 10.0 | 200.0 | g/dm^3 | Mass per volume of wood. Serves double duty: determines stem mass (and thus gravitational load) AND structural strength (`break_stress = wood_density * break_strength_factor`). Denser wood is heavier but stronger. |
| [`wood_flexibility`](../src/engine/genome.h) | 0.5 | 0.1 | 1.0 | fraction (0-1) | Droop threshold as a fraction of break stress. At 0.5, drooping begins at 50% of breaking capacity. Lower = stiffer (droop starts closer to break point). |
| [`stress_hormone_threshold`](../src/engine/genome.h) | 0.3 | 0.0 | 0.8 | ratio (0-1) | Stress ratio (`stress / break_stress`) below which no stress hormone is produced. Default 0.3 means branches below 30% of their breaking capacity produce no hormone — they're structurally fine. Only branches approaching their limits trigger a response. |
| [`stress_hormone_production_rate`](../src/engine/genome.h) | 0.5 | 0.01 | 2.0 | signal/tick | Hormone produced per unit of excess stress ratio (above threshold). The input is normalized 0-1: `(stress_ratio - threshold) / (1 - threshold)`. |
| [`stress_hormone_diffusion_rate`](../src/engine/genome.h) | 0.15 | 0.01 | 0.5 | fraction/tick | How fast stress hormone spreads to neighboring nodes. Moderate — keeps the signal somewhat local to the stressed region. |
| [`stress_hormone_decay_rate`](../src/engine/genome.h) | 0.2 | 0.01 | 0.5 | fraction/tick | How fast stress hormone fades. At 0.2, 80% remains each tick — signal clears within ~10 ticks once load is relieved. |
| [`stress_thickening_boost`](../src/engine/genome.h) | 1.0 | 0.0 | 5.0 | multiplier | How strongly stress hormone boosts radial thickening. Applied as `thickening_rate *= (1 + stress_hormone * boost)`. Stressed stems thicken faster to increase structural capacity. |
| [`stress_elongation_inhibition`](../src/engine/genome.h) | 1.0 | 0.0 | 5.0 | multiplier | How strongly stress hormone suppresses internode elongation. Applied as `elongation *= max(0, 1 - stress_hormone * inhibition)`. Stressed regions stop elongating and focus on structural reinforcement. |
| [`stress_gravitropism_boost`](../src/engine/genome.h) | 0.5 | 0.0 | 5.0 | multiplier | How strongly stress hormone pulls shoot growth direction toward vertical. Applied at the shoot apical meristem as a blend toward the up vector, capped at 50% pull. Only branches under genuine structural threat get redirected upward. |

**How stress is computed:** For each above-ground node: `torque = child_mass * gravity * lever_arm` (lateral offset of children's center of mass), then `stress = torque / (pi * radius^2)`. The stress ratio is `stress / (wood_density * break_strength_factor)`. If this ratio exceeds `stress_hormone_threshold`, hormone is produced proportional to the excess.

**Key interactions:** Stress hormone is **reset to zero every tick** then re-produced. Three downstream effects: thickening (reinforcement), elongation inhibition (stop growing longer), gravitropism (grow upward to reduce lever arm). The `stress_hormone_threshold` is critical — it prevents structurally sound branches from being forced vertical by the gravitropism response.
