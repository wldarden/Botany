# Stem Cuticle and Visual Expression

## The Idea

A stem cuticle parameter controls water loss through stem surfaces — a real physical property that varies enormously between plant types. A cactus has an extremely thick, waxy cuticle that nearly eliminates stem transpiration. A tropical rain-forest tree has a thin cuticle; water loss through bark is relatively high but rainfall keeps it irrelevant. A young herbaceous stem has almost no cuticle at all and would desiccate in an hour in a desert.

This parameter does double duty: it drives real sim behavior (water loss rate from STEM nodes) and it has a distinctive visual expression. A waxy cuticle looks waxy — specular highlights, slight blue-green sheen, smooth surface. A thin cuticle on a wet tropical stem looks wet and dark. A thick corky cuticle on an old tree looks dull and textured, like bark.

Even if the physics contribution were minimal, the visual contribution is load-bearing. The sim becomes legible when you can look at a plant and read its strategy from its appearance. A stem that looks waxy and bluish should be a drought-adapted plant. A stem that looks wet and dark green should be a water-rich environment plant. This visual feedback loop is how you tell at a glance whether evolution is doing something interesting.

---

## How the Parameter Works

**`stem_cuticle_thickness`** — a genome parameter on STEM (and ROOT?) nodes.

Mechanically: scales water loss rate from stem surfaces. In the unified water model, every STEM node loses a small amount of water each tick proportional to its surface area and the vapor pressure of the environment. Cuticle thickness is a multiplier on this loss rate:
- Very high cuticle (cactus: 0.9–1.0) → near-zero stem water loss.
- Medium cuticle (typical tree: 0.3–0.5) → moderate stem water loss, not a bottleneck when water is available.
- Very low cuticle (herbaceous/aquatic: 0.0–0.1) → high stem water loss; plant must be continuously supplied with water to stay hydrated.

The parameter is evolvable. In a dry world, high cuticle is strongly selected for because it reduces the continuous drain on water reserves. In a wet world, there's no cost to low cuticle and no benefit to high cuticle — selection is neutral and other traits dominate.

---

## Visual Expression

Three distinct visual regimes worth rendering differently:

**Waxy / succulent (high cuticle).** Stems appear lighter, slightly blue-green or gray-green, with specular highlight (shinier surface). In OpenGL terms: higher specular coefficient, slight color shift toward cool desaturated green or gray-green. This is what a cactus, aloe, or jade plant looks like.

**Barky / mature (medium-high cuticle, old wood).** Dull surface, brown-gray, low specular. Cuticle that has transitioned from waxy green to corky periderm. Can be approximated by blending toward brown and reducing specularity as stem age increases — the cuticle parameter sets the ceiling, but old stems trend toward matte bark appearance regardless.

**Wet / tropical (low cuticle).** Stems appear darker and more saturated — the wet-look of a plant that's fully hydrated and not waxy. Slight darkening and saturation boost to the base green. Low specular (wet surface is actually diffuse at scale, not specular unless it's running water).

These could also interact with the water content of the node, once water is modeled: a high-cuticle cactus that's drought-stressed (low water content) would look more matte/pale than a well-hydrated one, even with the same cuticle parameter. Visual hydration status = a direct readout of the water model.

---

## Other Genome Parameters That Deserve Visual Expression

The cuticle idea generalizes. Several genome parameters already affect plant behavior in ways that could be made visually legible:

**Wood density.** Already affects stress/droop. Dense wood could render darker and slightly cooler in color (think dense tropical hardwood vs pale fast-grown softwood). Old nodes = denser wood over time from thickening → color transition from young pale wood toward dark heartwood.

**Sugar content (current).** A node with high sugar could render with slightly more saturation or warmth — well-fed tissue is metabolically active and often greener. Starvation stress → pale, yellowish, pre-abscission appearance. The abscission yellowing for leaves already does this partially.

**Water content (future).** Turgor-high node = slightly swollen radius, darker/more saturated color. Turgor-low = slightly shrunken, pale, matte. This is the most powerful one — you could see a plant wilt in real time as its water level drops.

**Stress hormone level.** High stress → reddish tint (anthocyanins are produced under mechanical stress and drought in many real plants, giving a reddish or purplish cast). Not just autumn coloring — stress response mid-growth season.

**Auxin concentration.** The current color mode already shows auxin as a heatmap. In a "default" color mode, high-auxin zones (near shoot tips) might have subtly more yellow-green coloring (young, hormonally active tissue vs mature green).

**Leaf age.** Already does yellowing via the senescence system. Could extend: very young leaves could be lighter green or slightly reddish (many real plants have red/bronze young leaves due to anthocyanin before chlorophyll fully develops).

---

## Philosophy

Visual expressiveness of genome parameters is load-bearing for the sim, not decorative. When a parameter has both mechanical and visual expression, the sim becomes readable without the debug heatmaps. You should be able to watch evolution run and see the population shifting — stems getting waxier as drought pressure increases, leaf size shrinking, root mass deepening — just by watching the visual output. That's the version of this sim that's compelling to look at.

The rule of thumb: any genome parameter that describes a *material property* of plant tissue should have a visual expression. Rates and thresholds (auxin_production_rate, cytokinin_threshold) probably don't need visual encoding. But tissue properties — density, cuticle, water content, sugar content, age — all have real optical correlates in actual plant tissue.
