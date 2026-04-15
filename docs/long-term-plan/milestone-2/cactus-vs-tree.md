# Cactus vs Tree: Divergence Case Study

This is the validation exercise for milestone 2. The question it answers: given the tissue library and world physics built during this milestone, can the same engine produce a recognizable cactus-like strategy in a dry environment and a tree-like strategy in a wet, lit environment, without either being hand-designed?

If yes, milestone 2 is working. If no, something in the world physics or tissue library is still missing.

---

## What Each Strategy Actually Is

**Tree:** Compete primarily for light. Grow tall fast. Invest in secondary thickening to support height. Maintain a broad canopy that shades competitors below. Spend water freely because water is not the limiting resource.

**Cactus:** Survive primarily by retaining water. Minimize transpiring surface (no conventional leaves, or leaves vestigial). Store water in stem parenchyma as a buffer against long dry periods. Photosynthesize slowly through waxy stems. Grow extremely slowly because cell expansion is water-limited. Win by outlasting drought that kills competitors.

The strategies aren't just different morphologies — they're different solutions to which resource is scarce. In a tree-optimal environment, water is abundant and light is the bottleneck. In a cactus-optimal environment, light is abundant and water is the bottleneck. The world physics must create these two niches for the strategies to have ecological meaning.

---

## What Each Strategy Needs from the Tissue Library

### Cactus

| Need | Covered by | Status |
|---|---|---|
| No conventional leaves (or vestigial) | `max_leaf_size` → near 0 | Tunable, exists |
| Stem photosynthesis | `stem_photosynthesis_rate` parameter on STEM | New parameter, not yet added |
| Water storage in stem | Water-storage parenchyma (high water capacity, turgor model) | New tissue type, new code needed |
| Very slow growth | Low `growth_rate`, `internode_elongation_rate` | Tunable, exists |
| Waxy cuticle | `stem_cuticle_thickness` parameter | New parameter, not yet added |
| Wide shallow roots | Low `root_gravitropism_strength`, high `root_branch_angle` | Tunable, exists |
| CAM photosynthesis timing | Day/night cycle + stomatal aperture model | World physics + future layer |
| Spines | Not needed for selection pressure in animal-less world | Can skip |

**Blockers:** Water-storage parenchyma is the one non-negotiable. Without it, cactus plants are just slow-growing plants with small leaves — not plants that strategically store water to survive drought. The `stem_photosynthesis_rate` and `stem_cuticle_thickness` parameters are smaller additions but important for full cactus morphology.

### Tree

| Need | Covered by | Status |
|---|---|---|
| Tall fast growth | High `growth_rate`, `internode_elongation_rate` | Tunable, exists |
| Strong apical dominance | High `auxin_production_rate`, appropriate decay | Tunable, exists |
| Secondary thickening for structural support | `thickening_rate`, auxin-driven | Exists (approximate) |
| Broad canopy to shade competitors | `branch_angle`, phyllotaxis, `max_leaf_size` | Tunable, exists |
| Deep lateral roots | `root_gravitropism_strength`, `root_branch_angle`, `root_growth_rate` | Tunable, exists |
| Deciduous leaf cycle (optional) | Abscission exists; seasonal trigger missing | World needs seasonality |
| Bark (water retention from stem surfaces) | `stem_cuticle_thickness` at medium value | New parameter |

**Blockers:** The tree strategy is closer to what already works. The main additions needed from world physics are light occlusion (so tall trees actually shade competitors) and light competition (so broad canopy is rewarded with more sugar). Without light occlusion, there's no fitness advantage to height; the tree's fast-growing, height-seeking behavior is costly but gains nothing.

---

## What World Physics Each Strategy Needs

| Physics feature | Why cactus needs it | Why tree needs it |
|---|---|---|
| Soil water model | Water is the resource cactus hoards; drought makes cactus viable | Trees need reliable water supply; roots compete in wet zones |
| Evapotranspiration | Penalizes large leaves in dry environments — the key cactus pressure | Creates cost/benefit tradeoff for canopy size |
| Light occlusion | Minimal benefit (cactus is short) | Trees must shade competitors to win; without this, no benefit to height |
| Moisture gradient | Cactus needs a dry zone where its strategy wins | Trees need a wet zone where water isn't limiting |
| Root competition | Matters in dry zones where soil water is scarce | Matters in any zone; roots from multiple plants share soil cells |

The two strategies need the same physics but respond to different ends of the gradient.

---

## Success Criteria

Milestone 2's cactus↔tree validation passes when:

1. In the dry zone of a moisture-gradient world, a population under selection pressure converges on: small or absent leaves, stem photosynthesis, high cuticle, water-storage stem tissue, slow growth, wide shallow roots.
2. In the wet zone with light competition, a population under selection pressure converges on: large leaves, tall fast growth, strong apical dominance, deep or broad roots.
3. These are emergent outcomes from the genome evolution system, not hand-tuned starting points.
4. The two populations are visually distinct — the cactus zone looks like desert plants; the tree zone looks like a forest.

The third criterion is the hardest. It requires that the fitness function rewards ecological success (resource capture relative to neighbors) and that the genome mutation rates allow sufficient exploration. Neither of those is fully settled yet.

---

## Third Pole: Grass

If the case study were extended to a third strategy, grass would be the most structurally interesting — not because it's between cactus and tree on the water axis, but because it represents a completely different growth architecture (intercalary meristems, tiller branching from the base, no secondary thickening) that can occupy the same environmental conditions as a tree but with a different competitive strategy: low-growing, fire-resistant, grazing-adapted, fast recovery from disturbance.

In an animal-less world, the grass advantage is diminished (no grazing pressure that the basal meristem position specifically counters). But rhizome-connected clonal spread and the ability to colonize open ground are still real competitive advantages. Grass as a third validation case for milestone 2 is probably out of scope — it requires intercalary meristems (new code) and its selective advantage is weaker without herbivory. But it's the right next test after cactus↔tree divergence is working.
