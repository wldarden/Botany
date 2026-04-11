# Gibberellin & Ethylene Hormone Systems

Two new hormones that add realistic internode variation and canopy self-pruning to the plant growth simulator.

**Codebase context:** Meristems are organized in `src/engine/meristems/` with per-type files (shoot_apical, shoot_axillary, root_apical, root_axillary) plus helpers.h and a dispatch file (meristem.cpp). Sugar has dedicated functions in sugar.h/cpp (produce, diffuse, consume, grow_leaves, prune). Hormones (auxin/cytokinin) live in hormone.h/cpp. Construction costs are in WorldParams, not Genome.

## Gibberellin (GA)

### Biology

Gibberellin promotes cell elongation in stem internodes. In real plants, young expanding leaves are the primary production site. GA acts locally — it doesn't need whole-plant transport. The effect is concentration-dependent: more GA = faster elongation and longer final internode length. Plants in shade produce more GA (shade avoidance / etiolation), causing them to stretch toward light.

### Data Model

New field on `Node`:
- `float gibberellin` — GA concentration, reset to 0 each tick (signal model, same as auxin/cytokinin)

### Genome Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ga_production_rate` | 0.5 | GA produced per dm of leaf_size per tick |
| `ga_leaf_age_max` | 168 | Only leaves younger than 7 days produce GA (ticks) |
| `ga_elongation_sensitivity` | 2.0 | How strongly GA boosts elongation rate |
| `ga_length_sensitivity` | 1.5 | How strongly GA boosts target internode length |

### Computation

Reset-and-recompute each tick (no persistence):

1. Reset `gibberellin = 0` on all nodes
2. For each LEAF node where `age < ga_leaf_age_max`:
   - `parent.gibberellin += leaf_size * ga_production_rate`
   - `grandparent.gibberellin += leaf_size * ga_production_rate * 0.3` (if grandparent exists)

### Effect on Internode Elongation

Modifies the existing intercalary growth in `src/engine/meristems/meristem.cpp` (`tick_meristems()`):

```
ga_boost = 1.0 + node.gibberellin * ga_elongation_sensitivity
effective_rate = internode_elongation_rate * ga_boost

ga_length_boost = 1.0 + node.gibberellin * ga_length_sensitivity
effective_target = max_internode_length * ga_length_boost
```

Without nearby young leaves (GA ≈ 0), elongation works exactly as before (multiplier = 1.0). Near vigorous young foliage, internodes stretch longer and faster. A weak branch with small/old leaves stays compact.

### Tick Placement

After auxin + cytokinin, before sugar. GA only needs leaf age and size (always available).

## Ethylene

### Biology

Ethylene is a gaseous hormone that diffuses through air, not vascular tissue. It's produced under stress (starvation, shade, crowding, old age) and triggers leaf abscission and growth inhibition. The spatial gas-cloud behavior means ethylene from one node affects nearby nodes regardless of tree connectivity. This creates self-thinning cascades: shaded leaves produce ethylene, nearby leaves also senesce, the branch loses sugar sources, starvation produces more ethylene, and the branch dies. This is how real canopies self-prune.

### Data Model

New fields on `Node`:
- `float ethylene` — ethylene concentration, reset to 0 each tick
- `uint32_t senescence_ticks` — 0 = healthy, >0 = senescing (persists across ticks, irreversible)

### Genome Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ethylene_starvation_rate` | 0.3 | Ethylene production when sugar = 0 |
| `ethylene_shade_rate` | 0.2 | Ethylene production from low light |
| `ethylene_shade_threshold` | 0.3 | light_exposure below which shade-ethylene kicks in |
| `ethylene_age_rate` | 0.05 | Ethylene production ramp from old age |
| `ethylene_age_onset` | 720 | Tick age (30 days) when age-ethylene starts |
| `ethylene_crowding_rate` | 0.1 | Ethylene production per nearby node |
| `ethylene_crowding_radius` | 0.5 | Radius (dm) for crowding density count |
| `ethylene_diffusion_radius` | 1.0 | Gas cloud spread distance (dm) |
| `ethylene_abscission_threshold` | 0.5 | Ethylene level triggering leaf senescence |
| `ethylene_elongation_inhibition` | 1.0 | Strength of elongation suppression |
| `senescence_duration` | 48 | Ticks (2 days) from senescence start to leaf drop |

### Computation

Reset-and-recompute each tick. Two phases:

**Phase 1 — Local production:**

For each node, accumulate ethylene from applicable triggers:
1. Sugar starvation: if `sugar <= 0` → `ethylene += ethylene_starvation_rate`
2. Low light (LEAF only): if `light_exposure < ethylene_shade_threshold` → `ethylene += ethylene_shade_rate * (1.0 - light_exposure)`
3. Old age (LEAF only): if `age > ethylene_age_onset` → `ethylene += ethylene_age_rate * (age - onset) / onset`
4. Crowding: count nodes within `ethylene_crowding_radius` in 3D space → `ethylene += ethylene_crowding_rate * nearby_count`

**Phase 2 — Spatial gas diffusion:**

Compute-then-apply pattern (order-independent, like sugar diffusion):
1. For each node with ethylene > 0, find all nodes within `ethylene_diffusion_radius`
2. Each neighbor accumulates: `source.ethylene * (1.0 - distance / radius)` (linear falloff)
3. After computing all contributions, apply them atomically

Brute-force O(n²) spatial query. Typical plant sizes (hundreds to low thousands of nodes) make this acceptable. Can add spatial hashing later if profiling shows a bottleneck.

### Effects

**Leaf abscission (gradual):**
- When a LEAF node's ethylene > `ethylene_abscission_threshold` and `senescence_ticks == 0`: begin senescence
- Each tick during senescence: increment `senescence_ticks`
- Senescing leaves stop producing sugar (checked in `produce_sugar()`)
- Visual: color interpolates green → yellow → brown based on `senescence_ticks / senescence_duration`
- When `senescence_ticks >= senescence_duration`: remove the leaf node from the tree

**Elongation inhibition (opposes GA):**
```
ethylene_inhibition = max(0.0, 1.0 - node.ethylene * ethylene_elongation_inhibition)
effective_rate = internode_elongation_rate * ga_boost * ethylene_inhibition
```

High ethylene can completely suppress elongation even if GA is high. This creates the realistic dynamic where stressed internodes stay compact despite nearby young leaves.

### Tick Placement

After sugar phase (needs `light_exposure` and sugar levels). New abscission phase runs after ethylene, before meristems.

## Engine Tick Order

```
Previous: auxin → cytokinin → sugar → meristems → positions
New:      auxin → cytokinin → GA → sugar → ethylene → abscission → meristems → positions
```

- GA before sugar: only needs leaf age/size
- Ethylene after sugar: needs light_exposure and sugar levels as trigger inputs
- Abscission after ethylene: removes completed-senescence leaves before meristem tick
- Meristems last: elongation reads both GA boost and ethylene inhibition

## New Source Files

- `src/engine/gibberellin.h` / `gibberellin.cpp` — `compute_gibberellin(Plant&)`
- `src/engine/ethylene.h` / `ethylene.cpp` — `compute_ethylene(Plant&, WorldParams&)` + `process_abscission(Plant&)`

Follows the existing pattern of one header/source pair per system (hormone.h/cpp, sugar.h/cpp).

## Modifications to Existing Files

- `src/engine/node.h` — add `gibberellin`, `ethylene`, `senescence_ticks` fields to `Node`
- `src/engine/genome.h` — add GA and ethylene parameters to `Genome` + `default_genome()`
- `src/engine/engine.h/cpp` — add GA, ethylene, and abscission calls to tick loop
- `src/engine/meristems/meristem.cpp` — modify intercalary elongation to apply `ga_boost * ethylene_inhibition`
- `src/engine/sugar.cpp` — senescing leaves (`senescence_ticks > 0`) produce zero sugar in `produce_sugar()`
- `src/renderer/renderer.cpp` — new `gibberellin` and `ethylene` color modes; senescence leaf coloring in all modes
- `CMakeLists.txt` — add `gibberellin.cpp` and `ethylene.cpp` to `botany_engine` library sources
- `tests/` — new `test_gibberellin.cpp` and `test_ethylene.cpp` files; add to CMakeLists test sources

## Renderer Changes

**New color modes:**
- `--color gibberellin` — heatmap of GA concentration
- `--color ethylene` — heatmap of ethylene concentration

**Senescence visuals (all color modes):**
- `progress = senescence_ticks / senescence_duration` (0.0 → 1.0)
- Leaf color interpolates: green → yellow → brown
- This is a material change, not a debug overlay — visible in all modes including default

## Tests

### Unit Tests

1. **GA production:** Young leaf produces GA on parent node; old leaf (age > ga_leaf_age_max) produces none
2. **GA locality:** GA appears on parent and grandparent of producing leaf, nowhere else
3. **GA elongation effect:** Node with high GA elongates faster and to a longer target than node with zero GA
4. **Ethylene starvation trigger:** Node with sugar = 0 produces ethylene; node with sugar > 0 does not (from this trigger)
5. **Ethylene shade trigger:** Leaf with low light_exposure produces ethylene proportional to shade depth
6. **Ethylene age trigger:** Old leaf (past onset) produces ethylene; young leaf does not
7. **Ethylene crowding trigger:** Node with many neighbors within radius produces more ethylene than isolated node
8. **Spatial diffusion:** Ethylene from source reaches neighbor within radius with distance falloff; does not reach node outside radius
9. **Abscission lifecycle:** Leaf above threshold enters senescence → increments senescence_ticks → stops producing sugar → removed after duration
10. **GA-ethylene interaction:** Ethylene suppresses GA-driven elongation (combined modifier < GA-only modifier)

### Integration Test

11. **Self-thinning cascade:** A cluster of shaded leaves produces ethylene → nearby leaves senesce → sugar drops → more ethylene → inner branch prunes. Verify that after sufficient ticks, interior shaded leaves are removed while outer sun-exposed leaves survive.
