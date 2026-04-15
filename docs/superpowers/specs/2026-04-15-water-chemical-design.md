# Water Chemical — Step One Design

Add water as a capacity-based resource chemical. Roots absorb it, leaves lose it via transpiration and photosynthesis. Transport uses the existing unified transport system. No gameplay effects in step one — water deficit doesn't gate growth or cause wilting. This step establishes the chemical, tunes production/consumption/transport, and adds visualization.

## Water as a Chemical

Add `WATER` to `ChemicalID`. Water is a capacity-based resource, same model as sugar.

### Storage Capacity

Per-node capacity computed from tissue volume/area, mirroring sugar:

- **Stems and Roots**: `water_storage_density_stem * pi * r^2 * length`
- **Leaves**: `water_storage_density_leaf * leaf_area`
- **Meristems**: fixed `water_cap_meristem`
- **Minimum floor**: same as sugar (`sugar_cap_minimum` — shared floor for all resource chemicals)

### Transport

Plugs into the existing `transport_chemical()` function with its own parameter set:

- Concentration-based diffusion (level / capacity), same mechanism as sugar
- Slight positive bias to nudge water upward (models transpiration pull without explicit modeling)
- Water transports faster than sugar with less radius dependence — xylem transport is less bottlenecked than phloem
- No decay — water is a resource, not a signal

### No Step-One Effects

Water level reaching zero has no consequences. No wilting, no growth gating, no turgor mechanics. Future steps will add these.

## Production: Root Absorption

Every `RootNode` and `RootApicalNode` absorbs water in `produce()`.

```
water_absorbed = water_absorption_rate * surface_area * soil_moisture
```

- **Surface area**: `2 * pi * radius * length` for root segments; `2 * pi * radius * radius` (hemisphere approximation) for apical tips
- **soil_moisture**: new field on `WorldParams`, range 0-1, analogous to `light_level`
- Absorbed water is capped at the node's water capacity (cannot exceed cap)

## Consumption: Transpiration

Leaves lose water via transpiration in `LeafNode::produce()`.

```
water_lost = transpiration_rate * leaf_area * light_exposure
```

- `light_exposure` models stomata opening — shaded leaves transpire less
- Water level floors at zero (no negative water)
- No consequences for hitting zero in step one

## Consumption: Photosynthesis Water Cost

Small water deduction tied to sugar production in `LeafNode::produce()`.

```
water_cost = sugar_produced * photosynthesis_water_ratio
```

A fixed ratio linking the two resource economies. Deducted after sugar is produced. Water floors at zero.

## New Parameters

### Genome (evolvable)

| Parameter | Description | Units |
|-----------|-------------|-------|
| `water_absorption_rate` | Root absorption per unit surface area per unit soil moisture | ml / (dm² * hr) |
| `transpiration_rate` | Leaf water loss per unit area at full light | ml / (dm² * hr) |
| `photosynthesis_water_ratio` | Water consumed per sugar produced | ml / g |
| `water_storage_density_stem` | Water capacity per volume (StemNode + RootNode) | ml / dm³ |
| `water_storage_density_leaf` | Water capacity per leaf area | ml / dm² |
| `water_cap_meristem` | Fixed water capacity for meristem nodes | ml |
| `water_diffusion_rate` | Gradient responsiveness for transport | dimensionless |
| `water_bias` | Upward equilibrium shift (small positive value) | dimensionless |
| `water_base_transport` | Throughput floor | ml / hr |
| `water_transport_scale` | Radius scaling on throughput | ml / hr |

### WorldParams (physical constants)

| Parameter | Description | Units |
|-----------|-------------|-------|
| `soil_moisture` | Soil water availability | 0-1 range |

### Evolution Bridge

- New `water` linkage group for crossover
- All water genome params added to `build_genome_template()` with per-gene mutation config

### Renderer

- Add `water` as a color mode option in the heatmap visualization (blue gradient)

## Tick Integration

Water production and consumption happen inside existing `produce()` methods — no new global passes:

1. `RootNode::produce()` / `RootApicalNode::produce()` — absorb water
2. `LeafNode::produce()` — transpiration loss, then photosynthesis water cost (after sugar production)
3. `Node::transport_chemicals()` — water diffuses via existing unified transport (already handles any `ChemicalID` with capacity)

## What Step One Does NOT Include

- Turgor pressure / wilting
- Growth gating on water availability
- Water affecting transport capacity of other chemicals
- Transpiration pull as an explicit transport mechanism
- Humidity / spatial water vapor
- Root water-seeking tropism
