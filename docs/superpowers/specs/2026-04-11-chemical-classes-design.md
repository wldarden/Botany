# Chemical Class System Design

## Goal

Replace hardcoded per-chemical fields and free-function transport with a class-based chemical system that makes adding new chemicals trivial. Chemicals are passive data definitions. Nodes own all behavior.

## Three Categories

| Category | Chemicals | Key properties |
|----------|-----------|---------------|
| **Hormone** | auxin, cytokinin, gibberellin | Reset each tick, biased local transport (child<->parent), decay |
| **Resource** | sugar | Persistent across ticks, gradient diffusion, storage caps, production/consumption |
| **Volatile** | ethylene | Reset each tick, spatial position-based diffusion cloud (not graph-based) |

## Folder Layout

```
src/engine/chemical/
  chemical.h          // ChemicalID enum, ChemicalCategory enum, ChemicalDef base struct

  hormone/
    hormone.h         // HormoneDef : ChemicalDef
    auxin.h           // AuxinDef : HormoneDef
    cytokinin.h       // CytokininDef : HormoneDef
    gibberellin.h     // GibberellinDef : HormoneDef

  resource/
    resource.h        // ResourceDef : ChemicalDef
    sugar.h           // SugarDef : ResourceDef

  volatile/
    volatile_base.h   // VolatileDef : ChemicalDef
    ethylene.h        // EthyleneDef : VolatileDef
```

## Chemical Definitions — Passive Param Structs

Chemicals are data, not actors. No `.cpp` files unless truly needed. No behavior, no utility functions.

### `chemical.h`

- `ChemicalID` enum: `Auxin, Cytokinin, Gibberellin, Sugar, Ethylene`
- `ChemicalCategory` enum: `Hormone, Resource, Volatile`
- `ChemicalDef` base struct: `id`, `name`, `category`

### Category structs (shared default params)

- `HormoneDef : ChemicalDef` — `resets_each_tick = true`, default `decay_rate`, `transport_rate`, `directional_bias`
- `ResourceDef : ChemicalDef` — `resets_each_tick = false`, default `transport_conductance`
- `VolatileDef : ChemicalDef` — `resets_each_tick = true`, default `diffusion_radius`

### Concrete chemical structs

Each fills in its specific universal values. These are constants of chemistry/physics, not evolvable:

- `AuxinDef` — basipetal bias, specific decay rate
- `CytokininDef` — acropetal bias, specific decay rate
- `GibberellinDef` — local bias, specific decay rate
- `SugarDef` — conductance scaling
- `EthyleneDef` — diffusion radius, crowding radius, production costs

## Param Ownership: Universal vs Evolvable

**Universal params** (live on ChemicalDef structs): properties that are the same for every plant. Production costs, base decay rates, transport physics. These are chemistry.

**Evolvable params** (live on Genome): properties that vary between plants. Sensitivity thresholds, production rate modifiers, transport rate overrides. These are genetics.

Genome keeps its existing flat struct. Fields like `auxin_threshold` and `ga_production_rate` stay as named members — Genome is a fixed genetic blueprint, not a dynamic container. No map-based storage for Genome.

## Node Storage — Map Replaces Hardcoded Fields

`Node` replaces:
```cpp
float auxin;
float cytokinin;
float sugar;
float gibberellin;
float ethylene;
```

With:
```cpp
std::unordered_map<ChemicalID, float> chemicals;
```

All code reading `node.auxin` migrates to `node.chemicals[ChemicalID::Auxin]` (or a convenience accessor).

## Transport — Nodes Stay In Charge

### Hormone transport (local, on Node)

`Node::transport_chemicals()` generalizes from hardcoded auxin/cytokinin blocks to a loop over all chemicals in the map, dispatching by category:

- **Hormones**: existing `transport_chemical()` biased flow helper stays in `node.cpp`. Called in a loop — each hormone's `HormoneDef` provides default rates, Genome provides per-chemical overrides.
- **Resources**: sugar's cap-aware gradient diffusion called from the same loop.
- **Volatiles**: skipped (not graph-based).

Gibberellin gets migrated from its current global reset-and-deposit pattern into the local biased transport pattern, consistent with auxin and cytokinin.

### Volatile diffusion (plant-level exception)

Ethylene's spatial position-based cloud spread is the one behavior that can't live on a node — it requires cross-node spatial queries. This stays as a plant-level pass in `Plant::tick()`.

### Abscission moves to LeafNode

`process_abscission()` moves from `ethylene.cpp` into `LeafNode::tick()`. The leaf reads its own ethylene level and manages its own senescence. Ethylene's definition has no behavior.

## Files Deleted

| Current file | Fate |
|---|---|
| `hormone.h/.cpp` | Deleted (already nearly empty) |
| `sugar.h/.cpp` | Deleted — `sugar_cap()` moves to node helpers; `grow_leaves` is dead code; diffusion/production/consumption already handled by node ticks |
| `gibberellin.h/.cpp` | Deleted — becomes `HormoneDef` + local transport |
| `ethylene.h/.cpp` | Deleted — spatial diffusion moves to Plant-level pass; `process_abscission()` moves to `LeafNode::tick()` |

## Adding a New Chemical

To add a new hormone (e.g. abscisic acid):

1. Add `AbscisicAcid` to `ChemicalID` enum in `chemical.h`
2. Create `chemical/hormone/abscisic_acid.h` with `AbscisicAcidDef : HormoneDef` filling in universal params
3. Add evolvable fields to `Genome` (sensitivity, production modifier, etc.)
4. Update the node types that produce/consume it in their `tick()` methods

No base class changes. No transport code changes. No new `.cpp` files unless the chemical has truly novel behavior.
