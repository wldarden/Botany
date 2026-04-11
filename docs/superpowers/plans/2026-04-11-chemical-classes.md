# Chemical Class System Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace hardcoded per-chemical fields on Node with a map-based chemical system organized by category (Hormone/Resource/Volatile), making it trivial to add new chemicals.

**Architecture:** Chemical definitions are passive param structs in `src/engine/chemical/`. Node storage becomes `unordered_map<ChemicalID, float>`. Node::transport_chemicals() loops over a hormone registry. Gibberellin migrates from a global pass to local biased transport. Abscission moves from ethylene.cpp into LeafNode::tick().

**Tech Stack:** C++17, existing engine patterns

**Spec:** `docs/superpowers/specs/2026-04-11-chemical-classes-design.md`

---

### Task 1: Create chemical definition framework

**Files:**
- Create: `src/engine/chemical/chemical.h`
- Create: `src/engine/chemical/hormone/hormone.h`
- Create: `src/engine/chemical/hormone/auxin.h`
- Create: `src/engine/chemical/hormone/cytokinin.h`
- Create: `src/engine/chemical/hormone/gibberellin.h`
- Create: `src/engine/chemical/resource/resource.h`
- Create: `src/engine/chemical/resource/sugar.h`
- Create: `src/engine/chemical/volatile/volatile_base.h`
- Create: `src/engine/chemical/volatile/ethylene.h`
- Create: `src/engine/chemical/chemical_registry.h`

- [ ] **Step 1: Create `src/engine/chemical/chemical.h`**

```cpp
// src/engine/chemical/chemical.h
#pragma once

#include <cstdint>

namespace botany {

enum class ChemicalID : uint8_t {
    Auxin,
    Cytokinin,
    Gibberellin,
    Sugar,
    Ethylene,
};

enum class ChemicalCategory : uint8_t {
    Hormone,
    Resource,
    Volatile,
};

struct ChemicalDef {
    ChemicalID id;
    const char* name;
    ChemicalCategory category;
};

} // namespace botany
```

- [ ] **Step 2: Create `src/engine/chemical/hormone/hormone.h`**

```cpp
// src/engine/chemical/hormone/hormone.h
#pragma once

#include "engine/chemical/chemical.h"
#include "engine/genome.h"

namespace botany {

// Params extracted from a plant's genome for one hormone.
// Used by Node::transport_chemicals() to loop over all hormones generically.
struct HormoneTransportParams {
    ChemicalID id;
    float transport_rate;
    float directional_bias;
    float decay_rate;
};

} // namespace botany
```

- [ ] **Step 3: Create concrete hormone defs**

`src/engine/chemical/hormone/auxin.h`:
```cpp
#pragma once

#include "engine/chemical/chemical.h"

namespace botany {

inline constexpr ChemicalDef auxin_def = {
    ChemicalID::Auxin, "auxin", ChemicalCategory::Hormone
};

} // namespace botany
```

`src/engine/chemical/hormone/cytokinin.h`:
```cpp
#pragma once

#include "engine/chemical/chemical.h"

namespace botany {

inline constexpr ChemicalDef cytokinin_def = {
    ChemicalID::Cytokinin, "cytokinin", ChemicalCategory::Hormone
};

} // namespace botany
```

`src/engine/chemical/hormone/gibberellin.h`:
```cpp
#pragma once

#include "engine/chemical/chemical.h"

namespace botany {

inline constexpr ChemicalDef gibberellin_def = {
    ChemicalID::Gibberellin, "gibberellin", ChemicalCategory::Hormone
};

} // namespace botany
```

- [ ] **Step 4: Create resource and volatile defs**

`src/engine/chemical/resource/resource.h`:
```cpp
#pragma once

#include "engine/chemical/chemical.h"

namespace botany {

// Resources persist across ticks, use gradient diffusion.
// Transport logic is unique per resource (sugar has cap-aware diffusion).

} // namespace botany
```

`src/engine/chemical/resource/sugar.h`:
```cpp
#pragma once

#include "engine/chemical/chemical.h"

namespace botany {

inline constexpr ChemicalDef sugar_def = {
    ChemicalID::Sugar, "sugar", ChemicalCategory::Resource
};

} // namespace botany
```

`src/engine/chemical/volatile/volatile_base.h`:
```cpp
#pragma once

#include "engine/chemical/chemical.h"

namespace botany {

// Volatiles use spatial (position-based) diffusion, not graph transport.
// Diffused at the Plant level, not in Node::transport_chemicals().

} // namespace botany
```

`src/engine/chemical/volatile/ethylene.h`:
```cpp
#pragma once

#include "engine/chemical/chemical.h"

namespace botany {

inline constexpr ChemicalDef ethylene_def = {
    ChemicalID::Ethylene, "ethylene", ChemicalCategory::Volatile
};

} // namespace botany
```

- [ ] **Step 5: Create `src/engine/chemical/chemical_registry.h`**

```cpp
// src/engine/chemical/chemical_registry.h
#pragma once

#include "engine/chemical/chemical.h"
#include "engine/chemical/hormone/hormone.h"
#include "engine/genome.h"
#include <array>

namespace botany {

// All chemical IDs in the system. Add new chemicals here.
inline constexpr std::array<ChemicalID, 5> all_chemical_ids = {
    ChemicalID::Auxin,
    ChemicalID::Cytokinin,
    ChemicalID::Gibberellin,
    ChemicalID::Sugar,
    ChemicalID::Ethylene,
};

// Extract hormone transport params from a genome.
// Add new hormones here — one line per hormone.
inline std::array<HormoneTransportParams, 3> hormone_params(const Genome& g) {
    return {{
        {ChemicalID::Auxin, g.auxin_transport_rate, g.auxin_directional_bias, g.auxin_decay_rate},
        {ChemicalID::Cytokinin, g.cytokinin_transport_rate, g.cytokinin_directional_bias, g.cytokinin_decay_rate},
        {ChemicalID::Gibberellin, g.ga_transport_rate, g.ga_directional_bias, g.ga_decay_rate},
    }};
}

} // namespace botany
```

**Note:** The `ChemicalDef` instances (`auxin_def`, etc.) are metadata — nothing in the code references them yet. They serve as documentation and a canonical place to add universal (non-genome) params later (e.g., production costs, molar weights). The transport loop uses `hormone_params()` which extracts from Genome directly.

**Note:** The spec says sugar.h/.cpp and ethylene.h/.cpp are "Deleted." The plan deviates: we keep `sugar.cpp` (for `sugar_cap` + `compute_light_exposure`) and `ethylene.cpp` (for `compute_ethylene` spatial diffusion). These are plant-level utility functions, not chemical definitions. Moving them to plant.cpp would work too but adds no value.

- [ ] **Step 6: Build to verify new files compile**

Run: `/usr/local/bin/cmake --build build`
Expected: Success — new files are header-only, not yet included anywhere.

- [ ] **Step 7: Commit**

```bash
git add src/engine/chemical/
git commit -m "feat: chemical definition framework — enums, defs, and registry"
```

---

### Task 2: Add GA transport genome params

**Files:**
- Modify: `src/engine/genome.h`

Gibberellin needs transport params to work as a hormone. Currently it has `ga_production_rate`, `ga_leaf_age_max`, `ga_elongation_sensitivity`, `ga_length_sensitivity`. It needs `ga_transport_rate`, `ga_directional_bias`, `ga_decay_rate`.

- [ ] **Step 1: Add GA transport fields to Genome**

In `src/engine/genome.h`, after the existing GA fields:

```cpp
    // Gibberellin — promotes internode elongation, produced by young leaves
    float ga_production_rate;         // GA produced per dm of leaf_size per tick
    uint32_t ga_leaf_age_max;         // ticks — only leaves younger than this produce GA
    float ga_elongation_sensitivity;  // how strongly GA boosts elongation rate
    float ga_length_sensitivity;      // how strongly GA boosts target internode length
    float ga_transport_rate;          // fraction transported per tick (biased transport)
    float ga_directional_bias;        // -1=basipetal, 0=gradient, +1=acropetal
    float ga_decay_rate;              // exponential decay per tick
```

- [ ] **Step 2: Add defaults in `default_genome()`**

After the existing GA defaults:

```cpp
        // Gibberellin
        .ga_production_rate = 0.5f,
        .ga_leaf_age_max = 168,               // 7 days
        .ga_elongation_sensitivity = 2.0f,
        .ga_length_sensitivity = 1.5f,
        .ga_transport_rate = 0.2f,            // moderate transport
        .ga_directional_bias = -0.7f,         // mostly basipetal (leaf -> parent -> trunk)
        .ga_decay_rate = 0.15f,               // decays faster than auxin — short-range signal
```

Note: GA is basipetal (negative bias) because it's produced in leaves and acts on the internode below. Higher decay than auxin ensures it stays local to the producing leaf's neighborhood.

- [ ] **Step 3: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`
Expected: All tests pass — new fields are just added, not used yet.

- [ ] **Step 4: Commit**

```bash
git add src/engine/genome.h
git commit -m "feat: add gibberellin transport params to Genome"
```

---

### Task 3: Add chemicals map to Node + accessor

**Files:**
- Modify: `src/engine/node/node.h`
- Modify: `src/engine/node/node.cpp`

Add the `unordered_map<ChemicalID, float>` alongside the existing fields. Add `chemical()` accessors. Both storage systems coexist during migration.

- [ ] **Step 1: Add map and accessor to node.h**

Add includes at top of `src/engine/node/node.h`:
```cpp
#include <unordered_map>
#include "engine/chemical/chemical.h"
```

Add to the Node class, after the existing chemical fields:

```cpp
    // Chemical storage — map-based, replaces individual fields.
    // During migration, both map and old fields exist.
    std::unordered_map<ChemicalID, float> chemicals;

    float& chemical(ChemicalID id) { return chemicals[id]; }
    float chemical(ChemicalID id) const {
        auto it = chemicals.find(id);
        return it != chemicals.end() ? it->second : 0.0f;
    }
```

We need a hash for ChemicalID. Add before the Node class:
```cpp
} // namespace botany

// Hash for ChemicalID so it works as unordered_map key
template<>
struct std::hash<botany::ChemicalID> {
    std::size_t operator()(botany::ChemicalID id) const noexcept {
        return static_cast<std::size_t>(id);
    }
};

namespace botany {
```

- [ ] **Step 2: Initialize map in Node constructor**

In `src/engine/node/node.cpp`, in the Node constructor body, add after the initializer list:

```cpp
Node::Node(uint32_t id, NodeType type, glm::vec3 position, float radius)
    : id(id)
    , parent(nullptr)
    , offset(position)
    , position(position)
    , radius(radius)
    , type(type)
    , age(0)
    , auxin(0.0f)
    , cytokinin(0.0f)
{
    // Initialize all chemical map entries to zero
    chemicals[ChemicalID::Auxin] = 0.0f;
    chemicals[ChemicalID::Cytokinin] = 0.0f;
    chemicals[ChemicalID::Gibberellin] = 0.0f;
    chemicals[ChemicalID::Sugar] = 0.0f;
    chemicals[ChemicalID::Ethylene] = 0.0f;
}
```

Add include at top: `#include "engine/chemical/chemical.h"`

- [ ] **Step 3: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`
Expected: All tests pass — map exists alongside old fields, nothing reads it yet.

- [ ] **Step 4: Commit**

```bash
git add src/engine/node/node.h src/engine/node/node.cpp
git commit -m "feat: add chemicals map + accessor to Node (dual storage)"
```

---

### Task 4: Migrate all code from old fields to map

**Files:**
- Modify: `src/engine/node/node.cpp`
- Modify: `src/engine/node/stem_node.cpp`
- Modify: `src/engine/node/root_node.cpp`
- Modify: `src/engine/node/leaf_node.cpp`
- Modify: `src/engine/node/meristems/shoot_apical.cpp`
- Modify: `src/engine/node/meristems/shoot_axillary.cpp`
- Modify: `src/engine/node/meristems/root_apical.cpp`
- Modify: `src/engine/node/meristems/root_axillary.cpp`
- Modify: `src/engine/sugar.cpp`
- Modify: `src/engine/gibberellin.cpp`
- Modify: `src/engine/ethylene.cpp`
- Modify: `src/engine/plant.cpp`
- Modify: `src/app_realtime.cpp`
- Modify: `src/serialization/serializer.h`
- Modify: `src/serialization/serializer.cpp`
- Modify: All test files that reference chemical fields

This is a mechanical find-and-replace. Every occurrence of a direct chemical field access becomes a map accessor call. Both map and old fields are kept in sync by writing to both during this step.

- [ ] **Step 1: Add sync helpers to Node**

In `src/engine/node/node.cpp`, add a private helper that syncs map → old fields after every write. Actually, the simpler approach: for this task, every place that writes to a chemical writes to BOTH the map and the old field. Every place that reads, reads from the map.

The key replacements across all files:

**Reads** (switch to map):
- `node.auxin` → `node.chemical(ChemicalID::Auxin)`
- `node.cytokinin` → `node.chemical(ChemicalID::Cytokinin)`
- `node.sugar` → `node.chemical(ChemicalID::Sugar)`
- `node.gibberellin` → `node.chemical(ChemicalID::Gibberellin)`
- `node.ethylene` → `node.chemical(ChemicalID::Ethylene)`

Same patterns for `node->auxin`, `n.auxin`, `sel.auxin`, etc.

**Writes** (write to BOTH map and old field, until old fields removed):
- `node.auxin = X` → `node.chemical(ChemicalID::Auxin) = X; node.auxin = X;`
- `node.auxin += X` → `node.chemical(ChemicalID::Auxin) += X; node.auxin += X;`

This dual-write ensures nothing breaks during migration.

- [ ] **Step 2: Migrate `src/engine/node/node.cpp` transport_chemicals()**

This is the most critical migration — transport reads/writes chemicals on both `this` and `parent`.

Replace the hormone transport section:
```cpp
void Node::transport_chemicals(const Genome& g) {
    if (parent) {
        // Auxin: basipetal
        transport_chemical(chemical(ChemicalID::Auxin), parent->chemical(ChemicalID::Auxin),
            g.auxin_transport_rate, g.auxin_directional_bias, g.auxin_decay_rate);
        // Sync old fields
        auxin = chemical(ChemicalID::Auxin);
        parent->auxin = parent->chemical(ChemicalID::Auxin);

        // Cytokinin: acropetal
        transport_chemical(chemical(ChemicalID::Cytokinin), parent->chemical(ChemicalID::Cytokinin),
            g.cytokinin_transport_rate, g.cytokinin_directional_bias, g.cytokinin_decay_rate);
        cytokinin = chemical(ChemicalID::Cytokinin);
        parent->cytokinin = parent->chemical(ChemicalID::Cytokinin);

        // Sugar: gradient-based, cap-aware
        float my_cap = sugar_cap(*this, g);
        float parent_cap = sugar_cap(*parent, g);
        float diff = chemical(ChemicalID::Sugar) - parent->chemical(ChemicalID::Sugar);
        float min_r = std::min(radius, parent->radius);
        if (type == NodeType::LEAF || is_meristem()
            || parent->type == NodeType::LEAF || parent->is_meristem())
            min_r = std::max(min_r, 0.01f);
        float conductance = std::min(min_r * min_r * 3.14159f * g.sugar_transport_conductance, 0.25f);
        float flow = diff * conductance;
        if (flow > 0.0f) {
            float headroom = std::max(0.0f, parent_cap - parent->chemical(ChemicalID::Sugar));
            flow = std::min({flow, chemical(ChemicalID::Sugar), headroom});
        } else {
            float headroom = std::max(0.0f, my_cap - chemical(ChemicalID::Sugar));
            flow = std::max({flow, -parent->chemical(ChemicalID::Sugar), -headroom});
        }
        chemical(ChemicalID::Sugar) -= flow;
        parent->chemical(ChemicalID::Sugar) += flow;
        sugar = chemical(ChemicalID::Sugar);
        parent->sugar = parent->chemical(ChemicalID::Sugar);
    } else {
        // Seed: no parent, just decay hormones
        chemical(ChemicalID::Auxin) *= (1.0f - g.auxin_decay_rate);
        chemical(ChemicalID::Cytokinin) *= (1.0f - g.cytokinin_decay_rate);
        auxin = chemical(ChemicalID::Auxin);
        cytokinin = chemical(ChemicalID::Cytokinin);
    }
}
```

- [ ] **Step 3: Migrate Node::tick() in node.cpp**

Replace `sugar` references with map:
```cpp
void Node::tick(Plant& plant, const WorldParams& world) {
    age++;
    const Genome& g = plant.genome();

    float cost = maintenance_cost(g);
    chemical(ChemicalID::Sugar) = std::max(0.0f, chemical(ChemicalID::Sugar) - cost);
    sugar = chemical(ChemicalID::Sugar);

    float cap = sugar_cap(*this, g);
    chemical(ChemicalID::Sugar) = std::min(chemical(ChemicalID::Sugar), cap);
    sugar = chemical(ChemicalID::Sugar);

    if (chemical(ChemicalID::Sugar) <= 0.0f) starvation_ticks++;
    else starvation_ticks = 0;

    if (starvation_ticks >= world.starvation_ticks_max && parent) {
        die(plant);
        return;
    }

    transport_chemicals(g);
}
```

- [ ] **Step 4: Migrate meristem files**

`shoot_apical.cpp` — change `auxin +=` and `sugar` refs:
```cpp
void ShootApicalNode::tick(Plant& plant, const WorldParams& world) {
    chemical(ChemicalID::Auxin) += plant.genome().auxin_production_rate;
    auxin = chemical(ChemicalID::Auxin);
    // ... rest uses sugar via chemical(ChemicalID::Sugar) ...
```

`shoot_axillary.cpp`:
```cpp
bool ShootAxillaryNode::can_activate(const Genome& g, const WorldParams& world) const {
    float stem_auxin = parent ? parent->chemical(ChemicalID::Auxin) : chemical(ChemicalID::Auxin);
    if (stem_auxin >= g.auxin_threshold) return false;

    float parent_sugar = parent ? parent->chemical(ChemicalID::Sugar) : chemical(ChemicalID::Sugar);
    if (parent_sugar < g.sugar_activation_shoot) return false;
    if (chemical(ChemicalID::Sugar) < world.sugar_cost_activation) return false;
    return true;
}
```

`root_apical.cpp`:
```cpp
void RootApicalNode::tick(Plant& plant, const WorldParams& world) {
    chemical(ChemicalID::Cytokinin) += plant.genome().cytokinin_production_rate;
    cytokinin = chemical(ChemicalID::Cytokinin);
    // ...
```

`root_axillary.cpp`:
```cpp
bool RootAxillaryNode::can_activate(const Genome& g, const WorldParams& world) const {
    float stem_cytokinin = parent ? parent->chemical(ChemicalID::Cytokinin) : chemical(ChemicalID::Cytokinin);
    if (stem_cytokinin >= g.cytokinin_threshold) return false;

    float parent_sugar = parent ? parent->chemical(ChemicalID::Sugar) : chemical(ChemicalID::Sugar);
    if (parent_sugar < g.sugar_activation_root) return false;
    if (chemical(ChemicalID::Sugar) < world.sugar_cost_activation) return false;
    return true;
}
```

- [ ] **Step 5: Migrate stem_node.cpp and root_node.cpp**

In `StemNode::elongate()`, replace `gibberellin` and `ethylene`:
```cpp
    float ga_boost = 1.0f + chemical(ChemicalID::Gibberellin) * g.ga_elongation_sensitivity;
    float eth_inhibit = std::max(0.0f, 1.0f - chemical(ChemicalID::Ethylene) * g.ethylene_elongation_inhibition);
    // ...
    float max_len = g.max_internode_length * (1.0f + chemical(ChemicalID::Gibberellin) * g.ga_length_sensitivity);
```

Apply the same pattern to `RootNode::elongate()` — same fields, same replacements.

For `sugar` usage in `thicken()` and `elongate()`, replace `sugar` with `chemical(ChemicalID::Sugar)` and sync.

- [ ] **Step 6: Migrate leaf_node.cpp**

Replace all `sugar` references with `chemical(ChemicalID::Sugar)` and sync. This covers `photosynthesize()`, `phototropism()`, and `grow()`.

- [ ] **Step 7: Migrate sugar.cpp**

Replace `node.sugar` and `leaf.sugar` with map accessors in `sugar_cap()`, `compute_light_exposure()`, and the dead `grow_leaves()` code. For `sugar_cap()`, the function takes `const Node&` so use `node.chemical(ChemicalID::Sugar)` for any sugar reads (there are none — `sugar_cap` only reads geometry).

- [ ] **Step 8: Migrate gibberellin.cpp**

```cpp
void compute_gibberellin(Plant& plant) {
    const Genome& g = plant.genome();

    plant.for_each_node_mut([](Node& node) {
        node.chemical(ChemicalID::Gibberellin) = 0.0f;
        node.gibberellin = 0.0f;
    });

    plant.for_each_node_mut([&](Node& node) {
        auto* leaf = node.as_leaf();
        if (!leaf) return;
        if (node.age >= g.ga_leaf_age_max) return;
        if (leaf->leaf_size < 1e-6f) return;

        float production = leaf->leaf_size * g.ga_production_rate;
        if (node.parent) {
            node.parent->chemical(ChemicalID::Gibberellin) += production;
            node.parent->gibberellin = node.parent->chemical(ChemicalID::Gibberellin);
            if (node.parent->parent) {
                node.parent->parent->chemical(ChemicalID::Gibberellin) += production * 0.3f;
                node.parent->parent->gibberellin = node.parent->parent->chemical(ChemicalID::Gibberellin);
            }
        }
    });
}
```

- [ ] **Step 9: Migrate ethylene.cpp**

Replace all `node.ethylene` with `node.chemical(ChemicalID::Ethylene)` + sync, and `node.sugar` with `node.chemical(ChemicalID::Sugar)`. Same for pointer access patterns.

- [ ] **Step 10: Migrate app_realtime.cpp**

Replace chemical accessor lambdas:
```cpp
if (color_chemical == "auxin") {
    accessor = [](const Node& n) { return n.chemical(ChemicalID::Auxin); };
} else if (color_chemical == "cytokinin") {
    accessor = [](const Node& n) { return n.chemical(ChemicalID::Cytokinin); };
} else if (color_chemical == "sugar") {
    accessor = [](const Node& n) { return n.chemical(ChemicalID::Sugar); };
} else if (color_chemical == "gibberellin") {
    accessor = [](const Node& n) { return n.chemical(ChemicalID::Gibberellin); };
} else if (color_chemical == "ethylene") {
    accessor = [](const Node& n) { return n.chemical(ChemicalID::Ethylene); };
}
```

And all ImGui display code reading `sel.auxin`, `sel.sugar`, etc. Add include for `engine/chemical/chemical.h`.

- [ ] **Step 11: Migrate serializer**

`serializer.cpp` save_tick:
```cpp
write_val(out, node.chemical(ChemicalID::Auxin));
write_val(out, node.chemical(ChemicalID::Cytokinin));
write_val(out, node.chemical(ChemicalID::Sugar));
```

`serializer.h` NodeSnapshot and `serializer.cpp` load_tick can stay with named fields — NodeSnapshot is a serialization struct, not a Node. The field names `auxin`, `cytokinin`, `sugar` are fine there.

- [ ] **Step 12: Migrate all test files**

Every test that reads/writes `node.auxin`, `node->auxin`, `shoot->auxin`, `seed->auxin`, `plant.seed()->auxin`, etc. Replace with `chemical(ChemicalID::...)` or `->chemical(ChemicalID::...)`.

Key files:
- `tests/test_hormone.cpp` — `shoot->auxin`, `seed->auxin`, `root->cytokinin`, etc.
- `tests/test_sugar.cpp` — `node.sugar`, `seed->sugar`, `leaf->sugar`, etc.
- `tests/test_gibberellin.cpp` — `stem->gibberellin`, etc.
- `tests/test_ethylene.cpp` — `node.ethylene`, `stem->ethylene`, `leaf->ethylene`, etc.
- `tests/test_node.cpp`, `tests/test_plant.cpp`, `tests/test_engine.cpp` — any chemical field refs

Add `#include "engine/chemical/chemical.h"` to each test file that needs `ChemicalID`.

- [ ] **Step 13: Build and run ALL tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`
Expected: All tests pass. Both map and old fields are in sync.

- [ ] **Step 14: Commit**

```bash
git add -A
git commit -m "refactor: migrate all chemical access from fields to map"
```

---

### Task 5: Remove old chemical fields from Node

**Files:**
- Modify: `src/engine/node/node.h`
- Modify: `src/engine/node/node.cpp`
- Modify: All files that still reference old fields via sync code

- [ ] **Step 1: Remove old fields from node.h**

Remove these lines from the Node class:
```cpp
    float auxin;
    float cytokinin;
    float sugar = 0.0f;
    float gibberellin = 0.0f;
    float ethylene = 0.0f;
```

- [ ] **Step 2: Remove sync lines from node.cpp**

In the Node constructor, remove the `auxin(0.0f)` and `cytokinin(0.0f)` initializer list entries.

Across ALL files modified in Task 4, remove every "sync old field" line — i.e., every line that writes to `auxin =`, `cytokinin =`, `sugar =`, `gibberellin =`, `ethylene =` (the old field). There should be no remaining references to the old fields.

- [ ] **Step 3: Grep for any remaining old field references**

Run: `grep -rn '\.auxin\b\|->auxin\b\|\.cytokinin\b\|->cytokinin\b\|\.gibberellin\b\|->gibberellin\b\|\.ethylene\b\|->ethylene\b' src/ tests/`

Exclude: `serializer.h` (NodeSnapshot keeps named fields — that's fine).
Exclude: `chemical/` directory (def files mention the names as strings).
Exclude: comments.

Fix any remaining references.

For `sugar`: more careful grep needed since "sugar" appears in genome fields (`sugar_production_rate`, etc.) and variable names. Grep for `\.sugar\b` and `->sugar\b` on Node/pointer-to-Node contexts specifically.

- [ ] **Step 4: Build and run ALL tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`
Expected: All tests pass. Old fields are gone.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "refactor: remove old chemical fields from Node — map is sole storage"
```

---

### Task 6: Generalize hormone transport with registry

**Files:**
- Modify: `src/engine/node/node.cpp`

Replace the hardcoded per-hormone transport calls with a loop over the hormone registry. Initially includes only auxin and cytokinin (gibberellin still uses the global pass until Task 7).

- [ ] **Step 1: Add include**

Add to top of `src/engine/node/node.cpp`:
```cpp
#include "engine/chemical/chemical_registry.h"
```

- [ ] **Step 2: Replace hormone transport with registry loop**

In `Node::transport_chemicals()`, replace the explicit auxin and cytokinin blocks with:

```cpp
void Node::transport_chemicals(const Genome& g) {
    if (parent) {
        // Hormones: biased local transport via registry
        for (const auto& hp : hormone_params(g)) {
            // Skip gibberellin until it's migrated from the global pass
            if (hp.id == ChemicalID::Gibberellin) continue;
            transport_chemical(chemical(hp.id), parent->chemical(hp.id),
                hp.transport_rate, hp.directional_bias, hp.decay_rate);
        }

        // Sugar: cap-aware gradient diffusion (unique to this resource)
        float my_cap = sugar_cap(*this, g);
        float parent_cap = sugar_cap(*parent, g);
        float diff = chemical(ChemicalID::Sugar) - parent->chemical(ChemicalID::Sugar);
        float min_r = std::min(radius, parent->radius);
        if (type == NodeType::LEAF || is_meristem()
            || parent->type == NodeType::LEAF || parent->is_meristem())
            min_r = std::max(min_r, 0.01f);
        float conductance = std::min(min_r * min_r * 3.14159f * g.sugar_transport_conductance, 0.25f);
        float flow = diff * conductance;
        if (flow > 0.0f) {
            float headroom = std::max(0.0f, parent_cap - parent->chemical(ChemicalID::Sugar));
            flow = std::min({flow, chemical(ChemicalID::Sugar), headroom});
        } else {
            float headroom = std::max(0.0f, my_cap - chemical(ChemicalID::Sugar));
            flow = std::max({flow, -parent->chemical(ChemicalID::Sugar), -headroom});
        }
        chemical(ChemicalID::Sugar) -= flow;
        parent->chemical(ChemicalID::Sugar) += flow;
    } else {
        // Seed: just decay hormones
        for (const auto& hp : hormone_params(g)) {
            if (hp.id == ChemicalID::Gibberellin) continue;
            chemical(hp.id) *= (1.0f - hp.decay_rate);
        }
    }
}
```

- [ ] **Step 3: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`
Expected: All tests pass — behavior unchanged, just using registry loop.

- [ ] **Step 4: Commit**

```bash
git add src/engine/node/node.cpp
git commit -m "refactor: generalize hormone transport to use registry loop"
```

---

### Task 7: Migrate gibberellin to local hormone transport

**Files:**
- Modify: `src/engine/node/node.cpp` (remove gibberellin skip)
- Modify: `src/engine/node/leaf_node.h`
- Modify: `src/engine/node/leaf_node.cpp` (add GA production)
- Modify: `src/engine/plant.cpp` (remove compute_gibberellin call)
- Modify: `tests/test_gibberellin.cpp` (update tests for new behavior)

This is a **behavioral change**: GA goes from "deposit directly on parent + 0.3x grandparent each tick" to "produce on leaf, biased local transport spreads it." The transport params in Genome control how far it reaches. Tests need updating.

- [ ] **Step 1: Remove gibberellin skip in transport_chemicals**

In `src/engine/node/node.cpp`, in `transport_chemicals()`, remove the two lines:
```cpp
            if (hp.id == ChemicalID::Gibberellin) continue;
```
(One in the transport loop, one in the seed decay block.)

- [ ] **Step 2: Add GA production to LeafNode::tick()**

In `src/engine/node/leaf_node.cpp`, add GA production in `LeafNode::tick()` after the `Node::tick()` call:

```cpp
void LeafNode::tick(Plant& plant, const WorldParams& world) {
    Node::tick(plant, world);
    const Genome& g = plant.genome();

    // Gibberellin production: young leaves produce GA
    if (age < g.ga_leaf_age_max && leaf_size > 1e-6f && senescence_ticks == 0) {
        chemical(ChemicalID::Gibberellin) += leaf_size * g.ga_production_rate;
    }

    photosynthesize(g, world);
    phototropism(g, world);
    grow(g, world);
}
```

Add include: `#include "engine/chemical/chemical.h"`

- [ ] **Step 3: Remove compute_gibberellin call from Plant::tick()**

In `src/engine/plant.cpp`, in `Plant::tick()`:

Remove: `compute_gibberellin(*this);`
Remove: `#include "engine/gibberellin.h"`

- [ ] **Step 4: Rewrite gibberellin tests**

The old tests check `compute_gibberellin()` depositing GA on parent/grandparent. The new behavior is: leaf produces GA on itself, biased transport moves it toward parent. Tests should verify:

```cpp
TEST_CASE("Young leaf produces GA on itself", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->age = 10;
    stem->add_child(leaf);

    WorldParams wp = default_world_params();
    leaf->tick(plant, wp);

    // Leaf should have GA after its tick (produced on itself)
    REQUIRE(leaf->chemical(ChemicalID::Gibberellin) > 0.0f);
}

TEST_CASE("GA spreads to parent via transport", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->age = 10;
    stem->add_child(leaf);

    WorldParams wp = default_world_params();
    // Tick multiple times to let transport propagate
    for (int i = 0; i < 10; i++) {
        plant.tick(wp);
    }

    // Parent stem should have received GA via biased transport
    REQUIRE(stem->chemical(ChemicalID::Gibberellin) > 0.0f);
}

TEST_CASE("Old leaf produces no GA", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->age = g.ga_leaf_age_max + 1;
    stem->add_child(leaf);

    WorldParams wp = default_world_params();
    leaf->tick(plant, wp);

    REQUIRE(leaf->chemical(ChemicalID::Gibberellin) < 1e-6f);
}

TEST_CASE("GA decays without production", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);
    stem->chemical(ChemicalID::Gibberellin) = 1.0f;

    WorldParams wp = default_world_params();
    for (int i = 0; i < 20; i++) {
        plant.tick(wp);
    }

    // Should have decayed significantly
    REQUIRE(stem->chemical(ChemicalID::Gibberellin) < 0.1f);
}
```

- [ ] **Step 5: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`
Expected: All tests pass. GA now uses local transport.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/node.cpp src/engine/node/leaf_node.h src/engine/node/leaf_node.cpp src/engine/plant.cpp tests/test_gibberellin.cpp
git commit -m "feat: migrate gibberellin from global deposit to local biased transport"
```

---

### Task 8: Move abscission into LeafNode

**Files:**
- Modify: `src/engine/node/leaf_node.h`
- Modify: `src/engine/node/leaf_node.cpp`
- Modify: `src/engine/plant.cpp` (remove process_abscission call)
- Modify: `tests/test_ethylene.cpp` (update abscission tests)

- [ ] **Step 1: Add abscission logic to LeafNode::tick()**

In `src/engine/node/leaf_node.cpp`, add abscission handling in `LeafNode::tick()`, after growth:

```cpp
void LeafNode::tick(Plant& plant, const WorldParams& world) {
    Node::tick(plant, world);
    const Genome& g = plant.genome();

    // GA production
    if (age < g.ga_leaf_age_max && leaf_size > 1e-6f && senescence_ticks == 0) {
        chemical(ChemicalID::Gibberellin) += leaf_size * g.ga_production_rate;
    }

    photosynthesize(g, world);
    phototropism(g, world);
    grow(g, world);

    // Abscission: ethylene triggers senescence, then leaf drop
    if (senescence_ticks == 0 && chemical(ChemicalID::Ethylene) > g.ethylene_abscission_threshold) {
        senescence_ticks = 1;
    }
    if (senescence_ticks > 0) {
        senescence_ticks++;
        if (senescence_ticks >= g.senescence_duration) {
            die(plant);
            return;
        }
    }
}
```

- [ ] **Step 2: Remove process_abscission call from Plant::tick()**

In `src/engine/plant.cpp`:

Remove: `process_abscission(*this);`
Remove: `#include "engine/ethylene.h"` (if no longer needed — `compute_ethylene` is still called)

Check: `compute_ethylene(*this, world)` is still called in Plant::tick(). Keep that include if needed.

- [ ] **Step 3: Update abscission tests**

The tests for abscission previously called `process_abscission(plant)` directly. Now abscission happens during `leaf->tick()`. Update tests:

```cpp
TEST_CASE("Leaf above ethylene threshold begins senescence", "[ethylene][abscission]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->chemical(ChemicalID::Ethylene) = g.ethylene_abscission_threshold + 0.1f;
    leaf->chemical(ChemicalID::Sugar) = 1.0f;  // enough to survive tick
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    leaf->tick(plant, wp);

    REQUIRE(leaf->as_leaf()->senescence_ticks > 0);
}

TEST_CASE("Leaf below ethylene threshold stays healthy", "[ethylene][abscission]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->chemical(ChemicalID::Ethylene) = g.ethylene_abscission_threshold * 0.5f;
    leaf->chemical(ChemicalID::Sugar) = 1.0f;
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    leaf->tick(plant, wp);

    REQUIRE(leaf->as_leaf()->senescence_ticks == 0);
}

TEST_CASE("Senescing leaf is removed after senescence_duration", "[ethylene][abscission]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->chemical(ChemicalID::Ethylene) = g.ethylene_abscission_threshold + 0.1f;
    leaf->as_leaf()->senescence_ticks = g.senescence_duration - 2;
    leaf->chemical(ChemicalID::Sugar) = 1.0f;
    plant.seed_mut()->add_child(leaf);

    uint32_t count_before = plant.node_count();
    WorldParams wp = default_world_params();
    leaf->tick(plant, wp);

    REQUIRE(plant.node_count() < count_before);
}
```

- [ ] **Step 4: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`
Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/engine/node/leaf_node.h src/engine/node/leaf_node.cpp src/engine/plant.cpp tests/test_ethylene.cpp
git commit -m "refactor: move abscission logic from ethylene.cpp into LeafNode::tick()"
```

---

### Task 9: Delete old files and update CMakeLists

**Files:**
- Delete: `src/engine/hormone.h`, `src/engine/hormone.cpp`
- Delete: `src/engine/gibberellin.h`, `src/engine/gibberellin.cpp`
- Modify: `src/engine/sugar.h` (keep only `sugar_cap` and `compute_light_exposure`)
- Modify: `src/engine/sugar.cpp` (remove dead code: `grow_leaves`, `transport_sugar`)
- Modify: `src/engine/ethylene.h` (keep only `compute_ethylene`)
- Modify: `src/engine/ethylene.cpp` (remove `process_abscission`)
- Modify: `src/engine/plant.cpp` (remove dead includes)
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Delete hormone files**

```bash
git rm src/engine/hormone.h src/engine/hormone.cpp
```

Remove any `#include "engine/hormone.h"` from `src/engine/plant.cpp` or other files.

- [ ] **Step 2: Delete gibberellin files**

```bash
git rm src/engine/gibberellin.h src/engine/gibberellin.cpp
```

Remove any remaining `#include "engine/gibberellin.h"` (plant.cpp was already cleaned in Task 7).

- [ ] **Step 3: Clean up sugar files**

In `src/engine/sugar.h`, remove:
- `void transport_sugar(Plant& plant, const WorldParams& world);`

Keep:
- `float sugar_cap(const Node& node, const Genome& g);`
- `void compute_light_exposure(Plant& plant, const WorldParams& world);`

In `src/engine/sugar.cpp`, delete `grow_leaves()` entirely (dead code — LeafNode::grow() handles this). Delete `transport_sugar()` — it was just calling `compute_light_exposure()` and a comment.

In `src/engine/plant.cpp`, replace `transport_sugar(*this, world)` with `compute_light_exposure(*this, world)` (since that's all transport_sugar does now). Add `#include "engine/sugar.h"` if not already present.

- [ ] **Step 4: Clean up ethylene files**

In `src/engine/ethylene.h`, remove `void process_abscission(Plant& plant);`

In `src/engine/ethylene.cpp`, delete the `process_abscission()` function entirely.

- [ ] **Step 5: Update CMakeLists.txt**

In `CMakeLists.txt`, in the `botany_engine` source list:

Remove: `src/engine/hormone.cpp`
Remove: `src/engine/gibberellin.cpp`

Keep: `src/engine/sugar.cpp` (still has sugar_cap + compute_light_exposure)
Keep: `src/engine/ethylene.cpp` (still has compute_ethylene)

- [ ] **Step 6: Clean up remaining includes**

Grep for any remaining includes of deleted files:
```bash
grep -rn 'engine/hormone\.h\|engine/gibberellin\.h' src/ tests/
```
Remove any found.

- [ ] **Step 7: Build and run ALL tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests`
Expected: All tests pass. No dead code remains.

- [ ] **Step 8: Commit**

```bash
git add -A
git commit -m "cleanup: delete old chemical files, remove dead code"
```

---

### Task 10: Final verification and integration test

- [ ] **Step 1: Run full test suite**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests --reporter compact`
Expected: All tests pass.

- [ ] **Step 2: Run the realtime viewer**

Run: `./build/botany_realtime --color auxin`
Expected: Plants grow normally. Auxin heatmap displays correctly.

Run: `./build/botany_realtime --color sugar`
Expected: Sugar heatmap displays correctly.

Run: `./build/botany_realtime --color ethylene`
Expected: Ethylene heatmap displays correctly. Leaves still senesce and drop.

- [ ] **Step 3: Verify adding-a-chemical checklist**

Confirm that adding a new hormone (hypothetical) requires only:
1. One entry in `ChemicalID` enum
2. One def header in `chemical/hormone/`
3. One line in `chemical_registry.h` → `hormone_params()`
4. Genome fields for evolvable params
5. Node type `tick()` methods for production/consumption

No changes to transport code, Node base class, or chemical infrastructure.

- [ ] **Step 4: Final commit tag**

```bash
git tag chemical-class-system
```
