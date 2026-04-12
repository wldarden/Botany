# Mass & Stress System Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add gravitational physics — each above-ground node computes mass, accumulates subtree mass/moment, computes structural stress, and responds via a stress hormone that modulates thickening, elongation, and gravitropism. Branches droop under moderate load and snap under severe load.

**Architecture:** Mass and stress are computed locally during each node's tick (DFS pre-order, one tick stale for child values). A new stress hormone (`ChemicalID::Stress`) uses the existing chemical transport system. Droop/break are checked after stress computation. Nine new genome parameters in a `stress` linkage group. Seven new WorldParams constants.

**Tech Stack:** C++17, glm, Catch2

**Build:** `/usr/local/bin/cmake --build build`

**Test:** `./build/botany_tests`

---

## File Structure

**Modified files:**
- `src/engine/genome.h` — 9 new genome fields + defaults
- `src/engine/world_params.h` — 7 new WorldParams constants
- `src/engine/chemical/chemical.h` — add `ChemicalID::Stress`
- `src/engine/chemical/chemical_registry.h` — register stress in `all_chemical_ids` and `diffusion_params`
- `src/engine/node/node.h` — add `total_mass`, `mass_moment`, `stress` fields
- `src/engine/node/node.cpp` — mass/stress computation + stress hormone production + droop/break in `tick()`
- `src/engine/node/stem_node.cpp` — stress hormone effects on thickening/elongation
- `src/engine/node/meristems/shoot_apical.cpp` — gravitropism boost from stress hormone
- `src/evolution/genome_bridge.cpp` — register 9 new genes + stress linkage group
- `src/renderer/renderer.cpp` — stress heatmap color mode
- `src/app_realtime.cpp` — stress overlay button
- `tests/test_evolution.cpp` — update round-trip test for new genome fields

**No new files.**

---

### Task 1: Genome + WorldParams — add stress parameters

**Files:**
- Modify: `src/engine/genome.h`
- Modify: `src/engine/world_params.h`

- [ ] **Step 1: Add stress fields to Genome**

In `src/engine/genome.h`, add after the ethylene block (before the closing `};` of the struct):

```cpp
    // Stress — mechanical load response
    float wood_density;                   // g/dm³ — mass per volume, also determines strength
    float wood_flexibility;               // 0-1 — droop threshold as fraction of break threshold
    float stress_hormone_production_rate; // hormone per unit stress
    float stress_hormone_diffusion_rate;  // fraction diffused per tick
    float stress_hormone_decay_rate;      // fraction decayed per tick
    float stress_thickening_boost;        // thickening multiplier per unit stress hormone
    float stress_elongation_inhibition;   // elongation suppression per unit stress hormone
    float stress_gravitropism_boost;      // gravitropism pull per unit stress hormone
```

Add defaults in `default_genome()` after the ethylene block:

```cpp
        // Stress
        .wood_density = 50.0f,                    // g/dm³ — light deciduous wood
        .wood_flexibility = 0.5f,                 // droop starts at 50% of break stress
        .stress_hormone_production_rate = 0.1f,   // moderate signaling
        .stress_hormone_diffusion_rate = 0.15f,   // moderate local diffusion
        .stress_hormone_decay_rate = 0.2f,        // fades quickly — local signal
        .stress_thickening_boost = 1.0f,          // 1:1 hormone-to-thickening boost
        .stress_elongation_inhibition = 1.0f,     // 1:1 hormone-to-elongation suppression
        .stress_gravitropism_boost = 0.5f,        // moderate upward correction
```

- [ ] **Step 2: Add stress constants to WorldParams**

In `src/engine/world_params.h`, add to the `WorldParams` struct after `light_direction`:

```cpp
    // Stress physics
    float gravity = 9.81f;                  // m/s² — gravitational acceleration
    float break_strength_factor = 5.0f;     // stress units per (g/dm³) of wood density
    float reference_wood_density = 50.0f;   // g/dm³ — density at which sugar costs are calibrated
    float leaf_mass_density = 5.0f;         // g/dm² of leaf area
    float meristem_mass = 0.1f;             // g — fixed mass for meristem tips
    float ground_support_height = 0.5f;     // dm — below this Y, stress is zeroed
    float droop_rate = 0.01f;               // radians/tick — max angular droop when overstressed
```

Add matching values in `default_world_params()`.

- [ ] **Step 3: Build and verify**

Run: `/usr/local/bin/cmake --build build`
Expected: Build succeeds (new fields have defaults, no code references them yet).

- [ ] **Step 4: Commit**

```bash
git add src/engine/genome.h src/engine/world_params.h
git commit -m "feat(stress): add genome and WorldParams fields for mass/stress system"
```

---

### Task 2: Chemical registry — add Stress chemical

**Files:**
- Modify: `src/engine/chemical/chemical.h`
- Modify: `src/engine/chemical/chemical_registry.h`
- Modify: `src/engine/node/node.cpp`

- [ ] **Step 1: Add Stress to ChemicalID enum**

In `src/engine/chemical/chemical.h`, add `Stress` to the enum:

```cpp
enum class ChemicalID : uint8_t {
    Auxin,
    Cytokinin,
    Gibberellin,
    Sugar,
    Ethylene,
    Stress,
};
```

- [ ] **Step 2: Register in chemical_registry.h**

In `src/engine/chemical/chemical_registry.h`, update `all_chemical_ids` to include `Stress`:

```cpp
inline constexpr std::array<ChemicalID, 6> all_chemical_ids = {
    ChemicalID::Auxin,
    ChemicalID::Cytokinin,
    ChemicalID::Gibberellin,
    ChemicalID::Sugar,
    ChemicalID::Ethylene,
    ChemicalID::Stress,
};
```

Update `diffusion_params` to return 5 entries (add stress):

```cpp
inline std::array<ChemicalDiffusionParams, 5> diffusion_params(const Genome& g) {
    return {{
        {ChemicalID::Auxin,       g.auxin_diffusion_rate,     g.auxin_decay_rate},
        {ChemicalID::Cytokinin,   g.cytokinin_diffusion_rate, g.cytokinin_decay_rate},
        {ChemicalID::Gibberellin, g.ga_diffusion_rate,        g.ga_decay_rate},
        {ChemicalID::Sugar,       g.sugar_diffusion_rate,     0.0f},
        {ChemicalID::Stress,      g.stress_hormone_diffusion_rate, g.stress_hormone_decay_rate},
    }};
}
```

- [ ] **Step 3: Initialize Stress chemical in Node constructor**

In `src/engine/node/node.cpp`, in the `Node` constructor, add after the Ethylene init:

```cpp
    chemicals[ChemicalID::Stress] = 0.0f;
```

- [ ] **Step 4: Build and verify**

Run: `/usr/local/bin/cmake --build build`
Expected: Build succeeds. Stress hormone now diffuses and decays through the existing transport system.

- [ ] **Step 5: Commit**

```bash
git add src/engine/chemical/chemical.h src/engine/chemical/chemical_registry.h src/engine/node/node.cpp
git commit -m "feat(stress): register Stress chemical in transport system"
```

---

### Task 3: Node fields — add mass/moment/stress

**Files:**
- Modify: `src/engine/node/node.h`
- Modify: `src/engine/node/node.cpp`
- Test: `tests/test_evolution.cpp`

- [ ] **Step 1: Write the failing test**

In `tests/test_evolution.cpp`, add:

```cpp
TEST_CASE("Node has mass and stress fields", "[stress]") {
    botany::Engine engine;
    botany::Genome g = botany::default_genome();
    engine.create_plant(g, glm::vec3(0.0f));

    // Run a few ticks so nodes exist and have been processed
    for (int i = 0; i < 50; i++) engine.tick();

    const botany::Plant& plant = engine.get_plant(0);
    bool found_nonzero_mass = false;
    plant.for_each_node([&](const botany::Node& n) {
        // total_mass should be positive for any node that exists
        if (n.total_mass > 0.0f) found_nonzero_mass = true;
    });
    REQUIRE(found_nonzero_mass);
}
```

- [ ] **Step 2: Build and verify it fails**

Run: `/usr/local/bin/cmake --build build`
Expected: Compile error — `total_mass` not a member of `Node`.

- [ ] **Step 3: Add fields to Node**

In `src/engine/node/node.h`, add to the public section after `uint32_t starvation_ticks = 0;`:

```cpp
    // --- Mass / stress (computed each tick, children's values are one tick stale) ---
    float total_mass = 0.0f;       // self mass + all children's total_mass
    glm::vec3 mass_moment{0.0f};   // self_mass * position + Σ child.mass_moment
    float stress = 0.0f;           // torque / cross-section (structural load)
```

- [ ] **Step 4: Add mass/stress computation to Node::tick()**

In `src/engine/node/node.cpp`, in `Node::tick()`, add after the `grow(plant, world);` call and before `transport_chemicals(g);`:

```cpp
    // --- Mass & stress computation ---
    // Self mass depends on node type
    float self_mass = 0.0f;
    bool is_underground = (type == NodeType::ROOT || type == NodeType::ROOT_APICAL || type == NodeType::ROOT_AXILLARY);

    if (!is_underground) {
        if (type == NodeType::LEAF) {
            auto* leaf = as_leaf();
            self_mass = leaf ? (leaf->leaf_size * leaf->leaf_size * world.leaf_mass_density) : 0.0f;
        } else if (type == NodeType::STEM) {
            float length = std::max(glm::length(offset), 0.01f);
            self_mass = 3.14159f * radius * radius * length * g.wood_density;
        } else if (is_meristem()) {
            self_mass = world.meristem_mass;
        }
    }

    // Accumulate subtree mass from direct children (their values are one tick stale)
    total_mass = self_mass;
    mass_moment = self_mass * position;
    for (const Node* child : children) {
        total_mass += child->total_mass;
        mass_moment += child->mass_moment;
    }

    // Stress computation (above-ground only, skip if near ground)
    stress = 0.0f;
    if (!is_underground && position.y > world.ground_support_height) {
        float child_mass = total_mass - self_mass;
        if (child_mass > 1e-6f && radius > 1e-6f) {
            glm::vec3 child_com = (mass_moment - self_mass * position) / child_mass;
            float dx = child_com.x - position.x;
            float dz = child_com.z - position.z;
            float lever_arm = std::sqrt(dx * dx + dz * dz);
            float torque = child_mass * world.gravity * lever_arm;
            float cross_section = 3.14159f * radius * radius;
            stress = torque / cross_section;
        }
    }

    // Stress hormone production (above-ground only)
    if (!is_underground && stress > 0.0f) {
        chemical(ChemicalID::Stress) += stress * g.stress_hormone_production_rate;
    }
```

- [ ] **Step 5: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[stress]"`
Expected: Test passes — nodes have positive total_mass after ticks.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/node.h src/engine/node/node.cpp tests/test_evolution.cpp
git commit -m "feat(stress): compute mass, moment, and stress per node during tick"
```

---

### Task 4: Droop and break

**Files:**
- Modify: `src/engine/node/node.cpp`
- Test: `tests/test_evolution.cpp`

- [ ] **Step 1: Write the failing test**

In `tests/test_evolution.cpp`, add:

```cpp
TEST_CASE("Extreme stress causes branch break", "[stress]") {
    botany::Genome g = botany::default_genome();
    g.wood_density = 200.0f;         // very heavy wood
    g.wood_flexibility = 0.1f;       // very rigid (snaps easily)
    g.growth_rate = 0.05f;           // fast growth to build mass quickly
    g.branch_angle = 1.4f;           // nearly horizontal branches

    botany::Engine engine;
    botany::WorldParams world = botany::default_world_params();
    world.break_strength_factor = 0.5f;  // very weak wood — easy to break
    engine.world_params_mut() = world;
    engine.create_plant(g, glm::vec3(0.0f));

    uint32_t initial_nodes = 0;
    uint32_t final_nodes = 0;
    bool had_nodes = false;

    for (int i = 0; i < 2000; i++) {
        engine.tick();
        uint32_t count = engine.get_plant(0).node_count();
        if (count > 10 && !had_nodes) {
            initial_nodes = count;
            had_nodes = true;
        }
        final_nodes = count;
    }

    // With very heavy, rigid, horizontal branches and weak wood,
    // some branches should have broken off
    if (had_nodes) {
        REQUIRE(final_nodes < initial_nodes);
    }
}
```

- [ ] **Step 2: Build and verify it fails**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[stress]"`
Expected: Test fails — no branch breaking implemented yet.

- [ ] **Step 3: Add droop and break logic to Node::tick()**

In `src/engine/node/node.cpp`, in `Node::tick()`, add after the stress hormone production block (and before `transport_chemicals`):

```cpp
    // --- Droop and break (above-ground stems only) ---
    if (type == NodeType::STEM && !is_underground && stress > 0.0f) {
        float break_stress = g.wood_density * world.break_strength_factor;
        float droop_threshold = break_stress * g.wood_flexibility;

        if (stress >= break_stress) {
            // Branch snaps — remove this node and entire subtree
            die(plant);
            return;  // node is dead, skip transport
        }

        if (stress > droop_threshold) {
            // Branch droops toward gravity
            float excess = (stress - droop_threshold) / (break_stress - droop_threshold);
            float droop_angle = std::min(excess * world.droop_rate, world.droop_rate);
            float len = glm::length(offset);
            if (len > 1e-4f) {
                glm::vec3 dir = offset / len;
                glm::vec3 down(0.0f, -1.0f, 0.0f);
                // Rotate dir toward down by droop_angle
                glm::vec3 axis = glm::cross(dir, down);
                float axis_len = glm::length(axis);
                if (axis_len > 1e-6f) {
                    axis /= axis_len;
                    float c = std::cos(droop_angle);
                    float s = std::sin(droop_angle);
                    glm::vec3 new_dir = dir * c
                        + glm::cross(axis, dir) * s
                        + axis * glm::dot(axis, dir) * (1.0f - c);
                    offset = glm::normalize(new_dir) * len;
                }
            }
        }
    }
```

- [ ] **Step 4: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[stress]"`
Expected: Both stress tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/engine/node/node.cpp tests/test_evolution.cpp
git commit -m "feat(stress): add droop and break response to structural stress"
```

---

### Task 5: Stress hormone effects on growth

**Files:**
- Modify: `src/engine/node/stem_node.cpp`
- Modify: `src/engine/node/meristems/shoot_apical.cpp`

- [ ] **Step 1: Apply stress hormone to thickening in StemNode**

In `src/engine/node/stem_node.cpp`, in `StemNode::thicken()`, add stress hormone boost. Replace:

```cpp
    float actual_rate = g.thickening_rate * gf;
```

With:

```cpp
    float stress_boost = 1.0f + chemical(ChemicalID::Stress) * g.stress_thickening_boost;
    float actual_rate = g.thickening_rate * gf * stress_boost;
```

- [ ] **Step 2: Apply stress hormone to elongation in StemNode**

In `src/engine/node/stem_node.cpp`, in `StemNode::elongate()`, add stress hormone inhibition. After the existing `eth_inhibit` line:

```cpp
    float eth_inhibit = std::max(0.0f, 1.0f - chemical(ChemicalID::Ethylene) * g.ethylene_elongation_inhibition);
```

Add:

```cpp
    float stress_inhibit = std::max(0.0f, 1.0f - chemical(ChemicalID::Stress) * g.stress_elongation_inhibition);
```

And update the `effective_rate` line to include it:

```cpp
    float effective_rate = g.internode_elongation_rate * ga_boost * eth_inhibit * stress_inhibit;
```

- [ ] **Step 3: Apply stress hormone to shoot gravitropism**

In `src/engine/node/meristems/shoot_apical.cpp`, in `ShootApicalNode::grow_tip()`, after the direction roll, add a vertical pull from stress hormone. After:

```cpp
    if (glm::length(growth_dir) < 1e-4f) roll_direction(g);
```

Add:

```cpp
    // Stress hormone pulls growth direction toward vertical
    float stress_grav = chemical(ChemicalID::Stress) * g.stress_gravitropism_boost;
    if (stress_grav > 1e-6f) {
        glm::vec3 up(0.0f, 1.0f, 0.0f);
        float blend = std::min(stress_grav, 0.5f);  // cap at 50% pull toward vertical
        growth_dir = glm::normalize(glm::mix(growth_dir, up, blend));
    }
```

- [ ] **Step 4: Apply sugar cost scaling with wood density**

In `src/engine/node/stem_node.cpp`, in both `thicken()` and `elongate()`, scale costs by wood density. In `thicken()`, replace:

```cpp
    float max_cost = g.thickening_rate * world.sugar_cost_thickening;
```

With:

```cpp
    float density_scale = g.wood_density / world.reference_wood_density;
    float max_cost = g.thickening_rate * world.sugar_cost_thickening * density_scale;
```

In `elongate()`, replace:

```cpp
    float max_cost = effective_rate * world.sugar_cost_elongation;
```

With:

```cpp
    float density_scale = g.wood_density / world.reference_wood_density;
    float max_cost = effective_rate * world.sugar_cost_elongation * density_scale;
```

Also update the sugar deduction lines to use the scaled cost (they should already reference `max_cost` or `world.sugar_cost_*` — make sure the actual deduction uses the density-scaled value).

In `thicken()`, the deduction line:
```cpp
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_thickening;
```
becomes:
```cpp
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_thickening * density_scale;
```

In `elongate()`, the deduction line:
```cpp
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_elongation;
```
becomes:
```cpp
    chemical(ChemicalID::Sugar) -= actual_rate * world.sugar_cost_elongation * density_scale;
```

- [ ] **Step 5: Build and run all tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[stress]"`
Expected: All stress tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/stem_node.cpp src/engine/node/meristems/shoot_apical.cpp
git commit -m "feat(stress): stress hormone modulates thickening, elongation, and gravitropism"
```

---

### Task 6: Evolution genome bridge — register stress genes

**Files:**
- Modify: `src/evolution/genome_bridge.cpp`
- Test: `tests/test_evolution.cpp`

- [ ] **Step 1: Register stress genes in build_genome_template**

In `src/evolution/genome_bridge.cpp`, in `build_genome_template()`, add after the ethylene block and before the linkage groups:

```cpp
    // --- Stress group (8 genes) ---
    reg(sg, "wood_density",                    g.wood_density,                    r, 10.0f, 200.0f, p);
    reg(sg, "wood_flexibility",                g.wood_flexibility,                r, 0.1f,  1.0f, p);
    reg(sg, "stress_hormone_production_rate",  g.stress_hormone_production_rate,  r, 0.0f,  1.0f, p);
    reg(sg, "stress_hormone_diffusion_rate",   g.stress_hormone_diffusion_rate,   r, 0.01f, 0.5f, p);
    reg(sg, "stress_hormone_decay_rate",       g.stress_hormone_decay_rate,       r, 0.01f, 0.5f, p);
    reg(sg, "stress_thickening_boost",         g.stress_thickening_boost,         r, 0.0f,  5.0f, p);
    reg(sg, "stress_elongation_inhibition",    g.stress_elongation_inhibition,    r, 0.0f,  5.0f, p);
    reg(sg, "stress_gravitropism_boost",       g.stress_gravitropism_boost,       r, 0.0f,  5.0f, p);
```

Add the stress linkage group after the ethylene linkage group:

```cpp
    sg.add_linkage_group({"stress", {
        "wood_density", "wood_flexibility",
        "stress_hormone_production_rate", "stress_hormone_diffusion_rate",
        "stress_hormone_decay_rate",
        "stress_thickening_boost", "stress_elongation_inhibition", "stress_gravitropism_boost"
    }});
```

- [ ] **Step 2: Add from_structured readback**

In `src/evolution/genome_bridge.cpp`, in `from_structured()`, add after the ethylene block:

```cpp
    // Stress
    g.wood_density                    = sg.get("wood_density");
    g.wood_flexibility                = sg.get("wood_flexibility");
    g.stress_hormone_production_rate  = sg.get("stress_hormone_production_rate");
    g.stress_hormone_diffusion_rate   = sg.get("stress_hormone_diffusion_rate");
    g.stress_hormone_decay_rate       = sg.get("stress_hormone_decay_rate");
    g.stress_thickening_boost         = sg.get("stress_thickening_boost");
    g.stress_elongation_inhibition    = sg.get("stress_elongation_inhibition");
    g.stress_gravitropism_boost       = sg.get("stress_gravitropism_boost");
```

- [ ] **Step 3: Update round-trip test**

In `tests/test_evolution.cpp`, in the "Genome round-trips through StructuredGenome" test, add a spot-check:

```cpp
    REQUIRE_THAT(restored.wood_density,
                 WithinAbs(original.wood_density, 1e-6));
    REQUIRE_THAT(restored.stress_thickening_boost,
                 WithinAbs(original.stress_thickening_boost, 1e-6));
```

Also update the linkage group count test from 8 to 9:

```cpp
    REQUIRE(groups.size() == 9);
```

- [ ] **Step 4: Build and run tests**

Run: `/usr/local/bin/cmake --build build && ./build/botany_tests "[evolution]"`
Expected: All evolution tests pass with 9 linkage groups.

- [ ] **Step 5: Commit**

```bash
git add src/evolution/genome_bridge.cpp tests/test_evolution.cpp
git commit -m "feat(stress): register stress genome parameters in evolution bridge"
```

---

### Task 7: Renderer — stress heatmap overlay

**Files:**
- Modify: `src/app_realtime.cpp`

- [ ] **Step 1: Add Stress overlay button**

In `src/app_realtime.cpp`, in the Overlays collapsing header section, after the Ethylene button block, add:

```cpp
            ImGui::SameLine();
            if (ImGui::Button("Stress")) {
                renderer.set_color_by_type(false);
                renderer.set_color_mode([](const Node& n) { return n.stress; });
                active_overlay = Overlay::STRESS;
            }
```

Also add `STRESS` to the `Overlay` enum (after `ETHYLENE`):

```cpp
    enum class Overlay { NONE, NODE_TYPE, AUXIN, CYTOKININ, SUGAR, LIGHT, GIBBERELLIN, ETHYLENE, STRESS };
```

- [ ] **Step 2: Build and verify**

Run: `/usr/local/bin/cmake --build build`
Expected: Build succeeds.

- [ ] **Step 3: Commit**

```bash
git add src/app_realtime.cpp
git commit -m "feat(stress): add stress heatmap overlay to realtime viewer"
```

---

### Task 8: Export — add stress fields to genome export

**Files:**
- Modify: `src/app_evolve.cpp`

- [ ] **Step 1: Add stress fields to save_genome**

In `src/app_evolve.cpp`, in the `save_genome()` function, add after the senescence_duration line:

```cpp
    out << "wood_density=" << g.wood_density << "\n";
    out << "wood_flexibility=" << g.wood_flexibility << "\n";
    out << "stress_hormone_production_rate=" << g.stress_hormone_production_rate << "\n";
    out << "stress_hormone_diffusion_rate=" << g.stress_hormone_diffusion_rate << "\n";
    out << "stress_hormone_decay_rate=" << g.stress_hormone_decay_rate << "\n";
    out << "stress_thickening_boost=" << g.stress_thickening_boost << "\n";
    out << "stress_elongation_inhibition=" << g.stress_elongation_inhibition << "\n";
    out << "stress_gravitropism_boost=" << g.stress_gravitropism_boost << "\n";
```

- [ ] **Step 2: Build and verify**

Run: `/usr/local/bin/cmake --build build`
Expected: Build succeeds.

- [ ] **Step 3: Commit**

```bash
git add src/app_evolve.cpp
git commit -m "feat(stress): add stress genome fields to export"
```
