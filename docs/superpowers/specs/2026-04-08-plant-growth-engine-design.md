# Plant Growth Simulation Engine — Design Spec

A biologically-grounded plant growth simulator driven by hormone transport (auxin/cytokinin) through a graph-based plant model, with a separate OpenGL+GLFW 3D renderer.

## Goals

- Simulate a single plant growing from seed to maturity using hormone-driven meristem behavior
- Engine is instance-based and renderer-independent — ready for multi-plant evolution runs later
- Two runtime modes: real-time (engine + renderer) and headless precompute with playback
- C++ with CMake. OpenGL + GLFW for rendering.

## Core Data Model — The Plant Graph

The plant is a directed tree graph of **Node**s. Each node is a point in 3D space where biological activity happens — hormones accumulate, meristems grow, leaves attach. The physical "segment" between two nodes is not a first-class object in the engine; it's implied by the parent-child connection and derived by the renderer.

```cpp
struct Node {
    uint id;
    Node* parent;
    vector<Node*> children;

    // Geometry (engine-owned, renderer reads)
    vec3 position;       // world position of this node
    float radius;        // thickness at this point

    // Biology
    NodeType type;       // STEM, ROOT
    uint age;            // ticks since creation
    float auxin;         // local auxin concentration
    float cytokinin;     // local cytokinin concentration

    // Meristem (null if this node is not an active growth point)
    Meristem* meristem;

    // Leaf (null if no leaf at this node)
    Leaf* leaf;
};

enum NodeType { STEM, ROOT };

struct Meristem {
    MeristemType type;           // APICAL, AXILLARY, ROOT_APICAL, ROOT_AXILLARY
    bool active;                 // axillary meristems start dormant
    uint ticks_since_last_node;  // for internode spacing
};

struct Leaf {
    float size;
};
```

The node graph is both the plant's physical body and its hormone transport network. The renderer draws a tapered cylinder between each node and its parent — bottom radius from parent node, top radius from child node, length from the distance between their positions.

## Hormone Transport Model

Each tick, hormones are transported through the node graph in two passes.

### Pass 1 — Auxin (downward, tips to roots)

Post-order traversal (leaf nodes inward). Each node:

1. **Production:** If it has an active apical meristem, produce auxin at `genome.auxin_production_rate`
2. **Transport:** Send fraction of local auxin to parent node: `flow = auxin * genome.auxin_transport_rate`
3. **Decay:** `auxin *= (1.0 - genome.auxin_decay_rate)`

### Pass 2 — Cytokinin (upward, roots to tips)

Pre-order traversal (root nodes outward). Each node:

1. **Production:** If it has an active root apical meristem, produce cytokinin at `genome.cytokinin_production_rate`
2. **Transport:** Send fraction to each child node, split proportionally: `flow_per_child = cytokinin * genome.cytokinin_transport_rate / num_children`
3. **Decay:** `cytokinin *= (1.0 - genome.cytokinin_decay_rate)`

### Emergent Behaviors

- **Apical dominance:** Topmost apical meristem floods stem with auxin, keeping axillary meristems dormant
- **Decapitation response:** If apical meristem is removed/stops, auxin drops, nearest axillary meristems activate (like pruning)
- **Root-shoot balance:** More roots = more cytokinin = more shoot growth potential. More shoots = more auxin = suppression of lateral branching.
- **Branching gradient:** Axillary meristems far from apex get less auxin, so lower branches activate before upper ones

## Meristem Tick Behavior

### Shoot Apical Meristem (stem tips)

1. **Extend** — Move this node's position outward along its growth direction by `genome.growth_rate`. This increases the distance from its parent node.
2. **Thicken** — `node.radius += genome.thickening_rate`
3. **Produce auxin** — (handled in hormone pass)
4. **Node spacing check** — Increment `ticks_since_last_node`. When it reaches `genome.internode_spacing`:
   - Spawn a new child node branching off at `genome.branch_angle`
   - Attach a dormant axillary meristem to that child node
   - That child node also gets a leaf (size from `genome.leaf_size`)
   - Reset counter
5. **Chain growth** — When the distance between this node and its parent exceeds `genome.max_internode_length`, insert a new intermediate node: create a child node at this node's current position with a fresh apical meristem, and this node becomes a fixed interior node (meristem removed). The stem grows as a chain of nodes.

### Shoot Axillary Meristem (branch junctions)

1. **Check activation** — If dormant, test: `local_auxin < genome.auxin_threshold AND local_cytokinin > genome.cytokinin_threshold`
2. **On activation** — Convert to an apical meristem. It now follows the apical tick cycle, growing a new branch outward from this node.

### Root Apical Meristem (root tips)

Same architecture as shoot apical meristem, mirrored:

1. **Extend** — Move this node's position along its growth direction by `genome.root_growth_rate`
2. **Produce cytokinin** — (handled in hormone pass)
3. **Node spacing check** — At intervals (`genome.root_internode_spacing`), spawn a new child node with a dormant root axillary meristem

### Root Axillary Meristem (root branch junctions)

Same pattern as shoot axillary, driven by cytokinin instead of auxin:

1. **Check activation** — Dormant until `local_cytokinin < genome.cytokinin_threshold`
2. **On activation** — Converts to a root apical meristem, grows a new root branch from this node

### Symmetry

The plant runs the same branching algorithm twice:
- **Shoots:** Apical produces auxin → suppresses axillary → auxin drops → branch activates
- **Roots:** Root apical produces cytokinin → suppresses root axillary → cytokinin drops → root branch activates

## Genome

A struct of float parameters — the fractal constants that define a species. Every meristem reads from the same genome.

```cpp
struct Genome {
    // Hormone production & sensitivity
    float auxin_production_rate;
    float auxin_transport_rate;
    float auxin_decay_rate;
    float auxin_threshold;         // below this, shoot axillary meristems can activate

    float cytokinin_production_rate;
    float cytokinin_transport_rate;
    float cytokinin_decay_rate;
    float cytokinin_threshold;     // below this, root axillary meristems can activate

    // Shoot growth
    float growth_rate;
    float max_internode_length;    // max distance between nodes before chain extends
    uint internode_spacing;        // ticks between axillary node spawns
    float branch_angle;
    float thickening_rate;

    // Root growth
    float root_growth_rate;
    float root_max_internode_length;
    uint root_internode_spacing;
    float root_branch_angle;

    // Geometry
    float leaf_size;
    float initial_radius;
};
```

Morphological diversity comes from these parameters interacting with the transport simulation. High auxin production + high threshold = tall columnar plant. Low values = bushy. Short internode spacing = dense nodes. Wide branch angle = spreading canopy.

## Seed Initialization (Tick 0)

When `create_plant()` is called, the engine creates the initial graph:

1. A single **seed node** at the given position with `genome.initial_radius`
2. Two child nodes created from the seed:
   - A **shoot node** with an active shoot apical meristem (growth direction = up)
   - A **root node** with an active root apical meristem (growth direction = down)
3. Both start extending on tick 1. The plant grows upward and downward simultaneously from the seed point.

## Engine Architecture

The engine is a standalone library with no rendering dependencies.

```cpp
class Engine {
public:
    PlantID create_plant(const Genome& genome, vec3 position);
    void tick();
    const Plant& get_plant(PlantID id) const;
    uint get_tick() const;
    void save_state(Writer& writer);
    void load_state(Reader& reader);
};
```

### Mode 1 — Real-time

Engine + renderer in same process. Main loop: `tick()` → renderer reads plant graph → draws frame → repeat. Speed controlled by ticks-per-frame. Renderer has strictly read-only access to the plant graph.

### Mode 2 — Headless precompute

No renderer. Engine runs N ticks, serializing state at each tick (or at configurable intervals) to a binary file.

Binary format:
```
Header:   { num_ticks, num_plants, genome_data }
Per tick: { tick_number, num_nodes, [node_data...] }
Node:     { id, parent_id, type, position, radius,
            auxin, cytokinin, has_leaf, leaf_size }
```

### Mode 2 — Playback

Separate viewer app loads the binary file and scrubs through ticks. Slider to jump to any tick. Play/pause/speed controls. No simulation — just reading precomputed frames.

## Renderer

Reads the plant graph, generates geometry, draws it. Never modifies the engine.

### Segment geometry
For each parent-child node pair, generate a tapered cylinder:
- Bottom at parent node's position, top at child node's position
- Bottom radius from parent node, top radius from child node
- Length derived from distance between the two node positions
- Subdivisions configurable (8 for performance, 16 for quality)

### Leaves
Flat quad (two triangles) oriented outward from the branch, sized by leaf size.

### Camera
Orbit camera — mouse drag rotates, scroll zooms. Target starts at seed position, re-centerable as plant grows.

### Lighting
Directional light (sun from above) + ambient.

### Ground plane
Flat plane at y=0. Roots below, shoots above.

### Coloring
- Stems/branches: brown (varying by age)
- Roots: darker brown
- Leaves: green
- Optional toggle: color by hormone concentration (red = high auxin, blue = high cytokinin)

### Playback UI
Tick counter, play/pause (spacebar), speed (arrow keys), scrub slider. ImGui for quick dev UI.

## Project Structure

```
botany/
├── CMakeLists.txt
├── src/
│   ├── engine/           # Pure simulation, no graphics
│   │   ├── engine.h/cpp
│   │   ├── plant.h/cpp
│   │   ├── node.h/cpp
│   │   ├── meristem.h/cpp
│   │   ├── genome.h/cpp
│   │   └── hormone.h/cpp
│   ├── renderer/         # OpenGL + GLFW, reads engine data
│   │   ├── renderer.h/cpp
│   │   ├── camera.h/cpp
│   │   └── plant_mesh.h/cpp
│   ├── serialization/    # Binary read/write
│   │   └── serializer.h/cpp
│   ├── app_realtime.cpp  # Mode 1 entry point
│   ├── app_headless.cpp  # Mode 2 precompute entry point
│   └── app_playback.cpp  # Mode 2 viewer entry point
```

Engine is a static library. The three apps link against it. Renderer only links into `app_realtime` and `app_playback`. The evolution library can link against just the engine and ignore rendering.
