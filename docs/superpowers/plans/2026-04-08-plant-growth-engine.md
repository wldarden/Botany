# Plant Growth Engine Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a hormone-driven plant growth simulator with a node-based graph model and an OpenGL 3D renderer.

**Architecture:** The engine is a standalone C++ static library that simulates plant growth through a tree graph of nodes. Meristems at node tips grow the graph; auxin/cytokinin flow through node connections to drive branching decisions. A separate renderer reads node data and draws tapered cylinders between parent-child pairs. Three executables: real-time viewer, headless precompute, and playback viewer.

**Tech Stack:** C++17, CMake, GLM (math), Catch2 (tests), GLFW (windowing), glad (OpenGL loader), Dear ImGui (playback UI)

**Spec:** `docs/superpowers/specs/2026-04-08-plant-growth-engine-design.md`

---

## File Structure

```
botany/
├── CMakeLists.txt
├── src/
│   ├── engine/
│   │   ├── genome.h              # Genome struct + default_genome()
│   │   ├── node.h                # Node, Meristem, Leaf structs + NodeType/MeristemType enums
│   │   ├── node.cpp              # Node constructor/destructor, add_child, tree utils
│   │   ├── plant.h               # Plant class — owns the node graph
│   │   ├── plant.cpp             # Plant creation (seed init), node access, traversal helpers
│   │   ├── hormone.h             # transport_auxin(), transport_cytokinin()
│   │   ├── hormone.cpp           # Hormone transport passes
│   │   ├── meristem.h            # tick_meristems()
│   │   ├── meristem.cpp          # Meristem tick behavior (apical, axillary, root)
│   │   ├── engine.h              # Engine class — tick loop, plant management
│   │   └── engine.cpp            # Engine implementation
│   ├── serialization/
│   │   ├── serializer.h          # save_recording(), load_recording(), Recording struct
│   │   └── serializer.cpp        # Binary read/write
│   ├── renderer/
│   │   ├── shader.h              # Shader compile/link helper
│   │   ├── shader.cpp
│   │   ├── camera.h              # OrbitCamera class
│   │   ├── camera.cpp
│   │   ├── renderer.h            # Renderer class — init, draw_plant, draw_ground
│   │   └── renderer.cpp
│   ├── app_realtime.cpp          # Mode 1: engine + renderer main loop
│   ├── app_headless.cpp          # Mode 2: precompute to binary file
│   └── app_playback.cpp          # Mode 2: load binary, scrub through ticks
├── tests/
│   ├── test_genome.cpp
│   ├── test_node.cpp
│   ├── test_plant.cpp
│   ├── test_hormone.cpp
│   ├── test_meristem.cpp
│   ├── test_engine.cpp
│   └── test_serializer.cpp
├── shaders/
│   ├── plant.vert
│   └── plant.frag
└── extern/                       # Third-party deps (fetched by CMake)
```

---

### Task 1: Project Scaffolding

**Files:**
- Create: `CMakeLists.txt`
- Create: `src/engine/genome.h`
- Create: `tests/test_genome.cpp`

This task sets up the CMake build, fetches dependencies, and creates the first compilable + testable file (Genome) to verify the toolchain works.

- [ ] **Step 1: Create CMakeLists.txt**

```cmake
cmake_minimum_required(VERSION 3.20)
project(botany LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Dependencies
include(FetchContent)

FetchContent_Declare(
    glm
    GIT_REPOSITORY https://github.com/g-truc/glm.git
    GIT_TAG 1.0.1
)

FetchContent_Declare(
    catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.2
)

FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG 3.4
)

FetchContent_Declare(
    glad
    GIT_REPOSITORY https://github.com/Dav1dde/glad.git
    GIT_TAG v2.0.6
    SOURCE_SUBDIR cmake
)

FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG v1.91.8
)

FetchContent_MakeAvailable(glm catch2 glfw)

# glad
FetchContent_GetProperties(glad)
if(NOT glad_POPULATED)
    FetchContent_Populate(glad)
endif()
add_subdirectory("${glad_SOURCE_DIR}/cmake" glad_cmake)
glad_add_library(glad_gl REPRODUCIBLE API gl:core=4.1)

# imgui - build as a static library
FetchContent_GetProperties(imgui)
if(NOT imgui_POPULATED)
    FetchContent_Populate(imgui)
endif()
add_library(imgui STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui PUBLIC ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends)
target_link_libraries(imgui PUBLIC glfw glad_gl)

# Engine library (no graphics dependencies)
add_library(botany_engine STATIC
    src/engine/node.cpp
    src/engine/plant.cpp
    src/engine/hormone.cpp
    src/engine/meristem.cpp
    src/engine/engine.cpp
)
target_include_directories(botany_engine PUBLIC src)
target_link_libraries(botany_engine PUBLIC glm::glm)

# Serialization library
add_library(botany_serialization STATIC
    src/serialization/serializer.cpp
)
target_include_directories(botany_serialization PUBLIC src)
target_link_libraries(botany_serialization PUBLIC botany_engine)

# Renderer library
add_library(botany_renderer STATIC
    src/renderer/shader.cpp
    src/renderer/camera.cpp
    src/renderer/renderer.cpp
)
target_include_directories(botany_renderer PUBLIC src)
target_link_libraries(botany_renderer PUBLIC botany_engine glad_gl glfw glm::glm)

# Executables
add_executable(botany_realtime src/app_realtime.cpp)
target_link_libraries(botany_realtime PRIVATE botany_engine botany_renderer)

add_executable(botany_headless src/app_headless.cpp)
target_link_libraries(botany_headless PRIVATE botany_engine botany_serialization)

add_executable(botany_playback src/app_playback.cpp)
target_link_libraries(botany_playback PRIVATE botany_engine botany_renderer botany_serialization imgui)

# Tests
enable_testing()
add_executable(botany_tests
    tests/test_genome.cpp
    tests/test_node.cpp
    tests/test_plant.cpp
    tests/test_hormone.cpp
    tests/test_meristem.cpp
    tests/test_engine.cpp
    tests/test_serializer.cpp
)
target_link_libraries(botany_tests PRIVATE botany_engine botany_serialization Catch2::Catch2WithMain)
include(CTest)
include(Catch)
catch_discover_tests(botany_tests)
```

- [ ] **Step 2: Create genome.h**

```cpp
// src/engine/genome.h
#pragma once

#include <cstdint>

namespace botany {

struct Genome {
    // Hormone production & sensitivity
    float auxin_production_rate;
    float auxin_transport_rate;
    float auxin_decay_rate;
    float auxin_threshold;

    float cytokinin_production_rate;
    float cytokinin_transport_rate;
    float cytokinin_decay_rate;
    float cytokinin_threshold;

    // Shoot growth
    float growth_rate;
    float max_internode_length;
    uint32_t internode_spacing;
    float branch_angle;
    float thickening_rate;

    // Root growth
    float root_growth_rate;
    float root_max_internode_length;
    uint32_t root_internode_spacing;
    float root_branch_angle;

    // Geometry
    float leaf_size;
    float initial_radius;
};

inline Genome default_genome() {
    return Genome{
        .auxin_production_rate = 1.0f,
        .auxin_transport_rate = 0.3f,
        .auxin_decay_rate = 0.05f,
        .auxin_threshold = 0.5f,

        .cytokinin_production_rate = 1.0f,
        .cytokinin_transport_rate = 0.3f,
        .cytokinin_decay_rate = 0.05f,
        .cytokinin_threshold = 0.5f,

        .growth_rate = 0.1f,
        .max_internode_length = 1.0f,
        .internode_spacing = 5,
        .branch_angle = 0.785f,  // ~45 degrees
        .thickening_rate = 0.001f,

        .root_growth_rate = 0.08f,
        .root_max_internode_length = 0.8f,
        .root_internode_spacing = 4,
        .root_branch_angle = 0.925f,  // ~53 degrees

        .leaf_size = 0.3f,
        .initial_radius = 0.05f,
    };
}

} // namespace botany
```

- [ ] **Step 3: Create test_genome.cpp**

```cpp
// tests/test_genome.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/genome.h"

using namespace botany;

TEST_CASE("default_genome returns valid parameters", "[genome]") {
    Genome g = default_genome();

    SECTION("hormone rates are positive") {
        REQUIRE(g.auxin_production_rate > 0.0f);
        REQUIRE(g.auxin_transport_rate > 0.0f);
        REQUIRE(g.auxin_decay_rate > 0.0f);
        REQUIRE(g.cytokinin_production_rate > 0.0f);
        REQUIRE(g.cytokinin_transport_rate > 0.0f);
        REQUIRE(g.cytokinin_decay_rate > 0.0f);
    }

    SECTION("transport rates are fractions (0, 1]") {
        REQUIRE(g.auxin_transport_rate <= 1.0f);
        REQUIRE(g.cytokinin_transport_rate <= 1.0f);
    }

    SECTION("decay rates are fractions (0, 1)") {
        REQUIRE(g.auxin_decay_rate < 1.0f);
        REQUIRE(g.cytokinin_decay_rate < 1.0f);
    }

    SECTION("growth rates are positive") {
        REQUIRE(g.growth_rate > 0.0f);
        REQUIRE(g.root_growth_rate > 0.0f);
    }

    SECTION("internode spacing is at least 1") {
        REQUIRE(g.internode_spacing >= 1);
        REQUIRE(g.root_internode_spacing >= 1);
    }
}
```

- [ ] **Step 4: Create stub source files so CMake can link**

Create minimal stubs for every `.cpp` file referenced in CMakeLists.txt so the project compiles. Each stub is an empty file or has a minimal implementation. The real implementations come in later tasks.

`src/engine/node.cpp`:
```cpp
#include "engine/node.h"
```

`src/engine/node.h`:
```cpp
#pragma once
namespace botany {}
```

`src/engine/plant.h`:
```cpp
#pragma once
namespace botany {}
```

`src/engine/plant.cpp`:
```cpp
#include "engine/plant.h"
```

`src/engine/hormone.h`:
```cpp
#pragma once
namespace botany {}
```

`src/engine/hormone.cpp`:
```cpp
#include "engine/hormone.h"
```

`src/engine/meristem.h`:
```cpp
#pragma once
namespace botany {}
```

`src/engine/meristem.cpp`:
```cpp
#include "engine/meristem.h"
```

`src/engine/engine.h`:
```cpp
#pragma once
namespace botany {}
```

`src/engine/engine.cpp`:
```cpp
#include "engine/engine.h"
```

`src/serialization/serializer.h`:
```cpp
#pragma once
namespace botany {}
```

`src/serialization/serializer.cpp`:
```cpp
#include "serialization/serializer.h"
```

`src/renderer/shader.h`:
```cpp
#pragma once
namespace botany {}
```

`src/renderer/shader.cpp`:
```cpp
#include "renderer/shader.h"
```

`src/renderer/camera.h`:
```cpp
#pragma once
namespace botany {}
```

`src/renderer/camera.cpp`:
```cpp
#include "renderer/camera.h"
```

`src/renderer/renderer.h`:
```cpp
#pragma once
namespace botany {}
```

`src/renderer/renderer.cpp`:
```cpp
#include "renderer/renderer.h"
```

`src/app_realtime.cpp`:
```cpp
int main() { return 0; }
```

`src/app_headless.cpp`:
```cpp
int main() { return 0; }
```

`src/app_playback.cpp`:
```cpp
int main() { return 0; }
```

Stub test files (empty Catch2 files):

`tests/test_node.cpp`:
```cpp
#include <catch2/catch_test_macros.hpp>
```

`tests/test_plant.cpp`:
```cpp
#include <catch2/catch_test_macros.hpp>
```

`tests/test_hormone.cpp`:
```cpp
#include <catch2/catch_test_macros.hpp>
```

`tests/test_meristem.cpp`:
```cpp
#include <catch2/catch_test_macros.hpp>
```

`tests/test_engine.cpp`:
```cpp
#include <catch2/catch_test_macros.hpp>
```

`tests/test_serializer.cpp`:
```cpp
#include <catch2/catch_test_macros.hpp>
```

- [ ] **Step 5: Build and run tests**

```bash
cd botany && mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
cmake --build . -j$(sysctl -n hw.ncpu)
ctest --output-on-failure
```

Expected: All genome tests pass. Build succeeds for all targets.

- [ ] **Step 6: Commit**

```bash
git init
echo -e "build/\n.DS_Store\ncompile_commands.json\n.superpowers/" > .gitignore
git add CMakeLists.txt src/ tests/ .gitignore docs/
git commit -m "feat: project scaffolding with CMake, deps, and Genome struct"
```

---

### Task 2: Node and Meristem Data Structures

**Files:**
- Modify: `src/engine/node.h`
- Modify: `src/engine/node.cpp`
- Modify: `tests/test_node.cpp`

- [ ] **Step 1: Write the failing tests**

```cpp
// tests/test_node.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/node.h"

using namespace botany;

TEST_CASE("Node creation with default values", "[node]") {
    Node node(1, NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);

    REQUIRE(node.id == 1);
    REQUIRE(node.type == NodeType::STEM);
    REQUIRE(node.position.y == 1.0f);
    REQUIRE(node.radius == 0.05f);
    REQUIRE(node.parent == nullptr);
    REQUIRE(node.children.empty());
    REQUIRE(node.age == 0);
    REQUIRE(node.auxin == 0.0f);
    REQUIRE(node.cytokinin == 0.0f);
    REQUIRE(node.meristem == nullptr);
    REQUIRE(node.leaf == nullptr);
}

TEST_CASE("add_child establishes parent-child relationship", "[node]") {
    Node parent(1, NodeType::STEM, glm::vec3(0.0f), 0.05f);
    Node child(2, NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.04f);

    parent.add_child(&child);

    REQUIRE(parent.children.size() == 1);
    REQUIRE(parent.children[0] == &child);
    REQUIRE(child.parent == &parent);
}

TEST_CASE("Node can have a meristem attached", "[node]") {
    Node node(1, NodeType::STEM, glm::vec3(0.0f), 0.05f);
    Meristem m{MeristemType::APICAL, true, 0};
    node.meristem = &m;

    REQUIRE(node.meristem != nullptr);
    REQUIRE(node.meristem->type == MeristemType::APICAL);
    REQUIRE(node.meristem->active == true);
}

TEST_CASE("Node can have a leaf attached", "[node]") {
    Node node(1, NodeType::STEM, glm::vec3(0.0f), 0.05f);
    Leaf l{0.3f};
    node.leaf = &l;

    REQUIRE(node.leaf != nullptr);
    REQUIRE(node.leaf->size == 0.3f);
}

TEST_CASE("Meristem types cover all four variants", "[node]") {
    REQUIRE(static_cast<int>(MeristemType::APICAL) != static_cast<int>(MeristemType::AXILLARY));
    REQUIRE(static_cast<int>(MeristemType::ROOT_APICAL) != static_cast<int>(MeristemType::ROOT_AXILLARY));
}
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu) 2>&1
```

Expected: Compilation fails — `Node`, `Meristem`, `Leaf`, `NodeType`, `MeristemType` not defined.

- [ ] **Step 3: Implement node.h**

```cpp
// src/engine/node.h
#pragma once

#include <cstdint>
#include <vector>
#include <glm/vec3.hpp>

namespace botany {

enum class NodeType { STEM, ROOT };

enum class MeristemType { APICAL, AXILLARY, ROOT_APICAL, ROOT_AXILLARY };

struct Meristem {
    MeristemType type;
    bool active;
    uint32_t ticks_since_last_node;
};

struct Leaf {
    float size;
};

struct Node {
    uint32_t id;
    Node* parent;
    std::vector<Node*> children;

    glm::vec3 position;
    float radius;

    NodeType type;
    uint32_t age;
    float auxin;
    float cytokinin;

    Meristem* meristem;
    Leaf* leaf;

    Node(uint32_t id, NodeType type, glm::vec3 position, float radius);

    void add_child(Node* child);
};

} // namespace botany
```

- [ ] **Step 4: Implement node.cpp**

```cpp
// src/engine/node.cpp
#include "engine/node.h"

namespace botany {

Node::Node(uint32_t id, NodeType type, glm::vec3 position, float radius)
    : id(id)
    , parent(nullptr)
    , position(position)
    , radius(radius)
    , type(type)
    , age(0)
    , auxin(0.0f)
    , cytokinin(0.0f)
    , meristem(nullptr)
    , leaf(nullptr)
{}

void Node::add_child(Node* child) {
    children.push_back(child);
    child->parent = this;
}

} // namespace botany
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu) && ctest --output-on-failure -R "node"
```

Expected: All node tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node.h src/engine/node.cpp tests/test_node.cpp
git commit -m "feat: Node, Meristem, Leaf structs with parent-child relationships"
```

---

### Task 3: Plant Class — Graph Ownership and Seed Initialization

**Files:**
- Modify: `src/engine/plant.h`
- Modify: `src/engine/plant.cpp`
- Modify: `tests/test_plant.cpp`

The Plant class owns all nodes and meristems (manages their memory) and provides the seed initialization logic.

- [ ] **Step 1: Write the failing tests**

```cpp
// tests/test_plant.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/plant.h"

using namespace botany;

TEST_CASE("Plant seed initialization creates correct graph", "[plant]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f, 0.0f, 0.0f));

    SECTION("has a seed node at origin") {
        const Node* seed = plant.seed();
        REQUIRE(seed != nullptr);
        REQUIRE(seed->position == glm::vec3(0.0f));
        REQUIRE(seed->parent == nullptr);
    }

    SECTION("seed has two children: shoot and root") {
        const Node* seed = plant.seed();
        REQUIRE(seed->children.size() == 2);
    }

    SECTION("one child is a STEM with apical meristem") {
        const Node* seed = plant.seed();
        bool found_shoot = false;
        for (const Node* child : seed->children) {
            if (child->type == NodeType::STEM && child->meristem &&
                child->meristem->type == MeristemType::APICAL) {
                found_shoot = true;
                REQUIRE(child->meristem->active == true);
                REQUIRE(child->position.y >= 0.0f);
            }
        }
        REQUIRE(found_shoot);
    }

    SECTION("one child is a ROOT with root apical meristem") {
        const Node* seed = plant.seed();
        bool found_root = false;
        for (const Node* child : seed->children) {
            if (child->type == NodeType::ROOT && child->meristem &&
                child->meristem->type == MeristemType::ROOT_APICAL) {
                found_root = true;
                REQUIRE(child->meristem->active == true);
                REQUIRE(child->position.y <= 0.0f);
            }
        }
        REQUIRE(found_root);
    }

    SECTION("initial radius from genome") {
        const Node* seed = plant.seed();
        REQUIRE(seed->radius == g.initial_radius);
    }

    SECTION("node_count starts at 3") {
        REQUIRE(plant.node_count() == 3);
    }
}

TEST_CASE("Plant provides access to all nodes", "[plant]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    std::vector<const Node*> all;
    plant.for_each_node([&](const Node& n) {
        all.push_back(&n);
    });

    REQUIRE(all.size() == 3);
}

TEST_CASE("Plant can create new nodes", "[plant]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* parent = plant.seed_mut();
    Node* child = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 2.0f, 0.0f), g.initial_radius);
    parent->add_child(child);

    REQUIRE(plant.node_count() == 4);
    REQUIRE(child->parent == parent);
}

TEST_CASE("Plant can create meristems and leaves", "[plant]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* node = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    Meristem* m = plant.create_meristem(MeristemType::AXILLARY, false);
    node->meristem = m;
    Leaf* l = plant.create_leaf(0.3f);
    node->leaf = l;

    REQUIRE(node->meristem->type == MeristemType::AXILLARY);
    REQUIRE(node->meristem->active == false);
    REQUIRE(node->leaf->size == 0.3f);
}
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu) 2>&1
```

Expected: Compilation fails — `Plant` class not defined.

- [ ] **Step 3: Implement plant.h**

```cpp
// src/engine/plant.h
#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <vector>
#include <glm/vec3.hpp>
#include "engine/genome.h"
#include "engine/node.h"

namespace botany {

class Plant {
public:
    Plant(const Genome& genome, glm::vec3 position);

    const Genome& genome() const { return genome_; }

    const Node* seed() const { return nodes_[0].get(); }
    Node* seed_mut() { return nodes_[0].get(); }

    uint32_t node_count() const { return static_cast<uint32_t>(nodes_.size()); }

    Node* create_node(NodeType type, glm::vec3 position, float radius);
    Meristem* create_meristem(MeristemType type, bool active);
    Leaf* create_leaf(float size);

    void for_each_node(std::function<void(const Node&)> fn) const;
    void for_each_node_mut(std::function<void(Node&)> fn);

    uint32_t next_id() { return next_id_++; }

private:
    Genome genome_;
    uint32_t next_id_ = 0;
    std::vector<std::unique_ptr<Node>> nodes_;
    std::vector<std::unique_ptr<Meristem>> meristems_;
    std::vector<std::unique_ptr<Leaf>> leaves_;
};

} // namespace botany
```

- [ ] **Step 4: Implement plant.cpp**

```cpp
// src/engine/plant.cpp
#include "engine/plant.h"

namespace botany {

Plant::Plant(const Genome& genome, glm::vec3 position)
    : genome_(genome)
{
    // Seed node
    Node* seed = create_node(NodeType::STEM, position, genome.initial_radius);

    // Shoot apical meristem node (just above seed)
    Node* shoot = create_node(NodeType::STEM, position + glm::vec3(0.0f, 0.01f, 0.0f), genome.initial_radius);
    Meristem* shoot_m = create_meristem(MeristemType::APICAL, true);
    shoot->meristem = shoot_m;
    seed->add_child(shoot);

    // Root apical meristem node (just below seed)
    Node* root = create_node(NodeType::ROOT, position - glm::vec3(0.0f, 0.01f, 0.0f), genome.initial_radius);
    Meristem* root_m = create_meristem(MeristemType::ROOT_APICAL, true);
    root->meristem = root_m;
    seed->add_child(root);
}

Node* Plant::create_node(NodeType type, glm::vec3 position, float radius) {
    auto node = std::make_unique<Node>(next_id(), type, position, radius);
    Node* ptr = node.get();
    nodes_.push_back(std::move(node));
    return ptr;
}

Meristem* Plant::create_meristem(MeristemType type, bool active) {
    auto m = std::make_unique<Meristem>(Meristem{type, active, 0});
    Meristem* ptr = m.get();
    meristems_.push_back(std::move(m));
    return ptr;
}

Leaf* Plant::create_leaf(float size) {
    auto l = std::make_unique<Leaf>(Leaf{size});
    Leaf* ptr = l.get();
    leaves_.push_back(std::move(l));
    return ptr;
}

void Plant::for_each_node(std::function<void(const Node&)> fn) const {
    for (const auto& node : nodes_) {
        fn(*node);
    }
}

void Plant::for_each_node_mut(std::function<void(Node&)> fn) {
    for (auto& node : nodes_) {
        fn(*node);
    }
}

} // namespace botany
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu) && ctest --output-on-failure -R "plant"
```

Expected: All plant tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/engine/plant.h src/engine/plant.cpp tests/test_plant.cpp
git commit -m "feat: Plant class with seed initialization and node ownership"
```

---

### Task 4: Hormone Transport

**Files:**
- Modify: `src/engine/hormone.h`
- Modify: `src/engine/hormone.cpp`
- Modify: `tests/test_hormone.cpp`

Two free functions that operate on a Plant: `transport_auxin()` (post-order, tips→roots) and `transport_cytokinin()` (pre-order, roots→tips).

- [ ] **Step 1: Write the failing tests**

```cpp
// tests/test_hormone.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/plant.h"
#include "engine/hormone.h"

using namespace botany;
using Catch::Matchers::WithinAbs;

// Helper: build a simple 3-node linear stem (seed -> mid -> tip_with_apical)
static Plant make_linear_stem() {
    Genome g = default_genome();
    g.auxin_production_rate = 1.0f;
    g.auxin_transport_rate = 0.5f;
    g.auxin_decay_rate = 0.1f;
    Plant plant(g, glm::vec3(0.0f));
    // Plant already has: seed(0) -> shoot(1), seed(0) -> root(2)
    // shoot(1) has an APICAL meristem
    return plant;
}

TEST_CASE("transport_auxin: apical meristem produces auxin at its node", "[hormone]") {
    Plant plant = make_linear_stem();
    const Genome& g = plant.genome();

    transport_auxin(plant);

    // The shoot node (child 0 of seed, which is STEM type) should have auxin
    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::STEM) { shoot = c; break; }
    }
    REQUIRE(shoot != nullptr);
    // After production and decay: auxin = production_rate * (1 - decay_rate)
    // After transport some went to parent too
    REQUIRE(shoot->auxin > 0.0f);
}

TEST_CASE("transport_auxin: auxin flows from child to parent", "[hormone]") {
    Plant plant = make_linear_stem();

    transport_auxin(plant);

    const Node* seed = plant.seed();
    // Seed should have received auxin from the shoot child
    REQUIRE(seed->auxin > 0.0f);
}

TEST_CASE("transport_auxin: non-meristem nodes don't produce auxin", "[hormone]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Root node has ROOT_APICAL, not APICAL — should not produce auxin
    transport_auxin(plant);

    const Node* seed = plant.seed();
    const Node* root = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::ROOT) { root = c; break; }
    }
    REQUIRE(root != nullptr);
    REQUIRE(root->auxin == 0.0f);
}

TEST_CASE("transport_cytokinin: root apical produces cytokinin", "[hormone]") {
    Plant plant = make_linear_stem();

    transport_cytokinin(plant);

    const Node* seed = plant.seed();
    const Node* root = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::ROOT) { root = c; break; }
    }
    REQUIRE(root != nullptr);
    REQUIRE(root->cytokinin > 0.0f);
}

TEST_CASE("transport_cytokinin: cytokinin flows from parent to children", "[hormone]") {
    Genome g = default_genome();
    g.cytokinin_production_rate = 2.0f;
    Plant plant(g, glm::vec3(0.0f));

    // Run multiple passes so cytokinin propagates from root through seed to shoot
    for (int i = 0; i < 5; i++) {
        transport_cytokinin(plant);
    }

    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::STEM) { shoot = c; break; }
    }
    REQUIRE(shoot != nullptr);
    // After several passes, cytokinin should have reached the shoot through the seed
    REQUIRE(shoot->cytokinin > 0.0f);
}

TEST_CASE("transport_auxin: auxin decays each pass", "[hormone]") {
    Plant plant = make_linear_stem();

    // Run one pass
    transport_auxin(plant);
    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::STEM) { shoot = c; break; }
    }

    float after_one = shoot->auxin;

    // Run another pass — auxin should be higher (production continues) but
    // the per-node value demonstrates decay is happening (not accumulating without bound)
    transport_auxin(plant);
    transport_auxin(plant);
    transport_auxin(plant);

    // With decay, concentration should stabilize rather than grow linearly
    // After 4 passes with production=1, transport=0.5, decay=0.1:
    // It should converge, not keep growing proportional to pass count
    float after_four = shoot->auxin;
    REQUIRE(after_four < 4.0f * after_one);  // Not linear growth
}
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu) 2>&1
```

Expected: Compilation fails — `transport_auxin`, `transport_cytokinin` not defined.

- [ ] **Step 3: Implement hormone.h**

```cpp
// src/engine/hormone.h
#pragma once

namespace botany {

class Plant;

void transport_auxin(Plant& plant);
void transport_cytokinin(Plant& plant);

} // namespace botany
```

- [ ] **Step 4: Implement hormone.cpp**

```cpp
// src/engine/hormone.cpp
#include "engine/hormone.h"
#include "engine/plant.h"
#include "engine/node.h"

namespace botany {

// Post-order traversal: process children before parent
static void auxin_postorder(Node& node, const Genome& genome) {
    for (Node* child : node.children) {
        auxin_postorder(*child, genome);
    }

    // Production: only shoot apical meristems produce auxin
    if (node.meristem && node.meristem->active &&
        node.meristem->type == MeristemType::APICAL) {
        node.auxin += genome.auxin_production_rate;
    }

    // Transport: send fraction to parent
    if (node.parent) {
        float flow = node.auxin * genome.auxin_transport_rate;
        node.parent->auxin += flow;
        node.auxin -= flow;
    }

    // Decay
    node.auxin *= (1.0f - genome.auxin_decay_rate);
}

void transport_auxin(Plant& plant) {
    auxin_postorder(*plant.seed_mut(), plant.genome());
}

// Pre-order traversal: process parent before children
static void cytokinin_preorder(Node& node, const Genome& genome) {
    // Production: only root apical meristems produce cytokinin
    if (node.meristem && node.meristem->active &&
        node.meristem->type == MeristemType::ROOT_APICAL) {
        node.cytokinin += genome.cytokinin_production_rate;
    }

    // Transport: send fraction to children, split evenly
    if (!node.children.empty()) {
        float flow_total = node.cytokinin * genome.cytokinin_transport_rate;
        float flow_per_child = flow_total / static_cast<float>(node.children.size());
        for (Node* child : node.children) {
            child->cytokinin += flow_per_child;
        }
        node.cytokinin -= flow_total;
    }

    // Decay
    node.cytokinin *= (1.0f - genome.cytokinin_decay_rate);

    for (Node* child : node.children) {
        cytokinin_preorder(*child, genome);
    }
}

void transport_cytokinin(Plant& plant) {
    cytokinin_preorder(*plant.seed_mut(), plant.genome());
}

} // namespace botany
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu) && ctest --output-on-failure -R "hormone"
```

Expected: All hormone tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/engine/hormone.h src/engine/hormone.cpp tests/test_hormone.cpp
git commit -m "feat: auxin and cytokinin transport through node graph"
```

---

### Task 5: Meristem Tick Behavior

**Files:**
- Modify: `src/engine/meristem.h`
- Modify: `src/engine/meristem.cpp`
- Modify: `tests/test_meristem.cpp`

A single function `tick_meristems(Plant&)` that iterates all nodes, finds active meristems, and executes their tick behavior: extend, thicken, spawn axillary nodes, chain growth, and axillary activation.

- [ ] **Step 1: Write the failing tests**

```cpp
// tests/test_meristem.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/plant.h"
#include "engine/meristem.h"

using namespace botany;
using Catch::Matchers::WithinAbs;

TEST_CASE("Shoot apical meristem extends node position upward", "[meristem]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::STEM) { shoot = c; break; }
    }
    float y_before = shoot->position.y;

    tick_meristems(plant);

    REQUIRE(shoot->position.y > y_before);
    REQUIRE_THAT(shoot->position.y - y_before, WithinAbs(g.growth_rate, 0.001f));
}

TEST_CASE("Shoot apical meristem thickens its node", "[meristem]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::STEM) { shoot = c; break; }
    }
    float r_before = shoot->radius;

    tick_meristems(plant);

    REQUIRE(shoot->radius > r_before);
}

TEST_CASE("Shoot apical spawns axillary node at internode_spacing", "[meristem]") {
    Genome g = default_genome();
    g.internode_spacing = 3;
    Plant plant(g, glm::vec3(0.0f));

    uint32_t initial_count = plant.node_count();

    // Tick 3 times to hit internode_spacing
    for (int i = 0; i < 3; i++) {
        tick_meristems(plant);
    }

    // Should have spawned an axillary node with a dormant meristem and a leaf
    REQUIRE(plant.node_count() > initial_count);

    // Find the new axillary node
    bool found_axillary = false;
    plant.for_each_node([&](const Node& n) {
        if (n.meristem && n.meristem->type == MeristemType::AXILLARY) {
            found_axillary = true;
            REQUIRE(n.meristem->active == false);
            REQUIRE(n.leaf != nullptr);
            REQUIRE(n.leaf->size == g.leaf_size);
        }
    });
    REQUIRE(found_axillary);
}

TEST_CASE("Root apical meristem extends node position downward", "[meristem]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    const Node* seed = plant.seed();
    const Node* root = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::ROOT) { root = c; break; }
    }
    float y_before = root->position.y;

    tick_meristems(plant);

    REQUIRE(root->position.y < y_before);
}

TEST_CASE("Root apical spawns root axillary at root_internode_spacing", "[meristem]") {
    Genome g = default_genome();
    g.root_internode_spacing = 2;
    Plant plant(g, glm::vec3(0.0f));

    for (int i = 0; i < 2; i++) {
        tick_meristems(plant);
    }

    bool found_root_axillary = false;
    plant.for_each_node([&](const Node& n) {
        if (n.meristem && n.meristem->type == MeristemType::ROOT_AXILLARY) {
            found_root_axillary = true;
            REQUIRE(n.meristem->active == false);
        }
    });
    REQUIRE(found_root_axillary);
}

TEST_CASE("Chain growth: apical meristem transfers to new node when internode too long", "[meristem]") {
    Genome g = default_genome();
    g.growth_rate = 0.5f;
    g.max_internode_length = 0.6f;
    g.internode_spacing = 100; // high so axillary doesn't interfere
    Plant plant(g, glm::vec3(0.0f));

    // After 2 ticks the shoot moves 1.0 from its start, exceeding max_internode_length of 0.6
    tick_meristems(plant);
    tick_meristems(plant);

    // The original shoot node should no longer have a meristem (became interior)
    // A new node should have the apical meristem
    int apical_count = 0;
    plant.for_each_node([&](const Node& n) {
        if (n.meristem && n.meristem->type == MeristemType::APICAL && n.meristem->active) {
            apical_count++;
        }
    });
    REQUIRE(apical_count == 1); // Still exactly one apical

    // Node count should have increased (new chain node created)
    REQUIRE(plant.node_count() > 3);
}

TEST_CASE("Axillary meristem activates when auxin low and cytokinin high", "[meristem]") {
    Genome g = default_genome();
    g.auxin_threshold = 1.0f;
    g.cytokinin_threshold = 0.0f; // cytokinin just needs to be > 0
    g.internode_spacing = 1;
    Plant plant(g, glm::vec3(0.0f));

    // Tick once to spawn an axillary node
    tick_meristems(plant);

    // Manually set hormones on the axillary node to trigger activation
    Node* axillary_node = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.meristem && n.meristem->type == MeristemType::AXILLARY) {
            axillary_node = &n;
        }
    });
    REQUIRE(axillary_node != nullptr);
    axillary_node->auxin = 0.1f;      // below threshold of 1.0
    axillary_node->cytokinin = 0.5f;   // above threshold of 0.0

    tick_meristems(plant);

    // Should have converted to APICAL and be active
    REQUIRE(axillary_node->meristem->type == MeristemType::APICAL);
    REQUIRE(axillary_node->meristem->active == true);
}

TEST_CASE("Axillary meristem stays dormant when auxin is high", "[meristem]") {
    Genome g = default_genome();
    g.auxin_threshold = 1.0f;
    g.internode_spacing = 1;
    Plant plant(g, glm::vec3(0.0f));

    tick_meristems(plant);

    Node* axillary_node = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.meristem && n.meristem->type == MeristemType::AXILLARY) {
            axillary_node = &n;
        }
    });
    REQUIRE(axillary_node != nullptr);
    axillary_node->auxin = 5.0f; // way above threshold

    tick_meristems(plant);

    REQUIRE(axillary_node->meristem->type == MeristemType::AXILLARY);
    REQUIRE(axillary_node->meristem->active == false);
}

TEST_CASE("Node age increments each tick", "[meristem]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    tick_meristems(plant);
    tick_meristems(plant);

    const Node* seed = plant.seed();
    REQUIRE(seed->age == 2);
}
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu) 2>&1
```

Expected: Compilation fails — `tick_meristems` not defined.

- [ ] **Step 3: Implement meristem.h**

```cpp
// src/engine/meristem.h
#pragma once

namespace botany {

class Plant;

void tick_meristems(Plant& plant);

} // namespace botany
```

- [ ] **Step 4: Implement meristem.cpp**

```cpp
// src/engine/meristem.cpp
#include "engine/meristem.h"
#include "engine/plant.h"
#include "engine/node.h"
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>
#include <cmath>
#include <vector>

namespace botany {

// Compute growth direction for a node: normalized vector from parent to this node
// If no parent (seed), defaults to up for STEM, down for ROOT
static glm::vec3 growth_direction(const Node& node) {
    if (node.parent) {
        glm::vec3 dir = node.position - node.parent->position;
        float len = glm::length(dir);
        if (len > 0.0001f) {
            return dir / len;
        }
    }
    return (node.type == NodeType::STEM) ? glm::vec3(0.0f, 1.0f, 0.0f)
                                          : glm::vec3(0.0f, -1.0f, 0.0f);
}

// Compute a branch direction offset from the main growth direction
static glm::vec3 branch_direction(const glm::vec3& main_dir, float angle, uint32_t seed) {
    // Find a perpendicular vector
    glm::vec3 perp;
    if (std::abs(main_dir.y) < 0.9f) {
        perp = glm::normalize(glm::cross(main_dir, glm::vec3(0.0f, 1.0f, 0.0f)));
    } else {
        perp = glm::normalize(glm::cross(main_dir, glm::vec3(1.0f, 0.0f, 0.0f)));
    }

    // Rotate perp around main_dir by a pseudo-random angle based on seed
    float rotate = static_cast<float>(seed) * 2.399f; // golden angle in radians
    float cos_r = std::cos(rotate);
    float sin_r = std::sin(rotate);
    glm::vec3 perp2 = glm::normalize(glm::cross(main_dir, perp));
    glm::vec3 rotated_perp = perp * cos_r + perp2 * sin_r;

    // Tilt from main_dir toward rotated_perp by branch_angle
    float cos_a = std::cos(angle);
    float sin_a = std::sin(angle);
    return glm::normalize(main_dir * cos_a + rotated_perp * sin_a);
}

static void tick_shoot_apical(Node& node, Plant& plant) {
    const Genome& g = plant.genome();
    Meristem* m = node.meristem;

    // Extend
    glm::vec3 dir = growth_direction(node);
    node.position += dir * g.growth_rate;

    // Thicken
    node.radius += g.thickening_rate;

    // Node spacing check
    m->ticks_since_last_node++;
    if (m->ticks_since_last_node >= g.internode_spacing) {
        m->ticks_since_last_node = 0;

        // Spawn axillary node branching off
        glm::vec3 branch_dir = branch_direction(dir, g.branch_angle, node.id);
        glm::vec3 ax_pos = node.position + branch_dir * 0.01f;
        Node* axillary = plant.create_node(NodeType::STEM, ax_pos, g.initial_radius * 0.5f);
        Meristem* ax_m = plant.create_meristem(MeristemType::AXILLARY, false);
        axillary->meristem = ax_m;
        Leaf* leaf = plant.create_leaf(g.leaf_size);
        axillary->leaf = leaf;
        node.add_child(axillary);
    }

    // Chain growth: if distance to parent exceeds max, insert new node
    if (node.parent) {
        float dist = glm::length(node.position - node.parent->position);
        if (dist > g.max_internode_length) {
            // Create a new node at current position with the apical meristem
            Node* new_tip = plant.create_node(NodeType::STEM, node.position, node.radius);
            new_tip->meristem = node.meristem;
            node.meristem = nullptr;

            // Reparent: new_tip becomes child of node, inherits node's children
            // that were just added (axillary nodes stay on node)
            node.add_child(new_tip);

            // Reset the meristem counter on the transferred meristem
            new_tip->meristem->ticks_since_last_node = 0;
        }
    }
}

static void tick_root_apical(Node& node, Plant& plant) {
    const Genome& g = plant.genome();
    Meristem* m = node.meristem;

    // Extend downward
    glm::vec3 dir = growth_direction(node);
    node.position += dir * g.root_growth_rate;

    // Node spacing check — spawn root axillary
    m->ticks_since_last_node++;
    if (m->ticks_since_last_node >= g.root_internode_spacing) {
        m->ticks_since_last_node = 0;

        glm::vec3 branch_dir = branch_direction(dir, g.root_branch_angle, node.id);
        glm::vec3 ax_pos = node.position + branch_dir * 0.01f;
        Node* axillary = plant.create_node(NodeType::ROOT, ax_pos, g.initial_radius * 0.3f);
        Meristem* ax_m = plant.create_meristem(MeristemType::ROOT_AXILLARY, false);
        axillary->meristem = ax_m;
        node.add_child(axillary);
    }

    // Chain growth
    if (node.parent) {
        float dist = glm::length(node.position - node.parent->position);
        if (dist > g.root_max_internode_length) {
            Node* new_tip = plant.create_node(NodeType::ROOT, node.position, node.radius);
            new_tip->meristem = node.meristem;
            node.meristem = nullptr;
            node.add_child(new_tip);
            new_tip->meristem->ticks_since_last_node = 0;
        }
    }
}

static void tick_shoot_axillary(Node& node, const Genome& g) {
    if (node.auxin < g.auxin_threshold && node.cytokinin > g.cytokinin_threshold) {
        node.meristem->type = MeristemType::APICAL;
        node.meristem->active = true;
        node.meristem->ticks_since_last_node = 0;
    }
}

static void tick_root_axillary(Node& node, const Genome& g) {
    if (node.cytokinin < g.cytokinin_threshold) {
        node.meristem->type = MeristemType::ROOT_APICAL;
        node.meristem->active = true;
        node.meristem->ticks_since_last_node = 0;
    }
}

void tick_meristems(Plant& plant) {
    const Genome& g = plant.genome();

    // Collect nodes to tick first, since ticking may add new nodes
    std::vector<Node*> to_tick;
    plant.for_each_node_mut([&](Node& n) {
        n.age++;
        if (n.meristem) {
            to_tick.push_back(&n);
        }
    });

    for (Node* node : to_tick) {
        Meristem* m = node->meristem;
        if (!m) continue; // may have been transferred by chain growth

        switch (m->type) {
            case MeristemType::APICAL:
                if (m->active) tick_shoot_apical(*node, plant);
                break;
            case MeristemType::AXILLARY:
                if (!m->active) tick_shoot_axillary(*node, g);
                break;
            case MeristemType::ROOT_APICAL:
                if (m->active) tick_root_apical(*node, plant);
                break;
            case MeristemType::ROOT_AXILLARY:
                if (!m->active) tick_root_axillary(*node, g);
                break;
        }
    }
}

} // namespace botany
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu) && ctest --output-on-failure -R "meristem"
```

Expected: All meristem tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/engine/meristem.h src/engine/meristem.cpp tests/test_meristem.cpp
git commit -m "feat: meristem tick behavior — apical, axillary, root, chain growth"
```

---

### Task 6: Engine Class — Tick Loop and Plant Management

**Files:**
- Modify: `src/engine/engine.h`
- Modify: `src/engine/engine.cpp`
- Modify: `tests/test_engine.cpp`

The Engine owns Plants, runs the full tick cycle (hormones + meristems), and exposes read-only access for the renderer.

- [ ] **Step 1: Write the failing tests**

```cpp
// tests/test_engine.cpp
#include <catch2/catch_test_macros.hpp>
#include "engine/engine.h"

using namespace botany;

TEST_CASE("Engine starts at tick 0", "[engine]") {
    Engine engine;
    REQUIRE(engine.get_tick() == 0);
}

TEST_CASE("Engine creates a plant and returns valid ID", "[engine]") {
    Engine engine;
    Genome g = default_genome();
    PlantID id = engine.create_plant(g, glm::vec3(0.0f));

    const Plant& plant = engine.get_plant(id);
    REQUIRE(plant.node_count() == 3); // seed + shoot + root
}

TEST_CASE("Engine tick advances tick counter", "[engine]") {
    Engine engine;
    engine.tick();
    REQUIRE(engine.get_tick() == 1);
    engine.tick();
    REQUIRE(engine.get_tick() == 2);
}

TEST_CASE("Engine tick grows the plant", "[engine]") {
    Engine engine;
    Genome g = default_genome();
    PlantID id = engine.create_plant(g, glm::vec3(0.0f));

    const Node* shoot_before = nullptr;
    engine.get_plant(id).for_each_node([&](const Node& n) {
        if (n.meristem && n.meristem->type == MeristemType::APICAL) {
            shoot_before = &n;
        }
    });
    float y_before = shoot_before->position.y;

    engine.tick();

    // Shoot should have moved up
    float y_after = shoot_before->position.y;
    REQUIRE(y_after > y_before);
}

TEST_CASE("Engine runs multiple ticks and plant grows complex structure", "[engine]") {
    Engine engine;
    Genome g = default_genome();
    g.internode_spacing = 3;
    PlantID id = engine.create_plant(g, glm::vec3(0.0f));

    for (int i = 0; i < 20; i++) {
        engine.tick();
    }

    // After 20 ticks with internode_spacing=3, should have several axillary nodes
    const Plant& plant = engine.get_plant(id);
    REQUIRE(plant.node_count() > 3);

    // Verify hormones are non-zero somewhere
    bool found_auxin = false;
    plant.for_each_node([&](const Node& n) {
        if (n.auxin > 0.0f) found_auxin = true;
    });
    REQUIRE(found_auxin);
}

TEST_CASE("Engine supports multiple plants", "[engine]") {
    Engine engine;
    Genome g = default_genome();
    PlantID id1 = engine.create_plant(g, glm::vec3(-5.0f, 0.0f, 0.0f));
    PlantID id2 = engine.create_plant(g, glm::vec3(5.0f, 0.0f, 0.0f));

    REQUIRE(id1 != id2);

    engine.tick();

    REQUIRE(engine.get_plant(id1).seed()->position.x < 0.0f);
    REQUIRE(engine.get_plant(id2).seed()->position.x > 0.0f);
}
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu) 2>&1
```

Expected: Compilation fails — `Engine`, `PlantID` not defined.

- [ ] **Step 3: Implement engine.h**

```cpp
// src/engine/engine.h
#pragma once

#include <cstdint>
#include <memory>
#include <vector>
#include <glm/vec3.hpp>
#include "engine/genome.h"
#include "engine/plant.h"

namespace botany {

using PlantID = uint32_t;

class Engine {
public:
    Engine() = default;

    PlantID create_plant(const Genome& genome, glm::vec3 position);
    void tick();

    const Plant& get_plant(PlantID id) const;
    Plant& get_plant_mut(PlantID id);

    uint32_t get_tick() const { return tick_; }
    uint32_t plant_count() const { return static_cast<uint32_t>(plants_.size()); }

private:
    uint32_t tick_ = 0;
    std::vector<std::unique_ptr<Plant>> plants_;
};

} // namespace botany
```

- [ ] **Step 4: Implement engine.cpp**

```cpp
// src/engine/engine.cpp
#include "engine/engine.h"
#include "engine/hormone.h"
#include "engine/meristem.h"

namespace botany {

PlantID Engine::create_plant(const Genome& genome, glm::vec3 position) {
    PlantID id = static_cast<PlantID>(plants_.size());
    plants_.push_back(std::make_unique<Plant>(genome, position));
    return id;
}

void Engine::tick() {
    for (auto& plant : plants_) {
        transport_auxin(*plant);
        transport_cytokinin(*plant);
        tick_meristems(*plant);
    }
    tick_++;
}

const Plant& Engine::get_plant(PlantID id) const {
    return *plants_.at(id);
}

Plant& Engine::get_plant_mut(PlantID id) {
    return *plants_.at(id);
}

} // namespace botany
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu) && ctest --output-on-failure -R "engine"
```

Expected: All engine tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/engine/engine.h src/engine/engine.cpp tests/test_engine.cpp
git commit -m "feat: Engine class — tick loop with hormone transport and meristem growth"
```

---

### Task 7: Binary Serialization

**Files:**
- Modify: `src/serialization/serializer.h`
- Modify: `src/serialization/serializer.cpp`
- Modify: `tests/test_serializer.cpp`

Saves and loads plant state per-tick to a binary file for headless precompute and playback.

- [ ] **Step 1: Write the failing tests**

```cpp
// tests/test_serializer.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <sstream>
#include "engine/engine.h"
#include "serialization/serializer.h"

using namespace botany;

TEST_CASE("Round-trip: save and load single tick", "[serializer]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));
    engine.tick();

    // Save
    std::stringstream ss;
    Recording rec;
    rec.genome = g;
    save_tick(ss, engine, 0);

    // Load
    ss.seekg(0);
    TickSnapshot snap = load_tick(ss);

    REQUIRE(snap.tick_number == 1);
    REQUIRE(snap.nodes.size() == engine.get_plant(0).node_count());
}

TEST_CASE("Round-trip: node positions preserved", "[serializer]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(5.0f, 0.0f, -3.0f));
    engine.tick();
    engine.tick();

    std::stringstream ss;
    save_tick(ss, engine, 0);

    ss.seekg(0);
    TickSnapshot snap = load_tick(ss);

    // Find a node and verify position matches
    const Plant& plant = engine.get_plant(0);
    std::vector<const Node*> engine_nodes;
    plant.for_each_node([&](const Node& n) {
        engine_nodes.push_back(&n);
    });

    REQUIRE(snap.nodes.size() == engine_nodes.size());
    for (size_t i = 0; i < snap.nodes.size(); i++) {
        REQUIRE(snap.nodes[i].id == engine_nodes[i]->id);
        REQUIRE(snap.nodes[i].position.x == engine_nodes[i]->position.x);
        REQUIRE(snap.nodes[i].position.y == engine_nodes[i]->position.y);
        REQUIRE(snap.nodes[i].position.z == engine_nodes[i]->position.z);
        REQUIRE(snap.nodes[i].radius == engine_nodes[i]->radius);
    }
}

TEST_CASE("Round-trip: parent_id preserved for graph reconstruction", "[serializer]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));
    for (int i = 0; i < 5; i++) engine.tick();

    std::stringstream ss;
    save_tick(ss, engine, 0);

    ss.seekg(0);
    TickSnapshot snap = load_tick(ss);

    // Seed node should have parent_id of UINT32_MAX (no parent)
    REQUIRE(snap.nodes[0].parent_id == UINT32_MAX);
    // Other nodes should have valid parent ids
    for (size_t i = 1; i < snap.nodes.size(); i++) {
        REQUIRE(snap.nodes[i].parent_id != UINT32_MAX);
    }
}

TEST_CASE("Save recording header and multiple ticks", "[serializer]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));

    std::stringstream ss;
    save_recording_header(ss, g, 10);

    for (int i = 0; i < 10; i++) {
        engine.tick();
        save_tick(ss, engine, 0);
    }

    // Read it back
    ss.seekg(0);
    RecordingHeader header = load_recording_header(ss);
    REQUIRE(header.num_ticks == 10);

    // Load first tick
    TickSnapshot snap1 = load_tick(ss);
    REQUIRE(snap1.tick_number == 1);
    REQUIRE(snap1.nodes.size() >= 3);
}
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu) 2>&1
```

Expected: Compilation fails — serializer types not defined.

- [ ] **Step 3: Implement serializer.h**

```cpp
// src/serialization/serializer.h
#pragma once

#include <cstdint>
#include <iostream>
#include <vector>
#include <glm/vec3.hpp>
#include "engine/genome.h"
#include "engine/node.h"

namespace botany {

class Engine;

struct NodeSnapshot {
    uint32_t id;
    uint32_t parent_id;  // UINT32_MAX if no parent
    NodeType type;
    glm::vec3 position;
    float radius;
    float auxin;
    float cytokinin;
    bool has_leaf;
    float leaf_size;
};

struct TickSnapshot {
    uint32_t tick_number;
    std::vector<NodeSnapshot> nodes;
};

struct RecordingHeader {
    uint32_t num_ticks;
    Genome genome;
};

struct Recording {
    Genome genome;
};

void save_recording_header(std::ostream& out, const Genome& genome, uint32_t num_ticks);
RecordingHeader load_recording_header(std::istream& in);

void save_tick(std::ostream& out, const Engine& engine, uint32_t plant_id);
TickSnapshot load_tick(std::istream& in);

} // namespace botany
```

- [ ] **Step 4: Implement serializer.cpp**

```cpp
// src/serialization/serializer.cpp
#include "serialization/serializer.h"
#include "engine/engine.h"

namespace botany {

template<typename T>
static void write_val(std::ostream& out, const T& val) {
    out.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

template<typename T>
static T read_val(std::istream& in) {
    T val;
    in.read(reinterpret_cast<char*>(&val), sizeof(T));
    return val;
}

void save_recording_header(std::ostream& out, const Genome& genome, uint32_t num_ticks) {
    write_val(out, num_ticks);
    write_val(out, genome);
}

RecordingHeader load_recording_header(std::istream& in) {
    RecordingHeader header;
    header.num_ticks = read_val<uint32_t>(in);
    header.genome = read_val<Genome>(in);
    return header;
}

void save_tick(std::ostream& out, const Engine& engine, uint32_t plant_id) {
    const Plant& plant = engine.get_plant(plant_id);
    uint32_t tick = engine.get_tick();
    uint32_t count = plant.node_count();

    write_val(out, tick);
    write_val(out, count);

    plant.for_each_node([&](const Node& node) {
        write_val(out, node.id);
        uint32_t parent_id = node.parent ? node.parent->id : UINT32_MAX;
        write_val(out, parent_id);
        write_val(out, node.type);
        write_val(out, node.position);
        write_val(out, node.radius);
        write_val(out, node.auxin);
        write_val(out, node.cytokinin);
        bool has_leaf = (node.leaf != nullptr);
        write_val(out, has_leaf);
        float leaf_size = has_leaf ? node.leaf->size : 0.0f;
        write_val(out, leaf_size);
    });
}

TickSnapshot load_tick(std::istream& in) {
    TickSnapshot snap;
    snap.tick_number = read_val<uint32_t>(in);
    uint32_t count = read_val<uint32_t>(in);
    snap.nodes.resize(count);

    for (uint32_t i = 0; i < count; i++) {
        NodeSnapshot& ns = snap.nodes[i];
        ns.id = read_val<uint32_t>(in);
        ns.parent_id = read_val<uint32_t>(in);
        ns.type = read_val<NodeType>(in);
        ns.position = read_val<glm::vec3>(in);
        ns.radius = read_val<float>(in);
        ns.auxin = read_val<float>(in);
        ns.cytokinin = read_val<float>(in);
        ns.has_leaf = read_val<bool>(in);
        ns.leaf_size = read_val<float>(in);
    }

    return snap;
}

} // namespace botany
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu) && ctest --output-on-failure -R "serializer"
```

Expected: All serializer tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/serialization/serializer.h src/serialization/serializer.cpp tests/test_serializer.cpp
git commit -m "feat: binary serialization for tick snapshots and recordings"
```

---

### Task 8: OpenGL Renderer — Window, Shaders, Basic Pipeline

**Files:**
- Modify: `src/renderer/shader.h`
- Modify: `src/renderer/shader.cpp`
- Modify: `src/renderer/camera.h`
- Modify: `src/renderer/camera.cpp`
- Modify: `src/renderer/renderer.h`
- Modify: `src/renderer/renderer.cpp`
- Create: `shaders/plant.vert`
- Create: `shaders/plant.frag`

No automated tests for renderer (requires OpenGL context). Verified by visual inspection in Task 10.

- [ ] **Step 1: Create vertex shader**

```glsl
// shaders/plant.vert
#version 410 core

layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec3 aColor;

uniform mat4 uModel;
uniform mat4 uView;
uniform mat4 uProjection;

out vec3 vNormal;
out vec3 vFragPos;
out vec3 vColor;

void main() {
    vec4 worldPos = uModel * vec4(aPos, 1.0);
    vFragPos = worldPos.xyz;
    vNormal = mat3(transpose(inverse(uModel))) * aNormal;
    vColor = aColor;
    gl_Position = uProjection * uView * worldPos;
}
```

- [ ] **Step 2: Create fragment shader**

```glsl
// shaders/plant.frag
#version 410 core

in vec3 vNormal;
in vec3 vFragPos;
in vec3 vColor;

out vec4 FragColor;

uniform vec3 uLightDir;
uniform vec3 uAmbient;

void main() {
    vec3 norm = normalize(vNormal);
    vec3 lightDir = normalize(-uLightDir);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 result = (uAmbient + diff * vec3(1.0)) * vColor;
    FragColor = vec4(result, 1.0);
}
```

- [ ] **Step 3: Implement shader.h and shader.cpp**

```cpp
// src/renderer/shader.h
#pragma once

#include <string>
#include <cstdint>
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>

namespace botany {

class Shader {
public:
    Shader() = default;
    bool load(const std::string& vert_path, const std::string& frag_path);
    void use() const;
    void set_mat4(const std::string& name, const glm::mat4& mat) const;
    void set_vec3(const std::string& name, const glm::vec3& vec) const;
    uint32_t id() const { return program_; }

private:
    uint32_t program_ = 0;
};

} // namespace botany
```

```cpp
// src/renderer/shader.cpp
#include "renderer/shader.h"
#include <glad/gl.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <glm/gtc/type_ptr.hpp>

namespace botany {

static std::string read_file(const std::string& path) {
    std::ifstream file(path);
    std::stringstream ss;
    ss << file.rdbuf();
    return ss.str();
}

static uint32_t compile_shader(const std::string& source, GLenum type) {
    uint32_t shader = glCreateShader(type);
    const char* src = source.c_str();
    glShaderSource(shader, 1, &src, nullptr);
    glCompileShader(shader);

    int success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char log[512];
        glGetShaderInfoLog(shader, 512, nullptr, log);
        std::cerr << "Shader compilation error: " << log << std::endl;
    }
    return shader;
}

bool Shader::load(const std::string& vert_path, const std::string& frag_path) {
    std::string vert_src = read_file(vert_path);
    std::string frag_src = read_file(frag_path);

    uint32_t vert = compile_shader(vert_src, GL_VERTEX_SHADER);
    uint32_t frag = compile_shader(frag_src, GL_FRAGMENT_SHADER);

    program_ = glCreateProgram();
    glAttachShader(program_, vert);
    glAttachShader(program_, frag);
    glLinkProgram(program_);

    int success;
    glGetProgramiv(program_, GL_LINK_STATUS, &success);
    if (!success) {
        char log[512];
        glGetProgramInfoLog(program_, 512, nullptr, log);
        std::cerr << "Shader link error: " << log << std::endl;
        return false;
    }

    glDeleteShader(vert);
    glDeleteShader(frag);
    return true;
}

void Shader::use() const {
    glUseProgram(program_);
}

void Shader::set_mat4(const std::string& name, const glm::mat4& mat) const {
    glUniformMatrix4fv(glGetUniformLocation(program_, name.c_str()), 1, GL_FALSE, glm::value_ptr(mat));
}

void Shader::set_vec3(const std::string& name, const glm::vec3& vec) const {
    glUniform3fv(glGetUniformLocation(program_, name.c_str()), 1, glm::value_ptr(vec));
}

} // namespace botany
```

- [ ] **Step 4: Implement camera.h and camera.cpp**

```cpp
// src/renderer/camera.h
#pragma once

#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>

namespace botany {

class OrbitCamera {
public:
    OrbitCamera();

    glm::mat4 view_matrix() const;
    glm::mat4 projection_matrix(float aspect) const;

    void rotate(float dx, float dy);
    void zoom(float delta);
    void set_target(glm::vec3 target) { target_ = target; }

    glm::vec3 target() const { return target_; }

private:
    glm::vec3 target_;
    float distance_;
    float yaw_;    // radians
    float pitch_;  // radians
    float fov_;
    float near_;
    float far_;
};

} // namespace botany
```

```cpp
// src/renderer/camera.cpp
#include "renderer/camera.h"
#include <glm/gtc/matrix_transform.hpp>
#include <cmath>
#include <algorithm>

namespace botany {

OrbitCamera::OrbitCamera()
    : target_(0.0f, 0.0f, 0.0f)
    , distance_(5.0f)
    , yaw_(0.0f)
    , pitch_(0.3f)
    , fov_(45.0f)
    , near_(0.1f)
    , far_(100.0f)
{}

glm::mat4 OrbitCamera::view_matrix() const {
    float x = distance_ * std::cos(pitch_) * std::sin(yaw_);
    float y = distance_ * std::sin(pitch_);
    float z = distance_ * std::cos(pitch_) * std::cos(yaw_);
    glm::vec3 eye = target_ + glm::vec3(x, y, z);
    return glm::lookAt(eye, target_, glm::vec3(0.0f, 1.0f, 0.0f));
}

glm::mat4 OrbitCamera::projection_matrix(float aspect) const {
    return glm::perspective(glm::radians(fov_), aspect, near_, far_);
}

void OrbitCamera::rotate(float dx, float dy) {
    yaw_ += dx * 0.005f;
    pitch_ += dy * 0.005f;
    pitch_ = std::clamp(pitch_, -1.5f, 1.5f);
}

void OrbitCamera::zoom(float delta) {
    distance_ -= delta * 0.5f;
    distance_ = std::clamp(distance_, 0.5f, 50.0f);
}

} // namespace botany
```

- [ ] **Step 5: Implement renderer.h and renderer.cpp**

```cpp
// src/renderer/renderer.h
#pragma once

#include <string>
#include <cstdint>
#include "renderer/shader.h"
#include "renderer/camera.h"

struct GLFWwindow;

namespace botany {

class Plant;
struct TickSnapshot;

class Renderer {
public:
    bool init(int width, int height, const std::string& shader_dir);
    void shutdown();

    GLFWwindow* window() const { return window_; }
    OrbitCamera& camera() { return camera_; }

    void begin_frame();
    void draw_plant(const Plant& plant);
    void draw_snapshot(const TickSnapshot& snapshot);
    void draw_ground();
    void end_frame();

private:
    GLFWwindow* window_ = nullptr;
    Shader shader_;
    OrbitCamera camera_;
    int width_;
    int height_;

    uint32_t ground_vao_ = 0;
    uint32_t ground_vbo_ = 0;

    void setup_ground();
    void draw_cylinder(glm::vec3 start, glm::vec3 end,
                       float r_start, float r_end,
                       glm::vec3 color, int segments = 8);
    void draw_leaf(glm::vec3 position, glm::vec3 direction, float size);
};

} // namespace botany
```

```cpp
// src/renderer/renderer.cpp
#include "renderer/renderer.h"
#include <glad/gl.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/geometric.hpp>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <iostream>
#include "engine/plant.h"
#include "engine/node.h"
#include "serialization/serializer.h"

namespace botany {

bool Renderer::init(int width, int height, const std::string& shader_dir) {
    width_ = width;
    height_ = height;

    if (!glfwInit()) {
        std::cerr << "Failed to init GLFW" << std::endl;
        return false;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    window_ = glfwCreateWindow(width, height, "Botany", nullptr, nullptr);
    if (!window_) {
        std::cerr << "Failed to create window" << std::endl;
        glfwTerminate();
        return false;
    }
    glfwMakeContextCurrent(window_);

    int version = gladLoadGL(glfwGetProcAddress);
    if (!version) {
        std::cerr << "Failed to load OpenGL" << std::endl;
        return false;
    }

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.53f, 0.81f, 0.92f, 1.0f); // sky blue

    if (!shader_.load(shader_dir + "/plant.vert", shader_dir + "/plant.frag")) {
        return false;
    }

    setup_ground();
    return true;
}

void Renderer::shutdown() {
    if (ground_vao_) { glDeleteVertexArrays(1, &ground_vao_); ground_vao_ = 0; }
    if (ground_vbo_) { glDeleteBuffers(1, &ground_vbo_); ground_vbo_ = 0; }
    if (window_) { glfwDestroyWindow(window_); window_ = nullptr; }
    glfwTerminate();
}

void Renderer::begin_frame() {
    glfwGetFramebufferSize(window_, &width_, &height_);
    glViewport(0, 0, width_, height_);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    shader_.use();
    float aspect = static_cast<float>(width_) / static_cast<float>(height_);
    shader_.set_mat4("uView", camera_.view_matrix());
    shader_.set_mat4("uProjection", camera_.projection_matrix(aspect));
    shader_.set_mat4("uModel", glm::mat4(1.0f));
    shader_.set_vec3("uLightDir", glm::vec3(-0.3f, -1.0f, -0.5f));
    shader_.set_vec3("uAmbient", glm::vec3(0.3f, 0.3f, 0.3f));
}

void Renderer::draw_plant(const Plant& plant) {
    plant.for_each_node([&](const Node& node) {
        if (!node.parent) return; // seed has no segment to draw

        glm::vec3 color;
        if (node.type == NodeType::STEM) {
            color = glm::vec3(0.45f, 0.3f, 0.15f); // brown
        } else {
            color = glm::vec3(0.35f, 0.2f, 0.1f); // dark brown for roots
        }

        draw_cylinder(node.parent->position, node.position,
                      node.parent->radius, node.radius, color);

        if (node.leaf) {
            glm::vec3 dir = (node.parent)
                ? glm::normalize(node.position - node.parent->position)
                : glm::vec3(0.0f, 1.0f, 0.0f);
            draw_leaf(node.position, dir, node.leaf->size);
        }
    });
}

void Renderer::draw_snapshot(const TickSnapshot& snapshot) {
    // Build id -> index map for parent lookup
    std::unordered_map<uint32_t, size_t> id_to_idx;
    for (size_t i = 0; i < snapshot.nodes.size(); i++) {
        id_to_idx[snapshot.nodes[i].id] = i;
    }

    for (const auto& ns : snapshot.nodes) {
        if (ns.parent_id == UINT32_MAX) continue;

        auto it = id_to_idx.find(ns.parent_id);
        if (it == id_to_idx.end()) continue;
        const auto& parent = snapshot.nodes[it->second];

        glm::vec3 color;
        if (ns.type == NodeType::STEM) {
            color = glm::vec3(0.45f, 0.3f, 0.15f);
        } else {
            color = glm::vec3(0.35f, 0.2f, 0.1f);
        }

        draw_cylinder(parent.position, ns.position,
                      parent.radius, ns.radius, color);

        if (ns.has_leaf) {
            glm::vec3 dir = glm::normalize(ns.position - parent.position);
            draw_leaf(ns.position, dir, ns.leaf_size);
        }
    }
}

void Renderer::draw_cylinder(glm::vec3 start, glm::vec3 end,
                              float r_start, float r_end,
                              glm::vec3 color, int segments) {
    glm::vec3 axis = end - start;
    float height = glm::length(axis);
    if (height < 0.0001f) return;
    glm::vec3 axis_norm = axis / height;

    // Find perpendicular vectors
    glm::vec3 perp1;
    if (std::abs(axis_norm.y) < 0.9f) {
        perp1 = glm::normalize(glm::cross(axis_norm, glm::vec3(0.0f, 1.0f, 0.0f)));
    } else {
        perp1 = glm::normalize(glm::cross(axis_norm, glm::vec3(1.0f, 0.0f, 0.0f)));
    }
    glm::vec3 perp2 = glm::cross(axis_norm, perp1);

    // Generate vertices
    std::vector<float> vertices;
    for (int i = 0; i < segments; i++) {
        float a0 = (2.0f * 3.14159f * i) / segments;
        float a1 = (2.0f * 3.14159f * (i + 1)) / segments;
        float c0 = std::cos(a0), s0 = std::sin(a0);
        float c1 = std::cos(a1), s1 = std::sin(a1);

        glm::vec3 n0 = perp1 * c0 + perp2 * s0;
        glm::vec3 n1 = perp1 * c1 + perp2 * s1;

        glm::vec3 b0 = start + n0 * r_start;
        glm::vec3 b1 = start + n1 * r_start;
        glm::vec3 t0 = end + n0 * r_end;
        glm::vec3 t1 = end + n1 * r_end;

        // Two triangles per quad
        auto push = [&](glm::vec3 pos, glm::vec3 norm) {
            vertices.push_back(pos.x); vertices.push_back(pos.y); vertices.push_back(pos.z);
            vertices.push_back(norm.x); vertices.push_back(norm.y); vertices.push_back(norm.z);
            vertices.push_back(color.x); vertices.push_back(color.y); vertices.push_back(color.z);
        };

        push(b0, n0); push(t0, n0); push(b1, n1);
        push(b1, n1); push(t0, n0); push(t1, n1);
    }

    // Upload and draw
    uint32_t vao, vbo;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STREAM_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);

    glDrawArrays(GL_TRIANGLES, 0, static_cast<int>(vertices.size() / 9));

    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
}

void Renderer::draw_leaf(glm::vec3 position, glm::vec3 direction, float size) {
    // Flat quad oriented outward from branch
    glm::vec3 perp;
    if (std::abs(direction.y) < 0.9f) {
        perp = glm::normalize(glm::cross(direction, glm::vec3(0.0f, 1.0f, 0.0f)));
    } else {
        perp = glm::normalize(glm::cross(direction, glm::vec3(1.0f, 0.0f, 0.0f)));
    }
    glm::vec3 up = glm::normalize(glm::cross(perp, direction));

    float half = size * 0.5f;
    glm::vec3 p0 = position + perp * half + up * half;
    glm::vec3 p1 = position - perp * half + up * half;
    glm::vec3 p2 = position - perp * half - up * half;
    glm::vec3 p3 = position + perp * half - up * half;
    glm::vec3 normal = glm::normalize(glm::cross(p1 - p0, p3 - p0));
    glm::vec3 color(0.2f, 0.6f, 0.15f); // green

    float vertices[] = {
        p0.x, p0.y, p0.z, normal.x, normal.y, normal.z, color.x, color.y, color.z,
        p1.x, p1.y, p1.z, normal.x, normal.y, normal.z, color.x, color.y, color.z,
        p2.x, p2.y, p2.z, normal.x, normal.y, normal.z, color.x, color.y, color.z,
        p0.x, p0.y, p0.z, normal.x, normal.y, normal.z, color.x, color.y, color.z,
        p2.x, p2.y, p2.z, normal.x, normal.y, normal.z, color.x, color.y, color.z,
        p3.x, p3.y, p3.z, normal.x, normal.y, normal.z, color.x, color.y, color.z,
    };

    uint32_t vao, vbo;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STREAM_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);

    glDrawArrays(GL_TRIANGLES, 0, 6);

    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
}

void Renderer::setup_ground() {
    float size = 20.0f;
    glm::vec3 color(0.4f, 0.35f, 0.25f); // dirt
    glm::vec3 normal(0.0f, 1.0f, 0.0f);
    float verts[] = {
        -size, 0, -size, normal.x, normal.y, normal.z, color.x, color.y, color.z,
         size, 0, -size, normal.x, normal.y, normal.z, color.x, color.y, color.z,
         size, 0,  size, normal.x, normal.y, normal.z, color.x, color.y, color.z,
        -size, 0, -size, normal.x, normal.y, normal.z, color.x, color.y, color.z,
         size, 0,  size, normal.x, normal.y, normal.z, color.x, color.y, color.z,
        -size, 0,  size, normal.x, normal.y, normal.z, color.x, color.y, color.z,
    };

    glGenVertexArrays(1, &ground_vao_);
    glGenBuffers(1, &ground_vbo_);
    glBindVertexArray(ground_vao_);
    glBindBuffer(GL_ARRAY_BUFFER, ground_vbo_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(verts), verts, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);
}

void Renderer::draw_ground() {
    glBindVertexArray(ground_vao_);
    glDrawArrays(GL_TRIANGLES, 0, 6);
}

void Renderer::end_frame() {
    glfwSwapBuffers(window_);
    glfwPollEvents();
}

} // namespace botany
```

- [ ] **Step 5: Build to verify compilation**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu)
```

Expected: All targets compile. No runtime test — visual verification in Task 10.

- [ ] **Step 6: Commit**

```bash
git add src/renderer/ shaders/
git commit -m "feat: OpenGL renderer with shader pipeline, orbit camera, cylinder drawing"
```

---

### Task 9: Real-Time App (Mode 1)

**Files:**
- Modify: `src/app_realtime.cpp`

- [ ] **Step 1: Implement app_realtime.cpp**

```cpp
// src/app_realtime.cpp
#include <GLFW/glfw3.h>
#include <iostream>
#include "engine/engine.h"
#include "renderer/renderer.h"

using namespace botany;

static bool mouse_pressed = false;
static double last_mouse_x = 0, last_mouse_y = 0;
static Renderer* g_renderer = nullptr;

static void mouse_button_callback(GLFWwindow*, int button, int action, int) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        mouse_pressed = (action == GLFW_PRESS);
    }
}

static void cursor_pos_callback(GLFWwindow* window, double x, double y) {
    if (mouse_pressed && g_renderer) {
        float dx = static_cast<float>(x - last_mouse_x);
        float dy = static_cast<float>(y - last_mouse_y);
        g_renderer->camera().rotate(dx, dy);
    }
    last_mouse_x = x;
    last_mouse_y = y;
}

static void scroll_callback(GLFWwindow*, double, double yoffset) {
    if (g_renderer) {
        g_renderer->camera().zoom(static_cast<float>(yoffset));
    }
}

int main() {
    Engine engine;
    Genome g = default_genome();
    PlantID plant_id = engine.create_plant(g, glm::vec3(0.0f));

    Renderer renderer;
    if (!renderer.init(1280, 800, "shaders")) {
        std::cerr << "Failed to initialize renderer" << std::endl;
        return 1;
    }
    g_renderer = &renderer;

    GLFWwindow* window = renderer.window();
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_pos_callback);
    glfwSetScrollCallback(window, scroll_callback);

    int ticks_per_frame = 1;
    bool paused = false;

    while (!glfwWindowShouldClose(window)) {
        // Input
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, true);
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
            paused = !paused;
        if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
            ticks_per_frame = std::min(ticks_per_frame + 1, 100);
        if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
            ticks_per_frame = std::max(ticks_per_frame - 1, 1);

        // Simulate
        if (!paused) {
            for (int i = 0; i < ticks_per_frame; i++) {
                engine.tick();
            }
        }

        // Render
        renderer.begin_frame();
        renderer.draw_ground();
        renderer.draw_plant(engine.get_plant(plant_id));
        renderer.end_frame();
    }

    renderer.shutdown();
    return 0;
}
```

- [ ] **Step 2: Build and run**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu)
./botany_realtime
```

Expected: Window opens with sky-blue background, brown ground plane. A plant starts growing from the center — stem extends upward, root extends downward. Branches emerge from axillary nodes over time. Mouse drag rotates camera, scroll zooms.

- [ ] **Step 3: Commit**

```bash
git add src/app_realtime.cpp
git commit -m "feat: real-time app — engine + renderer main loop with camera controls"
```

---

### Task 10: Headless Precompute App (Mode 2)

**Files:**
- Modify: `src/app_headless.cpp`

- [ ] **Step 1: Implement app_headless.cpp**

```cpp
// src/app_headless.cpp
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "engine/engine.h"
#include "serialization/serializer.h"

using namespace botany;

int main(int argc, char* argv[]) {
    int num_ticks = 200;
    std::string output_path = "recording.bin";

    if (argc >= 2) num_ticks = std::atoi(argv[1]);
    if (argc >= 3) output_path = argv[2];

    std::cout << "Running " << num_ticks << " ticks, saving to " << output_path << std::endl;

    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));

    std::ofstream file(output_path, std::ios::binary);
    if (!file) {
        std::cerr << "Failed to open " << output_path << std::endl;
        return 1;
    }

    save_recording_header(file, g, static_cast<uint32_t>(num_ticks));

    for (int i = 0; i < num_ticks; i++) {
        engine.tick();
        save_tick(file, engine, 0);

        if ((i + 1) % 50 == 0) {
            std::cout << "Tick " << (i + 1) << "/" << num_ticks
                      << " (" << engine.get_plant(0).node_count() << " nodes)" << std::endl;
        }
    }

    std::cout << "Done. Final node count: " << engine.get_plant(0).node_count() << std::endl;
    return 0;
}
```

- [ ] **Step 2: Build and run**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu)
./botany_headless 100 test_recording.bin
```

Expected: Prints progress every 50 ticks, creates `test_recording.bin`.

- [ ] **Step 3: Commit**

```bash
git add src/app_headless.cpp
git commit -m "feat: headless precompute app — runs engine and saves binary recording"
```

---

### Task 11: Playback Viewer App (Mode 2)

**Files:**
- Modify: `src/app_playback.cpp`

- [ ] **Step 1: Implement app_playback.cpp**

```cpp
// src/app_playback.cpp
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "serialization/serializer.h"
#include "renderer/renderer.h"

using namespace botany;

static bool mouse_pressed = false;
static double last_mouse_x = 0, last_mouse_y = 0;
static Renderer* g_renderer = nullptr;
static bool imgui_wants_mouse = false;

static void mouse_button_callback(GLFWwindow*, int button, int action, int) {
    if (imgui_wants_mouse) return;
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        mouse_pressed = (action == GLFW_PRESS);
    }
}

static void cursor_pos_callback(GLFWwindow*, double x, double y) {
    if (imgui_wants_mouse) return;
    if (mouse_pressed && g_renderer) {
        float dx = static_cast<float>(x - last_mouse_x);
        float dy = static_cast<float>(y - last_mouse_y);
        g_renderer->camera().rotate(dx, dy);
    }
    last_mouse_x = x;
    last_mouse_y = y;
}

static void scroll_callback(GLFWwindow*, double, double yoffset) {
    if (imgui_wants_mouse) return;
    if (g_renderer) {
        g_renderer->camera().zoom(static_cast<float>(yoffset));
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: botany_playback <recording.bin>" << std::endl;
        return 1;
    }

    std::string input_path = argv[1];

    // Load entire recording into memory
    std::ifstream file(input_path, std::ios::binary);
    if (!file) {
        std::cerr << "Failed to open " << input_path << std::endl;
        return 1;
    }

    RecordingHeader header = load_recording_header(file);
    std::vector<TickSnapshot> ticks;
    ticks.reserve(header.num_ticks);
    for (uint32_t i = 0; i < header.num_ticks; i++) {
        ticks.push_back(load_tick(file));
    }
    std::cout << "Loaded " << ticks.size() << " ticks" << std::endl;

    // Init renderer
    Renderer renderer;
    if (!renderer.init(1280, 800, "shaders")) {
        std::cerr << "Failed to initialize renderer" << std::endl;
        return 1;
    }
    g_renderer = &renderer;

    GLFWwindow* window = renderer.window();
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_pos_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // Init ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 410");
    ImGui::StyleColorsDark();

    int current_tick = 0;
    bool playing = false;
    int speed = 1;
    float timer = 0.0f;
    float last_time = static_cast<float>(glfwGetTime());

    while (!glfwWindowShouldClose(window)) {
        float now = static_cast<float>(glfwGetTime());
        float dt = now - last_time;
        last_time = now;

        if (playing) {
            timer += dt * speed * 30.0f; // 30 ticks per second base rate
            while (timer >= 1.0f && current_tick < static_cast<int>(ticks.size()) - 1) {
                current_tick++;
                timer -= 1.0f;
            }
        }

        // Render
        renderer.begin_frame();
        renderer.draw_ground();
        if (!ticks.empty()) {
            renderer.draw_snapshot(ticks[current_tick]);
        }

        // ImGui
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        imgui_wants_mouse = ImGui::GetIO().WantCaptureMouse;

        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_Once);
        ImGui::Begin("Playback Controls");
        ImGui::Text("Tick: %d / %d", current_tick + 1, static_cast<int>(ticks.size()));
        ImGui::Text("Nodes: %d", static_cast<int>(ticks[current_tick].nodes.size()));
        ImGui::SliderInt("##tick", &current_tick, 0, static_cast<int>(ticks.size()) - 1);
        if (ImGui::Button(playing ? "Pause" : "Play")) {
            playing = !playing;
        }
        ImGui::SameLine();
        ImGui::SliderInt("Speed", &speed, 1, 10);
        if (ImGui::Button("Reset")) {
            current_tick = 0;
            playing = false;
        }
        ImGui::End();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        renderer.end_frame();
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    renderer.shutdown();
    return 0;
}
```

- [ ] **Step 2: Build and run**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu)
./botany_headless 200 recording.bin
./botany_playback recording.bin
```

Expected: Playback window opens showing the recorded plant growth. ImGui panel shows tick counter and slider. Play/pause and speed controls work. Scrubbing the slider jumps to any point in the recording.

- [ ] **Step 3: Commit**

```bash
git add src/app_playback.cpp
git commit -m "feat: playback viewer with ImGui controls and tick scrubbing"
```

---

### Task 12: Integration Test and Tuning

**Files:**
- No new files — this task validates the full system and tunes default genome values.

- [ ] **Step 1: Run the full test suite**

```bash
cd botany/build && cmake --build . -j$(sysctl -n hw.ncpu) && ctest --output-on-failure
```

Expected: All tests pass.

- [ ] **Step 2: Run the real-time app and observe plant growth**

```bash
cd botany/build && ./botany_realtime
```

Observe for ~30 seconds:
- Does the stem grow upward? Does it look like a plant?
- Do branches emerge from the stem at regular intervals?
- Do roots grow downward with their own branching?
- Does the plant eventually develop lateral branches (axillary activation)?

If the growth looks wrong (too fast, too slow, no branching, immediate branching), adjust the default genome values in `src/engine/genome.h`. Common adjustments:
- If no branches activate: lower `auxin_threshold` or raise `cytokinin_production_rate`
- If all branches activate immediately: raise `auxin_production_rate` or lower `auxin_threshold`
- If plant grows too fast/slow: adjust `growth_rate`
- If internodes are too long/short: adjust `max_internode_length` and `internode_spacing`

- [ ] **Step 3: Run headless + playback pipeline**

```bash
cd botany/build
./botany_headless 300 tuning.bin
./botany_playback tuning.bin
```

Scrub through the timeline. Verify the full growth history plays back correctly.

- [ ] **Step 4: Commit any genome tuning changes**

```bash
git add -A
git commit -m "tune: adjust default genome values for realistic growth behavior"
```
