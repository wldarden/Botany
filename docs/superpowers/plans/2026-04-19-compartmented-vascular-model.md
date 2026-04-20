# Compartmented Vascular Model Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the current single-pool chemical model and pairwise-Jacobi vascular transport with a compartmented LocalEnv + phloem + xylem structure driven by sub-stepped Münch pressure-flow.

**Architecture:** Every `Node` owns a `LocalEnv` (parenchyma chemicals). `StemNode` and `RootNode` additionally own a `TransportPool` for phloem (sugar) and one for xylem (water + cytokinin). Specialty nodes (leaves, apicals, root apicals) have only `local_env` and interact with their nearest ancestor's transport pools. A new `vascular_sub_stepped()` function replaces `vascular_transport()` — it runs a fixed N iterations of inject-radial-extract-Jacobi per tick. Plant tick order flips: metabolism first, vascular after.

**Tech Stack:** C++17, Catch2 tests, CMake build via `/usr/local/bin/cmake --build build`, tests via `./build/botany_tests`.

**Spec:** See [`../specs/2026-04-19-compartmented-vascular-model-design.md`](../specs/2026-04-19-compartmented-vascular-model-design.md) for full context and design rationale.

---

## File Structure

### New files

- `src/engine/compartments.h` — `LocalEnv` and `TransportPool` struct definitions. Small header, included by node types.
- `src/engine/vascular_sub_stepped.h` / `vascular_sub_stepped.cpp` — new vascular algorithm. During Phase D it coexists with the old `vascular.cpp`; after Phase F the old one goes away.
- `tests/test_compartments.cpp` — unit tests for `LocalEnv` and `TransportPool` structs, plus compartment-invariant tests.
- `tests/test_vascular_sub_stepped.cpp` — tests specific to the new algorithm (mass conservation, pressure propagation, distance-dependent supply, etc).

### Modified files

- `src/engine/node/node.h` / `node.cpp` — `local_env` replaces the raw `chemicals` map; add `local()` accessor; add virtual `phloem()` and `xylem()` returning `nullptr`; add `nearest_phloem_upstream()` and `nearest_xylem_upstream()` walk-up helpers.
- `src/engine/node/stem_node.h` / `stem_node.cpp` — add `TransportPool phloem_pool_`, `TransportPool xylem_pool_`; override `phloem()` and `xylem()`.
- `src/engine/node/root_node.h` / `root_node.cpp` — same as stem.
- `src/engine/genome.h` — add radial permeability and reserve/target params (see Task 13 for exact list).
- `src/engine/world_params.h` — add `vascular_substeps`.
- `src/engine/plant.cpp` — flip tick order in `Plant::tick()` to metabolism-first-vascular-after; eventually call new vascular function instead of old.
- `src/engine/vascular.h` / `vascular.cpp` — eventually trimmed down: `has_vasculature`, `phloem_resolve`, `xylem_resolve`, and the old `vascular_transport()` get deleted after cutover.
- **All files calling `node.chemical(id)`** — mechanical sweep to `node.local().chemical(id)`. Estimated ~30 files across `src/engine/`, `src/evolution/`, `src/serialization/`, and tests. Covered in Task 4.
- `CLAUDE.md` — Phase G rewrite of transport/tick sections.

### Phases

The plan is organized into 7 phases (A through G) mirroring the spec's migration plan. Each phase ends with tests green:

- **Phase A** (Tasks 1–3): Add compartment types, additive only.
- **Phase B** (Task 4): Mechanical sweep of chemical-access API.
- **Phase C** (Tasks 5–7): Pool members on StemNode/RootNode; walk-up helpers.
- **Phase D** (Tasks 8–19): Implement new vascular algorithm alongside old.
- **Phase E** (Task 20): Cutover — swap Plant::tick to use new algorithm.
- **Phase F** (Tasks 21–22): Remove dead code.
- **Phase G** (Task 23): Update CLAUDE.md.

---

## Phase A — Additive compartment types

### Task 1: Create `LocalEnv` and `TransportPool` types

**Files:**
- Create: `src/engine/compartments.h`
- Test: `tests/test_compartments.cpp`

- [ ] **Step 1: Write the failing test**

Create `tests/test_compartments.cpp`:

```cpp
#include <catch2/catch_test_macros.hpp>
#include "engine/compartments.h"

using namespace botany;

TEST_CASE("LocalEnv stores and reads chemicals by ChemicalID", "[compartments]") {
    LocalEnv env;
    env.chemical(ChemicalID::Sugar) = 1.5f;
    REQUIRE(env.chemical(ChemicalID::Sugar) == 1.5f);
    REQUIRE(env.chemical(ChemicalID::Water) == 0.0f);
}

TEST_CASE("TransportPool stores and reads chemicals by ChemicalID", "[compartments]") {
    TransportPool pool;
    pool.chemical(ChemicalID::Sugar) = 2.0f;
    REQUIRE(pool.chemical(ChemicalID::Sugar) == 2.0f);
    REQUIRE(pool.chemical(ChemicalID::Water) == 0.0f);
}

TEST_CASE("LocalEnv const accessor returns 0 for missing chemicals", "[compartments]") {
    const LocalEnv env;
    REQUIRE(env.chemical(ChemicalID::Sugar) == 0.0f);
}
```

Add the test to `CMakeLists.txt`'s `botany_tests` target (grep for existing `test_genome.cpp` entry and add `tests/test_compartments.cpp` alongside it).

- [ ] **Step 2: Run test — verify fail**

```bash
cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -5
```
Expected: compile error — `compartments.h` does not exist.

- [ ] **Step 3: Create `src/engine/compartments.h`**

```cpp
#pragma once

#include <unordered_map>
#include "engine/chemical/chemical.h"

// Hash for ChemicalID (duplicate of node.h's definition, but compartments.h
// needs to be includable standalone by headers/tests that don't pull in node.h)
namespace std {
    template<>
    struct hash<botany::ChemicalID> {
        std::size_t operator()(botany::ChemicalID id) const noexcept {
            return static_cast<std::size_t>(id);
        }
    };
}

namespace botany {

// Per-node local compartment.  Every node owns one.  Holds the chemicals the
// node itself uses for metabolism, growth, and signaling — i.e. everything
// that is NOT in the long-distance transport stream.
struct LocalEnv {
    std::unordered_map<ChemicalID, float> chemicals;

    float& chemical(ChemicalID id) { return chemicals[id]; }
    float chemical(ChemicalID id) const {
        auto it = chemicals.find(id);
        return it != chemicals.end() ? it->second : 0.0f;
    }
};

// Per-stem/root vascular conduit.  A StemNode or RootNode owns one phloem and
// one xylem TransportPool — representing the sieve tubes (phloem) and vessel
// elements (xylem) in that segment.  Specialty nodes (leaves, meristems) do
// not own any TransportPool.
struct TransportPool {
    std::unordered_map<ChemicalID, float> chemicals;

    float& chemical(ChemicalID id) { return chemicals[id]; }
    float chemical(ChemicalID id) const {
        auto it = chemicals.find(id);
        return it != chemicals.end() ? it->second : 0.0f;
    }
};

} // namespace botany
```

- [ ] **Step 4: Run test — verify pass**

```bash
cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[compartments]" 2>&1 | tail -5
```
Expected: PASS — three assertions pass.

- [ ] **Step 5: Run full suite to ensure nothing else broke**

```bash
./build/botany_tests 2>&1 | tail -3
```
Expected: `All tests passed`.

- [ ] **Step 6: Commit**

```bash
git add src/engine/compartments.h tests/test_compartments.cpp CMakeLists.txt
git commit -m "$(cat <<'EOF'
feat(vascular): introduce LocalEnv and TransportPool types

Additive-only change.  New compartment types with identical chemical()
accessor as Node.  Neither is wired into any node yet.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md
EOF
)"
```

---

### Task 2: Migrate `Node::chemicals` to `Node::local_env`

**Files:**
- Modify: `src/engine/node/node.h`
- Modify: `src/engine/node/node.cpp`

**Rationale:** Replace the raw `std::unordered_map<ChemicalID, float> chemicals` member with a `LocalEnv local_env` member. The existing `chemical(id)` method delegates to `local_env.chemical(id)`, so all external callers (which only use `chemical(id)`) are unaffected. This is behavior-preserving.

- [ ] **Step 1: Edit `src/engine/node/node.h`**

Add include near the top (after the existing `engine/chemical/chemical.h`):

```cpp
#include "engine/compartments.h"
```

Remove the `hash<ChemicalID>` template specialization at lines 10–15 (it now lives in `compartments.h`).

Replace the existing chemical-storage lines (currently `std::unordered_map<ChemicalID, float> chemicals;` and both `chemical(id)` methods) with:

```cpp
    // Local compartment — this node's parenchyma chemicals (sugar, water,
    // auxin, cytokinin, gibberellin, stress — anything NOT in the phloem or
    // xylem transport stream).  All chemical access goes through local().
    LocalEnv local_env;

    LocalEnv& local() { return local_env; }
    const LocalEnv& local() const { return local_env; }

    // Legacy accessor — forwards to local().chemical(id).  Retained for the
    // duration of Phase B; all call sites migrate to node.local().chemical(id)
    // in Task 4.  Remove after sweep is complete.
    float& chemical(ChemicalID id) { return local_env.chemical(id); }
    float chemical(ChemicalID id) const { return local_env.chemical(id); }
```

- [ ] **Step 2: Build and run full suite**

```bash
cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests 2>&1 | tail -3
```
Expected: builds clean, `All tests passed`. No behavior change because `chemical(id)` still works identically.

- [ ] **Step 3: Commit**

```bash
git add src/engine/node/node.h
git commit -m "$(cat <<'EOF'
refactor(vascular): Node.chemicals → Node.local_env (additive)

Replaces the raw chemical map on Node with a LocalEnv instance.  The
existing chemical(id) accessor delegates to local_env.chemical(id), so
all external call sites are unaffected — behavior-preserving refactor.

Phase B (Task 4) migrates every node.chemical(id) call site to
node.local().chemical(id) and then the legacy accessor is removed.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md
EOF
)"
```

---

### Task 3: Add virtual `phloem()` and `xylem()` to `Node`

**Files:**
- Modify: `src/engine/node/node.h`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_compartments.cpp`:

```cpp
#include "engine/node/node.h"
// Include concrete subclasses so we can instantiate them in the test
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"

TEST_CASE("Base Node returns nullptr from phloem() and xylem()", "[compartments]") {
    // We can't instantiate Node directly (abstract in spirit), but we can
    // verify the virtual default by checking a LeafNode — which inherits
    // the default nullptr behavior in this task (overrides arrive in Task 5/6).
    LeafNode leaf(/* id */ 1, glm::vec3(0), /* radius */ 0.01f);
    REQUIRE(leaf.phloem() == nullptr);
    REQUIRE(leaf.xylem() == nullptr);
}
```

If `LeafNode` constructor signature differs, adjust to match its actual signature — grep `src/engine/node/tissues/leaf.h` for the constructor.

- [ ] **Step 2: Run test — verify fail**

```bash
cd /Users/wldarden/learning/botany && /usr/local/bin/cmake --build build 2>&1 | tail -5
```
Expected: compile error — `leaf.phloem()` no such member.

- [ ] **Step 3: Add the virtuals to `Node` in `src/engine/node/node.h`**

Below the existing `virtual ~Node() = default;` line and the tree-operations block, add:

```cpp
    // --- Compartment access ---
    // Default: nullptr.  Stem and Root override to return their conduit pools.
    // Callers must null-check the return before dereferencing.
    virtual TransportPool* phloem() { return nullptr; }
    virtual const TransportPool* phloem() const { return nullptr; }
    virtual TransportPool* xylem()  { return nullptr; }
    virtual const TransportPool* xylem()  const { return nullptr; }
```

- [ ] **Step 4: Run test — verify pass**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[compartments]" 2>&1 | tail -5
```
Expected: PASS.

- [ ] **Step 5: Run full suite**

```bash
./build/botany_tests 2>&1 | tail -3
```
Expected: `All tests passed`.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/node.h tests/test_compartments.cpp
git commit -m "$(cat <<'EOF'
feat(vascular): Node virtual phloem()/xylem() accessors

Default return nullptr.  StemNode and RootNode will override in Task 5/6
to return their actual TransportPool members.  Leaves and meristems
inherit the nullptr default — a call to leaf.phloem() is a loud,
immediately-visible bug rather than a hidden type confusion.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md
EOF
)"
```

---
## Phase B — Chemical access API migration

### Task 4: Sweep all `node.chemical(id)` call sites to `node.local().chemical(id)`

**Files:**
- Modify: every `.h` and `.cpp` under `src/` and `tests/` that calls `chemical(` on a node-type object. Start with the grep below to get the exact list.

**Rationale:** Explicit compartment access at every call site is self-documenting — future readers see which compartment each piece of code is touching. Still additive-only in behavior (the legacy `chemical(id)` accessor continues to forward to `local().chemical(id)` until this task is complete, at which point the legacy accessor is removed in the final step).

- [ ] **Step 1: Inventory the call sites**

```bash
cd /Users/wldarden/learning/botany && grep -rn "\.chemical(\|->chemical(" src/ tests/ --include="*.cpp" --include="*.h" | wc -l
```
Expected: a number in the hundreds. Get the file list:

```bash
grep -rln "\.chemical(\|->chemical(" src/ tests/ --include="*.cpp" --include="*.h"
```

- [ ] **Step 2: Rewrite each call site**

The mechanical transformation is:
- `node.chemical(id)` → `node.local().chemical(id)`
- `node->chemical(id)` → `node->local().chemical(id)`
- `n.chemical(ChemicalID::X)` → `n.local().chemical(ChemicalID::X)`

Watch out for these patterns that look similar but aren't `Node::chemical`:
- `diffusion_params.chemical` (field access on a struct — leave as-is)
- `compartments.h`'s own `LocalEnv::chemical(id)` or `TransportPool::chemical(id)` (these are the new pool accessors, leave as-is)
- Any `ChemicalID` enum access (leave as-is)

Do the sweep in batches by directory, one commit per batch so diffs stay reviewable:

Batch 1: `src/engine/node/` (the core node code)
Batch 2: `src/engine/` (vascular, plant, engine, sugar, hormone, etc.)
Batch 3: `src/evolution/`, `src/serialization/`, `src/renderer/` if applicable
Batch 4: `src/app_*.cpp`
Batch 5: `tests/`

After each batch:
```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests 2>&1 | tail -3
```
Expected after every batch: build clean, `All tests passed`.

- [ ] **Step 3: Commit after each batch**

Example for Batch 1:
```bash
git add src/engine/node/
git commit -m "refactor(vascular): migrate src/engine/node/ to node.local().chemical(id)"
```

Replace the path and commit message for each subsequent batch.

- [ ] **Step 4: Remove the legacy accessor from `Node`**

Once the sweep is complete and all batches committed, edit `src/engine/node/node.h` to delete the legacy forwarding accessor added in Task 2:

```cpp
    // DELETE THESE LINES:
    float& chemical(ChemicalID id) { return local_env.chemical(id); }
    float chemical(ChemicalID id) const { return local_env.chemical(id); }
```

- [ ] **Step 5: Build — catch any missed call sites**

```bash
/usr/local/bin/cmake --build build 2>&1 | grep "error:" | head -10
```

Any errors like "no member named 'chemical' in 'botany::StemNode'" are missed call sites. Fix each in place (rewrite to `node.local().chemical(id)`). Rebuild until clean.

- [ ] **Step 6: Run full suite**

```bash
./build/botany_tests 2>&1 | tail -3
```
Expected: `All tests passed`. No behavior change — only access-API rewrite.

- [ ] **Step 7: Commit removal**

```bash
git add src/engine/node/node.h
git commit -m "$(cat <<'EOF'
refactor(vascular): remove legacy Node::chemical() accessor

All call sites now use node.local().chemical(id) explicitly, so the
forwarding accessor added in Task 2 is no longer needed.  Type checker
now enforces that chemical access goes through the local compartment.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md
EOF
)"
```

---


## Phase C — Pool members and walk-up helpers

### Task 5: Add `TransportPool` members to `StemNode` and override accessors

**Files:**
- Modify: `src/engine/node/stem_node.h`
- Modify: `src/engine/node/stem_node.cpp` (if any initialization needed)
- Modify: `tests/test_compartments.cpp`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_compartments.cpp`:

```cpp
TEST_CASE("StemNode exposes non-null phloem() and xylem()", "[compartments]") {
    StemNode stem(/* id */ 42, glm::vec3(0), /* radius */ 0.015f);
    REQUIRE(stem.phloem() != nullptr);
    REQUIRE(stem.xylem()  != nullptr);

    // The pools start empty.
    REQUIRE(stem.phloem()->chemical(ChemicalID::Sugar) == 0.0f);
    REQUIRE(stem.xylem()->chemical(ChemicalID::Water)  == 0.0f);
}

TEST_CASE("StemNode phloem and xylem are independent pools", "[compartments]") {
    StemNode stem(1, glm::vec3(0), 0.015f);
    stem.phloem()->chemical(ChemicalID::Sugar) = 5.0f;
    stem.xylem()->chemical(ChemicalID::Water)  = 3.0f;
    REQUIRE(stem.phloem()->chemical(ChemicalID::Sugar) == 5.0f);
    REQUIRE(stem.xylem()->chemical(ChemicalID::Water)  == 3.0f);
    REQUIRE(stem.phloem()->chemical(ChemicalID::Water) == 0.0f);  // water not in phloem
    REQUIRE(stem.xylem()->chemical(ChemicalID::Sugar)  == 0.0f);  // sugar not in xylem
}
```

If the `StemNode` constructor takes additional args, adjust — grep `stem_node.h` for its signature.

- [ ] **Step 2: Run test — verify fail**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -5
```
Expected: either compile error (no `phloem()` override) or test failure (pools return nullptr from base Node).

- [ ] **Step 3: Add pool members to `src/engine/node/stem_node.h`**

Include at the top (after existing includes):
```cpp
#include "engine/compartments.h"
```

In the class body, add private members and override public accessors:

```cpp
class StemNode : public Node {
public:
    StemNode(uint32_t id, glm::vec3 position, float radius);
    // ... existing methods ...

    // Compartment overrides — StemNode owns both pool types.
    TransportPool* phloem() override { return &phloem_pool_; }
    const TransportPool* phloem() const override { return &phloem_pool_; }
    TransportPool* xylem()  override { return &xylem_pool_;  }
    const TransportPool* xylem()  const override { return &xylem_pool_;  }

private:
    TransportPool phloem_pool_;
    TransportPool xylem_pool_;
};
```

- [ ] **Step 4: Run test — verify pass**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[compartments]" 2>&1 | tail -5
```
Expected: new StemNode tests PASS.

- [ ] **Step 5: Run full suite**

```bash
./build/botany_tests 2>&1 | tail -3
```
Expected: `All tests passed`. No behavior change — the pools exist but nothing reads or writes them yet.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/stem_node.h tests/test_compartments.cpp
git commit -m "$(cat <<'EOF'
feat(vascular): StemNode owns phloem and xylem TransportPools

Two TransportPool members on every StemNode, with phloem()/xylem()
overrides returning their addresses.  No code reads or writes these
pools yet — they're latent until the new vascular algorithm lands in
Phase D.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md
EOF
)"
```

---

### Task 6: Add `TransportPool` members to `RootNode` and override accessors

**Files:**
- Modify: `src/engine/node/root_node.h`
- Modify: `tests/test_compartments.cpp`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_compartments.cpp`:

```cpp
TEST_CASE("RootNode exposes non-null phloem() and xylem()", "[compartments]") {
    RootNode root(1, glm::vec3(0, -0.05f, 0), 0.015f);
    REQUIRE(root.phloem() != nullptr);
    REQUIRE(root.xylem()  != nullptr);
}

TEST_CASE("LeafNode/ApicalNode/RootApicalNode return nullptr pools", "[compartments]") {
    LeafNode leaf(1, glm::vec3(0), 0.01f);
    ApicalNode apical(2, glm::vec3(0), 0.01f);
    RootApicalNode root_apical(3, glm::vec3(0, -0.05f, 0), 0.01f);

    REQUIRE(leaf.phloem()        == nullptr);
    REQUIRE(leaf.xylem()         == nullptr);
    REQUIRE(apical.phloem()      == nullptr);
    REQUIRE(apical.xylem()       == nullptr);
    REQUIRE(root_apical.phloem() == nullptr);
    REQUIRE(root_apical.xylem()  == nullptr);
}
```

Adjust constructor signatures if necessary — grep each subclass header for its constructor.

- [ ] **Step 2: Run test — verify fail**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[compartments]" 2>&1 | tail -10
```
Expected: the RootNode test fails because RootNode doesn't override phloem()/xylem() yet.

- [ ] **Step 3: Add pool members to `src/engine/node/root_node.h`**

Same pattern as StemNode — include `engine/compartments.h`, add `TransportPool phloem_pool_` and `TransportPool xylem_pool_` as private members, override the four accessors (2 mutable + 2 const).

- [ ] **Step 4: Run tests — verify pass**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[compartments]" 2>&1 | tail -5
```
Expected: all compartment tests PASS.

- [ ] **Step 5: Run full suite**

```bash
./build/botany_tests 2>&1 | tail -3
```
Expected: `All tests passed`.

- [ ] **Step 6: Commit**

```bash
git add src/engine/node/root_node.h tests/test_compartments.cpp
git commit -m "$(cat <<'EOF'
feat(vascular): RootNode owns phloem and xylem TransportPools

Mirrors the StemNode change — roots are the source side of xylem and
have their own section of the phloem network, so they need their own
pool members.  LeafNode, ApicalNode, and RootApicalNode correctly
inherit the nullptr default (verified by new test).

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md
EOF
)"
```

---

### Task 7: Walk-up helpers `nearest_phloem_upstream()` and `nearest_xylem_upstream()`

**Files:**
- Modify: `src/engine/node/node.h`
- Modify: `src/engine/node/node.cpp`
- Modify: `tests/test_compartments.cpp`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_compartments.cpp`:

```cpp
TEST_CASE("nearest_phloem_upstream finds parent stem", "[compartments][walkup]") {
    StemNode stem(1, glm::vec3(0), 0.015f);
    LeafNode leaf(2, glm::vec3(0.05f, 0, 0), 0.01f);
    stem.add_child(&leaf);
    leaf.parent = &stem;
    REQUIRE(leaf.nearest_phloem_upstream() == stem.phloem());
}

TEST_CASE("nearest_phloem_upstream walks past non-conduit ancestors", "[compartments][walkup]") {
    StemNode stem(1, glm::vec3(0), 0.015f);
    ApicalNode apical(2, glm::vec3(0.05f, 0, 0), 0.01f);
    LeafNode leaf(3, glm::vec3(0.08f, 0, 0), 0.01f);
    stem.add_child(&apical);
    apical.parent = &stem;
    apical.add_child(&leaf);
    leaf.parent = &apical;
    // Leaf's direct parent is apical (no phloem).  Walk-up must skip
    // past the apical and find the stem.
    REQUIRE(leaf.nearest_phloem_upstream() == stem.phloem());
}

TEST_CASE("nearest_phloem_upstream returns nullptr when no upstream conduit", "[compartments][walkup]") {
    LeafNode orphan(1, glm::vec3(0), 0.01f);
    // No parent set.
    REQUIRE(orphan.nearest_phloem_upstream() == nullptr);
}

TEST_CASE("nearest_xylem_upstream finds parent root", "[compartments][walkup]") {
    RootNode root(1, glm::vec3(0, -0.05f, 0), 0.015f);
    RootApicalNode ra(2, glm::vec3(0, -0.10f, 0), 0.01f);
    root.add_child(&ra);
    ra.parent = &root;
    REQUIRE(ra.nearest_xylem_upstream() == root.xylem());
}
```

- [ ] **Step 2: Run tests — verify fail**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -5
```
Expected: compile error — no such method.

- [ ] **Step 3: Declare in `src/engine/node/node.h`** (add to public section):

```cpp
    // Walk-up: find the nearest ancestor that owns a TransportPool of the
    // given kind.  Used by leaves, apicals, and root apicals (which have no
    // pool of their own) to locate the parent conduit they load into or
    // unload from.  Returns nullptr if no ancestor has a matching pool.
    TransportPool* nearest_phloem_upstream();
    const TransportPool* nearest_phloem_upstream() const;
    TransportPool* nearest_xylem_upstream();
    const TransportPool* nearest_xylem_upstream() const;
```

- [ ] **Step 4: Implement in `src/engine/node/node.cpp`** (add near the end of the file):

```cpp
TransportPool* Node::nearest_phloem_upstream() {
    for (Node* n = parent; n != nullptr; n = n->parent) {
        if (auto p = n->phloem()) return p;
    }
    return nullptr;
}

const TransportPool* Node::nearest_phloem_upstream() const {
    for (const Node* n = parent; n != nullptr; n = n->parent) {
        if (auto p = n->phloem()) return p;
    }
    return nullptr;
}

TransportPool* Node::nearest_xylem_upstream() {
    for (Node* n = parent; n != nullptr; n = n->parent) {
        if (auto p = n->xylem()) return p;
    }
    return nullptr;
}

const TransportPool* Node::nearest_xylem_upstream() const {
    for (const Node* n = parent; n != nullptr; n = n->parent) {
        if (auto p = n->xylem()) return p;
    }
    return nullptr;
}
```

- [ ] **Step 5: Run tests — verify pass**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[walkup]" 2>&1 | tail -5
```
Expected: four walk-up tests PASS.

- [ ] **Step 6: Run full suite**

```bash
./build/botany_tests 2>&1 | tail -3
```
Expected: `All tests passed`.

- [ ] **Step 7: Commit**

```bash
git add src/engine/node/node.h src/engine/node/node.cpp tests/test_compartments.cpp
git commit -m "$(cat <<'EOF'
feat(vascular): walk-up helpers for finding upstream conduits

Adds Node::nearest_phloem_upstream() and nearest_xylem_upstream() that
walk the parent chain until finding the first ancestor with the matching
TransportPool.  Used by specialty nodes (leaves, apicals, root apicals)
to locate the conduit they should load into or unload from, since their
direct parent may also be a specialty node (e.g., leaf parented to an
apical).

O(depth) per call, recomputed on demand.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md
EOF
)"
```

---

## Phase D — New vascular algorithm (alongside old)

Phase D builds `vascular_sub_stepped()` from scratch. During this phase, the old `vascular_transport()` is still wired into `Plant::tick()` and continues running — the old algorithm keeps the tests green while the new algorithm is built, tested, and verified. At the end of Phase D, the new function exists and has its own passing test suite, but isn't called from the main tick yet. Phase E is the cutover.

### Task 8: Add genome and world params for the new model

**Files:**
- Modify: `src/engine/genome.h`
- Modify: `src/engine/world_params.h`
- Modify: `src/evolution/genome_bridge.cpp` (register new genes for evolution)

- [ ] **Step 1: Add genome fields** in `src/engine/genome.h` under a new "Compartmented vascular" section (find an existing section boundary and add after canalization or vascular params):

```cpp
    // --- Compartmented vascular model (2026-04-19) ---
    // Radial flow permeability curve (per chemical class): perm(r) = base ×
    // (floor + (1 - floor) / (1 + (r / half_radius)²)).  Young thin stems
    // leak freely; mature thick trunks asymptote to `floor × base`.  See
    // spec section 6 for full rationale.

    // Phloem radial (sugar between stem's own local_env and own phloem):
    float base_radial_permeability_sugar;   // max permeability at r = 0
    float radial_floor_fraction_sugar;      // asymptote as r → ∞, fraction of base
    float radial_half_radius_sugar;         // dm — curve inflection

    // Xylem radial (water and cytokinin between stem's own local_env and own xylem):
    float base_radial_permeability_water;
    float radial_floor_fraction_water;
    float radial_half_radius_water;

    // Pipe cross-section fractions (what portion of stem cross-section is
    // sieve tubes vs vessel elements):
    float phloem_fraction;                  // 0–1, portion of π·r² that's phloem
    float xylem_fraction;                   // 0–1, portion of π·r² that's xylem

    // Source/sink targets:
    float leaf_reserve_fraction_sugar;      // leaf keeps this fraction of its
                                            // local sugar before loading to phloem
    float meristem_sink_target_fraction;    // meristem refills to this fraction
                                            // of sugar_cap per tick
    float leaf_turgor_target_fraction;      // leaves refill to this fraction of
                                            // water_cap via xylem pull
```

- [ ] **Step 2: Add defaults** in `default_genome()`:

```cpp
        // Compartmented vascular defaults — starting values from spec section 6.
        .base_radial_permeability_sugar = 1.0f,
        .radial_floor_fraction_sugar    = 0.1f,
        .radial_half_radius_sugar       = 0.3f,   // dm
        .base_radial_permeability_water = 1.0f,
        .radial_floor_fraction_water    = 0.1f,
        .radial_half_radius_water       = 0.3f,
        .phloem_fraction                = 0.05f,  // 5% of cross-section is sieve tubes
        .xylem_fraction                 = 0.2f,   // 20% of cross-section is vessels
        .leaf_reserve_fraction_sugar    = 0.3f,   // matches existing phloem_reserve_fraction value
        .meristem_sink_target_fraction  = 0.05f,  // matches existing meristem_sink_fraction
        .leaf_turgor_target_fraction    = 0.7f,   // leaves aim to fill to 70% of water_cap
```

- [ ] **Step 3: Add world param** in `src/engine/world_params.h`:

```cpp
    uint32_t vascular_substeps = 25;   // N — number of sub-steps in the
                                       // sub-stepped vascular loop.  Each
                                       // sub-step propagates the pressure
                                       // wave by ~1 hop.  Plants with chains
                                       // longer than this will show
                                       // distance-dependent apical supply
                                       // (intentional — real hydraulic limit).
```

- [ ] **Step 4: Register new genes in `src/evolution/genome_bridge.cpp`**

Find the existing `build_genome_template()` and add entries for each new field, following the pattern of existing canalization entries. Use reasonable ranges — e.g., `radial_floor_fraction_sugar` range [0.01, 0.5], `radial_half_radius_sugar` range [0.05, 2.0], etc. Match the mutation-strength pattern used for nearby params.

- [ ] **Step 5: Build and run full suite**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests 2>&1 | tail -3
```
Expected: build clean, `All tests passed`. New params exist but aren't read anywhere yet — no behavior change.

- [ ] **Step 6: Commit**

```bash
git add src/engine/genome.h src/engine/world_params.h src/evolution/genome_bridge.cpp
git commit -m "$(cat <<'EOF'
feat(vascular): add genome and world params for compartmented model

Adds the radial-permeability curve params (per chemical class), pipe
cross-section fractions, source/sink target fractions, and the
vascular_substeps world param.  None of these are read by any code
yet — they're wired up in Tasks 9-19.  Registering them now lets
evolution experiments pick up the new genes in any serialized genome.

Defaults match spec section 6 starting values.  Ranges for mutation
chosen to match neighboring params in the genome bridge.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md
EOF
)"
```

---

### Task 9: Create `vascular_sub_stepped.h` with empty function stub

**Files:**
- Create: `src/engine/vascular_sub_stepped.h`
- Create: `src/engine/vascular_sub_stepped.cpp`
- Create: `tests/test_vascular_sub_stepped.cpp`
- Modify: `CMakeLists.txt` (add new files to botany_engine target and the new test file to botany_tests target)

- [ ] **Step 1: Write the failing test**

`tests/test_vascular_sub_stepped.cpp`:

```cpp
#include <catch2/catch_test_macros.hpp>
#include "engine/vascular_sub_stepped.h"
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/world_params.h"

using namespace botany;

TEST_CASE("vascular_sub_stepped exists and can be called", "[vascular_sub_stepped]") {
    Genome g = default_genome();
    WorldParams world;
    Plant plant(g, glm::vec3(0.0f));

    // Should not crash on an empty plant (just a seed).
    vascular_sub_stepped(plant, g, world);
    SUCCEED();
}
```

- [ ] **Step 2: Run — verify fail**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3
```
Expected: compile error — header does not exist.

- [ ] **Step 3: Create `src/engine/vascular_sub_stepped.h`**

```cpp
#pragma once

namespace botany {

class Plant;
struct Genome;
struct WorldParams;

// Sub-stepped vascular transport.  Replaces the pairwise-Jacobi
// vascular_transport in Phase E.  For the duration of Phase D this function
// coexists with the old one but is not wired into Plant::tick() yet.
//
// Algorithm: N = world.vascular_substeps iterations of
//   1. Inject — sources push budget/N into parent's conduit
//   2. Radial flow — stem/root local_env ⇄ own phloem/xylem, radius-dependent
//      permeability
//   3. Extract — sinks pull budget/N from parent's conduit
//   4. Longitudinal Jacobi — one pass of neighbor pressure equalization
// See spec: docs/superpowers/specs/2026-04-19-compartmented-vascular-model-design.md
void vascular_sub_stepped(Plant& plant, const Genome& g, const WorldParams& world);

} // namespace botany
```

- [ ] **Step 4: Create `src/engine/vascular_sub_stepped.cpp`**

```cpp
#include "engine/vascular_sub_stepped.h"
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/world_params.h"

namespace botany {

void vascular_sub_stepped(Plant& /*plant*/, const Genome& /*g*/, const WorldParams& /*world*/) {
    // Empty stub — per-sub-step logic is added in Tasks 10-18.
}

} // namespace botany
```

- [ ] **Step 5: Add files to `CMakeLists.txt`**

Find the existing `botany_engine` target or equivalent (grep for `vascular.cpp`). Add `src/engine/vascular_sub_stepped.cpp` in the same list. Add `tests/test_vascular_sub_stepped.cpp` to the test target alongside existing test cpp files.

- [ ] **Step 6: Build and run test**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[vascular_sub_stepped]" 2>&1 | tail -5
```
Expected: PASS (trivial — just checks the function can be called).

- [ ] **Step 7: Commit**

```bash
git add src/engine/vascular_sub_stepped.h src/engine/vascular_sub_stepped.cpp tests/test_vascular_sub_stepped.cpp CMakeLists.txt
git commit -m "$(cat <<'EOF'
feat(vascular): scaffold for vascular_sub_stepped

Empty function stub plus trivial test that it can be called without
crashing.  Subsequent tasks flesh out the algorithm step by step.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md
EOF
)"
```

---

### Task 10: Radial permeability curve function

**Files:**
- Modify: `src/engine/vascular_sub_stepped.cpp`
- Modify: `tests/test_vascular_sub_stepped.cpp`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_vascular_sub_stepped.cpp`:

```cpp
#include "engine/genome.h"

// Expose the internal helper for testing via a forward declaration.
// (Declared inside anonymous namespace in the cpp; we re-declare here.)
namespace botany {
    float radial_permeability_sugar(float radius, const Genome& g);
    float radial_permeability_water(float radius, const Genome& g);
}

TEST_CASE("radial_permeability approaches base at r=0", "[vascular_sub_stepped][radial]") {
    Genome g = default_genome();
    float perm = radial_permeability_sugar(0.0f, g);
    REQUIRE(perm == g.base_radial_permeability_sugar);
}

TEST_CASE("radial_permeability at half-radius is between base and floor", "[vascular_sub_stepped][radial]") {
    Genome g = default_genome();
    float base   = g.base_radial_permeability_sugar;
    float floor  = g.radial_floor_fraction_sugar;
    float r_half = g.radial_half_radius_sugar;
    float perm = radial_permeability_sugar(r_half, g);
    // At r = r_half, denominator = 1 + 1 = 2, so (1 - floor)/2 remains.
    // perm = base × (floor + (1 - floor) / 2)
    float expected = base * (floor + (1.0f - floor) * 0.5f);
    REQUIRE(perm == Catch::Approx(expected).margin(1e-5f));
}

TEST_CASE("radial_permeability asymptotes to floor at large r", "[vascular_sub_stepped][radial]") {
    Genome g = default_genome();
    float perm = radial_permeability_sugar(100.0f, g);
    float expected_floor = g.base_radial_permeability_sugar * g.radial_floor_fraction_sugar;
    REQUIRE(perm == Catch::Approx(expected_floor).margin(0.01f));
}
```

- [ ] **Step 2: Run — verify fail**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -5
```
Expected: compile error — `radial_permeability_sugar` not defined.

- [ ] **Step 3: Implement in `src/engine/vascular_sub_stepped.cpp`**

Add at top of file, inside the `namespace botany`:

```cpp
#include "engine/genome.h"

namespace botany {

float radial_permeability_sugar(float radius, const Genome& g) {
    const float base   = g.base_radial_permeability_sugar;
    const float floor  = g.radial_floor_fraction_sugar;
    const float r_half = g.radial_half_radius_sugar;
    const float ratio  = radius / r_half;
    return base * (floor + (1.0f - floor) / (1.0f + ratio * ratio));
}

float radial_permeability_water(float radius, const Genome& g) {
    const float base   = g.base_radial_permeability_water;
    const float floor  = g.radial_floor_fraction_water;
    const float r_half = g.radial_half_radius_water;
    const float ratio  = radius / r_half;
    return base * (floor + (1.0f - floor) / (1.0f + ratio * ratio));
}
```

- [ ] **Step 4: Run test — verify pass**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[radial]" 2>&1 | tail -5
```
Expected: PASS.

- [ ] **Step 5: Run full suite + commit**

```bash
./build/botany_tests 2>&1 | tail -3
git add src/engine/vascular_sub_stepped.cpp tests/test_vascular_sub_stepped.cpp
git commit -m "feat(vascular): radial_permeability curve (inverse-saturation with floor)

Implements the radius-dependent radial permeability formula:
  perm(r) = base × (floor + (1 - floor) / (1 + (r / half_radius)²))

Young thin stems → full base permeability (leaky, grow fast).
Mature thick trunks → asymptote to base × floor (act as pipes).

Two variants (sugar and water) so each chemical class can evolve its
own curve.  Three evolvable genome params per variant.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md"
```

---

### Task 11: Phloem and xylem capacity helpers

**Files:**
- Modify: `src/engine/vascular_sub_stepped.cpp`
- Modify: `tests/test_vascular_sub_stepped.cpp`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_vascular_sub_stepped.cpp`:

```cpp
namespace botany {
    float phloem_capacity(const Node& n, const Genome& g);
    float xylem_capacity(const Node& n, const Genome& g);
}

TEST_CASE("phloem_capacity scales with r² × length × phloem_fraction", "[vascular_sub_stepped][capacity]") {
    Genome g = default_genome();
    StemNode stem(1, glm::vec3(0), 0.1f);  // radius 0.1 dm
    // We need to give the stem a length via offset magnitude.  Set offset to
    // a unit vector along Y so length = 1 dm.
    stem.offset = glm::vec3(0.0f, 1.0f, 0.0f);

    float expected = 3.14159f * 0.1f * 0.1f * 1.0f * g.phloem_fraction;
    REQUIRE(phloem_capacity(stem, g) == Catch::Approx(expected).margin(1e-4f));
}

TEST_CASE("phloem_capacity returns zero for a leaf (no phloem)", "[vascular_sub_stepped][capacity]") {
    Genome g = default_genome();
    LeafNode leaf(1, glm::vec3(0), 0.01f);
    REQUIRE(phloem_capacity(leaf, g) == 0.0f);
}
```

- [ ] **Step 2: Run — verify fail**

Expected: compile error, helpers not defined.

- [ ] **Step 3: Implement in `src/engine/vascular_sub_stepped.cpp`**

```cpp
#include "engine/node/node.h"
#include <glm/geometric.hpp>  // for glm::length
#include <cmath>

// Returns the length of the internode (distance from parent to this node).
// For the seed (no parent), returns 1.0 so the seed has a meaningful capacity.
static float node_length(const Node& n) {
    if (!n.parent) return 1.0f;  // seed: use unit length
    return glm::length(n.offset);
}

float phloem_capacity(const Node& n, const Genome& g) {
    if (!n.phloem()) return 0.0f;   // only stems/roots have phloem
    const float L = node_length(n);
    return 3.14159265358979f * n.radius * n.radius * L * g.phloem_fraction;
}

float xylem_capacity(const Node& n, const Genome& g) {
    if (!n.xylem()) return 0.0f;
    const float L = node_length(n);
    return 3.14159265358979f * n.radius * n.radius * L * g.xylem_fraction;
}
```

- [ ] **Step 4: Run test — verify pass**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[capacity]" 2>&1 | tail -5
```
Expected: PASS.

- [ ] **Step 5: Run full suite + commit**

```bash
./build/botany_tests 2>&1 | tail -3
git add src/engine/vascular_sub_stepped.cpp tests/test_vascular_sub_stepped.cpp
git commit -m "feat(vascular): phloem_capacity and xylem_capacity helpers

Capacity = π · r² · L · fraction.  Returns 0 for nodes without the
matching pool (leaves, meristems).

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md"
```

---

### Task 12: Plant-wide flattening helper

**Files:**
- Modify: `src/engine/vascular_sub_stepped.cpp`

**Rationale:** The inject/extract/radial/Jacobi steps all need to walk every node in the plant. Doing a recursive walk every sub-step is wasteful. Build a flat node list once per vascular call and reuse across all N sub-steps.

- [ ] **Step 1: Add flattening helper**

In `src/engine/vascular_sub_stepped.cpp`, add (below the helpers from Tasks 10–11):

```cpp
#include <vector>

namespace {

// Flat view of the plant for efficient sub-step iteration.  Populated once
// at the start of vascular_sub_stepped(), then reused across all N iterations.
struct FlatNodes {
    std::vector<botany::Node*> all;             // every node, DFS order
    std::vector<botany::Node*> conduits;        // only nodes with phloem/xylem
    std::vector<botany::Node*> leaves;          // LeafNode only
    std::vector<botany::Node*> roots;           // RootNode only
    std::vector<botany::Node*> apicals;         // ApicalNode only (includes SA)
    std::vector<botany::Node*> root_apicals;    // RootApicalNode only
};

void collect_dfs(botany::Node* n, FlatNodes& flat) {
    if (!n) return;
    flat.all.push_back(n);
    if (n->phloem() || n->xylem()) flat.conduits.push_back(n);
    switch (n->type) {
        case botany::NodeType::LEAF:        flat.leaves.push_back(n); break;
        case botany::NodeType::ROOT:        flat.roots.push_back(n); break;
        case botany::NodeType::APICAL:      flat.apicals.push_back(n); break;
        case botany::NodeType::ROOT_APICAL: flat.root_apicals.push_back(n); break;
        default: break;
    }
    for (botany::Node* child : n->children) collect_dfs(child, flat);
}

FlatNodes flatten(botany::Plant& plant) {
    FlatNodes flat;
    collect_dfs(plant.seed_mut(), flat);
    return flat;
}

} // anonymous namespace
```

- [ ] **Step 2: Build and run full suite**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests 2>&1 | tail -3
```
Expected: build clean, `All tests passed`. The helper isn't called anywhere yet — just a compilation check.

- [ ] **Step 3: Commit**

```bash
git add src/engine/vascular_sub_stepped.cpp
git commit -m "feat(vascular): flattening helper for sub-step iteration

Single DFS walk produces buckets for nodes with conduits, leaves, roots,
apicals, and root apicals.  Built once per vascular call and reused
across N sub-steps to avoid N×tree-walks.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md"
```

---

### Task 13: Budget computation

**Files:**
- Modify: `src/engine/vascular_sub_stepped.cpp`
- Modify: `tests/test_vascular_sub_stepped.cpp`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_vascular_sub_stepped.cpp`:

```cpp
namespace botany {
    struct VascularBudget {
        float sugar_supply;   // for this node, how much to inject into parent phloem this tick
        float sugar_demand;   // for this node, how much to pull from parent phloem this tick
        float water_demand;   // for this node, how much to pull from parent xylem this tick
        float cytokinin_supply; // for root apicals, how much to push into parent xylem
    };
    VascularBudget compute_budget(Node& n, const Genome& g, const WorldParams& world);
}

TEST_CASE("leaf budget: sugar_supply is local sugar above reserve", "[vascular_sub_stepped][budget]") {
    Genome g = default_genome();
    WorldParams world;
    LeafNode leaf(1, glm::vec3(0), 0.02f);
    // Assume leaf has 10 g sugar in local_env and sugar_cap is computed internally.
    // The reserve is leaf_reserve_fraction_sugar × sugar_cap.
    leaf.local().chemical(ChemicalID::Sugar) = 10.0f;
    VascularBudget b = compute_budget(leaf, g, world);
    REQUIRE(b.sugar_supply > 0.0f);
    REQUIRE(b.sugar_demand == 0.0f);
}

TEST_CASE("dormant apical budget: zero demand", "[vascular_sub_stepped][budget]") {
    Genome g = default_genome();
    WorldParams world;
    ApicalNode apical(1, glm::vec3(0), 0.01f);
    apical.active = false;  // dormant
    VascularBudget b = compute_budget(apical, g, world);
    REQUIRE(b.sugar_demand == 0.0f);
    REQUIRE(b.water_demand == 0.0f);
}

TEST_CASE("active apical below target has positive sugar_demand", "[vascular_sub_stepped][budget]") {
    Genome g = default_genome();
    WorldParams world;
    ApicalNode apical(1, glm::vec3(0), 0.01f);
    apical.active = true;
    apical.local().chemical(ChemicalID::Sugar) = 0.0f;
    VascularBudget b = compute_budget(apical, g, world);
    REQUIRE(b.sugar_demand > 0.0f);
}
```

(Field name `apical.active` assumed — grep `apical.h` for the actual activation flag name and adjust.)

- [ ] **Step 2: Run — verify fail**

Expected: compile error, `VascularBudget` undefined.

- [ ] **Step 3: Implement in `src/engine/vascular_sub_stepped.cpp`**

Define the struct in the header (`vascular_sub_stepped.h`) so it can be tested, and implement `compute_budget` in the cpp. Include `engine/sugar.h` and whatever header exposes `sugar_cap` / `water_cap`.

```cpp
// in vascular_sub_stepped.h, inside namespace botany:
struct VascularBudget {
    float sugar_supply   = 0.0f;
    float sugar_demand   = 0.0f;
    float water_demand   = 0.0f;
    float cytokinin_supply = 0.0f;
};

VascularBudget compute_budget(Node& n, const Genome& g, const WorldParams& world);
```

```cpp
// in vascular_sub_stepped.cpp, inside namespace botany:
#include "engine/sugar.h"     // sugar_cap
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"

VascularBudget compute_budget(Node& n, const Genome& g, const WorldParams& world) {
    VascularBudget b;
    const float sugar = n.local().chemical(ChemicalID::Sugar);
    const float water = n.local().chemical(ChemicalID::Water);
    const float cap_s = sugar_cap(n, g);
    const float cap_w = water_cap(n, g);

    switch (n.type) {
        case NodeType::LEAF: {
            const float reserve = g.leaf_reserve_fraction_sugar * cap_s;
            b.sugar_supply = std::max(0.0f, sugar - reserve);
            const float turgor_target = g.leaf_turgor_target_fraction * cap_w;
            b.water_demand = std::max(0.0f, turgor_target - water);
            break;
        }
        case NodeType::APICAL: {
            if (n.as_apical() && n.as_apical()->active) {
                const float target = g.meristem_sink_target_fraction * cap_s;
                b.sugar_demand = std::max(0.0f, target - sugar);
                const float turgor_target = g.leaf_turgor_target_fraction * cap_w;
                b.water_demand = std::max(0.0f, turgor_target - water);
            }
            break;
        }
        case NodeType::ROOT_APICAL: {
            if (n.as_root_apical() && n.as_root_apical()->active) {
                const float target = g.meristem_sink_target_fraction * cap_s;
                b.sugar_demand = std::max(0.0f, target - sugar);
                const float turgor_target = g.leaf_turgor_target_fraction * cap_w;
                b.water_demand = std::max(0.0f, turgor_target - water);
            }
            // Root apicals also push cytokinin into parent xylem.  They produce
            // cytokinin in their local_env during the tick phase; here we treat
            // all local cytokinin as surplus (they don't use it themselves).
            b.cytokinin_supply = n.local().chemical(ChemicalID::Cytokinin);
            break;
        }
        // STEM and ROOT get their sugar via radial flow, not direct injection.
        // No inject/extract budget needed for them.
        default: break;
    }
    return b;
}
```

- [ ] **Step 4: Run tests — verify pass**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[budget]" 2>&1 | tail -5
```
Expected: PASS.

- [ ] **Step 5: Run full suite + commit**

```bash
./build/botany_tests 2>&1 | tail -3
git add src/engine/vascular_sub_stepped.h src/engine/vascular_sub_stepped.cpp tests/test_vascular_sub_stepped.cpp
git commit -m "feat(vascular): compute_budget for per-node supply/demand snapshot

Classifies each node as source/sink/neither for each chemical this tick.
Leaves are sugar sources and water sinks.  Active meristems are sugar
and water sinks.  Root apicals push cytokinin into parent xylem.

Budgets are computed once at the start of vascular_sub_stepped so that
amortization across N sub-steps is pure division.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md"
```

---

### Task 14: Inject step (sources push into parent conduit)

**Files:**
- Modify: `src/engine/vascular_sub_stepped.cpp`
- Modify: `tests/test_vascular_sub_stepped.cpp`

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("inject step: leaf pushes sugar into parent stem phloem", "[vascular_sub_stepped][inject]") {
    Genome g = default_genome();
    WorldParams world;
    StemNode stem(1, glm::vec3(0), 0.02f);
    LeafNode leaf(2, glm::vec3(0.05f, 0, 0), 0.01f);
    stem.add_child(&leaf);
    leaf.parent = &stem;

    leaf.local().chemical(ChemicalID::Sugar) = 10.0f;
    VascularBudget b = compute_budget(leaf, g, world);
    float budget_slice = b.sugar_supply / 10.0f;  // N = 10 sub-steps

    float leaf_before   = leaf.local().chemical(ChemicalID::Sugar);
    float phloem_before = stem.phloem()->chemical(ChemicalID::Sugar);

    inject_step(leaf, b, /* N */ 10, g);

    float leaf_after   = leaf.local().chemical(ChemicalID::Sugar);
    float phloem_after = stem.phloem()->chemical(ChemicalID::Sugar);

    REQUIRE(leaf_after   == Catch::Approx(leaf_before   - budget_slice).margin(1e-5f));
    REQUIRE(phloem_after == Catch::Approx(phloem_before + budget_slice).margin(1e-5f));
}
```

- [ ] **Step 2: Run — verify fail**

Expected: compile error, `inject_step` undefined.

- [ ] **Step 3: Implement**

In `vascular_sub_stepped.h` declare:

```cpp
void inject_step(Node& source, const VascularBudget& b, uint32_t N, const Genome& g);
```

In `vascular_sub_stepped.cpp`:

```cpp
void inject_step(Node& source, const VascularBudget& b, uint32_t N, const Genome& /*g*/) {
    if (N == 0) return;

    // Sugar injection: leaves push into parent phloem (active pump, uphill allowed).
    if (b.sugar_supply > 0.0f) {
        if (auto* target_pool = source.nearest_phloem_upstream()) {
            const float slice = b.sugar_supply / static_cast<float>(N);
            const float actual = std::min(slice, source.local().chemical(ChemicalID::Sugar));
            source.local().chemical(ChemicalID::Sugar)   -= actual;
            target_pool->chemical(ChemicalID::Sugar)     += actual;
        }
    }

    // Cytokinin injection: root apicals push into parent xylem.
    if (b.cytokinin_supply > 0.0f) {
        if (auto* target_pool = source.nearest_xylem_upstream()) {
            const float slice = b.cytokinin_supply / static_cast<float>(N);
            const float actual = std::min(slice, source.local().chemical(ChemicalID::Cytokinin));
            source.local().chemical(ChemicalID::Cytokinin) -= actual;
            target_pool->chemical(ChemicalID::Cytokinin)   += actual;
        }
    }
}
```

- [ ] **Step 4: Run tests — verify pass**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[inject]" 2>&1 | tail -5
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
./build/botany_tests 2>&1 | tail -3
git add src/engine/vascular_sub_stepped.h src/engine/vascular_sub_stepped.cpp tests/test_vascular_sub_stepped.cpp
git commit -m "feat(vascular): inject_step — sources push budget/N into parent conduit

Leaves push sugar into parent.phloem (active pump).  Root apicals push
cytokinin into parent.xylem (gradient-based but capped by budget slice).
Capped by actual local availability so we never go negative.

Walk-up via nearest_phloem_upstream()/nearest_xylem_upstream() handles
the case where the direct parent is another specialty node.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md"
```

---

### Task 15: Extract step (sinks pull from parent conduit)

**Files:**
- Modify: `src/engine/vascular_sub_stepped.cpp`
- Modify: `tests/test_vascular_sub_stepped.cpp`

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("extract step: meristem pulls sugar from parent phloem", "[vascular_sub_stepped][extract]") {
    Genome g = default_genome();
    WorldParams world;
    StemNode stem(1, glm::vec3(0), 0.02f);
    ApicalNode apical(2, glm::vec3(0.05f, 0, 0), 0.01f);
    stem.add_child(&apical);
    apical.parent = &stem;
    apical.active = true;

    stem.phloem()->chemical(ChemicalID::Sugar) = 5.0f;
    apical.local().chemical(ChemicalID::Sugar) = 0.0f;

    VascularBudget b = compute_budget(apical, g, world);
    REQUIRE(b.sugar_demand > 0.0f);
    float slice = b.sugar_demand / 10.0f;

    extract_step(apical, b, /* N */ 10, g);

    // Sugar moved from phloem to meristem local, up to slice.
    REQUIRE(apical.local().chemical(ChemicalID::Sugar) == Catch::Approx(slice).margin(1e-5f));
    REQUIRE(stem.phloem()->chemical(ChemicalID::Sugar) == Catch::Approx(5.0f - slice).margin(1e-5f));
}

TEST_CASE("extract step capped by available phloem sugar", "[vascular_sub_stepped][extract]") {
    Genome g = default_genome();
    WorldParams world;
    StemNode stem(1, glm::vec3(0), 0.02f);
    ApicalNode apical(2, glm::vec3(0.05f, 0, 0), 0.01f);
    stem.add_child(&apical);
    apical.parent = &stem;
    apical.active = true;

    stem.phloem()->chemical(ChemicalID::Sugar) = 0.01f;  // almost dry
    apical.local().chemical(ChemicalID::Sugar) = 0.0f;

    VascularBudget b = compute_budget(apical, g, world);
    extract_step(apical, b, /* N */ 10, g);

    // Capped by phloem availability — no negative phloem.
    REQUIRE(stem.phloem()->chemical(ChemicalID::Sugar) >= 0.0f);
    REQUIRE(apical.local().chemical(ChemicalID::Sugar) <= 0.01f);
}
```

- [ ] **Step 2: Implement in `vascular_sub_stepped.cpp`**

Declaration in header:
```cpp
void extract_step(Node& sink, const VascularBudget& b, uint32_t N, const Genome& g);
```

Implementation:
```cpp
void extract_step(Node& sink, const VascularBudget& b, uint32_t N, const Genome& /*g*/) {
    if (N == 0) return;

    if (b.sugar_demand > 0.0f) {
        if (auto* source_pool = sink.nearest_phloem_upstream()) {
            const float slice = b.sugar_demand / static_cast<float>(N);
            const float actual = std::min(slice, source_pool->chemical(ChemicalID::Sugar));
            source_pool->chemical(ChemicalID::Sugar)   -= actual;
            sink.local().chemical(ChemicalID::Sugar)   += actual;
        }
    }

    if (b.water_demand > 0.0f) {
        if (auto* source_pool = sink.nearest_xylem_upstream()) {
            const float slice = b.water_demand / static_cast<float>(N);
            const float actual = std::min(slice, source_pool->chemical(ChemicalID::Water));
            source_pool->chemical(ChemicalID::Water)   -= actual;
            sink.local().chemical(ChemicalID::Water)   += actual;
            // Cytokinin rides along passively — it's carried in xylem water.
            // Move cytokinin proportional to water drawn from this pool.
            const float water_after = source_pool->chemical(ChemicalID::Water);
            if (water_after + actual > 1e-8f) {  // avoid div-by-zero
                const float cyto_ratio = actual / (water_after + actual);
                const float cyto_move  = source_pool->chemical(ChemicalID::Cytokinin) * cyto_ratio;
                source_pool->chemical(ChemicalID::Cytokinin) -= cyto_move;
                sink.local().chemical(ChemicalID::Cytokinin) += cyto_move;
            }
        }
    }
}
```

- [ ] **Step 3: Run tests — verify pass**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[extract]" 2>&1 | tail -5
```
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
./build/botany_tests 2>&1 | tail -3
git add src/engine/vascular_sub_stepped.h src/engine/vascular_sub_stepped.cpp tests/test_vascular_sub_stepped.cpp
git commit -m "feat(vascular): extract_step — sinks pull budget/N from parent conduit

Capped by actual pool availability (never drives pool negative).  Water
extraction also carries cytokinin proportionally — cytokinin rides in
xylem solution, so it follows water mass.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md"
```

---

### Task 16: Radial flow step

**Files:**
- Modify: `src/engine/vascular_sub_stepped.cpp`
- Modify: `tests/test_vascular_sub_stepped.cpp`

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("radial flow equilibrates stem local and phloem", "[vascular_sub_stepped][radial_flow]") {
    Genome g = default_genome();
    StemNode stem(1, glm::vec3(0), 0.015f);
    stem.offset = glm::vec3(0.0f, 0.05f, 0.0f);

    // Start with high phloem, empty local.
    stem.phloem()->chemical(ChemicalID::Sugar) = 1.0f;
    stem.local().chemical(ChemicalID::Sugar)   = 0.0f;

    float total_before = stem.phloem()->chemical(ChemicalID::Sugar)
                       + stem.local().chemical(ChemicalID::Sugar);

    // Run radial flow many times — local should approach phloem concentration.
    for (int i = 0; i < 100; ++i) {
        radial_flow_step(stem, /* N */ 1, g);
    }

    float total_after = stem.phloem()->chemical(ChemicalID::Sugar)
                      + stem.local().chemical(ChemicalID::Sugar);

    // Conservation.
    REQUIRE(total_after == Catch::Approx(total_before).margin(1e-4f));
    // After many iterations, concentrations should be close (concentration =
    // chemical / capacity).  Both compartments share the same chemical here;
    // we're just checking that some flowed from phloem to local.
    REQUIRE(stem.local().chemical(ChemicalID::Sugar) > 0.0f);
    REQUIRE(stem.phloem()->chemical(ChemicalID::Sugar) < 1.0f);
}

TEST_CASE("radial permeability is smaller for thicker stems", "[vascular_sub_stepped][radial_flow]") {
    Genome g = default_genome();
    float perm_young = radial_permeability_sugar(0.01f, g);
    float perm_mature = radial_permeability_sugar(2.0f, g);
    REQUIRE(perm_mature < perm_young);
}
```

- [ ] **Step 2: Implement**

Declaration in header:
```cpp
void radial_flow_step(Node& conduit, uint32_t N, const Genome& g);
```

Implementation in `.cpp`. The formula moves chemical between `local` and `phloem`/`xylem` toward gradient equality, rate-limited by the radius-dependent permeability. Use concentrations (chemical / capacity) for the gradient to keep behavior sensible as pool sizes change.

```cpp
#include "engine/sugar.h"  // sugar_cap
// water_cap should already be imported via earlier tasks

void radial_flow_step(Node& n, uint32_t N, const Genome& g) {
    if (N == 0) return;
    const float inv_N = 1.0f / static_cast<float>(N);

    // --- Sugar: phloem <-> local ---
    if (auto* phl = n.phloem()) {
        const float perm     = radial_permeability_sugar(n.radius, g);
        const float cap_phl  = phloem_capacity(n, g);
        const float cap_loc  = sugar_cap(n, g);
        if (cap_phl > 1e-8f && cap_loc > 1e-8f) {
            const float conc_phl = phl->chemical(ChemicalID::Sugar) / cap_phl;
            const float conc_loc = n.local().chemical(ChemicalID::Sugar) / cap_loc;
            // Flow toward equal concentration; positive = phloem → local.
            const float dconc = conc_phl - conc_loc;
            float flow = perm * dconc * inv_N;
            // Cap the flow to avoid overshoot in one step.
            const float max_equalize_volume = std::min(cap_phl, cap_loc) * 0.5f;
            flow = std::clamp(flow, -max_equalize_volume, max_equalize_volume);
            // Don't draw more than what's present in either side.
            if (flow > 0.0f) flow = std::min(flow, phl->chemical(ChemicalID::Sugar));
            else             flow = -std::min(-flow, n.local().chemical(ChemicalID::Sugar));

            phl->chemical(ChemicalID::Sugar)       -= flow;
            n.local().chemical(ChemicalID::Sugar)  += flow;
        }
    }

    // --- Water + cytokinin: xylem <-> local ---
    if (auto* xyl = n.xylem()) {
        const float perm     = radial_permeability_water(n.radius, g);
        const float cap_xyl  = xylem_capacity(n, g);
        const float cap_loc  = water_cap(n, g);
        if (cap_xyl > 1e-8f && cap_loc > 1e-8f) {
            const float conc_xyl = xyl->chemical(ChemicalID::Water) / cap_xyl;
            const float conc_loc = n.local().chemical(ChemicalID::Water) / cap_loc;
            const float dconc = conc_xyl - conc_loc;
            float flow = perm * dconc * inv_N;
            const float max_equalize_volume = std::min(cap_xyl, cap_loc) * 0.5f;
            flow = std::clamp(flow, -max_equalize_volume, max_equalize_volume);
            if (flow > 0.0f) flow = std::min(flow, xyl->chemical(ChemicalID::Water));
            else             flow = -std::min(-flow, n.local().chemical(ChemicalID::Water));

            xyl->chemical(ChemicalID::Water)       -= flow;
            n.local().chemical(ChemicalID::Water)  += flow;

            // Cytokinin rides with water proportionally.
            const float water_before = xyl->chemical(ChemicalID::Water) + flow;
            if (water_before > 1e-8f && flow != 0.0f) {
                const float cyto_ratio = flow / water_before;
                const float cyto_flow  = xyl->chemical(ChemicalID::Cytokinin) * cyto_ratio;
                xyl->chemical(ChemicalID::Cytokinin)      -= cyto_flow;
                n.local().chemical(ChemicalID::Cytokinin) += cyto_flow;
            }
        }
    }
}
```

- [ ] **Step 3: Run tests — verify pass**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[radial_flow]" 2>&1 | tail -5
```
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
./build/botany_tests 2>&1 | tail -3
git add src/engine/vascular_sub_stepped.h src/engine/vascular_sub_stepped.cpp tests/test_vascular_sub_stepped.cpp
git commit -m "feat(vascular): radial_flow_step — radius-dependent local↔conduit exchange

Bidirectional gradient flow between a stem/root's own local_env and own
phloem/xylem.  Rate scales with radial_permeability(radius) per spec
section 6 — young stems leak freely, mature trunks approach the floor.

Mass-conservative by construction (paired ±flow).  Water radial flow
carries cytokinin proportionally.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md"
```

---

### Task 17: Longitudinal Jacobi step

**Files:**
- Modify: `src/engine/vascular_sub_stepped.cpp`
- Modify: `tests/test_vascular_sub_stepped.cpp`

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("Jacobi equalizes two connected phloem pools", "[vascular_sub_stepped][jacobi]") {
    Genome g = default_genome();
    StemNode parent_stem(1, glm::vec3(0), 0.02f);
    parent_stem.offset = glm::vec3(0.0f, 0.05f, 0.0f);
    StemNode child_stem(2, glm::vec3(0, 0.05f, 0), 0.02f);
    child_stem.offset = glm::vec3(0.0f, 0.05f, 0.0f);
    parent_stem.add_child(&child_stem);
    child_stem.parent = &parent_stem;

    parent_stem.phloem()->chemical(ChemicalID::Sugar) = 1.0f;
    child_stem.phloem()->chemical(ChemicalID::Sugar)  = 0.0f;

    float total_before = parent_stem.phloem()->chemical(ChemicalID::Sugar)
                       + child_stem.phloem()->chemical(ChemicalID::Sugar);

    // Many Jacobi iterations should bring pressures close together.
    for (int i = 0; i < 100; ++i) {
        jacobi_step(parent_stem, child_stem, g);
    }

    float total_after = parent_stem.phloem()->chemical(ChemicalID::Sugar)
                      + child_stem.phloem()->chemical(ChemicalID::Sugar);

    REQUIRE(total_after == Catch::Approx(total_before).margin(1e-4f));
    // Some sugar moved from parent to child.
    REQUIRE(child_stem.phloem()->chemical(ChemicalID::Sugar) > 0.0f);
    REQUIRE(parent_stem.phloem()->chemical(ChemicalID::Sugar) < 1.0f);
}
```

- [ ] **Step 2: Implement**

Declaration in header:
```cpp
// One Jacobi pass between two adjacent conduit nodes.  Equalizes both
// phloem (sugar) and xylem (water + cytokinin) pressures between them.
void jacobi_step(Node& parent, Node& child, const Genome& g);
```

Implementation:
```cpp
void jacobi_step(Node& parent, Node& child, const Genome& g) {
    // Phloem sugar.
    if (auto* pp = parent.phloem()) {
        if (auto* cp = child.phloem()) {
            const float cap_p = phloem_capacity(parent, g);
            const float cap_c = phloem_capacity(child,  g);
            if (cap_p > 1e-8f && cap_c > 1e-8f) {
                const float pressure_p = pp->chemical(ChemicalID::Sugar) / cap_p;
                const float pressure_c = cp->chemical(ChemicalID::Sugar) / cap_c;
                // Conductance: use smaller-pipe cross-section × canalization.
                const float base_cond = std::min(cap_p, cap_c);
                const float bias      = parent.get_bias_multiplier(&child, g);
                const float conductance = base_cond * bias;
                // Flow, clamped to half the pressure difference × smaller cap
                // to avoid overshoot.
                float flow = conductance * (pressure_p - pressure_c);
                const float max_move = 0.5f * std::min(cap_p, cap_c);
                flow = std::clamp(flow, -max_move, max_move);
                // Don't drive either pool negative.
                if (flow > 0.0f) flow = std::min(flow, pp->chemical(ChemicalID::Sugar));
                else             flow = -std::min(-flow, cp->chemical(ChemicalID::Sugar));

                pp->chemical(ChemicalID::Sugar) -= flow;
                cp->chemical(ChemicalID::Sugar) += flow;
            }
        }
    }

    // Xylem water + cytokinin.
    if (auto* pp = parent.xylem()) {
        if (auto* cp = child.xylem()) {
            const float cap_p = xylem_capacity(parent, g);
            const float cap_c = xylem_capacity(child,  g);
            if (cap_p > 1e-8f && cap_c > 1e-8f) {
                const float pressure_p = pp->chemical(ChemicalID::Water) / cap_p;
                const float pressure_c = cp->chemical(ChemicalID::Water) / cap_c;
                const float base_cond  = std::min(cap_p, cap_c);
                const float bias       = parent.get_bias_multiplier(&child, g);
                const float conductance = base_cond * bias;
                float flow = conductance * (pressure_p - pressure_c);
                const float max_move = 0.5f * std::min(cap_p, cap_c);
                flow = std::clamp(flow, -max_move, max_move);
                if (flow > 0.0f) flow = std::min(flow, pp->chemical(ChemicalID::Water));
                else             flow = -std::min(-flow, cp->chemical(ChemicalID::Water));

                // Water move.
                pp->chemical(ChemicalID::Water) -= flow;
                cp->chemical(ChemicalID::Water) += flow;

                // Cytokinin rides with water.
                const float water_source_before = (flow > 0.0f)
                    ? pp->chemical(ChemicalID::Water) + flow
                    : cp->chemical(ChemicalID::Water) - flow;
                if (water_source_before > 1e-8f && std::abs(flow) > 1e-8f) {
                    TransportPool* src = (flow > 0.0f) ? pp : cp;
                    TransportPool* dst = (flow > 0.0f) ? cp : pp;
                    const float cyto_ratio = std::abs(flow) / water_source_before;
                    const float cyto_move  = src->chemical(ChemicalID::Cytokinin) * cyto_ratio;
                    src->chemical(ChemicalID::Cytokinin) -= cyto_move;
                    dst->chemical(ChemicalID::Cytokinin) += cyto_move;
                }
            }
        }
    }
}
```

- [ ] **Step 3: Run tests — verify pass**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[jacobi]" 2>&1 | tail -5
```
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
./build/botany_tests 2>&1 | tail -3
git add src/engine/vascular_sub_stepped.h src/engine/vascular_sub_stepped.cpp tests/test_vascular_sub_stepped.cpp
git commit -m "feat(vascular): jacobi_step — neighbor pressure equalization

Pure neighbor equalizer: looks at two adjacent conduits' pressures
(chemical / capacity) and moves chemical to reduce the gradient.  Has no
awareness of sources or sinks — pressure fields created by inject/extract
at node locations drive routing automatically.

Conductance scales with smaller pipe's capacity × canalization bias.
Clamped to half-equilibrium volume to avoid overshoot.  Water Jacobi
carries cytokinin proportionally.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md"
```

---

### Task 18: Wire up the full sub-step loop in `vascular_sub_stepped()`

**Files:**
- Modify: `src/engine/vascular_sub_stepped.cpp`
- Modify: `tests/test_vascular_sub_stepped.cpp`

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("full vascular_sub_stepped delivers sugar to apex within N hops", "[vascular_sub_stepped][integration]") {
    Genome g = default_genome();
    WorldParams world;
    world.vascular_substeps = 20;
    Plant plant(g, glm::vec3(0.0f));

    // Build a 10-stem chain above the seed.
    Node* tip_stem = plant.seed_mut();
    for (int i = 0; i < 10; ++i) {
        Node* stem = plant.create_node(NodeType::STEM,
            glm::vec3(0.0f, 0.05f * (i + 1), 0.0f), 0.015f);
        tip_stem->add_child(stem);
        tip_stem = stem;
    }

    // Zero all sugar; pile 10 g at the seed.
    plant.for_each_node_mut([&](Node& n) {
        n.local().chemical(ChemicalID::Sugar) = 0.0f;
        if (auto* p = n.phloem()) p->chemical(ChemicalID::Sugar) = 0.0f;
    });
    // Seed is a STEM with phloem.  Put sugar directly into the phloem so
    // Jacobi can propagate it without waiting on leaf loading.
    plant.seed_mut()->phloem()->chemical(ChemicalID::Sugar) = 10.0f;

    vascular_sub_stepped(plant, g, world);

    // The tip stem is 10 hops from seed.  N=20 sub-steps means 20 Jacobi
    // iterations, so pressure wave reaches the tip.
    float tip_phloem = tip_stem->phloem()->chemical(ChemicalID::Sugar);
    REQUIRE(tip_phloem > 0.001f);
}

TEST_CASE("vascular_sub_stepped conserves mass", "[vascular_sub_stepped][integration][conservation]") {
    Genome g = default_genome();
    WorldParams world;
    world.vascular_substeps = 20;
    Plant plant(g, glm::vec3(0.0f));

    // Build a tiny plant: seed + one stem + one leaf.
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0, 0.05f, 0), 0.015f);
    plant.seed_mut()->add_child(stem);
    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.02f, 0.05f, 0), 0.01f);
    stem->add_child(leaf);

    // Seed sugar and water in various compartments.
    plant.seed_mut()->phloem()->chemical(ChemicalID::Sugar) = 1.0f;
    stem->phloem()->chemical(ChemicalID::Sugar)             = 0.5f;
    leaf->local().chemical(ChemicalID::Sugar)               = 2.0f;
    plant.seed_mut()->xylem()->chemical(ChemicalID::Water)  = 1.0f;
    stem->xylem()->chemical(ChemicalID::Water)              = 0.5f;

    // Pre-sum across all pools.
    float sugar_before = 0.0f, water_before = 0.0f;
    plant.for_each_node([&](const Node& n) {
        sugar_before += n.local().chemical(ChemicalID::Sugar);
        water_before += n.local().chemical(ChemicalID::Water);
        if (auto* p = n.phloem()) sugar_before += p->chemical(ChemicalID::Sugar);
        if (auto* x = n.xylem())  water_before += x->chemical(ChemicalID::Water);
    });

    vascular_sub_stepped(plant, g, world);

    float sugar_after = 0.0f, water_after = 0.0f;
    plant.for_each_node([&](const Node& n) {
        sugar_after += n.local().chemical(ChemicalID::Sugar);
        water_after += n.local().chemical(ChemicalID::Water);
        if (auto* p = n.phloem()) sugar_after += p->chemical(ChemicalID::Sugar);
        if (auto* x = n.xylem())  water_after += x->chemical(ChemicalID::Water);
    });

    REQUIRE(sugar_after == Catch::Approx(sugar_before).margin(1e-4f));
    REQUIRE(water_after == Catch::Approx(water_before).margin(1e-4f));
}
```

- [ ] **Step 2: Implement the full loop**

Replace the stub body of `vascular_sub_stepped()` in `.cpp`:

```cpp
void vascular_sub_stepped(Plant& plant, const Genome& g, const WorldParams& world) {
    const uint32_t N = std::max<uint32_t>(1, world.vascular_substeps);

    FlatNodes flat = flatten(plant);

    // --- Part A: Budget snapshot ---
    // One budget per node.  Indexed by position in flat.all so we can reuse.
    std::vector<VascularBudget> budgets;
    budgets.reserve(flat.all.size());
    for (Node* n : flat.all) {
        budgets.push_back(compute_budget(*n, g, world));
    }

    // --- Part B: Sub-step loop ---
    for (uint32_t iter = 0; iter < N; ++iter) {
        // Step 1: Inject at sources.
        for (size_t i = 0; i < flat.all.size(); ++i) {
            const VascularBudget& b = budgets[i];
            if (b.sugar_supply > 0.0f || b.cytokinin_supply > 0.0f) {
                inject_step(*flat.all[i], b, N, g);
            }
        }

        // Step 2: Radial flow on every conduit.
        for (Node* n : flat.conduits) {
            radial_flow_step(*n, N, g);
        }

        // Step 3: Extract at sinks.
        for (size_t i = 0; i < flat.all.size(); ++i) {
            const VascularBudget& b = budgets[i];
            if (b.sugar_demand > 0.0f || b.water_demand > 0.0f) {
                extract_step(*flat.all[i], b, N, g);
            }
        }

        // Step 4: Longitudinal Jacobi across every conduit edge.
        for (Node* n : flat.conduits) {
            for (Node* child : n->children) {
                if (child->phloem() || child->xylem()) {
                    jacobi_step(*n, *child, g);
                }
            }
        }
    }
}
```

- [ ] **Step 3: Run tests — verify pass**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[vascular_sub_stepped][integration]" 2>&1 | tail -5
```
Expected: PASS.

- [ ] **Step 4: Run full suite**

```bash
./build/botany_tests 2>&1 | tail -3
```
Expected: `All tests passed`. The old vascular_transport is still running; our new function isn't wired in yet.

- [ ] **Step 5: Commit**

```bash
git add src/engine/vascular_sub_stepped.cpp tests/test_vascular_sub_stepped.cpp
git commit -m "feat(vascular): full vascular_sub_stepped loop

Integrates the helpers built in Tasks 10-17 into the full algorithm:
N iterations of [inject sources → radial flow → extract sinks → Jacobi
neighbors].  Budgets snapshot once at start, then amortized across N.

Integration tests verify:
  - Sugar propagates from seed phloem to tip of 10-stem chain within one
    tick when N >= chain length.
  - Mass conservation: total sugar and water unchanged across the
    vascular pass.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md"
```

---

### Task 19: Compartment invariant and distance-dependent supply tests

**Files:**
- Modify: `tests/test_vascular_sub_stepped.cpp`

- [ ] **Step 1: Write the invariant tests**

```cpp
TEST_CASE("compartment residency invariants hold after vascular", "[vascular_sub_stepped][invariants]") {
    Genome g = default_genome();
    WorldParams world;
    Plant plant(g, glm::vec3(0.0f));

    // Populate with a typical tick's chemicals.
    plant.for_each_node_mut([&](Node& n) {
        n.local().chemical(ChemicalID::Sugar)     = 1.0f;
        n.local().chemical(ChemicalID::Water)     = 0.5f;
        n.local().chemical(ChemicalID::Auxin)     = 0.1f;
        if (auto* p = n.phloem()) p->chemical(ChemicalID::Sugar) = 0.3f;
        if (auto* x = n.xylem())  x->chemical(ChemicalID::Water) = 0.2f;
    });

    vascular_sub_stepped(plant, g, world);

    // Sugar must never appear in any xylem pool.
    // Water/cytokinin must never appear in any phloem pool.
    // Auxin/gibberellin/stress must never appear in any transport pool.
    plant.for_each_node([&](const Node& n) {
        if (auto* p = n.phloem()) {
            REQUIRE(p->chemical(ChemicalID::Water)     == 0.0f);
            REQUIRE(p->chemical(ChemicalID::Cytokinin) == 0.0f);
            REQUIRE(p->chemical(ChemicalID::Auxin)     == 0.0f);
        }
        if (auto* x = n.xylem()) {
            REQUIRE(x->chemical(ChemicalID::Sugar)     == 0.0f);
            REQUIRE(x->chemical(ChemicalID::Auxin)     == 0.0f);
        }
    });
}

TEST_CASE("distal apex under-supplied when chain exceeds N", "[vascular_sub_stepped][hydraulic_limit]") {
    Genome g = default_genome();
    WorldParams world;
    world.vascular_substeps = 5;  // deliberately small
    Plant plant(g, glm::vec3(0.0f));

    // Build a 15-stem chain — longer than N.
    Node* tip = plant.seed_mut();
    for (int i = 0; i < 15; ++i) {
        Node* stem = plant.create_node(NodeType::STEM,
            glm::vec3(0.0f, 0.05f * (i + 1), 0.0f), 0.015f);
        tip->add_child(stem);
        tip = stem;
    }

    // Pressure at seed only.
    plant.for_each_node_mut([&](Node& n) {
        if (auto* p = n.phloem()) p->chemical(ChemicalID::Sugar) = 0.0f;
    });
    plant.seed_mut()->phloem()->chemical(ChemicalID::Sugar) = 10.0f;

    vascular_sub_stepped(plant, g, world);

    // With N=5 and 15-stem chain, sugar should fall off with distance.  Check
    // that the tip phloem has less than the stem 5 hops from seed.
    Node* mid_stem = plant.seed_mut()->children[0];
    for (int i = 0; i < 4; ++i) mid_stem = mid_stem->children[0];  // 5 hops
    REQUIRE(tip->phloem()->chemical(ChemicalID::Sugar)
          < mid_stem->phloem()->chemical(ChemicalID::Sugar));
}
```

- [ ] **Step 2: Run tests — verify pass**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests "[invariants],[hydraulic_limit]" 2>&1 | tail -5
```
Expected: PASS.

- [ ] **Step 3: Run full suite + commit**

```bash
./build/botany_tests 2>&1 | tail -3
git add tests/test_vascular_sub_stepped.cpp
git commit -m "test(vascular): compartment invariants + hydraulic-limit regression

- Chemicals stay in their allowed compartments after a full vascular
  pass (no water in phloem, no sugar in xylem, etc).
- With N smaller than chain length, distal apex receives less sugar
  than mid-chain — confirms hydraulic limitation behaves as designed.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md"
```

---

## Phase E — Cutover

### Task 20: Switch `Plant::tick()` to new algorithm + tick-then-vascular ordering

**Files:**
- Modify: `src/engine/plant.cpp`
- Modify: possibly several existing test files that asserted old-algorithm behavior

**Rationale:** This is the behavior-changing commit. The metabolic tick now runs before vascular (Ordering B from the spec), and `vascular_sub_stepped` replaces `vascular_transport`. Existing tests that asserted pairwise-Jacobi-specific end states will need updates.

- [ ] **Step 1: Modify `Plant::tick()`**

Find the current tick body. It probably looks like:

```cpp
void Plant::tick(const WorldParams& world) {
    vascular_transport(*this, genome_, world);   // OLD
    tick_tree(seed_, world);
    // ...
}
```

Change to:

```cpp
#include "engine/vascular_sub_stepped.h"

void Plant::tick(const WorldParams& world) {
    tick_tree(seed_, world);                      // Phase 1: per-node metabolism
    vascular_sub_stepped(*this, genome_, world);  // Phase 2: vascular transport
    // ... any remaining post-tick maintenance ...
}
```

Remove the `#include "engine/vascular.h"` line if it's only present for the `vascular_transport` call (if other things in `plant.cpp` still use functions from vascular.h, leave the include).

- [ ] **Step 2: Build — expect some test failures**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests 2>&1 | tail -10
```

Failure patterns to expect:
- Tests in `test_vascularization.cpp` that assert pairwise-Jacobi behavior (canalization ratchet, conductance-weighted distribution). These need to be restated in terms of the new algorithm.
- Tests that assumed vascular runs before tick (i.e., that rely on sugar being delivered to sinks *before* the node tick in the same tick). Under Ordering B this is now a 1-tick lag; tests need to iterate one extra tick.
- Sugar-economy tests with specific equilibrium-state assertions may fail due to different tuning dynamics.

- [ ] **Step 3: Fix failing tests one at a time**

For each failing test:
1. Read its assertion carefully.
2. Is it testing a mass-conservation / plant-survival property (still valid, fix the test)? Or is it testing a mechanism-specific outcome (pairwise Jacobi edge behavior — needs restating)?
3. For mass-conservation / survival tests: usually just add `plant.tick(world)` once more so the new 1-tick lag is accommodated.
4. For mechanism-specific tests: either restate the assertion in terms of the new algorithm's equivalent invariant, or delete the test if it was testing a property that no longer exists.

Keep a running list of which tests were changed and why, in the commit message.

- [ ] **Step 4: Run full suite iteratively until green**

```bash
./build/botany_tests 2>&1 | grep -E "FAILED|failed" | head -20
```

Expected: incrementally fewer failures as you fix them. Final state: all green.

- [ ] **Step 5: Run the realtime viewer for smoke test**

```bash
./build/botany_realtime &
# Let it run for ~30 seconds, watch for:
#   - Plant grows continuously (not dying)
#   - No obvious visual glitches
#   - Realtime perf acceptable (no stutter)
# Kill with Ctrl+C.
```

- [ ] **Step 6: Commit**

```bash
git add src/engine/plant.cpp tests/
git commit -m "$(cat <<'EOF'
refactor(vascular): cutover to vascular_sub_stepped + tick-then-vascular

BEHAVIOR CHANGE.  This is the Phase E cutover in the compartmented
vascular refactor:

  - Plant::tick() now runs tick_tree() first, then vascular_sub_stepped().
    Previously ran vascular_transport() first.
  - The old pairwise-Jacobi vascular_transport() is no longer called from
    the main tick path.  (File is not yet deleted — Phase F.)

Tests updated:
  [list the specific tests that needed changes, with one-line reasons]

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md
EOF
)"
```

---

## Phase F — Dead code removal

### Task 21: Delete old vascular code

**Files:**
- Modify: `src/engine/vascular.cpp` — delete `has_vasculature`, `phloem_resolve`, `xylem_resolve`, `vascular_transport` and any private helpers used only by them
- Modify: `src/engine/vascular.h` — delete the corresponding declarations
- Modify: Any caller of `has_vasculature()` that isn't Plant::tick (there's one in `node.cpp` around line 382 — `transport_with_children`)

**Rationale:** The old algorithm is dead code after the Phase E cutover. Delete it so future readers don't confuse it with the new one.

- [ ] **Step 1: Check for remaining callers of old functions**

```bash
cd /Users/wldarden/learning/botany
grep -rn "has_vasculature\|phloem_resolve\|xylem_resolve\|vascular_transport" src/ --include="*.cpp" --include="*.h"
```

Expected: matches only in the files we're about to edit (`vascular.h`, `vascular.cpp`), and possibly in `node.cpp`'s `transport_with_children` function. Any other caller needs handling.

- [ ] **Step 2: Remove `has_vasculature()` gate from `transport_with_children()`**

In `src/engine/node/node.cpp` around line 380, the block that skips vascular chemicals on mature-to-mature edges:

```cpp
bool parent_vascular = has_vasculature(*this, g);
// ...
if (is_vascular_chemical(dp.id) && parent_vascular && has_vasculature(*child, g)) {
    continue;
}
```

Under the new model, `transport_with_children` only handles non-vascular chemicals (auxin, gibberellin, stress). Replace the skip logic with a direct type-based skip:

```cpp
// Vascular chemicals are now handled entirely by vascular_sub_stepped();
// per-node local diffusion only runs for signaling chemicals.
if (is_vascular_chemical(dp.id)) continue;
```

Remove the `parent_vascular` variable and both `has_vasculature()` calls.

- [ ] **Step 3: Delete `vascular.cpp` contents**

Open `src/engine/vascular.cpp`. Delete:
- The `has_vasculature` function
- `phloem_resolve`
- `xylem_resolve`
- `vascular_transport`
- All private static helpers that are no longer used: `pipe_capacity`, `build_flat`, `chem_name`, `node_type_name`, `phloem_ring_area`, `edge_pipe_radius`, `node_volume`, `compute_phloem_pressure`, `unloading_permeability`, etc.
- The debug logging functions if they only serve the old algorithm

If `vascular.cpp` ends up empty, delete the file entirely and remove it from `CMakeLists.txt`.

- [ ] **Step 4: Clean up `vascular.h`**

Delete the corresponding declarations. If `vascular.h` ends up empty, delete it and remove any `#include "engine/vascular.h"` lines from other files.

- [ ] **Step 5: Build**

```bash
/usr/local/bin/cmake --build build 2>&1 | grep "error:" | head -10
```

Fix any errors (missing includes, etc) in place.

- [ ] **Step 6: Run full suite**

```bash
./build/botany_tests 2>&1 | tail -3
```
Expected: `All tests passed`.

- [ ] **Step 7: Commit**

```bash
git add src/engine/vascular.cpp src/engine/vascular.h src/engine/node/node.cpp CMakeLists.txt
git commit -m "$(cat <<'EOF'
refactor(vascular): delete old pairwise-Jacobi vascular code

Phase F of the compartmented vascular refactor.  vascular_sub_stepped
has been in production since Phase E cutover; the old vascular_transport
+ phloem_resolve + xylem_resolve + has_vasculature are dead code.

Also simplified transport_with_children() in node.cpp — it no longer
needs to gate on has_vasculature() because vascular chemicals (sugar,
water, cytokinin) are now handled exclusively by vascular_sub_stepped,
and transport_with_children only runs for signaling chemicals.

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md
EOF
)"
```

---

### Task 22: Remove unused genome/world params

**Files:**
- Modify: `src/engine/genome.h` — remove `vascular_radius_threshold` (no longer used after has_vasculature is gone)
- Modify: `src/engine/world_params.h` — remove `phloem_iterations` (was already deprecated; the old vascular_transport was its only user)
- Modify: `src/evolution/genome_bridge.cpp` — remove the corresponding gene entries

- [ ] **Step 1: Confirm no remaining usage**

```bash
grep -rn "vascular_radius_threshold\|phloem_iterations" src/ tests/ --include="*.cpp" --include="*.h"
```
Expected: only in genome.h / world_params.h / genome_bridge.cpp (the declarations/registrations). If other references exist, investigate — they may be obsolete.

- [ ] **Step 2: Delete the fields**

In `genome.h`, remove the `vascular_radius_threshold` line and its default in `default_genome()`.
In `world_params.h`, remove the `phloem_iterations` line.
In `genome_bridge.cpp`, remove the corresponding gene entries.

- [ ] **Step 3: Build and run full suite**

```bash
/usr/local/bin/cmake --build build 2>&1 | tail -3 && ./build/botany_tests 2>&1 | tail -3
```
Expected: build clean, `All tests passed`.

- [ ] **Step 4: Commit**

```bash
git add src/engine/genome.h src/engine/world_params.h src/evolution/genome_bridge.cpp
git commit -m "$(cat <<'EOF'
refactor(vascular): remove obsolete vascular_radius_threshold and phloem_iterations

vascular_radius_threshold was the gate for has_vasculature(), which was
deleted in Task 21.  phloem_iterations was the iteration count for the
old pairwise-Jacobi, also gone.  Serialized genomes that contain these
fields will ignore them (genome_bridge no longer registers them).

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md
EOF
)"
```

---

## Phase G — Documentation

### Task 23: Rewrite CLAUDE.md transport and tick sections

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Update the Chemical Transport Model section**

Find the `## Chemical Transport Model` section in CLAUDE.md. Replace its contents to describe the new model:

```markdown
## Chemical Transport Model

The plant uses a **compartmented dual transport system** matching real plant biology:

**Three compartments per node class:**
- `local_env` — every node owns one. Holds the node's own parenchyma chemicals: sugar, water, auxin, cytokinin, gibberellin, stress (anything the node itself uses for metabolism, growth, or signaling).
- `phloem` — StemNode and RootNode only. Sieve-tube content, carries sugar longitudinally.
- `xylem` — StemNode and RootNode only. Vessel content, carries water and cytokinin longitudinally.

Specialty nodes (LeafNode, ApicalNode, RootApicalNode) have only `local_env`. They interact with their nearest ancestor's phloem/xylem via `nearest_phloem_upstream()` and `nearest_xylem_upstream()` walk-up helpers.

All chemical access goes through explicit compartment accessors:
- `node.local().chemical(id)` for the local compartment
- `stem.phloem()->chemical(id)` and `stem.xylem()->chemical(id)` for transport pools

**Two transport pathways:**

1. **Sub-stepped vascular bulk flow** (`vascular_sub_stepped.cpp`) — N iterations per tick of:
   1. Inject — sources (leaves push sugar, root apicals push cytokinin) transfer `budget/N` into their walk-up parent's conduit
   2. Radial flow — stem/root local_env ⇄ own phloem/xylem, rate-limited by radius-dependent `radial_permeability(r)` (see Radial Permeability section)
   3. Extract — sinks (meristems pull sugar, leaves/meristems pull water) transfer `budget/N` from their walk-up parent's conduit
   4. Longitudinal Jacobi — one pass of neighbor pressure equalization across every (parent, child) edge that has matching conduit pools

   Jacobi is a pure neighbor equalizer — no awareness of sources or sinks. Pressure fields created by local injection and extraction drive routing automatically.

2. **Local diffusion** (`Node::transport_with_children()`) — per-node during DFS walk, handles only signaling chemicals: auxin, gibberellin, stress. Short-range cell-to-cell transport. Uses `transport_received` buffer for anti-teleportation.

| Chemical | Pathway | Notes |
|----------|---------|-------|
| Auxin | Local diffusion | Polar cell-to-cell, basipetal bias |
| Gibberellin | Local diffusion | Short-range, local to producing leaf |
| Stress | Local diffusion | Local mechanical alarm signal |
| Ethylene | Spatial gas diffusion | Global 3D pass (unchanged) |
| Sugar | Vascular (phloem) + radial | Leaves (source) → phloem → radial into stem local, direct pull at meristems |
| Water | Vascular (xylem) + radial | Roots (source, via radial) → xylem → direct pull at leaves/meristems |
| Cytokinin | Vascular (xylem) | Rides the xylem stream proportionally with water |

**Budgets are frozen at the start of each vascular pass.** `compute_budget(node, g, world)` classifies each node as source/sink for each chemical using the state at the moment vascular starts. Amortization across N sub-steps is then pure division — no dynamic re-evaluation within the loop.

**Fixed N (`world.vascular_substeps`).** Each Jacobi iteration propagates pressure by ~1 hop. Plants whose longest source-to-sink chain exceeds N will show distance-dependent apical supply — distal apices get less sugar than proximal ones. This is intentional and biologically realistic (real tall plants have hydraulic limitations that cap practical height). Evolution can select for genomes whose radial-permeability curve sacrifices trunk-tissue access for apical delivery under a given N budget.
```

- [ ] **Step 2: Add Radial Permeability section**

Insert after the Chemical Transport Model section:

```markdown
## Radial Permeability (radius-dependent)

Radial flow between a stem/root's own `local_env` and its own `phloem`/`xylem` is rate-limited by a radius-dependent permeability. Young thin stems are leaky (get plenty of nutrients, grow fast). Mature thick trunks asymptote to a floor (enough to maintain themselves, not enough to siphon flow away from distal apices).

```
radial_permeability(r) = base × (floor + (1 - floor) / (1 + (r / half_radius)²))
```

Independent curves for phloem-radial (sugar) and xylem-radial (water + cytokinin), each with its own genome params: `base_radial_permeability_*`, `radial_floor_fraction_*`, `radial_half_radius_*`.

Combined with fixed N, the curve shape is what makes tall plants viable in the sim — a mature trunk with low radial permeability acts as a hydraulic highway delivering water/sugar to the canopy rather than bleeding it into intermediate stem tissue.
```

- [ ] **Step 3: Update the Tick Control Flow section**

Find `## Tick Control Flow`. Rewrite:

```markdown
## Tick Control Flow

`Plant::tick()` runs two phases in order:

**Phase 1 — Per-node metabolism** (`tick_tree()`): DFS walk from seed. Each `Node::tick()`:
1. `age++`, sync position
2. Produce (photosynthesis into leaf.local, water absorption into root.local, hormone production)
3. Pay maintenance from `local_env`; starvation check
4. `update_tissue()` — virtual, type-specific growth; consumes from `local_env`
5. Sync position, update physics
6. Local signaling diffusion (auxin, gibberellin, stress between neighboring `local_env`s)
7. Decay, ethylene spatial pass

**Phase 2 — Vascular transport** (`vascular_sub_stepped()`): sub-stepped loop of inject / radial / extract / Jacobi — see Chemical Transport Model above.

Ordering is **metabolism first, vascular after**. This gives the model a 1-tick lag between a leaf producing sugar (during Phase 1 of tick N) and that sugar being loaded into phloem (during Phase 2 of tick N+1). At 1-hour tick granularity this is biologically defensible — real plants buffer an hour of sugar in mesophyll cells easily.
```

- [ ] **Step 4: Update Key Design Decisions**

Replace any lines mentioning the old vascular model, `has_vasculature`, `vascular_radius_threshold`, `phloem_iterations`, etc. Add bullets for the new model:

```markdown
- **Compartmented chemical model** — every node has a `local_env`; stems/roots additionally have typed `phloem` and `xylem` TransportPools. Specialty nodes use walk-up helpers to reach their nearest ancestor's conduit.
- **Sub-stepped vascular with fixed N** — N iterations per tick of inject → radial → extract → Jacobi. Fixed N produces realistic hydraulic limitation on tall plants.
- **Tick-then-vascular ordering** — metabolism runs first, vascular after. 1-tick lag accepted at 1-hour granularity.
```

- [ ] **Step 5: Update Tuning Parameters section**

Remove any entries for `vascular_radius_threshold`, `phloem_iterations`, etc. Add entries for the new radial permeability and compartment params, with brief one-line descriptions.

- [ ] **Step 6: Commit**

```bash
git add CLAUDE.md
git commit -m "$(cat <<'EOF'
docs: rewrite CLAUDE.md transport/tick sections for compartmented model

Phase G of the compartmented vascular refactor.  Replaces the old
"Chemical Transport Model" and "Tick Control Flow" sections with
descriptions of the new architecture:

  - Three compartments per node class (local_env, phloem, xylem)
  - Sub-stepped vascular loop (inject → radial → extract → Jacobi)
  - Fixed N vascular_substeps (hydraulic limitation)
  - Radial permeability curve (young-leaky → mature-pipe asymmetry)
  - Tick-then-vascular ordering with 1-tick lag

Plan: docs/superpowers/plans/2026-04-19-compartmented-vascular-model.md
EOF
)"
```

---
