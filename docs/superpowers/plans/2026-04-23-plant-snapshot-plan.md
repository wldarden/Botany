# Plant Snapshot Save/Load — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `.tree` binary save/load so a live plant in `botany_realtime` can be snapshotted to disk and later reloaded to continue simulation from exactly where it left off.

**Architecture:** A new `src/serialization/plant_snapshot.{h,cpp}` module owns the `.tree` format, deliberately separate from the playback-oriented `serializer.*`. Loading constructs a Plant via a new `from_empty` factory + `install_node` API, then rebuilds parent/child pointers and canalization maps. `botany_realtime` gains `--load-plant` and `--genome-override` CLI flags plus `Cmd/Ctrl+S` hotkey for in-viewer saving. A shared `src/engine/genome_io.{h,cpp}` is hoisted out of `app_realtime.cpp` so the snapshot module can use the override-genome parser.

**Tech Stack:** C++17, CMake (FetchContent). Catch2 v3 for tests. GLM for vec3. GLFW for key input. ImGui for the toast.

---

## File structure

**Create:**
- `src/engine/genome_io.h` / `genome_io.cpp` — text-format genome I/O hoisted from `app_realtime.cpp`.
- `src/serialization/plant_snapshot.h` / `plant_snapshot.cpp` — `.tree` binary format, `save_plant_snapshot()`, `load_plant_snapshot()`.
- `tests/test_plant_snapshot.cpp` — round-trip, continuation, override, corrupt-file, canalization tests.

**Modify:**
- `src/engine/plant.h` / `plant.cpp` — add `Plant::from_empty()` factory + `install_node()` + `set_next_id()`.
- `src/engine/engine.h` / `engine.cpp` — add `adopt_plant(std::unique_ptr<Plant>)` + `set_tick(uint32_t)`.
- `src/app_realtime.cpp` — replace local `load_genome_file` with `genome_io` call; parse `--load-plant` + `--genome-override`; install loaded plant; add `Cmd/Ctrl+S` hotkey + toast.
- `CMakeLists.txt` — add the new `.cpp` files to `botany_engine` (genome_io), to a new `botany_serialization` or into `botany_engine` (plant_snapshot), and to `botany_tests`.
- `CLAUDE.md` — add a short "Save/load" section under Apps.

---

## Task 1: Hoist `load_genome_file` into `src/engine/genome_io`

**Files:**
- Create: `src/engine/genome_io.h`
- Create: `src/engine/genome_io.cpp`
- Modify: `src/app_realtime.cpp:34-127` — remove local `load_genome_file`, include new header, call shared function.
- Modify: `CMakeLists.txt` — add `src/engine/genome_io.cpp` to `botany_engine`.
- Test: `tests/test_genome_io.cpp`

- [ ] **Step 1.1: Write the failing test**

Create `tests/test_genome_io.cpp`:
```cpp
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <fstream>
#include "engine/genome.h"
#include "engine/genome_io.h"

using namespace botany;

TEST_CASE("load_genome_file parses key=value pairs", "[genome_io]") {
    auto tmp = std::filesystem::temp_directory_path() / "botany_genome_test.txt";
    {
        std::ofstream f(tmp);
        f << "growth_rate=0.012345\n";
        f << "shoot_plastochron=77\n";
        f << "max_leaf_size=1.5\n";
    }

    Genome g = load_genome_file(tmp.string());
    REQUIRE(g.growth_rate == 0.012345f);
    REQUIRE(g.shoot_plastochron == 77u);
    REQUIRE(g.max_leaf_size == 1.5f);
    // Unspecified fields stay at default.
    REQUIRE(g.branch_angle == default_genome().branch_angle);

    std::filesystem::remove(tmp);
}

TEST_CASE("load_genome_file returns default on missing file", "[genome_io]") {
    Genome g = load_genome_file("/nonexistent/path/never_here.txt");
    REQUIRE(g.growth_rate == default_genome().growth_rate);
}
```

- [ ] **Step 1.2: Add the new cpp to CMake**

Edit `CMakeLists.txt`, add `src/engine/genome_io.cpp` to the `botany_engine` source list (alphabetical by current convention, so after `engine.cpp`). Add `tests/test_genome_io.cpp` to the `botany_tests` source list.

- [ ] **Step 1.3: Run test to verify it fails**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5
```
Expected: compile error — `engine/genome_io.h` not found.

- [ ] **Step 1.4: Create the header**

Create `src/engine/genome_io.h`:
```cpp
#pragma once

#include <string>
#include "engine/genome.h"

namespace botany {

// Parse a ".env"-style text genome file: one "key=value" per line.
// Unrecognized keys are silently ignored (so old files survive field additions).
// Missing file: returns default_genome().
Genome load_genome_file(const std::string& path);

// Write a genome as "key=value" lines.  Writes every field currently parsed by
// load_genome_file so round-trip is lossless for those fields.
// Returns true on success.
bool save_genome_file(const Genome& g, const std::string& path);

} // namespace botany
```

- [ ] **Step 1.5: Create the implementation**

Create `src/engine/genome_io.cpp`. Copy the full body of `load_genome_file` from `src/app_realtime.cpp:34-127` (all ~70 `get_f` / `get_u` calls). Add the symmetric `save_genome_file`:
```cpp
#include "engine/genome_io.h"
#include <fstream>
#include <map>
#include <string>

namespace botany {

Genome load_genome_file(const std::string& path) {
    Genome g = default_genome();
    std::ifstream file(path);
    if (!file) return g;

    std::map<std::string, std::string> fields;
    std::string line;
    while (std::getline(file, line)) {
        auto eq = line.find('=');
        if (eq != std::string::npos) {
            fields[line.substr(0, eq)] = line.substr(eq + 1);
        }
    }

    auto get_f = [&](const char* name, float& out) {
        auto it = fields.find(name);
        if (it != fields.end()) out = std::stof(it->second);
    };
    auto get_u = [&](const char* name, uint32_t& out) {
        auto it = fields.find(name);
        if (it != fields.end()) out = static_cast<uint32_t>(std::stoul(it->second));
    };

    // --- all get_f / get_u calls copied verbatim from app_realtime.cpp:57-124 ---
    // (Keep the order identical so diffs stay small when new genome fields
    //  are added.)
    get_f("apical_auxin_baseline", g.apical_auxin_baseline);
    // ... (every line from app_realtime.cpp lines 57-124 verbatim)
    get_f("vascular_radius_threshold", g.vascular_radius_threshold);

    return g;
}

bool save_genome_file(const Genome& g, const std::string& path) {
    std::ofstream file(path);
    if (!file) return false;

    file << "apical_auxin_baseline=" << g.apical_auxin_baseline << "\n";
    // ... mirror every get_f / get_u from load_genome_file in the same order.
    file << "vascular_radius_threshold=" << g.vascular_radius_threshold << "\n";

    return static_cast<bool>(file);
}

} // namespace botany
```

NOTE: Implementer must copy every field line from `app_realtime.cpp:57-124` into both functions. Use `std::fixed << std::setprecision(9)` or default `operator<<` — default is acceptable for this format since `std::stof` round-trips `operator<<` output for normal floats.

- [ ] **Step 1.6: Update app_realtime to use the shared header**

In `src/app_realtime.cpp`:
1. Add `#include "engine/genome_io.h"` near the other engine includes.
2. Delete the local `static Genome load_genome_file(...)` at lines 34-127.
3. Existing call sites (`g = load_genome_file(...)`) already resolve correctly via the new namespace-scoped function.

- [ ] **Step 1.7: Run test to verify it passes**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[genome_io]" 2>&1 | tail -5
```
Expected: PASS, 2 assertions.

- [ ] **Step 1.8: Run full test suite to verify nothing regressed**

```
./build/botany_tests 2>&1 | tail -3
```
Expected: all tests pass (count rises by 2).

- [ ] **Step 1.9: Commit**

```bash
git add src/engine/genome_io.h src/engine/genome_io.cpp src/app_realtime.cpp CMakeLists.txt tests/test_genome_io.cpp
git commit -m "$(cat <<'EOF'
genome_io: hoist text genome parser into shared module

Moves load_genome_file out of app_realtime.cpp into a new
src/engine/genome_io.h/cpp so plant_snapshot can reuse it for
the upcoming --genome-override flag.  Adds a symmetric
save_genome_file() and a round-trip test.

EOF
)"
```

---

## Task 2: `plant_snapshot` module skeleton + binary I/O helpers

**Files:**
- Create: `src/serialization/plant_snapshot.h`
- Create: `src/serialization/plant_snapshot.cpp`
- Test: `tests/test_plant_snapshot.cpp`
- Modify: `CMakeLists.txt`

- [ ] **Step 2.1: Write a failing test for the binary helpers**

Create `tests/test_plant_snapshot.cpp`:
```cpp
#include <catch2/catch_test_macros.hpp>
#include <sstream>
#include "serialization/plant_snapshot.h"

using namespace botany;

TEST_CASE("plant_snapshot header magic roundtrips", "[plant_snapshot][header]") {
    std::stringstream ss;
    plant_snapshot_write_magic(ss);
    bool ok = plant_snapshot_check_magic(ss);
    REQUIRE(ok);
}

TEST_CASE("plant_snapshot check_magic rejects wrong magic", "[plant_snapshot][header]") {
    std::stringstream ss;
    ss.write("XXXX", 4);
    bool ok = plant_snapshot_check_magic(ss);
    REQUIRE_FALSE(ok);
}
```

- [ ] **Step 2.2: Create the header**

Create `src/serialization/plant_snapshot.h`:
```cpp
#pragma once

#include <cstdint>
#include <iosfwd>
#include <memory>
#include <optional>
#include <string>
#include "engine/genome.h"

namespace botany {

class Plant;

constexpr uint32_t PLANT_SNAPSHOT_VERSION = 1;
// Magic "BTNT" (Botany Tree snapshot).
constexpr char PLANT_SNAPSHOT_MAGIC[4] = {'B','T','N','T'};

// --- Public save/load API (implemented in later tasks) ---

struct SaveResult {
    bool ok = false;
    std::string path;   // populated on success
    std::string error;  // populated on failure
};

// Writes a .tree snapshot of `p` to `dir/plant_YYYYMMDD_HHMMSS_tick<N>.tree`.
// Creates `dir` if missing.  The tick number shown in the filename is engine_tick.
SaveResult save_plant_snapshot(const Plant& p,
                               uint64_t engine_tick,
                               const std::string& dir = "saves");

struct LoadedPlant {
    std::unique_ptr<Plant> plant;
    Genome   genome;       // the genome the plant is running under (file's or override)
    uint64_t engine_tick;  // tick count read from the file header
};

// Reads a .tree snapshot from `path`.  If `genome_override` is set it replaces
// the embedded genome; otherwise the embedded one is used.  Throws std::runtime_error
// with a human-readable message on any failure.
LoadedPlant load_plant_snapshot(const std::string& path,
                                const std::optional<Genome>& genome_override);

// --- Header helpers (exposed for testing; implementation details) ---

void plant_snapshot_write_magic(std::ostream& out);
bool plant_snapshot_check_magic(std::istream& in);

} // namespace botany
```

- [ ] **Step 2.3: Create the implementation with magic helpers + binary helpers**

Create `src/serialization/plant_snapshot.cpp`:
```cpp
#include "serialization/plant_snapshot.h"
#include "engine/plant.h"
#include <cstring>
#include <iostream>
#include <stdexcept>

namespace botany {

// -- Binary helpers (file-local). Values are written host-endian; this format
// is not expected to be portable across differently-endian machines, same as
// the existing serializer.cpp convention.
namespace {
template<typename T>
void write_val(std::ostream& out, const T& v) {
    out.write(reinterpret_cast<const char*>(&v), sizeof(T));
}

template<typename T>
T read_val(std::istream& in) {
    T v;
    in.read(reinterpret_cast<char*>(&v), sizeof(T));
    if (!in) throw std::runtime_error("plant_snapshot: unexpected EOF");
    return v;
}
} // namespace

void plant_snapshot_write_magic(std::ostream& out) {
    out.write(PLANT_SNAPSHOT_MAGIC, 4);
}

bool plant_snapshot_check_magic(std::istream& in) {
    char buf[4];
    in.read(buf, 4);
    if (!in) return false;
    return std::memcmp(buf, PLANT_SNAPSHOT_MAGIC, 4) == 0;
}

// Stubs — filled in by later tasks.
SaveResult save_plant_snapshot(const Plant&, uint64_t, const std::string&) {
    return SaveResult{false, "", "save_plant_snapshot not implemented yet"};
}

LoadedPlant load_plant_snapshot(const std::string&, const std::optional<Genome>&) {
    throw std::runtime_error("load_plant_snapshot not implemented yet");
}

} // namespace botany
```

- [ ] **Step 2.4: Add to CMake**

Edit `CMakeLists.txt`:
1. Add `src/serialization/plant_snapshot.cpp` to the `botany_engine` source list (right after `src/serialization/serializer.cpp` if that's in there; otherwise add to a sensible alphabetical spot).
2. Add `tests/test_plant_snapshot.cpp` to the `botany_tests` source list.

(Check first whether `src/serialization/serializer.cpp` is part of `botany_engine` or a separate lib — match that location.)

- [ ] **Step 2.5: Run test to verify it passes**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[plant_snapshot]" 2>&1 | tail -5
```
Expected: 2 test cases pass.

- [ ] **Step 2.6: Run full test suite**

```
./build/botany_tests 2>&1 | tail -3
```
Expected: all tests pass.

- [ ] **Step 2.7: Commit**

```bash
git add src/serialization/plant_snapshot.h src/serialization/plant_snapshot.cpp tests/test_plant_snapshot.cpp CMakeLists.txt
git commit -m "$(cat <<'EOF'
plant_snapshot: module skeleton with magic-header round-trip

New src/serialization/plant_snapshot.{h,cpp} provides the public
save/load API (stubbed) plus write/check magic helpers and
internal binary read/write templates.  save/load stubs throw
for now; later tasks fill them in.

EOF
)"
```

---

## Task 3: Genome binary serialization

**Files:**
- Modify: `src/serialization/plant_snapshot.cpp` — add `write_genome_binary` / `read_genome_binary`.
- Modify: `src/serialization/plant_snapshot.h` — expose them (internal header OK) for tests.
- Test: `tests/test_plant_snapshot.cpp`

- [ ] **Step 3.1: Write the failing test**

Append to `tests/test_plant_snapshot.cpp`:
```cpp
#include "serialization/plant_snapshot.h"

TEST_CASE("genome binary round-trip is identity", "[plant_snapshot][genome]") {
    Genome original = default_genome();
    original.growth_rate = 0.01234f;
    original.shoot_plastochron = 123u;
    original.meristem_dormancy_death_ticks = 4567u;
    original.max_internode_length = 2.7f;

    std::stringstream ss;
    write_genome_binary(ss, original);

    std::stringstream rs(ss.str());
    Genome loaded = read_genome_binary(rs);

    REQUIRE(loaded.growth_rate == original.growth_rate);
    REQUIRE(loaded.shoot_plastochron == original.shoot_plastochron);
    REQUIRE(loaded.meristem_dormancy_death_ticks == original.meristem_dormancy_death_ticks);
    REQUIRE(loaded.max_internode_length == original.max_internode_length);
    // Canary: default field we didn't touch.
    REQUIRE(loaded.branch_angle == original.branch_angle);
}
```

- [ ] **Step 3.2: Expose the functions in the header**

Add to `src/serialization/plant_snapshot.h` (bottom, before `} // namespace botany`):
```cpp
// --- Genome binary I/O (internal; exposed for tests) ---
// Writes every Genome field in declared order.  Adding or reordering a
// Genome field requires editing these two functions together AND bumping
// PLANT_SNAPSHOT_VERSION.
void   write_genome_binary(std::ostream& out, const Genome& g);
Genome read_genome_binary(std::istream& in);
```

- [ ] **Step 3.3: Implement the functions**

In `src/serialization/plant_snapshot.cpp` (after the magic helpers), add:
```cpp
void write_genome_binary(std::ostream& out, const Genome& g) {
    // Written in declaration order from genome.h.
    // Every field in the Genome struct must appear here AND in read_genome_binary.
    // When a new Genome field is added, add a line to both and bump
    // PLANT_SNAPSHOT_VERSION in the header.
    write_val(out, g.apical_auxin_baseline);
    write_val(out, g.apical_growth_auxin_multiplier);
    write_val(out, g.auxin_diffusion_rate);
    write_val(out, g.auxin_decay_rate);
    write_val(out, g.auxin_threshold);
    write_val(out, g.auxin_shade_boost);
    write_val(out, g.auxin_sugar_half_saturation);
    write_val(out, g.auxin_water_half_saturation);
    write_val(out, g.cytokinin_sugar_half_saturation);
    write_val(out, g.cytokinin_water_half_saturation);
    write_val(out, g.quiescence_threshold);
    write_val(out, g.meristem_dormancy_death_ticks);
    write_val(out, g.starvation_ticks_max_stem);
    write_val(out, g.starvation_ticks_max_root);
    write_val(out, g.auxin_bias);
    write_val(out, g.leaf_auxin_baseline);
    write_val(out, g.leaf_growth_auxin_multiplier);
    write_val(out, g.stem_auxin_max_boost);
    write_val(out, g.stem_auxin_half_saturation);
    write_val(out, g.root_auxin_max_boost);
    write_val(out, g.root_auxin_half_saturation);
    write_val(out, g.leaf_auxin_max_boost);
    write_val(out, g.leaf_auxin_half_saturation);
    write_val(out, g.apical_auxin_max_boost);
    write_val(out, g.apical_auxin_half_saturation);
    write_val(out, g.root_apical_auxin_max_boost);
    write_val(out, g.root_apical_auxin_half_saturation);
    write_val(out, g.cytokinin_production_rate);
    write_val(out, g.cytokinin_diffusion_rate);
    write_val(out, g.cytokinin_decay_rate);
    write_val(out, g.cytokinin_threshold);
    write_val(out, g.cytokinin_growth_threshold);
    write_val(out, g.cytokinin_bias);
    write_val(out, g.hormone_base_transport);
    write_val(out, g.hormone_transport_scale);
    write_val(out, g.sugar_base_transport);
    write_val(out, g.sugar_transport_scale);
    write_val(out, g.stem_photosynthesis_rate);
    write_val(out, g.stem_green_radius_threshold);
    write_val(out, g.growth_rate);
    write_val(out, g.shoot_plastochron);
    write_val(out, g.branch_angle);
    write_val(out, g.cambium_responsiveness);
    write_val(out, g.internode_elongation_rate);
    write_val(out, g.max_internode_length);
    write_val(out, g.internode_maturation_ticks);
    write_val(out, g.root_growth_rate);
    write_val(out, g.root_plastochron);
    write_val(out, g.root_branch_angle);
    write_val(out, g.root_internode_elongation_rate);
    write_val(out, g.root_internode_maturation_ticks);
    write_val(out, g.root_gravitropism_strength);
    write_val(out, g.root_gravitropism_depth);
    write_val(out, g.root_cytokinin_production_rate);
    write_val(out, g.root_tip_auxin_production_rate);
    write_val(out, g.root_auxin_growth_threshold);
    write_val(out, g.root_ck_growth_floor);
    write_val(out, g.root_auxin_activation_threshold);
    write_val(out, g.root_cytokinin_inhibition_threshold);
    write_val(out, g.primary_root_lateral_delay_internodes);
    write_val(out, g.max_leaf_size);
    write_val(out, g.leaf_growth_rate);
    write_val(out, g.leaf_bud_size);
    write_val(out, g.leaf_petiole_length);
    write_val(out, g.leaf_opacity);
    write_val(out, g.initial_radius);
    write_val(out, g.root_initial_radius);
    write_val(out, g.tip_offset);
    write_val(out, g.growth_noise);
    write_val(out, g.sugar_diffusion_rate);
    write_val(out, g.seed_sugar);
    write_val(out, g.sugar_storage_density_wood);
    write_val(out, g.sugar_storage_density_leaf);
    write_val(out, g.sugar_cap_minimum);
    write_val(out, g.sugar_cap_meristem);
    write_val(out, g.water_absorption_rate);
    write_val(out, g.transpiration_rate);
    write_val(out, g.photosynthesis_water_ratio);
    write_val(out, g.water_storage_density_stem);
    write_val(out, g.water_storage_density_leaf);
    write_val(out, g.water_cap_meristem);
    write_val(out, g.water_diffusion_rate);
    write_val(out, g.water_bias);
    write_val(out, g.water_base_transport);
    write_val(out, g.water_transport_scale);
    write_val(out, g.leaf_phototropism_rate);
    write_val(out, g.meristem_gravitropism_rate);
    write_val(out, g.meristem_phototropism_rate);
    write_val(out, g.ga_production_rate);
    write_val(out, g.ga_leaf_age_max);
    write_val(out, g.ga_elongation_sensitivity);
    write_val(out, g.ga_length_sensitivity);
    write_val(out, g.ga_diffusion_rate);
    write_val(out, g.ga_decay_rate);
    write_val(out, g.leaf_abscission_ticks);
    write_val(out, g.min_leaf_age_before_abscission);
    write_val(out, g.ethylene_starvation_rate);
    write_val(out, g.ethylene_starvation_tick_threshold);
    write_val(out, g.ethylene_shade_rate);
    write_val(out, g.ethylene_shade_threshold);
    write_val(out, g.ethylene_age_rate);
    write_val(out, g.ethylene_age_onset);
    write_val(out, g.ethylene_crowding_rate);
    write_val(out, g.ethylene_crowding_radius);
    write_val(out, g.ethylene_diffusion_radius);
    write_val(out, g.ethylene_abscission_threshold);
    write_val(out, g.ethylene_elongation_inhibition);
    write_val(out, g.senescence_duration);
    write_val(out, g.wood_density);
    write_val(out, g.wood_flexibility);
    write_val(out, g.stress_hormone_threshold);
    write_val(out, g.stress_hormone_production_rate);
    write_val(out, g.stress_hormone_diffusion_rate);
    write_val(out, g.stress_hormone_decay_rate);
    write_val(out, g.stress_thickening_boost);
    write_val(out, g.stress_elongation_inhibition);
    write_val(out, g.stress_gravitropism_boost);
    write_val(out, g.elastic_recovery_rate);
    write_val(out, g.smoothing_rate);
    write_val(out, g.canalization_weight);
    write_val(out, g.pin_capacity_per_area);
    write_val(out, g.pin_base_efficiency);
    write_val(out, g.meristem_sink_fraction);
    write_val(out, g.vascular_radius_threshold);
    write_val(out, g.base_radial_permeability_sugar);
    write_val(out, g.radial_floor_fraction_sugar);
    write_val(out, g.radial_half_radius_sugar);
    write_val(out, g.base_radial_permeability_water);
    write_val(out, g.radial_floor_fraction_water);
    write_val(out, g.radial_half_radius_water);
    write_val(out, g.phloem_fraction);
    write_val(out, g.xylem_fraction);
    write_val(out, g.leaf_reserve_fraction_sugar);
    write_val(out, g.meristem_sink_target_fraction);
    write_val(out, g.leaf_turgor_target_fraction);
    write_val(out, g.root_water_reserve_fraction);
}

Genome read_genome_binary(std::istream& in) {
    Genome g = default_genome();
    // Mirror write_genome_binary exactly.  Keep the order identical.
    g.apical_auxin_baseline = read_val<float>(in);
    g.apical_growth_auxin_multiplier = read_val<float>(in);
    g.auxin_diffusion_rate = read_val<float>(in);
    // ... (mirror every field in the same order)
    g.root_water_reserve_fraction = read_val<float>(in);
    return g;
}
```

NOTE for implementer: the write list above is the authoritative field list for v1. If the current `genome.h` has ADDED a field since this plan was written, add it to BOTH functions AND bump `PLANT_SNAPSHOT_VERSION`. Check with:
```
grep -cE "^    (float|uint32_t|int)" src/engine/genome.h
```
— the count should equal the number of `write_val` lines.

- [ ] **Step 3.4: Verify field count**

```
cd /Users/wldarden/learning/botany
echo "Genome scalar fields: $(grep -cE '^    (float|uint32_t|int)' src/engine/genome.h)"
echo "write_val lines:       $(grep -c 'write_val' src/serialization/plant_snapshot.cpp)"
```
The numbers must match (both ~127). Mismatch means a field was missed.

- [ ] **Step 3.5: Build + run tests**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[plant_snapshot]" 2>&1 | tail -5
```
Expected: 3 test cases pass.

- [ ] **Step 3.6: Commit**

```bash
git add src/serialization/plant_snapshot.h src/serialization/plant_snapshot.cpp tests/test_plant_snapshot.cpp
git commit -m "$(cat <<'EOF'
plant_snapshot: field-by-field Genome binary I/O

write_genome_binary / read_genome_binary serialize every Genome
field in declared order, avoiding cross-platform struct padding
issues.  Adding a new Genome field requires editing both
functions in tandem and bumping PLANT_SNAPSHOT_VERSION.

EOF
)"
```

---

## Task 4: Node schema — common fields + local_chemicals + auxin_flow_bias

**Files:**
- Modify: `src/serialization/plant_snapshot.cpp` — add `write_node_common`, `read_node_common`, `write_chem_map`, `read_chem_map`, `write_bias_map`, `read_bias_map`.

- [ ] **Step 4.1: Write failing test — common fields round-trip**

Append to `tests/test_plant_snapshot.cpp`:
```cpp
#include "engine/plant.h"
#include "engine/node/node.h"
#include "engine/node/stem_node.h"

TEST_CASE("node common fields round-trip through binary", "[plant_snapshot][node]") {
    // Build one StemNode and populate every common field.
    StemNode s(42, glm::vec3(1.0f, 2.0f, 3.0f), 0.025f);
    s.parent = nullptr;
    s.offset      = glm::vec3(0.1f, 0.2f, 0.3f);
    s.rest_offset = glm::vec3(0.4f, 0.5f, 0.6f);
    s.position    = glm::vec3(1.0f, 2.0f, 3.0f);
    s.age = 77;
    s.starvation_ticks = 11;
    s.dormant_ticks = 0;
    s.ever_active = true;
    s.local().chemical(ChemicalID::Sugar)  = 1.5f;
    s.local().chemical(ChemicalID::Auxin)  = 0.25f;
    s.local().chemical(ChemicalID::Water)  = 2.0f;

    std::stringstream ss;
    write_node_common(ss, s, /*parent_id=*/ UINT32_MAX);

    std::stringstream rs(ss.str());
    NodeCommonRecord rec = read_node_common(rs);

    REQUIRE(rec.id == 42u);
    REQUIRE(rec.parent_id == UINT32_MAX);
    REQUIRE(rec.type == NodeType::STEM);
    REQUIRE(rec.age == 77u);
    REQUIRE(rec.starvation_ticks == 11u);
    REQUIRE(rec.dormant_ticks == 0u);
    REQUIRE(rec.radius == 0.025f);
    REQUIRE(rec.position == glm::vec3(1.0f, 2.0f, 3.0f));
    REQUIRE(rec.offset == glm::vec3(0.1f, 0.2f, 0.3f));
    REQUIRE(rec.rest_offset == glm::vec3(0.4f, 0.5f, 0.6f));
    REQUIRE(rec.ever_active == true);
    REQUIRE(rec.local_chemicals.at(ChemicalID::Sugar) == 1.5f);
    REQUIRE(rec.local_chemicals.at(ChemicalID::Auxin) == 0.25f);
    REQUIRE(rec.local_chemicals.at(ChemicalID::Water) == 2.0f);
}
```

- [ ] **Step 4.2: Add `NodeCommonRecord` and function signatures to header**

In `src/serialization/plant_snapshot.h`, above the genome-binary section:
```cpp
#include <unordered_map>
#include "engine/node/node.h"

// Per-node common-field record returned by read_node_common.  The snapshot
// loader uses this to populate the allocated subclass, then reads any
// type-specific trailer.
struct NodeCommonRecord {
    uint32_t  id;
    uint32_t  parent_id;
    NodeType  type;
    uint32_t  age;
    uint32_t  starvation_ticks;
    uint32_t  dormant_ticks;
    glm::vec3 position;
    glm::vec3 offset;
    glm::vec3 rest_offset;
    float     radius;
    bool      ever_active;
    std::unordered_map<ChemicalID, float> local_chemicals;
    std::unordered_map<uint32_t,   float> auxin_flow_bias; // keyed by child_id
};

void             write_node_common(std::ostream& out, const class Node& n, uint32_t parent_id);
NodeCommonRecord read_node_common(std::istream& in);
```

- [ ] **Step 4.3: Implement the helpers in `plant_snapshot.cpp`**

Add (after `read_genome_binary`):
```cpp
namespace {

void write_chem_map(std::ostream& out, const std::unordered_map<ChemicalID, float>& m) {
    uint16_t count = static_cast<uint16_t>(m.size());
    write_val(out, count);
    for (const auto& kv : m) {
        write_val(out, static_cast<uint8_t>(kv.first));
        write_val(out, kv.second);
    }
}

std::unordered_map<ChemicalID, float> read_chem_map(std::istream& in) {
    uint16_t count = read_val<uint16_t>(in);
    std::unordered_map<ChemicalID, float> m;
    m.reserve(count);
    for (uint16_t i = 0; i < count; i++) {
        auto id = static_cast<ChemicalID>(read_val<uint8_t>(in));
        float v = read_val<float>(in);
        m[id] = v;
    }
    return m;
}

void write_bias_map(std::ostream& out, const std::unordered_map<Node*, float>& m) {
    uint16_t count = static_cast<uint16_t>(m.size());
    write_val(out, count);
    for (const auto& kv : m) {
        // Safe: by the time we serialize, every node has a valid id.
        uint32_t child_id = kv.first ? kv.first->id : UINT32_MAX;
        write_val(out, child_id);
        write_val(out, kv.second);
    }
}

std::unordered_map<uint32_t, float> read_bias_map(std::istream& in) {
    uint16_t count = read_val<uint16_t>(in);
    std::unordered_map<uint32_t, float> m;
    m.reserve(count);
    for (uint16_t i = 0; i < count; i++) {
        uint32_t id = read_val<uint32_t>(in);
        float v    = read_val<float>(in);
        m[id] = v;
    }
    return m;
}

} // namespace

void write_node_common(std::ostream& out, const Node& n, uint32_t parent_id) {
    write_val(out, n.id);
    write_val(out, parent_id);
    write_val(out, static_cast<uint8_t>(n.type));
    write_val(out, n.age);
    write_val(out, n.starvation_ticks);
    write_val(out, n.dormant_ticks);
    write_val(out, n.position);
    write_val(out, n.offset);
    write_val(out, n.rest_offset);
    write_val(out, n.radius);
    write_val(out, static_cast<uint8_t>(n.ever_active ? 1 : 0));
    write_chem_map(out, n.local().chemicals);
    write_bias_map(out, n.auxin_flow_bias);
}

NodeCommonRecord read_node_common(std::istream& in) {
    NodeCommonRecord r;
    r.id               = read_val<uint32_t>(in);
    r.parent_id        = read_val<uint32_t>(in);
    r.type             = static_cast<NodeType>(read_val<uint8_t>(in));
    r.age              = read_val<uint32_t>(in);
    r.starvation_ticks = read_val<uint32_t>(in);
    r.dormant_ticks    = read_val<uint32_t>(in);
    r.position         = read_val<glm::vec3>(in);
    r.offset           = read_val<glm::vec3>(in);
    r.rest_offset      = read_val<glm::vec3>(in);
    r.radius           = read_val<float>(in);
    r.ever_active      = read_val<uint8_t>(in) != 0;
    r.local_chemicals  = read_chem_map(in);
    r.auxin_flow_bias  = read_bias_map(in);
    return r;
}
```

- [ ] **Step 4.4: Build + test**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[plant_snapshot]" 2>&1 | tail -5
```
Expected: 4 test cases pass.

- [ ] **Step 4.5: Commit**

```bash
git add src/serialization/plant_snapshot.h src/serialization/plant_snapshot.cpp tests/test_plant_snapshot.cpp
git commit -m "plant_snapshot: common node-field binary I/O"
```

---

## Task 5: Node schema — conduit pools + type-specific trailers

**Files:**
- Modify: `src/serialization/plant_snapshot.h` — type-specific record structs + function signatures.
- Modify: `src/serialization/plant_snapshot.cpp` — write/read for each type trailer and conduit pools.

- [ ] **Step 5.1: Write failing tests for each type trailer**

Append to `tests/test_plant_snapshot.cpp`:
```cpp
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"

TEST_CASE("leaf trailer round-trips", "[plant_snapshot][node]") {
    LeafNode l(10, glm::vec3(0.0f), 0.01f);
    l.leaf_size = 0.5f;
    l.light_exposure = 0.75f;
    l.senescence_ticks = 12;
    l.deficit_ticks = 5;
    l.facing = glm::vec3(0.0f, 0.707f, 0.707f);

    std::stringstream ss;
    write_leaf_trailer(ss, l);

    std::stringstream rs(ss.str());
    LeafTrailer t = read_leaf_trailer(rs);
    REQUIRE(t.leaf_size == 0.5f);
    REQUIRE(t.light_exposure == 0.75f);
    REQUIRE(t.senescence_ticks == 12u);
    REQUIRE(t.deficit_ticks == 5u);
    REQUIRE(t.facing == glm::vec3(0.0f, 0.707f, 0.707f));
}

TEST_CASE("apical trailer round-trips", "[plant_snapshot][node]") {
    ApicalNode a(11, glm::vec3(0.0f), 0.01f);
    a.active = false;
    a.is_primary = true;
    a.growth_dir = glm::vec3(0.1f, 0.9f, 0.0f);
    a.ticks_since_last_node = 7;

    std::stringstream ss;
    write_apical_trailer(ss, a);

    std::stringstream rs(ss.str());
    ApicalTrailer t = read_apical_trailer(rs);
    REQUIRE(t.active == false);
    REQUIRE(t.is_primary == true);
    REQUIRE(t.growth_dir == glm::vec3(0.1f, 0.9f, 0.0f));
    REQUIRE(t.ticks_since_last_node == 7u);
}

TEST_CASE("root_apical trailer round-trips", "[plant_snapshot][node]") {
    RootApicalNode r(12, glm::vec3(0.0f), 0.005f);
    r.active = true;
    r.is_primary = true;
    r.growth_dir = glm::vec3(0.0f, -1.0f, 0.0f);
    r.ticks_since_last_node = 3;
    r.internodes_spawned = 8;

    std::stringstream ss;
    write_root_apical_trailer(ss, r);

    std::stringstream rs(ss.str());
    RootApicalTrailer t = read_root_apical_trailer(rs);
    REQUIRE(t.active == true);
    REQUIRE(t.is_primary == true);
    REQUIRE(t.growth_dir == glm::vec3(0.0f, -1.0f, 0.0f));
    REQUIRE(t.ticks_since_last_node == 3u);
    REQUIRE(t.internodes_spawned == 8u);
}

TEST_CASE("conduit pool chemicals round-trip", "[plant_snapshot][node]") {
    StemNode s(42, glm::vec3(0.0f), 0.02f);
    s.phloem()->chemical(ChemicalID::Sugar) = 0.05f;
    s.xylem()->chemical(ChemicalID::Water) = 0.3f;
    s.xylem()->chemical(ChemicalID::Cytokinin) = 0.015f;

    std::stringstream ss;
    write_conduit_pools(ss, s);

    std::stringstream rs(ss.str());
    ConduitPools p = read_conduit_pools(rs);
    REQUIRE(p.phloem.at(ChemicalID::Sugar) == 0.05f);
    REQUIRE(p.xylem.at(ChemicalID::Water)  == 0.3f);
    REQUIRE(p.xylem.at(ChemicalID::Cytokinin) == 0.015f);
}
```

- [ ] **Step 5.2: Add record structs + signatures to header**

In `src/serialization/plant_snapshot.h`:
```cpp
struct LeafTrailer {
    float     leaf_size;
    float     light_exposure;
    uint32_t  senescence_ticks;
    uint32_t  deficit_ticks;
    glm::vec3 facing;
};

struct ApicalTrailer {
    bool      active;
    bool      is_primary;
    glm::vec3 growth_dir;
    uint32_t  ticks_since_last_node;
};

struct RootApicalTrailer {
    bool      active;
    bool      is_primary;
    glm::vec3 growth_dir;
    uint32_t  ticks_since_last_node;
    uint32_t  internodes_spawned;
};

struct ConduitPools {
    std::unordered_map<ChemicalID, float> phloem;
    std::unordered_map<ChemicalID, float> xylem;
};

class LeafNode;
class ApicalNode;
class RootApicalNode;

void write_leaf_trailer(std::ostream& out, const LeafNode& l);
LeafTrailer read_leaf_trailer(std::istream& in);

void write_apical_trailer(std::ostream& out, const ApicalNode& a);
ApicalTrailer read_apical_trailer(std::istream& in);

void write_root_apical_trailer(std::ostream& out, const RootApicalNode& r);
RootApicalTrailer read_root_apical_trailer(std::istream& in);

void write_conduit_pools(std::ostream& out, const Node& n);
ConduitPools read_conduit_pools(std::istream& in);
```

- [ ] **Step 5.3: Implement in `plant_snapshot.cpp`**

Add includes at the top:
```cpp
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"
```

Add (after read_node_common):
```cpp
void write_leaf_trailer(std::ostream& out, const LeafNode& l) {
    write_val(out, l.leaf_size);
    write_val(out, l.light_exposure);
    write_val(out, l.senescence_ticks);
    write_val(out, l.deficit_ticks);
    write_val(out, l.facing);
}

LeafTrailer read_leaf_trailer(std::istream& in) {
    LeafTrailer t;
    t.leaf_size        = read_val<float>(in);
    t.light_exposure   = read_val<float>(in);
    t.senescence_ticks = read_val<uint32_t>(in);
    t.deficit_ticks    = read_val<uint32_t>(in);
    t.facing           = read_val<glm::vec3>(in);
    return t;
}

void write_apical_trailer(std::ostream& out, const ApicalNode& a) {
    write_val(out, static_cast<uint8_t>(a.active ? 1 : 0));
    write_val(out, static_cast<uint8_t>(a.is_primary ? 1 : 0));
    write_val(out, a.growth_dir);
    write_val(out, a.ticks_since_last_node);
}

ApicalTrailer read_apical_trailer(std::istream& in) {
    ApicalTrailer t;
    t.active                = read_val<uint8_t>(in) != 0;
    t.is_primary            = read_val<uint8_t>(in) != 0;
    t.growth_dir            = read_val<glm::vec3>(in);
    t.ticks_since_last_node = read_val<uint32_t>(in);
    return t;
}

void write_root_apical_trailer(std::ostream& out, const RootApicalNode& r) {
    write_val(out, static_cast<uint8_t>(r.active ? 1 : 0));
    write_val(out, static_cast<uint8_t>(r.is_primary ? 1 : 0));
    write_val(out, r.growth_dir);
    write_val(out, r.ticks_since_last_node);
    write_val(out, r.internodes_spawned);
}

RootApicalTrailer read_root_apical_trailer(std::istream& in) {
    RootApicalTrailer t;
    t.active                = read_val<uint8_t>(in) != 0;
    t.is_primary            = read_val<uint8_t>(in) != 0;
    t.growth_dir            = read_val<glm::vec3>(in);
    t.ticks_since_last_node = read_val<uint32_t>(in);
    t.internodes_spawned    = read_val<uint32_t>(in);
    return t;
}

void write_conduit_pools(std::ostream& out, const Node& n) {
    const TransportPool* p = n.phloem();
    const TransportPool* x = n.xylem();
    // Serialize empty maps when the pool is absent — caller already gated on
    // node type, but be defensive.
    static const std::unordered_map<ChemicalID, float> kEmpty;
    write_chem_map(out, p ? p->chemicals : kEmpty);
    write_chem_map(out, x ? x->chemicals : kEmpty);
}

ConduitPools read_conduit_pools(std::istream& in) {
    ConduitPools p;
    p.phloem = read_chem_map(in);
    p.xylem  = read_chem_map(in);
    return p;
}
```

- [ ] **Step 5.4: Build + test**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[plant_snapshot][node]" 2>&1 | tail -5
```
Expected: all 4 node trailer tests pass.

- [ ] **Step 5.5: Commit**

```bash
git add src/serialization/plant_snapshot.h src/serialization/plant_snapshot.cpp tests/test_plant_snapshot.cpp
git commit -m "plant_snapshot: conduit pools + per-type node trailers"
```

---

## Task 6: Plant API — `from_empty` factory + `install_node` + `set_next_id`

**Files:**
- Modify: `src/engine/plant.h` — add the new public methods.
- Modify: `src/engine/plant.cpp` — implement them.
- Test: `tests/test_plant_snapshot.cpp`

- [ ] **Step 6.1: Write failing test**

Append to `tests/test_plant_snapshot.cpp`:
```cpp
TEST_CASE("Plant::from_empty yields an empty plant", "[plant_snapshot][plant]") {
    auto p = Plant::from_empty(default_genome());
    REQUIRE(p != nullptr);
    REQUIRE(p->node_count() == 0);
}

TEST_CASE("Plant::install_node + set_next_id populate nodes manually", "[plant_snapshot][plant]") {
    auto p = Plant::from_empty(default_genome());
    auto s = std::make_unique<StemNode>(99, glm::vec3(0.0f), 0.02f);
    Node* raw = s.get();
    p->install_node(std::move(s));
    p->set_next_id(100);
    REQUIRE(p->node_count() == 1);
    REQUIRE(p->seed() == raw);
    REQUIRE(p->next_id() == 100u); // next_id() increments; now 101
}
```

- [ ] **Step 6.2: Add declarations to `plant.h`**

Insert into the `public:` section (after `Plant(const Genome& genome, glm::vec3 position);`):
```cpp
    // Factory for loading from a snapshot: constructs a Plant with the given
    // genome but no default seed/SA/RA.  Caller must install_node() for every
    // node from the snapshot and then set_next_id() to match the snapshot's
    // next-id counter before the plant is ticked.
    static std::unique_ptr<Plant> from_empty(const Genome& genome);

    // Install a pre-constructed node.  Takes ownership.  Caller is responsible
    // for setting parent / children links, chemicals, and any subclass state
    // before this plant is ticked.  Used by the snapshot loader.
    void install_node(std::unique_ptr<Node> node);

    // Override the id counter so future create_node() calls don't collide with
    // ids already loaded from a snapshot.
    void set_next_id(uint32_t next);
```

- [ ] **Step 6.3: Add the private tag-dispatch constructor to `plant.h`**

In the `private:` section of `plant.h` (near the existing ctor helpers):
```cpp
    // Snapshot-load constructor: tag-dispatched so it doesn't collide with the
    // normal (Genome, vec3) seeding constructor.  Leaves nodes_ empty.  The
    // loader populates nodes via install_node() and then calls set_next_id().
    Plant(const Genome& genome, bool /*empty_tag*/);
```

- [ ] **Step 6.4: Implement the empty ctor, factory, and helpers in `plant.cpp`**

Append to `plant.cpp` (after the existing public constructor body):
```cpp
Plant::Plant(const Genome& genome, bool /*empty_tag*/)
    : genome_(genome)
{
    // Intentionally empty — snapshot loader populates nodes_ and sets next_id_.
}

std::unique_ptr<Plant> Plant::from_empty(const Genome& genome) {
    return std::unique_ptr<Plant>(new Plant(genome, /*empty_tag=*/true));
}

void Plant::install_node(std::unique_ptr<Node> node) {
    nodes_.push_back(std::move(node));
}

void Plant::set_next_id(uint32_t next) {
    next_id_ = next;
}
```

- [ ] **Step 6.5: Build + test**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[plant_snapshot][plant]" 2>&1 | tail -5
```
Expected: both new tests pass.

- [ ] **Step 6.6: Full-suite regression check**

```
./build/botany_tests 2>&1 | tail -3
```
Expected: all pre-existing tests still pass.

- [ ] **Step 6.7: Commit**

```bash
git add src/engine/plant.h src/engine/plant.cpp tests/test_plant_snapshot.cpp
git commit -m "plant: from_empty factory + install_node/set_next_id for snapshot load"
```

---

## Task 7: `save_plant_snapshot()` — header + per-node walk

**Files:**
- Modify: `src/serialization/plant_snapshot.cpp`.
- Test: `tests/test_plant_snapshot.cpp` (intermediate "file exists and header checks out").

- [ ] **Step 7.1: Write failing test — save produces a file with a valid header**

Append:
```cpp
#include <filesystem>

TEST_CASE("save_plant_snapshot writes a file with valid header", "[plant_snapshot][save]") {
    Plant plant(default_genome(), glm::vec3(0.0f));
    auto tmp = std::filesystem::temp_directory_path() / "botany_snap_test";
    std::filesystem::remove_all(tmp);

    SaveResult r = save_plant_snapshot(plant, /*engine_tick=*/42, tmp.string());
    REQUIRE(r.ok);
    REQUIRE(!r.path.empty());
    REQUIRE(std::filesystem::exists(r.path));

    // Re-open and check magic + version + tick.
    std::ifstream in(r.path, std::ios::binary);
    REQUIRE(plant_snapshot_check_magic(in));
    uint32_t version = 0;
    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    REQUIRE(version == PLANT_SNAPSHOT_VERSION);
    uint64_t engine_tick = 0;
    in.read(reinterpret_cast<char*>(&engine_tick), sizeof(engine_tick));
    REQUIRE(engine_tick == 42u);

    std::filesystem::remove_all(tmp);
}
```

- [ ] **Step 7.2: Implement `save_plant_snapshot`**

Replace the stub in `plant_snapshot.cpp`:
```cpp
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace {

std::string make_timestamp_now() {
    auto t = std::chrono::system_clock::now();
    std::time_t tt = std::chrono::system_clock::to_time_t(t);
    std::tm lt{};
#if defined(_WIN32)
    localtime_s(&lt, &tt);
#else
    localtime_r(&tt, &lt);
#endif
    std::ostringstream os;
    os << std::put_time(&lt, "%Y%m%d_%H%M%S");
    return os.str();
}

// DFS from seed so parent always precedes children in the file.
void dfs_collect(const Node* n, std::vector<const Node*>& out) {
    if (!n) return;
    out.push_back(n);
    for (const Node* c : n->children) dfs_collect(c, out);
}

void write_node_full(std::ostream& out, const Node& n) {
    uint32_t parent_id = n.parent ? n.parent->id : UINT32_MAX;
    write_node_common(out, n, parent_id);
    if (n.type == NodeType::STEM || n.type == NodeType::ROOT) {
        write_conduit_pools(out, n);
    }
    switch (n.type) {
        case NodeType::LEAF:        write_leaf_trailer(out, *n.as_leaf()); break;
        case NodeType::APICAL:      write_apical_trailer(out, *n.as_apical()); break;
        case NodeType::ROOT_APICAL: write_root_apical_trailer(out, *n.as_root_apical()); break;
        case NodeType::STEM:
        case NodeType::ROOT:        /* no trailer beyond conduit pools */ break;
    }
}

} // namespace

SaveResult save_plant_snapshot(const Plant& p, uint64_t engine_tick, const std::string& dir) {
    SaveResult r;
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) { r.error = "cannot create dir: " + ec.message(); return r; }

    std::ostringstream name;
    name << "plant_" << make_timestamp_now() << "_tick" << engine_tick << ".tree";
    auto full = std::filesystem::path(dir) / name.str();

    std::ofstream out(full, std::ios::binary);
    if (!out) { r.error = "cannot open output: " + full.string(); return r; }

    // Header
    plant_snapshot_write_magic(out);
    write_val(out, PLANT_SNAPSHOT_VERSION);
    write_val(out, engine_tick);
    write_genome_binary(out, p.genome());

    // Collect nodes in DFS order (parents precede children).
    std::vector<const Node*> order;
    dfs_collect(p.seed(), order);

    // next_node_id: one past the max id we saw, consistent with create_node contract.
    uint32_t max_id = 0;
    for (const Node* n : order) if (n->id >= max_id) max_id = n->id;
    uint32_t next_id = max_id + 1;

    uint32_t count = static_cast<uint32_t>(order.size());
    write_val(out, count);
    write_val(out, next_id);

    for (const Node* n : order) write_node_full(out, *n);

    out.flush();
    if (!out) { r.error = "write failed"; return r; }

    r.ok = true;
    r.path = full.string();
    return r;
}
```

- [ ] **Step 7.3: Build + test**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[plant_snapshot][save]" 2>&1 | tail -5
```
Expected: save test passes.

- [ ] **Step 7.4: Commit**

```bash
git add src/serialization/plant_snapshot.cpp tests/test_plant_snapshot.cpp
git commit -m "plant_snapshot: save writes magic+version+tick+genome+nodes"
```

---

## Task 8: `load_plant_snapshot()` — three-pass reconstruction

**Files:**
- Modify: `src/serialization/plant_snapshot.cpp`.
- Test: `tests/test_plant_snapshot.cpp` — round-trip identity test.

- [ ] **Step 8.1: Write failing test — round-trip identity on a tiny grown plant**

Append:
```cpp
#include "engine/engine.h"

static void deep_compare_nodes(const Node& a, const Node& b) {
    REQUIRE(a.id == b.id);
    REQUIRE(a.type == b.type);
    REQUIRE(a.age == b.age);
    REQUIRE(a.starvation_ticks == b.starvation_ticks);
    REQUIRE(a.dormant_ticks == b.dormant_ticks);
    REQUIRE(a.radius == b.radius);
    REQUIRE(a.position == b.position);
    REQUIRE(a.offset == b.offset);
    REQUIRE(a.rest_offset == b.rest_offset);
    REQUIRE(a.ever_active == b.ever_active);
    REQUIRE(a.local().chemicals.size() == b.local().chemicals.size());
    for (const auto& kv : a.local().chemicals) {
        REQUIRE(b.local().chemical(kv.first) == kv.second);
    }
    // Parent id match (NOT pointer)
    uint32_t ap = a.parent ? a.parent->id : UINT32_MAX;
    uint32_t bp = b.parent ? b.parent->id : UINT32_MAX;
    REQUIRE(ap == bp);
}

TEST_CASE("save+load round-trip preserves tiny plant", "[plant_snapshot][roundtrip]") {
    Engine engine;
    engine.create_plant(default_genome(), glm::vec3(0.0f));
    for (int i = 0; i < 20; i++) engine.tick();

    const Plant& original = engine.get_plant(0);
    auto tmp = std::filesystem::temp_directory_path() / "botany_rt_test";
    std::filesystem::remove_all(tmp);

    SaveResult sr = save_plant_snapshot(original, engine.get_tick(), tmp.string());
    REQUIRE(sr.ok);

    LoadedPlant lp = load_plant_snapshot(sr.path, std::nullopt);
    REQUIRE(lp.plant != nullptr);
    REQUIRE(lp.plant->node_count() == original.node_count());
    REQUIRE(lp.engine_tick == engine.get_tick());

    // Compare every node by id.
    std::unordered_map<uint32_t, const Node*> orig_by_id;
    original.for_each_node([&](const Node& n) { orig_by_id[n.id] = &n; });
    lp.plant->for_each_node([&](const Node& n) {
        auto it = orig_by_id.find(n.id);
        REQUIRE(it != orig_by_id.end());
        deep_compare_nodes(*it->second, n);
    });

    std::filesystem::remove_all(tmp);
}
```

- [ ] **Step 8.2: Implement `load_plant_snapshot`**

Replace the stub with:
```cpp
LoadedPlant load_plant_snapshot(const std::string& path,
                                const std::optional<Genome>& genome_override) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("plant_snapshot: cannot open " + path);

    if (!plant_snapshot_check_magic(in))
        throw std::runtime_error("plant_snapshot: bad magic in " + path);

    uint32_t version = read_val<uint32_t>(in);
    if (version != PLANT_SNAPSHOT_VERSION)
        throw std::runtime_error("plant_snapshot: unsupported version " + std::to_string(version));

    uint64_t engine_tick = read_val<uint64_t>(in);
    Genome file_genome   = read_genome_binary(in);
    uint32_t node_count  = read_val<uint32_t>(in);
    uint32_t next_id     = read_val<uint32_t>(in);

    Genome active = genome_override ? *genome_override : file_genome;

    auto plant = Plant::from_empty(active);

    struct LoadedRecord {
        NodeCommonRecord common;
        std::optional<ConduitPools> pools;
        std::optional<LeafTrailer>  leaf;
        std::optional<ApicalTrailer> apical;
        std::optional<RootApicalTrailer> root_apical;
    };

    std::vector<LoadedRecord> records;
    records.reserve(node_count);

    // Pass 1: read records + allocate subclass + populate common state + chemicals.
    std::unordered_map<uint32_t, Node*> by_id;
    for (uint32_t i = 0; i < node_count; i++) {
        LoadedRecord rec;
        rec.common = read_node_common(in);
        if (rec.common.type == NodeType::STEM || rec.common.type == NodeType::ROOT)
            rec.pools = read_conduit_pools(in);
        switch (rec.common.type) {
            case NodeType::LEAF:        rec.leaf = read_leaf_trailer(in); break;
            case NodeType::APICAL:      rec.apical = read_apical_trailer(in); break;
            case NodeType::ROOT_APICAL: rec.root_apical = read_root_apical_trailer(in); break;
            case NodeType::STEM:
            case NodeType::ROOT:        break;
        }
        records.push_back(std::move(rec));
    }

    // Allocate + fill common state.
    for (const LoadedRecord& rec : records) {
        std::unique_ptr<Node> n;
        switch (rec.common.type) {
            case NodeType::STEM:        n = std::make_unique<StemNode>(rec.common.id, rec.common.position, rec.common.radius); break;
            case NodeType::ROOT:        n = std::make_unique<RootNode>(rec.common.id, rec.common.position, rec.common.radius); break;
            case NodeType::LEAF:        n = std::make_unique<LeafNode>(rec.common.id, rec.common.position, rec.common.radius); break;
            case NodeType::APICAL:      n = std::make_unique<ApicalNode>(rec.common.id, rec.common.position, rec.common.radius); break;
            case NodeType::ROOT_APICAL: n = std::make_unique<RootApicalNode>(rec.common.id, rec.common.position, rec.common.radius); break;
        }
        n->age              = rec.common.age;
        n->starvation_ticks = rec.common.starvation_ticks;
        n->dormant_ticks    = rec.common.dormant_ticks;
        n->offset           = rec.common.offset;
        n->rest_offset      = rec.common.rest_offset;
        n->position         = rec.common.position;
        n->ever_active      = rec.common.ever_active;
        n->local().chemicals = rec.common.local_chemicals;

        if (rec.pools) {
            if (auto* p = n->phloem()) p->chemicals = rec.pools->phloem;
            if (auto* x = n->xylem())  x->chemicals = rec.pools->xylem;
        }
        if (rec.leaf) {
            auto* l = n->as_leaf();
            l->leaf_size        = rec.leaf->leaf_size;
            l->light_exposure   = rec.leaf->light_exposure;
            l->senescence_ticks = rec.leaf->senescence_ticks;
            l->deficit_ticks    = rec.leaf->deficit_ticks;
            l->facing           = rec.leaf->facing;
        }
        if (rec.apical) {
            auto* a = n->as_apical();
            a->active                = rec.apical->active;
            a->is_primary            = rec.apical->is_primary;
            a->growth_dir            = rec.apical->growth_dir;
            a->ticks_since_last_node = rec.apical->ticks_since_last_node;
        }
        if (rec.root_apical) {
            auto* r = n->as_root_apical();
            r->active                = rec.root_apical->active;
            r->is_primary            = rec.root_apical->is_primary;
            r->growth_dir            = rec.root_apical->growth_dir;
            r->ticks_since_last_node = rec.root_apical->ticks_since_last_node;
            r->internodes_spawned    = rec.root_apical->internodes_spawned;
        }

        by_id[rec.common.id] = n.get();
        plant->install_node(std::move(n));
    }

    // Pass 2: wire parent/children pointers.
    for (const LoadedRecord& rec : records) {
        Node* self = by_id.at(rec.common.id);
        if (rec.common.parent_id == UINT32_MAX) continue;
        auto pit = by_id.find(rec.common.parent_id);
        if (pit == by_id.end())
            throw std::runtime_error("plant_snapshot: node " + std::to_string(rec.common.id)
                                   + " references unknown parent " + std::to_string(rec.common.parent_id));
        self->parent = pit->second;
        pit->second->children.push_back(self);
    }

    // Pass 3: re-key auxin_flow_bias from child_id → child ptr.
    for (const LoadedRecord& rec : records) {
        Node* self = by_id.at(rec.common.id);
        for (const auto& kv : rec.common.auxin_flow_bias) {
            auto cit = by_id.find(kv.first);
            if (cit != by_id.end()) self->auxin_flow_bias[cit->second] = kv.second;
            // else: child was removed between save and load — skip.
        }
    }

    plant->set_next_id(next_id);

    LoadedPlant out;
    out.plant = std::move(plant);
    out.genome = active;
    out.engine_tick = engine_tick;
    return out;
}
```

Note the bias map reconstruction: the key is the child `Node*`, value is the bias. The saved map has `child_id → bias`. Look up `child_id` in `by_id` to get the pointer.

- [ ] **Step 8.3: Build + test**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[plant_snapshot][roundtrip]" 2>&1 | tail -5
```
Expected: round-trip identity passes.

- [ ] **Step 8.4: Full-suite regression**

```
./build/botany_tests 2>&1 | tail -3
```
Expected: all tests pass.

- [ ] **Step 8.5: Commit**

```bash
git add src/serialization/plant_snapshot.cpp tests/test_plant_snapshot.cpp
git commit -m "plant_snapshot: three-pass load reconstructs tree + canalization"
```

---

## Task 9: Continuation-equivalence test

**Files:**
- Test: `tests/test_plant_snapshot.cpp`.

- [ ] **Step 9.1: Write the test**

Append:
```cpp
TEST_CASE("continuation equivalence: load + run vs direct-run match", "[plant_snapshot][continuation]") {
    // A: grow 30 ticks
    Engine engineA;
    engineA.create_plant(default_genome(), glm::vec3(0.0f));
    for (int i = 0; i < 30; i++) engineA.tick();

    // Save A at tick 30, continue to tick 60.
    auto tmp = std::filesystem::temp_directory_path() / "botany_cont_test";
    std::filesystem::remove_all(tmp);
    SaveResult sr = save_plant_snapshot(engineA.get_plant(0), engineA.get_tick(), tmp.string());
    REQUIRE(sr.ok);
    for (int i = 0; i < 30; i++) engineA.tick();

    // B: load A's snapshot at tick 30 and continue for 30 ticks.
    Engine engineB;
    LoadedPlant lp = load_plant_snapshot(sr.path, std::nullopt);
    engineB.adopt_plant(std::move(lp.plant));
    engineB.set_tick(static_cast<uint32_t>(lp.engine_tick));
    for (int i = 0; i < 30; i++) engineB.tick();

    // Compare terminal state.
    REQUIRE(engineA.get_plant(0).node_count() == engineB.get_plant(0).node_count());
    std::unordered_map<uint32_t, const Node*> a_by_id;
    engineA.get_plant(0).for_each_node([&](const Node& n) { a_by_id[n.id] = &n; });
    engineB.get_plant(0).for_each_node([&](const Node& n) {
        auto it = a_by_id.find(n.id);
        REQUIRE(it != a_by_id.end());
        deep_compare_nodes(*it->second, n);
    });

    std::filesystem::remove_all(tmp);
}
```

This test depends on Task 10 (Engine::adopt_plant / set_tick). It is written now, and Task 10 makes it pass.

- [ ] **Step 9.2: Verify it fails to compile**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5
```
Expected: error about `adopt_plant` or `set_tick` not being a member of Engine.

Move on to Task 10 without committing.

---

## Task 10: Engine API — `adopt_plant` + `set_tick`

**Files:**
- Modify: `src/engine/engine.h`, `src/engine/engine.cpp`.

- [ ] **Step 10.1: Add declarations to `engine.h`**

In the `public:` section:
```cpp
    // Install a plant constructed outside the engine (e.g., from a snapshot).
    // Returns the assigned PlantID.  Takes ownership.
    PlantID adopt_plant(std::unique_ptr<Plant> plant);

    // Override the current tick counter (used when loading a snapshot so HUDs
    // and age-dependent signals see the correct elapsed sim time).
    void set_tick(uint32_t tick);
```

- [ ] **Step 10.2: Implement in `engine.cpp`**

Append:
```cpp
PlantID Engine::adopt_plant(std::unique_ptr<Plant> plant) {
    PlantID id = static_cast<PlantID>(plants_.size());
    plants_.push_back(std::move(plant));
    return id;
}

void Engine::set_tick(uint32_t tick) {
    tick_ = tick;
}
```

- [ ] **Step 10.3: Build + run continuation test**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[plant_snapshot][continuation]" 2>&1 | tail -5
```
Expected: continuation-equivalence passes.

- [ ] **Step 10.4: Full-suite regression**

```
./build/botany_tests 2>&1 | tail -3
```
Expected: all tests pass.

- [ ] **Step 10.5: Commit**

```bash
git add src/engine/engine.h src/engine/engine.cpp tests/test_plant_snapshot.cpp
git commit -m "engine: adopt_plant + set_tick for snapshot-driven plant install"
```

---

## Task 11: Genome override + canalization + corrupt-file tests

**Files:**
- Test: `tests/test_plant_snapshot.cpp`.

- [ ] **Step 11.1: Write three tests**

Append:
```cpp
TEST_CASE("genome override replaces the embedded genome on load", "[plant_snapshot][override]") {
    Genome original = default_genome();
    original.growth_rate = 0.9999f;  // distinctive
    Plant plant(original, glm::vec3(0.0f));

    auto tmp = std::filesystem::temp_directory_path() / "botany_override_test";
    std::filesystem::remove_all(tmp);
    SaveResult sr = save_plant_snapshot(plant, 5, tmp.string());
    REQUIRE(sr.ok);

    Genome override_g = default_genome();
    override_g.growth_rate = 0.00001f; // different, distinctive

    LoadedPlant lp = load_plant_snapshot(sr.path, override_g);
    REQUIRE(lp.genome.growth_rate == 0.00001f);
    REQUIRE(lp.plant->genome().growth_rate == 0.00001f);

    std::filesystem::remove_all(tmp);
}

TEST_CASE("load rejects bad magic", "[plant_snapshot][corrupt]") {
    auto tmp = std::filesystem::temp_directory_path() / "botany_bad_magic.tree";
    { std::ofstream f(tmp, std::ios::binary); f.write("XXXX", 4); f.write("\x00\x00\x00\x00", 4); }
    REQUIRE_THROWS_AS(load_plant_snapshot(tmp.string(), std::nullopt), std::runtime_error);
    std::filesystem::remove(tmp);
}

TEST_CASE("load rejects unknown version", "[plant_snapshot][corrupt]") {
    auto tmp = std::filesystem::temp_directory_path() / "botany_bad_version.tree";
    { std::ofstream f(tmp, std::ios::binary);
      f.write("BTNT", 4);
      uint32_t bogus = 999999u;
      f.write(reinterpret_cast<const char*>(&bogus), sizeof(bogus));
    }
    REQUIRE_THROWS_AS(load_plant_snapshot(tmp.string(), std::nullopt), std::runtime_error);
    std::filesystem::remove(tmp);
}

TEST_CASE("load rejects truncated body", "[plant_snapshot][corrupt]") {
    // Save a valid file, then truncate it to just the header.
    Plant plant(default_genome(), glm::vec3(0.0f));
    auto tmp = std::filesystem::temp_directory_path() / "botany_trunc_test";
    std::filesystem::remove_all(tmp);
    SaveResult sr = save_plant_snapshot(plant, 0, tmp.string());
    REQUIRE(sr.ok);
    std::filesystem::resize_file(sr.path, 64);
    REQUIRE_THROWS_AS(load_plant_snapshot(sr.path, std::nullopt), std::runtime_error);
    std::filesystem::remove_all(tmp);
}

TEST_CASE("canalization bias preserved across save+load", "[plant_snapshot][canalization]") {
    // Construct a small tree by hand with deterministic bias values (avoids
    // relying on multi-tick sim dynamics for test stability).
    Engine engine;
    engine.create_plant(default_genome(), glm::vec3(0.0f));
    for (int i = 0; i < 40; i++) engine.tick();
    Plant& pl = engine.get_plant_mut(0);

    // Inject known bias on the seed so we have a non-empty map guaranteed.
    Node* seed = pl.seed_mut();
    REQUIRE(!seed->children.empty());
    Node* child0 = seed->children[0];
    seed->auxin_flow_bias[child0] = 0.73f;

    auto tmp = std::filesystem::temp_directory_path() / "botany_cana_test";
    std::filesystem::remove_all(tmp);
    SaveResult sr = save_plant_snapshot(pl, engine.get_tick(), tmp.string());
    REQUIRE(sr.ok);

    LoadedPlant lp = load_plant_snapshot(sr.path, std::nullopt);
    Node* loaded_seed = lp.plant->seed_mut();
    // Find the loaded-side child that has the same id as child0.
    Node* loaded_child0 = nullptr;
    for (Node* c : loaded_seed->children) {
        if (c->id == child0->id) { loaded_child0 = c; break; }
    }
    REQUIRE(loaded_child0 != nullptr);
    auto it = loaded_seed->auxin_flow_bias.find(loaded_child0);
    REQUIRE(it != loaded_seed->auxin_flow_bias.end());
    REQUIRE(it->second == 0.73f);

    std::filesystem::remove_all(tmp);
}
```

- [ ] **Step 11.2: Build + run**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5 && ./build/botany_tests "[plant_snapshot]" 2>&1 | tail -5
```
Expected: all plant_snapshot tests pass (10 total).

- [ ] **Step 11.3: Full-suite regression**

```
./build/botany_tests 2>&1 | tail -3
```
Expected: all tests pass.

- [ ] **Step 11.4: Commit**

```bash
git add tests/test_plant_snapshot.cpp
git commit -m "plant_snapshot: genome override + corrupt-file + canalization tests"
```

---

## Task 12: CLI integration in `botany_realtime`

**Files:**
- Modify: `src/app_realtime.cpp`.

- [ ] **Step 12.1: Parse the new flags**

Find the existing CLI/argv handling in `app_realtime.cpp` (near the top of `main`). Add:
```cpp
// Near top of main() before genome initialization:
std::string load_plant_path;
std::string genome_override_path;
for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    if (a == "--load-plant" && i + 1 < argc) load_plant_path = argv[++i];
    else if (a == "--genome-override" && i + 1 < argc) genome_override_path = argv[++i];
    // ... keep any existing flag handling that's already there
}
```

- [ ] **Step 12.2: Wire load into engine bootstrap**

Where the engine currently calls `engine.create_plant(g, ...)` (existing "fresh plant" path), wrap it so that a `--load-plant` path takes precedence:
```cpp
#include "serialization/plant_snapshot.h"
#include "engine/genome_io.h"

// ...
if (!load_plant_path.empty()) {
    std::optional<Genome> override;
    if (!genome_override_path.empty()) {
        override = load_genome_file(genome_override_path);
    }
    try {
        LoadedPlant lp = load_plant_snapshot(load_plant_path, override);
        engine.adopt_plant(std::move(lp.plant));
        engine.set_tick(static_cast<uint32_t>(lp.engine_tick));
        std::cerr << "Loaded plant from " << load_plant_path
                  << " at tick " << lp.engine_tick << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Failed to load snapshot: " << e.what() << "\n";
        return 1;
    }
} else {
    engine.create_plant(g, glm::vec3(0.0f));
    if (!genome_override_path.empty()) {
        std::cerr << "--genome-override ignored without --load-plant\n";
    }
}
```

- [ ] **Step 12.3: Build + sanity smoke**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5
```
Expected: clean build.

Manual sanity (no test, this is a viewer path):
```
./build/botany_realtime --load-plant /nonexistent.tree
```
Expected: exits with non-zero and prints "Failed to load snapshot: plant_snapshot: cannot open /nonexistent.tree".

- [ ] **Step 12.4: Full-suite regression**

```
./build/botany_tests 2>&1 | tail -3
```
Expected: all tests pass.

- [ ] **Step 12.5: Commit**

```bash
git add src/app_realtime.cpp
git commit -m "app_realtime: --load-plant and --genome-override CLI flags"
```

---

## Task 13: `Cmd/Ctrl+S` save hotkey + ImGui toast

**Files:**
- Modify: `src/app_realtime.cpp`.

- [ ] **Step 13.1: Add state for the save toast**

Near the top of `main()` (alongside other viewer state):
```cpp
// Save-toast state.
std::string save_toast_msg;
double save_toast_expires = 0.0;
```

- [ ] **Step 13.2: Detect Cmd/Ctrl+S once per press**

In the main loop, next to the existing `glfwGetKey(window, GLFW_KEY_...)` block:
```cpp
// Cmd+S (macOS) / Ctrl+S (Linux/Win) saves the primary plant.
static bool s_prev = false;
bool s_down = glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS;
bool mod    = glfwGetKey(window, GLFW_KEY_LEFT_SUPER)   == GLFW_PRESS
           || glfwGetKey(window, GLFW_KEY_RIGHT_SUPER)  == GLFW_PRESS
           || glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS
           || glfwGetKey(window, GLFW_KEY_RIGHT_CONTROL)== GLFW_PRESS;
if (s_down && mod && !s_prev) {
    if (engine.plant_count() > 0) {
        SaveResult r = save_plant_snapshot(engine.get_plant(0), engine.get_tick());
        save_toast_msg    = r.ok ? ("Saved \u2192 " + r.path) : ("Save failed: " + r.error);
        save_toast_expires = glfwGetTime() + 2.0;
    }
}
s_prev = s_down;
```

- [ ] **Step 13.3: Render the toast in the Controls panel**

Inside the existing ImGui Controls panel block, after other `ImGui::Text(...)` lines:
```cpp
if (glfwGetTime() < save_toast_expires && !save_toast_msg.empty()) {
    ImGui::Separator();
    ImGui::TextWrapped("%s", save_toast_msg.c_str());
}
```

- [ ] **Step 13.4: Build + manual smoke**

```
/usr/local/bin/cmake --build build 2>&1 | tail -5
```
Expected: clean build.

Manual: launch `./build/botany_realtime`, wait for plant to spawn, press Cmd+S (or Ctrl+S on Linux). Expected: toast appears for 2 s, `saves/` directory created in project root, file `plant_<timestamp>_tick<N>.tree` exists.

- [ ] **Step 13.5: Commit**

```bash
git add src/app_realtime.cpp
git commit -m "app_realtime: Cmd/Ctrl+S saves a .tree snapshot with ImGui toast"
```

---

## Task 14: Manual end-to-end smoke + CLAUDE.md update

**Files:**
- Modify: `CLAUDE.md`.

- [ ] **Step 14.1: Manual smoke walk-through**

1. `./build/botany_realtime` — grow a plant for a couple of minutes of real time.
2. Press Cmd/Ctrl+S, note the path in the toast.
3. Close the app.
4. Relaunch: `./build/botany_realtime --load-plant <path_from_step_2>`.
5. Verify the tree appears identical to what was saved.
6. Press play, verify growth continues without visual glitches.

If any step misbehaves, root-cause and fix before moving on.

- [ ] **Step 14.2: Update CLAUDE.md**

In the `## Build & Run` section, add an entry:
```
# Load a previously saved snapshot and continue simulating it
./build/botany_realtime --load-plant saves/plant_<timestamp>_tick<N>.tree [--genome-override some_genome.txt]

# Save the current plant from inside botany_realtime
Cmd+S (macOS) / Ctrl+S (Linux/Windows) writes saves/plant_<timestamp>_tick<N>.tree
```

In the `### Apps` section (just after `app_sugar_test.cpp` description), add:
```
Snapshot format: `src/serialization/plant_snapshot.{h,cpp}` writes `.tree`
binary files (magic "BTNT", single plant, full node + conduit + meristem
state). Playback recordings (`botany_headless --recording-interval N` → `.bin`)
are a separate, playback-only format owned by `serializer.*`.
```

- [ ] **Step 14.3: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: snapshot save/load in CLAUDE.md"
```

---

## Task 15: Push to origin

- [ ] **Step 15.1: Verify tree health**

```
git status
git log --oneline -15
./build/botany_tests 2>&1 | tail -3
```

- [ ] **Step 15.2: Push**

```
git push origin main
```
