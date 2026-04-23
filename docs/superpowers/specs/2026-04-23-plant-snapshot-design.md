# Plant Snapshot Save/Load — Design

**Date:** 2026-04-23
**Status:** Draft (approved in brainstorming)

## Problem

After 36k ticks of iteration, the sim produces a recognizable open-crown tree. There is no way to save that tree, close the app, reopen later, and keep growing from the same state. The existing `src/serialization/serializer.*` captures playback-only snapshots (position, radius, 3 chemicals) intended for `botany_playback` replay, which is structurally unsuitable for resuming live simulation.

## Goal

Snapshot a live plant to disk such that loading it into `botany_realtime` restores full internal state and the next tick is indistinguishable from what would have been tick N+1 in the original run.

## Non-goals

- Multi-plant / whole-engine snapshots (single plant only in v1)
- ImGui file browser / drag-drop loading (CLI + hotkey only in v1)
- Automatic schema migration between file format versions
- Evolve-app integration (treating snapshots as evolution seeds)

## Decisions locked in brainstorming

1. **Use case:** live continuation, not visualization-only.
2. **Genome handling:** baked into the file with an optional `--genome-override <path.txt>` CLI flag for experiment-on-old-trees workflows.
3. **UX:** `--load-plant <path>` on startup; `Cmd/Ctrl+S` in viewer saves to `saves/plant_<timestamp>_tick<N>.tree`.
4. **Format:** new `.tree` binary file, deliberately separate from the playback `.bin`.
5. **Scope:** single plant per file.

## File format

Little-endian binary. All integers fixed-width. Floats IEEE-754 32-bit.

### Header

```
magic        : char[4]       = "BTNT"
version      : uint32        (bumped on any schema change; no migration)
engine_tick  : uint64        (elapsed ticks in the saved world)
genome       : Genome struct (binary dump, size tied to version)
node_count   : uint32
next_node_id : uint32        (so post-load spawns don't collide with existing IDs)
```

### Per-node record (repeated node_count times, DFS order so parent precedes children)

**Common header for every node type:**
- `id: uint32`
- `parent_id: uint32` (`UINT32_MAX` for seed)
- `type: uint8` (NodeType enum)
- `age: uint32`
- `starvation_ticks: uint32`
- `dormant_ticks: uint32`
- `position: vec3`
- `offset: vec3`
- `rest_offset: vec3`
- `radius: float32`
- `ever_active: uint8`
- `local_chemicals_count: uint16` then that many `(ChemicalID: uint8, value: float32)` pairs
- `auxin_flow_bias_count: uint16` then that many `(child_id: uint32, value: float32)` pairs

**If type is STEM or ROOT (conduit nodes), additionally:**
- `phloem_chemicals_count: uint16` then pairs
- `xylem_chemicals_count: uint16` then pairs

**Type-specific trailer:**
- **LEAF:** `leaf_size: float32`, `light_exposure: float32`, `senescence_ticks: uint32`, `deficit_ticks: uint32`, `facing: vec3`
- **APICAL:** `active: uint8`, `is_primary: uint8`, `growth_dir: vec3`, `ticks_since_last_node: uint32`
- **ROOT_APICAL:** same as APICAL + `internodes_spawned: uint32`
- **STEM / ROOT:** no extra trailer beyond conduit pools (covered above)

### Fields intentionally not persisted

Regenerated on first post-load tick:
- `total_mass`, `mass_moment`, `stress` — recomputed in `update_physics`
- All `tick_chem_*`, `tick_edge_*`, `tick_sugar_*`, `tick_auxin_produced`, `tick_cytokinin_produced`, `tick_sugar_transport` — per-tick flux counters, zeroed each tick anyway
- `transport_received`, `seed_collected`, `last_auxin_flux` — transient one-hop buffers

## Modules

### `src/serialization/plant_snapshot.h` / `.cpp`

Public API:
```cpp
// Save
struct SaveResult { bool ok; std::string path; std::string error; };
SaveResult save_plant_snapshot(const Plant& p, uint64_t engine_tick, const std::string& dir = "saves");

// Load
struct LoadedPlant {
    std::unique_ptr<Plant> plant;
    Genome genome;
    uint64_t engine_tick;
};
LoadedPlant load_plant_snapshot(const std::string& path, const std::optional<Genome>& genome_override);
```

Deliberately a separate translation unit from `serializer.cpp`. Playback and snapshot serve different purposes and shouldn't co-evolve.

### `src/app_realtime.cpp` changes

- Parse two new CLI flags: `--load-plant <path>`, `--genome-override <path.txt>`.
- On startup:
  - If `--load-plant` given: call `load_plant_snapshot`. Install the resulting plant + tick count into the Engine instead of spawning a fresh default plant.
  - If `--genome-override` given without `--load-plant`: warn and ignore (it only makes sense with load).
- In the GLFW key handler: on `S + Cmd/Ctrl`, call `save_plant_snapshot`. Store the `SaveResult` in a `save_toast_t` struct `{ message, expires_at }` and display it in the Controls panel while `now < expires_at`.

### `Plant` and `Engine` API additions

The current `Plant` constructor creates the seed + primary shoot/root apicals automatically. To install a loaded plant we need to construct one without that default initialization. Two equally clean options; implementation plan will pick one:

- **A.** Add `Plant::from_snapshot(Genome, std::vector<std::unique_ptr<Node>>&&, uint32_t next_id)` factory that skips the default-seed path.
- **B.** Add a private `Plant::Plant(Genome, empty_tag)` constructor + a friend `load_plant_snapshot` that populates `nodes_` and `id_counter_` directly.

Either way `Engine` needs a way to install an already-constructed plant:
- Add `Engine::adopt_plant(std::unique_ptr<Plant>, uint64_t tick_count)`.
- Alternatively, `Engine::set_tick_count(uint64_t)` + `Engine::add_plant(std::unique_ptr<Plant>)` as separate calls if that's more consistent with existing patterns.

The plan will pick the combination most consistent with the current `Engine` / `Plant` API.

### Reuse of existing code

- Text genome parsing: `load_genome_file` in `app_realtime.cpp:34` already parses the `.env`-style format. Hoist it into a small shared header `src/engine/genome_io.h` / `.cpp` so `plant_snapshot.cpp` can call it too.

## Load algorithm (detail)

1. Open file, verify `"BTNT"` magic. On mismatch: error, return failure.
2. Read `version`. If unrecognized: error.
3. Read `engine_tick`, `genome`, `node_count`, `next_node_id`.
4. If `genome_override` supplied, replace `genome` with it.
5. Construct an empty `Plant` with the chosen genome (no default seed — we will populate it).
6. **Pass 1 — allocate:** for each node record, allocate the correct subclass (`StemNode`, `RootNode`, `LeafNode`, `ApicalNode`, `RootApicalNode`), populate common fields + type-specific trailer + chemicals. Leave `parent = nullptr` and `children = {}`. Add to an `id → Node*` map and to the plant's node vector.
7. **Pass 2 — hook up tree:** for each node with a non-UINT32_MAX parent_id, look up the parent in the id map, set `node.parent`, append to `parent->children`.
8. **Pass 3 — resolve bias maps:** for each node, convert the saved `{child_id → bias}` list back to `{child_ptr → bias}` using the id map.
9. Set `plant.next_node_id()` from the saved header.
10. Return `LoadedPlant { plant, genome, engine_tick }` to the caller.

Failure at any step: return a meaningful error string; caller decides whether to fall back to default or abort.

## Save algorithm (detail)

1. Create `saves/` if it doesn't exist. On failure: error toast.
2. Format filename: `saves/plant_YYYYMMDD_HHMMSS_tick<N>.tree` from current local time.
3. Enumerate nodes via DFS from seed so parents appear before children.
4. For each node, translate `auxin_flow_bias` (`child_ptr → bias`) to `(child_id → bias)` using `child->id`.
5. Write header + node records via a small `BinaryWriter` helper.
6. Close + flush. Return the path.

## Testing strategy

### `tests/test_plant_snapshot.cpp`

1. **Round-trip identity.** Grow a small plant for 100 ticks with a fixed genome. Save. Load into a fresh Engine. Deep-compare every node: ids match, types match, every `local_chemical` value within 1e-6, every `auxin_flow_bias` entry within 1e-6, all conduit pools, all meristem flags. Parent pointers match. Children vectors match.
2. **Continuation equivalence.** Grow plant A for 100 ticks. Save. Continue A to tick 200. Load snapshot into plant B and run 100 more ticks. A@200 and B@200 are bit-identical: same node count, same IDs, same geometry, same chemical values.
3. **Genome override.** Save under genome_X (tweak a distinctive field like `growth_rate`). Load with override genome_Y. Verify the running plant's genome reports Y's value, not X's.
4. **Corrupt-file rejection.** Separate cases: bad magic, unknown version, truncated body, `parent_id` pointing at a nonexistent node, a non-seed with UINT32_MAX parent. Each should fail with a specific error string (no crashes).
5. **Canalization preserved.** Grow a plant long enough that the trunk has nonzero bias to at least 3 children. Save, load, assert the bias map on the loaded trunk has the same entries.

### Manual smoke

Grow a visible tree in `botany_realtime`. Press Cmd+S. Quit. Relaunch with `--load-plant <saved_path>`. Tree looks identical. Press play. Growth continues organically.

## Open risks

- **Schema drift.** Every new field we add on `Node`, `StemNode`, etc. silently invalidates old snapshot files. Mitigation: the `version` counter. If schema drift happens a lot, a later migration layer may be worth building. Not in v1.
- **Genome struct padding differences across platforms.** If `Genome` has implicit padding, a file saved on one build won't load on another. Mitigation: **do not** dump the struct raw. Serialize field-by-field in a stable declared order. The implementation plan will include a `write_genome(writer, g)` / `read_genome(reader)` pair that lists every field explicitly — the same way `evolve`'s `.txt` genome format lists them, just in binary. Adding a genome field means adding a line to both functions plus bumping the file `version`.
- **The auxin_flow_bias pointer↔id round-trip relies on IDs being unique.** If IDs ever get reused or wrap, this breaks. Mitigation: IDs are currently 32-bit monotonic from a plant-level counter; reasonable for a plant that never exceeds ~4B nodes over its lifetime. Not a v1 concern.
