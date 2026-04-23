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

// --- Genome binary I/O (internal; exposed for tests) ---
// Writes every Genome field in declared order.  Adding or reordering a
// Genome field requires editing these two functions together AND bumping
// PLANT_SNAPSHOT_VERSION.
void   write_genome_binary(std::ostream& out, const Genome& g);
Genome read_genome_binary(std::istream& in);

} // namespace botany
