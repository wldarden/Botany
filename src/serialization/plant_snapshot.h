#pragma once

#include <cstdint>
#include <iosfwd>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <glm/vec3.hpp>
#include "engine/genome.h"
#include "engine/node/node.h"

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

void             write_node_common(std::ostream& out, const Node& n, uint32_t parent_id);
NodeCommonRecord read_node_common(std::istream& in);

// --- Per-type node trailer records ---

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

void        write_leaf_trailer(std::ostream& out, const LeafNode& l);
LeafTrailer read_leaf_trailer(std::istream& in);

void          write_apical_trailer(std::ostream& out, const ApicalNode& a);
ApicalTrailer read_apical_trailer(std::istream& in);

void              write_root_apical_trailer(std::ostream& out, const RootApicalNode& r);
RootApicalTrailer read_root_apical_trailer(std::istream& in);

void         write_conduit_pools(std::ostream& out, const Node& n);
ConduitPools read_conduit_pools(std::istream& in);

// --- Genome binary I/O (internal; exposed for tests) ---
// Writes every Genome field in declared order.  Adding or reordering a
// Genome field requires editing these two functions together AND bumping
// PLANT_SNAPSHOT_VERSION.
void   write_genome_binary(std::ostream& out, const Genome& g);
Genome read_genome_binary(std::istream& in);

} // namespace botany
