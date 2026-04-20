#pragma once

#include <fstream>
#include <string>
#include <unordered_map>
#include "engine/chemical/chemical.h"

namespace botany {

class Plant;
struct WorldParams;

// Per-tick global economy tracker.  For each chemical of interest, records:
//   storage  — total across all pools (local_env + phloem + xylem) at tick end
//   produced — total created this tick (photosynthesis, absorption,
//              hormone synthesis, etc.) from known per-node trackers
//   consumed — total destroyed this tick (maintenance, decay, transpiration,
//              growth use) — partly measured, partly derived from mass balance
//   delta    — storage_now - storage_prev (positive = net gain, negative = loss)
//
// One row per (tick, chemical) pair.  Designed to make it obvious when a
// chemical is being massively over-produced (storage ballooning) or
// under-supplied (storage trending to zero).
class GlobalEconomyLog {
public:
    void open(const std::string& path);
    void close();
    bool is_open() const { return file_.is_open(); }

    void log_tick(uint32_t tick, const Plant& plant, const WorldParams& world);

private:
    std::ofstream file_;
    bool header_written_ = false;
    // Last-tick storage per chemical, for delta computation.
    std::unordered_map<ChemicalID, float> prev_storage_;
};

} // namespace botany
