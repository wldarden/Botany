// src/engine/perf_log.h
// Lightweight performance logger — tracks time spent in each phase of the tick.
#pragma once

#include <chrono>
#include <fstream>
#include <string>
#include <cstdint>

namespace botany {

struct PerfStats {
    // Top-level engine buckets.
    double light_ms = 0;        // CPU light/shadow pass (0 when GPU LightSystem drives it, or on skipped intervals)
    double plant_tick_ms = 0;   // total across all plants; phase0+phase1+phase2 should sum to ~this
    double debug_log_ms = 0;    // --debug-log CSV write (not plant-sim)

    // Plant::tick_tree phase breakdown (all summed across all plants in a tick).
    double phase0_auxin_transport_ms = 0;  // diffuse_auxin_across_seed_junction + pin_transport
    double phase1_metabolism_ms      = 0;  // per-node DFS metabolism walk
    double phase2_vascular_ms        = 0;  // full vascular_sub_stepped pass

    // Vascular sub-phase breakdown (summed across N sub-steps × all plants).
    // These should sum to ~phase2_vascular_ms (remainder is loop/budget overhead).
    double vascular_inject_ms  = 0;
    double vascular_jacobi_ms  = 0;
    double vascular_radial_ms  = 0;
    double vascular_extract_ms = 0;

    // Per-tick sugar accounting (g glucose, summed across all nodes in all plants)
    float sugar_spent_maintenance = 0;  // consumed by maintenance this tick
    float sugar_spent_growth      = 0;  // consumed by growth / phototropism this tick
    float sugar_spent_transport   = 0;  // moved between nodes this tick (= total exported = total imported)

    uint32_t node_count = 0;

    void reset() { *this = PerfStats{}; }
};

class PerfLog {
public:
    void open(const std::string& path) {
        file_.open(path);
        if (file_) {
            file_ << "tick,nodes,total_ms,light_ms,plant_ms,debug_ms,"
                  << "phase0_auxin_ms,phase1_metabolism_ms,phase2_vascular_ms,"
                  << "vascular_inject_ms,vascular_jacobi_ms,vascular_radial_ms,vascular_extract_ms,"
                  << "sugar_spent_maintenance,sugar_spent_growth,sugar_spent_transport\n";
        }
    }

    bool is_open() const { return file_.is_open(); }
    void close() { file_.close(); }

    PerfStats& stats() { return stats_; }

    void flush(uint32_t tick) {
        if (!file_) return;
        double total = stats_.light_ms + stats_.plant_tick_ms + stats_.debug_log_ms;
        file_ << tick << ","
              << stats_.node_count << ","
              << total << ","
              << stats_.light_ms << ","
              << stats_.plant_tick_ms << ","
              << stats_.debug_log_ms << ","
              << stats_.phase0_auxin_transport_ms << ","
              << stats_.phase1_metabolism_ms << ","
              << stats_.phase2_vascular_ms << ","
              << stats_.vascular_inject_ms << ","
              << stats_.vascular_jacobi_ms << ","
              << stats_.vascular_radial_ms << ","
              << stats_.vascular_extract_ms << ","
              << stats_.sugar_spent_maintenance << ","
              << stats_.sugar_spent_growth << ","
              << stats_.sugar_spent_transport << "\n";
        stats_.reset();
    }

private:
    std::ofstream file_;
    PerfStats stats_;
};

// Scoped timer — adds elapsed time to a double on destruction
struct ScopedTimer {
    double& target;
    std::chrono::high_resolution_clock::time_point start;
    static inline double dummy = 0.0;  // sink for when perf logging is off

    explicit ScopedTimer(double& t) : target(t), start(std::chrono::high_resolution_clock::now()) {}
    ~ScopedTimer() {
        auto end = std::chrono::high_resolution_clock::now();
        target += std::chrono::duration<double, std::milli>(end - start).count();
    }
};

} // namespace botany
