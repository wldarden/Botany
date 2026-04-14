// src/engine/perf_log.h
// Lightweight performance logger — tracks time spent in each phase of the tick.
#pragma once

#include <chrono>
#include <fstream>
#include <string>
#include <cstdint>

namespace botany {

struct PerfStats {
    double light_ms = 0;
    double plant_tick_ms = 0;
    double debug_log_ms = 0;

    // Node-level breakdown (accumulated across all nodes in a tick)
    double node_position_ms = 0;
    double node_maintenance_ms = 0;
    double node_grow_ms = 0;
    double node_mass_stress_ms = 0;
    double node_droop_break_ms = 0;
    double node_transport_ms = 0;

    uint32_t node_count = 0;
    uint32_t shoot_apical_count = 0;

    void reset() { *this = PerfStats{}; }
};

class PerfLog {
public:
    void open(const std::string& path) {
        file_.open(path);
        if (file_) {
            file_ << "tick,nodes,SA,total_ms,light_ms,plant_ms,debug_ms,"
                  << "position_ms,maintenance_ms,grow_ms,mass_stress_ms,droop_break_ms,transport_ms\n";
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
              << stats_.shoot_apical_count << ","
              << total << ","
              << stats_.light_ms << ","
              << stats_.plant_tick_ms << ","
              << stats_.debug_log_ms << ","
              << stats_.node_position_ms << ","
              << stats_.node_maintenance_ms << ","
              << stats_.node_grow_ms << ","
              << stats_.node_mass_stress_ms << ","
              << stats_.node_droop_break_ms << ","
              << stats_.node_transport_ms << "\n";
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
