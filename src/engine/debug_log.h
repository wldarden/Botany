#pragma once

#include <fstream>
#include <string>

namespace botany {

class Plant;
struct WorldParams;

class DebugLog {
public:
    void open(const std::string& path);
    void close();
    bool is_open() const { return file_.is_open(); }

    // Log one line per leaf per tick: sugar balance, production, maintenance, exposure, etc.
    void log_tick(uint32_t tick, const Plant& plant, const WorldParams& world);

private:
    std::ofstream file_;
    bool header_written_ = false;
};

} // namespace botany
