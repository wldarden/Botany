#include "engine/world_params.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

namespace botany {

WorldParams load_world_params(const std::string& path) {
    WorldParams wp = default_world_params();
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Warning: could not open " << path << ", using defaults" << std::endl;
        return wp;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Simple key-value parsing: look for "key": value
        if (line.find("light_level") != std::string::npos) {
            auto colon = line.find(':');
            if (colon != std::string::npos) {
                wp.light_level = std::stof(line.substr(colon + 1));
            }
        } else if (line.find("starvation_ticks_max") != std::string::npos) {
            auto colon = line.find(':');
            if (colon != std::string::npos) {
                wp.starvation_ticks_max = static_cast<uint32_t>(std::stoi(line.substr(colon + 1)));
            }
        }
    }
    return wp;
}

} // namespace botany
