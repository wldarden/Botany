// src/app_headless.cpp
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "engine/engine.h"
#include "serialization/serializer.h"

using namespace botany;

int main(int argc, char* argv[]) {
    int num_ticks = 200;
    std::string output_path = "recording.bin";
    std::string debug_log_path;    // optional; empty = no debug log
    std::string economy_log_path;  // optional; empty = no global economy log

    for (int a = 1; a < argc; ++a) {
        std::string arg = argv[a];
        if (arg == "--debug-log" && a + 1 < argc) {
            debug_log_path = argv[++a];
        } else if (arg == "--economy-log" && a + 1 < argc) {
            economy_log_path = argv[++a];
        } else if (num_ticks == 200 && arg[0] != '-') {
            num_ticks = std::atoi(arg.c_str());
        } else if (output_path == "recording.bin" && arg[0] != '-') {
            output_path = arg;
        }
    }

    std::cout << "Running " << num_ticks << " ticks, saving to " << output_path << std::endl;
    if (!debug_log_path.empty()) {
        std::cout << "Debug log: " << debug_log_path << std::endl;
    }

    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));
    if (!debug_log_path.empty()) {
        engine.debug_log().open(debug_log_path);
    }
    if (!economy_log_path.empty()) {
        engine.global_economy_log().open(economy_log_path);
        std::cout << "Global economy log: " << economy_log_path << std::endl;
    }

    std::ofstream file(output_path, std::ios::binary);
    if (!file) {
        std::cerr << "Failed to open " << output_path << std::endl;
        return 1;
    }

    save_recording_header(file, g, static_cast<uint32_t>(num_ticks));

    for (int i = 0; i < num_ticks; i++) {
        engine.tick();
        save_tick(file, engine, 0);

        if ((i + 1) % 50 == 0) {
            std::cout << "Tick " << (i + 1) << "/" << num_ticks
                      << " (" << engine.get_plant(0).node_count() << " nodes)" << std::endl;
        }
    }

    std::cout << "Done. Final node count: " << engine.get_plant(0).node_count() << std::endl;
    return 0;
}
