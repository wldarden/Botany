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

    if (argc >= 2) num_ticks = std::atoi(argv[1]);
    if (argc >= 3) output_path = argv[2];

    std::cout << "Running " << num_ticks << " ticks, saving to " << output_path << std::endl;

    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));

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
