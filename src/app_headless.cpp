// src/app_headless.cpp
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include "engine/engine.h"
#include "engine/plant.h"
#include "engine/node/node.h"
#include "engine/chemical/chemical.h"
#include "serialization/serializer.h"

using namespace botany;

namespace {

// Walk from `start` along a strict parent→child chain, greedily picking the
// child whose subtree contains the largest number of nodes of `accepted_type`
// (or of `accepted_type_alt`).  Ties broken by child->id for determinism.
// Returns the nodes visited, in order from start to tip (start is NOT included).
std::vector<const Node*> walk_longest_chain(
    const Node* start,
    NodeType accepted_type,
    NodeType accepted_type_alt,
    NodeType tip_type)
{
    // Recursive subtree-depth measurement (number of nodes through a given child).
    auto measure = [&](const Node* root, auto&& self) -> int {
        int n = 1;
        int best = 0;
        for (const Node* c : root->children) {
            if (c->type != accepted_type && c->type != accepted_type_alt && c->type != tip_type) continue;
            best = std::max(best, self(c, self));
        }
        return n + best;
    };

    std::vector<const Node*> chain;
    const Node* cur = start;
    while (cur) {
        const Node* next = nullptr;
        int best_depth = -1;
        for (const Node* c : cur->children) {
            if (c->type != accepted_type && c->type != accepted_type_alt && c->type != tip_type) continue;
            int d = measure(c, measure);
            if (d > best_depth || (d == best_depth && next && c->id < next->id)) {
                best_depth = d;
                next = c;
            }
        }
        if (!next) break;
        chain.push_back(next);
        cur = next;
    }
    return chain;
}

void write_chain_profile(const std::string& path, const Plant& plant) {
    std::ofstream csv(path);
    if (!csv.is_open()) {
        std::cerr << "Failed to open chain-profile file " << path << std::endl;
        return;
    }
    csv << "chain_side,depth,node_id,type,y,radius,auxin,cytokinin,sugar,water,"
           "xyl_cytokinin,xyl_water,phloem_sugar\n";

    const Node* seed = plant.seed();
    auto type_str = [](NodeType t) {
        switch (t) {
            case NodeType::STEM: return "STEM";
            case NodeType::ROOT: return "ROOT";
            case NodeType::LEAF: return "LEAF";
            case NodeType::APICAL: return "SA";
            case NodeType::ROOT_APICAL: return "RA";
        }
        return "?";
    };
    auto emit = [&](const char* side, int depth, const Node* n) {
        const auto* xyl = n->xylem();
        const auto* phl = n->phloem();
        csv << side << ',' << depth << ',' << n->id << ',' << type_str(n->type) << ','
            << n->position.y << ',' << n->radius << ','
            << n->local().chemical(ChemicalID::Auxin) << ','
            << n->local().chemical(ChemicalID::Cytokinin) << ','
            << n->local().chemical(ChemicalID::Sugar) << ','
            << n->local().chemical(ChemicalID::Water) << ','
            << (xyl ? xyl->chemical(ChemicalID::Cytokinin) : 0.0f) << ','
            << (xyl ? xyl->chemical(ChemicalID::Water) : 0.0f) << ','
            << (phl ? phl->chemical(ChemicalID::Sugar) : 0.0f) << '\n';
    };

    // Row 0: seed on both sides.
    emit("shoot", 0, seed);
    emit("root",  0, seed);

    auto shoot_chain = walk_longest_chain(seed,
        NodeType::STEM, NodeType::STEM, NodeType::APICAL);
    int d = 1;
    for (const Node* n : shoot_chain) emit("shoot", d++, n);

    auto root_chain = walk_longest_chain(seed,
        NodeType::ROOT, NodeType::ROOT, NodeType::ROOT_APICAL);
    d = 1;
    for (const Node* n : root_chain) emit("root", d++, n);

    std::cerr << "Chain profile written: " << path
              << " (shoot chain " << shoot_chain.size()
              << " nodes, root chain " << root_chain.size() << " nodes)" << std::endl;
}

} // namespace

int main(int argc, char* argv[]) {
    int num_ticks = 200;
    int recording_interval = 1;   // save every Nth tick; 0 disables recording entirely
    std::string output_path = "recording.bin";
    std::string debug_log_path;    // optional; empty = no debug log
    std::string economy_log_path;  // optional; empty = no global economy log
    std::string chain_profile_path; // optional; empty = no chain-profile dump
    bool num_ticks_set = false;
    bool output_path_set = false;

    for (int a = 1; a < argc; ++a) {
        std::string arg = argv[a];
        if (arg == "--debug-log" && a + 1 < argc) {
            debug_log_path = argv[++a];
        } else if (arg == "--economy-log" && a + 1 < argc) {
            economy_log_path = argv[++a];
        } else if (arg == "--chain-profile" && a + 1 < argc) {
            chain_profile_path = argv[++a];
        } else if (arg == "--recording-interval" && a + 1 < argc) {
            recording_interval = std::atoi(argv[++a]);
            if (recording_interval < 0) recording_interval = 0;
        } else if (!num_ticks_set && arg[0] != '-') {
            num_ticks = std::atoi(arg.c_str());
            num_ticks_set = true;
        } else if (!output_path_set && arg[0] != '-') {
            output_path = arg;
            output_path_set = true;
        }
    }

    // Count how many ticks will actually be saved: ticks at indices
    // 0, interval, 2*interval, ... up to < num_ticks.  This matches the
    // tick-loop condition `i % interval == 0` below.
    uint32_t saved_tick_count = 0;
    if (recording_interval > 0) {
        saved_tick_count = (num_ticks + recording_interval - 1) / recording_interval;
    }

    if (recording_interval == 0) {
        std::cout << "Running " << num_ticks << " ticks, recording DISABLED (--recording-interval 0)" << std::endl;
    } else if (recording_interval == 1) {
        std::cout << "Running " << num_ticks << " ticks, saving every tick to " << output_path << std::endl;
    } else {
        std::cout << "Running " << num_ticks << " ticks, saving every "
                  << recording_interval << " ticks (" << saved_tick_count
                  << " snapshots total) to " << output_path << std::endl;
    }
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

    // Recording file: opened only if recording_interval > 0.  Header records
    // the *actual* number of saved ticks so playback can consume it directly.
    // Sparse recording (interval > 1) keeps disk usage flat for long sims
    // while still allowing periodic inspection of trajectory.
    std::ofstream file;
    if (recording_interval > 0) {
        file.open(output_path, std::ios::binary);
        if (!file) {
            std::cerr << "Failed to open " << output_path << std::endl;
            return 1;
        }
        save_recording_header(file, g, saved_tick_count);
    }

    for (int i = 0; i < num_ticks; i++) {
        engine.tick();
        if (recording_interval > 0 && (i % recording_interval) == 0) {
            save_tick(file, engine, 0);
        }

        if ((i + 1) % 50 == 0) {
            std::cout << "Tick " << (i + 1) << "/" << num_ticks
                      << " (" << engine.get_plant(0).node_count() << " nodes)" << std::endl;
        }
    }

    if (!chain_profile_path.empty()) {
        write_chain_profile(chain_profile_path, engine.get_plant(0));
    }

    std::cout << "Done. Final node count: " << engine.get_plant(0).node_count() << std::endl;
    return 0;
}
