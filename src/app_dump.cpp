// src/app_dump.cpp
// Reads a binary recording and dumps it as human-readable text for debugging.
//
// Usage:
//   ./build/botany_dump recording.bin                  # dump all ticks
//   ./build/botany_dump recording.bin 10               # dump only tick 10
//   ./build/botany_dump recording.bin 0 5              # dump ticks 0 through 5
//   ./build/botany_dump recording.bin -tree 10         # dump tick 10 as indented tree
//   ./build/botany_dump recording.bin -stats           # summary stats per tick

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <glm/geometric.hpp>
#include "serialization/serializer.h"

using namespace botany;

static const char* node_type_str(NodeType t) {
    switch (t) {
        case NodeType::STEM:           return "STEM";
        case NodeType::ROOT:           return "ROOT";
        case NodeType::LEAF:           return "LEAF";
        case NodeType::APICAL:         return "APICAL";
        case NodeType::ROOT_APICAL:    return "ROOT_APICAL";
    }
    return "UNKNOWN";
}

static void print_genome(const Genome& g) {
    std::cout << "=== GENOME ===" << std::endl;
    std::cout << "  apical_auxin_baseline:    " << g.apical_auxin_baseline << std::endl;
    std::cout << "  auxin_diffusion_rate:     " << g.auxin_diffusion_rate << std::endl;
    std::cout << "  auxin_decay_rate:         " << g.auxin_decay_rate << std::endl;
    std::cout << "  auxin_threshold:          " << g.auxin_threshold << std::endl;
    std::cout << "  cytokinin_production_rate:" << g.cytokinin_production_rate << std::endl;
    std::cout << "  cytokinin_diffusion_rate: " << g.cytokinin_diffusion_rate << std::endl;
    std::cout << "  cytokinin_decay_rate:     " << g.cytokinin_decay_rate << std::endl;
    std::cout << "  cytokinin_threshold:      " << g.cytokinin_threshold << std::endl;
    std::cout << "  growth_rate:              " << g.growth_rate << std::endl;
    std::cout << "  max_internode_length:     " << g.max_internode_length << std::endl;
    std::cout << "  branch_angle:             " << g.branch_angle << " rad (" << (g.branch_angle * 180.0f / 3.14159f) << " deg)" << std::endl;
    std::cout << "  cambium_responsiveness:        " << g.cambium_responsiveness << std::endl;
    std::cout << "  vascular_radius_threshold:     " << g.vascular_radius_threshold << std::endl;
    std::cout << "  root_growth_rate:         " << g.root_growth_rate << std::endl;
    std::cout << "  max_internode_length:     " << g.max_internode_length << std::endl;
    std::cout << "  root_branch_angle:        " << g.root_branch_angle << " rad (" << (g.root_branch_angle * 180.0f / 3.14159f) << " deg)" << std::endl;
    std::cout << "  max_leaf_size:            " << g.max_leaf_size << std::endl;
    std::cout << "  initial_radius:           " << g.initial_radius << std::endl;
    std::cout << std::endl;
}

static void print_tick_full(const TickSnapshot& snap) {
    std::cout << "--- TICK " << snap.tick_number << " (" << snap.nodes.size() << " nodes) ---" << std::endl;

    // Build parent->children map and id lookup
    std::unordered_map<uint32_t, std::vector<uint32_t>> children_map;
    std::unordered_map<uint32_t, size_t> id_to_idx;
    for (size_t i = 0; i < snap.nodes.size(); i++) {
        const auto& n = snap.nodes[i];
        id_to_idx[n.id] = i;
        if (n.parent_id != UINT32_MAX) {
            children_map[n.parent_id].push_back(n.id);
        }
    }

    for (const auto& n : snap.nodes) {
        float dist_to_parent = 0.0f;
        if (n.parent_id != UINT32_MAX) {
            auto it = id_to_idx.find(n.parent_id);
            if (it != id_to_idx.end()) {
                dist_to_parent = glm::length(n.position - snap.nodes[it->second].position);
            }
        }

        int num_children = 0;
        auto cit = children_map.find(n.id);
        if (cit != children_map.end()) {
            num_children = static_cast<int>(cit->second.size());
        }

        std::cout << "  node[" << n.id << "]"
                  << " type=" << node_type_str(n.type)
                  << " parent=" << (n.parent_id == UINT32_MAX ? -1 : static_cast<int>(n.parent_id))
                  << " children=" << num_children
                  << " pos=(" << n.position.x << ", " << n.position.y << ", " << n.position.z << ")"
                  << " r=" << n.radius
                  << " dist=" << dist_to_parent
                  << " auxin=" << n.auxin
                  << " cyto=" << n.cytokinin
                  << " sugar=" << n.sugar;
        if (n.type == NodeType::LEAF) {
            std::cout << " LEAF(" << n.leaf_size << ")";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

static void print_tick_tree(const TickSnapshot& snap) {
    std::cout << "--- TICK " << snap.tick_number << " (" << snap.nodes.size() << " nodes) TREE ---" << std::endl;

    // Build lookup structures
    std::unordered_map<uint32_t, size_t> id_to_idx;
    std::unordered_map<uint32_t, std::vector<uint32_t>> children_map;
    uint32_t root_id = UINT32_MAX;

    for (size_t i = 0; i < snap.nodes.size(); i++) {
        const auto& n = snap.nodes[i];
        id_to_idx[n.id] = i;
        if (n.parent_id == UINT32_MAX) {
            root_id = n.id;
        } else {
            children_map[n.parent_id].push_back(n.id);
        }
    }

    // Recursive tree print
    struct TreePrinter {
        const TickSnapshot& snap;
        const std::unordered_map<uint32_t, size_t>& id_to_idx;
        const std::unordered_map<uint32_t, std::vector<uint32_t>>& children_map;

        void print(uint32_t id, const std::string& prefix, bool is_last) {
            auto it = id_to_idx.find(id);
            if (it == id_to_idx.end()) return;
            const auto& n = snap.nodes[it->second];

            float dist = 0.0f;
            if (n.parent_id != UINT32_MAX) {
                auto pit = id_to_idx.find(n.parent_id);
                if (pit != id_to_idx.end()) {
                    dist = glm::length(n.position - snap.nodes[pit->second].position);
                }
            }

            std::cout << prefix << (is_last ? "└── " : "├── ");
            std::cout << "[" << n.id << "] " << node_type_str(n.type)
                      << " r=" << n.radius
                      << " dist=" << dist
                      << " y=" << n.position.y
                      << " aux=" << n.auxin
                      << " cyt=" << n.cytokinin
                      << " sugar=" << n.sugar;
            if (n.type == NodeType::LEAF) std::cout << " LEAF";

            // Count children
            auto cit = children_map.find(n.id);
            int nc = cit != children_map.end() ? static_cast<int>(cit->second.size()) : 0;
            if (nc > 0) std::cout << " (" << nc << " children)";
            std::cout << std::endl;

            if (cit != children_map.end()) {
                const auto& kids = cit->second;
                for (size_t i = 0; i < kids.size(); i++) {
                    print(kids[i], prefix + (is_last ? "    " : "│   "), i == kids.size() - 1);
                }
            }
        }
    };

    if (root_id != UINT32_MAX) {
        TreePrinter tp{snap, id_to_idx, children_map};
        tp.print(root_id, "", true);
    }
    std::cout << std::endl;
}

static void print_tick_stats(const TickSnapshot& snap) {
    int stem_count = 0, root_count = 0, leaf_count = 0, leaf_node_count = 0;
    float min_radius = 1e9f, max_radius = 0.0f;
    float min_auxin = 1e9f, max_auxin = 0.0f, sum_auxin = 0.0f;
    float min_cyto = 1e9f, max_cyto = 0.0f, sum_cyto = 0.0f;
    float max_y = -1e9f, min_y = 1e9f;
    int max_children = 0;

    std::unordered_map<uint32_t, int> children_count;
    for (const auto& n : snap.nodes) {
        if (n.parent_id != UINT32_MAX) {
            children_count[n.parent_id]++;
        }
    }

    for (const auto& n : snap.nodes) {
        if (n.type == NodeType::STEM) stem_count++;
        else if (n.type == NodeType::ROOT) root_count++;
        else if (n.type == NodeType::LEAF) leaf_node_count++;
        if (n.type == NodeType::LEAF) leaf_count++;

        min_radius = std::min(min_radius, n.radius);
        max_radius = std::max(max_radius, n.radius);
        min_auxin = std::min(min_auxin, n.auxin);
        max_auxin = std::max(max_auxin, n.auxin);
        sum_auxin += n.auxin;
        min_cyto = std::min(min_cyto, n.cytokinin);
        max_cyto = std::max(max_cyto, n.cytokinin);
        sum_cyto += n.cytokinin;
        max_y = std::max(max_y, n.position.y);
        min_y = std::min(min_y, n.position.y);

        auto it = children_count.find(n.id);
        if (it != children_count.end()) {
            max_children = std::max(max_children, it->second);
        }
    }

    float n = static_cast<float>(snap.nodes.size());
    std::cout << "tick=" << snap.tick_number
              << " nodes=" << snap.nodes.size()
              << " stem=" << stem_count
              << " root=" << root_count
              << " leaves=" << leaf_count
              << " | radius[" << min_radius << ".." << max_radius << "]"
              << " | auxin[" << min_auxin << ".." << max_auxin << " avg=" << (sum_auxin / n) << "]"
              << " | cyto[" << min_cyto << ".." << max_cyto << " avg=" << (sum_cyto / n) << "]"
              << " | y[" << min_y << ".." << max_y << "]"
              << " | max_children=" << max_children
              << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage:" << std::endl;
        std::cerr << "  botany_dump <file.bin>                   # dump all ticks" << std::endl;
        std::cerr << "  botany_dump <file.bin> <tick>             # dump one tick" << std::endl;
        std::cerr << "  botany_dump <file.bin> <from> <to>        # dump tick range" << std::endl;
        std::cerr << "  botany_dump <file.bin> -tree <tick>        # tree view of one tick" << std::endl;
        std::cerr << "  botany_dump <file.bin> -stats              # one-line stats per tick" << std::endl;
        return 1;
    }

    std::string path = argv[1];
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        std::cerr << "Failed to open " << path << std::endl;
        return 1;
    }

    RecordingHeader header = load_recording_header(file);
    print_genome(header.genome);
    std::cout << "Recording: " << header.num_ticks << " ticks" << std::endl << std::endl;

    // Parse mode
    enum Mode { FULL, TREE, STATS } mode = FULL;
    int tick_from = -1, tick_to = -1;
    int tick_single = -1;

    if (argc >= 3) {
        std::string arg2 = argv[2];
        if (arg2 == "-stats") {
            mode = STATS;
        } else if (arg2 == "-tree") {
            mode = TREE;
            if (argc >= 4) tick_single = std::atoi(argv[3]);
            else tick_single = static_cast<int>(header.num_ticks) - 1;
        } else {
            tick_from = std::atoi(argv[2]);
            tick_to = tick_from;
            if (argc >= 4) {
                std::string arg3 = argv[3];
                if (arg3 != "-tree" && arg3 != "-stats") {
                    tick_to = std::atoi(argv[3]);
                }
            }
        }
    }

    for (uint32_t i = 0; i < header.num_ticks; i++) {
        TickSnapshot snap = load_tick(file);

        if (mode == STATS) {
            print_tick_stats(snap);
            continue;
        }

        if (mode == TREE) {
            if (tick_single == -1 || static_cast<int>(snap.tick_number) == tick_single) {
                print_tick_tree(snap);
            }
            continue;
        }

        // FULL mode
        if (tick_from == -1) {
            // dump all
            print_tick_full(snap);
        } else if (static_cast<int>(snap.tick_number) >= tick_from &&
                   static_cast<int>(snap.tick_number) <= tick_to) {
            print_tick_full(snap);
        }
    }

    return 0;
}
