// app_sugar_test.cpp — Static tree sugar economy tester.
// Builds hardcoded trees, freezes growth, runs N ticks of production/maintenance/transport.
// Usage: ./botany_sugar_test [--ticks N] [--csv] [--tree seedling|medium|large]

#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>
#include "engine/plant.h"
#include "engine/node/node.h"
#include "engine/node/leaf_node.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/sugar.h"
#include "engine/chemical/chemical.h"

using namespace botany;

// ── Frozen genome: all growth disabled, only sugar economy active ──

static Genome make_frozen_genome() {
    Genome g = default_genome();
    g.growth_rate = 0.0f;
    g.root_growth_rate = 0.0f;
    g.thickening_rate = 0.0f;
    g.internode_elongation_rate = 0.0f;
    g.root_internode_elongation_rate = 0.0f;
    g.leaf_growth_rate = 0.0f;
    g.shoot_plastochron = 999999;
    g.root_plastochron = 999999;
    g.auxin_production_rate = 0.0f;
    g.cytokinin_production_rate = 0.0f;
    g.ga_production_rate = 0.0f;
    g.leaf_phototropism_rate = 0.0f;
    // Keep starvation death off — we want to observe the economy, not kill nodes
    g.ethylene_starvation_rate = 0.0f;
    g.ethylene_shade_rate = 0.0f;
    g.ethylene_age_rate = 0.0f;
    g.ethylene_crowding_rate = 0.0f;
    return g;
}

// ── Per-tick stats ──

struct TickStats {
    float total_production = 0.0f;
    float total_maintenance = 0.0f;
    float net_balance = 0.0f;
    float total_sugar = 0.0f;
    float min_sugar = FLT_MAX;
    float max_sugar = 0.0f;
    int starving_count = 0;
    int node_count = 0;
};

static TickStats gather_stats(const Plant& plant, const WorldParams& world) {
    TickStats s{};
    plant.for_each_node([&](const Node& n) {
        s.node_count++;
        float sugar = n.chemical(ChemicalID::Sugar);
        s.total_sugar += sugar;
        s.min_sugar = std::min(s.min_sugar, sugar);
        s.max_sugar = std::max(s.max_sugar, sugar);
        if (sugar <= 0.0f) s.starving_count++;

        s.total_maintenance += n.maintenance_cost(world);

        if (auto* leaf = n.as_leaf()) {
            if (leaf->leaf_size > 1e-6f && leaf->senescence_ticks == 0) {
                float angle_eff = 1.0f;
                float flen = glm::length(leaf->facing);
                if (flen > 1e-4f) {
                    glm::vec3 normal = leaf->facing / flen;
                    angle_eff = std::max(0.0f, normal.y);
                }
                float area = leaf->leaf_size * leaf->leaf_size;
                s.total_production += leaf->light_exposure * angle_eff
                    * world.light_level * area * world.sugar_production_rate;
            }
        }
    });
    s.net_balance = s.total_production - s.total_maintenance;
    if (s.node_count == 0) s.min_sugar = 0.0f;
    return s;
}

// ── Tree description ──

struct TreeDesc {
    int node_count = 0;
    int leaf_count = 0;
    float total_leaf_area = 0.0f;
    float total_stem_volume = 0.0f;
    float total_root_volume = 0.0f;
    float initial_sugar = 0.0f;
};

static TreeDesc describe_tree(const Plant& plant) {
    TreeDesc d{};
    plant.for_each_node([&](const Node& n) {
        d.node_count++;
        d.initial_sugar += n.chemical(ChemicalID::Sugar);
        if (auto* leaf = n.as_leaf()) {
            d.leaf_count++;
            d.total_leaf_area += leaf->leaf_size * leaf->leaf_size;
        }
        float length = std::max(glm::length(n.offset), 0.01f);
        float volume = 3.14159f * n.radius * n.radius * length;
        if (n.type == NodeType::STEM) d.total_stem_volume += volume;
        if (n.type == NodeType::ROOT) d.total_root_volume += volume;
    });
    return d;
}

// ── Helpers ──

// Add a leaf to a parent node with given size and light exposure.
static Node* add_leaf(Plant& plant, Node* parent, glm::vec3 offset,
                       float leaf_size, float light_exposure) {
    Node* leaf = plant.create_node(NodeType::LEAF, offset, 0.0f);
    leaf->as_leaf()->leaf_size = leaf_size;
    leaf->as_leaf()->light_exposure = light_exposure;
    leaf->as_leaf()->facing = glm::vec3(0.0f, 1.0f, 0.0f);
    parent->add_child(leaf);
    return leaf;
}

// Add a stem internode to a parent.
static Node* add_stem(Plant& plant, Node* parent, glm::vec3 offset, float radius) {
    Node* stem = plant.create_node(NodeType::STEM, offset, radius);
    parent->add_child(stem);
    return stem;
}

// Add a root internode to a parent.
static Node* add_root(Plant& plant, Node* parent, glm::vec3 offset, float radius) {
    Node* root = plant.create_node(NodeType::ROOT, offset, radius);
    parent->add_child(root);
    return root;
}

// Re-assert light exposure on all leaves before each tick.
// The engine's ethylene/abscission system can trigger senescence which stops production.
// We also need to ensure light_exposure persists since there's no shadow map.
static void refresh_leaf_light(Plant& plant) {
    plant.for_each_node_mut([](Node& n) {
        if (auto* leaf = n.as_leaf()) {
            // Prevent senescence from shutting down production
            leaf->senescence_ticks = 0;
            leaf->deficit_ticks = 0;
        }
    });
}

// Set initial sugar distribution.
// Seed gets seed_sugar. Stems/roots get a small buffer. Leaves start empty
// (they need headroom to produce — if full, photosynthesize bails out).
static void distribute_initial_sugar(Plant& plant, const Genome& g) {
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::STEM || n.type == NodeType::ROOT) {
            float cap = sugar_cap(n, g);
            n.chemical(ChemicalID::Sugar) = cap * 0.2f;
        } else if (n.type == NodeType::LEAF) {
            n.chemical(ChemicalID::Sugar) = 0.0f; // empty — needs headroom to produce
        } else if (n.is_meristem()) {
            n.chemical(ChemicalID::Sugar) = 0.3f;
        }
    });
    // Seed keeps its seed_sugar
    plant.seed_mut()->chemical(ChemicalID::Sugar) = g.seed_sugar;
}

// ── Tree builders ──

static Plant build_seedling(const Genome& g) {
    Plant plant(g, glm::vec3(0.0f));
    Node* seed = plant.seed_mut();

    // Find auto-created meristems
    Node* shoot = nullptr;
    Node* root_tip = nullptr;
    for (Node* c : seed->children) {
        if (c->type == NodeType::SHOOT_APICAL) shoot = c;
        if (c->type == NodeType::ROOT_APICAL) root_tip = c;
    }

    // Insert 2 stem internodes between seed and shoot apical
    Node* stem1 = add_stem(plant, seed, glm::vec3(0.0f, 0.3f, 0.0f), g.initial_radius);
    Node* stem2 = add_stem(plant, stem1, glm::vec3(0.0f, 0.3f, 0.0f), g.initial_radius);

    // Move shoot apical to top of stem chain
    seed->replace_child(shoot, stem1);  // already added stem1 above, but need to remove shoot from seed
    // Oops — replace_child replaces in-place. Let me restructure.
    // Actually: seed already has stem1 as child (from add_stem).
    // And seed also has shoot from the constructor.
    // I need to remove shoot from seed's children and add it to stem2.

    // Remove shoot from seed (it was added by constructor)
    auto& sc = seed->children;
    sc.erase(std::remove(sc.begin(), sc.end(), shoot), sc.end());
    shoot->parent = nullptr;
    stem2->add_child(shoot);

    // Add leaves on each stem internode
    add_leaf(plant, stem1, glm::vec3(0.15f, 0.0f, 0.0f), 0.3f, 0.9f);
    add_leaf(plant, stem2, glm::vec3(-0.15f, 0.0f, 0.1f), 0.25f, 1.0f);

    distribute_initial_sugar(plant, g);
    return plant;
}

static Plant build_medium(const Genome& g) {
    Plant plant(g, glm::vec3(0.0f));
    Node* seed = plant.seed_mut();

    // Find auto-created meristems
    Node* shoot = nullptr;
    Node* root_tip = nullptr;
    for (Node* c : seed->children) {
        if (c->type == NodeType::SHOOT_APICAL) shoot = c;
        if (c->type == NodeType::ROOT_APICAL) root_tip = c;
    }

    // Remove shoot from seed — we'll place it at the top of the trunk
    auto& sc = seed->children;
    sc.erase(std::remove(sc.begin(), sc.end(), shoot), sc.end());
    shoot->parent = nullptr;

    // Remove root_tip from seed — we'll place it at the bottom of the root chain
    sc.erase(std::remove(sc.begin(), sc.end(), root_tip), sc.end());
    root_tip->parent = nullptr;

    // Build trunk: 5 stem internodes with increasing thickness
    float trunk_radii[] = {0.08f, 0.07f, 0.06f, 0.055f, 0.05f};
    float trunk_heights[] = {0.8f, 0.7f, 0.6f, 0.5f, 0.4f};
    float leaf_sizes[] = {1.2f, 1.0f, 0.9f, 0.8f, 0.6f};
    float leaf_light[] = {0.5f, 0.6f, 0.7f, 0.85f, 1.0f};  // lower leaves shaded

    Node* prev = seed;
    Node* trunk_nodes[5];
    for (int i = 0; i < 5; i++) {
        trunk_nodes[i] = add_stem(plant, prev, glm::vec3(0.0f, trunk_heights[i], 0.0f), trunk_radii[i]);
        // Leaf on each internode
        float angle = 2.4f * i;  // golden angle spiral
        glm::vec3 leaf_off(0.15f * std::cos(angle), 0.0f, 0.15f * std::sin(angle));
        add_leaf(plant, trunk_nodes[i], leaf_off, leaf_sizes[i], leaf_light[i]);
        prev = trunk_nodes[i];
    }

    // Shoot apical at trunk top
    prev->add_child(shoot);

    // Branch off internode 2: 3 stem internodes + 2 leaves
    Node* br1 = add_stem(plant, trunk_nodes[1], glm::vec3(0.4f, 0.2f, 0.0f), 0.04f);
    Node* br2 = add_stem(plant, br1, glm::vec3(0.3f, 0.15f, 0.0f), 0.035f);
    Node* br3 = add_stem(plant, br2, glm::vec3(0.2f, 0.1f, 0.0f), 0.03f);
    add_leaf(plant, br2, glm::vec3(0.1f, 0.05f, 0.1f), 0.9f, 0.7f);
    add_leaf(plant, br3, glm::vec3(0.1f, 0.05f, -0.1f), 0.7f, 0.85f);

    // Root system: 5 root internodes going down
    prev = seed;
    for (int i = 0; i < 5; i++) {
        float depth = -0.3f - 0.1f * i;
        prev = add_root(plant, prev, glm::vec3(0.05f * (i % 2 ? 1 : -1), depth, 0.0f), 0.025f);
    }
    prev->add_child(root_tip);

    distribute_initial_sugar(plant, g);
    return plant;
}

static Plant build_large(const Genome& g) {
    Plant plant(g, glm::vec3(0.0f));
    Node* seed = plant.seed_mut();

    // Find and detach auto-created meristems
    Node* shoot = nullptr;
    Node* root_tip = nullptr;
    for (Node* c : seed->children) {
        if (c->type == NodeType::SHOOT_APICAL) shoot = c;
        if (c->type == NodeType::ROOT_APICAL) root_tip = c;
    }
    auto& sc = seed->children;
    sc.erase(std::remove(sc.begin(), sc.end(), shoot), sc.end());
    shoot->parent = nullptr;
    sc.erase(std::remove(sc.begin(), sc.end(), root_tip), sc.end());
    root_tip->parent = nullptr;

    // Thick trunk: 5 internodes
    float trunk_radii[] = {0.15f, 0.13f, 0.11f, 0.09f, 0.07f};
    Node* prev = seed;
    Node* trunk[5];
    for (int i = 0; i < 5; i++) {
        trunk[i] = add_stem(plant, prev, glm::vec3(0.0f, 1.0f, 0.0f), trunk_radii[i]);
        prev = trunk[i];
    }
    prev->add_child(shoot);

    // Trunk leaves on each internode
    for (int i = 0; i < 5; i++) {
        float angle = 2.4f * i;
        glm::vec3 leaf_off(0.15f * std::cos(angle), 0.0f, 0.15f * std::sin(angle));
        add_leaf(plant, trunk[i], leaf_off, 1.2f, 0.5f + 0.1f * i);
    }

    // 3 branches, each off different trunk internodes
    struct BranchDef {
        int trunk_idx;
        glm::vec3 dir;
        int internodes;
        float base_radius;
    };
    BranchDef branches[] = {
        {1, glm::vec3(0.5f, 0.3f, 0.0f),  4, 0.06f},
        {2, glm::vec3(-0.3f, 0.3f, 0.4f), 3, 0.05f},
        {3, glm::vec3(0.0f, 0.3f, -0.5f), 4, 0.055f},
    };

    float branch_leaf_sizes[] = {1.3f, 1.1f, 1.0f, 0.8f};
    float branch_leaf_light[] = {0.6f, 0.75f, 0.85f, 1.0f};

    for (auto& bd : branches) {
        prev = trunk[bd.trunk_idx];
        for (int j = 0; j < bd.internodes; j++) {
            float r = bd.base_radius * (1.0f - 0.15f * j);
            glm::vec3 off = bd.dir * (0.4f + 0.1f * j);
            Node* br = add_stem(plant, prev, off, r);

            // Leaf on each branch internode
            float angle = 2.4f * j;
            glm::vec3 leaf_off(0.12f * std::cos(angle), 0.05f, 0.12f * std::sin(angle));
            float ls = (j < 4) ? branch_leaf_sizes[j] : 0.8f;
            float ll = (j < 4) ? branch_leaf_light[j] : 0.9f;
            add_leaf(plant, br, leaf_off, ls, ll);

            prev = br;
        }
    }

    // Root system: 10 internodes, branching
    prev = seed;
    Node* root_chain[10];
    for (int i = 0; i < 10; i++) {
        float r = 0.05f - 0.003f * i;
        glm::vec3 off(0.03f * ((i % 3) - 1), -0.3f, 0.03f * ((i % 2) ? 1 : -1));
        root_chain[i] = add_root(plant, prev, off, std::max(r, 0.01f));
        prev = root_chain[i];
    }
    prev->add_child(root_tip);

    distribute_initial_sugar(plant, g);
    return plant;
}

// ── Run simulation ──

static void print_header() {
    std::cout << std::fixed << std::setprecision(4);
}

static void print_tree_description(const char* name, const TreeDesc& d) {
    std::cout << "\n=== " << name << " ===\n";
    std::cout << "  Nodes: " << d.node_count
              << "  Leaves: " << d.leaf_count
              << "  Leaf area: " << d.total_leaf_area << " dm²\n";
    std::cout << "  Stem volume: " << d.total_stem_volume << " dm³"
              << "  Root volume: " << d.total_root_volume << " dm³\n";
    std::cout << "  Initial sugar: " << d.initial_sugar << " g\n";
}

static void run_tree(const char* name, Plant& plant, const WorldParams& world,
                     int num_ticks, bool csv_mode) {
    const Genome& g = plant.genome();
    TreeDesc desc = describe_tree(plant);

    if (!csv_mode) {
        print_tree_description(name, desc);
    }

    // Accumulators for summary
    float sum_production = 0.0f;
    float sum_maintenance = 0.0f;
    int peak_starving = 0;

    for (int tick = 0; tick < num_ticks; tick++) {
        // Re-assert leaf state (prevent senescence in static test)
        refresh_leaf_light(plant);

        // Simulate growth consumption: in a real sim, meristems consume sugar
        // for growth, keeping leaves below cap. With growth frozen, leaves fill
        // to cap and photosynthesize bails out. Drain leaf sugar before each tick
        // so they have headroom to produce.
        plant.for_each_node_mut([&](Node& n) {
            if (n.type == NodeType::LEAF) {
                // In a real sim, growth drains sugar from the tree so leaves always
                // have headroom to produce. With growth frozen, we must drain leaves
                // manually to prevent photosynthesize from bailing at cap.


                n.chemical(ChemicalID::Sugar) = 0.0f;
            }
        });

        // Gather pre-tick stats (maintenance estimate + tree state)
        TickStats pre = gather_stats(plant, world);
        float sugar_produced_before = plant.total_sugar_produced();

        plant.tick(world);

        // Gather post-tick stats, use actual production from engine
        TickStats post = gather_stats(plant, world);
        float actual_production = plant.total_sugar_produced() - sugar_produced_before;
        pre.total_production = actual_production;
        pre.net_balance = actual_production - pre.total_maintenance;

        sum_production += pre.total_production;
        sum_maintenance += pre.total_maintenance;
        peak_starving = std::max(peak_starving, post.starving_count);

        if (csv_mode) {
            std::cout << name << ","
                      << tick << ","
                      << pre.total_production << ","
                      << pre.total_maintenance << ","
                      << pre.net_balance << ","
                      << post.total_sugar << ","
                      << post.min_sugar << ","
                      << post.max_sugar << ","
                      << post.starving_count << "\n";
        }
    }

    if (!csv_mode) {
        float avg_prod = sum_production / num_ticks;
        float avg_maint = sum_maintenance / num_ticks;
        TickStats final_stats = gather_stats(plant, world);

        std::cout << "  --- After " << num_ticks << " ticks ---\n";
        std::cout << "  Avg production:  " << avg_prod << " g/hr\n";
        std::cout << "  Avg maintenance: " << avg_maint << " g/hr\n";
        std::cout << "  Net balance:     " << (avg_prod - avg_maint) << " g/hr";
        if (avg_prod - avg_maint >= 0) std::cout << " (surplus)";
        else                            std::cout << " (DEFICIT)";
        std::cout << "\n";
        std::cout << "  Prod/maint ratio: " << (avg_maint > 0 ? avg_prod / avg_maint : 0) << "x\n";
        std::cout << "  Final total sugar: " << final_stats.total_sugar << " g\n";
        std::cout << "  Sugar range: [" << final_stats.min_sugar << ", " << final_stats.max_sugar << "]\n";
        std::cout << "  Peak starving nodes: " << peak_starving << "\n";
    }
}

// ── Main ──

int main(int argc, char* argv[]) {
    int num_ticks = 500;
    bool csv_mode = false;
    const char* tree_filter = nullptr;

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--ticks") == 0 && i + 1 < argc) {
            num_ticks = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "--csv") == 0) {
            csv_mode = true;
        } else if (std::strcmp(argv[i], "--tree") == 0 && i + 1 < argc) {
            tree_filter = argv[++i];
        }
    }

    Genome g = make_frozen_genome();
    WorldParams world = default_world_params();

    // Disable starvation death — we want to observe, not prune
    world.starvation_ticks_max = 999999;

    if (csv_mode) {
        std::cout << "tree,tick,production,maintenance,net,total_sugar,min_sugar,max_sugar,starving\n";
    } else {
        std::cout << "Sugar Economy Test — " << num_ticks << " ticks\n";
        std::cout << "Production rate: " << world.sugar_production_rate << " g/(dm²·hr)\n";
        std::cout << "Maintenance: leaf=" << world.sugar_maintenance_leaf
                  << " stem=" << world.sugar_maintenance_stem
                  << " root=" << world.sugar_maintenance_root
                  << " meristem=" << world.sugar_maintenance_meristem << "\n";
    }

    struct TreeDef {
        const char* name;
        Plant (*builder)(const Genome&);
    };
    TreeDef trees[] = {
        {"seedling", build_seedling},
        {"medium",   build_medium},
        {"large",    build_large},
    };

    for (auto& def : trees) {
        if (tree_filter && std::strcmp(tree_filter, def.name) != 0) continue;
        Plant plant = def.builder(g);
        run_tree(def.name, plant, world, num_ticks, csv_mode);
    }

    return 0;
}
