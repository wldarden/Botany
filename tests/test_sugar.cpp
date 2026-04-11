#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <glm/geometric.hpp>
#include "engine/plant.h"
#include "engine/node/leaf_node.h"
#include "engine/node/meristem_node.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/sugar.h"

using namespace botany;
using Catch::Matchers::WithinAbs;

// === Production tests ===

TEST_CASE("LEAF nodes produce sugar proportional to light and leaf_size", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 2.0f;

    produce_sugar(plant, wp);

    float expected = wp.light_level * leaf->as_leaf()->leaf_size * g.sugar_production_rate;
    REQUIRE_THAT(leaf->sugar, WithinAbs(expected, 1e-6));
}

TEST_CASE("Zero light produces zero sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 0.0f;

    produce_sugar(plant, wp);

    REQUIRE(leaf->sugar == 0.0f);
}

TEST_CASE("Non-LEAF nodes do not produce sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Zero out any initial sugar (e.g. seed_sugar) so we only measure production
    plant.for_each_node_mut([](Node& n) { n.sugar = 0.0f; });

    WorldParams wp = default_world_params();
    produce_sugar(plant, wp);

    plant.for_each_node([](const Node& n) {
        if (n.type != NodeType::LEAF) {
            REQUIRE(n.sugar == 0.0f);
        }
    });
}

// === Consumption tests ===

TEST_CASE("consume_sugar deducts volume-based maintenance cost", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Create a stem with known dimensions for predictable cost
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.1f);
    plant.seed_mut()->add_child(stem);

    // Set sugar within cap so safety clamp doesn't interfere
    float cap = sugar_cap(*stem, g);
    float initial_sugar = cap * 0.9f;
    stem->sugar = initial_sugar;

    // Zero out other nodes' sugar so we only check our stem
    plant.seed_mut()->sugar = sugar_cap(*plant.seed_mut(), g);

    consume_sugar(plant);

    float length = glm::length(stem->offset);  // 1.0
    float volume = 3.14159f * 0.1f * 0.1f * length;
    float expected_cost = g.sugar_maintenance_stem * volume;
    REQUIRE_THAT(stem->sugar, WithinAbs(initial_sugar - expected_cost, 1e-5));
}

TEST_CASE("LEAF maintenance cost uses leaf area", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    plant.seed_mut()->add_child(leaf);

    // Set sugar within cap so safety clamp doesn't interfere
    float cap = sugar_cap(*leaf, g);
    float initial_sugar = cap * 0.9f;
    leaf->sugar = initial_sugar;

    consume_sugar(plant);

    float expected_cost = g.sugar_maintenance_leaf * 0.5f * 0.5f;
    REQUIRE_THAT(leaf->sugar, WithinAbs(initial_sugar - expected_cost, 1e-6));
}

TEST_CASE("Sugar cannot go below zero after consumption", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    seed->sugar = 0.0001f;

    consume_sugar(plant);

    REQUIRE(seed->sugar >= 0.0f);
}

TEST_CASE("Active meristem tips have additional maintenance cost", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* shoot_tip = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.is_meristem() && n.as_meristem()->is_tip() && n.as_meristem()->active &&
            n.type == NodeType::SHOOT_APICAL) {
            shoot_tip = &n;
        }
    });
    REQUIRE(shoot_tip != nullptr);

    // Set sugar within cap so safety clamp doesn't interfere
    float cap = sugar_cap(*shoot_tip, g);
    float initial_sugar = cap * 0.9f;
    shoot_tip->sugar = initial_sugar;

    consume_sugar(plant);

    // Meristem nodes have no volume-based maintenance, only the per-tip cost
    float expected_cost = g.sugar_maintenance_meristem;
    REQUIRE_THAT(shoot_tip->sugar, WithinAbs(initial_sugar - expected_cost, 1e-5));
}

// === Diffusion tests ===

TEST_CASE("Sugar diffuses from high to low concentration", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    seed->sugar = 100.0f;

    float seed_before = seed->sugar;
    // Transport chemicals locally on each child (pulls sugar from seed)
    for (Node* child : seed->children) {
        child->transport_chemicals(g);
    }

    REQUIRE(seed->sugar < seed_before);
    for (const Node* child : seed->children) {
        REQUIRE(child->sugar > 0.0f);
    }
}

TEST_CASE("Sugar diffusion preserves total sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    plant.seed_mut()->sugar = 100.0f;

    float total_before = 0.0f;
    plant.for_each_node([&](const Node& n) { total_before += n.sugar; });

    // Run several rounds of local transport (one pass per round)
    for (int i = 0; i < 5; i++) {
        plant.for_each_node_mut([&](Node& n) {
            n.transport_chemicals(g);
        });
    }

    float total_after = 0.0f;
    plant.for_each_node([&](const Node& n) { total_after += n.sugar; });

    REQUIRE_THAT(total_after, WithinAbs(total_before, 1e-4));
}

TEST_CASE("More transport ticks produce smoother distribution", "[sugar]") {
    Genome g = default_genome();

    // Create plants with large-cap children so cap clamping doesn't interfere
    Plant plant1(g, glm::vec3(0.0f));
    for (int i = 0; i < 3; i++) {
        Node* c = plant1.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
        plant1.seed_mut()->add_child(c);
    }
    plant1.seed_mut()->sugar = 10.0f;
    // 1 round of local transport
    plant1.for_each_node_mut([&](Node& n) { n.transport_chemicals(g); });

    Plant plant2(g, glm::vec3(0.0f));
    for (int i = 0; i < 3; i++) {
        Node* c = plant2.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
        plant2.seed_mut()->add_child(c);
    }
    plant2.seed_mut()->sugar = 10.0f;
    // 10 rounds of local transport
    for (int i = 0; i < 10; i++) {
        plant2.for_each_node_mut([&](Node& n) { n.transport_chemicals(g); });
    }

    float child_sugar_1tick = 0.0f;
    float child_sugar_10tick = 0.0f;
    for (const Node* c : plant1.seed()->children) { child_sugar_1tick += c->sugar; }
    for (const Node* c : plant2.seed()->children) { child_sugar_10tick += c->sugar; }

    REQUIRE(child_sugar_10tick > child_sugar_1tick);
}

TEST_CASE("Sugar diffusion conserves total across many ticks", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Build a bigger tree: seed -> 5 children, each with 2 grandchildren
    for (int i = 0; i < 5; i++) {
        Node* child = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
        plant.seed_mut()->add_child(child);
        for (int j = 0; j < 2; j++) {
            Node* gc = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.03f);
            child->add_child(gc);
        }
    }

    // Give sugar to just a few nodes
    plant.seed_mut()->sugar = 50.0f;
    plant.seed_mut()->children[0]->sugar = 20.0f;
    plant.seed_mut()->children[2]->sugar = 10.0f;

    float total_before = 0.0f;
    plant.for_each_node([&](const Node& n) { total_before += n.sugar; });

    for (int ticks : {1, 5, 15, 50}) {
        // Reset to initial distribution
        plant.for_each_node_mut([](Node& n) { n.sugar = 0.0f; });
        plant.seed_mut()->sugar = 50.0f;
        plant.seed_mut()->children[0]->sugar = 20.0f;
        plant.seed_mut()->children[2]->sugar = 10.0f;

        for (int t = 0; t < ticks; t++) {
            plant.for_each_node_mut([&](Node& n) {
                n.transport_chemicals(g);
            });
        }

        float total_after = 0.0f;
        plant.for_each_node([&](const Node& n) { total_after += n.sugar; });

        REQUIRE_THAT(total_after, WithinAbs(total_before, 1e-3));
    }
}

// === Starvation tests ===

TEST_CASE("Starvation ticks increment when sugar is zero", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    seed->sugar = 0.0f;

    consume_sugar(plant);

    REQUIRE(seed->starvation_ticks > 0);
}

TEST_CASE("Starvation ticks reset when sugar is available", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    seed->starvation_ticks = 10;
    seed->sugar = 100.0f;

    consume_sugar(plant);

    REQUIRE(seed->starvation_ticks == 0);
}

TEST_CASE("Starved non-seed nodes are pruned after max starvation ticks", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Create a child node and starve it
    Node* child = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(child);

    uint32_t count_before = plant.node_count();
    child->starvation_ticks = 100; // way past max

    WorldParams wp = default_world_params();
    wp.starvation_ticks_max = 50;

    prune_starved_nodes(plant, wp);

    REQUIRE(plant.node_count() < count_before);
}

TEST_CASE("Seed node is never pruned even when starved", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    plant.seed_mut()->starvation_ticks = 1000; // extremely starved

    WorldParams wp = default_world_params();
    wp.starvation_ticks_max = 50;

    uint32_t count_before = plant.node_count();
    prune_starved_nodes(plant, wp);

    // Seed should still exist
    REQUIRE(plant.seed() != nullptr);
}

TEST_CASE("Subtree removal cleans up all descendants", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Create a chain: seed -> A -> B -> C
    Node* a = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    Node* b = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 2.0f, 0.0f), 0.05f);
    Node* c = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 3.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(a);
    a->add_child(b);
    b->add_child(c);

    uint32_t a_id = a->id;

    // Starve A — should remove A, B, and C
    a->starvation_ticks = 100;

    WorldParams wp = default_world_params();
    wp.starvation_ticks_max = 50;

    prune_starved_nodes(plant, wp);

    // A, B, C should all be gone
    bool found_a = false;
    plant.for_each_node([&](const Node& n) {
        if (n.id == a_id) found_a = true;
    });
    REQUIRE_FALSE(found_a);
}

// === Integration test ===

TEST_CASE("transport_sugar runs full produce-diffuse-consume cycle", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 1.0f;
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 5.0f;

    transport_sugar(plant, wp);

    REQUIRE(leaf->sugar > 0.0f);
    REQUIRE(plant.seed()->sugar > 0.0f);
}

// === Sugar cap tests ===

TEST_CASE("sugar_cap scales with stem volume", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Create a stem with known dimensions: radius 0.1, internode length 1.0
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.1f);
    plant.seed_mut()->add_child(stem);

    float volume = 3.14159f * 0.1f * 0.1f * 1.0f;  // π * r² * length
    float expected = volume * g.sugar_storage_density_wood;
    REQUIRE_THAT(sugar_cap(*stem, g), WithinAbs(expected, 1e-4));
}

TEST_CASE("sugar_cap scales with leaf area", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;  // large enough that area cap exceeds minimum
    plant.seed_mut()->add_child(leaf);

    float area = 0.5f * 0.5f;
    float expected = area * g.sugar_storage_density_leaf;
    REQUIRE(expected > g.sugar_cap_minimum);  // precondition: area cap wins over floor
    REQUIRE_THAT(sugar_cap(*leaf, g), WithinAbs(expected, 1e-6));
}

TEST_CASE("sugar_cap returns minimum for tiny nodes", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Tiny node: radius 0.001, offset 0.001 — volume-based cap is below minimum
    Node* tiny = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.001f, 0.0f), 0.001f);
    plant.seed_mut()->add_child(tiny);

    REQUIRE_THAT(sugar_cap(*tiny, g), WithinAbs(g.sugar_cap_minimum, 1e-6));
}

TEST_CASE("sugar_cap for seed node covers seed_sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Seed node (no parent) should have cap >= seed_sugar
    REQUIRE(sugar_cap(*plant.seed(), g) >= g.seed_sugar);
}

TEST_CASE("sugar_cap for root scales with volume", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* root = plant.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.8f, 0.0f), 0.025f);
    plant.seed_mut()->add_child(root);

    float length = 0.8f;  // glm::length of offset
    float volume = 3.14159f * 0.025f * 0.025f * length;
    float expected = volume * g.sugar_storage_density_wood;
    // Volume cap (0.0785) > minimum (0.01), so volume wins
    REQUIRE_THAT(sugar_cap(*root, g), WithinAbs(expected, 1e-5));
}

TEST_CASE("Production skipped when leaf sugar is at cap", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.3f;
    plant.seed_mut()->add_child(leaf);

    // Fill leaf to its cap
    float cap = sugar_cap(*leaf, g);
    leaf->sugar = cap;

    WorldParams wp = default_world_params();
    produce_sugar(plant, wp);

    // Sugar should not have increased
    REQUIRE_THAT(leaf->sugar, WithinAbs(cap, 1e-6));
}

TEST_CASE("Production works normally when leaf is below cap", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    leaf->sugar = 0.0f;
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 2.0f;

    produce_sugar(plant, wp);

    REQUIRE(leaf->sugar > 0.0f);
}

TEST_CASE("Diffusion does not push sugar past receiver cap", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Create a child node with small cap (tiny radius, short offset)
    Node* child = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.1f, 0.0f), 0.02f);
    plant.seed_mut()->add_child(child);

    float child_cap = sugar_cap(*child, g);

    // Give seed lots of sugar, child already at cap
    plant.seed_mut()->sugar = 100.0f;
    child->sugar = child_cap;

    // Run several rounds of local transport
    for (int i = 0; i < 10; i++) {
        child->transport_chemicals(g);
    }

    // Child should not exceed its cap
    REQUIRE(child->sugar <= child_cap + 1e-6f);
}

TEST_CASE("Safety clamp caps sugar after consumption", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.1f);
    plant.seed_mut()->add_child(stem);

    float cap = sugar_cap(*stem, g);
    // Set sugar way above cap (simulating an edge case)
    stem->sugar = cap * 10.0f;

    consume_sugar(plant);

    REQUIRE(stem->sugar <= cap + 1e-6f);
}

TEST_CASE("Senescing leaf produces no sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    leaf->as_leaf()->senescence_ticks = 1; // senescing
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 2.0f;

    produce_sugar(plant, wp);

    REQUIRE(leaf->sugar == 0.0f);
}

TEST_CASE("Diffusion still conserves sugar when caps are not hit", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Build a tree with nodes that have large caps (big radii, long offsets)
    for (int i = 0; i < 3; i++) {
        Node* child = plant.create_node(NodeType::STEM,
            glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
        plant.seed_mut()->add_child(child);
    }

    // Give modest sugar amounts (well below caps)
    plant.seed_mut()->sugar = 5.0f;
    plant.seed_mut()->children[0]->sugar = 2.0f;

    float total_before = 0.0f;
    plant.for_each_node([&](const Node& n) { total_before += n.sugar; });

    // Run several rounds of local transport
    for (int i = 0; i < 10; i++) {
        plant.for_each_node_mut([&](Node& n) {
            n.transport_chemicals(g);
        });
    }

    float total_after = 0.0f;
    plant.for_each_node([&](const Node& n) { total_after += n.sugar; });

    REQUIRE_THAT(total_after, WithinAbs(total_before, 1e-4));
}
