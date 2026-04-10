#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/plant.h"
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
    leaf->leaf_size = 0.5f;
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 2.0f;

    produce_sugar(plant, wp);

    float expected = wp.light_level * leaf->leaf_size * g.sugar_production_rate;
    REQUIRE_THAT(leaf->sugar, WithinAbs(expected, 1e-6));
}

TEST_CASE("Zero light produces zero sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.5f;
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

TEST_CASE("consume_sugar deducts maintenance cost by node type", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    seed->sugar = 10.0f;

    consume_sugar(plant);

    float expected_cost = g.sugar_maintenance_stem * seed->radius;
    REQUIRE_THAT(seed->sugar, WithinAbs(10.0f - expected_cost, 1e-6));
}

TEST_CASE("LEAF maintenance cost uses leaf_size", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->leaf_size = 0.5f;
    leaf->sugar = 10.0f;
    plant.seed_mut()->add_child(leaf);

    consume_sugar(plant);

    float expected_cost = g.sugar_maintenance_leaf * leaf->leaf_size;
    REQUIRE_THAT(leaf->sugar, WithinAbs(10.0f - expected_cost, 1e-6));
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
        if (n.meristem && n.meristem->is_tip() && n.meristem->active &&
            n.type == NodeType::STEM) {
            shoot_tip = &n;
        }
    });
    REQUIRE(shoot_tip != nullptr);
    shoot_tip->sugar = 10.0f;

    consume_sugar(plant);

    float expected_cost = g.sugar_maintenance_stem * shoot_tip->radius
                        + g.sugar_maintenance_meristem;
    REQUIRE_THAT(shoot_tip->sugar, WithinAbs(10.0f - expected_cost, 1e-6));
}

// === Diffusion tests ===

TEST_CASE("Sugar diffuses from high to low concentration", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    seed->sugar = 100.0f;

    WorldParams wp = default_world_params();
    wp.sugar_diffusion_iterations = 1;

    float seed_before = seed->sugar;
    diffuse_sugar(plant, wp);

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

    WorldParams wp = default_world_params();
    wp.sugar_diffusion_iterations = 5;
    diffuse_sugar(plant, wp);

    float total_after = 0.0f;
    plant.for_each_node([&](const Node& n) { total_after += n.sugar; });

    REQUIRE_THAT(total_after, WithinAbs(total_before, 1e-4));
}

TEST_CASE("Multiple diffusion iterations produce smoother distribution", "[sugar]") {
    Genome g = default_genome();

    Plant plant1(g, glm::vec3(0.0f));
    plant1.seed_mut()->sugar = 100.0f;
    WorldParams wp1 = default_world_params();
    wp1.sugar_diffusion_iterations = 1;
    diffuse_sugar(plant1, wp1);

    Plant plant2(g, glm::vec3(0.0f));
    plant2.seed_mut()->sugar = 100.0f;
    WorldParams wp2 = default_world_params();
    wp2.sugar_diffusion_iterations = 10;
    diffuse_sugar(plant2, wp2);

    float child_sugar_1iter = 0.0f;
    float child_sugar_10iter = 0.0f;
    for (const Node* c : plant1.seed()->children) { child_sugar_1iter += c->sugar; }
    for (const Node* c : plant2.seed()->children) { child_sugar_10iter += c->sugar; }

    REQUIRE(child_sugar_10iter > child_sugar_1iter);
}

TEST_CASE("Sugar diffusion conserves total across many iterations", "[sugar]") {
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

    WorldParams wp = default_world_params();
    for (int iters : {1, 5, 15, 50}) {
        // Reset to initial distribution
        plant.for_each_node_mut([](Node& n) { n.sugar = 0.0f; });
        plant.seed_mut()->sugar = 50.0f;
        plant.seed_mut()->children[0]->sugar = 20.0f;
        plant.seed_mut()->children[2]->sugar = 10.0f;

        wp.sugar_diffusion_iterations = iters;
        diffuse_sugar(plant, wp);

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
    leaf->leaf_size = 1.0f;
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 5.0f;
    wp.sugar_diffusion_iterations = 20;

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
    leaf->leaf_size = 0.3f;
    plant.seed_mut()->add_child(leaf);

    float area = 0.3f * 0.3f;
    float expected = area * g.sugar_storage_density_leaf;
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
