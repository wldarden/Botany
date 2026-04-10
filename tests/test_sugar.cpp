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
