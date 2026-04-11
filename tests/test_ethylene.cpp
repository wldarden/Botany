#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/ethylene.h"

using namespace botany;
using Catch::Matchers::WithinAbs;

TEST_CASE("Starvation produces ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);
    stem->sugar = 0.0f;

    // Place far from other nodes so spatial diffusion doesn't interfere
    stem->position = glm::vec3(100.0f, 100.0f, 100.0f);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    REQUIRE_THAT(stem->ethylene, WithinAbs(g.ethylene_starvation_rate, 0.1));
}

TEST_CASE("Fed node produces no starvation ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);
    stem->sugar = 5.0f;
    stem->position = glm::vec3(100.0f, 100.0f, 100.0f);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    REQUIRE(stem->ethylene < 0.01f);
}

TEST_CASE("Shaded leaf produces ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->sugar = 1.0f;
    leaf->as_leaf()->light_exposure = 0.1f;
    leaf->position = glm::vec3(100.0f, 100.0f, 100.0f);
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    float expected = g.ethylene_shade_rate * (1.0f - 0.1f);
    REQUIRE(leaf->ethylene > expected * 0.5f);
}

TEST_CASE("Well-lit leaf produces no shade ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->sugar = 1.0f;
    leaf->as_leaf()->light_exposure = 0.8f;
    leaf->position = glm::vec3(100.0f, 100.0f, 100.0f);
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    REQUIRE(leaf->ethylene < 0.01f);
}

TEST_CASE("Old leaf produces age ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->sugar = 1.0f;
    leaf->as_leaf()->light_exposure = 1.0f;
    leaf->age = g.ethylene_age_onset + 360;
    leaf->position = glm::vec3(100.0f, 100.0f, 100.0f);
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    float expected = g.ethylene_age_rate * 360.0f / static_cast<float>(g.ethylene_age_onset);
    REQUIRE(leaf->ethylene > expected * 0.5f);
}

TEST_CASE("Crowded nodes produce ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    glm::vec3 center(5.0f, 5.0f, 5.0f);
    Node* target = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
    target->sugar = 1.0f;
    target->position = center;
    plant.seed_mut()->add_child(target);

    for (int i = 0; i < 5; i++) {
        Node* n = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
        n->sugar = 1.0f;
        n->position = center + glm::vec3(0.1f * i, 0.0f, 0.0f);
        plant.seed_mut()->add_child(n);
    }

    plant.recompute_world_positions();

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    REQUIRE(target->ethylene > 0.1f);
}

TEST_CASE("Spatial diffusion spreads ethylene to nearby nodes", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* source = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
    source->sugar = 0.0f;
    source->position = glm::vec3(10.0f, 10.0f, 10.0f);
    plant.seed_mut()->add_child(source);

    Node* nearby = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
    nearby->sugar = 1.0f;
    nearby->position = glm::vec3(10.5f, 10.0f, 10.0f);
    plant.seed_mut()->add_child(nearby);

    Node* far = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
    far->sugar = 1.0f;
    far->position = glm::vec3(15.0f, 10.0f, 10.0f);
    plant.seed_mut()->add_child(far);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    REQUIRE(nearby->ethylene > 0.0f);
    REQUIRE(far->ethylene < 0.01f);
}

TEST_CASE("Ethylene resets to zero before recomputing", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    plant.seed_mut()->ethylene = 999.0f;
    plant.seed_mut()->sugar = 1.0f;
    plant.seed_mut()->position = glm::vec3(100.0f, 100.0f, 100.0f);

    plant.for_each_node_mut([](Node& n) {
        n.sugar = 1.0f;
    });

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    REQUIRE(plant.seed()->ethylene < 999.0f);
}

TEST_CASE("Leaf above ethylene threshold begins senescence", "[ethylene][abscission]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->ethylene = g.ethylene_abscission_threshold + 0.1f;
    plant.seed_mut()->add_child(leaf);

    process_abscission(plant);

    REQUIRE(leaf->as_leaf()->senescence_ticks > 0);
}

TEST_CASE("Leaf below ethylene threshold stays healthy", "[ethylene][abscission]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->ethylene = g.ethylene_abscission_threshold * 0.5f;
    plant.seed_mut()->add_child(leaf);

    process_abscission(plant);

    REQUIRE(leaf->as_leaf()->senescence_ticks == 0);
}

TEST_CASE("Senescing leaf is removed after senescence_duration", "[ethylene][abscission]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->ethylene = g.ethylene_abscission_threshold + 0.1f;
    leaf->as_leaf()->senescence_ticks = g.senescence_duration - 1; // almost done
    plant.seed_mut()->add_child(leaf);

    uint32_t count_before = plant.node_count();
    process_abscission(plant);

    REQUIRE(plant.node_count() < count_before);
}

TEST_CASE("Non-leaf nodes do not senesce", "[ethylene][abscission]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    stem->ethylene = g.ethylene_abscission_threshold + 10.0f; // very high
    plant.seed_mut()->add_child(stem);

    uint32_t count_before = plant.node_count();
    process_abscission(plant);

    // senescence_ticks only exists on LeafNode, so stem nodes can't senesce
    REQUIRE(plant.node_count() == count_before);
}

TEST_CASE("Self-thinning cascade prunes shaded interior leaves", "[ethylene][integration]") {
    Genome g = default_genome();
    // Make abscission fast for test
    g.senescence_duration = 5;
    g.ethylene_shade_rate = 1.0f;      // strong shade response
    g.ethylene_shade_threshold = 0.5f;
    g.ethylene_abscission_threshold = 0.3f; // easy to trigger
    // Disable crowding and diffusion so only direct shade signal matters
    g.ethylene_crowding_rate = 0.0f;
    g.ethylene_diffusion_radius = 0.05f; // very short-range — no cross-leaf spread

    Plant plant(g, glm::vec3(0.0f));

    // Build a small branch with shaded interior leaves and sunlit outer leaves.
    // seed -> stem -> inner_leaf (shaded) + outer_leaf (sunlit)
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);
    plant.seed_mut()->sugar = 100.0f; // plenty of sugar

    Node* inner_leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    inner_leaf->as_leaf()->leaf_size = 0.2f;
    inner_leaf->as_leaf()->light_exposure = 0.1f; // heavily shaded
    inner_leaf->position = glm::vec3(0.1f, 1.1f, 0.0f);
    inner_leaf->sugar = 1.0f; // fed — only shade triggers ethylene
    stem->add_child(inner_leaf);

    Node* outer_leaf = plant.create_node(NodeType::LEAF, glm::vec3(-0.1f, 0.1f, 0.0f), 0.0f);
    outer_leaf->as_leaf()->leaf_size = 0.2f;
    outer_leaf->as_leaf()->light_exposure = 0.9f; // well-lit
    outer_leaf->position = glm::vec3(-0.1f, 1.1f, 0.0f);
    outer_leaf->sugar = 1.0f; // fed — no starvation ethylene
    stem->add_child(outer_leaf);

    uint32_t inner_id = inner_leaf->id;
    uint32_t outer_id = outer_leaf->id;

    WorldParams wp = default_world_params();

    // Run ethylene + abscission for enough ticks to trigger senescence and removal
    for (int i = 0; i < 20; i++) {
        // Keep light_exposure fixed (simulate persistent shade)
        plant.for_each_node_mut([&](Node& n) {
            if (n.id == inner_id) n.as_leaf()->light_exposure = 0.1f;
            if (n.id == outer_id) n.as_leaf()->light_exposure = 0.9f;
        });
        compute_ethylene(plant, wp);
        process_abscission(plant);
    }

    // Inner (shaded) leaf should have been removed
    bool inner_exists = false;
    bool outer_exists = false;
    plant.for_each_node([&](const Node& n) {
        if (n.id == inner_id) inner_exists = true;
        if (n.id == outer_id) outer_exists = true;
    });

    REQUIRE_FALSE(inner_exists); // shaded leaf pruned
    REQUIRE(outer_exists);       // sunlit leaf survived
}
