#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/plant.h"
#include "engine/node/tissues/leaf.h"
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
    stem->local().chemical(ChemicalID::Sugar) = 0.0f;

    // Place far from other nodes so spatial diffusion doesn't interfere
    stem->position = glm::vec3(100.0f, 100.0f, 100.0f);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    REQUIRE_THAT(stem->local().chemical(ChemicalID::Ethylene), WithinAbs(g.ethylene_starvation_rate, 0.1));
}

TEST_CASE("Fed node produces no starvation ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);
    stem->local().chemical(ChemicalID::Sugar) = 5.0f;
    stem->position = glm::vec3(100.0f, 100.0f, 100.0f);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    REQUIRE(stem->local().chemical(ChemicalID::Ethylene) < 0.01f);
}

TEST_CASE("Shaded leaf produces ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->local().chemical(ChemicalID::Sugar) = 1.0f;
    leaf->as_leaf()->light_exposure = 0.1f;
    leaf->position = glm::vec3(100.0f, 100.0f, 100.0f);
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    float expected = g.ethylene_shade_rate * (1.0f - 0.1f);
    REQUIRE(leaf->local().chemical(ChemicalID::Ethylene) > expected * 0.5f);
}

TEST_CASE("Well-lit leaf produces no shade ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->local().chemical(ChemicalID::Sugar) = 1.0f;
    leaf->as_leaf()->light_exposure = 0.8f;
    leaf->position = glm::vec3(100.0f, 100.0f, 100.0f);
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    REQUIRE(leaf->local().chemical(ChemicalID::Ethylene) < 0.01f);
}

TEST_CASE("Old leaf produces age ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->local().chemical(ChemicalID::Sugar) = 1.0f;
    leaf->as_leaf()->light_exposure = 1.0f;
    leaf->age = g.ethylene_age_onset + 360;
    leaf->position = glm::vec3(100.0f, 100.0f, 100.0f);
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    float expected = g.ethylene_age_rate * 360.0f / static_cast<float>(g.ethylene_age_onset);
    REQUIRE(leaf->local().chemical(ChemicalID::Ethylene) > expected * 0.5f);
}

TEST_CASE("Crowded nodes produce ethylene", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    glm::vec3 center(5.0f, 5.0f, 5.0f);
    Node* target = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
    target->local().chemical(ChemicalID::Sugar) = 1.0f;
    target->position = center;
    plant.seed_mut()->add_child(target);

    for (int i = 0; i < 5; i++) {
        Node* n = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
        n->local().chemical(ChemicalID::Sugar) = 1.0f;
        n->position = center + glm::vec3(0.1f * i, 0.0f, 0.0f);
        plant.seed_mut()->add_child(n);
    }


    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    REQUIRE(target->local().chemical(ChemicalID::Ethylene) > 0.1f);
}

TEST_CASE("Spatial diffusion spreads ethylene to nearby nodes", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* source = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
    source->local().chemical(ChemicalID::Sugar) = 0.0f;
    source->position = glm::vec3(10.0f, 10.0f, 10.0f);
    plant.seed_mut()->add_child(source);

    Node* nearby = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
    nearby->local().chemical(ChemicalID::Sugar) = 1.0f;
    nearby->position = glm::vec3(10.5f, 10.0f, 10.0f);
    plant.seed_mut()->add_child(nearby);

    Node* far = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
    far->local().chemical(ChemicalID::Sugar) = 1.0f;
    far->position = glm::vec3(15.0f, 10.0f, 10.0f);
    plant.seed_mut()->add_child(far);

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    REQUIRE(nearby->local().chemical(ChemicalID::Ethylene) > 0.0f);
    REQUIRE(far->local().chemical(ChemicalID::Ethylene) < 0.01f);
}

TEST_CASE("Ethylene resets to zero before recomputing", "[ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    plant.seed_mut()->local().chemical(ChemicalID::Ethylene) = 999.0f;
    plant.seed_mut()->local().chemical(ChemicalID::Sugar) = 1.0f;
    plant.seed_mut()->position = glm::vec3(100.0f, 100.0f, 100.0f);

    plant.for_each_node_mut([](Node& n) {
        n.local().chemical(ChemicalID::Sugar) = 1.0f;
    });

    WorldParams wp = default_world_params();
    compute_ethylene(plant, wp);

    REQUIRE(plant.seed()->local().chemical(ChemicalID::Ethylene) < 999.0f);
}

TEST_CASE("Leaf above ethylene threshold begins senescence", "[ethylene][abscission]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->local().chemical(ChemicalID::Ethylene) = g.ethylene_abscission_threshold + 0.1f;
    leaf->local().chemical(ChemicalID::Sugar) = 1.0f;
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    leaf->tick(plant, wp);

    REQUIRE(leaf->as_leaf()->senescence_ticks > 0);
}

TEST_CASE("Leaf below ethylene threshold stays healthy", "[ethylene][abscission]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->local().chemical(ChemicalID::Ethylene) = g.ethylene_abscission_threshold * 0.5f;
    leaf->local().chemical(ChemicalID::Sugar) = 1.0f;
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    leaf->tick(plant, wp);

    REQUIRE(leaf->as_leaf()->senescence_ticks == 0);
}

TEST_CASE("Senescing leaf is removed after senescence_duration", "[ethylene][abscission]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->local().chemical(ChemicalID::Ethylene) = g.ethylene_abscission_threshold + 0.1f;
    leaf->as_leaf()->senescence_ticks = g.senescence_duration - 1; // almost done
    leaf->local().chemical(ChemicalID::Sugar) = 1.0f;
    plant.seed_mut()->add_child(leaf);

    uint32_t count_before = plant.node_count();
    WorldParams wp = default_world_params();
    leaf->tick(plant, wp);
    plant.flush_removals();

    REQUIRE(plant.node_count() < count_before);
}

// Disabled: compute_ethylene is currently disabled in Plant::tick() (crowding kills leaves instantly)
