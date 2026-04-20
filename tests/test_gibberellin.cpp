#include <catch2/catch_test_macros.hpp>
#include "engine/plant.h"
#include "engine/node/tissues/leaf.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/chemical/chemical.h"

using namespace botany;

TEST_CASE("Young leaf produces GA on itself", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->age = 10;
    leaf->local().chemical(ChemicalID::Sugar) = 1.0f;
    stem->add_child(leaf);

    WorldParams wp = default_world_params();
    leaf->tick(plant, wp);

    REQUIRE(leaf->local().chemical(ChemicalID::Gibberellin) > 0.0f);
}

TEST_CASE("GA spreads to parent via transport", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);
    stem->local().chemical(ChemicalID::Sugar) = 5.0f;

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->age = 10;
    leaf->local().chemical(ChemicalID::Sugar) = 50.0f;
    stem->add_child(leaf);

    WorldParams wp = default_world_params();

    // Tick the leaf to produce GA, then tick stem to transport it
    leaf->tick(plant, wp);
    stem->tick(plant, wp);  // stem transports with leaf child

    REQUIRE(stem->local().chemical(ChemicalID::Gibberellin) > 0.0f);
}

TEST_CASE("Old leaf produces no GA", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->age = g.ga_leaf_age_max + 1;
    leaf->local().chemical(ChemicalID::Sugar) = 1.0f;
    stem->add_child(leaf);

    WorldParams wp = default_world_params();
    leaf->tick(plant, wp);

    // No GA production (age exceeds max), and no prior GA in system
    REQUIRE(leaf->local().chemical(ChemicalID::Gibberellin) < 1e-6f);
}

TEST_CASE("GA decays without production", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);
    stem->local().chemical(ChemicalID::Gibberellin) = 1.0f;
    stem->local().chemical(ChemicalID::Sugar) = 5.0f;

    WorldParams wp = default_world_params();
    for (int i = 0; i < 20; i++) {
        stem->tick(plant, wp);
    }

    REQUIRE(stem->local().chemical(ChemicalID::Gibberellin) < 0.1f);
}

TEST_CASE("Senescing leaf produces no GA", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->age = 10;
    leaf->as_leaf()->senescence_ticks = 1;
    leaf->local().chemical(ChemicalID::Sugar) = 1.0f;
    stem->add_child(leaf);

    WorldParams wp = default_world_params();
    leaf->tick(plant, wp);

    REQUIRE(leaf->local().chemical(ChemicalID::Gibberellin) < 1e-6f);
}
