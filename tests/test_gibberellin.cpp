#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/plant.h"
#include "engine/node/leaf_node.h"
#include "engine/genome.h"
#include "engine/gibberellin.h"

using namespace botany;
using Catch::Matchers::WithinAbs;

TEST_CASE("Young leaf produces GA on parent node", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Create: seed -> stem -> leaf
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->age = 10; // young
    stem->add_child(leaf);

    compute_gibberellin(plant);

    float expected = leaf->as_leaf()->leaf_size * g.ga_production_rate;
    REQUIRE_THAT(stem->chemical(ChemicalID::Gibberellin), WithinAbs(expected, 1e-6));
}

TEST_CASE("Old leaf produces no GA", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->age = g.ga_leaf_age_max + 1; // too old
    stem->add_child(leaf);

    compute_gibberellin(plant);

    REQUIRE(stem->chemical(ChemicalID::Gibberellin) == 0.0f);
}

TEST_CASE("GA reaches grandparent at reduced fraction", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // seed -> grandparent -> parent -> leaf
    Node* grandparent = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(grandparent);

    Node* parent = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    grandparent->add_child(parent);

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->age = 10;
    parent->add_child(leaf);

    compute_gibberellin(plant);

    float parent_expected = leaf->as_leaf()->leaf_size * g.ga_production_rate;
    float grandparent_expected = leaf->as_leaf()->leaf_size * g.ga_production_rate * 0.3f;

    REQUIRE_THAT(parent->chemical(ChemicalID::Gibberellin), WithinAbs(parent_expected, 1e-6));
    REQUIRE_THAT(grandparent->chemical(ChemicalID::Gibberellin), WithinAbs(grandparent_expected, 1e-6));
}

TEST_CASE("GA does not spread beyond grandparent", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // seed -> great_grandparent -> grandparent -> parent -> leaf
    Node* ggp = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(ggp);

    Node* gp = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    ggp->add_child(gp);

    Node* parent = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    gp->add_child(parent);

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.1f, 0.1f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.2f;
    leaf->age = 10;
    parent->add_child(leaf);

    compute_gibberellin(plant);

    // Great-grandparent and seed should have no GA
    REQUIRE(ggp->chemical(ChemicalID::Gibberellin) == 0.0f);
    REQUIRE(plant.seed()->chemical(ChemicalID::Gibberellin) == 0.0f);
}

TEST_CASE("GA resets to zero before recomputing", "[gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(stem);

    // Manually set GA to something nonzero
    stem->chemical(ChemicalID::Gibberellin) = 999.0f; stem->gibberellin = 999.0f;

    // No leaves — compute should reset to 0
    compute_gibberellin(plant);

    REQUIRE(stem->chemical(ChemicalID::Gibberellin) == 0.0f);
}
