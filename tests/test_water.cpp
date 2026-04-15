#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <glm/geometric.hpp>
#include "engine/plant.h"
#include "engine/node/tissues/leaf.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/sugar.h"
#include "engine/chemical/chemical.h"

using namespace botany;
using Catch::Matchers::WithinAbs;

// === Water cap tests ===

TEST_CASE("water_cap scales with stem volume", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.1f);
    plant.seed_mut()->add_child(stem);

    float volume = 3.14159f * 0.1f * 0.1f * 1.0f;
    float expected = volume * g.water_storage_density_stem;
    REQUIRE_THAT(water_cap(*stem, g), WithinAbs(expected, 1e-4));
}

TEST_CASE("water_cap scales with leaf area", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    plant.seed_mut()->add_child(leaf);

    float area = 0.5f * 0.5f;
    float expected = area * g.water_storage_density_leaf;
    REQUIRE_THAT(water_cap(*leaf, g), WithinAbs(expected, 1e-6));
}

TEST_CASE("water_cap for root scales with volume", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* root = plant.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.8f, 0.0f), 0.025f);
    plant.seed_mut()->add_child(root);

    float length = 0.8f;
    float volume = 3.14159f * 0.025f * 0.025f * length;
    float expected = volume * g.water_storage_density_stem;
    REQUIRE_THAT(water_cap(*root, g), WithinAbs(expected, 1e-5));
}

TEST_CASE("water_cap for meristem is fixed", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* apical = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::APICAL) apical = &n;
    });
    REQUIRE(apical != nullptr);
    REQUIRE_THAT(water_cap(*apical, g), WithinAbs(g.water_cap_meristem, 1e-6));
}

TEST_CASE("water_cap returns minimum for tiny nodes", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* tiny = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.001f, 0.0f), 0.001f);
    plant.seed_mut()->add_child(tiny);

    REQUIRE_THAT(water_cap(*tiny, g), WithinAbs(g.sugar_cap_minimum, 1e-6));
}

// === Water transport tests ===

TEST_CASE("Water diffuses from high to low concentration", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* child = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.1f, 0.0f), 0.02f);
    plant.seed_mut()->add_child(child);

    plant.seed_mut()->chemical(ChemicalID::Water) = 50.0f;
    child->chemical(ChemicalID::Water) = 0.0f;

    plant.seed_mut()->transport_with_children(g);

    REQUIRE(child->chemical(ChemicalID::Water) > 0.0f);
    REQUIRE(plant.seed()->chemical(ChemicalID::Water) < 50.0f);
}

TEST_CASE("Water transport conserves total water", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    for (int i = 0; i < 3; i++) {
        Node* child = plant.create_node(NodeType::STEM,
            glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
        plant.seed_mut()->add_child(child);
    }

    plant.seed_mut()->chemical(ChemicalID::Water) = 20.0f;
    plant.seed_mut()->children[0]->chemical(ChemicalID::Water) = 5.0f;

    float total_before = 0.0f;
    plant.for_each_node([&](const Node& n) { total_before += n.chemical(ChemicalID::Water); });

    for (int i = 0; i < 10; i++) {
        plant.for_each_node_mut([&](Node& n) {
            n.transport_with_children(g);
        });
    }

    float total_after = 0.0f;
    plant.for_each_node([&](const Node& n) { total_after += n.chemical(ChemicalID::Water); });

    REQUIRE_THAT(total_after, WithinAbs(total_before, 1e-4));
}
