// tests/test_plant.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/plant.h"
#include "engine/meristem_types.h"

using namespace botany;

TEST_CASE("Plant seed initialization creates correct graph", "[plant]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f, 0.0f, 0.0f));

    SECTION("has a seed node at origin") {
        const Node* seed = plant.seed();
        REQUIRE(seed != nullptr);
        REQUIRE(seed->position == glm::vec3(0.0f));
        REQUIRE(seed->parent == nullptr);
    }

    SECTION("seed has two children: shoot and root") {
        const Node* seed = plant.seed();
        REQUIRE(seed->children.size() == 2);
    }

    SECTION("one child is a STEM with apical meristem") {
        const Node* seed = plant.seed();
        bool found_shoot = false;
        for (const Node* child : seed->children) {
            if (child->type == NodeType::STEM && child->meristem &&
                child->meristem->type() == MeristemType::APICAL) {
                found_shoot = true;
                REQUIRE(child->meristem->active == true);
                REQUIRE(child->position.y >= 0.0f);
            }
        }
        REQUIRE(found_shoot);
    }

    SECTION("one child is a ROOT with root apical meristem") {
        const Node* seed = plant.seed();
        bool found_root = false;
        for (const Node* child : seed->children) {
            if (child->type == NodeType::ROOT && child->meristem &&
                child->meristem->type() == MeristemType::ROOT_APICAL) {
                found_root = true;
                REQUIRE(child->meristem->active == true);
                REQUIRE(child->position.y <= 0.0f);
            }
        }
        REQUIRE(found_root);
    }

    SECTION("initial radius from genome") {
        const Node* seed = plant.seed();
        REQUIRE(seed->radius == g.initial_radius);
    }

    SECTION("node_count starts at 3") {
        REQUIRE(plant.node_count() == 3);
    }
}

TEST_CASE("Plant provides access to all nodes", "[plant]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    std::vector<const Node*> all;
    plant.for_each_node([&](const Node& n) {
        all.push_back(&n);
    });

    REQUIRE(all.size() == 3);
}

TEST_CASE("Plant can create new nodes", "[plant]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* parent = plant.seed_mut();
    Node* child = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 2.0f, 0.0f), g.initial_radius);
    parent->add_child(child);

    REQUIRE(plant.node_count() == 4);
    REQUIRE(child->parent == parent);
}

TEST_CASE("Plant can create meristems and leaves", "[plant]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* node = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    Meristem* m = plant.create_meristem<ShootAxillaryMeristem>();
    node->meristem = m;
    Leaf* l = plant.create_leaf(0.3f);
    node->leaf = l;

    REQUIRE(node->meristem->type() == MeristemType::AXILLARY);
    REQUIRE(node->meristem->active == false);
    REQUIRE(node->leaf->size == 0.3f);
}
