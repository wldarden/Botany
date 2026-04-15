// tests/test_plant.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/plant.h"
#include "engine/node/tissues/leaf.h"

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

    SECTION("one child is a APICAL meristem node") {
        const Node* seed = plant.seed();
        bool found_shoot = false;
        for (const Node* child : seed->children) {
            if (child->type == NodeType::APICAL) {
                found_shoot = true;
                REQUIRE(child->is_meristem());
                REQUIRE(child->position.y >= 0.0f);
            }
        }
        REQUIRE(found_shoot);
    }

    SECTION("one child is a ROOT_APICAL meristem node") {
        const Node* seed = plant.seed();
        bool found_root = false;
        for (const Node* child : seed->children) {
            if (child->type == NodeType::ROOT_APICAL) {
                found_root = true;
                REQUIRE(child->is_meristem());
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

TEST_CASE("Plant can create LEAF nodes", "[plant]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 1.0f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.3f;
    REQUIRE(leaf->type == NodeType::LEAF);
    REQUIRE(leaf->as_leaf()->leaf_size == 0.3f);
}

TEST_CASE("Plant::remove_subtree removes node and descendants", "[plant]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* child = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    Node* grandchild = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 2.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(child);
    child->add_child(grandchild);

    uint32_t id_child = child->id;
    uint32_t count_before = plant.node_count();

    plant.remove_subtree(child);

    REQUIRE(plant.node_count() == count_before - 2);

    // Child should no longer be in seed's children
    bool found = false;
    for (const Node* c : plant.seed()->children) {
        if (c->id == id_child) found = true;
    }
    REQUIRE_FALSE(found);
}
