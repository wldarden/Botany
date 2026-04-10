#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/node.h"
#include "engine/meristem_types.h"

using namespace botany;

TEST_CASE("Node creation with default values", "[node]") {
    Node node(1, NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);

    REQUIRE(node.id == 1);
    REQUIRE(node.type == NodeType::STEM);
    REQUIRE(node.position.y == 1.0f);
    REQUIRE(node.radius == 0.05f);
    REQUIRE(node.parent == nullptr);
    REQUIRE(node.children.empty());
    REQUIRE(node.age == 0);
    REQUIRE(node.auxin == 0.0f);
    REQUIRE(node.cytokinin == 0.0f);
    REQUIRE(node.meristem == nullptr);
    REQUIRE(node.sugar == 0.0f);
    REQUIRE(node.leaf_size == 0.0f);
}

TEST_CASE("add_child establishes parent-child relationship", "[node]") {
    Node parent(1, NodeType::STEM, glm::vec3(0.0f), 0.05f);
    Node child(2, NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.04f);

    parent.add_child(&child);

    REQUIRE(parent.children.size() == 1);
    REQUIRE(parent.children[0] == &child);
    REQUIRE(child.parent == &parent);
}

TEST_CASE("Node can have a meristem attached", "[node]") {
    Node node(1, NodeType::STEM, glm::vec3(0.0f), 0.05f);
    ShootApicalMeristem m;
    node.meristem = &m;

    REQUIRE(node.meristem != nullptr);
    REQUIRE(node.meristem->type() == MeristemType::APICAL);
    REQUIRE(node.meristem->active == true);
}

TEST_CASE("LEAF node stores leaf_size", "[node]") {
    Node node(1, NodeType::LEAF, glm::vec3(0.0f), 0.0f);
    node.leaf_size = 0.3f;
    REQUIRE(node.type == NodeType::LEAF);
    REQUIRE(node.leaf_size == 0.3f);
    REQUIRE(node.sugar == 0.0f);
}

TEST_CASE("Meristem types cover all four variants", "[node]") {
    REQUIRE(static_cast<int>(MeristemType::APICAL) != static_cast<int>(MeristemType::AXILLARY));
    REQUIRE(static_cast<int>(MeristemType::ROOT_APICAL) != static_cast<int>(MeristemType::ROOT_AXILLARY));
}
