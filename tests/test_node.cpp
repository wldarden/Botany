#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/node/node.h"
#include "engine/node/stem_node.h"
#include "engine/node/leaf_node.h"
#include "engine/node/meristems/shoot_apical.h"

using namespace botany;

TEST_CASE("Node creation with default values", "[node]") {
    StemNode node(1, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);

    REQUIRE(node.id == 1);
    REQUIRE(node.type == NodeType::STEM);
    REQUIRE(node.position.y == 1.0f);
    REQUIRE(node.radius == 0.05f);
    REQUIRE(node.parent == nullptr);
    REQUIRE(node.children.empty());
    REQUIRE(node.age == 0);
    REQUIRE(node.chemical(ChemicalID::Auxin) == 0.0f);
    REQUIRE(node.chemical(ChemicalID::Cytokinin) == 0.0f);
    REQUIRE(node.is_meristem() == false);
    REQUIRE(node.chemical(ChemicalID::Sugar) == 0.0f);
}

TEST_CASE("add_child establishes parent-child relationship", "[node]") {
    StemNode parent(1, glm::vec3(0.0f), 0.05f);
    StemNode child(2, glm::vec3(0.0f, 1.0f, 0.0f), 0.04f);

    parent.add_child(&child);

    REQUIRE(parent.children.size() == 1);
    REQUIRE(parent.children[0] == &child);
    REQUIRE(child.parent == &parent);
}

TEST_CASE("ShootApicalNode has correct meristem properties", "[node]") {
    ShootApicalNode node(1, glm::vec3(0.0f), 0.05f);

    REQUIRE(node.type == NodeType::SHOOT_APICAL);
    REQUIRE(node.is_meristem() == true);
    REQUIRE(node.as_shoot_apical() != nullptr);
}

TEST_CASE("LEAF node stores leaf_size", "[node]") {
    LeafNode node(1, glm::vec3(0.0f), 0.0f);
    node.leaf_size = 0.3f;
    REQUIRE(node.type == NodeType::LEAF);
    REQUIRE(node.leaf_size == 0.3f);
    REQUIRE(node.chemical(ChemicalID::Sugar) == 0.0f);
}

TEST_CASE("NodeType covers all four meristem variants", "[node]") {
    REQUIRE(static_cast<int>(NodeType::SHOOT_APICAL) != static_cast<int>(NodeType::SHOOT_AXILLARY));
    REQUIRE(static_cast<int>(NodeType::ROOT_APICAL) != static_cast<int>(NodeType::ROOT_AXILLARY));
}
