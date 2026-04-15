#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/node/node.h"
#include "engine/node/stem_node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"
#include "engine/plant.h"
#include "engine/world_params.h"

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

TEST_CASE("ApicalNode has correct meristem properties", "[node]") {
    ApicalNode node(1, glm::vec3(0.0f), 0.05f);

    REQUIRE(node.type == NodeType::APICAL);
    REQUIRE(node.is_meristem() == true);
    REQUIRE(node.as_apical() != nullptr);
}

TEST_CASE("LEAF node stores leaf_size", "[node]") {
    LeafNode node(1, glm::vec3(0.0f), 0.0f);
    node.leaf_size = 0.3f;
    REQUIRE(node.type == NodeType::LEAF);
    REQUIRE(node.leaf_size == 0.3f);
    REQUIRE(node.chemical(ChemicalID::Sugar) == 0.0f);
}

TEST_CASE("NodeType covers both meristem types", "[node]") {
    REQUIRE(static_cast<int>(NodeType::APICAL) != static_cast<int>(NodeType::ROOT_APICAL));
}

TEST_CASE("Canalization: replace_child transfers biases", "[canalization]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* parent_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* old_child = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* new_child = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);

    plant.seed_mut()->add_child(parent_node);
    parent_node->add_child(old_child);

    parent_node->auxin_flow_bias[old_child] = 0.5f;
    parent_node->structural_flow_bias[old_child] = 0.3f;

    parent_node->replace_child(old_child, new_child);

    REQUIRE(parent_node->auxin_flow_bias.count(new_child) == 1);
    REQUIRE(parent_node->auxin_flow_bias[new_child] == 0.5f);
    REQUIRE(parent_node->structural_flow_bias[new_child] == 0.3f);
    REQUIRE(parent_node->auxin_flow_bias.count(old_child) == 0);
    REQUIRE(parent_node->structural_flow_bias.count(old_child) == 0);
}

TEST_CASE("Canalization: die cleans up parent bias entries", "[canalization]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* parent_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* child_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);

    plant.seed_mut()->add_child(parent_node);
    parent_node->add_child(child_node);

    parent_node->auxin_flow_bias[child_node] = 0.7f;
    parent_node->structural_flow_bias[child_node] = 0.4f;

    child_node->die(plant);

    REQUIRE(parent_node->auxin_flow_bias.count(child_node) == 0);
    REQUIRE(parent_node->structural_flow_bias.count(child_node) == 0);
}

TEST_CASE("Canalization: get_bias_multiplier returns 1.0 with no bias", "[canalization]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* parent_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* child_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    plant.seed_mut()->add_child(parent_node);
    parent_node->add_child(child_node);

    REQUIRE(parent_node->get_bias_multiplier(child_node, g) == 1.0f);
}

TEST_CASE("Canalization: get_bias_multiplier scales with canalization_weight", "[canalization]") {
    Genome g = default_genome();
    g.canalization_weight = 2.0f;
    Plant plant(g, glm::vec3(0.0f));

    Node* parent_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* child_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    plant.seed_mut()->add_child(parent_node);
    parent_node->add_child(child_node);

    parent_node->auxin_flow_bias[child_node] = 0.5f;
    parent_node->structural_flow_bias[child_node] = 0.3f;

    float expected = 1.0f + 2.0f * (0.5f + 0.3f);
    REQUIRE(parent_node->get_bias_multiplier(child_node, g) == expected);
}
