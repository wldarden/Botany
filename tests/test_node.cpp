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

TEST_CASE("Canalization: biased child gets larger sugar share", "[canalization]") {
    Genome g = default_genome();
    g.canalization_weight = 1.0f;
    Plant plant(g, glm::vec3(0.0f));

    // Use small parent sugar so children don't hit headroom caps.
    // Both children are identical STEM nodes with cap ~0.39g.
    // Parent sugar = 0.3g → budget < headroom, so bias split shows through.
    Node* parent_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* childA = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* childB = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    plant.seed_mut()->add_child(parent_node);
    parent_node->add_child(childA);
    parent_node->add_child(childB);

    parent_node->chemical(ChemicalID::Sugar) = 0.3f;
    childA->chemical(ChemicalID::Sugar) = 0.0f;
    childB->chemical(ChemicalID::Sugar) = 0.0f;

    // bias_mult for childA = 1 + 1*(1+0) = 2.0, for childB = 1.0
    // childA should get 2/3 of the budget, childB 1/3
    parent_node->auxin_flow_bias[childA] = 1.0f;
    parent_node->auxin_flow_bias[childB] = 0.0f;

    parent_node->transport_with_children(g);
    // Flush received buffers (normally done by tick())
    for (Node* c : parent_node->children) {
        for (auto& [id, amt] : c->transport_received) c->chemical(id) += amt;
        c->transport_received.clear();
    }

    REQUIRE(childA->chemical(ChemicalID::Sugar) > childB->chemical(ChemicalID::Sugar) * 1.5f);
}

TEST_CASE("Canalization: zero canalization_weight disables bias", "[canalization]") {
    Genome g = default_genome();
    g.canalization_weight = 0.0f;  // disabled
    Plant plant(g, glm::vec3(0.0f));

    // Same small-sugar setup as the bias test
    Node* parent_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* childA = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* childB = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    plant.seed_mut()->add_child(parent_node);
    parent_node->add_child(childA);
    parent_node->add_child(childB);

    parent_node->chemical(ChemicalID::Sugar) = 0.3f;
    childA->chemical(ChemicalID::Sugar) = 0.0f;
    childB->chemical(ChemicalID::Sugar) = 0.0f;

    // Big biases — but weight=0 should make them irrelevant
    parent_node->auxin_flow_bias[childA] = 5.0f;
    parent_node->auxin_flow_bias[childB] = 0.0f;

    parent_node->transport_with_children(g);
    // Flush received buffers (normally done by tick())
    for (Node* c : parent_node->children) {
        for (auto& [id, amt] : c->transport_received) c->chemical(id) += amt;
        c->transport_received.clear();
    }

    // Both should get roughly equal sugar
    float ratio = childA->chemical(ChemicalID::Sugar) /
                  std::max(childB->chemical(ChemicalID::Sugar), 1e-8f);
    REQUIRE(ratio > 0.8f);
    REQUIRE(ratio < 1.2f);
}

TEST_CASE("Canalization: transport with biases doesn't break child-to-parent flow", "[canalization]") {
    Genome g = default_genome();
    g.canalization_weight = 1.0f;
    Plant plant(g, glm::vec3(0.0f));

    // Verify that setting biases doesn't break normal child→parent transport.
    // Children have auxin, parent has none — auxin flows basipetally.
    Node* parent_node = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* childA = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    Node* childB = plant.create_node(NodeType::STEM, glm::vec3(0, 0.1f, 0), 0.05f);
    plant.seed_mut()->add_child(parent_node);
    parent_node->add_child(childA);
    parent_node->add_child(childB);

    childA->chemical(ChemicalID::Auxin) = 5.0f;
    childB->chemical(ChemicalID::Auxin) = 5.0f;
    parent_node->chemical(ChemicalID::Auxin) = 0.0f;

    parent_node->auxin_flow_bias[childA] = 2.0f;
    parent_node->auxin_flow_bias[childB] = 0.0f;

    parent_node->transport_with_children(g);

    // Parent should have received auxin from both children
    REQUIRE(parent_node->chemical(ChemicalID::Auxin) > 0.0f);
    // Both children should have given some auxin
    REQUIRE(childA->chemical(ChemicalID::Auxin) < 5.0f);
    REQUIRE(childB->chemical(ChemicalID::Auxin) < 5.0f);
}

TEST_CASE("Canalization: transient bias builds from auxin flux", "[canalization]") {
    Genome g = default_genome();
    g.transient_gain = 2.0f;
    g.transient_rate = 0.5f;  // fast for testing
    g.canalization_weight = 1.0f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world_params();

    // Give shoot apical lots of sugar so it produces auxin and grows
    for (int i = 0; i < 10; i++) {
        plant.for_each_node_mut([](Node& n) {
            n.chemical(ChemicalID::Sugar) = 100.0f;
            n.chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant.tick(world);
    }

    // The seed should have non-zero transient bias toward the shoot apical child
    bool found_nonzero_bias = false;
    Node* seed = plant.seed_mut();
    for (auto& [child, bias] : seed->auxin_flow_bias) {
        if (bias > 0.01f) found_nonzero_bias = true;
    }
    REQUIRE(found_nonzero_bias);
}

TEST_CASE("Canalization: transient bias decays without flux", "[canalization]") {
    Genome g = default_genome();
    g.transient_gain = 2.0f;
    g.transient_rate = 0.5f;
    g.apical_auxin_baseline = 0.0f;  // no auxin production
    g.canalization_weight = 1.0f;
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    Node* child = nullptr;
    for (Node* c : seed->children) {
        if (c->type == NodeType::APICAL) { child = c; break; }
    }
    REQUIRE(child != nullptr);

    // Manually set a transient bias
    seed->auxin_flow_bias[child] = 1.0f;

    // Tick with no auxin production — bias should decay toward 0
    WorldParams world = default_world_params();
    for (int i = 0; i < 10; i++) {
        plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; });
        plant.tick(world);
    }

    // Bias should have decayed significantly (target=0, rate=0.5)
    REQUIRE(seed->auxin_flow_bias[child] < 0.1f);
}

TEST_CASE("Canalization: structural bias builds with sustained flux", "[canalization]") {
    Genome g = default_genome();
    g.structural_threshold = 0.01f;  // low for testing
    g.structural_growth_rate = 0.1f; // fast for testing
    g.structural_max = 5.0f;
    g.canalization_weight = 1.0f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world_params();

    for (int i = 0; i < 10; i++) {
        plant.for_each_node_mut([](Node& n) {
            n.chemical(ChemicalID::Sugar) = 100.0f;
            n.chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant.tick(world);
    }

    bool found_structural = false;
    Node* seed = plant.seed_mut();
    for (auto& [child, bias] : seed->structural_flow_bias) {
        if (bias > 0.01f) found_structural = true;
    }
    REQUIRE(found_structural);
}

TEST_CASE("Canalization: structural bias does NOT decay", "[canalization]") {
    Genome g = default_genome();
    g.structural_threshold = 0.01f;
    g.structural_growth_rate = 0.1f;
    g.structural_max = 5.0f;
    g.apical_auxin_baseline = 0.0f;  // disable auxin production for decay test
    g.canalization_weight = 1.0f;
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    Node* shoot = nullptr;
    for (Node* c : seed->children) {
        if (c->type == NodeType::APICAL) { shoot = c; break; }
    }
    REQUIRE(shoot != nullptr);

    // Manually set structural bias
    seed->structural_flow_bias[shoot] = 1.5f;
    float bias_before = seed->structural_flow_bias[shoot];

    // Tick with zero auxin — no flux, structural should NOT decay
    WorldParams world = default_world_params();
    for (int i = 0; i < 20; i++) {
        plant.for_each_node_mut([](Node& n) {
            n.chemical(ChemicalID::Sugar) = 100.0f;
            n.chemical(ChemicalID::Auxin) = 0.0f;
        });
        plant.tick(world);
    }

    REQUIRE(seed->structural_flow_bias[shoot] >= bias_before);
}

TEST_CASE("Canalization: structural bias capped at structural_max", "[canalization]") {
    Genome g = default_genome();
    g.structural_max = 1.0f;
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    Node* child = nullptr;
    for (Node* c : seed->children) { child = c; break; }
    REQUIRE(child != nullptr);

    // Set bias above max
    seed->structural_flow_bias[child] = 5.0f;

    // update_canalization should clamp to max
    seed->update_canalization(plant.genome());

    REQUIRE(seed->structural_flow_bias[child] <= g.structural_max);
}
