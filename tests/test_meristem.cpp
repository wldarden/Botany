// tests/test_meristem.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/plant.h"
#include "engine/meristem.h"

using namespace botany;
using Catch::Matchers::WithinAbs;

TEST_CASE("Shoot apical meristem extends node position upward", "[meristem]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::STEM) { shoot = c; break; }
    }
    float y_before = shoot->position.y;

    plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
    tick_meristems(plant);

    REQUIRE(shoot->position.y > y_before);
    REQUIRE_THAT(shoot->position.y - y_before, WithinAbs(g.growth_rate, 0.001f));
}

TEST_CASE("Secondary growth thickens interior nodes, not tips", "[meristem]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::STEM) { shoot = c; break; }
    }
    float seed_r_before = seed->radius;
    float shoot_r_before = shoot->radius;

    plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
    tick_meristems(plant);

    // Seed (interior node, no active meristem) should thicken
    REQUIRE(seed->radius > seed_r_before);
    // Shoot tip (active apical meristem) should NOT thicken
    REQUIRE(shoot->radius == shoot_r_before);
}

TEST_CASE("Chain growth spawns axillary node and LEAF child on interior node", "[meristem]") {
    Genome g = default_genome();
    g.growth_rate = 0.5f;
    g.max_internode_length = 0.6f;
    Plant plant(g, glm::vec3(0.0f));

    uint32_t initial_count = plant.node_count();

    // Tick until chain growth fires (distance exceeds max_internode_length)
    for (int i = 0; i < 2; i++) {
        plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
        tick_meristems(plant);
    }

    // Should have spawned an axillary node with a dormant meristem and a LEAF node
    REQUIRE(plant.node_count() > initial_count);

    bool found_axillary = false;
    bool found_leaf = false;
    plant.for_each_node([&](const Node& n) {
        if (n.meristem && n.meristem->type() == MeristemType::AXILLARY) {
            found_axillary = true;
            REQUIRE(n.meristem->active == false);
        }
        if (n.type == NodeType::LEAF) {
            found_leaf = true;
            REQUIRE(n.leaf_size == g.leaf_size);
        }
    });
    REQUIRE(found_axillary);
    REQUIRE(found_leaf);
}

TEST_CASE("Interior STEM nodes have at most 3 children", "[meristem]") {
    Genome g = default_genome();
    g.growth_rate = 0.2f;
    g.max_internode_length = 0.3f;
    Plant plant(g, glm::vec3(0.0f));

    // Run enough ticks to create several interior nodes
    for (int i = 0; i < 20; i++) {
        plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
        tick_meristems(plant);
    }

    // Every STEM node should have at most 3 children: continuation tip + axillary + LEAF
    plant.for_each_node([&](const Node& n) {
        if (n.type == NodeType::STEM) {
            REQUIRE(n.children.size() <= 3);
        }
    });
}

TEST_CASE("Root apical meristem extends node position downward", "[meristem]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    const Node* seed = plant.seed();
    const Node* root = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::ROOT) { root = c; break; }
    }
    float y_before = root->position.y;

    plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
    tick_meristems(plant);

    REQUIRE(root->position.y < y_before);
}

TEST_CASE("Root chain growth spawns root axillary on interior node", "[meristem]") {
    Genome g = default_genome();
    g.root_growth_rate = 0.5f;
    g.root_max_internode_length = 0.6f;
    Plant plant(g, glm::vec3(0.0f));

    for (int i = 0; i < 2; i++) {
        plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
        tick_meristems(plant);
    }

    bool found_root_axillary = false;
    plant.for_each_node([&](const Node& n) {
        if (n.meristem && n.meristem->type() == MeristemType::ROOT_AXILLARY) {
            found_root_axillary = true;
            REQUIRE(n.meristem->active == false);
        }
    });
    REQUIRE(found_root_axillary);
}

TEST_CASE("Chain growth: apical meristem transfers to new node when internode too long", "[meristem]") {
    Genome g = default_genome();
    g.growth_rate = 0.5f;
    g.max_internode_length = 0.6f;
    Plant plant(g, glm::vec3(0.0f));

    // After 2 ticks the shoot moves 1.0 from its start, exceeding max_internode_length of 0.6
    plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
    tick_meristems(plant);
    plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
    tick_meristems(plant);

    // Should still have exactly one active apical meristem
    int apical_count = 0;
    plant.for_each_node([&](const Node& n) {
        if (n.meristem && n.meristem->type() == MeristemType::APICAL && n.meristem->active) {
            apical_count++;
        }
    });
    REQUIRE(apical_count == 1);

    // Node count should have increased (new chain + axillary nodes created)
    REQUIRE(plant.node_count() > 3);
}

TEST_CASE("Axillary meristem activates when auxin low and cytokinin high", "[meristem]") {
    Genome g = default_genome();
    g.auxin_threshold = 1.0f;
    g.cytokinin_threshold = 0.0f; // cytokinin just needs to be > 0
    g.growth_rate = 0.5f;
    g.max_internode_length = 0.4f; // chain growth fires quickly to create axillary
    Plant plant(g, glm::vec3(0.0f));

    // Tick once — chain growth fires, creating an interior node with axillary
    plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
    tick_meristems(plant);

    // Manually set hormones on the axillary node to trigger activation
    Node* axillary_node = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.meristem && n.meristem->type() == MeristemType::AXILLARY) {
            axillary_node = &n;
        }
    });
    REQUIRE(axillary_node != nullptr);
    axillary_node->auxin = 0.1f;      // below threshold of 1.0
    axillary_node->cytokinin = 0.5f;   // above threshold of 0.0

    plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
    tick_meristems(plant);

    // Should have converted to APICAL and be active
    REQUIRE(axillary_node->meristem->type() == MeristemType::APICAL);
    REQUIRE(axillary_node->meristem->active == true);
}

TEST_CASE("Axillary meristem stays dormant when auxin is high", "[meristem]") {
    Genome g = default_genome();
    g.auxin_threshold = 1.0f;
    g.growth_rate = 0.5f;
    g.max_internode_length = 0.4f;
    Plant plant(g, glm::vec3(0.0f));

    plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
    tick_meristems(plant);

    Node* axillary_node = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.meristem && n.meristem->type() == MeristemType::AXILLARY) {
            axillary_node = &n;
        }
    });
    REQUIRE(axillary_node != nullptr);
    // Shoot axillaries sense auxin on their parent stem node
    axillary_node->auxin = 5.0f; // way above threshold
    if (axillary_node->parent) {
        axillary_node->parent->auxin = 5.0f;
    }

    plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
    tick_meristems(plant);

    REQUIRE(axillary_node->meristem->type() == MeristemType::AXILLARY);
    REQUIRE(axillary_node->meristem->active == false);
}

TEST_CASE("Node age increments each tick", "[meristem]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
    tick_meristems(plant);
    plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
    tick_meristems(plant);

    const Node* seed = plant.seed();
    REQUIRE(seed->age == 2);
}

TEST_CASE("Shoot apical meristem does not grow without sugar", "[meristem][sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Find shoot tip and record position
    Node* shoot_tip = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.meristem && n.meristem->type() == MeristemType::APICAL && n.meristem->active) {
            shoot_tip = &n;
        }
    });
    REQUIRE(shoot_tip != nullptr);

    // Zero sugar — should not grow
    shoot_tip->sugar = 0.0f;
    glm::vec3 pos_before = shoot_tip->position;
    tick_meristems(plant);
    REQUIRE(shoot_tip->position.y == pos_before.y);
}

TEST_CASE("Shoot apical meristem grows and deducts sugar", "[meristem][sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* shoot_tip = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.meristem && n.meristem->type() == MeristemType::APICAL && n.meristem->active) {
            shoot_tip = &n;
        }
    });
    REQUIRE(shoot_tip != nullptr);

    shoot_tip->sugar = 100.0f;
    float sugar_before = shoot_tip->sugar;
    tick_meristems(plant);
    REQUIRE(shoot_tip->sugar < sugar_before);
}

TEST_CASE("Shoot axillary does not activate without sugar", "[meristem][sugar]") {
    Genome g = default_genome();
    g.growth_rate = 0.5f;
    g.max_internode_length = 0.6f;
    Plant plant(g, glm::vec3(0.0f));

    // Run until axillary buds exist
    // Give all nodes plenty of sugar for growth
    for (int i = 0; i < 5; i++) {
        plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
        tick_meristems(plant);
    }

    // Find a dormant axillary and zero its sugar
    Node* axillary_node = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.meristem && n.meristem->type() == MeristemType::AXILLARY) {
            axillary_node = &n;
        }
    });

    if (axillary_node) {
        axillary_node->sugar = 0.0f;
        // Set parent auxin low so hormone condition passes
        if (axillary_node->parent) {
            axillary_node->parent->auxin = 0.0f;
        }
        tick_meristems(plant);
        // Should still be axillary (not enough sugar to activate)
        REQUIRE(axillary_node->meristem != nullptr);
        REQUIRE(axillary_node->meristem->type() == MeristemType::AXILLARY);
    }
}

TEST_CASE("Thickening does not occur without sugar", "[meristem][sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    float radius_before = seed->radius;
    seed->sugar = 0.0f;

    tick_meristems(plant);
    REQUIRE(seed->radius == radius_before);
}

TEST_CASE("Thickening deducts sugar", "[meristem][sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    seed->sugar = 100.0f;
    float sugar_before = seed->sugar;

    tick_meristems(plant);
    REQUIRE(seed->sugar < sugar_before);
    REQUIRE(seed->radius > g.initial_radius);
}
