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

    tick_meristems(plant);

    REQUIRE(shoot->position.y > y_before);
    REQUIRE_THAT(shoot->position.y - y_before, WithinAbs(g.growth_rate, 0.001f));
}

TEST_CASE("Shoot apical meristem thickens its node", "[meristem]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::STEM) { shoot = c; break; }
    }
    float r_before = shoot->radius;

    tick_meristems(plant);

    REQUIRE(shoot->radius > r_before);
}

TEST_CASE("Shoot apical spawns axillary node at internode_spacing", "[meristem]") {
    Genome g = default_genome();
    g.internode_spacing = 3;
    Plant plant(g, glm::vec3(0.0f));

    uint32_t initial_count = plant.node_count();

    // Tick 3 times to hit internode_spacing
    for (int i = 0; i < 3; i++) {
        tick_meristems(plant);
    }

    // Should have spawned an axillary node with a dormant meristem and a leaf
    REQUIRE(plant.node_count() > initial_count);

    // Find the new axillary node
    bool found_axillary = false;
    plant.for_each_node([&](const Node& n) {
        if (n.meristem && n.meristem->type == MeristemType::AXILLARY) {
            found_axillary = true;
            REQUIRE(n.meristem->active == false);
            REQUIRE(n.leaf != nullptr);
            REQUIRE(n.leaf->size == g.leaf_size);
        }
    });
    REQUIRE(found_axillary);
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

    tick_meristems(plant);

    REQUIRE(root->position.y < y_before);
}

TEST_CASE("Root apical spawns root axillary at root_internode_spacing", "[meristem]") {
    Genome g = default_genome();
    g.root_internode_spacing = 2;
    Plant plant(g, glm::vec3(0.0f));

    for (int i = 0; i < 2; i++) {
        tick_meristems(plant);
    }

    bool found_root_axillary = false;
    plant.for_each_node([&](const Node& n) {
        if (n.meristem && n.meristem->type == MeristemType::ROOT_AXILLARY) {
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
    g.internode_spacing = 100; // high so axillary doesn't interfere
    Plant plant(g, glm::vec3(0.0f));

    // After 2 ticks the shoot moves 1.0 from its start, exceeding max_internode_length of 0.6
    tick_meristems(plant);
    tick_meristems(plant);

    // The original shoot node should no longer have a meristem (became interior)
    // A new node should have the apical meristem
    int apical_count = 0;
    plant.for_each_node([&](const Node& n) {
        if (n.meristem && n.meristem->type == MeristemType::APICAL && n.meristem->active) {
            apical_count++;
        }
    });
    REQUIRE(apical_count == 1); // Still exactly one apical

    // Node count should have increased (new chain node created)
    REQUIRE(plant.node_count() > 3);
}

TEST_CASE("Axillary meristem activates when auxin low and cytokinin high", "[meristem]") {
    Genome g = default_genome();
    g.auxin_threshold = 1.0f;
    g.cytokinin_threshold = 0.0f; // cytokinin just needs to be > 0
    g.internode_spacing = 1;
    Plant plant(g, glm::vec3(0.0f));

    // Tick once to spawn an axillary node
    tick_meristems(plant);

    // Manually set hormones on the axillary node to trigger activation
    Node* axillary_node = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.meristem && n.meristem->type == MeristemType::AXILLARY) {
            axillary_node = &n;
        }
    });
    REQUIRE(axillary_node != nullptr);
    axillary_node->auxin = 0.1f;      // below threshold of 1.0
    axillary_node->cytokinin = 0.5f;   // above threshold of 0.0

    tick_meristems(plant);

    // Should have converted to APICAL and be active
    REQUIRE(axillary_node->meristem->type == MeristemType::APICAL);
    REQUIRE(axillary_node->meristem->active == true);
}

TEST_CASE("Axillary meristem stays dormant when auxin is high", "[meristem]") {
    Genome g = default_genome();
    g.auxin_threshold = 1.0f;
    g.internode_spacing = 1;
    Plant plant(g, glm::vec3(0.0f));

    tick_meristems(plant);

    Node* axillary_node = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.meristem && n.meristem->type == MeristemType::AXILLARY) {
            axillary_node = &n;
        }
    });
    REQUIRE(axillary_node != nullptr);
    axillary_node->auxin = 5.0f; // way above threshold

    tick_meristems(plant);

    REQUIRE(axillary_node->meristem->type == MeristemType::AXILLARY);
    REQUIRE(axillary_node->meristem->active == false);
}

TEST_CASE("Node age increments each tick", "[meristem]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    tick_meristems(plant);
    tick_meristems(plant);

    const Node* seed = plant.seed();
    REQUIRE(seed->age == 2);
}
