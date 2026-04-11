// tests/test_meristem.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <glm/geometric.hpp>
#include "engine/plant.h"
#include "engine/meristems/meristem.h"
#include "engine/world_params.h"

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
    tick_meristems(plant, default_world_params());
    plant.recompute_world_positions();

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
    tick_meristems(plant, default_world_params());

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
    // Keep auxin high so axillary buds stay dormant
    for (int i = 0; i < 2; i++) {
        plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; n.auxin = 1.0f; });
        tick_meristems(plant, default_world_params());
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
            REQUIRE(n.as_leaf()->leaf_size == g.leaf_bud_size);
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
        tick_meristems(plant, default_world_params());
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
    tick_meristems(plant, default_world_params());
    plant.recompute_world_positions();

    REQUIRE(root->position.y < y_before);
}

TEST_CASE("Root chain growth spawns root axillary on interior node", "[meristem]") {
    Genome g = default_genome();
    g.root_growth_rate = 0.5f;
    g.root_max_internode_length = 0.6f;
    Plant plant(g, glm::vec3(0.0f));

    for (int i = 0; i < 2; i++) {
        plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
        tick_meristems(plant, default_world_params());
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
    // Keep auxin high so axillary buds stay dormant (don't spawn extra apicals)
    plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; n.auxin = 1.0f; });
    tick_meristems(plant, default_world_params());
    plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; n.auxin = 1.0f; });
    tick_meristems(plant, default_world_params());

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
    tick_meristems(plant, default_world_params());

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
    tick_meristems(plant, default_world_params());

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
    tick_meristems(plant, default_world_params());

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
    tick_meristems(plant, default_world_params());

    REQUIRE(axillary_node->meristem->type() == MeristemType::AXILLARY);
    REQUIRE(axillary_node->meristem->active == false);
}

TEST_CASE("Node age increments each tick", "[meristem]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
    tick_meristems(plant, default_world_params());
    plant.for_each_node_mut([](Node& n) { n.sugar = 100.0f; });
    tick_meristems(plant, default_world_params());

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
    tick_meristems(plant, default_world_params());
    plant.recompute_world_positions();
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
    tick_meristems(plant, default_world_params());
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
        tick_meristems(plant, default_world_params());
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
        tick_meristems(plant, default_world_params());
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

    tick_meristems(plant, default_world_params());
    REQUIRE(seed->radius == radius_before);
}

TEST_CASE("Thickening deducts sugar", "[meristem][sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    seed->sugar = 100.0f;
    float sugar_before = seed->sugar;

    tick_meristems(plant, default_world_params());
    REQUIRE(seed->sugar < sugar_before);
    REQUIRE(seed->radius > g.initial_radius);
}

TEST_CASE("Shoot growth scales with sugar level", "[meristem][sugar]") {
    Genome g = default_genome();
    Plant plant1(g, glm::vec3(0.0f));
    Plant plant2(g, glm::vec3(0.0f));

    // Find shoot tips in both plants
    Node* tip1 = nullptr;
    Node* tip2 = nullptr;
    plant1.for_each_node_mut([&](Node& n) {
        if (n.meristem && n.meristem->type() == MeristemType::APICAL && n.meristem->active)
            tip1 = &n;
    });
    plant2.for_each_node_mut([&](Node& n) {
        if (n.meristem && n.meristem->type() == MeristemType::APICAL && n.meristem->active)
            tip2 = &n;
    });
    REQUIRE(tip1 != nullptr);
    REQUIRE(tip2 != nullptr);

    // Low sugar = slow growth, high sugar = full growth
    WorldParams w = default_world_params();
    float max_cost = g.growth_rate * w.sugar_cost_shoot_growth;
    tip1->sugar = g.sugar_save_shoot + max_cost * 0.5f; // half growth
    tip2->sugar = g.sugar_save_shoot + max_cost * 2.0f; // full growth (capped at 1.0)

    glm::vec3 pos1_before = tip1->position;
    glm::vec3 pos2_before = tip2->position;

    tick_meristems(plant1, default_world_params());
    plant1.recompute_world_positions();
    tick_meristems(plant2, default_world_params());
    plant2.recompute_world_positions();

    float dist1 = glm::length(tip1->position - pos1_before);
    float dist2 = glm::length(tip2->position - pos2_before);

    // Half sugar should produce less growth than full sugar
    REQUIRE(dist1 > 0.0f);
    REQUIRE(dist2 > dist1);
}

TEST_CASE("Growth at save_threshold produces zero growth", "[meristem][sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* shoot_tip = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.meristem && n.meristem->type() == MeristemType::APICAL && n.meristem->active)
            shoot_tip = &n;
    });
    REQUIRE(shoot_tip != nullptr);

    // Exactly at save threshold — no growth
    shoot_tip->sugar = g.sugar_save_shoot;
    glm::vec3 pos_before = shoot_tip->position;
    tick_meristems(plant, default_world_params());
    plant.recompute_world_positions();
    REQUIRE(shoot_tip->position.y == pos_before.y);
}

TEST_CASE("GA boosts intercalary elongation rate", "[meristem][gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
    WorldParams wp = default_world_params();

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    stem->age = 1;
    stem->sugar = 5.0f;
    plant.seed_mut()->add_child(stem);

    // Run without GA
    stem->gibberellin = 0.0f;
    float offset_before = glm::length(stem->offset);
    tick_meristems(plant, wp);
    float growth_no_ga = glm::length(stem->offset) - offset_before;

    // Reset
    stem->offset = glm::vec3(0.0f, 0.5f, 0.0f);
    stem->age = 1;
    stem->sugar = 5.0f;

    // Run with GA
    stem->gibberellin = 1.0f;
    offset_before = glm::length(stem->offset);
    tick_meristems(plant, wp);
    float growth_with_ga = glm::length(stem->offset) - offset_before;

    REQUIRE(growth_with_ga > growth_no_ga);
}

TEST_CASE("Ethylene inhibits elongation", "[meristem][ethylene]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
    WorldParams wp = default_world_params();

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    stem->age = 1;
    stem->sugar = 5.0f;
    plant.seed_mut()->add_child(stem);

    // Run without ethylene
    stem->ethylene = 0.0f;
    float offset_before = glm::length(stem->offset);
    tick_meristems(plant, wp);
    float growth_no_eth = glm::length(stem->offset) - offset_before;

    // Reset
    stem->offset = glm::vec3(0.0f, 0.5f, 0.0f);
    stem->age = 1;
    stem->sugar = 5.0f;

    // Run with high ethylene
    stem->ethylene = 2.0f;
    offset_before = glm::length(stem->offset);
    tick_meristems(plant, wp);
    float growth_with_eth = glm::length(stem->offset) - offset_before;

    REQUIRE(growth_with_eth < growth_no_eth);
}

TEST_CASE("Thickening scales with sugar level", "[meristem][sugar]") {
    Genome g = default_genome();
    Plant plant1(g, glm::vec3(0.0f));
    Plant plant2(g, glm::vec3(0.0f));

    Node* seed1 = plant1.seed_mut();
    Node* seed2 = plant2.seed_mut();

    WorldParams w = default_world_params();
    float max_cost = g.thickening_rate * w.sugar_cost_thickening;
    seed1->sugar = g.sugar_save_stem + max_cost * 0.5f;
    seed2->sugar = g.sugar_save_stem + max_cost * 2.0f;

    float r1_before = seed1->radius;
    float r2_before = seed2->radius;

    // Give shoot tips sugar so they don't interfere
    plant1.for_each_node_mut([&](Node& n) {
        if (n.meristem && n.meristem->is_tip()) n.sugar = 100.0f;
    });
    plant2.for_each_node_mut([&](Node& n) {
        if (n.meristem && n.meristem->is_tip()) n.sugar = 100.0f;
    });

    tick_meristems(plant1, default_world_params());
    tick_meristems(plant2, default_world_params());

    float thicken1 = seed1->radius - r1_before;
    float thicken2 = seed2->radius - r2_before;

    REQUIRE(thicken1 > 0.0f);
    REQUIRE(thicken2 > thicken1);
}
