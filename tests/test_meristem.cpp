// tests/test_meristem.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <glm/geometric.hpp>
#include "engine/plant.h"
#include "engine/node/leaf_node.h"
#include "engine/node/meristems/shoot_axillary.h"
#include "engine/node/meristems/root_axillary.h"
#include "engine/world_params.h"

using namespace botany;
using Catch::Matchers::WithinAbs;

TEST_CASE("Shoot apical meristem extends node position upward", "[meristem]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::SHOOT_APICAL) { shoot = c; break; }
    }
    REQUIRE(shoot != nullptr);
    float y_before = shoot->position.y;

    plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; });
    plant.tick(default_world_params());

    REQUIRE(shoot->position.y > y_before);
    REQUIRE_THAT(shoot->position.y - y_before, WithinAbs(g.growth_rate, 0.001f));
}

TEST_CASE("Secondary growth thickens interior nodes, not tips", "[meristem]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::SHOOT_APICAL) { shoot = c; break; }
    }
    REQUIRE(shoot != nullptr);
    float seed_r_before = seed->radius;
    float shoot_r_before = shoot->radius;

    plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; });
    plant.tick(default_world_params());

    // Seed (interior node, no active meristem) should thicken
    REQUIRE(seed->radius > seed_r_before);
    // Shoot tip (active apical meristem) should NOT thicken
    REQUIRE(shoot->radius == shoot_r_before);
}

TEST_CASE("Chain growth spawns axillary node and LEAF child on interior node", "[meristem]") {
    Genome g = default_genome();
    g.growth_rate = 0.5f;
    g.shoot_plastochron = 1; // spawn internode every tick for fast testing
    Plant plant(g, glm::vec3(0.0f));

    uint32_t initial_count = plant.node_count();

    // Tick until chain growth fires (plastochron-based: 1 tick)
    // Keep auxin high so axillary buds stay dormant
    for (int i = 0; i < 2; i++) {
        plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; n.chemical(ChemicalID::Auxin) = 1.0f; n.chemical(ChemicalID::Cytokinin) = 1.0f; });
        plant.tick(default_world_params());
    }

    // Should have spawned an axillary node with a dormant meristem and a LEAF node
    REQUIRE(plant.node_count() > initial_count);

    bool found_axillary = false;
    bool found_leaf = false;
    plant.for_each_node([&](const Node& n) {
        if (n.type == NodeType::SHOOT_AXILLARY) {
            found_axillary = true;
            REQUIRE(n.as_shoot_axillary()->active == false);
        }
        if (n.type == NodeType::LEAF) {
            found_leaf = true;
        }
    });
    REQUIRE(found_axillary);
    REQUIRE(found_leaf);
}

TEST_CASE("Interior STEM nodes have at most 3 children", "[meristem]") {
    Genome g = default_genome();
    g.growth_rate = 0.2f;
    g.shoot_plastochron = 1; // spawn internode every tick for fast testing
    Plant plant(g, glm::vec3(0.0f));

    // Run enough ticks to create several interior nodes
    for (int i = 0; i < 20; i++) {
        plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; n.chemical(ChemicalID::Cytokinin) = 1.0f; });
        plant.tick(default_world_params());
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
        if (c->type == NodeType::ROOT_APICAL) { root = c; break; }
    }
    REQUIRE(root != nullptr);
    float y_before = root->position.y;

    plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; });
    plant.tick(default_world_params());

    REQUIRE(root->position.y < y_before);
}

TEST_CASE("Root chain growth spawns root axillary on interior node", "[meristem]") {
    Genome g = default_genome();
    g.root_growth_rate = 0.5f;
    g.root_plastochron = 1; // spawn root node every tick for fast testing
    Plant plant(g, glm::vec3(0.0f));

    for (int i = 0; i < 3; i++) {
        plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; });
        plant.tick(default_world_params());
    }

    bool found_root_axillary = false;
    plant.for_each_node([&](const Node& n) {
        if (n.type == NodeType::ROOT_AXILLARY) {
            found_root_axillary = true;
            REQUIRE(n.as_root_axillary()->active == false);
        }
    });
    REQUIRE(found_root_axillary);
}

TEST_CASE("Chain growth: apical meristem transfers to new node after plastochron", "[meristem]") {
    Genome g = default_genome();
    g.growth_rate = 0.5f;
    g.shoot_plastochron = 1; // spawn internode every tick for fast testing
    Plant plant(g, glm::vec3(0.0f));

    // After 1 tick the plastochron fires and spawns an internode
    // Keep auxin high so axillary buds stay dormant (don't spawn extra apicals)
    plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; n.chemical(ChemicalID::Auxin) = 1.0f; n.chemical(ChemicalID::Cytokinin) = 1.0f; });
    plant.tick(default_world_params());
    plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; n.chemical(ChemicalID::Auxin) = 1.0f; n.chemical(ChemicalID::Cytokinin) = 1.0f; });
    plant.tick(default_world_params());

    // Should still have exactly one active shoot apical meristem
    int apical_count = 0;
    plant.for_each_node([&](const Node& n) {
        if (n.type == NodeType::SHOOT_APICAL) {
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
    g.auxin_production_rate = 0.0f; // disable auxin production so manual values stick
    g.cytokinin_threshold = 0.0f;   // low threshold so activation isn't blocked by cytokinin
    g.growth_rate = 0.5f;
    g.shoot_plastochron = 1; // spawn internode every tick for fast testing
    Plant plant(g, glm::vec3(0.0f));

    // Tick twice — plastochron fires on the second tick, creating an interior node with axillary
    for (int i = 0; i < 2; i++) {
        plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; n.chemical(ChemicalID::Cytokinin) = 1.0f; });
        plant.tick(default_world_params());
    }

    // Count shoot apicals before activation
    int apical_before = 0;
    int axillary_before = 0;
    plant.for_each_node([&](const Node& n) {
        if (n.type == NodeType::SHOOT_APICAL) apical_before++;
        if (n.type == NodeType::SHOOT_AXILLARY) axillary_before++;
    });
    REQUIRE(axillary_before > 0);

    // Set parent auxin low on all axillary nodes to trigger activation
    plant.for_each_node_mut([&](Node& n) {
        n.chemical(ChemicalID::Auxin) = 0.0f;  // clear all auxin (no production to interfere)
        n.chemical(ChemicalID::Sugar) = 100.0f;
    });

    plant.tick(default_world_params());

    // Activation replaces SHOOT_AXILLARY with a new SHOOT_APICAL, so apical count should increase
    int apical_after = 0;
    plant.for_each_node([&](const Node& n) {
        if (n.type == NodeType::SHOOT_APICAL) apical_after++;
    });
    REQUIRE(apical_after > apical_before);
}

TEST_CASE("Axillary meristem stays dormant when auxin is high", "[meristem]") {
    Genome g = default_genome();
    g.auxin_threshold = 1.0f;
    g.growth_rate = 0.5f;
    g.shoot_plastochron = 1; // spawn internode every tick for fast testing
    Plant plant(g, glm::vec3(0.0f));

    // Tick twice — plastochron fires on the second tick, creating axillary buds
    for (int i = 0; i < 2; i++) {
        plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; n.chemical(ChemicalID::Cytokinin) = 1.0f; });
        plant.tick(default_world_params());
    }

    Node* axillary_node = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::SHOOT_AXILLARY) {
            axillary_node = &n;
        }
    });
    REQUIRE(axillary_node != nullptr);
    // Shoot axillaries sense auxin on their parent stem node
    axillary_node->chemical(ChemicalID::Auxin) = 5.0f; // way above threshold
    if (axillary_node->parent) {
        axillary_node->parent->chemical(ChemicalID::Auxin) = 5.0f;
    }

    plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; });
    plant.tick(default_world_params());

    REQUIRE(axillary_node->type == NodeType::SHOOT_AXILLARY);
    REQUIRE(axillary_node->as_shoot_axillary()->active == false);
}

TEST_CASE("Node age increments each tick", "[meristem]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; });
    plant.tick(default_world_params());
    plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; });
    plant.tick(default_world_params());

    const Node* seed = plant.seed();
    REQUIRE(seed->age == 2);
}

TEST_CASE("Shoot apical meristem does not grow without sugar", "[meristem][sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Find shoot tip and record position
    Node* shoot_tip = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::SHOOT_APICAL) {
            shoot_tip = &n;
        }
    });
    REQUIRE(shoot_tip != nullptr);

    // Zero sugar everywhere — should not grow
    plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 0.0f; });
    glm::vec3 pos_before = shoot_tip->position;
    plant.tick(default_world_params());
    REQUIRE(shoot_tip->position.y == pos_before.y);
}

TEST_CASE("Shoot apical meristem grows and deducts sugar", "[meristem][sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* shoot_tip = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::SHOOT_APICAL) {
            shoot_tip = &n;
        }
    });
    REQUIRE(shoot_tip != nullptr);

    shoot_tip->chemical(ChemicalID::Sugar) = 100.0f;
    float sugar_before = shoot_tip->chemical(ChemicalID::Sugar);
    plant.tick(default_world_params());
    REQUIRE(shoot_tip->chemical(ChemicalID::Sugar) < sugar_before);
}

TEST_CASE("Shoot axillary does not activate without sugar", "[meristem][sugar]") {
    Genome g = default_genome();
    g.growth_rate = 0.5f;
    g.shoot_plastochron = 1; // spawn internode every tick for fast testing
    Plant plant(g, glm::vec3(0.0f));

    // Run until axillary buds exist
    // Give all nodes plenty of sugar and cytokinin for growth
    for (int i = 0; i < 5; i++) {
        plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 100.0f; n.chemical(ChemicalID::Cytokinin) = 1.0f; });
        plant.tick(default_world_params());
    }

    // Find a dormant axillary and zero its sugar
    Node* axillary_node = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::SHOOT_AXILLARY) {
            axillary_node = &n;
        }
    });

    if (axillary_node) {
        // Zero sugar and cytokinin everywhere so:
        // - axillary can't activate (no sugar)
        // - meristems can't grow/spawn (no sugar + no cytokinin), preventing vector reallocation
        plant.for_each_node_mut([](Node& n) {
            n.chemical(ChemicalID::Sugar) = 0.0f;
            n.chemical(ChemicalID::Cytokinin) = 0.0f;
            n.chemical(ChemicalID::Auxin) = 0.0f;
        });
        plant.tick(default_world_params());
        // Should still be axillary (not enough sugar to activate)
        REQUIRE(axillary_node->type == NodeType::SHOOT_AXILLARY);
        REQUIRE(axillary_node->as_shoot_axillary()->active == false);
    }
}

TEST_CASE("Thickening does not occur without sugar", "[meristem][sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    float radius_before = seed->radius;
    seed->chemical(ChemicalID::Sugar) = 0.0f;

    plant.tick(default_world_params());
    REQUIRE(seed->radius == radius_before);
}

TEST_CASE("Thickening deducts sugar", "[meristem][sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    seed->chemical(ChemicalID::Sugar) = 100.0f;
    float sugar_before = seed->chemical(ChemicalID::Sugar);

    plant.tick(default_world_params());
    REQUIRE(seed->chemical(ChemicalID::Sugar) < sugar_before);
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
        if (n.type == NodeType::SHOOT_APICAL)
            tip1 = &n;
    });
    plant2.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::SHOOT_APICAL)
            tip2 = &n;
    });
    REQUIRE(tip1 != nullptr);
    REQUIRE(tip2 != nullptr);

    // Zero all nodes so transport doesn't interfere with sugar setup
    plant1.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 0.0f; });
    plant2.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 0.0f; });

    // Low sugar = slow growth, high sugar = full growth
    WorldParams w = default_world_params();
    float max_cost = g.growth_rate * w.sugar_cost_shoot_growth;
    tip1->chemical(ChemicalID::Sugar) = max_cost * 0.5f; // half growth
    tip2->chemical(ChemicalID::Sugar) = max_cost * 2.0f; // full growth (capped at 1.0)
    // Saturate cytokinin so it doesn't gate this test
    tip1->chemical(ChemicalID::Cytokinin) = g.cytokinin_growth_threshold * 2.0f;
    tip2->chemical(ChemicalID::Cytokinin) = g.cytokinin_growth_threshold * 2.0f;

    glm::vec3 pos1_before = tip1->position;
    glm::vec3 pos2_before = tip2->position;

    plant1.tick(default_world_params());
    plant1.recompute_world_positions();
    plant2.tick(default_world_params());
    plant2.recompute_world_positions();

    float dist1 = glm::length(tip1->position - pos1_before);
    float dist2 = glm::length(tip2->position - pos2_before);

    // Half sugar should produce less growth than full sugar
    REQUIRE(dist1 > 0.0f);
    REQUIRE(dist2 > dist1);
}

TEST_CASE("Zero sugar produces zero growth", "[meristem][sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* shoot_tip = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::SHOOT_APICAL)
            shoot_tip = &n;
    });
    REQUIRE(shoot_tip != nullptr);

    // Zero sugar everywhere — no growth possible
    plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 0.0f; });
    glm::vec3 pos_before = shoot_tip->position;
    plant.tick(default_world_params());
    REQUIRE(shoot_tip->position.y == pos_before.y);
}

TEST_CASE("GA boosts intercalary elongation rate", "[meristem][gibberellin]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
    WorldParams wp = default_world_params();

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    stem->age = 1;
    stem->chemical(ChemicalID::Sugar) = 5.0f;
    plant.seed_mut()->add_child(stem);

    // Run without GA
    stem->chemical(ChemicalID::Gibberellin) = 0.0f;
    float offset_before = glm::length(stem->offset);
    plant.tick(wp);
    float growth_no_ga = glm::length(stem->offset) - offset_before;

    // Reset
    stem->offset = glm::vec3(0.0f, 0.5f, 0.0f);
    stem->age = 1;
    stem->chemical(ChemicalID::Sugar) = 5.0f;

    // Run with GA
    stem->chemical(ChemicalID::Gibberellin) = 1.0f;
    offset_before = glm::length(stem->offset);
    plant.tick(wp);
    float growth_with_ga = glm::length(stem->offset) - offset_before;

    REQUIRE(growth_with_ga > growth_no_ga);
}

TEST_CASE("Ethylene inhibits elongation", "[meristem][ethylene]") {
    WorldParams wp = default_world_params();

    // Plant 1: no ethylene on stem
    Genome g = default_genome();
    Plant plant1(g, glm::vec3(0.0f));

    Node* stem1 = plant1.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    stem1->age = 1;
    stem1->chemical(ChemicalID::Sugar) = 5.0f;
    plant1.seed_mut()->chemical(ChemicalID::Sugar) = 5.0f;
    plant1.seed_mut()->add_child(stem1);

    // Plant 2: inject ethylene directly onto stem
    Plant plant2(g, glm::vec3(0.0f));

    Node* stem2 = plant2.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.05f);
    stem2->age = 1;
    stem2->chemical(ChemicalID::Sugar) = 5.0f;
    stem2->chemical(ChemicalID::Ethylene) = 1.0f;  // direct injection
    plant2.seed_mut()->chemical(ChemicalID::Sugar) = 5.0f;
    plant2.seed_mut()->add_child(stem2);

    float offset1_before = glm::length(stem1->offset);
    float offset2_before = glm::length(stem2->offset);

    plant1.tick(wp);
    plant2.tick(wp);

    float growth1 = glm::length(stem1->offset) - offset1_before;
    float growth2 = glm::length(stem2->offset) - offset2_before;

    REQUIRE(growth1 > 0.0f);
    REQUIRE(growth2 < growth1);
}

TEST_CASE("Thickening scales with sugar level", "[meristem][sugar]") {
    Genome g = default_genome();
    Plant plant1(g, glm::vec3(0.0f));
    Plant plant2(g, glm::vec3(0.0f));

    Node* seed1 = plant1.seed_mut();
    Node* seed2 = plant2.seed_mut();

    WorldParams w = default_world_params();
    // With reserve_fraction, a node with more sugar has more available for growth
    // (available = sugar * (1 - reserve_fraction))
    // Thickening cost is tiny (~0.00002g). Need very low sugar so
    // reserve fraction actually limits growth.
    // Zero ALL sugar in plant1, then give seed barely any — prevents diffusion from children
    plant1.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 0.0f; });
    seed1->chemical(ChemicalID::Sugar) = 0.00005f; // barely any

    // Plant2 gets plenty everywhere so diffusion doesn't drain the seed
    plant2.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 1.0f; });

    float r1_before = seed1->radius;
    float r2_before = seed2->radius;

    plant1.tick(default_world_params());
    plant2.tick(default_world_params());

    float thicken1 = seed1->radius - r1_before;
    float thicken2 = seed2->radius - r2_before;

    REQUIRE(thicken1 > 0.0f);
    REQUIRE(thicken2 > thicken1);
}
