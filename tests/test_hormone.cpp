#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/plant.h"
#include "engine/hormone.h"

using namespace botany;
using Catch::Matchers::WithinAbs;

// Helper: build a simple 3-node linear stem (seed -> mid -> tip_with_apical)
static Plant make_linear_stem() {
    Genome g = default_genome();
    g.auxin_production_rate = 1.0f;
    g.auxin_transport_rate = 0.5f;
    g.auxin_decay_rate = 0.1f;
    Plant plant(g, glm::vec3(0.0f));
    // Plant already has: seed(0) -> shoot(1), seed(0) -> root(2)
    // shoot(1) has an APICAL meristem
    return plant;
}

TEST_CASE("transport_auxin: apical meristem produces auxin at its node", "[hormone]") {
    Plant plant = make_linear_stem();
    const Genome& g = plant.genome();

    transport_auxin(plant);

    // The shoot node (child 0 of seed, which is STEM type) should have auxin
    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::STEM) { shoot = c; break; }
    }
    REQUIRE(shoot != nullptr);
    REQUIRE(shoot->auxin > 0.0f);
}

TEST_CASE("transport_auxin: auxin flows from child to parent", "[hormone]") {
    Plant plant = make_linear_stem();

    transport_auxin(plant);

    const Node* seed = plant.seed();
    // Seed should have received auxin from the shoot child
    REQUIRE(seed->auxin > 0.0f);
}

TEST_CASE("transport_auxin: non-meristem nodes don't produce auxin", "[hormone]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Root node has ROOT_APICAL, not APICAL — should not produce auxin
    transport_auxin(plant);

    const Node* seed = plant.seed();
    const Node* root = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::ROOT) { root = c; break; }
    }
    REQUIRE(root != nullptr);
    // Root doesn't produce auxin, but may receive a small amount via
    // spillback from the seed node. Verify it's much less than what
    // the shoot tip's side of the tree gets.
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::STEM) { shoot = c; break; }
    }
    REQUIRE(root->auxin < shoot->auxin);
}

TEST_CASE("transport_cytokinin: root apical produces cytokinin", "[hormone]") {
    Plant plant = make_linear_stem();

    transport_cytokinin(plant);

    const Node* seed = plant.seed();
    const Node* root = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::ROOT) { root = c; break; }
    }
    REQUIRE(root != nullptr);
    REQUIRE(root->cytokinin > 0.0f);
}

TEST_CASE("transport_cytokinin: cytokinin flows from parent to children", "[hormone]") {
    Genome g = default_genome();
    g.cytokinin_production_rate = 2.0f;
    Plant plant(g, glm::vec3(0.0f));

    // Run multiple passes so cytokinin propagates from root through seed to shoot
    for (int i = 0; i < 5; i++) {
        transport_cytokinin(plant);
    }

    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::STEM) { shoot = c; break; }
    }
    REQUIRE(shoot != nullptr);
    REQUIRE(shoot->cytokinin > 0.0f);
}

TEST_CASE("transport_auxin: auxin decays each pass", "[hormone]") {
    Plant plant = make_linear_stem();

    transport_auxin(plant);
    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::STEM) { shoot = c; break; }
    }

    float after_one = shoot->auxin;

    transport_auxin(plant);
    transport_auxin(plant);
    transport_auxin(plant);

    float after_four = shoot->auxin;
    REQUIRE(after_four < 4.0f * after_one);  // Not linear growth
}
