#include <catch2/catch_test_macros.hpp>
#include "engine/plant.h"
#include "engine/node/meristem_node.h"
#include "engine/world_params.h"

using namespace botany;

static WorldParams default_world() {
    WorldParams w;
    return w;
}

// Helper: tick just enough for auxin to propagate (no chain growth)
static void tick_n(Plant& plant, const WorldParams& world, int n) {
    for (int i = 0; i < n; i++) {
        plant.tick(world);
    }
}

TEST_CASE("Auxin: shoot apical produces auxin", "[hormone]") {
    Genome g = default_genome();
    g.auxin_production_rate = 1.0f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world();

    tick_n(plant, world, 1);

    // The shoot apical node should have auxin (produced - transported fraction)
    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::SHOOT_APICAL) { shoot = c; break; }
    }
    REQUIRE(shoot != nullptr);
    REQUIRE(shoot->auxin > 0.0f);
}

TEST_CASE("Auxin: basipetal transport reaches parent", "[hormone]") {
    Genome g = default_genome();
    g.auxin_production_rate = 1.0f;
    g.auxin_transport_rate = 0.5f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world();

    tick_n(plant, world, 1);

    // Seed should have received auxin from shoot apical child
    REQUIRE(plant.seed()->auxin > 0.0f);
}

TEST_CASE("Auxin: root apical doesn't produce auxin", "[hormone]") {
    Genome g = default_genome();
    g.auxin_production_rate = 1.0f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world();

    tick_n(plant, world, 3);

    const Node* seed = plant.seed();
    const Node* root = nullptr;
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::ROOT_APICAL) root = c;
        if (c->type == NodeType::SHOOT_APICAL) shoot = c;
    }
    REQUIRE(root != nullptr);
    REQUIRE(shoot != nullptr);
    // Root should have much less auxin than shoot (only receives via seed spillover)
    REQUIRE(root->auxin < shoot->auxin);
}

TEST_CASE("Auxin: decays over time without production", "[hormone]") {
    Genome g = default_genome();
    g.auxin_production_rate = 1.0f;
    g.auxin_decay_rate = 0.1f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world();

    // Produce for a few ticks
    tick_n(plant, world, 5);

    const Node* seed = plant.seed();
    float auxin_at_5 = seed->auxin;

    // Keep ticking — production continues but transport + decay should prevent
    // linear accumulation
    tick_n(plant, world, 20);
    float auxin_at_25 = seed->auxin;

    // Should not be 5x the earlier amount (decay limits accumulation)
    REQUIRE(auxin_at_25 < 5.0f * auxin_at_5);
}

TEST_CASE("Cytokinin: root apical produces cytokinin", "[hormone]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world();

    tick_n(plant, world, 1);

    const Node* seed = plant.seed();
    const Node* root = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::ROOT_APICAL) { root = c; break; }
    }
    REQUIRE(root != nullptr);
    REQUIRE(root->cytokinin > 0.0f);
}

TEST_CASE("Cytokinin: cytokinin flows from parent to children", "[hormone]") {
    Genome g = default_genome();
    g.cytokinin_production_rate = 2.0f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world();

    tick_n(plant, world, 5);

    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::SHOOT_APICAL) { shoot = c; break; }
    }
    REQUIRE(shoot != nullptr);
    REQUIRE(shoot->cytokinin > 0.0f);
}
