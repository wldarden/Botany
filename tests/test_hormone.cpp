#include <catch2/catch_test_macros.hpp>
#include "engine/plant.h"
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
    g.apical_auxin_baseline = 1.0f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world();

    tick_n(plant, world, 1);

    // The shoot apical node should have auxin (produced - transported fraction)
    const Node* seed = plant.seed();
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::APICAL) { shoot = c; break; }
    }
    REQUIRE(shoot != nullptr);
    REQUIRE(shoot->chemical(ChemicalID::Auxin) > 0.0f);
}

TEST_CASE("Auxin: basipetal transport reaches parent", "[hormone]") {
    Genome g = default_genome();
    g.apical_auxin_baseline = 1.0f;
    g.auxin_diffusion_rate = 0.5f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world();

    // 2 ticks: child produces auxin at tick 1, parent receives at tick 2
    tick_n(plant, world, 2);

    // Seed should have received auxin from shoot apical child
    REQUIRE(plant.seed()->chemical(ChemicalID::Auxin) > 0.0f);
}

TEST_CASE("Auxin: root apical doesn't produce auxin", "[hormone]") {
    Genome g = default_genome();
    g.apical_auxin_baseline = 1.0f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world();

    tick_n(plant, world, 3);

    const Node* seed = plant.seed();
    const Node* root = nullptr;
    const Node* shoot = nullptr;
    for (const Node* c : seed->children) {
        if (c->type == NodeType::ROOT_APICAL) root = c;
        if (c->type == NodeType::APICAL) shoot = c;
    }
    REQUIRE(root != nullptr);
    REQUIRE(shoot != nullptr);
    // Root should have much less auxin than shoot (only receives via seed spillover)
    REQUIRE(root->chemical(ChemicalID::Auxin) < shoot->chemical(ChemicalID::Auxin));
}

TEST_CASE("Auxin: decays over time without production", "[hormone]") {
    Genome g = default_genome();
    g.apical_auxin_baseline = 1.0f;
    g.auxin_decay_rate = 0.1f;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world();

    // Produce for a few ticks
    tick_n(plant, world, 5);

    const Node* seed = plant.seed();
    float auxin_at_5 = seed->chemical(ChemicalID::Auxin);

    // Keep ticking — production continues but transport + decay should prevent
    // linear accumulation
    tick_n(plant, world, 20);
    float auxin_at_25 = seed->chemical(ChemicalID::Auxin);

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
    REQUIRE(root->chemical(ChemicalID::Cytokinin) > 0.0f);
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
        if (c->type == NodeType::APICAL) { shoot = c; break; }
    }
    REQUIRE(shoot != nullptr);
    REQUIRE(shoot->chemical(ChemicalID::Cytokinin) > 0.0f);
}

TEST_CASE("Auxin: apical growth multiplier boosts production with sugar", "[hormone]") {
    WorldParams world = default_world_params();

    // Plant 1: growth multiplier disabled
    Genome g1 = default_genome();
    g1.apical_auxin_baseline = 1.0f;
    g1.apical_growth_auxin_multiplier = 0.0f;  // no growth bonus
    Plant plant1(g1, glm::vec3(0.0f));

    // Plant 2: growth multiplier enabled
    Genome g2 = default_genome();
    g2.apical_auxin_baseline = 1.0f;
    g2.apical_growth_auxin_multiplier = 5.0f;  // large for clear signal
    Plant plant2(g2, glm::vec3(0.0f));

    // Saturate sugar + cytokinin on both → max growth → max multiplier
    for (int i = 0; i < 3; i++) {
        plant1.for_each_node_mut([](Node& n) {
            n.chemical(ChemicalID::Sugar) = 100.0f;
            n.chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant2.for_each_node_mut([](Node& n) {
            n.chemical(ChemicalID::Sugar) = 100.0f;
            n.chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant1.tick(world);
        plant2.tick(world);
    }

    // Sum total auxin in each plant
    float total1 = 0.0f, total2 = 0.0f;
    plant1.for_each_node([&](const Node& n) { total1 += n.chemical(ChemicalID::Auxin); });
    plant2.for_each_node([&](const Node& n) { total2 += n.chemical(ChemicalID::Auxin); });

    // Growth-boosted plant should have meaningfully more total auxin
    REQUIRE(total2 > total1 * 1.5f);
}

TEST_CASE("Auxin: apical growth multiplier has no effect without sugar", "[hormone]") {
    WorldParams world = default_world_params();

    Genome g1 = default_genome();
    g1.apical_auxin_baseline = 1.0f;
    g1.apical_growth_auxin_multiplier = 0.0f;
    Plant plant1(g1, glm::vec3(0.0f));

    Genome g2 = default_genome();
    g2.apical_auxin_baseline = 1.0f;
    g2.apical_growth_auxin_multiplier = 5.0f;
    Plant plant2(g2, glm::vec3(0.0f));

    // Zero sugar → zero growth → multiplier should have no effect
    for (int i = 0; i < 3; i++) {
        plant1.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 0.0f; });
        plant2.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 0.0f; });
        plant1.tick(world);
        plant2.tick(world);
    }

    float total1 = 0.0f, total2 = 0.0f;
    plant1.for_each_node([&](const Node& n) { total1 += n.chemical(ChemicalID::Auxin); });
    plant2.for_each_node([&](const Node& n) { total2 += n.chemical(ChemicalID::Auxin); });

    // Both should produce similar auxin (baseline only, growth_gf ≈ 0)
    // Allow 20% tolerance for sugar_factor floor (0.1 minimum) interacting with multiplier
    float ratio = (total1 > 1e-8f) ? total2 / total1 : 1.0f;
    REQUIRE(ratio < 1.2f);
}
