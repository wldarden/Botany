#include <catch2/catch_test_macros.hpp>
#include "engine/plant.h"
#include "engine/sugar.h"
#include "engine/world_params.h"
#include "engine/node/tissues/leaf.h"

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
    REQUIRE(shoot->local().chemical(ChemicalID::Auxin) > 0.0f);
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
    REQUIRE(plant.seed()->local().chemical(ChemicalID::Auxin) > 0.0f);
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
    REQUIRE(root->local().chemical(ChemicalID::Auxin) < shoot->local().chemical(ChemicalID::Auxin));
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
    float auxin_at_5 = seed->local().chemical(ChemicalID::Auxin);

    // Keep ticking — production continues but transport + decay should prevent
    // linear accumulation
    tick_n(plant, world, 20);
    float auxin_at_25 = seed->local().chemical(ChemicalID::Auxin);

    // Should not be 8x the earlier amount (decay limits accumulation).
    // With slower diffusion (0.05) auxin disperses more slowly so the local
    // concentration at tick 5 is still low; by tick 25 it has climbed but
    // the system is still far from the 20x that linear production would give.
    REQUIRE(auxin_at_25 < 8.0f * auxin_at_5);
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
    // Production is tracked per-tick in tick_cytokinin_produced — diagnostic counter
    // that is NOT affected by transport, diffusion, or decay.  The RA's local pool
    // is drained each tick by vascular inject + local diffusion to the rest of the
    // plant (seed, shoot), so checking the local pool alone understates production.
    REQUIRE(root->tick_cytokinin_produced > 0.0f);
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
    REQUIRE(shoot->local().chemical(ChemicalID::Cytokinin) > 0.0f);
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
            n.local().chemical(ChemicalID::Sugar) = 100.0f;
            n.local().chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant2.for_each_node_mut([](Node& n) {
            n.local().chemical(ChemicalID::Sugar) = 100.0f;
            n.local().chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant1.tick(world);
        plant2.tick(world);
    }

    // Sum total auxin in each plant
    float total1 = 0.0f, total2 = 0.0f;
    plant1.for_each_node([&](const Node& n) { total1 += n.local().chemical(ChemicalID::Auxin); });
    plant2.for_each_node([&](const Node& n) { total2 += n.local().chemical(ChemicalID::Auxin); });

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
        plant1.for_each_node_mut([](Node& n) { n.local().chemical(ChemicalID::Sugar) = 0.0f; });
        plant2.for_each_node_mut([](Node& n) { n.local().chemical(ChemicalID::Sugar) = 0.0f; });
        plant1.tick(world);
        plant2.tick(world);
    }

    float total1 = 0.0f, total2 = 0.0f;
    plant1.for_each_node([&](const Node& n) { total1 += n.local().chemical(ChemicalID::Auxin); });
    plant2.for_each_node([&](const Node& n) { total2 += n.local().chemical(ChemicalID::Auxin); });

    // Both should produce similar auxin (baseline only, growth_gf ≈ 0)
    // Allow 20% tolerance for sugar_factor floor (0.1 minimum) interacting with multiplier
    float ratio = (total1 > 1e-8f) ? total2 / total1 : 1.0f;
    REQUIRE(ratio < 1.2f);
}

TEST_CASE("Auxin: growing leaf produces auxin", "[hormone]") {
    Genome g = default_genome();
    g.apical_auxin_baseline = 0.0f;          // disable meristem auxin
    g.leaf_auxin_baseline = 1.0f;             // high for easy detection
    g.leaf_growth_auxin_multiplier = 0.5f;
    g.growth_rate = 0.5f;
    g.shoot_plastochron = 1;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world_params();

    // Grow plant until leaves exist
    for (int i = 0; i < 5; i++) {
        plant.for_each_node_mut([](Node& n) {
            n.local().chemical(ChemicalID::Sugar) = 100.0f;
            n.local().chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant.tick(world);
    }

    // Find a growing leaf
    LeafNode* growing_leaf = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (auto* leaf = n.as_leaf()) {
            if (leaf->leaf_size < g.max_leaf_size) growing_leaf = leaf;
        }
    });
    REQUIRE(growing_leaf != nullptr);

    // Zero all auxin, give sugar + water for leaf growth
    plant.for_each_node_mut([&](Node& n) {
        n.local().chemical(ChemicalID::Auxin) = 0.0f;
        n.local().chemical(ChemicalID::Sugar) = 100.0f;
        n.local().chemical(ChemicalID::Water) = water_cap(n, g);
    });

    plant.tick(world);

    // Growing leaf should have produced auxin (some may have transported, but should retain some)
    REQUIRE(growing_leaf->local().chemical(ChemicalID::Auxin) > 0.0f);
}

TEST_CASE("Auxin: full-size leaf produces zero auxin", "[hormone]") {
    Genome g = default_genome();
    g.apical_auxin_baseline = 0.0f;          // disable meristem auxin
    g.root_tip_auxin_production_rate = 0.0f;  // disable root tip auxin
    g.leaf_auxin_baseline = 1.0f;
    g.leaf_growth_auxin_multiplier = 0.5f;
    g.growth_rate = 0.5f;
    g.shoot_plastochron = 1;
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world_params();

    // Grow plant until leaves exist
    for (int i = 0; i < 5; i++) {
        plant.for_each_node_mut([](Node& n) {
            n.local().chemical(ChemicalID::Sugar) = 100.0f;
            n.local().chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        plant.tick(world);
    }

    // Force ALL leaves to max size
    plant.for_each_node_mut([&](Node& n) {
        if (auto* leaf = n.as_leaf()) {
            leaf->leaf_size = g.max_leaf_size;
        }
    });

    // Zero all auxin everywhere
    plant.for_each_node_mut([](Node& n) {
        n.local().chemical(ChemicalID::Auxin) = 0.0f;
        n.local().chemical(ChemicalID::Sugar) = 100.0f;
    });

    plant.tick(world);

    // No auxin source anywhere → all auxin should be zero
    float total_auxin = 0.0f;
    plant.for_each_node([&](const Node& n) { total_auxin += n.local().chemical(ChemicalID::Auxin); });
    REQUIRE(total_auxin < 1e-6f);
}

TEST_CASE("Auxin: leaf auxin scales with growth amount", "[hormone]") {
    WorldParams world = default_world_params();

    // Plant 1: low leaf multiplier
    Genome g1 = default_genome();
    g1.apical_auxin_baseline = 0.0f;
    g1.leaf_auxin_baseline = 1.0f;
    g1.leaf_growth_auxin_multiplier = 0.1f;
    g1.growth_rate = 0.5f;
    g1.shoot_plastochron = 1;
    Plant plant1(g1, glm::vec3(0.0f));

    // Plant 2: high leaf multiplier
    Genome g2 = default_genome();
    g2.apical_auxin_baseline = 0.0f;
    g2.leaf_auxin_baseline = 1.0f;
    g2.leaf_growth_auxin_multiplier = 0.9f;
    g2.growth_rate = 0.5f;
    g2.shoot_plastochron = 1;
    Plant plant2(g2, glm::vec3(0.0f));

    // Grow both until leaves exist, then measure auxin.
    // Set water high so leaves are fully hydrated — this test is about auxin, not water balance.
    for (int i = 0; i < 5; i++) {
        plant1.for_each_node_mut([](Node& n) {
            n.local().chemical(ChemicalID::Sugar) = 100.0f;
            n.local().chemical(ChemicalID::Cytokinin) = 1.0f;
            n.local().chemical(ChemicalID::Water) = 100.0f;
        });
        plant2.for_each_node_mut([](Node& n) {
            n.local().chemical(ChemicalID::Sugar) = 100.0f;
            n.local().chemical(ChemicalID::Cytokinin) = 1.0f;
            n.local().chemical(ChemicalID::Water) = 100.0f;
        });
        plant1.tick(world);
        plant2.tick(world);
    }

    // Zero auxin, give sugar + water, tick once more
    plant1.for_each_node_mut([](Node& n) {
        n.local().chemical(ChemicalID::Auxin) = 0.0f;
        n.local().chemical(ChemicalID::Sugar) = 100.0f;
        n.local().chemical(ChemicalID::Water) = 100.0f;
    });
    plant2.for_each_node_mut([](Node& n) {
        n.local().chemical(ChemicalID::Auxin) = 0.0f;
        n.local().chemical(ChemicalID::Sugar) = 100.0f;
        n.local().chemical(ChemicalID::Water) = 100.0f;
    });
    plant1.tick(world);
    plant2.tick(world);

    float total1 = 0.0f, total2 = 0.0f;
    plant1.for_each_node([&](const Node& n) { total1 += n.local().chemical(ChemicalID::Auxin); });
    plant2.for_each_node([&](const Node& n) { total2 += n.local().chemical(ChemicalID::Auxin); });

    // Higher multiplier → more leaf auxin in the system
    REQUIRE(total2 > total1 * 2.0f);
}

// --- Cross-seed transport tests ---
// These tests verify that chemicals produced on one side of the seed node
// actually cross it and reach the other side.  Bootstrap cytokinin/auxin
// is zeroed out so only newly-produced chemical drives the measurement.

TEST_CASE("Transport: auxin crosses seed node from shoot to root", "[hormone][transport]") {
    Genome g = default_genome();
    g.apical_auxin_baseline = 1.0f;      // strong signal for clear detection
    g.auxin_diffusion_rate  = 0.3f;      // faster than default so it spreads in fewer ticks
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world_params();

    // Zero all auxin so we are measuring only shoot-produced auxin.
    // (Bootstrap auxin in seed/shoot would otherwise mask whether transport works.)
    plant.for_each_node_mut([](Node& n) { n.local().chemical(ChemicalID::Auxin) = 0.0f; });

    for (int i = 0; i < 50; i++) {
        plant.for_each_node_mut([](Node& n) {
            n.local().chemical(ChemicalID::Sugar) = 10.0f;   // keep nodes alive & growing
        });
        plant.tick(world);
    }

    // After growth the apical has chained upward — search all nodes.
    const Node* shoot = nullptr;
    const Node* root  = nullptr;
    plant.for_each_node([&](const Node& n) {
        if (n.type == NodeType::APICAL     && shoot == nullptr) shoot = &n;
        if (n.type == NodeType::ROOT_APICAL && root  == nullptr) root  = &n;
    });
    REQUIRE(shoot != nullptr);
    REQUIRE(root  != nullptr);

    // Shoot produces auxin — must have some remaining after transport+decay
    REQUIRE(shoot->local().chemical(ChemicalID::Auxin) > 0.0f);

    // Auxin must have crossed the seed and reached the root apical.
    // Auxin accumulates at its source (shoot tip), so root level is lower than
    // shoot level — we only require root > 0, i.e. transport DID cross the seed.
    REQUIRE(root->local().chemical(ChemicalID::Auxin) > 0.0f);
}

TEST_CASE("Transport: cytokinin crosses seed node from root to shoot", "[hormone][transport]") {
    Genome g = default_genome();
    g.root_cytokinin_production_rate = 1.0f; // strong signal for clear detection
    g.cytokinin_diffusion_rate       = 0.3f; // faster than default
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world = default_world_params();

    // Zero all cytokinin so we measure only root-produced cytokinin.
    plant.for_each_node_mut([](Node& n) { n.local().chemical(ChemicalID::Cytokinin) = 0.0f; });

    for (int i = 0; i < 50; i++) {
        plant.for_each_node_mut([](Node& n) {
            n.local().chemical(ChemicalID::Sugar) = 10.0f;
        });
        plant.tick(world);
    }

    // After growth the apical has chained upward — search all nodes.
    const Node* shoot = nullptr;
    const Node* root  = nullptr;
    plant.for_each_node([&](const Node& n) {
        if (n.type == NodeType::APICAL     && shoot == nullptr) shoot = &n;
        if (n.type == NodeType::ROOT_APICAL && root  == nullptr) root  = &n;
    });
    REQUIRE(shoot != nullptr);
    REQUIRE(root  != nullptr);

    // Root produces cytokinin — verify by checking total plant cytokinin > 0.
    // Under tick-then-vascular ordering, the RA's local cytokinin is moved into the
    // xylem stream each tick (by vascular_sub_stepped), so root->local() will be near 0.
    // Instead, verify that cytokinin is present somewhere in the plant.
    float cyto_total = 0.0f;
    plant.for_each_node([&](const Node& n) {
        cyto_total += n.local().chemical(ChemicalID::Cytokinin);
        if (auto* xyl = n.xylem()) cyto_total += xyl->chemical(ChemicalID::Cytokinin);
    });
    REQUIRE(cyto_total > 0.0f);

    // Cytokinin must have crossed the seed and reached the shoot apical.
    // Cytokinin accumulates at its source (root tip), so root level stays higher
    // than shoot level — we only require shoot > 0, i.e. transport DID cross the seed.
    REQUIRE(shoot->local().chemical(ChemicalID::Cytokinin) > 0.0f);
}
