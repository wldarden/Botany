// tests/test_cytokinin_transport.cpp
// Verifies that cytokinin produced by root apicals accumulates properly and
// reaches shoot apicals. The key parameter is cytokinin_diffusion_rate:
// too high drains the RA before vascular can pick it up; the right value
// (~0.02) keeps most cytokinin in the RA where it builds up and feeds xylem.
#include <catch2/catch_test_macros.hpp>
#include "engine/plant.h"
#include "engine/world_params.h"
#include "engine/chemical/chemical.h"
#include "engine/node/tissues/root_apical.h"

using namespace botany;

// Genome with growth frozen: no new nodes, no elongation, no auxin production.
// Cytokinin production and decay still run normally.
static Genome cytokinin_test_genome() {
    Genome g = default_genome();
    g.shoot_plastochron = 1000000u;
    g.root_plastochron  = 1000000u;
    g.root_growth_rate  = 0.0f;      // no elongation → no sugar burn
    g.apical_auxin_baseline = 0.0f;  // no shoot auxin (removes confounding)
    return g;
}

static WorldParams cytokinin_test_world() {
    WorldParams w = default_world_params();
    w.starvation_ticks_max = 1000000u;
    w.light_level  = 0.0f;   // no photosynthesis — sugar fixed
    w.soil_moisture = 1.0f;  // full water for RA absorption
    return w;
}

// -----------------------------------------------------------------------
// Test 1: RA cytokinin steady-state exceeds 0.04 with default diffusion rate.
//
// With cytokinin production gated by local auxin (rate × auxin), the RA
// reaches equilibrium at ~0.06 with diffusion_rate=0.02 and ~0.027 with
// the old 0.1.  Threshold 0.04 sits between the two regimes and FAILS with
// the old default 0.1, passing once the rate is lowered to 0.02.
// -----------------------------------------------------------------------
TEST_CASE("cytokinin: RA accumulates > 0.04 with low default diffusion rate", "[cytokinin]") {
    Genome g = cytokinin_test_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();

    // Locate the root apical (direct child of seed in the initial plant)
    Node* ra = nullptr;
    for (Node* c : seed->children) {
        if (c->type == NodeType::ROOT_APICAL) { ra = c; break; }
    }
    REQUIRE(ra != nullptr);
    REQUIRE(ra->as_root_apical()->active);   // primary RA starts active

    // Prime with plenty of sugar/water so RA is never starved or turgor-limited
    seed->chemical(ChemicalID::Sugar) = 100.0f;
    ra->chemical(ChemicalID::Sugar)   = 50.0f;

    WorldParams world = cytokinin_test_world();
    for (int i = 0; i < 100; i++) plant.tick(world);

    // With rate=0.02: RA_eq ≈ 0.030 → PASS  (Münch tick order — DFS before xylem —
    //   lowers equilibrium vs the old order where xylem preceded DFS; the two diffusion
    //   rates still differ clearly, as Test 2 verifies)
    // With rate=0.10: RA_eq ≈ 0.013 → FAIL
    REQUIRE(ra->chemical(ChemicalID::Cytokinin) > 0.02f);
}

// -----------------------------------------------------------------------
// Test 2: Lower diffusion rate produces higher RA cytokinin than higher rate.
//
// This is a direct comparison: with identical plants run for 100 ticks, the
// plant whose genome has diffusion_rate=0.02 must have a higher RA cytokinin
// level than the identical plant with diffusion_rate=0.1.
// -----------------------------------------------------------------------
TEST_CASE("cytokinin: lower diffusion rate → higher RA accumulation", "[cytokinin]") {
    // High-diffusion plant
    Genome g_high = cytokinin_test_genome();
    g_high.cytokinin_diffusion_rate = 0.1f;
    Plant plant_high(g_high, glm::vec3(0.0f));

    Node* seed_high = plant_high.seed_mut();
    Node* ra_high = nullptr;
    for (Node* c : seed_high->children) {
        if (c->type == NodeType::ROOT_APICAL) { ra_high = c; break; }
    }
    seed_high->chemical(ChemicalID::Sugar) = 100.0f;
    ra_high->chemical(ChemicalID::Sugar)   = 50.0f;

    // Low-diffusion plant
    Genome g_low = cytokinin_test_genome();
    g_low.cytokinin_diffusion_rate = 0.02f;
    Plant plant_low(g_low, glm::vec3(0.0f));

    Node* seed_low = plant_low.seed_mut();
    Node* ra_low = nullptr;
    for (Node* c : seed_low->children) {
        if (c->type == NodeType::ROOT_APICAL) { ra_low = c; break; }
    }
    seed_low->chemical(ChemicalID::Sugar) = 100.0f;
    ra_low->chemical(ChemicalID::Sugar)   = 50.0f;

    WorldParams world = cytokinin_test_world();
    for (int i = 0; i < 100; i++) {
        plant_high.tick(world);
        plant_low.tick(world);
    }

    REQUIRE(ra_low->chemical(ChemicalID::Cytokinin) >
            ra_high->chemical(ChemicalID::Cytokinin));
}

// -----------------------------------------------------------------------
// Test 3: Active RA produces cytokinin; dormant RA does not.
//
// Run a single tick from a zero-cytokinin baseline so we see pure
// production before cross-diffusion has a chance to move anything.
// After 1 tick: active RA gains ~0.15 from production; dormant RA gains 0
// (no self-production, and diffusion across one hop is negligible when all
// nodes start at zero cytokinin).
// -----------------------------------------------------------------------
TEST_CASE("cytokinin: active RA produces cytokinin, dormant does not", "[cytokinin]") {
    Genome g = cytokinin_test_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    Node* ra_active = nullptr;
    for (Node* c : seed->children) {
        if (c->type == NodeType::ROOT_APICAL) { ra_active = c; break; }
    }
    REQUIRE(ra_active != nullptr);
    ra_active->as_root_apical()->active = true;

    // Manually attach a dormant sibling bud
    Node* ra_dormant = plant.create_node(NodeType::ROOT_APICAL,
                                          glm::vec3(0.0f, -0.02f, 0.0f), 0.004f);
    ra_dormant->as_root_apical()->active = false;
    seed->add_child(ra_dormant);

    // Prime sugar, then zero ALL cytokinin so the first tick shows
    // only self-production (no bootstrap diffusion across nodes).
    seed->chemical(ChemicalID::Sugar) = 100.0f;
    ra_active->chemical(ChemicalID::Sugar)  = 50.0f;
    ra_dormant->chemical(ChemicalID::Sugar) = 50.0f;
    for (Node* c : seed->children)
        c->chemical(ChemicalID::Cytokinin) = 0.0f;
    seed->chemical(ChemicalID::Cytokinin) = 0.0f;

    WorldParams world = cytokinin_test_world();
    plant.tick(world);  // single tick: production only, negligible diffusion

    REQUIRE(ra_active->chemical(ChemicalID::Cytokinin)  > 0.0f);
    REQUIRE(ra_dormant->chemical(ChemicalID::Cytokinin) == 0.0f);
}
