// tests/test_vascularization.cpp — Integration tests for vascular-driven thickening,
// canalization ratchet, and bias-weighted local distribution.
#include <catch2/catch_test_macros.hpp>
#include "engine/plant.h"
#include "engine/world_params.h"
#include "engine/sugar.h"
#include "engine/node/tissues/leaf.h"

using namespace botany;

// Genome with all autonomous growth disabled so tests control the topology and
// chemicals directly.  cambium_responsiveness is kept at default so thickening
// can happen when structural bias IS present.
static Genome frozen_genome() {
    Genome g = default_genome();
    g.apical_auxin_baseline          = 0.0f;  // no shoot auxin production
    g.root_tip_auxin_production_rate = 0.0f;  // no root auxin production
    g.root_cytokinin_production_rate = 0.0f;  // no root cytokinin production
    g.shoot_plastochron              = 1000000u;  // meristems never spawn
    g.root_plastochron               = 1000000u;
    return g;
}

// World with photosynthesis and starvation death both disabled so sugar levels
// only change via transport, and nodes never die during the test.
static WorldParams static_world() {
    WorldParams w = default_world_params();
    w.starvation_ticks_max = 1000000u;
    w.light_level          = 0.0f;  // no photosynthesis — sugar moves only via transport
    return w;
}

// -----------------------------------------------------------------------
// Test 1: Zero structural_flow_bias → zero thickening
//
// A manually-attached STEM with no bias entry in seed's map should never
// thicken regardless of sugar supply.  With frozen auxin production the
// structural bias cannot accumulate either, so the radius must stay exactly
// at its initial value after 50 ticks.
// -----------------------------------------------------------------------
TEST_CASE("Vascularization: no structural bias means no thickening", "[vascularization]") {
    Genome g = frozen_genome();
    // cambium_responsiveness is non-zero by default — thickening WOULD happen
    // if there were any structural bias.

    Plant plant(g, glm::vec3(0.0f));
    Node* seed = plant.seed_mut();
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.1f, 0.0f), 0.05f);
    seed->add_child(stem);

    // No structural_flow_bias entry set — seed->structural_flow_bias[stem] = 0.
    // With no auxin production, update_canalization() sees zero flux each tick,
    // so the structural bias never rises above the 1e-6 threshold in thicken().
    stem->chemical(ChemicalID::Sugar) = 10.0f;
    seed->chemical(ChemicalID::Sugar) = 10.0f;

    float radius_before = stem->radius;

    WorldParams world = static_world();
    for (int i = 0; i < 50; i++) plant.tick(world);

    REQUIRE(stem->radius == radius_before);
}

// -----------------------------------------------------------------------
// Test 2: Higher structural bias → faster thickening
//
// Two sibling STEM nodes start with the same radius and ample sugar but
// different structural_flow_bias values on the parent (seed) side.  After
// 50 ticks the high-bias stem must have a larger radius gain.
// -----------------------------------------------------------------------
TEST_CASE("Vascularization: higher bias produces faster thickening", "[vascularization]") {
    Genome g = frozen_genome();

    Plant plant(g, glm::vec3(0.0f));
    Node* seed   = plant.seed_mut();
    Node* stem_a = plant.create_node(NodeType::STEM, glm::vec3( 0.1f, 0.1f, 0.0f), 0.05f);
    Node* stem_b = plant.create_node(NodeType::STEM, glm::vec3(-0.1f, 0.1f, 0.0f), 0.05f);
    seed->add_child(stem_a);
    seed->add_child(stem_b);

    // 5:1 bias ratio gives a clear, measurable thickening difference.
    // With no auxin production the biases are stable across all 50 ticks.
    seed->structural_flow_bias[stem_a] = 0.5f;
    seed->structural_flow_bias[stem_b] = 0.1f;

    // Sugar is far above the maintenance threshold so sugar_gf ≈ 1 for both.
    stem_a->chemical(ChemicalID::Sugar) = 20.0f;
    stem_b->chemical(ChemicalID::Sugar) = 20.0f;
    seed->chemical(ChemicalID::Sugar)   = 20.0f;

    float a_initial = stem_a->radius;
    float b_initial = stem_b->radius;

    WorldParams world = static_world();
    for (int i = 0; i < 50; i++) plant.tick(world);

    float a_growth = stem_a->radius - a_initial;
    float b_growth = stem_b->radius - b_initial;

    REQUIRE(b_growth > 0.0f);     // low-bias stem also thickens (bias=0.1 > 1e-6)
    REQUIRE(a_growth > b_growth); // high-bias stem thickens proportionally more
}

// -----------------------------------------------------------------------
// Test 3: Canalization ratchet: auxin flux builds structural_flow_bias
//
// Manually seeded auxin on a stem flows basipetally toward the seed.
// Each tick the flow is recorded in seed->last_auxin_flux[stem] and if it
// exceeds structural_threshold the permanent bias ratchets up by
// structural_growth_rate.  After 100 ticks the bias must be > 0.
// -----------------------------------------------------------------------
TEST_CASE("Vascularization: canalization ratchet builds from auxin flow", "[vascularization]") {
    Genome g = frozen_genome();
    g.cambium_responsiveness = 0.0f;  // thickening off — don't consume test sugar

    Plant plant(g, glm::vec3(0.0f));
    Node* seed = plant.seed_mut();
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.1f, 0.0f), 0.05f);
    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.2f, 0.0f), 0.05f);
    seed->add_child(stem);
    stem->add_child(leaf);

    // No bias at start — canalization must build it from scratch.
    REQUIRE(seed->structural_flow_bias.count(stem) == 0);

    // Seed auxin at the stem end.  With basipetal bias (auxin_bias = -0.1)
    // this flows toward the seed, registering as flux on the seed↔stem edge.
    // A large initial amount sustains flux above structural_threshold long enough
    // to build measurable bias before the auxin decays away.
    stem->chemical(ChemicalID::Auxin) = 5.0f;

    // Ample sugar and a full-sized leaf prevent starvation interfering.
    seed->chemical(ChemicalID::Sugar) = 10.0f;
    stem->chemical(ChemicalID::Sugar) = 10.0f;
    leaf->chemical(ChemicalID::Sugar) = 10.0f;
    leaf->as_leaf()->leaf_size        = 0.3f;

    WorldParams world = static_world();
    for (int i = 0; i < 100; i++) plant.tick(world);

    // The structural bias on the seed→stem connection must have grown from auxin flow.
    REQUIRE(seed->structural_flow_bias.count(stem) > 0);
    REQUIRE(seed->structural_flow_bias.at(stem) > 0.0f);
}

// -----------------------------------------------------------------------
// Test 4: Bias-weighted distribution favors the high-bias branch's leaf
//
// Two seed→stem→leaf branches with 2:1 structural_flow_bias.  Setting
// vascular_conductance_threshold above the bias values forces all sugar
// transport through bias-weighted local diffusion rather than the vascular
// bulk pass.  After 20 ticks the high-bias branch's leaf must hold more
// sugar than the low-bias one.
// -----------------------------------------------------------------------
TEST_CASE("Vascularization: conductance-weighted distribution favors high-bias branch", "[vascularization]") {
    Genome g = frozen_genome();
    g.cambium_responsiveness        = 0.0f;
    // Raise the vascular admission threshold so our manually-set biases (≤ 2.0)
    // keep the stems out of the vascular bulk pass.  All sugar then moves via
    // local diffusion, which weights each child's share by bias_mult =
    // 1 + canalization_weight * structural_flow_bias.
    g.vascular_conductance_threshold = 100.0f;

    Plant plant(g, glm::vec3(0.0f));
    Node* seed   = plant.seed_mut();
    Node* stem_a = plant.create_node(NodeType::STEM, glm::vec3( 0.1f, 0.1f, 0.0f), 0.05f);
    Node* stem_b = plant.create_node(NodeType::STEM, glm::vec3(-0.1f, 0.1f, 0.0f), 0.05f);
    Node* leaf_a = plant.create_node(NodeType::LEAF, glm::vec3( 0.1f, 0.2f, 0.0f), 0.05f);
    Node* leaf_b = plant.create_node(NodeType::LEAF, glm::vec3(-0.1f, 0.2f, 0.0f), 0.05f);
    seed->add_child(stem_a);
    seed->add_child(stem_b);
    stem_a->add_child(leaf_a);
    stem_b->add_child(leaf_b);

    // 2:1 bias ratio → bias_mult_a = 3.0, bias_mult_b = 2.0.
    // stem_a wins a disproportionate share of seed's sugar every tick.
    seed->structural_flow_bias[stem_a] = 2.0f;
    seed->structural_flow_bias[stem_b] = 1.0f;

    // Sugar source: seed.  Both leaves start empty and act as distal sinks.
    seed->chemical(ChemicalID::Sugar) = 10.0f;
    leaf_a->as_leaf()->leaf_size      = 0.3f;
    leaf_b->as_leaf()->leaf_size      = 0.3f;

    WorldParams world = static_world();
    for (int i = 0; i < 20; i++) plant.tick(world);

    REQUIRE(leaf_a->chemical(ChemicalID::Sugar) > leaf_b->chemical(ChemicalID::Sugar));
}
