// tests/test_vascularization.cpp — Integration tests for vascular-driven thickening,
// canalization ratchet, and bias-weighted vascular distribution.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "engine/plant.h"
#include "engine/world_params.h"
#include "engine/sugar.h"
#include "engine/vascular.h"
#include "engine/node/tissues/leaf.h"
#include "engine/chemical/chemical.h"

using namespace botany;
using Catch::Approx;

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
// Test 4: Conductance-first vascular distribution favors the high-bias branch
//
// Topology: seed → leaf_src (phloem source) + stem_high→apical_high (sink) +
//                                               stem_low→apical_low  (sink).
// stem_high has 2× the structural_flow_bias of stem_low (both above
// vascular_conductance_threshold so both are vascular conduits).
// Both apicals start empty, so they have equal demand.  leaf_src provides
// limited sugar (surplus < combined demand), forcing the water-filling loop
// to ration by conductance weight: weight = pipe_cap × (1 + bias).
// The high-bias branch (weight 3×) receives more than the low-bias (weight 2×).
// -----------------------------------------------------------------------
TEST_CASE("Vascularization: conductance-weighted vascular pass favors high-bias branch", "[vascularization]") {
    Genome g = frozen_genome();
    g.cambium_responsiveness = 0.0f;  // no thickening during this pass

    Plant plant(g, glm::vec3(0.0f));
    Node* seed = plant.seed_mut();

    // Leaf source: leaf_size = 1.0 → sugar_cap = 2.0, reserve = 0.6.
    // Setting sugar = 0.7 → surplus = 0.1, which is less than the combined
    // pipe-capped demand of the two stems (~0.063 each), so rationing applies.
    Node* leaf_src = plant.create_node(NodeType::LEAF,
                                       glm::vec3(0.1f, 0.0f, 0.0f), 0.01f);
    leaf_src->as_leaf()->leaf_size = 1.0f;
    seed->add_child(leaf_src);

    // Two branches with same stem radius but different structural bias.
    float stem_r = 0.05f;
    Node* stem_high   = plant.create_node(NodeType::STEM,
                                          glm::vec3(-0.3f, 0.5f, 0.0f), stem_r);
    Node* apical_high = plant.create_node(NodeType::APICAL,
                                          glm::vec3(0.0f,  0.1f, 0.0f), 0.02f);
    seed->add_child(stem_high);
    stem_high->add_child(apical_high);
    seed->structural_flow_bias[stem_high] = 2.0f;  // above vascular_conductance_threshold

    Node* stem_low   = plant.create_node(NodeType::STEM,
                                         glm::vec3( 0.3f, 0.5f, 0.0f), stem_r);
    Node* apical_low = plant.create_node(NodeType::APICAL,
                                         glm::vec3(0.0f,  0.1f, 0.0f), 0.02f);
    seed->add_child(stem_low);
    stem_low->add_child(apical_low);
    seed->structural_flow_bias[stem_low] = 1.0f;   // above vascular_conductance_threshold

    // Pre-fill every default node to its cap so their sugar demand = 0.
    // This isolates our two test apicals as the only active phloem sinks.
    plant.for_each_node_mut([&](Node& n) {
        if (&n != leaf_src && &n != apical_high && &n != apical_low)
            n.chemical(ChemicalID::Sugar) = sugar_cap(n, g);
    });
    // Set the leaf source and empty sinks after the fill.
    leaf_src->chemical(ChemicalID::Sugar)   = 0.7f;
    apical_high->chemical(ChemicalID::Sugar) = 0.0f;
    apical_low->chemical(ChemicalID::Sugar)  = 0.0f;

    // Single vascular transport pass — tests the distribution algorithm directly.
    vascular_transport(plant, g, static_world());

    float got_high = apical_high->chemical(ChemicalID::Sugar);
    float got_low  = apical_low->chemical(ChemicalID::Sugar);

    // Both branches received sugar (budget was non-zero), but the high-bias
    // branch (weight 3×) received proportionally more than the low-bias (2×).
    REQUIRE(got_high > 0.0f);
    REQUIRE(got_low  > 0.0f);
    REQUIRE(got_high > got_low);
}

// -----------------------------------------------------------------------
// Test 5: get_parent_structural_bias() accessor used by the Vascular overlay
//
// The Vascular color overlay uses `node.get_parent_structural_bias()` as its
// ChemicalAccessor.  This test verifies the three cases that accessor will hit:
//   (a) a regular node returns its parent's recorded bias for it
//   (b) a node with no parent entry yet returns 0
//   (c) the seed (no parent) returns the max of all its children's biases,
//       giving it a color proportional to the busiest branch through it.
// -----------------------------------------------------------------------
TEST_CASE("Vascularization overlay: get_parent_structural_bias returns correct values", "[vascularization]") {
    Genome g = frozen_genome();
    Plant plant(g, glm::vec3(0.0f));
    Node* seed   = plant.seed_mut();
    Node* stem_a = plant.create_node(NodeType::STEM, glm::vec3( 0.1f, 0.1f, 0.0f), 0.05f);
    Node* stem_b = plant.create_node(NodeType::STEM, glm::vec3(-0.1f, 0.1f, 0.0f), 0.05f);
    Node* stem_c = plant.create_node(NodeType::STEM, glm::vec3( 0.0f, 0.2f, 0.0f), 0.05f);
    seed->add_child(stem_a);
    seed->add_child(stem_b);
    seed->add_child(stem_c);

    seed->structural_flow_bias[stem_a] = 0.4f;
    seed->structural_flow_bias[stem_b] = 1.8f;
    // stem_c has no entry — should return 0

    // (a) node with a recorded entry returns that value
    REQUIRE(stem_a->get_parent_structural_bias() == Approx(0.4f));
    REQUIRE(stem_b->get_parent_structural_bias() == Approx(1.8f));

    // (b) node with no entry in parent's map returns 0
    REQUIRE(stem_c->get_parent_structural_bias() == Approx(0.0f));

    // (c) seed (no parent) returns max of its children's biases
    REQUIRE(seed->get_parent_structural_bias() == Approx(1.8f));
}
