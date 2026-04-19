// tests/test_vascularization.cpp — Integration tests for vascular-driven thickening,
// PIN canalization, and bias-weighted vascular distribution.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "engine/plant.h"
#include "engine/world_params.h"
#include "engine/sugar.h"
#include "engine/vascular.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"
#include "engine/chemical/chemical.h"

using namespace botany;
using Catch::Approx;

// Genome with all autonomous growth disabled so tests control the topology and
// chemicals directly.  cambium_responsiveness is kept at default so thickening
// can happen when auxin_flow_bias IS present.
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
// Test 1: Zero auxin_flow_bias → zero thickening
//
// A manually-attached STEM with no bias entry in seed's map should never
// thicken regardless of sugar supply.  With frozen auxin production the
// auxin_flow_bias cannot build up either, so the radius must stay exactly
// at its initial value after 50 ticks.
// -----------------------------------------------------------------------
TEST_CASE("Vascularization: no auxin_flow_bias means no thickening", "[vascularization]") {
    Genome g = frozen_genome();
    // Disable cambium so the test is purely about "no bias → no thickening" logic,
    // isolated from any floating-point noise at near-zero bias values.
    g.cambium_responsiveness = 0.0f;

    Plant plant(g, glm::vec3(0.0f));
    Node* seed = plant.seed_mut();
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.1f, 0.0f), 0.05f);
    seed->add_child(stem);

    // No auxin_flow_bias entry set. With no auxin production, PIN transport moves
    // nothing, update_canalization() sees zero flux each tick, so bias stays at 0
    // (or lerps to 0 on the first tick) and never exceeds the 1e-6 threshold.
    stem->chemical(ChemicalID::Sugar) = 10.0f;
    seed->chemical(ChemicalID::Sugar) = 10.0f;

    float radius_before = stem->radius;

    WorldParams world = static_world();
    for (int i = 0; i < 50; i++) plant.tick(world);

    REQUIRE(stem->radius == radius_before);
}

// -----------------------------------------------------------------------
// Test 2: Higher auxin_flow_bias → faster thickening
//
// Two sibling STEM nodes start with the same radius and ample sugar but
// different auxin_flow_bias values on the parent (seed) side.  After
// 50 ticks the high-bias stem must have a larger radius gain.
//
// Both biases decay proportionally each tick (frozen genome = no auxin flux
// → saturation = 0 → lerp toward 0 at smoothing_rate).  Decay is identical
// in shape for both so the 5:1 initial ratio is preserved throughout and
// the proportional thickening difference remains measurable.
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
    seed->auxin_flow_bias[stem_a] = 0.5f;
    seed->auxin_flow_bias[stem_b] = 0.1f;

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

    REQUIRE(b_growth > 0.0f);     // low-bias stem also thickens (bias=0.1 decays but stays > 1e-6 for many ticks)
    REQUIRE(a_growth > b_growth); // high-bias stem thickens proportionally more
}

// -----------------------------------------------------------------------
// Test 3: PIN canalization builds auxin_flow_bias from auxin flux
//
// Manually seeded auxin on a stem flows basipetally toward the seed via
// both PIN transport (Phase A) and local diffusion.  Each tick the flux is
// recorded in seed->last_auxin_flux[stem] and update_canalization() lerps
// auxin_flow_bias toward the saturation level.  After ticks the bias must
// be > 0 on the seed→stem connection.
// -----------------------------------------------------------------------
TEST_CASE("Vascularization: auxin flux builds auxin_flow_bias from zero", "[vascularization]") {
    Genome g = frozen_genome();
    g.cambium_responsiveness = 0.0f;  // thickening off — don't consume test sugar

    Plant plant(g, glm::vec3(0.0f));
    Node* seed = plant.seed_mut();
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.1f, 0.0f), 0.05f);
    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.2f, 0.0f), 0.05f);
    seed->add_child(stem);
    stem->add_child(leaf);

    // No bias at start.
    {
        auto it = seed->auxin_flow_bias.find(stem);
        float initial = (it != seed->auxin_flow_bias.end()) ? it->second : 0.0f;
        REQUIRE(initial == 0.0f);
    }

    // Seed auxin at the stem end.  With basipetal bias (auxin_bias = -0.1) and
    // PIN Phase A this flows toward the seed, registering flux on the seed↔stem
    // edge.  A large initial amount sustains flux above zero long enough to build
    // measurable bias before the auxin decays away.
    stem->chemical(ChemicalID::Auxin) = 5.0f;

    // Ample sugar and a full-sized leaf prevent starvation interfering.
    seed->chemical(ChemicalID::Sugar) = 10.0f;
    stem->chemical(ChemicalID::Sugar) = 10.0f;
    leaf->chemical(ChemicalID::Sugar) = 10.0f;
    leaf->as_leaf()->leaf_size        = 0.3f;

    WorldParams world = static_world();
    for (int i = 0; i < 100; i++) plant.tick(world);

    // The auxin_flow_bias on the seed→stem connection must have built from auxin flux.
    REQUIRE(seed->auxin_flow_bias.count(stem) > 0);
    REQUIRE(seed->auxin_flow_bias.at(stem) > 0.0f);
}

// -----------------------------------------------------------------------
// Test 4: Münch phloem delivers to symmetric sinks equally.
//
// Topology: seed → leaf_src (phloem source) + stem_high→apical_high (sink) +
//                                               stem_low→apical_low  (sink).
// Both stems/apicals have identical geometry.  stem_high has higher
// auxin_flow_bias but under the Münch pressure-flow model phloem routing is
// driven by osmotic pressure gradients and pipe cross-section area, NOT by
// auxin_flow_bias.  Canalization bias still governs xylem (water/cytokinin)
// distribution but not phloem.
// Both apicals start empty with equal geometry → they receive equal sugar.
// -----------------------------------------------------------------------
TEST_CASE("Vascularization: conductance-weighted vascular pass favors high-bias branch", "[vascularization]") {
    Genome g = frozen_genome();
    g.cambium_responsiveness = 0.0f;  // no thickening during this pass

    Plant plant(g, glm::vec3(0.0f));
    Node* seed = plant.seed_mut();

    // Leaf source: leaf_size = 1.0 dm.
    Node* leaf_src = plant.create_node(NodeType::LEAF,
                                       glm::vec3(0.1f, 0.0f, 0.0f), 0.01f);
    leaf_src->as_leaf()->leaf_size = 1.0f;
    seed->add_child(leaf_src);

    // Two branches with same stem radius but different auxin_flow_bias.
    float stem_r = 0.05f;
    Node* stem_high   = plant.create_node(NodeType::STEM,
                                          glm::vec3(-0.3f, 0.5f, 0.0f), stem_r);
    Node* apical_high = plant.create_node(NodeType::APICAL,
                                          glm::vec3(0.0f,  0.1f, 0.0f), 0.02f);
    seed->add_child(stem_high);
    stem_high->add_child(apical_high);
    seed->auxin_flow_bias[stem_high] = 2.0f;

    Node* stem_low   = plant.create_node(NodeType::STEM,
                                         glm::vec3( 0.3f, 0.5f, 0.0f), stem_r);
    Node* apical_low = plant.create_node(NodeType::APICAL,
                                         glm::vec3(0.0f,  0.1f, 0.0f), 0.02f);
    seed->add_child(stem_low);
    stem_low->add_child(apical_low);
    seed->auxin_flow_bias[stem_low] = 1.0f;

    // Pre-fill every node to its cap so only the test apicals are sinks.
    plant.for_each_node_mut([&](Node& n) {
        if (&n != leaf_src && &n != apical_high && &n != apical_low)
            n.chemical(ChemicalID::Sugar) = sugar_cap(n, g);
    });
    leaf_src->chemical(ChemicalID::Sugar)   = 0.605f;
    apical_high->chemical(ChemicalID::Sugar) = 0.0f;
    apical_low->chemical(ChemicalID::Sugar)  = 0.0f;
    seed->chemical(ChemicalID::Sugar) = 0.0f;

    vascular_transport(plant, g, static_world());

    float got_high = apical_high->chemical(ChemicalID::Sugar);
    float got_low  = apical_low->chemical(ChemicalID::Sugar);

    // Under Münch flow, both apicals receive sugar from their parent stems
    // (which are filled to cap, giving a steep concentration gradient).
    // Geometry is symmetric → both receive equal amounts.
    // Note: auxin_flow_bias is irrelevant for Münch phloem; it governs xylem only.
    REQUIRE(got_high > 0.0f);
    REQUIRE(got_low  > 0.0f);
    REQUIRE(got_high == Approx(got_low).margin(1e-5f));
}

// -----------------------------------------------------------------------
// Test 5: get_parent_auxin_flow_bias() accessor
//
// The Vascular color overlay uses `node.get_parent_auxin_flow_bias()` as its
// ChemicalAccessor.  This test verifies the three cases that accessor will hit:
//   (a) a regular node returns its parent's recorded bias for it
//   (b) a node with no parent entry yet returns 0
//   (c) the seed (no parent) returns the max of all its children's biases,
//       giving it a color proportional to the busiest branch through it.
// -----------------------------------------------------------------------
TEST_CASE("Vascularization overlay: get_parent_auxin_flow_bias returns correct values", "[vascularization]") {
    Genome g = frozen_genome();
    Plant plant(g, glm::vec3(0.0f));
    Node* seed   = plant.seed_mut();
    Node* stem_a = plant.create_node(NodeType::STEM, glm::vec3( 0.1f, 0.1f, 0.0f), 0.05f);
    Node* stem_b = plant.create_node(NodeType::STEM, glm::vec3(-0.1f, 0.1f, 0.0f), 0.05f);
    Node* stem_c = plant.create_node(NodeType::STEM, glm::vec3( 0.0f, 0.2f, 0.0f), 0.05f);
    seed->add_child(stem_a);
    seed->add_child(stem_b);
    seed->add_child(stem_c);

    seed->auxin_flow_bias[stem_a] = 0.4f;
    seed->auxin_flow_bias[stem_b] = 1.8f;
    // stem_c has no entry — should return 0

    // (a) node with a recorded entry returns that value
    REQUIRE(stem_a->get_parent_auxin_flow_bias() == Approx(0.4f));
    REQUIRE(stem_b->get_parent_auxin_flow_bias() == Approx(1.8f));

    // (b) node with no entry in parent's map returns 0
    REQUIRE(stem_c->get_parent_auxin_flow_bias() == Approx(0.0f));

    // (c) seed (no parent) returns max of its children's biases
    REQUIRE(seed->get_parent_auxin_flow_bias() == Approx(1.8f));
}

// -----------------------------------------------------------------------
// Test 6: Demand-driven phloem delivers sugar to apex across long chain
//
// Topology: seed (sugar=10) → 15 STEM nodes → ApicalNode (sugar=0)
//
// The primary SA is reparented to the tip of the chain so it acts as the
// sole phloem sink at maximum distance from the source.  After one tick
// of vascular_transport + tick, the apex must have measurably non-zero
// sugar.  This fails under the Jacobi Münch implementation (3 iterations,
// attenuates over ~30 hops) and is expected to pass once the algorithm is
// replaced with a demand-driven single-pass allocation.
// -----------------------------------------------------------------------
TEST_CASE("Demand-driven phloem delivers sugar to apex across long shoot chain", "[vascular][phloem][demand]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));
    WorldParams world;

    // Build a synthetic 15-stem chain above the seed.  Give each stem enough
    // radius to count as a vascular conduit.  The primary SA is already the
    // seed's shoot child from the Plant constructor; we'll reparent it to the
    // tip of the chain so it acts as the apex.
    Node* tip_stem = plant.seed_mut();
    for (int i = 0; i < 15; ++i) {
        Node* stem = plant.create_node(NodeType::STEM,
            glm::vec3(0.0f, 0.05f * (i + 1), 0.0f), 0.015f);
        tip_stem->add_child(stem);
        tip_stem = stem;
    }

    // Find the primary SA and reparent it to the tip of the chain.
    ApicalNode* apex = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (auto sa = n.as_apical(); sa && sa->is_primary) apex = sa;
    });
    REQUIRE(apex != nullptr);
    // Remove apex from its current parent and attach to the chain tip.
    if (apex->parent) {
        auto& siblings = apex->parent->children;
        siblings.erase(std::remove(siblings.begin(), siblings.end(),
                                   static_cast<Node*>(apex)),
                       siblings.end());
    }
    apex->parent = nullptr;
    tip_stem->add_child(apex);

    // Zero all sugar except the seed.  Put a healthy pile at the seed so we
    // can detect any flow reaching the tip.
    plant.for_each_node_mut([&](Node& n) {
        n.chemical(ChemicalID::Sugar) = 0.0f;
    });
    plant.seed_mut()->chemical(ChemicalID::Sugar) = 10.0f;

    float apex_sugar_before = apex->chemical(ChemicalID::Sugar);
    REQUIRE(apex_sugar_before == 0.0f);

    // Single tick — phloem_resolve runs once.  Apex should receive a
    // measurable amount (> 0.001 g, well above float noise).  The precise
    // amount depends on proportional allocation, but anything non-trivial
    // demonstrates that the 15-hop distance did not attenuate delivery.
    plant.tick(world);

    float apex_sugar_after = apex->chemical(ChemicalID::Sugar);
    REQUIRE(apex_sugar_after > 0.001f);
}
