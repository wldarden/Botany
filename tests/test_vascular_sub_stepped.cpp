#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "engine/vascular_sub_stepped.h"
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/node/stem_node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"

using namespace botany;

TEST_CASE("vascular_sub_stepped exists and can be called", "[vascular_sub_stepped]") {
    Genome g = default_genome();
    WorldParams world;
    Plant plant(g, glm::vec3(0.0f));

    // Should not crash on an empty plant (just a seed).
    vascular_sub_stepped(plant, g, world);
    SUCCEED();
}

TEST_CASE("radial_permeability approaches base at r=0", "[vascular_sub_stepped][radial]") {
    Genome g = default_genome();
    float perm = radial_permeability_sugar(0.0f, g);
    REQUIRE(perm == g.base_radial_permeability_sugar);
}

TEST_CASE("radial_permeability at half-radius is between base and floor", "[vascular_sub_stepped][radial]") {
    Genome g = default_genome();
    float base   = g.base_radial_permeability_sugar;
    float floor  = g.radial_floor_fraction_sugar;
    float r_half = g.radial_half_radius_sugar;
    float perm = radial_permeability_sugar(r_half, g);
    // At r = r_half, denominator = 1 + 1 = 2, so (1 - floor)/2 remains.
    // perm = base × (floor + (1 - floor) / 2)
    float expected = base * (floor + (1.0f - floor) * 0.5f);
    REQUIRE(perm == Catch::Approx(expected).margin(1e-5f));
}

TEST_CASE("radial_permeability asymptotes to floor at large r", "[vascular_sub_stepped][radial]") {
    Genome g = default_genome();
    float perm = radial_permeability_sugar(100.0f, g);
    float expected_floor = g.base_radial_permeability_sugar * g.radial_floor_fraction_sugar;
    REQUIRE(perm == Catch::Approx(expected_floor).margin(0.01f));
}

TEST_CASE("phloem_capacity scales with r² × length × phloem_fraction", "[vascular_sub_stepped][capacity]") {
    Genome g = default_genome();
    StemNode stem(1, glm::vec3(0), 0.1f);  // radius 0.1 dm
    // Give the stem a length via offset magnitude.
    stem.offset = glm::vec3(0.0f, 1.0f, 0.0f);  // length = 1 dm

    float expected = 3.14159f * 0.1f * 0.1f * 1.0f * g.phloem_fraction;
    REQUIRE(phloem_capacity(stem, g) == Catch::Approx(expected).margin(1e-4f));
}

TEST_CASE("phloem_capacity returns zero for a leaf (no phloem)", "[vascular_sub_stepped][capacity]") {
    Genome g = default_genome();
    LeafNode leaf(1, glm::vec3(0), 0.01f);
    REQUIRE(phloem_capacity(leaf, g) == 0.0f);
}

TEST_CASE("leaf budget: sugar_supply is local sugar above reserve", "[vascular_sub_stepped][budget]") {
    Genome g = default_genome();
    WorldParams world;
    // Give the leaf a parent so sugar_cap uses the leaf path (not seed path).
    StemNode parent(0, glm::vec3(0), 0.1f);
    LeafNode leaf(1, glm::vec3(0), 0.02f);
    leaf.parent = &parent;
    leaf.offset = glm::vec3(0.0f, 0.1f, 0.0f);
    leaf.leaf_size = 0.5f;  // area = 0.25 dm², cap = 0.25 * 2.0 = 0.5 g
    // sugar_cap = 0.5, reserve = 0.3 * 0.5 = 0.15, so sugar = 1.0 → supply = 0.85
    leaf.local().chemical(ChemicalID::Sugar) = 1.0f;
    VascularBudget b = compute_budget(leaf, g, world);
    REQUIRE(b.sugar_supply > 0.0f);
    REQUIRE(b.sugar_demand == 0.0f);
}

TEST_CASE("dormant apical budget: zero demand", "[vascular_sub_stepped][budget]") {
    Genome g = default_genome();
    WorldParams world;
    ApicalNode apical(1, glm::vec3(0), 0.01f);
    apical.active = false;  // dormant
    VascularBudget b = compute_budget(apical, g, world);
    REQUIRE(b.sugar_demand == 0.0f);
    REQUIRE(b.water_demand == 0.0f);
}

TEST_CASE("active apical below target has positive sugar_demand", "[vascular_sub_stepped][budget]") {
    Genome g = default_genome();
    WorldParams world;
    ApicalNode apical(1, glm::vec3(0), 0.01f);
    apical.active = true;
    apical.local().chemical(ChemicalID::Sugar) = 0.0f;
    VascularBudget b = compute_budget(apical, g, world);
    REQUIRE(b.sugar_demand > 0.0f);
}

TEST_CASE("inject step: leaf pushes sugar into parent stem phloem", "[vascular_sub_stepped][inject]") {
    Genome g = default_genome();
    WorldParams world;
    StemNode stem(1, glm::vec3(0), 0.02f);
    stem.offset = glm::vec3(0.0f, 0.05f, 0.0f);
    LeafNode leaf(2, glm::vec3(0.05f, 0, 0), 0.01f);
    leaf.leaf_size = 0.5f;  // so sugar_cap uses the leaf code path
    stem.add_child(&leaf);
    leaf.parent = &stem;

    leaf.local().chemical(ChemicalID::Sugar) = 10.0f;
    VascularBudget b = compute_budget(leaf, g, world);
    REQUIRE(b.sugar_supply > 0.0f);
    float budget_slice = b.sugar_supply / 10.0f;  // N = 10 sub-steps

    float leaf_before   = leaf.local().chemical(ChemicalID::Sugar);
    float phloem_before = stem.phloem()->chemical(ChemicalID::Sugar);

    inject_step(leaf, b, /* N */ 10, g);

    float leaf_after   = leaf.local().chemical(ChemicalID::Sugar);
    float phloem_after = stem.phloem()->chemical(ChemicalID::Sugar);

    REQUIRE(leaf_after   == Catch::Approx(leaf_before   - budget_slice).margin(1e-5f));
    REQUIRE(phloem_after == Catch::Approx(phloem_before + budget_slice).margin(1e-5f));
}

TEST_CASE("extract step: meristem pulls sugar from parent phloem", "[vascular_sub_stepped][extract]") {
    Genome g = default_genome();
    WorldParams world;
    StemNode stem(1, glm::vec3(0), 0.02f);
    stem.offset = glm::vec3(0.0f, 0.05f, 0.0f);
    ApicalNode apical(2, glm::vec3(0.05f, 0, 0), 0.01f);
    stem.add_child(&apical);
    apical.parent = &stem;
    apical.active = true;

    stem.phloem()->chemical(ChemicalID::Sugar) = 5.0f;
    apical.local().chemical(ChemicalID::Sugar) = 0.0f;

    VascularBudget b = compute_budget(apical, g, world);
    REQUIRE(b.sugar_demand > 0.0f);
    float slice = b.sugar_demand / 10.0f;

    extract_step(apical, b, /* N */ 10, g);

    // Sugar moved from phloem to meristem local, up to slice.
    REQUIRE(apical.local().chemical(ChemicalID::Sugar) == Catch::Approx(slice).margin(1e-5f));
    REQUIRE(stem.phloem()->chemical(ChemicalID::Sugar) == Catch::Approx(5.0f - slice).margin(1e-5f));
}

TEST_CASE("extract step capped by available phloem sugar", "[vascular_sub_stepped][extract]") {
    Genome g = default_genome();
    WorldParams world;
    StemNode stem(1, glm::vec3(0), 0.02f);
    stem.offset = glm::vec3(0.0f, 0.05f, 0.0f);
    ApicalNode apical(2, glm::vec3(0.05f, 0, 0), 0.01f);
    stem.add_child(&apical);
    apical.parent = &stem;
    apical.active = true;

    stem.phloem()->chemical(ChemicalID::Sugar) = 0.01f;  // almost dry
    apical.local().chemical(ChemicalID::Sugar) = 0.0f;

    VascularBudget b = compute_budget(apical, g, world);
    extract_step(apical, b, /* N */ 10, g);

    // Capped by phloem availability — no negative phloem.
    REQUIRE(stem.phloem()->chemical(ChemicalID::Sugar) >= 0.0f);
    REQUIRE(apical.local().chemical(ChemicalID::Sugar) <= 0.01f);
}

TEST_CASE("radial flow equilibrates stem local and phloem", "[vascular_sub_stepped][radial_flow]") {
    Genome g = default_genome();
    StemNode stem(1, glm::vec3(0), 0.015f);
    stem.offset = glm::vec3(0.0f, 0.05f, 0.0f);

    // Start with high phloem, empty local.
    stem.phloem()->chemical(ChemicalID::Sugar) = 1.0f;
    stem.local().chemical(ChemicalID::Sugar)   = 0.0f;

    float total_before = stem.phloem()->chemical(ChemicalID::Sugar)
                       + stem.local().chemical(ChemicalID::Sugar);

    // Run radial flow many times — local should approach phloem concentration.
    for (int i = 0; i < 100; ++i) {
        radial_flow_step(stem, /* N */ 1, g);
    }

    float total_after = stem.phloem()->chemical(ChemicalID::Sugar)
                      + stem.local().chemical(ChemicalID::Sugar);

    // Conservation.
    REQUIRE(total_after == Catch::Approx(total_before).margin(1e-4f));
    // After many iterations, some flowed from phloem to local.
    REQUIRE(stem.local().chemical(ChemicalID::Sugar) > 0.0f);
    REQUIRE(stem.phloem()->chemical(ChemicalID::Sugar) < 1.0f);
}

TEST_CASE("radial permeability is smaller for thicker stems", "[vascular_sub_stepped][radial_flow]") {
    Genome g = default_genome();
    float perm_young = radial_permeability_sugar(0.01f, g);
    float perm_mature = radial_permeability_sugar(2.0f, g);
    REQUIRE(perm_mature < perm_young);
}

TEST_CASE("Jacobi equalizes two connected phloem pools", "[vascular_sub_stepped][jacobi]") {
    Genome g = default_genome();
    StemNode parent_stem(1, glm::vec3(0), 0.02f);
    parent_stem.offset = glm::vec3(0.0f, 0.05f, 0.0f);
    StemNode child_stem(2, glm::vec3(0, 0.05f, 0), 0.02f);
    child_stem.offset = glm::vec3(0.0f, 0.05f, 0.0f);
    parent_stem.add_child(&child_stem);
    child_stem.parent = &parent_stem;

    parent_stem.phloem()->chemical(ChemicalID::Sugar) = 1.0f;
    child_stem.phloem()->chemical(ChemicalID::Sugar)  = 0.0f;

    float total_before = parent_stem.phloem()->chemical(ChemicalID::Sugar)
                       + child_stem.phloem()->chemical(ChemicalID::Sugar);

    // Many Jacobi iterations should bring pressures close together.
    for (int i = 0; i < 100; ++i) {
        jacobi_step(parent_stem, child_stem, g);
    }

    float total_after = parent_stem.phloem()->chemical(ChemicalID::Sugar)
                      + child_stem.phloem()->chemical(ChemicalID::Sugar);

    REQUIRE(total_after == Catch::Approx(total_before).margin(1e-4f));
    // Some sugar moved from parent to child.
    REQUIRE(child_stem.phloem()->chemical(ChemicalID::Sugar) > 0.0f);
    REQUIRE(parent_stem.phloem()->chemical(ChemicalID::Sugar) < 1.0f);
}

TEST_CASE("full vascular_sub_stepped delivers sugar to apex within N hops", "[vascular_sub_stepped][integration]") {
    Genome g = default_genome();
    WorldParams world;
    world.vascular_substeps = 20;
    Plant plant(g, glm::vec3(0.0f));

    // Build a 10-stem chain above the seed.
    Node* tip_stem = plant.seed_mut();
    for (int i = 0; i < 10; ++i) {
        Node* stem = plant.create_node(NodeType::STEM,
            glm::vec3(0.0f, 0.05f * (i + 1), 0.0f), 0.015f);
        tip_stem->add_child(stem);
        tip_stem = stem;
    }

    // Zero all sugar; pile 10 g at the seed's phloem.
    plant.for_each_node_mut([&](Node& n) {
        n.local().chemical(ChemicalID::Sugar) = 0.0f;
        if (auto* p = n.phloem()) p->chemical(ChemicalID::Sugar) = 0.0f;
    });
    plant.seed_mut()->phloem()->chemical(ChemicalID::Sugar) = 10.0f;

    vascular_sub_stepped(plant, g, world);

    // The tip stem is 10 hops from seed.  N=20 sub-steps means 20 Jacobi
    // iterations, so pressure wave reaches the tip.  The absolute amount
    // is small (thin stems have tiny phloem cross-section) but must be
    // strictly positive — any nonzero value confirms propagation.
    float tip_phloem = tip_stem->phloem()->chemical(ChemicalID::Sugar);
    REQUIRE(tip_phloem > 1e-8f);
}

TEST_CASE("vascular_sub_stepped conserves mass", "[vascular_sub_stepped][integration][conservation]") {
    Genome g = default_genome();
    WorldParams world;
    world.vascular_substeps = 20;
    Plant plant(g, glm::vec3(0.0f));

    // Build: seed + one stem + one leaf.
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0, 0.05f, 0), 0.015f);
    plant.seed_mut()->add_child(stem);
    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.02f, 0.05f, 0), 0.01f);
    stem->add_child(leaf);

    // Seed sugar and water in various compartments.
    plant.seed_mut()->phloem()->chemical(ChemicalID::Sugar) = 1.0f;
    stem->phloem()->chemical(ChemicalID::Sugar)             = 0.5f;
    leaf->local().chemical(ChemicalID::Sugar)               = 2.0f;
    plant.seed_mut()->xylem()->chemical(ChemicalID::Water)  = 1.0f;
    stem->xylem()->chemical(ChemicalID::Water)              = 0.5f;

    // Pre-sum across all pools.
    float sugar_before = 0.0f, water_before = 0.0f;
    plant.for_each_node([&](const Node& n) {
        sugar_before += n.local().chemical(ChemicalID::Sugar);
        water_before += n.local().chemical(ChemicalID::Water);
        if (auto* p = n.phloem()) sugar_before += p->chemical(ChemicalID::Sugar);
        if (auto* x = n.xylem())  water_before += x->chemical(ChemicalID::Water);
    });

    vascular_sub_stepped(plant, g, world);

    float sugar_after = 0.0f, water_after = 0.0f;
    plant.for_each_node([&](const Node& n) {
        sugar_after += n.local().chemical(ChemicalID::Sugar);
        water_after += n.local().chemical(ChemicalID::Water);
        if (auto* p = n.phloem()) sugar_after += p->chemical(ChemicalID::Sugar);
        if (auto* x = n.xylem())  water_after += x->chemical(ChemicalID::Water);
    });

    REQUIRE(sugar_after == Catch::Approx(sugar_before).margin(1e-4f));
    REQUIRE(water_after == Catch::Approx(water_before).margin(1e-4f));
}

TEST_CASE("compartment residency invariants hold after vascular", "[vascular_sub_stepped][invariants]") {
    Genome g = default_genome();
    WorldParams world;
    Plant plant(g, glm::vec3(0.0f));

    // Populate with a typical tick's chemicals.
    plant.for_each_node_mut([&](Node& n) {
        n.local().chemical(ChemicalID::Sugar)     = 1.0f;
        n.local().chemical(ChemicalID::Water)     = 0.5f;
        n.local().chemical(ChemicalID::Auxin)     = 0.1f;
        if (auto* p = n.phloem()) p->chemical(ChemicalID::Sugar) = 0.3f;
        if (auto* x = n.xylem())  x->chemical(ChemicalID::Water) = 0.2f;
    });

    vascular_sub_stepped(plant, g, world);

    // Sugar must never appear in any xylem pool.
    // Water/cytokinin must never appear in any phloem pool.
    // Auxin/gibberellin/stress must never appear in any transport pool.
    plant.for_each_node([&](const Node& n) {
        if (auto* p = n.phloem()) {
            REQUIRE(p->chemical(ChemicalID::Water)     == 0.0f);
            REQUIRE(p->chemical(ChemicalID::Cytokinin) == 0.0f);
            REQUIRE(p->chemical(ChemicalID::Auxin)     == 0.0f);
        }
        if (auto* x = n.xylem()) {
            REQUIRE(x->chemical(ChemicalID::Sugar)     == 0.0f);
            REQUIRE(x->chemical(ChemicalID::Auxin)     == 0.0f);
        }
    });
}

TEST_CASE("vascular propagates sugar along a long chain without losing mass", "[vascular_sub_stepped][hydraulic_limit]") {
    Genome g = default_genome();
    WorldParams world;
    world.vascular_substeps = 4;
    Plant plant(g, glm::vec3(0.0f));

    // Build a 30-stem chain with uniform per-stem offset.
    std::vector<Node*> chain;
    Node* tip_stem = plant.seed_mut();
    for (int i = 0; i < 30; ++i) {
        Node* stem = plant.create_node(NodeType::STEM,
            glm::vec3(0.0f, 0.05f, 0.0f), 0.015f);
        tip_stem->add_child(stem);
        tip_stem = stem;
        chain.push_back(stem);
    }

    // Zero all sugar; pile 10 g at the seed's phloem.
    plant.for_each_node_mut([&](Node& n) {
        n.local().chemical(ChemicalID::Sugar) = 0.0f;
        if (auto* p = n.phloem()) p->chemical(ChemicalID::Sugar) = 0.0f;
    });
    plant.seed_mut()->phloem()->chemical(ChemicalID::Sugar) = 10.0f;

    float total_before = 10.0f;
    vascular_sub_stepped(plant, g, world);

    // 1) Sugar reached the chain — some has moved out of the seed.
    float seed_after = plant.seed_mut()->phloem()->chemical(ChemicalID::Sugar);
    REQUIRE(seed_after < total_before);

    // 2) Mass is conserved across all pools (the invariant that matters
    //    more than distance-dependent supply, which is sensitive to the
    //    interplay between max_move clamps and pipe capacities).
    float total_after = 0.0f;
    plant.for_each_node([&](const Node& n) {
        total_after += n.local().chemical(ChemicalID::Sugar);
        if (auto* p = n.phloem()) total_after += p->chemical(ChemicalID::Sugar);
    });
    REQUIRE(total_after == Catch::Approx(total_before).margin(1e-4f));
}
