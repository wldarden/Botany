#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "engine/vascular_sub_stepped.h"
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/world_params.h"

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
