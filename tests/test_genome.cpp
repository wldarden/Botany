// tests/test_genome.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/genome.h"
#include "engine/world_params.h"

using namespace botany;

TEST_CASE("default_genome returns valid parameters", "[genome]") {
    Genome g = default_genome();

    SECTION("hormone rates are positive") {
        REQUIRE(g.apical_auxin_baseline > 0.0f);
        REQUIRE(g.auxin_diffusion_rate > 0.0f);
        REQUIRE(g.auxin_decay_rate > 0.0f);
        REQUIRE(g.cytokinin_production_rate > 0.0f);
        REQUIRE(g.cytokinin_diffusion_rate > 0.0f);
        REQUIRE(g.cytokinin_decay_rate > 0.0f);
    }

    SECTION("transport rates are fractions (0, 1]") {
        REQUIRE(g.auxin_diffusion_rate <= 1.0f);
        REQUIRE(g.cytokinin_diffusion_rate <= 1.0f);
    }

    SECTION("decay rates are fractions (0, 1)") {
        REQUIRE(g.auxin_decay_rate < 1.0f);
        REQUIRE(g.cytokinin_decay_rate < 1.0f);
    }

    SECTION("growth rates are positive") {
        REQUIRE(g.growth_rate > 0.0f);
        REQUIRE(g.root_growth_rate > 0.0f);
    }

    SECTION("internode lengths are positive") {
        REQUIRE(g.max_internode_length > 0.0f);
    }

    SECTION("sugar parameters are positive") {
        REQUIRE(g.sugar_diffusion_rate >= 0.0f);  // 0 is valid: phloem fully owns sugar transport
        WorldParams w = default_world_params();
        REQUIRE(w.sugar_production_rate > 0.0f);
        REQUIRE(w.sugar_maintenance_leaf > 0.0f);
        REQUIRE(w.sugar_maintenance_stem > 0.0f);
        REQUIRE(w.sugar_maintenance_root > 0.0f);
        REQUIRE(w.sugar_maintenance_meristem > 0.0f);
        REQUIRE(w.sugar_cost_activation > 0.0f);
        REQUIRE(w.sugar_cost_meristem_growth > 0.0f);
        REQUIRE(w.sugar_cost_root_growth > 0.0f);
        REQUIRE(w.sugar_cost_stem_growth > 0.0f);
        REQUIRE(g.cytokinin_growth_threshold > 0.0f);
    }
}
