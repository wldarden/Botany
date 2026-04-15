#include <catch2/catch_test_macros.hpp>
#include "engine/node/meristems/helpers.h"
#include <cmath>

using namespace botany::meristem_helpers;

TEST_CASE("auxin_growth_factor: zero auxin returns 1.0", "[auxin_sensitivity]") {
    REQUIRE(auxin_growth_factor(0.0f, 0.5f, 0.2f) == 1.0f);
    REQUIRE(auxin_growth_factor(0.0f, -0.3f, 0.1f) == 1.0f);
}

TEST_CASE("auxin_growth_factor: positive max_boost promotes growth", "[auxin_sensitivity]") {
    float factor = auxin_growth_factor(0.2f, 0.5f, 0.2f);
    // At half-saturation, factor should be 1.0 + 0.5 * 0.5 = 1.25
    REQUIRE(factor > 1.0f);
    REQUIRE(factor < 1.5f);
    REQUIRE(std::abs(factor - 1.25f) < 0.01f);
}

TEST_CASE("auxin_growth_factor: negative max_boost inhibits growth", "[auxin_sensitivity]") {
    float factor = auxin_growth_factor(0.1f, -0.3f, 0.1f);
    // At half-saturation, factor should be 1.0 + (-0.3) * 0.5 = 0.85
    REQUIRE(factor < 1.0f);
    REQUIRE(factor > 0.7f);
    REQUIRE(std::abs(factor - 0.85f) < 0.01f);
}

TEST_CASE("auxin_growth_factor: high auxin asymptotes at 1.0 + max_boost", "[auxin_sensitivity]") {
    float factor = auxin_growth_factor(1000.0f, 0.5f, 0.2f);
    REQUIRE(std::abs(factor - 1.5f) < 0.01f);

    float factor_neg = auxin_growth_factor(1000.0f, -0.3f, 0.1f);
    REQUIRE(std::abs(factor_neg - 0.7f) < 0.01f);
}

TEST_CASE("auxin_growth_factor: zero max_boost always returns 1.0", "[auxin_sensitivity]") {
    REQUIRE(auxin_growth_factor(0.5f, 0.0f, 0.2f) == 1.0f);
    REQUIRE(auxin_growth_factor(100.0f, 0.0f, 0.1f) == 1.0f);
}

TEST_CASE("auxin_growth_factor: monotonically increases with auxin for positive boost", "[auxin_sensitivity]") {
    float f1 = auxin_growth_factor(0.05f, 0.5f, 0.2f);
    float f2 = auxin_growth_factor(0.1f, 0.5f, 0.2f);
    float f3 = auxin_growth_factor(0.5f, 0.5f, 0.2f);
    REQUIRE(f1 < f2);
    REQUIRE(f2 < f3);
}

TEST_CASE("auxin_growth_factor: monotonically decreases with auxin for negative boost", "[auxin_sensitivity]") {
    float f1 = auxin_growth_factor(0.05f, -0.3f, 0.1f);
    float f2 = auxin_growth_factor(0.1f, -0.3f, 0.1f);
    float f3 = auxin_growth_factor(0.5f, -0.3f, 0.1f);
    REQUIRE(f1 > f2);
    REQUIRE(f2 > f3);
}
