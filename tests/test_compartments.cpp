#include <catch2/catch_test_macros.hpp>
#include "engine/compartments.h"

using namespace botany;

TEST_CASE("LocalEnv stores and reads chemicals by ChemicalID", "[compartments]") {
    LocalEnv env;
    env.chemical(ChemicalID::Sugar) = 1.5f;
    REQUIRE(env.chemical(ChemicalID::Sugar) == 1.5f);
    REQUIRE(env.chemical(ChemicalID::Water) == 0.0f);
}

TEST_CASE("TransportPool stores and reads chemicals by ChemicalID", "[compartments]") {
    TransportPool pool;
    pool.chemical(ChemicalID::Sugar) = 2.0f;
    REQUIRE(pool.chemical(ChemicalID::Sugar) == 2.0f);
    REQUIRE(pool.chemical(ChemicalID::Water) == 0.0f);
}

TEST_CASE("LocalEnv const accessor returns 0 for missing chemicals", "[compartments]") {
    const LocalEnv env;
    REQUIRE(env.chemical(ChemicalID::Sugar) == 0.0f);
}
