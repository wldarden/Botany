#include <catch2/catch_test_macros.hpp>
#include "engine/compartments.h"
#include "engine/node/node.h"
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"

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

TEST_CASE("Base Node returns nullptr from phloem() and xylem()", "[compartments]") {
    // We can't instantiate Node directly (abstract in spirit), but we can
    // verify the virtual default by checking a LeafNode — which inherits
    // the default nullptr behavior in this task (overrides arrive in Task 5/6).
    LeafNode leaf(/* id */ 1, glm::vec3(0), /* radius */ 0.01f);
    REQUIRE(leaf.phloem() == nullptr);
    REQUIRE(leaf.xylem() == nullptr);
}
