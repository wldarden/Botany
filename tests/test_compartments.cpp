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

TEST_CASE("StemNode exposes non-null phloem() and xylem()", "[compartments]") {
    StemNode stem(/* id */ 42, glm::vec3(0), /* radius */ 0.015f);
    REQUIRE(stem.phloem() != nullptr);
    REQUIRE(stem.xylem()  != nullptr);

    // The pools start empty.
    REQUIRE(stem.phloem()->chemical(ChemicalID::Sugar) == 0.0f);
    REQUIRE(stem.xylem()->chemical(ChemicalID::Water)  == 0.0f);
}

TEST_CASE("StemNode phloem and xylem are independent pools", "[compartments]") {
    StemNode stem(1, glm::vec3(0), 0.015f);
    stem.phloem()->chemical(ChemicalID::Sugar) = 5.0f;
    stem.xylem()->chemical(ChemicalID::Water)  = 3.0f;
    REQUIRE(stem.phloem()->chemical(ChemicalID::Sugar) == 5.0f);
    REQUIRE(stem.xylem()->chemical(ChemicalID::Water)  == 3.0f);
    REQUIRE(stem.phloem()->chemical(ChemicalID::Water) == 0.0f);  // water not in phloem
    REQUIRE(stem.xylem()->chemical(ChemicalID::Sugar)  == 0.0f);  // sugar not in xylem
}

TEST_CASE("RootNode exposes non-null phloem() and xylem()", "[compartments]") {
    RootNode root(1, glm::vec3(0, -0.05f, 0), 0.015f);
    REQUIRE(root.phloem() != nullptr);
    REQUIRE(root.xylem()  != nullptr);
}

TEST_CASE("LeafNode/ApicalNode/RootApicalNode return nullptr pools", "[compartments]") {
    LeafNode leaf(1, glm::vec3(0), 0.01f);
    ApicalNode apical(2, glm::vec3(0), 0.01f);
    RootApicalNode root_apical(3, glm::vec3(0, -0.05f, 0), 0.01f);

    REQUIRE(leaf.phloem()        == nullptr);
    REQUIRE(leaf.xylem()         == nullptr);
    REQUIRE(apical.phloem()      == nullptr);
    REQUIRE(apical.xylem()       == nullptr);
    REQUIRE(root_apical.phloem() == nullptr);
    REQUIRE(root_apical.xylem()  == nullptr);
}

TEST_CASE("nearest_phloem_upstream finds parent stem", "[compartments][walkup]") {
    StemNode stem(1, glm::vec3(0), 0.015f);
    LeafNode leaf(2, glm::vec3(0.05f, 0, 0), 0.01f);
    stem.add_child(&leaf);
    leaf.parent = &stem;
    REQUIRE(leaf.nearest_phloem_upstream() == stem.phloem());
}

TEST_CASE("nearest_phloem_upstream walks past non-conduit ancestors", "[compartments][walkup]") {
    StemNode stem(1, glm::vec3(0), 0.015f);
    ApicalNode apical(2, glm::vec3(0.05f, 0, 0), 0.01f);
    LeafNode leaf(3, glm::vec3(0.08f, 0, 0), 0.01f);
    stem.add_child(&apical);
    apical.parent = &stem;
    apical.add_child(&leaf);
    leaf.parent = &apical;
    // Leaf's direct parent is apical (no phloem).  Walk-up must skip
    // past the apical and find the stem.
    REQUIRE(leaf.nearest_phloem_upstream() == stem.phloem());
}

TEST_CASE("nearest_phloem_upstream returns nullptr when no upstream conduit", "[compartments][walkup]") {
    LeafNode orphan(1, glm::vec3(0), 0.01f);
    // No parent set.
    REQUIRE(orphan.nearest_phloem_upstream() == nullptr);
}

TEST_CASE("nearest_xylem_upstream finds parent root", "[compartments][walkup]") {
    RootNode root(1, glm::vec3(0, -0.05f, 0), 0.015f);
    RootApicalNode ra(2, glm::vec3(0, -0.10f, 0), 0.01f);
    root.add_child(&ra);
    ra.parent = &root;
    REQUIRE(ra.nearest_xylem_upstream() == root.xylem());
}
