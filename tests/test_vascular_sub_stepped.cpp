#include <catch2/catch_test_macros.hpp>
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
