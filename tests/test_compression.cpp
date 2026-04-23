#include <catch2/catch_test_macros.hpp>
#include <glm/vec3.hpp>

#include "engine/compression.h"
#include "engine/plant.h"
#include "engine/genome.h"

using namespace botany;

TEST_CASE("compress_plant on a fresh young plant reports no merges", "[compression][stub]") {
    Plant plant(default_genome(), glm::vec3(0.0f));
    CompressionParams params;
    CompressionResult r = compress_plant(plant, params);
    REQUIRE(r.merges_performed == 0u);
    REQUIRE(r.passes_run == 0u);
}
