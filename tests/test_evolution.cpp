#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/engine.h"
#include "engine/genome.h"

using Catch::Matchers::WithinAbs;

TEST_CASE("evolution stub", "[evolution]") {
    REQUIRE(true);
}

TEST_CASE("Plant tracks total sugar produced", "[evolution]") {
    botany::Engine engine;
    botany::Genome g = botany::default_genome();
    engine.create_plant(g, glm::vec3(0.0f));

    // Run enough ticks for leaves to grow and photosynthesize
    for (int i = 0; i < 200; i++) {
        engine.tick();
    }

    float produced = engine.get_plant(0).total_sugar_produced();
    REQUIRE(produced > 0.0f);
}
