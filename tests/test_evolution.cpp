#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/engine.h"
#include "engine/genome.h"
#include "evolution/genome_bridge.h"

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

TEST_CASE("Genome round-trips through StructuredGenome", "[evolution]") {
    botany::Genome original = botany::default_genome();
    auto sg = botany::to_structured(original);
    botany::Genome restored = botany::from_structured(sg);

    // Spot-check representative fields from different groups
    REQUIRE_THAT(restored.auxin_production_rate,
                 WithinAbs(original.auxin_production_rate, 1e-6));
    REQUIRE_THAT(restored.branch_angle,
                 WithinAbs(original.branch_angle, 1e-6));
    REQUIRE_THAT(restored.sugar_production_rate,
                 WithinAbs(original.sugar_production_rate, 1e-6));
    REQUIRE_THAT(restored.ga_elongation_sensitivity,
                 WithinAbs(original.ga_elongation_sensitivity, 1e-6));
    REQUIRE_THAT(restored.ethylene_abscission_threshold,
                 WithinAbs(original.ethylene_abscission_threshold, 1e-6));
    REQUIRE(restored.internode_maturation_ticks == original.internode_maturation_ticks);
    REQUIRE(restored.senescence_duration == original.senescence_duration);
}

TEST_CASE("Genome template has linkage groups", "[evolution]") {
    auto tmpl = botany::build_genome_template(botany::default_genome());
    auto& groups = tmpl.linkage_groups();
    REQUIRE(groups.size() == 8);

    // Verify auxin group has 5 genes
    bool found_auxin = false;
    for (auto& g : groups) {
        if (g.name == "auxin") {
            found_auxin = true;
            REQUIRE(g.gene_tags.size() == 5);
        }
    }
    REQUIRE(found_auxin);
}
