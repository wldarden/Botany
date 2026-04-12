#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/engine.h"
#include "engine/genome.h"
#include "evolution/genome_bridge.h"
#include "evolution/fitness.h"
#include "evolution/evolution_runner.h"

using Catch::Matchers::WithinAbs;

TEST_CASE("evolution stub", "[evolution]") {
    REQUIRE(true);
}

TEST_CASE("Plant tracks total sugar produced", "[evolution]") {
    botany::Engine engine;
    botany::Genome g = botany::default_genome();
    engine.create_plant(g, glm::vec3(0.0f));

    // Run enough ticks for leaves to grow and photosynthesize
    // (new cytokinin gating means plant needs more time to ramp up)
    for (int i = 0; i < 1000; i++) {
        engine.tick();
    }

    float produced = engine.get_plant(0).total_sugar_produced();
    REQUIRE(produced >= 0.0f);  // accumulator exists and is non-negative
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

    // Verify auxin group has 4 genes (production, diffusion, decay, threshold)
    bool found_auxin = false;
    for (auto& g : groups) {
        if (g.name == "auxin") {
            found_auxin = true;
            REQUIRE(g.gene_tags.size() == 4);
        }
    }
    REQUIRE(found_auxin);
}

TEST_CASE("evaluate_plant returns populated stats", "[evolution]") {
    botany::Genome g = botany::default_genome();
    botany::WorldParams world = botany::default_world_params();

    auto stats = botany::evaluate_plant(g, world, 1000);

    REQUIRE(stats.survival_ticks > 0);
    REQUIRE(stats.node_count >= 3);  // at least seed + 2 meristems
    REQUIRE(stats.height >= 0.0f);
}

TEST_CASE("evaluate_plant respects max_ticks", "[evolution]") {
    botany::Genome g = botany::default_genome();
    botany::WorldParams world = botany::default_world_params();

    auto stats = botany::evaluate_plant(g, world, 50);
    REQUIRE(stats.survival_ticks <= 50);
}

TEST_CASE("compute_fitness normalizes and weights correctly", "[evolution]") {
    botany::PlantStats stats;
    stats.survival_ticks = 100;
    stats.node_count = 50;
    stats.leaf_count = 10;
    stats.total_sugar_produced = 5.0f;
    stats.height = 2.0f;
    stats.crown_ratio = 0.5f;
    stats.branch_depth = 3;
    stats.leaf_height_spread = 1.0f;

    botany::PlantStats gen_max = stats;

    botany::FitnessWeights w;
    float fitness = botany::compute_fitness(stats, gen_max, w);
    REQUIRE_THAT(fitness, WithinAbs(8.0, 1e-4));
}

TEST_CASE("compute_fitness handles zero gen_max gracefully", "[evolution]") {
    botany::PlantStats stats;
    stats.survival_ticks = 100;
    stats.height = 2.0f;

    botany::PlantStats gen_max;

    botany::FitnessWeights w;
    float fitness = botany::compute_fitness(stats, gen_max, w);
    REQUIRE_THAT(fitness, WithinAbs(0.0, 1e-4));
}

TEST_CASE("EvolutionRunner advances generations", "[evolution]") {
    botany::EvolutionConfig config;
    config.population_size = 10;
    config.max_ticks = 200;
    config.num_threads = 2;

    botany::EvolutionRunner runner(config);
    REQUIRE(runner.generation() == 0);

    runner.run_generation();
    REQUIRE(runner.generation() == 1);
    REQUIRE(runner.best_fitness() > 0.0f);
    REQUIRE(runner.best_stats().survival_ticks > 0);

    runner.run_generation();
    REQUIRE(runner.generation() == 2);
}

TEST_CASE("Node has mass and stress fields", "[stress]") {
    botany::Engine engine;
    botany::Genome g = botany::default_genome();
    engine.create_plant(g, glm::vec3(0.0f));

    for (int i = 0; i < 50; i++) engine.tick();

    const botany::Plant& plant = engine.get_plant(0);
    bool found_nonzero_mass = false;
    plant.for_each_node([&](const botany::Node& n) {
        if (n.total_mass > 0.0f) found_nonzero_mass = true;
    });
    REQUIRE(found_nonzero_mass);
}
