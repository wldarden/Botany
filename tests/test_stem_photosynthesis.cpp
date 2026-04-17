// tests/test_stem_photosynthesis.cpp — corticular (green stem) photosynthesis tests.
// Biology: young stems with chlorenchyma under the epidermis photosynthesize until
// radial thickening covers them with bark. This is corticular photosynthesis —
// real, significant (10-15% of branch carbon in young tissue).
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/chemical/chemical.h"
#include <glm/geometric.hpp>

using namespace botany;
using Catch::Matchers::WithinAbs;

TEST_CASE("Young stem produces sugar when exposed to light", "[stem_photosynthesis]") {
    Genome g = default_genome();
    g.stem_photosynthesis_rate = 0.01f;          // elevated for clear signal
    g.stem_green_radius_threshold = 0.05f;        // threshold above initial_radius (0.015)
    g.cambium_responsiveness = 0.0f;              // disable thickening

    Plant plant(g, glm::vec3(0.0f));
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.02f);
    plant.seed_mut()->add_child(stem);
    stem->chemical(ChemicalID::Sugar) = 100.0f;  // plenty for maintenance

    WorldParams wp = default_world_params();
    wp.light_level = 1.0f;

    float produced_before = plant.total_sugar_produced();
    stem->tick(plant, wp);
    float produced_after = plant.total_sugar_produced();

    REQUIRE(produced_after > produced_before);
}

TEST_CASE("Mature stem (radius >= threshold) produces no sugar", "[stem_photosynthesis]") {
    Genome g = default_genome();
    g.stem_photosynthesis_rate = 0.01f;
    g.stem_green_radius_threshold = 0.03f;        // threshold below our stem radius

    Plant plant(g, glm::vec3(0.0f));
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.1f);  // radius=0.1 >> threshold
    plant.seed_mut()->add_child(stem);
    stem->chemical(ChemicalID::Sugar) = 100.0f;

    WorldParams wp = default_world_params();
    wp.light_level = 1.0f;

    float produced_before = plant.total_sugar_produced();
    stem->tick(plant, wp);
    float produced_after = plant.total_sugar_produced();

    REQUIRE_THAT(produced_after, WithinAbs(produced_before, 1e-8f));
}

TEST_CASE("Zero light stops green stem photosynthesis", "[stem_photosynthesis]") {
    Genome g = default_genome();
    g.stem_photosynthesis_rate = 0.01f;
    g.stem_green_radius_threshold = 0.05f;
    g.cambium_responsiveness = 0.0f;

    Plant plant(g, glm::vec3(0.0f));
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.02f);
    plant.seed_mut()->add_child(stem);
    stem->chemical(ChemicalID::Sugar) = 100.0f;

    WorldParams wp = default_world_params();
    wp.light_level = 0.0f;

    float produced_before = plant.total_sugar_produced();
    stem->tick(plant, wp);
    float produced_after = plant.total_sugar_produced();

    REQUIRE_THAT(produced_after, WithinAbs(produced_before, 1e-8f));
}

TEST_CASE("Stem photosynthesis proportional to surface area (longer = more)", "[stem_photosynthesis]") {
    Genome g = default_genome();
    g.stem_photosynthesis_rate = 0.01f;
    g.stem_green_radius_threshold = 0.05f;
    g.cambium_responsiveness = 0.0f;

    // Short stem: length ~0.3 dm
    Plant plant_short(g, glm::vec3(0.0f));
    Node* short_stem = plant_short.create_node(NodeType::STEM, glm::vec3(0.0f, 0.3f, 0.0f), 0.02f);
    plant_short.seed_mut()->add_child(short_stem);
    short_stem->chemical(ChemicalID::Sugar) = 100.0f;

    // Long stem: length ~1.0 dm
    Plant plant_long(g, glm::vec3(0.0f));
    Node* long_stem = plant_long.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.02f);
    plant_long.seed_mut()->add_child(long_stem);
    long_stem->chemical(ChemicalID::Sugar) = 100.0f;

    WorldParams wp = default_world_params();
    wp.light_level = 1.0f;

    float before_short = plant_short.total_sugar_produced();
    short_stem->tick(plant_short, wp);
    float delta_short = plant_short.total_sugar_produced() - before_short;

    float before_long = plant_long.total_sugar_produced();
    long_stem->tick(plant_long, wp);
    float delta_long = plant_long.total_sugar_produced() - before_long;

    REQUIRE(delta_long > delta_short);  // more surface area = more production
}

TEST_CASE("Stem photosynthesis credited correctly to plant sugar_produced", "[stem_photosynthesis]") {
    Genome g = default_genome();
    g.stem_photosynthesis_rate = 0.01f;
    g.stem_green_radius_threshold = 0.05f;
    g.cambium_responsiveness = 0.0f;

    Plant plant(g, glm::vec3(0.0f));
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.02f);
    plant.seed_mut()->add_child(stem);
    stem->chemical(ChemicalID::Sugar) = 100.0f;
    // Freeze elongation so length stays exactly at 0.5
    stem->age = g.internode_maturation_ticks;

    WorldParams wp = default_world_params();
    wp.light_level = 1.0f;

    float produced_before = plant.total_sugar_produced();
    stem->tick(plant, wp);
    float delta = plant.total_sugar_produced() - produced_before;

    // Expected: 2 * PI * r * L * light_level * rate
    float length = glm::length(stem->offset);
    float expected = 2.0f * 3.14159f * 0.02f * length * 1.0f * 0.01f;
    REQUIRE_THAT(delta, WithinAbs(expected, 1e-5f));
}

TEST_CASE("Stem photosynthesis rate zero disables production", "[stem_photosynthesis]") {
    Genome g = default_genome();
    g.stem_photosynthesis_rate = 0.0f;           // explicitly disabled
    g.stem_green_radius_threshold = 0.05f;

    Plant plant(g, glm::vec3(0.0f));
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.5f, 0.0f), 0.02f);
    plant.seed_mut()->add_child(stem);
    stem->chemical(ChemicalID::Sugar) = 100.0f;

    WorldParams wp = default_world_params();
    wp.light_level = 1.0f;

    float produced_before = plant.total_sugar_produced();
    stem->tick(plant, wp);
    float produced_after = plant.total_sugar_produced();

    REQUIRE_THAT(produced_after, WithinAbs(produced_before, 1e-8f));
}
