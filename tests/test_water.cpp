#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <glm/geometric.hpp>
#include "engine/plant.h"
#include "engine/node/root_node.h"
#include "engine/node/tissues/root_apical.h"
#include "engine/node/tissues/leaf.h"
#include "engine/genome.h"
#include "engine/light.h"
#include "engine/world_params.h"
#include "engine/sugar.h"
#include "engine/chemical/chemical.h"

using namespace botany;

static void flush_buffers(Plant& plant) {
    plant.for_each_node_mut([](Node& n) {
        for (auto& [id, amount] : n.transport_received) {
            n.chemical(id) += amount;
        }
        n.transport_received.clear();
    });
}
using Catch::Matchers::WithinAbs;

TEST_CASE("Seed node starts with water", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    REQUIRE(plant.seed()->chemical(ChemicalID::Water) > 0.0f);
}

// === Water cap tests ===

TEST_CASE("water_cap scales with stem volume", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.1f);
    plant.seed_mut()->add_child(stem);

    float volume = 3.14159f * 0.1f * 0.1f * 1.0f;
    float expected = volume * g.water_storage_density_stem;
    REQUIRE_THAT(water_cap(*stem, g), WithinAbs(expected, 1e-4));
}

TEST_CASE("water_cap scales with leaf area", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    plant.seed_mut()->add_child(leaf);

    float area = 0.5f * 0.5f;
    float expected = area * g.water_storage_density_leaf;
    REQUIRE_THAT(water_cap(*leaf, g), WithinAbs(expected, 1e-6));
}

TEST_CASE("water_cap for root scales with volume", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* root = plant.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.8f, 0.0f), 0.025f);
    plant.seed_mut()->add_child(root);

    float length = 0.8f;
    float volume = 3.14159f * 0.025f * 0.025f * length;
    float expected = volume * g.water_storage_density_stem;
    REQUIRE_THAT(water_cap(*root, g), WithinAbs(expected, 1e-5));
}

TEST_CASE("water_cap for meristem is fixed", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* apical = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::APICAL) apical = &n;
    });
    REQUIRE(apical != nullptr);
    REQUIRE_THAT(water_cap(*apical, g), WithinAbs(g.water_cap_meristem, 1e-6));
}

TEST_CASE("water_cap returns minimum for tiny nodes", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* tiny = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.001f, 0.0f), 0.001f);
    plant.seed_mut()->add_child(tiny);

    REQUIRE_THAT(water_cap(*tiny, g), WithinAbs(g.sugar_cap_minimum, 1e-6));
}

// === Water transport tests ===

// Water is handled exclusively by xylem_resolve — transport_with_children() skips it.
// Verifying that invariant here.
TEST_CASE("Water does not move via local transport (xylem-only)", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* child = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.1f, 0.0f), 0.02f);
    plant.seed_mut()->add_child(child);

    plant.seed_mut()->chemical(ChemicalID::Water) = 50.0f;
    child->chemical(ChemicalID::Water) = 0.0f;

    plant.seed_mut()->transport_with_children(g);
    flush_buffers(plant);

    // Water must NOT move via local transport — xylem_resolve handles it.
    REQUIRE(child->chemical(ChemicalID::Water) == 0.0f);
    REQUIRE(plant.seed()->chemical(ChemicalID::Water) == 50.0f);
}

TEST_CASE("Water transport conserves total water", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    for (int i = 0; i < 3; i++) {
        Node* child = plant.create_node(NodeType::STEM,
            glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
        plant.seed_mut()->add_child(child);
    }

    plant.seed_mut()->chemical(ChemicalID::Water) = 20.0f;
    plant.seed_mut()->children[0]->chemical(ChemicalID::Water) = 5.0f;

    float total_before = 0.0f;
    plant.for_each_node([&](const Node& n) { total_before += n.chemical(ChemicalID::Water); });

    for (int i = 0; i < 10; i++) {
        plant.for_each_node_mut([&](Node& n) {
            n.transport_with_children(g);
        });
        flush_buffers(plant);
    }

    float total_after = 0.0f;
    plant.for_each_node([&](const Node& n) { total_after += n.chemical(ChemicalID::Water); });

    REQUIRE_THAT(total_after, WithinAbs(total_before, 1e-4));
}

// === Task 6: Root water absorption tests ===

TEST_CASE("Root nodes absorb water proportional to surface area", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* root = plant.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.5f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(root);
    root->chemical(ChemicalID::Water) = 0.0f;
    root->chemical(ChemicalID::Sugar) = 10.0f;

    WorldParams wp = default_world_params();
    wp.soil_moisture = 1.0f;

    root->tick(plant, wp);

    REQUIRE(root->chemical(ChemicalID::Water) > 0.0f);
}

TEST_CASE("Root absorption scales with soil_moisture", "[water]") {
    Genome g = default_genome();

    Plant plant_wet(g, glm::vec3(0.0f));
    Node* root_wet = plant_wet.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.5f, 0.0f), 0.05f);
    plant_wet.seed_mut()->add_child(root_wet);
    root_wet->chemical(ChemicalID::Water) = 0.0f;
    root_wet->chemical(ChemicalID::Sugar) = 10.0f;

    WorldParams wp_wet = default_world_params();
    wp_wet.soil_moisture = 1.0f;
    root_wet->tick(plant_wet, wp_wet);

    Plant plant_dry(g, glm::vec3(0.0f));
    Node* root_dry = plant_dry.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.5f, 0.0f), 0.05f);
    plant_dry.seed_mut()->add_child(root_dry);
    root_dry->chemical(ChemicalID::Water) = 0.0f;
    root_dry->chemical(ChemicalID::Sugar) = 10.0f;

    WorldParams wp_dry = default_world_params();
    wp_dry.soil_moisture = 0.2f;
    root_dry->tick(plant_dry, wp_dry);

    REQUIRE(root_wet->chemical(ChemicalID::Water) > root_dry->chemical(ChemicalID::Water));
}

TEST_CASE("Zero soil_moisture means zero water absorption", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* root = plant.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.5f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(root);
    root->chemical(ChemicalID::Water) = 0.0f;
    root->chemical(ChemicalID::Sugar) = 10.0f;

    WorldParams wp = default_world_params();
    wp.soil_moisture = 0.0f;

    root->tick(plant, wp);

    REQUIRE(root->chemical(ChemicalID::Water) < 1e-6f);
}

TEST_CASE("Root apical tips also absorb water", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* root_apical = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::ROOT_APICAL) root_apical = &n;
    });
    REQUIRE(root_apical != nullptr);

    root_apical->chemical(ChemicalID::Water) = 0.0f;
    root_apical->chemical(ChemicalID::Sugar) = 10.0f;

    WorldParams wp = default_world_params();
    wp.soil_moisture = 1.0f;

    root_apical->tick(plant, wp);

    REQUIRE(root_apical->chemical(ChemicalID::Water) > 0.0f);
}

// === Task 7: Leaf transpiration and photosynthesis water cost tests ===

TEST_CASE("Leaves lose water via transpiration", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    plant.seed_mut()->add_child(leaf);

    float initial_water = 5.0f;
    leaf->chemical(ChemicalID::Water) = initial_water;
    leaf->chemical(ChemicalID::Sugar) = 10.0f;

    WorldParams wp = default_world_params();
    wp.light_level = 2.0f;
    compute_light_exposure(plant, wp);
    leaf->tick(plant, wp);

    REQUIRE(leaf->chemical(ChemicalID::Water) < initial_water);
}

TEST_CASE("Transpiration scales with light exposure", "[water]") {
    Genome g = default_genome();

    // light_exposure (set by shadow computation) drives transpiration, not world.light_level.
    // Set it directly on each leaf to isolate the transpiration scaling.

    Plant plant_bright(g, glm::vec3(0.0f));
    Node* leaf_bright = plant_bright.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf_bright->as_leaf()->leaf_size = 0.5f;
    leaf_bright->as_leaf()->light_exposure = 1.0f;
    plant_bright.seed_mut()->add_child(leaf_bright);
    leaf_bright->chemical(ChemicalID::Water) = 5.0f;
    leaf_bright->chemical(ChemicalID::Sugar) = 10.0f;

    WorldParams wp = default_world_params();
    leaf_bright->tick(plant_bright, wp);

    Plant plant_dim(g, glm::vec3(0.0f));
    Node* leaf_dim = plant_dim.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf_dim->as_leaf()->leaf_size = 0.5f;
    leaf_dim->as_leaf()->light_exposure = 0.1f;
    plant_dim.seed_mut()->add_child(leaf_dim);
    leaf_dim->chemical(ChemicalID::Water) = 5.0f;
    leaf_dim->chemical(ChemicalID::Sugar) = 10.0f;

    leaf_dim->tick(plant_dim, wp);

    REQUIRE(leaf_bright->chemical(ChemicalID::Water) < leaf_dim->chemical(ChemicalID::Water));
}

TEST_CASE("Water does not go below zero from transpiration", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 1.0f;
    plant.seed_mut()->add_child(leaf);

    leaf->chemical(ChemicalID::Water) = 0.0001f;
    leaf->chemical(ChemicalID::Sugar) = 10.0f;

    WorldParams wp = default_world_params();
    wp.light_level = 10.0f;
    compute_light_exposure(plant, wp);
    leaf->tick(plant, wp);

    REQUIRE(leaf->chemical(ChemicalID::Water) >= 0.0f);
}

TEST_CASE("Photosynthesis consumes water proportional to sugar produced", "[water]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    plant.seed_mut()->add_child(leaf);

    leaf->chemical(ChemicalID::Water) = 100.0f;
    leaf->chemical(ChemicalID::Sugar) = 0.0f;

    WorldParams wp = default_world_params();
    wp.light_level = 2.0f;
    compute_light_exposure(plant, wp);

    float water_before = leaf->chemical(ChemicalID::Water);
    leaf->tick(plant, wp);
    float water_after = leaf->chemical(ChemicalID::Water);

    float leaf_area = 0.5f * 0.5f;
    float transpiration_only = g.transpiration_rate * leaf_area * leaf->as_leaf()->light_exposure;
    float total_loss = water_before - water_after;

    REQUIRE(total_loss > transpiration_only - 1e-6f);
}

// === Gradient-based self-limiting absorption tests ===

TEST_CASE("Root absorption self-limits as root fills up", "[water]") {
    // Gradient formula: absorbed ∝ (soil_moisture - fill_fraction)
    // A root at 90% capacity should absorb ~9× less than one at 10% capacity.
    // Ticking only the root node (no seed tick) means water_delta == absorbed from soil.
    Genome g = default_genome();
    WorldParams wp = default_world_params();
    wp.soil_moisture = 1.0f;

    // Root at 10% fill
    Plant plant_low(g, glm::vec3(0.0f));
    Node* root_low = plant_low.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.5f, 0.0f), 0.05f);
    plant_low.seed_mut()->add_child(root_low);
    float cap = water_cap(*root_low, g);
    root_low->chemical(ChemicalID::Water) = 0.1f * cap;
    root_low->chemical(ChemicalID::Sugar) = 10.0f;

    float before_low = root_low->chemical(ChemicalID::Water);
    root_low->tick(plant_low, wp);
    float absorbed_low = root_low->chemical(ChemicalID::Water) - before_low;

    // Root at 90% fill — same geometry, just more water inside
    Plant plant_high(g, glm::vec3(0.0f));
    Node* root_high = plant_high.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.5f, 0.0f), 0.05f);
    plant_high.seed_mut()->add_child(root_high);
    root_high->chemical(ChemicalID::Water) = 0.9f * cap;
    root_high->chemical(ChemicalID::Sugar) = 10.0f;

    float before_high = root_high->chemical(ChemicalID::Water);
    root_high->tick(plant_high, wp);
    float absorbed_high = root_high->chemical(ChemicalID::Water) - before_high;

    // 10%-full root has gradient 0.9, 90%-full has gradient 0.1 → ratio ~9×
    REQUIRE(absorbed_low > 0.0f);                       // 10%-full root absorbs
    REQUIRE(absorbed_high < absorbed_low * 0.5f);       // 90%-full absorbs much less
}

// === Seed initialization regression test ===

TEST_CASE("Seed water initializes within capacity (not to seed_sugar)", "[water]") {
    // Regression: seed was incorrectly initialized with genome.seed_sugar (48 ml)
    // which is far over the ~2 ml water capacity of the two meristem children.
    // Overcap blocks outward transport and masks the root-to-leaf gradient.
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    float seed_water = plant.seed()->chemical(ChemicalID::Water);
    float seed_cap   = water_cap(*plant.seed(), g);

    // Water must be at or below capacity (no overcap)
    REQUIRE(seed_water <= seed_cap + 1e-4f);

    // Capacity should be at least 2 × water_cap_meristem (shoot + root meristem children)
    REQUIRE(seed_cap >= 2.0f * g.water_cap_meristem - 1e-4f);

    // seed_sugar (48 g-glucose) must not be used as the water amount
    REQUIRE(seed_water < g.seed_sugar - 1.0f);
}
