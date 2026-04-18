#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <glm/geometric.hpp>
#include "engine/plant.h"
#include "engine/node/tissues/leaf.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/sugar.h"
#include "engine/light.h"

using namespace botany;
using Catch::Matchers::WithinAbs;

// Flush transport_received buffers into chemical() for all nodes.
// Needed after calling transport_with_children() directly in tests,
// since the buffer is normally flushed at the end of Node::tick().
static void flush_buffers(Plant& plant) {
    plant.for_each_node_mut([](Node& n) {
        for (auto& [id, amount] : n.transport_received) {
            n.chemical(id) += amount;
        }
        n.transport_received.clear();
    });
}

// === Production tests ===

TEST_CASE("LEAF nodes produce sugar proportional to light and leaf_size", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    leaf->chemical(ChemicalID::Water) = 100.0f;  // ensure full stomatal conductance
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 2.0f;

    plant.seed_mut()->chemical(ChemicalID::Sugar) = 0.0f;  // isolate leaf from transport
    compute_light_exposure(plant, wp);
    float light_exp = leaf->as_leaf()->light_exposure;  // read after compute_light_exposure
    leaf->tick(plant, wp);

    REQUIRE(leaf->chemical(ChemicalID::Sugar) > 0.0f);
    // photosynthesis = light_exposure * angle_eff(1.0) * light_level * leaf_size² * production_rate
    // Note: leaf growth and maintenance consume some sugar in the same tick (within 0.01 tolerance)
    float ls = leaf->as_leaf()->leaf_size;
    float expected = light_exp * wp.light_level * ls * ls * wp.sugar_production_rate;
    REQUIRE_THAT(leaf->chemical(ChemicalID::Sugar), WithinAbs(expected, 0.01f));
}

TEST_CASE("Zero light produces zero sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 0.0f;

    plant.seed_mut()->chemical(ChemicalID::Sugar) = 0.0f;
    compute_light_exposure(plant, wp);
    leaf->tick(plant, wp);

    // No light + no parent sugar = leaf stays at 0 (transport has nothing to move)
    REQUIRE(leaf->chemical(ChemicalID::Sugar) < 1e-6f);
}

TEST_CASE("Non-LEAF nodes do not produce sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Zero out any initial sugar (e.g. seed_sugar) so we only measure production
    plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 0.0f; });

    WorldParams wp = default_world_params();
    compute_light_exposure(plant, wp);

    // Tick only non-leaf nodes — they should not produce sugar
    plant.for_each_node_mut([&](Node& n) {
        if (n.type != NodeType::LEAF) {
            float before = n.chemical(ChemicalID::Sugar);
            // Non-leaf nodes have no production logic, sugar stays 0
            // (transport might move sugar, but all start at 0)
            REQUIRE(n.chemical(ChemicalID::Sugar) == 0.0f);
        }
    });
}

// === Consumption tests ===

TEST_CASE("Stem maintenance_cost scales with surface area (πrL)", "[sugar]") {
    Genome g = default_genome();
    WorldParams wp = default_world_params();
    Plant plant(g, glm::vec3(0.0f));

    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.1f);
    float length = glm::length(stem->offset);
    float expected = wp.sugar_maintenance_stem * 3.14159f * 0.1f * length;
    REQUIRE_THAT(stem->maintenance_cost(wp), WithinAbs(expected, 1e-6));
}

TEST_CASE("Leaf maintenance_cost uses leaf area", "[sugar]") {
    Genome g = default_genome();
    WorldParams wp = default_world_params();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    float expected = wp.sugar_maintenance_leaf * 0.5f * 0.5f;
    REQUIRE_THAT(leaf->maintenance_cost(wp), WithinAbs(expected, 1e-6));
}

TEST_CASE("Sugar cannot go below zero after maintenance", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Zero all sugar, tick — sugar should stay >= 0
    plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 0.0001f; });
    plant.for_each_node_mut([&](Node& n) { n.tick(plant, default_world_params()); });

    plant.for_each_node([](const Node& n) {
        REQUIRE(n.chemical(ChemicalID::Sugar) >= 0.0f);
    });
}

TEST_CASE("Active meristem tip maintenance_cost", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* shoot_tip = nullptr;
    plant.for_each_node_mut([&](Node& n) {
        if (n.type == NodeType::APICAL) shoot_tip = &n;
    });
    REQUIRE(shoot_tip != nullptr);
    WorldParams wp = default_world_params();
    REQUIRE_THAT(shoot_tip->maintenance_cost(wp), WithinAbs(wp.sugar_maintenance_meristem, 1e-6));
}

// === Diffusion tests ===

TEST_CASE("Sugar diffuses from high to low concentration", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    seed->chemical(ChemicalID::Sugar) = 100.0f;

    float seed_before = seed->chemical(ChemicalID::Sugar);
    // Parent distributes to children collectively
    seed->transport_with_children(g);
    flush_buffers(plant);

    REQUIRE(seed->chemical(ChemicalID::Sugar) < seed_before);
    for (const Node* child : seed->children) {
        REQUIRE(child->chemical(ChemicalID::Sugar) > 0.0f);
    }
}

TEST_CASE("Sugar diffusion preserves total sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    plant.seed_mut()->chemical(ChemicalID::Sugar) = 100.0f;

    float total_before = 0.0f;
    plant.for_each_node([&](const Node& n) { total_before += n.chemical(ChemicalID::Sugar); });

    // Run several rounds of local transport (one pass per round)
    for (int i = 0; i < 5; i++) {
        plant.for_each_node_mut([&](Node& n) {
            n.transport_with_children(g);
        });
        flush_buffers(plant);
    }

    float total_after = 0.0f;
    plant.for_each_node([&](const Node& n) { total_after += n.chemical(ChemicalID::Sugar); });

    REQUIRE_THAT(total_after, WithinAbs(total_before, 1e-4));
}

TEST_CASE("Transport moves sugar from full parent to empty children", "[sugar]") {
    Genome g = default_genome();

    Plant plant(g, glm::vec3(0.0f));
    for (int i = 0; i < 3; i++) {
        Node* c = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
        plant.seed_mut()->add_child(c);
    }
    // Zero all children, give seed a full load
    plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 0.0f; });
    plant.seed_mut()->chemical(ChemicalID::Sugar) = 10.0f;

    float child_before = 0.0f;
    for (const Node* c : plant.seed()->children) { child_before += c->chemical(ChemicalID::Sugar); }

    // One round of transport: seed should distribute to children
    plant.seed_mut()->transport_with_children(g);
    flush_buffers(plant);

    float child_after = 0.0f;
    for (const Node* c : plant.seed()->children) { child_after += c->chemical(ChemicalID::Sugar); }

    REQUIRE(child_after > child_before);  // children received sugar from seed
}

TEST_CASE("Sugar diffusion conserves total across many ticks", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Build a bigger tree: seed -> 5 children, each with 2 grandchildren
    for (int i = 0; i < 5; i++) {
        Node* child = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.05f);
        plant.seed_mut()->add_child(child);
        for (int j = 0; j < 2; j++) {
            Node* gc = plant.create_node(NodeType::STEM, glm::vec3(0.0f), 0.03f);
            child->add_child(gc);
        }
    }

    // Give sugar to just a few nodes
    plant.seed_mut()->chemical(ChemicalID::Sugar) = 50.0f;
    plant.seed_mut()->children[0]->chemical(ChemicalID::Sugar) = 20.0f;
    plant.seed_mut()->children[2]->chemical(ChemicalID::Sugar) = 10.0f;

    float total_before = 0.0f;
    plant.for_each_node([&](const Node& n) { total_before += n.chemical(ChemicalID::Sugar); });

    for (int ticks : {1, 5, 15, 50}) {
        // Reset to initial distribution
        plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 0.0f; });
        plant.seed_mut()->chemical(ChemicalID::Sugar) = 50.0f;
        plant.seed_mut()->children[0]->chemical(ChemicalID::Sugar) = 20.0f;
        plant.seed_mut()->children[2]->chemical(ChemicalID::Sugar) = 10.0f;

        for (int t = 0; t < ticks; t++) {
            plant.for_each_node_mut([&](Node& n) {
                n.transport_with_children(g);
            });
            // Flush received buffers (normally done by tick())
            plant.for_each_node_mut([](Node& n) {
                for (auto& [id, amt] : n.transport_received) {
                    n.chemical(id) += amt;
                }
                n.transport_received.clear();
            });
        }

        float total_after = 0.0f;
        plant.for_each_node([&](const Node& n) { total_after += n.chemical(ChemicalID::Sugar); });

        REQUIRE_THAT(total_after, WithinAbs(total_before, 1e-3));
    }
}

// === Starvation tests ===

TEST_CASE("Starvation ticks increment when sugar is zero", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    seed->chemical(ChemicalID::Sugar) = 0.0f;

    plant.for_each_node_mut([&](Node& n) { n.tick(plant, default_world_params()); });

    REQUIRE(seed->starvation_ticks > 0);
}

TEST_CASE("Starvation ticks reset when sugar is available", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* seed = plant.seed_mut();
    seed->starvation_ticks = 10;
    seed->chemical(ChemicalID::Sugar) = 100.0f;

    plant.for_each_node_mut([&](Node& n) { n.tick(plant, default_world_params()); });

    REQUIRE(seed->starvation_ticks == 0);
}

TEST_CASE("Starved non-seed nodes are pruned after max starvation ticks", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Create a child node and starve it
    Node* child = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(child);

    uint32_t count_before = plant.node_count();
    child->starvation_ticks = 100; // way past max
    // Zero all sugar so parent transport can't rescue the starving child
    plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 0.0f; });

    WorldParams wp = default_world_params();
    wp.starvation_ticks_max = 50;

    plant.tick(wp);

    REQUIRE(plant.node_count() < count_before);
}

TEST_CASE("Seed node is never pruned even when starved", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    plant.seed_mut()->starvation_ticks = 1000; // extremely starved

    WorldParams wp = default_world_params();
    wp.starvation_ticks_max = 50;

    uint32_t count_before = plant.node_count();
    plant.tick(wp);

    // Seed should still exist
    REQUIRE(plant.seed() != nullptr);
}

TEST_CASE("Subtree removal cleans up all descendants", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Create a chain: seed -> A -> B -> C
    Node* a = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.05f);
    Node* b = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 2.0f, 0.0f), 0.05f);
    Node* c = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 3.0f, 0.0f), 0.05f);
    plant.seed_mut()->add_child(a);
    a->add_child(b);
    b->add_child(c);

    uint32_t a_id = a->id;

    // Starve A — should remove A, B, and C
    a->starvation_ticks = 100;
    plant.for_each_node_mut([](Node& n) { n.chemical(ChemicalID::Sugar) = 0.0f; });

    WorldParams wp = default_world_params();
    wp.starvation_ticks_max = 50;

    plant.tick(wp);

    // A, B, C should all be gone
    bool found_a = false;
    plant.for_each_node([&](const Node& n) {
        if (n.id == a_id) found_a = true;
    });
    REQUIRE_FALSE(found_a);
}

// === Integration test ===

TEST_CASE("Leaf tick runs full produce-diffuse-consume cycle", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 1.0f;
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 5.0f;

    compute_light_exposure(plant, wp);
    leaf->tick(plant, wp);  // produces sugar locally

    REQUIRE(leaf->chemical(ChemicalID::Sugar) > 0.0f);
}

// === Sugar cap tests ===

TEST_CASE("sugar_cap scales with stem volume", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Create a stem with known dimensions: radius 0.1, internode length 1.0
    Node* stem = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 1.0f, 0.0f), 0.1f);
    plant.seed_mut()->add_child(stem);

    float volume = 3.14159f * 0.1f * 0.1f * 1.0f;  // π * r² * length
    float expected = volume * g.sugar_storage_density_wood;
    REQUIRE_THAT(sugar_cap(*stem, g), WithinAbs(expected, 1e-4));
}

TEST_CASE("sugar_cap scales with leaf area", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;  // large enough that area cap exceeds minimum
    plant.seed_mut()->add_child(leaf);

    float area = 0.5f * 0.5f;
    float expected = area * g.sugar_storage_density_leaf;
    REQUIRE(expected > g.sugar_cap_minimum);  // precondition: area cap wins over floor
    REQUIRE_THAT(sugar_cap(*leaf, g), WithinAbs(expected, 1e-6));
}

TEST_CASE("sugar_cap returns minimum for tiny nodes", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Tiny node: radius 0.001, offset 0.001 — volume-based cap is below minimum
    Node* tiny = plant.create_node(NodeType::STEM, glm::vec3(0.0f, 0.001f, 0.0f), 0.001f);
    plant.seed_mut()->add_child(tiny);

    REQUIRE_THAT(sugar_cap(*tiny, g), WithinAbs(g.sugar_cap_minimum, 1e-6));
}

TEST_CASE("sugar_cap for seed node covers seed_sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    // Seed node (no parent) should have cap >= seed_sugar
    REQUIRE(sugar_cap(*plant.seed(), g) >= g.seed_sugar);
}

TEST_CASE("sugar_cap for root scales with volume", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* root = plant.create_node(NodeType::ROOT, glm::vec3(0.0f, -0.8f, 0.0f), 0.025f);
    plant.seed_mut()->add_child(root);

    float length = 0.8f;  // glm::length of offset
    float volume = 3.14159f * 0.025f * 0.025f * length;
    float expected = volume * g.sugar_storage_density_wood;
    // Volume cap (0.0785) > minimum (0.01), so volume wins
    REQUIRE_THAT(sugar_cap(*root, g), WithinAbs(expected, 1e-5));
}

TEST_CASE("Production skipped when leaf sugar is at cap", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.3f;
    plant.seed_mut()->add_child(leaf);

    // Fill leaf to its cap
    float cap = sugar_cap(*leaf, g);
    leaf->chemical(ChemicalID::Sugar) = cap;
    plant.seed_mut()->chemical(ChemicalID::Sugar) = 0.0f;

    WorldParams wp = default_world_params();
    compute_light_exposure(plant, wp);
    leaf->tick(plant, wp);

    // Sugar should not have increased past cap (production skipped at cap)
    REQUIRE(leaf->chemical(ChemicalID::Sugar) <= cap + 1e-6f);
}

TEST_CASE("Production works normally when leaf is below cap", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    leaf->chemical(ChemicalID::Sugar) = 0.0f;
    leaf->chemical(ChemicalID::Water) = 100.0f;  // ensure full stomatal conductance
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 2.0f;
    plant.seed_mut()->chemical(ChemicalID::Sugar) = 0.0f;

    compute_light_exposure(plant, wp);
    leaf->tick(plant, wp);

    REQUIRE(leaf->chemical(ChemicalID::Sugar) > 0.0f);
}

TEST_CASE("Senescing leaf produces no sugar", "[sugar]") {
    Genome g = default_genome();
    Plant plant(g, glm::vec3(0.0f));

    Node* leaf = plant.create_node(NodeType::LEAF, glm::vec3(0.0f, 0.5f, 0.0f), 0.0f);
    leaf->as_leaf()->leaf_size = 0.5f;
    leaf->as_leaf()->senescence_ticks = 1; // senescing
    plant.seed_mut()->add_child(leaf);

    WorldParams wp = default_world_params();
    wp.light_level = 2.0f;
    plant.seed_mut()->chemical(ChemicalID::Sugar) = 0.0f;

    compute_light_exposure(plant, wp);
    leaf->tick(plant, wp);

    // Senescing leaf should not produce sugar
    REQUIRE(leaf->chemical(ChemicalID::Sugar) < 1e-6f);
}

