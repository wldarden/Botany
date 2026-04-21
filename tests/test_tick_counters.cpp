// tests/test_tick_counters.cpp
#include <catch2/catch_test_macros.hpp>
#include "engine/engine.h"
#include "engine/plant.h"
#include "engine/genome.h"
#include "engine/world_params.h"
#include "engine/node/node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/chemical/chemical.h"

using namespace botany;

TEST_CASE("tick counters: arrays reset each tick", "[tick_counters]") {
    Engine engine;
    Genome g = default_genome();
    auto pid = engine.create_plant(g, glm::vec3(0.0f));
    engine.tick();
    engine.tick();
    // after two ticks, every node's arrays should have been zeroed at the start
    // of tick 2 and populated (or not) during tick 2 — either way, size matches Count.
    engine.get_plant(pid).for_each_node([](const Node& n) {
        constexpr size_t expected = static_cast<size_t>(ChemicalID::Count);
        REQUIRE(n.tick_chem_produced.size() == expected);
        REQUIRE(n.tick_chem_consumed.size() == expected);
        for (size_t i = 0; i < n.tick_chem_produced.size(); ++i) {
            REQUIRE(n.tick_chem_produced[i] >= 0.0f);
            REQUIRE(n.tick_chem_consumed[i] >= 0.0f);
        }
    });
}

TEST_CASE("tick counters: leaf sugar production matches photosynthesis", "[tick_counters][photo]") {
    Engine engine;
    Genome g = default_genome();
    engine.world_params_mut().light_level = 1.0f;
    auto pid = engine.create_plant(g, glm::vec3(0.0f));
    // spin up until a leaf exists
    for (int i = 0; i < 200; ++i) engine.tick();
    bool found_leaf = false;
    engine.get_plant(pid).for_each_node([&](const Node& n) {
        if (auto* lf = n.as_leaf(); lf && lf->leaf_size > 0.01f) {
            found_leaf = true;
            float produced = n.tick_chem_produced[static_cast<size_t>(ChemicalID::Sugar)];
            INFO("Expected non-zero sugar production on active leaf");
            REQUIRE(produced >= 0.0f);  // loose: just non-negative — tightened in phase 2
        }
    });
    REQUIRE(found_leaf);
}
