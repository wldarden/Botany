// tests/test_engine.cpp
#include <catch2/catch_test_macros.hpp>
#include "engine/engine.h"
#include "engine/sugar.h"

using namespace botany;

TEST_CASE("Engine starts at tick 0", "[engine]") {
    Engine engine;
    REQUIRE(engine.get_tick() == 0);
}

TEST_CASE("Engine creates a plant and returns valid ID", "[engine]") {
    Engine engine;
    Genome g = default_genome();
    PlantID id = engine.create_plant(g, glm::vec3(0.0f));

    const Plant& plant = engine.get_plant(id);
    REQUIRE(plant.node_count() == 3);
}

TEST_CASE("Engine tick advances tick counter", "[engine]") {
    Engine engine;
    engine.tick();
    REQUIRE(engine.get_tick() == 1);
    engine.tick();
    REQUIRE(engine.get_tick() == 2);
}

TEST_CASE("Engine tick grows the plant", "[engine]") {
    Engine engine;
    Genome g = default_genome();
    // Raise minimum cap so tip nodes can store enough sugar to grow
    g.sugar_cap_minimum = 1.0f;
    PlantID id = engine.create_plant(g, glm::vec3(0.0f));

    const Node* shoot_before = nullptr;
    engine.get_plant(id).for_each_node([&](const Node& n) {
        if (n.type == NodeType::APICAL) {
            shoot_before = &n;
        }
    });
    float y_before = shoot_before->position.y;

    const Genome& genome = engine.get_plant(id).genome();
    engine.get_plant_mut(id).for_each_node_mut([&genome](Node& n) {
        n.local().chemical(ChemicalID::Sugar) = sugar_cap(n, genome);
        n.local().chemical(ChemicalID::Water) = water_cap(n, genome);
    });
    engine.tick();

    float y_after = shoot_before->position.y;
    REQUIRE(y_after > y_before);
}

TEST_CASE("Engine runs multiple ticks and plant grows complex structure", "[engine]") {
    Engine engine;
    Genome g = default_genome();
    g.growth_rate = 0.2f;
    g.shoot_plastochron = 1; // spawn internode every tick for fast testing
    g.root_plastochron = 1;
    // Raise minimum cap so tip nodes can store enough sugar to grow
    g.sugar_cap_minimum = 1.0f;
    PlantID id = engine.create_plant(g, glm::vec3(0.0f));

    for (int i = 0; i < 5; i++) {
        const Genome& genome = engine.get_plant(id).genome();
        engine.get_plant_mut(id).for_each_node_mut([&genome](Node& n) {
            n.local().chemical(ChemicalID::Sugar) = sugar_cap(n, genome);
            n.local().chemical(ChemicalID::Water) = water_cap(n, genome);
            n.local().chemical(ChemicalID::Cytokinin) = 1.0f;
        });
        engine.tick();
    }

    const Plant& plant = engine.get_plant(id);
    REQUIRE(plant.node_count() > 3);

    bool found_auxin = false;
    plant.for_each_node([&](const Node& n) {
        if (n.local().chemical(ChemicalID::Auxin) > 0.0f) found_auxin = true;
    });
    REQUIRE(found_auxin);
}

TEST_CASE("Engine supports multiple plants", "[engine]") {
    Engine engine;
    Genome g = default_genome();
    PlantID id1 = engine.create_plant(g, glm::vec3(-5.0f, 0.0f, 0.0f));
    PlantID id2 = engine.create_plant(g, glm::vec3(5.0f, 0.0f, 0.0f));

    REQUIRE(id1 != id2);

    engine.tick();

    REQUIRE(engine.get_plant(id1).seed()->position.x < 0.0f);
    REQUIRE(engine.get_plant(id2).seed()->position.x > 0.0f);
}
