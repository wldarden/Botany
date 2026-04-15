#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <sstream>
#include "engine/engine.h"
#include "engine/node/tissues/leaf.h"
#include "serialization/serializer.h"

using namespace botany;

TEST_CASE("Round-trip: save and load single tick", "[serializer]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));
    engine.tick();

    std::stringstream ss;
    Recording rec;
    rec.genome = g;
    save_tick(ss, engine, 0);

    ss.seekg(0);
    TickSnapshot snap = load_tick(ss);

    REQUIRE(snap.tick_number == 1);
    REQUIRE(snap.nodes.size() == engine.get_plant(0).node_count());
}

TEST_CASE("Round-trip: node positions preserved", "[serializer]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(5.0f, 0.0f, -3.0f));
    engine.tick();
    engine.tick();

    std::stringstream ss;
    save_tick(ss, engine, 0);

    ss.seekg(0);
    TickSnapshot snap = load_tick(ss);

    const Plant& plant = engine.get_plant(0);
    std::vector<const Node*> engine_nodes;
    plant.for_each_node([&](const Node& n) {
        engine_nodes.push_back(&n);
    });

    REQUIRE(snap.nodes.size() == engine_nodes.size());
    for (size_t i = 0; i < snap.nodes.size(); i++) {
        REQUIRE(snap.nodes[i].id == engine_nodes[i]->id);
        REQUIRE(snap.nodes[i].position.x == engine_nodes[i]->position.x);
        REQUIRE(snap.nodes[i].position.y == engine_nodes[i]->position.y);
        REQUIRE(snap.nodes[i].position.z == engine_nodes[i]->position.z);
        REQUIRE(snap.nodes[i].radius == engine_nodes[i]->radius);
        REQUIRE(snap.nodes[i].sugar == engine_nodes[i]->chemical(ChemicalID::Sugar));
        REQUIRE(snap.nodes[i].leaf_size == (engine_nodes[i]->type == NodeType::LEAF ? engine_nodes[i]->as_leaf()->leaf_size : 0.0f));
    }
}

TEST_CASE("Round-trip: parent_id preserved for graph reconstruction", "[serializer]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));
    for (int i = 0; i < 5; i++) engine.tick();

    std::stringstream ss;
    save_tick(ss, engine, 0);

    ss.seekg(0);
    TickSnapshot snap = load_tick(ss);

    REQUIRE(snap.nodes[0].parent_id == UINT32_MAX);
    for (size_t i = 1; i < snap.nodes.size(); i++) {
        REQUIRE(snap.nodes[i].parent_id != UINT32_MAX);
    }
}

TEST_CASE("Save recording header and multiple ticks", "[serializer]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));

    std::stringstream ss;
    save_recording_header(ss, g, 10);

    for (int i = 0; i < 10; i++) {
        engine.tick();
        save_tick(ss, engine, 0);
    }

    ss.seekg(0);
    RecordingHeader header = load_recording_header(ss);
    REQUIRE(header.num_ticks == 10);

    TickSnapshot snap1 = load_tick(ss);
    REQUIRE(snap1.tick_number == 1);
    REQUIRE(snap1.nodes.size() >= 3);
}
