#include <catch2/catch_test_macros.hpp>
#include <sstream>
#include "serialization/plant_snapshot.h"
#include "engine/genome.h"
#include "engine/plant.h"
#include "engine/node/node.h"
#include "engine/node/stem_node.h"

using namespace botany;

TEST_CASE("plant_snapshot header magic roundtrips", "[plant_snapshot][header]") {
    std::stringstream ss;
    plant_snapshot_write_magic(ss);
    bool ok = plant_snapshot_check_magic(ss);
    REQUIRE(ok);
}

TEST_CASE("plant_snapshot check_magic rejects wrong magic", "[plant_snapshot][header]") {
    std::stringstream ss;
    ss.write("XXXX", 4);
    bool ok = plant_snapshot_check_magic(ss);
    REQUIRE_FALSE(ok);
}

TEST_CASE("genome binary round-trip is identity", "[plant_snapshot][genome]") {
    Genome original = default_genome();
    original.growth_rate = 0.01234f;
    original.shoot_plastochron = 123u;
    original.meristem_dormancy_death_ticks = 4567u;
    original.max_internode_length = 2.7f;

    std::stringstream ss;
    write_genome_binary(ss, original);

    std::stringstream rs(ss.str());
    Genome loaded = read_genome_binary(rs);

    REQUIRE(loaded.growth_rate == original.growth_rate);
    REQUIRE(loaded.shoot_plastochron == original.shoot_plastochron);
    REQUIRE(loaded.meristem_dormancy_death_ticks == original.meristem_dormancy_death_ticks);
    REQUIRE(loaded.max_internode_length == original.max_internode_length);
    // Canary: default field we didn't touch.
    REQUIRE(loaded.branch_angle == original.branch_angle);
}

TEST_CASE("node common fields round-trip through binary", "[plant_snapshot][node]") {
    // Build one StemNode and populate every common field.
    StemNode s(42, glm::vec3(1.0f, 2.0f, 3.0f), 0.025f);
    s.parent = nullptr;
    s.offset      = glm::vec3(0.1f, 0.2f, 0.3f);
    s.rest_offset = glm::vec3(0.4f, 0.5f, 0.6f);
    s.position    = glm::vec3(1.0f, 2.0f, 3.0f);
    s.age = 77;
    s.starvation_ticks = 11;
    s.dormant_ticks = 0;
    s.ever_active = true;
    s.local().chemical(ChemicalID::Sugar)  = 1.5f;
    s.local().chemical(ChemicalID::Auxin)  = 0.25f;
    s.local().chemical(ChemicalID::Water)  = 2.0f;

    std::stringstream ss;
    write_node_common(ss, s, /*parent_id=*/ UINT32_MAX);

    std::stringstream rs(ss.str());
    NodeCommonRecord rec = read_node_common(rs);

    REQUIRE(rec.id == 42u);
    REQUIRE(rec.parent_id == UINT32_MAX);
    REQUIRE(rec.type == NodeType::STEM);
    REQUIRE(rec.age == 77u);
    REQUIRE(rec.starvation_ticks == 11u);
    REQUIRE(rec.dormant_ticks == 0u);
    REQUIRE(rec.radius == 0.025f);
    REQUIRE(rec.position == glm::vec3(1.0f, 2.0f, 3.0f));
    REQUIRE(rec.offset == glm::vec3(0.1f, 0.2f, 0.3f));
    REQUIRE(rec.rest_offset == glm::vec3(0.4f, 0.5f, 0.6f));
    REQUIRE(rec.ever_active == true);
    REQUIRE(rec.local_chemicals.at(ChemicalID::Sugar) == 1.5f);
    REQUIRE(rec.local_chemicals.at(ChemicalID::Auxin) == 0.25f);
    REQUIRE(rec.local_chemicals.at(ChemicalID::Water) == 2.0f);
}
