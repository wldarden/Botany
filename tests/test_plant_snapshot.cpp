#include <catch2/catch_test_macros.hpp>
#include <sstream>
#include "serialization/plant_snapshot.h"
#include "engine/genome.h"

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
