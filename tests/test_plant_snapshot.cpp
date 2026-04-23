#include <catch2/catch_test_macros.hpp>
#include <sstream>
#include "serialization/plant_snapshot.h"

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
