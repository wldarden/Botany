#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <fstream>
#include "engine/genome.h"
#include "engine/genome_io.h"

using namespace botany;

TEST_CASE("load_genome_file parses key=value pairs", "[genome_io]") {
    auto tmp = std::filesystem::temp_directory_path() / "botany_genome_test.txt";
    {
        std::ofstream f(tmp);
        f << "growth_rate=0.012345\n";
        f << "shoot_plastochron=77\n";
        f << "max_leaf_size=1.5\n";
    }

    Genome g = load_genome_file(tmp.string());
    REQUIRE(g.growth_rate == 0.012345f);
    REQUIRE(g.shoot_plastochron == 77u);
    REQUIRE(g.max_leaf_size == 1.5f);
    // Unspecified fields stay at default.
    REQUIRE(g.branch_angle == default_genome().branch_angle);

    std::filesystem::remove(tmp);
}

TEST_CASE("load_genome_file returns default on missing file", "[genome_io]") {
    Genome g = load_genome_file("/nonexistent/path/never_here.txt");
    REQUIRE(g.growth_rate == default_genome().growth_rate);
}
