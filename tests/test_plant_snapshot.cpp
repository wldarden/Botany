#include <catch2/catch_test_macros.hpp>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <unordered_map>
#include "serialization/plant_snapshot.h"
#include "engine/genome.h"
#include "engine/plant.h"
#include "engine/engine.h"
#include <memory>
#include "engine/node/node.h"
#include "engine/node/stem_node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"

using namespace botany;

static void deep_compare_nodes(const Node& a, const Node& b) {
    REQUIRE(a.id == b.id);
    REQUIRE(a.type == b.type);
    REQUIRE(a.age == b.age);
    REQUIRE(a.starvation_ticks == b.starvation_ticks);
    REQUIRE(a.dormant_ticks == b.dormant_ticks);
    REQUIRE(a.radius == b.radius);
    REQUIRE(a.position == b.position);
    REQUIRE(a.offset == b.offset);
    REQUIRE(a.rest_offset == b.rest_offset);
    REQUIRE(a.ever_active == b.ever_active);
    REQUIRE(a.local().chemicals.size() == b.local().chemicals.size());
    for (const auto& kv : a.local().chemicals) {
        REQUIRE(b.local().chemical(kv.first) == kv.second);
    }
    // Parent id match (NOT pointer)
    uint32_t ap = a.parent ? a.parent->id : UINT32_MAX;
    uint32_t bp = b.parent ? b.parent->id : UINT32_MAX;
    REQUIRE(ap == bp);
}

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

TEST_CASE("leaf trailer round-trips", "[plant_snapshot][node]") {
    LeafNode l(10, glm::vec3(0.0f), 0.01f);
    l.leaf_size = 0.5f;
    l.light_exposure = 0.75f;
    l.senescence_ticks = 12;
    l.deficit_ticks = 5;
    l.facing = glm::vec3(0.0f, 0.707f, 0.707f);

    std::stringstream ss;
    write_leaf_trailer(ss, l);

    std::stringstream rs(ss.str());
    LeafTrailer t = read_leaf_trailer(rs);
    REQUIRE(t.leaf_size == 0.5f);
    REQUIRE(t.light_exposure == 0.75f);
    REQUIRE(t.senescence_ticks == 12u);
    REQUIRE(t.deficit_ticks == 5u);
    REQUIRE(t.facing == glm::vec3(0.0f, 0.707f, 0.707f));
}

TEST_CASE("apical trailer round-trips", "[plant_snapshot][node]") {
    ApicalNode a(11, glm::vec3(0.0f), 0.01f);
    a.active = false;
    a.is_primary = true;
    a.growth_dir = glm::vec3(0.1f, 0.9f, 0.0f);
    a.ticks_since_last_node = 7;

    std::stringstream ss;
    write_apical_trailer(ss, a);

    std::stringstream rs(ss.str());
    ApicalTrailer t = read_apical_trailer(rs);
    REQUIRE(t.active == false);
    REQUIRE(t.is_primary == true);
    REQUIRE(t.growth_dir == glm::vec3(0.1f, 0.9f, 0.0f));
    REQUIRE(t.ticks_since_last_node == 7u);
}

TEST_CASE("root_apical trailer round-trips", "[plant_snapshot][node]") {
    RootApicalNode r(12, glm::vec3(0.0f), 0.005f);
    r.active = true;
    r.is_primary = true;
    r.growth_dir = glm::vec3(0.0f, -1.0f, 0.0f);
    r.ticks_since_last_node = 3;
    r.internodes_spawned = 8;

    std::stringstream ss;
    write_root_apical_trailer(ss, r);

    std::stringstream rs(ss.str());
    RootApicalTrailer t = read_root_apical_trailer(rs);
    REQUIRE(t.active == true);
    REQUIRE(t.is_primary == true);
    REQUIRE(t.growth_dir == glm::vec3(0.0f, -1.0f, 0.0f));
    REQUIRE(t.ticks_since_last_node == 3u);
    REQUIRE(t.internodes_spawned == 8u);
}

TEST_CASE("conduit pool chemicals round-trip", "[plant_snapshot][node]") {
    StemNode s(42, glm::vec3(0.0f), 0.02f);
    s.phloem()->chemical(ChemicalID::Sugar) = 0.05f;
    s.xylem()->chemical(ChemicalID::Water) = 0.3f;
    s.xylem()->chemical(ChemicalID::Cytokinin) = 0.015f;

    std::stringstream ss;
    write_conduit_pools(ss, s);

    std::stringstream rs(ss.str());
    ConduitPools p = read_conduit_pools(rs);
    REQUIRE(p.phloem.at(ChemicalID::Sugar) == 0.05f);
    REQUIRE(p.xylem.at(ChemicalID::Water)  == 0.3f);
    REQUIRE(p.xylem.at(ChemicalID::Cytokinin) == 0.015f);
}

TEST_CASE("Plant::from_empty yields an empty plant", "[plant_snapshot][plant]") {
    auto p = Plant::from_empty(default_genome());
    REQUIRE(p != nullptr);
    REQUIRE(p->node_count() == 0);
}

TEST_CASE("Plant::install_node + set_next_id populate nodes manually", "[plant_snapshot][plant]") {
    auto p = Plant::from_empty(default_genome());
    auto s = std::make_unique<StemNode>(99, glm::vec3(0.0f), 0.02f);
    Node* raw = s.get();
    p->install_node(std::move(s));
    p->set_next_id(100);
    REQUIRE(p->node_count() == 1);
    REQUIRE(p->seed() == raw);
    REQUIRE(p->next_id() == 100u); // next_id() increments; now 101
}

TEST_CASE("save_plant_snapshot writes a file with valid header", "[plant_snapshot][save]") {
    Plant plant(default_genome(), glm::vec3(0.0f));
    auto tmp = std::filesystem::temp_directory_path() / "botany_snap_test";
    std::filesystem::remove_all(tmp);

    SaveResult r = save_plant_snapshot(plant, /*engine_tick=*/42, tmp.string());
    REQUIRE(r.ok);
    REQUIRE(!r.path.empty());
    REQUIRE(std::filesystem::exists(r.path));

    // Re-open and check magic + version + tick.
    std::ifstream in(r.path, std::ios::binary);
    REQUIRE(plant_snapshot_check_magic(in));
    uint32_t version = 0;
    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    REQUIRE(version == PLANT_SNAPSHOT_VERSION);
    uint64_t engine_tick = 0;
    in.read(reinterpret_cast<char*>(&engine_tick), sizeof(engine_tick));
    REQUIRE(engine_tick == 42u);

    std::filesystem::remove_all(tmp);
}

TEST_CASE("save+load round-trip preserves tiny plant", "[plant_snapshot][roundtrip]") {
    Engine engine;
    engine.create_plant(default_genome(), glm::vec3(0.0f));
    for (int i = 0; i < 20; i++) engine.tick();

    const Plant& original = engine.get_plant(0);
    auto tmp = std::filesystem::temp_directory_path() / "botany_rt_test";
    std::filesystem::remove_all(tmp);

    SaveResult sr = save_plant_snapshot(original, engine.get_tick(), tmp.string());
    REQUIRE(sr.ok);

    LoadedPlant lp = load_plant_snapshot(sr.path, std::nullopt);
    REQUIRE(lp.plant != nullptr);
    REQUIRE(lp.plant->node_count() == original.node_count());
    REQUIRE(lp.engine_tick == engine.get_tick());

    std::unordered_map<uint32_t, const Node*> orig_by_id;
    original.for_each_node([&](const Node& n) { orig_by_id[n.id] = &n; });
    lp.plant->for_each_node([&](const Node& n) {
        auto it = orig_by_id.find(n.id);
        REQUIRE(it != orig_by_id.end());
        deep_compare_nodes(*it->second, n);
    });

    std::filesystem::remove_all(tmp);
}

TEST_CASE("continuation equivalence: load + run vs direct-run match", "[plant_snapshot][continuation]") {
    // A: grow 30 ticks
    Engine engineA;
    engineA.create_plant(default_genome(), glm::vec3(0.0f));
    for (int i = 0; i < 30; i++) engineA.tick();

    // Save A at tick 30, continue to tick 60.
    auto tmp = std::filesystem::temp_directory_path() / "botany_cont_test";
    std::filesystem::remove_all(tmp);
    SaveResult sr = save_plant_snapshot(engineA.get_plant(0), engineA.get_tick(), tmp.string());
    REQUIRE(sr.ok);
    for (int i = 0; i < 30; i++) engineA.tick();

    // B: load A's snapshot at tick 30 and continue for 30 ticks.
    Engine engineB;
    LoadedPlant lp = load_plant_snapshot(sr.path, std::nullopt);
    engineB.adopt_plant(std::move(lp.plant));
    engineB.set_tick(static_cast<uint32_t>(lp.engine_tick));
    for (int i = 0; i < 30; i++) engineB.tick();

    // Compare terminal state.
    REQUIRE(engineA.get_plant(0).node_count() == engineB.get_plant(0).node_count());
    std::unordered_map<uint32_t, const Node*> a_by_id;
    engineA.get_plant(0).for_each_node([&](const Node& n) { a_by_id[n.id] = &n; });
    engineB.get_plant(0).for_each_node([&](const Node& n) {
        auto it = a_by_id.find(n.id);
        REQUIRE(it != a_by_id.end());
        deep_compare_nodes(*it->second, n);
    });

    std::filesystem::remove_all(tmp);
}
