#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <glm/vec3.hpp>
#include <cmath>

#include "engine/compression.h"
#include "engine/plant.h"
#include "engine/genome.h"

using namespace botany;
using Catch::Approx;

TEST_CASE("compress_plant on a fresh young plant reports no merges", "[compression][stub]") {
    Plant plant(default_genome(), glm::vec3(0.0f));
    CompressionParams params;
    CompressionResult r = compress_plant(plant, params);
    REQUIRE(r.merges_performed == 0u);
    // One detecting pass runs; nothing mergeable → immediate bail-out.
    REQUIRE(r.passes_run == 1u);
}

#include <memory>
#include "engine/node/stem_node.h"
#include "engine/node/root_node.h"
#include "engine/node/tissues/leaf.h"
#include "engine/node/tissues/apical.h"

// Helpers: build a 3-node chain seed → A → B so we can test predicates on the
// (A, B) pair.  Returns raw pointers to the three nodes; plant owns them.
struct Chain3 {
    std::unique_ptr<Plant> plant;
    Node* seed;
    Node* A;
    Node* B;
};

static Chain3 build_stem_chain(const Genome& g,
                               float radiusA, float radiusB,
                               const glm::vec3& offsetA, const glm::vec3& offsetB,
                               uint32_t ageA, uint32_t ageB) {
    Chain3 out;
    out.plant = Plant::from_empty(g);

    auto seed = std::make_unique<StemNode>(0, glm::vec3(0.0f), g.initial_radius);
    out.seed = seed.get();
    seed->offset = glm::vec3(0.0f);

    auto a = std::make_unique<StemNode>(1, glm::vec3(0.0f), radiusA);
    out.A = a.get();
    a->parent = out.seed;
    a->offset = offsetA;
    a->age = ageA;
    out.seed->children.push_back(out.A);

    auto b = std::make_unique<StemNode>(2, glm::vec3(0.0f), radiusB);
    out.B = b.get();
    b->parent = out.A;
    b->offset = offsetB;
    b->age = ageB;
    out.A->children.push_back(out.B);

    out.plant->install_node(std::move(seed));
    out.plant->install_node(std::move(a));
    out.plant->install_node(std::move(b));
    out.plant->set_next_id(3);
    return out;
}

TEST_CASE("can_merge: same-type + mature + similar radius + straight passes", "[compression][predicate]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g,
        /*rA=*/0.030f, /*rB=*/0.030f,
        /*offA=*/glm::vec3(0.0f, 0.2f, 0.0f),
        /*offB=*/glm::vec3(0.0f, 0.2f, 0.0f),
        /*ageA=*/g.internode_maturation_ticks + 10,
        /*ageB=*/g.internode_maturation_ticks + 10);
    CompressionParams params;
    REQUIRE(can_merge(*c.A, *c.B, g, params));
}

TEST_CASE("can_merge: seed (no grandparent) is rejected", "[compression][predicate]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.03f, 0.03f,
        glm::vec3(0, 0.2f, 0), glm::vec3(0, 0.2f, 0),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    CompressionParams params;
    // Try to merge seed + A (seed has no parent → must be refused)
    REQUIRE_FALSE(can_merge(*c.seed, *c.A, g, params));
}

TEST_CASE("can_merge: bent chain (angle too large) is rejected", "[compression][predicate]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.03f, 0.03f,
        glm::vec3(0.0f, 0.2f, 0.0f),
        glm::vec3(0.2f, 0.0f, 0.0f), // 90° bend
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    CompressionParams params;
    REQUIRE_FALSE(can_merge(*c.A, *c.B, g, params));
}

TEST_CASE("can_merge: radius ratio above threshold is rejected", "[compression][predicate]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g,
        /*rA=*/0.030f, /*rB=*/0.060f, // 2× mismatch
        glm::vec3(0, 0.2f, 0), glm::vec3(0, 0.2f, 0),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    CompressionParams params;
    REQUIRE_FALSE(can_merge(*c.A, *c.B, g, params));
}

TEST_CASE("can_merge: still-elongating (young) nodes are rejected", "[compression][predicate]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.03f, 0.03f,
        glm::vec3(0, 0.2f, 0), glm::vec3(0, 0.2f, 0),
        /*ageA=*/10, /*ageB=*/10);
    CompressionParams params;
    REQUIRE_FALSE(can_merge(*c.A, *c.B, g, params));
}

TEST_CASE("can_merge: combined length above threshold is rejected", "[compression][predicate]") {
    Genome g = default_genome();
    // Two segments each at max_internode_length → combined = 2× max; threshold auto is 2×.
    // Make each slightly over max/2 so the sum exceeds the cap.
    float L = g.max_internode_length * 1.1f;
    auto c = build_stem_chain(g, 0.03f, 0.03f,
        glm::vec3(0, L, 0), glm::vec3(0, L, 0),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    CompressionParams params;
    // Default resolves to 2 * max_internode_length; 2 * 1.1 * max > 2 * max.
    REQUIRE_FALSE(can_merge(*c.A, *c.B, g, params));
}

TEST_CASE("can_merge: parent with 2 structural children is rejected", "[compression][predicate]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.03f, 0.03f,
        glm::vec3(0, 0.2f, 0), glm::vec3(0, 0.2f, 0),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    // Add a SECOND structural child to A — another stem branch.
    auto branch = std::make_unique<StemNode>(99, glm::vec3(0.0f), 0.03f);
    branch->parent = c.A;
    branch->offset = glm::vec3(0.1f, 0.1f, 0.0f);
    branch->age = g.internode_maturation_ticks + 10;
    c.A->children.push_back(branch.get());
    c.plant->install_node(std::move(branch));

    CompressionParams params;
    REQUIRE_FALSE(can_merge(*c.A, *c.B, g, params));
}

TEST_CASE("can_merge: parent with leaf + bud + one structural child passes", "[compression][predicate]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.03f, 0.03f,
        glm::vec3(0, 0.2f, 0), glm::vec3(0, 0.2f, 0),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    // Attach a leaf and a dormant axillary bud to A — these must NOT count as
    // structural children (leaves + meristems don't block compression).
    auto leaf = std::make_unique<LeafNode>(50, glm::vec3(0.0f), 0.0f);
    leaf->parent = c.A;
    leaf->offset = glm::vec3(0.1f, 0.0f, 0.0f);
    c.A->children.push_back(leaf.get());
    c.plant->install_node(std::move(leaf));

    auto bud = std::make_unique<ApicalNode>(51, glm::vec3(0.0f), 0.01f);
    bud->parent = c.A;
    bud->offset = glm::vec3(-0.1f, 0.0f, 0.0f);
    if (auto* a = bud->as_apical()) a->active = false;
    c.A->children.push_back(bud.get());
    c.plant->install_node(std::move(bud));

    CompressionParams params;
    REQUIRE(can_merge(*c.A, *c.B, g, params));
}

TEST_CASE("can_merge: stem cannot absorb root-type child", "[compression][predicate]") {
    Genome g = default_genome();
    auto plant = Plant::from_empty(g);
    auto seed = std::make_unique<StemNode>(0, glm::vec3(0.0f), g.initial_radius);
    auto stem = std::make_unique<StemNode>(1, glm::vec3(0.0f), 0.03f);
    stem->parent = seed.get();
    stem->offset = glm::vec3(0.0f, 0.2f, 0.0f);
    stem->age = g.internode_maturation_ticks + 10;
    seed->children.push_back(stem.get());
    auto root = std::make_unique<RootNode>(2, glm::vec3(0.0f), 0.03f);
    root->parent = stem.get();
    root->offset = glm::vec3(0.0f, 0.2f, 0.0f);
    root->age = g.internode_maturation_ticks + 10;
    stem->children.push_back(root.get());

    Node* stem_ptr = stem.get();
    Node* root_ptr = root.get();
    plant->install_node(std::move(seed));
    plant->install_node(std::move(stem));
    plant->install_node(std::move(root));
    plant->set_next_id(3);

    CompressionParams params;
    REQUIRE_FALSE(can_merge(*stem_ptr, *root_ptr, g, params));
}

// Forward-declaration of the test-visible helper.  compress_plant calls this
// internally but the tests exercise it directly for precise assertions.
namespace botany {
void merge_pair(Plant& plant, Node& parent, Node& child, CompressionResult& result);
}

TEST_CASE("merge_pair absorbs C's geometry into P", "[compression][merge]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.04f, 0.04f,
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);

    CompressionResult result;
    merge_pair(*c.plant, *c.A, *c.B, result);
    c.plant->flush_removals();

    // Post-merge A's offset is the grandparent→C vector.
    REQUIRE(c.A->offset == glm::vec3(0.0f, 2.0f, 0.0f));
    // Volume-preserving radius of two equal-radius, equal-length segments is the same radius.
    REQUIRE(c.A->radius == Approx(0.04f).epsilon(1e-4));
    REQUIRE(c.plant->node_count() == 2u); // seed + A; C removed
    REQUIRE(result.merges_performed == 1u);
}

TEST_CASE("merge_pair uses volume-preserving radius on unequal r/L", "[compression][merge]") {
    Genome g = default_genome();
    // r1=0.05 L1=1.0, r2=0.03 L2=1.0 → r_merged = sqrt((0.0025 + 0.0009) / 2) = sqrt(0.0017)
    auto c = build_stem_chain(g, 0.05f, 0.03f,
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);

    // Temporarily loosen the radius ratio gate so the test can proceed.
    // (can_merge defaults reject |0.05-0.03|/0.05 = 0.4 > 0.20.)
    // merge_pair doesn't re-check; it's the caller's responsibility.
    CompressionResult result;
    merge_pair(*c.plant, *c.A, *c.B, result);
    c.plant->flush_removals();

    float expected = std::sqrt((0.05f*0.05f + 0.03f*0.03f) / 2.0f);
    REQUIRE(c.A->radius == Approx(expected).epsilon(1e-4));
}

TEST_CASE("merge_pair sums local chemicals", "[compression][merge]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.04f, 0.04f,
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    c.A->local().chemical(ChemicalID::Sugar) = 0.003f;
    c.B->local().chemical(ChemicalID::Sugar) = 0.005f;
    c.A->local().chemical(ChemicalID::Water) = 2.0f;
    c.B->local().chemical(ChemicalID::Water) = 3.0f;

    CompressionResult result;
    merge_pair(*c.plant, *c.A, *c.B, result);
    c.plant->flush_removals();

    // Caps on a 0.04 dm radius, L=2 dm stem are large enough to not clamp these values.
    REQUIRE(c.A->local().chemical(ChemicalID::Sugar) == Approx(0.008f).epsilon(1e-5));
    REQUIRE(c.A->local().chemical(ChemicalID::Water) == Approx(5.0f).epsilon(1e-5));
    REQUIRE(result.delta_sugar == 0.0f); // no clamping
    REQUIRE(result.delta_water == 0.0f);
}

TEST_CASE("merge_pair sums conduit pool chemicals", "[compression][merge]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.04f, 0.04f,
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    c.A->phloem()->chemical(ChemicalID::Sugar) = 0.001f;
    c.B->phloem()->chemical(ChemicalID::Sugar) = 0.002f;
    c.A->xylem()->chemical(ChemicalID::Water) = 0.5f;
    c.B->xylem()->chemical(ChemicalID::Water) = 0.7f;

    CompressionResult result;
    merge_pair(*c.plant, *c.A, *c.B, result);
    c.plant->flush_removals();

    REQUIRE(c.A->phloem()->chemical(ChemicalID::Sugar) == Approx(0.003f).epsilon(1e-5));
    REQUIRE(c.A->xylem()->chemical(ChemicalID::Water)  == Approx(1.2f).epsilon(1e-5));
}

TEST_CASE("merge_pair reparents children of C to P with same offsets", "[compression][merge]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.04f, 0.04f,
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    // Attach a leaf and a dormant bud to C — both should migrate to A.
    auto leaf = std::make_unique<LeafNode>(100, glm::vec3(0.0f), 0.0f);
    leaf->parent = c.B;
    leaf->offset = glm::vec3(0.1f, 0.0f, 0.0f);
    Node* leaf_ptr = leaf.get();
    c.B->children.push_back(leaf_ptr);
    c.plant->install_node(std::move(leaf));

    auto bud = std::make_unique<ApicalNode>(101, glm::vec3(0.0f), 0.01f);
    bud->parent = c.B;
    bud->offset = glm::vec3(-0.1f, 0.0f, 0.0f);
    Node* bud_ptr = bud.get();
    if (auto* a = bud_ptr->as_apical()) a->active = false;
    c.B->children.push_back(bud_ptr);
    c.plant->install_node(std::move(bud));

    CompressionResult result;
    merge_pair(*c.plant, *c.A, *c.B, result);
    c.plant->flush_removals();

    // Both leaf and bud now have A as parent; their offsets unchanged.
    REQUIRE(leaf_ptr->parent == c.A);
    REQUIRE(leaf_ptr->offset == glm::vec3(0.1f, 0.0f, 0.0f));
    REQUIRE(bud_ptr->parent == c.A);
    REQUIRE(bud_ptr->offset == glm::vec3(-0.1f, 0.0f, 0.0f));

    // A's children now include leaf + bud (and no B).
    bool found_leaf = false, found_bud = false;
    for (Node* cc : c.A->children) {
        if (cc == leaf_ptr) found_leaf = true;
        if (cc == bud_ptr)  found_bud = true;
    }
    REQUIRE(found_leaf);
    REQUIRE(found_bud);
    REQUIRE(c.A->children.size() == 2u); // leaf + bud; B removed
}

TEST_CASE("merge_pair carries auxin_flow_bias from B→child into A→child", "[compression][merge]") {
    Genome g = default_genome();
    auto c = build_stem_chain(g, 0.04f, 0.04f,
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        g.internode_maturation_ticks + 10, g.internode_maturation_ticks + 10);
    // Attach a grandchild to B so bias transfer has somewhere to land.
    auto gc = std::make_unique<StemNode>(200, glm::vec3(0.0f), 0.04f);
    gc->parent = c.B;
    gc->offset = glm::vec3(0.0f, 1.0f, 0.0f);
    gc->age    = g.internode_maturation_ticks + 10;
    Node* gc_ptr = gc.get();
    c.B->children.push_back(gc_ptr);
    c.plant->install_node(std::move(gc));

    c.B->auxin_flow_bias[gc_ptr] = 0.73f;
    c.A->auxin_flow_bias[c.B]    = 0.44f; // will be destroyed

    CompressionResult result;
    merge_pair(*c.plant, *c.A, *c.B, result);
    c.plant->flush_removals();

    auto it = c.A->auxin_flow_bias.find(gc_ptr);
    REQUIRE(it != c.A->auxin_flow_bias.end());
    REQUIRE(it->second == 0.73f);
    // A→B entry is gone.
    REQUIRE(c.A->auxin_flow_bias.find(c.B) == c.A->auxin_flow_bias.end());
}

TEST_CASE("compress_plant merges a straight mature chain", "[compression][scan]") {
    // Build seed → A → B → C → D, all mature same radius, all straight.
    Genome g = default_genome();
    auto plant = Plant::from_empty(g);
    auto seed = std::make_unique<StemNode>(0, glm::vec3(0.0f), g.initial_radius);
    Node* seed_ptr = seed.get();
    plant->install_node(std::move(seed));

    auto add_child = [&](Node* parent, uint32_t id, const glm::vec3& offset) -> Node* {
        auto n = std::make_unique<StemNode>(id, glm::vec3(0.0f), 0.04f);
        n->parent = parent;
        n->offset = offset;
        n->age    = g.internode_maturation_ticks + 50;
        Node* raw = n.get();
        parent->children.push_back(raw);
        plant->install_node(std::move(n));
        return raw;
    };
    Node* A = add_child(seed_ptr, 1, glm::vec3(0.0f, 0.1f, 0.0f));
    Node* B = add_child(A,        2, glm::vec3(0.0f, 0.1f, 0.0f));
    Node* C = add_child(B,        3, glm::vec3(0.0f, 0.1f, 0.0f));
    Node* D = add_child(C,        4, glm::vec3(0.0f, 0.1f, 0.0f));
    (void)B; (void)D;
    plant->set_next_id(5);
    REQUIRE(plant->node_count() == 5u); // seed + 4

    CompressionParams params;
    CompressionResult r = compress_plant(*plant, params);
    REQUIRE(r.merges_performed > 0u);
    REQUIRE(plant->node_count() < 5u);
}

TEST_CASE("compress_plant is idempotent — second run merges nothing", "[compression][scan]") {
    Genome g = default_genome();
    auto plant = Plant::from_empty(g);
    auto seed = std::make_unique<StemNode>(0, glm::vec3(0.0f), g.initial_radius);
    Node* seed_ptr = seed.get();
    plant->install_node(std::move(seed));
    auto add = [&](Node* p, uint32_t id) {
        auto n = std::make_unique<StemNode>(id, glm::vec3(0.0f), 0.04f);
        n->parent = p; n->offset = glm::vec3(0.0f, 0.1f, 0.0f);
        n->age = g.internode_maturation_ticks + 50;
        Node* raw = n.get();
        p->children.push_back(raw);
        plant->install_node(std::move(n));
        return raw;
    };
    Node* A = add(seed_ptr, 1);
    Node* B = add(A, 2);
    (void)B;
    plant->set_next_id(3);

    CompressionParams params;
    CompressionResult r1 = compress_plant(*plant, params);
    CompressionResult r2 = compress_plant(*plant, params);
    REQUIRE(r1.merges_performed >= 1u);
    REQUIRE(r2.merges_performed == 0u);
}

TEST_CASE("compress_plant respects max_passes cap", "[compression][scan]") {
    Genome g = default_genome();
    auto plant = Plant::from_empty(g);
    auto seed = std::make_unique<StemNode>(0, glm::vec3(0.0f), g.initial_radius);
    Node* seed_ptr = seed.get();
    plant->install_node(std::move(seed));
    Node* last = seed_ptr;
    for (uint32_t id = 1; id < 20; ++id) {
        auto n = std::make_unique<StemNode>(id, glm::vec3(0.0f), 0.04f);
        n->parent = last; n->offset = glm::vec3(0.0f, 0.05f, 0.0f);
        n->age = g.internode_maturation_ticks + 50;
        Node* raw = n.get();
        last->children.push_back(raw);
        plant->install_node(std::move(n));
        last = raw;
    }
    plant->set_next_id(20);

    CompressionParams params;
    params.max_passes = 1;
    CompressionResult r = compress_plant(*plant, params);
    REQUIRE(r.passes_run == 1u);
}

#include "engine/engine.h"

TEST_CASE("Engine::trigger_compression runs immediately", "[compression][engine]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));
    for (int i = 0; i < 200; ++i) engine.tick(); // grow enough internodes
    uint32_t before = engine.get_plant(0).node_count();

    engine.trigger_compression();
    uint32_t after = engine.get_plant(0).node_count();
    // Either the plant has mergeable pairs (after < before) or it doesn't
    // (after == before).  Both are acceptable; we're testing the mechanism.
    REQUIRE(after <= before);
    REQUIRE(engine.last_compression().passes_run >= 1u);
}

TEST_CASE("Engine autocompress fires on interval", "[compression][engine]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));
    engine.enable_autocompress(true);
    engine.set_compression_interval(5);

    // Tick 20 times — interval 5 should trigger 4 compressions (ticks 5,10,15,20).
    for (int i = 0; i < 20; ++i) engine.tick();
    REQUIRE(engine.last_compression().passes_run >= 1u);
}

TEST_CASE("Engine autocompress disabled does not run compression", "[compression][engine]") {
    Engine engine;
    Genome g = default_genome();
    engine.create_plant(g, glm::vec3(0.0f));
    engine.enable_autocompress(false);
    engine.set_compression_interval(5);
    for (int i = 0; i < 20; ++i) engine.tick();
    REQUIRE(engine.last_compression().passes_run == 0u);
}
