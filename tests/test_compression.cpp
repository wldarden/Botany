#include <catch2/catch_test_macros.hpp>
#include <glm/vec3.hpp>

#include "engine/compression.h"
#include "engine/plant.h"
#include "engine/genome.h"

using namespace botany;

TEST_CASE("compress_plant on a fresh young plant reports no merges", "[compression][stub]") {
    Plant plant(default_genome(), glm::vec3(0.0f));
    CompressionParams params;
    CompressionResult r = compress_plant(plant, params);
    REQUIRE(r.merges_performed == 0u);
    REQUIRE(r.passes_run == 0u);
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
