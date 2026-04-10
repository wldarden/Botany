// src/engine/meristem.cpp
#include "engine/meristem.h"
#include "engine/plant.h"
#include "engine/node.h"
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>
#include <cmath>
#include <random>
#include <vector>

namespace botany {

// Compute growth direction for a node: normalized vector from parent to this node
static glm::vec3 growth_direction(const Node& node) {
    if (node.parent) {
        glm::vec3 dir = node.position - node.parent->position;
        float len = glm::length(dir);
        if (len > 0.0001f) {
            return dir / len;
        }
    }
    return (node.type == NodeType::STEM) ? glm::vec3(0.0f, 1.0f, 0.0f)
                                          : glm::vec3(0.0f, -1.0f, 0.0f);
}

// Compute a branch direction in the node's local reference frame.
static glm::vec3 branch_direction(const glm::vec3& main_dir, float angle, uint32_t seed) {
    glm::vec3 gravity(0.0f, -1.0f, 0.0f);

    glm::vec3 perp;
    float alignment = std::abs(glm::dot(main_dir, gravity));
    if (alignment < 0.99f) {
        perp = glm::normalize(glm::cross(main_dir, gravity));
    } else {
        perp = glm::normalize(glm::cross(main_dir, glm::vec3(1.0f, 0.0f, 0.0f)));
    }

    glm::vec3 perp2 = glm::cross(main_dir, perp);

    float rotate = static_cast<float>(seed) * 2.399f; // golden angle in radians
    float cos_r = std::cos(rotate);
    float sin_r = std::sin(rotate);
    glm::vec3 radial = perp * cos_r + perp2 * sin_r;

    float cos_a = std::cos(angle);
    float sin_a = std::sin(angle);
    return glm::normalize(main_dir * cos_a + radial * sin_a);
}

static std::mt19937 rng{std::random_device{}()};

// Perturb a direction by a random angle up to max_angle radians
static glm::vec3 perturb(const glm::vec3& dir, float max_angle) {
    if (max_angle < 1e-6f) return dir;

    std::uniform_real_distribution<float> angle_dist(0.0f, max_angle);
    std::uniform_real_distribution<float> rot_dist(0.0f, 2.0f * 3.14159f);

    float tilt = angle_dist(rng);
    float rot = rot_dist(rng);

    glm::vec3 perp;
    if (std::abs(dir.y) < 0.9f) {
        perp = glm::normalize(glm::cross(dir, glm::vec3(0.0f, 1.0f, 0.0f)));
    } else {
        perp = glm::normalize(glm::cross(dir, glm::vec3(1.0f, 0.0f, 0.0f)));
    }
    glm::vec3 perp2 = glm::cross(dir, perp);

    glm::vec3 radial = perp * std::cos(rot) + perp2 * std::sin(rot);
    return glm::normalize(dir * std::cos(tilt) + radial * std::sin(tilt));
}

// Compute growth fraction based on available sugar.
// Returns 0.0 if sugar is at or below save threshold.
// Returns 1.0 if sugar is sufficient for max growth.
// Linearly interpolates in between.
static float sugar_growth_fraction(float sugar, float save_threshold, float max_cost) {
    if (max_cost < 1e-6f) return 1.0f; // zero-cost growth always happens
    float available = std::max(sugar - save_threshold, 0.0f);
    return std::min(available / max_cost, 1.0f);
}

// ---------------------------------------------------------------------------
// Subclass tick() implementations
// ---------------------------------------------------------------------------

void ShootApicalMeristem::tick(Node& node, Plant& plant) {
    if (!active) return;

    const Genome& g = plant.genome();

    // Extend with slight angular noise
    glm::vec3 dir = perturb(growth_direction(node), g.growth_noise);

    // Sugar-scaled growth
    float max_cost = g.growth_rate * g.sugar_cost_shoot_growth;
    float gf = sugar_growth_fraction(node.sugar, g.sugar_save_shoot, max_cost);
    if (gf < 1e-6f) return; // no sugar for growth

    float actual_rate = g.growth_rate * gf;
    float actual_cost = actual_rate * g.sugar_cost_shoot_growth;
    node.sugar -= actual_cost;

    node.position += dir * actual_rate;

    // Chain growth: when distance to parent exceeds max, this node becomes
    // an interior node with one axillary meristem + leaf, and a new tip node
    // carries the apical meristem forward.
    if (node.parent) {
        float dist = glm::length(node.position - node.parent->position);
        if (dist > g.max_internode_length) {
            // This node becomes interior — give it an axillary meristem + leaf
            glm::vec3 branch_dir = branch_direction(dir, g.branch_angle, node.id);
            glm::vec3 ax_pos = node.position + branch_dir * g.tip_offset;
            Node* axillary = plant.create_node(NodeType::STEM, ax_pos, g.initial_radius * 0.5f);
            axillary->meristem = plant.create_meristem<ShootAxillaryMeristem>();
            node.add_child(axillary);

            // Create leaf as a separate LEAF node on this interior node
            glm::vec3 leaf_dir = branch_direction(dir, g.branch_angle, node.id + 1000);
            glm::vec3 leaf_pos = node.position + leaf_dir * g.tip_offset;
            Node* leaf_node = plant.create_node(NodeType::LEAF, leaf_pos, 0.0f);
            leaf_node->leaf_size = g.leaf_size;
            node.add_child(leaf_node);

            // Create new tip node slightly ahead along current direction
            // so it inherits the branch's growth direction
            Node* new_tip = plant.create_node(NodeType::STEM, node.position + dir * g.tip_offset, node.radius);
            new_tip->meristem = node.meristem;
            node.meristem = nullptr;
            node.add_child(new_tip);
        }
    }
}

void ShootAxillaryMeristem::tick(Node& node, Plant& plant) {
    if (active) return; // already activated — nothing to do

    const Genome& g = plant.genome();
    // Sense auxin on the parent stem node (where auxin actually flows),
    // not on this side-branch node which never receives basipetal auxin.
    // Shoot axillaries activate when auxin drops low enough —
    // meaning they're far from any active shoot tip.
    float stem_auxin = node.parent ? node.parent->auxin : node.auxin;
    if (stem_auxin < g.auxin_threshold) {
        // Activation costs sugar
        if (node.sugar < g.sugar_cost_activation) return;
        node.sugar -= g.sugar_cost_activation;

        // Replace this dormant axillary with an active shoot apical
        auto* apical = plant.create_meristem<ShootApicalMeristem>();
        apical->ticks_since_last_node = 0;
        node.meristem = apical;
    }
}

void RootApicalMeristem::tick(Node& node, Plant& plant) {
    if (!active) return;

    const Genome& g = plant.genome();

    glm::vec3 dir = perturb(growth_direction(node), g.growth_noise);

    // Gravitropism: roots sense gravity and turn downward when near or above
    // the soil surface. The closer to the surface, the stronger the correction.
    if (node.position.y > -g.root_gravitropism_depth) {
        float exposure = (node.position.y + g.root_gravitropism_depth)
                       / g.root_gravitropism_depth; // 0 at depth, 1+ at surface
        exposure = glm::clamp(exposure, 0.0f, 1.0f);
        float strength = exposure * g.root_gravitropism_strength;
        glm::vec3 down(0.0f, -1.0f, 0.0f);
        dir = glm::normalize(dir + down * strength);
    }

    // Sugar-scaled growth
    float max_cost = g.root_growth_rate * g.sugar_cost_root_growth;
    float gf = sugar_growth_fraction(node.sugar, g.sugar_save_root, max_cost);
    if (gf < 1e-6f) return;

    float actual_rate = g.root_growth_rate * gf;
    float actual_cost = actual_rate * g.sugar_cost_root_growth;
    node.sugar -= actual_cost;

    node.position += dir * actual_rate;

    // Chain growth: same logic as shoot — interior node gets one root axillary
    if (node.parent) {
        float dist = glm::length(node.position - node.parent->position);
        if (dist > g.root_max_internode_length) {
            // This node becomes interior — give it a root axillary meristem
            // (only if we haven't hit the root meristem cap)
            if (!plant.root_meristems_at_cap()) {
                glm::vec3 branch_dir = branch_direction(dir, g.root_branch_angle, node.id);
                glm::vec3 ax_pos = node.position + branch_dir * g.tip_offset;
                Node* axillary = plant.create_node(NodeType::ROOT, ax_pos, g.root_initial_radius * 0.5f);
                axillary->meristem = plant.create_meristem<RootAxillaryMeristem>();
                node.add_child(axillary);
            }

            // Create new tip slightly ahead so it inherits the root's growth direction
            Node* new_tip = plant.create_node(NodeType::ROOT, node.position + dir * g.tip_offset, node.radius);
            new_tip->meristem = node.meristem;
            node.meristem = nullptr;
            node.add_child(new_tip);
        }
    }
}

void RootAxillaryMeristem::tick(Node& node, Plant& plant) {
    if (active) return;

    const Genome& g = plant.genome();
    // Sense cytokinin on the parent root node, not this side-branch node
    float stem_cytokinin = node.parent ? node.parent->cytokinin : node.cytokinin;
    if (stem_cytokinin < g.cytokinin_threshold) {
        // Activation costs sugar
        if (node.sugar < g.sugar_cost_activation) return;
        node.sugar -= g.sugar_cost_activation;

        // Replace this dormant root axillary with an active root apical
        auto* apical = plant.create_meristem<RootApicalMeristem>();
        apical->ticks_since_last_node = 0;
        node.meristem = apical;
    }
}

// ---------------------------------------------------------------------------
// Top-level tick dispatch
// ---------------------------------------------------------------------------

void tick_meristems(Plant& plant) {
    const Genome& g = plant.genome();

    // Collect nodes to tick first, since ticking may add new nodes
    std::vector<Node*> to_tick;
    plant.for_each_node_mut([&](Node& n) {
        n.age++;

        // Secondary growth: interior nodes (no active tip meristem) thicken.
        bool is_active_tip = n.meristem && n.meristem->is_tip() && n.meristem->active;
        if (!is_active_tip && n.type != NodeType::LEAF) {
            float max_cost = g.thickening_rate * g.sugar_cost_thickening;
            float gf = sugar_growth_fraction(n.sugar, g.sugar_save_stem, max_cost);
            if (gf > 1e-6f) {
                float actual_rate = g.thickening_rate * gf;
                float actual_cost = actual_rate * g.sugar_cost_thickening;
                n.sugar -= actual_cost;
                n.radius += actual_rate;
            }
        }
        if (n.meristem) {
            to_tick.push_back(&n);
        }
    });

    for (Node* node : to_tick) {
        if (node->meristem) {
            node->meristem->tick(*node, plant);
        }
    }
}

} // namespace botany
