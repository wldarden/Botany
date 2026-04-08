// src/engine/meristem.cpp
#include "engine/meristem.h"
#include "engine/plant.h"
#include "engine/node.h"
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>
#include <cmath>
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

// Compute a branch direction offset from the main growth direction
static glm::vec3 branch_direction(const glm::vec3& main_dir, float angle, uint32_t seed) {
    glm::vec3 perp;
    if (std::abs(main_dir.y) < 0.9f) {
        perp = glm::normalize(glm::cross(main_dir, glm::vec3(0.0f, 1.0f, 0.0f)));
    } else {
        perp = glm::normalize(glm::cross(main_dir, glm::vec3(1.0f, 0.0f, 0.0f)));
    }

    float rotate = static_cast<float>(seed) * 2.399f; // golden angle
    float cos_r = std::cos(rotate);
    float sin_r = std::sin(rotate);
    glm::vec3 perp2 = glm::cross(main_dir, perp);
    glm::vec3 rotated_perp = perp * cos_r + perp2 * sin_r;

    float cos_a = std::cos(angle);
    float sin_a = std::sin(angle);
    return glm::normalize(main_dir * cos_a + rotated_perp * sin_a);
}

static void tick_shoot_apical(Node& node, Plant& plant) {
    const Genome& g = plant.genome();
    Meristem* m = node.meristem;

    // Extend
    glm::vec3 dir = growth_direction(node);
    node.position += dir * g.growth_rate;

    // Thicken
    node.radius += g.thickening_rate;

    // Node spacing check
    m->ticks_since_last_node++;
    if (m->ticks_since_last_node >= g.internode_spacing) {
        m->ticks_since_last_node = 0;

        glm::vec3 branch_dir = branch_direction(dir, g.branch_angle, node.id);
        glm::vec3 ax_pos = node.position + branch_dir * 0.01f;
        Node* axillary = plant.create_node(NodeType::STEM, ax_pos, g.initial_radius * 0.5f);
        Meristem* ax_m = plant.create_meristem(MeristemType::AXILLARY, false);
        axillary->meristem = ax_m;
        Leaf* leaf = plant.create_leaf(g.leaf_size);
        axillary->leaf = leaf;
        node.add_child(axillary);
    }

    // Chain growth
    if (node.parent) {
        float dist = glm::length(node.position - node.parent->position);
        if (dist > g.max_internode_length) {
            Node* new_tip = plant.create_node(NodeType::STEM, node.position, node.radius);
            new_tip->meristem = node.meristem;
            node.meristem = nullptr;
            node.add_child(new_tip);
            new_tip->meristem->ticks_since_last_node = 0;
        }
    }
}

static void tick_root_apical(Node& node, Plant& plant) {
    const Genome& g = plant.genome();
    Meristem* m = node.meristem;

    glm::vec3 dir = growth_direction(node);
    node.position += dir * g.root_growth_rate;

    m->ticks_since_last_node++;
    if (m->ticks_since_last_node >= g.root_internode_spacing) {
        m->ticks_since_last_node = 0;

        glm::vec3 branch_dir = branch_direction(dir, g.root_branch_angle, node.id);
        glm::vec3 ax_pos = node.position + branch_dir * 0.01f;
        Node* axillary = plant.create_node(NodeType::ROOT, ax_pos, g.initial_radius * 0.3f);
        Meristem* ax_m = plant.create_meristem(MeristemType::ROOT_AXILLARY, false);
        axillary->meristem = ax_m;
        node.add_child(axillary);
    }

    if (node.parent) {
        float dist = glm::length(node.position - node.parent->position);
        if (dist > g.root_max_internode_length) {
            Node* new_tip = plant.create_node(NodeType::ROOT, node.position, node.radius);
            new_tip->meristem = node.meristem;
            node.meristem = nullptr;
            node.add_child(new_tip);
            new_tip->meristem->ticks_since_last_node = 0;
        }
    }
}

static void tick_shoot_axillary(Node& node, const Genome& g) {
    if (node.auxin < g.auxin_threshold && node.cytokinin > g.cytokinin_threshold) {
        node.meristem->type = MeristemType::APICAL;
        node.meristem->active = true;
        node.meristem->ticks_since_last_node = 0;
    }
}

static void tick_root_axillary(Node& node, const Genome& g) {
    if (node.cytokinin < g.cytokinin_threshold) {
        node.meristem->type = MeristemType::ROOT_APICAL;
        node.meristem->active = true;
        node.meristem->ticks_since_last_node = 0;
    }
}

void tick_meristems(Plant& plant) {
    const Genome& g = plant.genome();

    // Collect nodes to tick first, since ticking may add new nodes
    std::vector<Node*> to_tick;
    plant.for_each_node_mut([&](Node& n) {
        n.age++;
        if (n.meristem) {
            to_tick.push_back(&n);
        }
    });

    for (Node* node : to_tick) {
        Meristem* m = node->meristem;
        if (!m) continue;

        switch (m->type) {
            case MeristemType::APICAL:
                if (m->active) tick_shoot_apical(*node, plant);
                break;
            case MeristemType::AXILLARY:
                if (!m->active) tick_shoot_axillary(*node, g);
                break;
            case MeristemType::ROOT_APICAL:
                if (m->active) tick_root_apical(*node, plant);
                break;
            case MeristemType::ROOT_AXILLARY:
                if (!m->active) tick_root_axillary(*node, g);
                break;
        }
    }
}

} // namespace botany
