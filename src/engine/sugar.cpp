#include "engine/sugar.h"
#include "engine/plant.h"
#include "engine/world_params.h"
#include "engine/node/node.h"
#include "engine/node/leaf_node.h"
#include "engine/node/meristem_node.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <glm/geometric.hpp>
#include <unordered_map>

namespace botany {

float sugar_cap(const Node& node, const Genome& g) {
    float length = std::max(glm::length(node.offset), 0.01f);

    switch (node.type) {
        case NodeType::LEAF: {
            float ls = node.as_leaf()->leaf_size;
            float area = ls * ls;
            return std::max(g.sugar_cap_minimum, area * g.sugar_storage_density_leaf);
        }
        case NodeType::STEM:
        case NodeType::ROOT: {
            float volume = 3.14159f * node.radius * node.radius * length;
            float cap = volume * g.sugar_storage_density_wood;
            // Seed node (no parent): cap must hold initial reserves
            if (!node.parent) {
                cap = std::max(cap, g.seed_sugar);
            }
            return std::max(g.sugar_cap_minimum, cap);
        }
        case NodeType::SHOOT_APICAL:
        case NodeType::SHOOT_AXILLARY:
        case NodeType::ROOT_APICAL:
        case NodeType::ROOT_AXILLARY:
            return g.sugar_cap_minimum;
    }
    return g.sugar_cap_minimum;
}

void produce_sugar(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    // Collect all above-ground nodes as shadow casters, and leaves as targets.
    struct Caster { float x, y, z, shadow_radius; };
    struct LeafTarget { Node* node; float x, y, z, size; };

    std::vector<Caster> casters;
    std::vector<LeafTarget> leaves;

    plant.for_each_node_mut([&](Node& node) {
        // Only above-ground physical nodes cast shadows
        // (roots are underground, meristems are virtual growth points)
        if (node.type != NodeType::ROOT && !node.is_meristem()) {
            float sr = node.radius;
            if (auto* leaf = node.as_leaf()) {
                leaf->light_exposure = 1.0f;  // reset
                sr = leaf->leaf_size;
            }
            if (sr > 1e-6f) {
                casters.push_back({node.position.x, node.position.y,
                                   node.position.z, sr});
            }
        }
        if (auto* leaf = node.as_leaf()) {
            if (leaf->leaf_size > 1e-6f && leaf->senescence_ticks == 0) {
                leaves.push_back({&node, node.position.x, node.position.y,
                                  node.position.z, leaf->leaf_size});
            }
        }
    });

    // For each leaf, accumulate shade from all nodes above it.
    // Beer-Lambert: exposure = exp(-k * accumulated_shade_area)
    float k = world.light_extinction_coeff;
    for (auto& leaf : leaves) {
        float shade = 0.0f;
        for (const auto& c : casters) {
            if (c.y <= leaf.y) continue;  // only shaded by things above
            float dx = c.x - leaf.x;
            float dz = c.z - leaf.z;
            float horiz2 = dx * dx + dz * dz;
            float reach = c.shadow_radius + leaf.size;
            if (horiz2 >= reach * reach) continue;
            float overlap = 1.0f - std::sqrt(horiz2) / reach;
            shade += c.shadow_radius * overlap;
        }
        leaf.node->as_leaf()->light_exposure = std::exp(-k * shade);

        // Leaf angle efficiency: a leaf facing straight up captures full light,
        // one edge-on captures very little. Use dot product of leaf normal
        // with sky direction, clamped to [0,1].
        float angle_efficiency = 1.0f;
        float offset_len = glm::length(leaf.node->offset);
        if (offset_len > 1e-4f) {
            glm::vec3 leaf_normal = leaf.node->offset / offset_len;
            angle_efficiency = std::max(0.0f, leaf_normal.y); // dot with (0,1,0)
        }

        // Feedback inhibition: full storage stops photosynthesis
        float cap = sugar_cap(*leaf.node, g);
        if (leaf.node->sugar >= cap) continue;

        auto* lf = leaf.node->as_leaf();
        leaf.node->sugar += lf->light_exposure * angle_efficiency
                          * world.light_level * lf->leaf_size
                          * g.sugar_production_rate;
    }
}

void consume_sugar(Plant& plant) {
    const Genome& g = plant.genome();
    plant.for_each_node_mut([&](Node& node) {
        float cost = 0.0f;
        float length = std::max(glm::length(node.offset), 0.01f);

        switch (node.type) {
            case NodeType::LEAF: {
                float ls = node.as_leaf()->leaf_size;
                cost = g.sugar_maintenance_leaf * ls * ls;
                break;
            }
            case NodeType::STEM: {
                float volume = 3.14159f * node.radius * node.radius * length;
                cost = g.sugar_maintenance_stem * volume;
                break;
            }
            case NodeType::ROOT: {
                float volume = 3.14159f * node.radius * node.radius * length;
                cost = g.sugar_maintenance_root * volume;
                break;
            }
            case NodeType::SHOOT_APICAL:
            case NodeType::SHOOT_AXILLARY:
            case NodeType::ROOT_APICAL:
            case NodeType::ROOT_AXILLARY:
                break; // meristem tip cost handled below
        }
        if (auto* m = node.as_meristem()) {
            if (m->is_tip() && m->active) {
                cost += g.sugar_maintenance_meristem;
            }
        }
        node.sugar = std::max(0.0f, node.sugar - cost);

        // Safety clamp: cap sugar to node's storage capacity
        float cap = sugar_cap(node, g);
        node.sugar = std::min(node.sugar, cap);

        // Track starvation
        if (node.sugar <= 0.0f) {
            node.starvation_ticks++;
        } else {
            node.starvation_ticks = 0;
        }
    });
}

void prune_starved_nodes(Plant& plant, const WorldParams& world) {
    std::vector<Node*> to_prune;
    plant.for_each_node_mut([&](Node& node) {
        if (node.starvation_ticks >= world.starvation_ticks_max) {
            // Only prune if parent is NOT also starved (prune from the top)
            bool parent_starved = node.parent &&
                node.parent->starvation_ticks >= world.starvation_ticks_max;
            if (!parent_starved && node.parent != nullptr) {  // never prune seed (parent == nullptr)
                to_prune.push_back(&node);
            }
        }
    });

    for (Node* n : to_prune) {
        plant.remove_subtree(n);
    }
}

void grow_leaves(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();
    glm::vec3 up(0.0f, 1.0f, 0.0f);

    plant.for_each_node_mut([&](Node& node) {
        auto* leaf = node.as_leaf();
        if (!leaf) return;

        // Phototropism: rotate leaf offset toward sky
        if (g.leaf_phototropism_rate > 1e-8f) {
            float offset_len = glm::length(node.offset);
            if (offset_len > 1e-4f) {
                glm::vec3 dir = node.offset / offset_len;
                float cos_angle = glm::dot(dir, up);
                if (cos_angle < 0.999f) { // not already pointing up
                    glm::vec3 axis = glm::cross(dir, up);
                    float axis_len = glm::length(axis);
                    if (axis_len > 1e-6f) {
                        axis /= axis_len;
                        float angle_to_up = std::acos(std::min(cos_angle, 1.0f));
                        float turn = std::min(g.leaf_phototropism_rate, angle_to_up);

                        float cost = turn * world.sugar_cost_phototropism;
                        if (node.sugar >= cost) {
                            node.sugar -= cost;
                            float c = std::cos(turn);
                            float s = std::sin(turn);
                            glm::vec3 new_dir = dir * c
                                + glm::cross(axis, dir) * s
                                + axis * glm::dot(axis, dir) * (1.0f - c);
                            node.offset = glm::normalize(new_dir) * offset_len;
                        }
                    }
                }
            }
        }

        // Size growth
        if (leaf->leaf_size >= g.max_leaf_size) return;

        float max_growth = g.leaf_growth_rate;
        float remaining = g.max_leaf_size - leaf->leaf_size;
        float growth = std::min(max_growth, remaining);

        float cost = growth * world.sugar_cost_leaf_growth;
        if (node.sugar < cost) {
            growth = node.sugar / world.sugar_cost_leaf_growth;
            cost = node.sugar;
        }
        if (growth < 1e-7f) return;

        leaf->leaf_size += growth;
        node.sugar -= cost;
    });
}

void transport_sugar(Plant& plant, const WorldParams& world) {
    produce_sugar(plant, world);
    // Sugar diffusion now happens locally in Node::transport_chemicals()
    consume_sugar(plant);
    prune_starved_nodes(plant, world);
}

} // namespace botany
