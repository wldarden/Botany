// src/engine/meristems/meristem.cpp
#include "engine/meristems/meristem.h"
#include "engine/meristems/helpers.h"
#include "engine/plant.h"
#include "engine/node.h"
#include "engine/world_params.h"
#include <algorithm>
#include <glm/geometric.hpp>
#include <vector>

namespace botany {

using namespace meristem_helpers;

void tick_meristems(Plant& plant, const WorldParams& world) {
    const Genome& g = plant.genome();

    // Collect nodes to tick first, since ticking may add new nodes
    std::vector<Node*> to_tick;
    plant.for_each_node_mut([&](Node& n) {
        n.age++;

        // Secondary growth: interior nodes (no active tip meristem) thicken.
        bool is_active_tip = n.meristem && n.meristem->is_tip() && n.meristem->active;
        if (!is_active_tip && n.type != NodeType::LEAF) {
            float max_cost = g.thickening_rate * world.sugar_cost_thickening;
            float gf = sugar_growth_fraction(n.sugar, g.sugar_save_stem, max_cost);
            if (gf > 1e-6f) {
                float actual_rate = g.thickening_rate * gf;
                float actual_cost = actual_rate * world.sugar_cost_thickening;
                n.sugar -= actual_cost;
                n.radius += actual_rate;
            }
        }
        // Intercalary growth: young interior nodes elongate their internode.
        // Eligible = no meristem (not a tip or dormant bud), not a leaf, has parent, young enough.
        if (!n.meristem && n.type != NodeType::LEAF && n.parent) {
            float elong_rate = (n.type == NodeType::STEM) ? g.internode_elongation_rate
                                                          : g.root_internode_elongation_rate;
            uint32_t mat_ticks = (n.type == NodeType::STEM) ? g.internode_maturation_ticks
                                                            : g.root_internode_maturation_ticks;
            float save = (n.type == NodeType::STEM) ? g.sugar_save_stem : g.sugar_save_root;

            if (n.age < mat_ticks && elong_rate > 1e-8f) {
                // GA boosts elongation rate
                float ga_boost = 1.0f + n.gibberellin * g.ga_elongation_sensitivity;
                // Ethylene inhibits elongation
                float eth_inhibit = std::max(0.0f, 1.0f - n.ethylene * g.ethylene_elongation_inhibition);
                float effective_rate = elong_rate * ga_boost * eth_inhibit;

                // GA-modulated max internode length
                float max_len = (n.type == NodeType::STEM) ? g.max_internode_length
                                                           : g.root_max_internode_length;
                max_len *= (1.0f + n.gibberellin * g.ga_length_sensitivity);
                float current_len = glm::length(n.offset);
                if (current_len >= max_len) return;

                float max_cost = effective_rate * world.sugar_cost_elongation;
                float gf = sugar_growth_fraction(n.sugar, save, max_cost);
                if (gf > 1e-6f) {
                    float actual_rate = effective_rate * gf;
                    float actual_cost = actual_rate * world.sugar_cost_elongation;
                    n.sugar -= actual_cost;
                    float len = glm::length(n.offset);
                    if (len > 1e-4f) {
                        n.offset += (n.offset / len) * actual_rate;
                    }
                }
            }
        }

        if (n.meristem) {
            to_tick.push_back(&n);
        }
    });

    for (Node* node : to_tick) {
        if (node->meristem) {
            node->meristem->tick(*node, plant, world);
        }
    }
}

} // namespace botany
