#include "engine/hormone.h"
#include "engine/plant.h"
#include "engine/node.h"

namespace botany {

// Post-order traversal: process children before parent
static void auxin_postorder(Node& node, const Genome& genome) {
    for (Node* child : node.children) {
        auxin_postorder(*child, genome);
    }

    // Production: only shoot apical meristems produce auxin
    if (node.meristem && node.meristem->active &&
        node.meristem->type == MeristemType::APICAL) {
        node.auxin += genome.auxin_production_rate;
    }

    // Transport: send fraction to parent
    if (node.parent) {
        float flow = node.auxin * genome.auxin_transport_rate;
        node.parent->auxin += flow;
        node.auxin -= flow;
    }

    // Decay
    node.auxin *= (1.0f - genome.auxin_decay_rate);
}

void transport_auxin(Plant& plant) {
    auxin_postorder(*plant.seed_mut(), plant.genome());
}

// Pre-order traversal for cytokinin: roots -> tips
// Phase 1 (post-order): collect cytokinin upward from root tips to seed
// Phase 2 (pre-order): distribute cytokinin downward from seed to shoot tips
static void cytokinin_collect(Node& node, const Genome& genome) {
    for (Node* child : node.children) {
        cytokinin_collect(*child, genome);
    }

    // Production: only root apical meristems produce cytokinin
    if (node.meristem && node.meristem->active &&
        node.meristem->type == MeristemType::ROOT_APICAL) {
        node.cytokinin += genome.cytokinin_production_rate;
    }

    // Transport upward: send fraction to parent
    if (node.parent) {
        float flow = node.cytokinin * genome.cytokinin_transport_rate;
        node.parent->cytokinin += flow;
        node.cytokinin -= flow;
    }
}

static void cytokinin_distribute(Node& node, const Genome& genome) {
    // Transport downward: send fraction to children, split evenly
    if (!node.children.empty()) {
        float flow_total = node.cytokinin * genome.cytokinin_transport_rate;
        float flow_per_child = flow_total / static_cast<float>(node.children.size());
        for (Node* child : node.children) {
            child->cytokinin += flow_per_child;
        }
        node.cytokinin -= flow_total;
    }

    // Decay
    node.cytokinin *= (1.0f - genome.cytokinin_decay_rate);

    for (Node* child : node.children) {
        cytokinin_distribute(*child, genome);
    }
}

void transport_cytokinin(Plant& plant) {
    cytokinin_collect(*plant.seed_mut(), plant.genome());
    cytokinin_distribute(*plant.seed_mut(), plant.genome());
}

} // namespace botany
