#include "engine/hormone.h"
#include "engine/plant.h"
#include "engine/node.h"

namespace botany {

// Phase 1 (post-order): produce auxin at shoot tips and flow basipetally
// toward the seed. Junction nodes accumulate auxin from all child branches.
static void auxin_collect(Node& node, const Genome& genome) {
    for (Node* child : node.children) {
        auxin_collect(*child, genome);
    }

    // Production: only active shoot apical meristems produce auxin
    if (node.meristem && node.meristem->active &&
        node.meristem->type() == MeristemType::APICAL) {
        node.auxin += genome.auxin_production_rate;
    }

    // Basipetal transport: send fraction toward parent (seed-ward)
    if (node.parent) {
        float flow = node.auxin * genome.auxin_transport_rate;
        node.parent->auxin += flow;
        node.auxin -= flow;
    }
}

// Phase 2 (pre-order): a small fraction of accumulated auxin spills back
// into child branches — the "traffic jam" at junctions where multiple
// branches feed auxin into one node pushes some back up into branches.
static void auxin_spillback(Node& node, const Genome& genome) {
    if (!node.children.empty()) {
        float flow_total = node.auxin * genome.auxin_spillback_rate;
        float flow_per_child = flow_total / static_cast<float>(node.children.size());
        for (Node* child : node.children) {
            child->auxin += flow_per_child;
        }
        node.auxin -= flow_total;
    }

    // Decay
    node.auxin *= (1.0f - genome.auxin_decay_rate);

    for (Node* child : node.children) {
        auxin_spillback(*child, genome);
    }
}

static void reset_auxin(Node& node) {
    node.auxin = 0.0f;
    for (Node* child : node.children) {
        reset_auxin(*child);
    }
}

void transport_auxin(Plant& plant) {
    reset_auxin(*plant.seed_mut());
    auxin_collect(*plant.seed_mut(), plant.genome());
    auxin_spillback(*plant.seed_mut(), plant.genome());
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
        node.meristem->type() == MeristemType::ROOT_APICAL) {
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

static void reset_cytokinin(Node& node) {
    node.cytokinin = 0.0f;
    for (Node* child : node.children) {
        reset_cytokinin(*child);
    }
}

void transport_cytokinin(Plant& plant) {
    reset_cytokinin(*plant.seed_mut());
    cytokinin_collect(*plant.seed_mut(), plant.genome());
    cytokinin_distribute(*plant.seed_mut(), plant.genome());
}

} // namespace botany
