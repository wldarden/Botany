#pragma once

#include "engine/node/node.h"

namespace botany {

class LeafNode : public Node {
public:
    // --- Leaf state ---
    float leaf_size = 0.0f;
    float light_exposure = 1.0f;               // average of 5 GPU samples (used by sim)
    float sample_exposure[5] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f};  // per-corner GPU values (used by renderer)
    uint32_t senescence_ticks = 0;
    uint32_t deficit_ticks = 0;  // consecutive ticks where maintenance > production
    glm::vec3 facing = glm::vec3(0.0f, 1.0f, 0.0f);  // blade orientation (independent of attachment point)

    // --- Node overrides ---
    LeafNode(uint32_t id, glm::vec3 position, float radius);
    float maintenance_cost(const WorldParams& world) const override;
    void update_tissue(Plant& plant, const WorldParams& world) override;
    void compute_growth_reserve(const Genome& g, const WorldParams& world) override;

private:
    // --- Leaf behavior ---
    void produce_gibberellin(const Genome& g);
    float photosynthesize(Plant& plant, const Genome& g, const WorldParams& world);
    void transpire(const Genome& g, const WorldParams& world);
    void check_carbon_balance(const Genome& g, const WorldParams& world, float net_sugar);
    bool advance_senescence(Plant& plant, const Genome& g);
    void phototropism(const Genome& g, const WorldParams& world);
    void expand(const Genome& g, const WorldParams& world);
};

} // namespace botany
