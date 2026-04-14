#pragma once

#include "engine/node/node.h"

namespace botany {

class LeafNode : public Node {
public:
    // --- Leaf state ---
    float leaf_size = 0.0f;
    float light_exposure = 1.0f;
    uint32_t senescence_ticks = 0;
    uint32_t deficit_ticks = 0;  // consecutive ticks where maintenance > production
    glm::vec3 facing = glm::vec3(0.0f, 1.0f, 0.0f);  // blade orientation (independent of attachment point)

    // --- Node overrides ---
    LeafNode(uint32_t id, glm::vec3 position, float radius);
    float maintenance_cost(const WorldParams& world) const override;
    void produce(Plant& plant, const WorldParams& world) override;
    void grow(Plant& plant, const WorldParams& world) override;

private:
    // --- Leaf behavior ---
    void photosynthesize(const Genome& g, const WorldParams& world);
    void phototropism(const Genome& g, const WorldParams& world);
    void grow_size(const Genome& g, const WorldParams& world);
};

} // namespace botany
