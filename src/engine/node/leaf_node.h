#pragma once

#include "engine/node/node.h"

namespace botany {

class LeafNode : public Node {
public:
    float leaf_size = 0.0f;
    float light_exposure = 1.0f;
    uint32_t senescence_ticks = 0;

    LeafNode(uint32_t id, glm::vec3 position, float radius);
    void tick(Plant& plant, const WorldParams& world) override;
    float maintenance_cost(const Genome& g) const override;

private:
    void photosynthesize(const Genome& g, const WorldParams& world);
    void phototropism(const Genome& g, const WorldParams& world);
    void grow(const Genome& g, const WorldParams& world);
};

} // namespace botany
