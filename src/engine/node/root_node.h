#pragma once

#include "engine/node/node.h"

namespace botany {

class RootNode : public Node {
public:
    // --- Node overrides ---
    RootNode(uint32_t id, glm::vec3 position, float radius);
    float maintenance_cost(const WorldParams& world) const override;
    void update_tissue(Plant& plant, const WorldParams& world) override;

private:
    // --- Root growth ---
    void absorb_water(const Genome& g, const WorldParams& world);
    void thicken(const Genome& g, const WorldParams& world);
    void elongate(const Genome& g, const WorldParams& world);
};

} // namespace botany
