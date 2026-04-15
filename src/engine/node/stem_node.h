#pragma once

#include "engine/node/node.h"

namespace botany {

class StemNode : public Node {
public:
    // --- Node overrides ---
    StemNode(uint32_t id, glm::vec3 position, float radius);
    float maintenance_cost(const WorldParams& world) const override;
    void tissue_tick(Plant& plant, const WorldParams& world) override;

private:
    // --- Stem growth ---
    void thicken(const Genome& g, const WorldParams& world);
    void elongate(const Genome& g, const WorldParams& world);
};

} // namespace botany
