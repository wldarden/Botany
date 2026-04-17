#pragma once

#include "engine/node/node.h"

namespace botany {

class StemNode : public Node {
public:
    // --- Node overrides ---
    StemNode(uint32_t id, glm::vec3 position, float radius);
    float maintenance_cost(const WorldParams& world) const override;
    void update_tissue(Plant& plant, const WorldParams& world) override;
    void compute_growth_reserve(const Genome& g, const WorldParams& world) override;

private:
    // --- Corticular photosynthesis ---
    void photosynthesize(Plant& plant, const Genome& g, const WorldParams& world);

    // --- Stem growth ---
    void thicken(const Genome& g, const WorldParams& world);
    void elongate(const Genome& g, const WorldParams& world);
};

} // namespace botany
