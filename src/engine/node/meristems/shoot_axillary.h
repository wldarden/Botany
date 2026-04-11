#pragma once

#include "engine/node/meristem_node.h"

namespace botany {

class ShootAxillaryNode : public MeristemNode {
public:
    ShootAxillaryNode(uint32_t id, glm::vec3 position, float radius);
    void tick(Plant& plant, const WorldParams& world) override;
};

} // namespace botany
