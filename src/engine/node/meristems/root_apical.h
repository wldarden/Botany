#pragma once

#include "engine/node/meristem_node.h"

namespace botany {

class RootApicalNode : public MeristemNode {
public:
    RootApicalNode(uint32_t id, glm::vec3 position, float radius);
    void tick(Plant& plant, const WorldParams& world) override;

    float target_internode_length = 0.0f;
};

} // namespace botany
