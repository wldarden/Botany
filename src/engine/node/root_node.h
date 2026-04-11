#pragma once

#include "engine/node/node.h"

namespace botany {

class RootNode : public Node {
public:
    RootNode(uint32_t id, glm::vec3 position, float radius);
    void tick(Plant& plant, const WorldParams& world) override;
};

} // namespace botany
