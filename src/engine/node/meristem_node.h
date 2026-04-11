#pragma once

#include "engine/node/node.h"

namespace botany {

class MeristemNode : public Node {
public:
    bool active;
    uint32_t ticks_since_last_node = 0;

    bool is_tip() const;

    void tick(Plant& plant, const WorldParams& world) override;

protected:
    MeristemNode(uint32_t id, NodeType type, glm::vec3 position, float radius, bool active);
};

} // namespace botany
