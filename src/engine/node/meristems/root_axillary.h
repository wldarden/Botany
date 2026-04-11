#pragma once

#include "engine/node/node.h"

namespace botany {

class RootAxillaryNode : public Node {
public:
    bool active = false;

    RootAxillaryNode(uint32_t id, glm::vec3 position, float radius);
    void tick(Plant& plant, const WorldParams& world) override;

private:
    bool can_activate(const Genome& g, const WorldParams& world) const;
    void activate(Plant& plant, const Genome& g, const WorldParams& world);
};

} // namespace botany
