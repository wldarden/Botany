#pragma once

#include "engine/node/meristem_node.h"

namespace botany {

class ShootApicalNode : public MeristemNode {
public:
    ShootApicalNode(uint32_t id, glm::vec3 position, float radius);
    void tick(Plant& plant, const WorldParams& world) override;

    float target_internode_length = 0.0f;
    uint32_t phyllotaxis_index = 0;

private:
    bool grow(const struct Genome& g, const WorldParams& world, const glm::vec3& dir);
    void split_internode(Plant& plant, const struct Genome& g, const glm::vec3& dir);
};

} // namespace botany
