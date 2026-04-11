#pragma once

#include "engine/node/node.h"

namespace botany {

class RootApicalNode : public Node {
public:
    RootApicalNode(uint32_t id, glm::vec3 position, float radius);
    void tick(Plant& plant, const WorldParams& world) override;
    float maintenance_cost(const Genome& g) const override;

    float target_internode_length = 0.0f;
    uint32_t ticks_since_last_node = 0;

private:
    glm::vec3 apply_gravitropism(const glm::vec3& dir, const Genome& g) const;
    bool grow(const Genome& g, const WorldParams& world, const glm::vec3& dir);
    void split_internode(Plant& plant, const Genome& g, const glm::vec3& dir);
};

} // namespace botany
