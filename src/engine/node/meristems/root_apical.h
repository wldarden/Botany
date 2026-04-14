#pragma once

#include "engine/node/node.h"

namespace botany {

class RootApicalNode : public Node {
public:
    // --- Node overrides ---
    RootApicalNode(uint32_t id, glm::vec3 position, float radius);
    void tick(Plant& plant, const WorldParams& world) override;
    float maintenance_cost(const WorldParams& world) const override;
    void grow(Plant& plant, const WorldParams& world) override;

    // --- Meristem state ---
    glm::vec3 growth_dir = glm::vec3(0.0f);
    uint32_t ticks_since_last_node = 0;

private:
    // --- Chain growth ---
    void roll_direction(const Genome& g);
    glm::vec3 apply_gravitropism(const glm::vec3& dir, const Genome& g) const;
    void grow_tip(const Genome& g, const WorldParams& world);
    void spawn_internode(Plant& plant, const Genome& g);
    void spawn_axillary(Plant& plant, Node* internode, const Genome& g, const glm::vec3& lateral_offset);
};

} // namespace botany
