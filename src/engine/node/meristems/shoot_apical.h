#pragma once

#include "engine/node/node.h"

namespace botany {

class ShootApicalNode : public Node {
public:
    // --- Node overrides ---
    ShootApicalNode(uint32_t id, glm::vec3 position, float radius);
    void tick(Plant& plant, const WorldParams& world) override;
    void produce(Plant& plant, const WorldParams& world) override;
    float maintenance_cost(const WorldParams& world) const override;
    void grow(Plant& plant, const WorldParams& world) override;

    // --- Meristem state ---
    glm::vec3 growth_dir = glm::vec3(0.0f);
    glm::vec3 set_point_dir = glm::vec3(0.0f, 1.0f, 0.0f); // preferred angle (vertical for trunk, branch angle for branches)
    uint32_t phyllotaxis_index = 0;
    uint32_t ticks_since_last_node = 0;

private:
    // --- Auxin production (environmentally modulated) ---
    float produce_auxin(const Plant& plant) const;
    float estimate_local_light() const;

    // --- Chain growth ---
    void roll_direction(const Genome& g, const WorldParams& world);
    void grow_tip(Plant& plant, const Genome& g, const WorldParams& world);
    void spawn_internode(Plant& plant, const Genome& g);
    void spawn_axillary(Plant& plant, Node* internode, const Genome& g, const glm::vec3& lateral_offset);
    void spawn_leaf(Plant& plant, Node* internode, const Genome& g, const glm::vec3& lateral_offset);
};

} // namespace botany
