#pragma once

#include "engine/node/node.h"

namespace botany {

class ApicalNode : public Node {
public:
    // --- Node overrides ---
    ApicalNode(uint32_t id, glm::vec3 position, float radius);
    float maintenance_cost(const WorldParams& world) const override;
    void update_tissue(Plant& plant, const WorldParams& world) override;

    // --- Meristem state ---
    bool active = true; // false = dormant bud, true = actively growing
    glm::vec3 growth_dir = glm::vec3(0.0f);
    glm::vec3 set_point_dir = glm::vec3(0.0f, 1.0f, 0.0f); // preferred angle (vertical for trunk, branch angle for branches)
    uint32_t phyllotaxis_index = 0;
    uint32_t ticks_since_last_node = 0;

private:
    // --- Dormant bud activation ---
    bool can_activate(const Genome& g, const WorldParams& world) const;
    void activate(const Genome& g, const WorldParams& world);

    // --- Tip production ---
    void produce_auxin(Plant& plant, const Genome& g, const WorldParams& world);
    void photosynthesize(Plant& plant, const Genome& g, const WorldParams& world);
    float estimate_local_light() const;

    // --- Chain growth ---
    void roll_direction(const Genome& g, const WorldParams& world);
    void elongate(Plant& plant, const Genome& g, const WorldParams& world);
    void check_spawn(Plant& plant, const Genome& g);
    void spawn_internode(Plant& plant, const Genome& g);
    void spawn_axillary(Plant& plant, Node* internode, const Genome& g, const glm::vec3& lateral_offset);
    void spawn_leaf(Plant& plant, Node* internode, const Genome& g, const glm::vec3& lateral_offset);
};

} // namespace botany
