#pragma once

#include "engine/node/node.h"

namespace botany {

class RootApicalNode : public Node {
public:
    // --- Node overrides ---
    RootApicalNode(uint32_t id, glm::vec3 position, float radius);
    float maintenance_cost(const WorldParams& world) const override;
    void update_tissue(Plant& plant, const WorldParams& world) override;

    // --- Meristem state ---
    bool active = true; // false = dormant bud, true = actively growing
    bool is_primary = false; // true = plant's primary root apex (never quiesces; can be re-promoted from lateral if original dies)
    glm::vec3 growth_dir = glm::vec3(0.0f);
    uint32_t ticks_since_last_node = 0;
    uint32_t internodes_spawned = 0; // count of internodes this RA has laid down (used by primary RA to suppress initial lateral buds)

private:
    // --- Dormant bud activation ---
    bool can_activate(const Genome& g, const WorldParams& world) const;
    void activate(const Genome& g, const WorldParams& world);

    // --- Water uptake ---
    void absorb_water(const Genome& g, const WorldParams& world);

    // --- Chain growth ---
    void roll_direction(const Genome& g);
    glm::vec3 apply_gravitropism(const glm::vec3& dir, const Genome& g) const;
    void elongate(const Genome& g, const WorldParams& world);
    void check_spawn(Plant& plant, const Genome& g);
    void spawn_internode(Plant& plant, const Genome& g);
    void spawn_axillary(Plant& plant, Node* internode, const Genome& g, const glm::vec3& lateral_offset);
};

} // namespace botany
