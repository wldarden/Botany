#include "engine/node/meristem_node.h"

namespace botany {

MeristemNode::MeristemNode(uint32_t id, NodeType type, glm::vec3 position, float radius, bool active)
    : Node(id, type, position, radius)
    , active(active)
{}

bool MeristemNode::is_tip() const {
    return type == NodeType::SHOOT_APICAL || type == NodeType::ROOT_APICAL;
}

void MeristemNode::tick(Plant& plant, const WorldParams& world) {
    Node::tick(plant, world);
    ticks_since_last_node++;
}

} // namespace botany
