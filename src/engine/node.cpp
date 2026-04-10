#include "engine/node.h"

namespace botany {

Node::Node(uint32_t id, NodeType type, glm::vec3 position, float radius)
    : id(id)
    , parent(nullptr)
    , position(position)
    , radius(radius)
    , type(type)
    , age(0)
    , auxin(0.0f)
    , cytokinin(0.0f)
    , meristem(nullptr)
{}

void Node::add_child(Node* child) {
    children.push_back(child);
    child->parent = this;
}

} // namespace botany
