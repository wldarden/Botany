#pragma once

#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>

namespace botany {

class OrbitCamera {
public:
    OrbitCamera();

    glm::mat4 view_matrix() const;
    glm::mat4 projection_matrix(float aspect) const;

    void rotate(float dx, float dy);
    void zoom(float delta);
    void set_target(glm::vec3 target) { target_ = target; }

    glm::vec3 target() const { return target_; }

private:
    glm::vec3 target_;
    float distance_;
    float yaw_;
    float pitch_;
    float fov_;
    float near_;
    float far_;
};

} // namespace botany
