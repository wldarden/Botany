#include "renderer/camera.h"
#include <glm/gtc/matrix_transform.hpp>
#include <cmath>
#include <algorithm>

namespace botany {

OrbitCamera::OrbitCamera()
    : target_(0.0f, 0.0f, 0.0f)
    , distance_(5.0f)
    , yaw_(0.0f)
    , pitch_(0.3f)
    , fov_(45.0f)
    , near_(0.1f)
    , far_(100.0f)
{}

glm::mat4 OrbitCamera::view_matrix() const {
    float x = distance_ * std::cos(pitch_) * std::sin(yaw_);
    float y = distance_ * std::sin(pitch_);
    float z = distance_ * std::cos(pitch_) * std::cos(yaw_);
    glm::vec3 eye = target_ + glm::vec3(x, y, z);
    return glm::lookAt(eye, target_, glm::vec3(0.0f, 1.0f, 0.0f));
}

glm::mat4 OrbitCamera::projection_matrix(float aspect) const {
    return glm::perspective(glm::radians(fov_), aspect, near_, far_);
}

void OrbitCamera::rotate(float dx, float dy) {
    yaw_ += dx * 0.005f;
    pitch_ += dy * 0.005f;
    pitch_ = std::clamp(pitch_, -1.5f, 1.5f);
}

void OrbitCamera::zoom(float delta) {
    distance_ -= delta * 0.5f;
    distance_ = std::clamp(distance_, 0.5f, 50.0f);
}

} // namespace botany
