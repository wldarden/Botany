#pragma once

#include <string>
#include <cstdint>
#include <unordered_map>
#include "renderer/shader.h"
#include "renderer/camera.h"

struct GLFWwindow;

namespace botany {

class Plant;
struct TickSnapshot;

class Renderer {
public:
    bool init(int width, int height, const std::string& shader_dir);
    void shutdown();

    GLFWwindow* window() const { return window_; }
    OrbitCamera& camera() { return camera_; }

    void begin_frame();
    void draw_plant(const Plant& plant);
    void draw_snapshot(const TickSnapshot& snapshot);
    void draw_ground();
    void end_frame();

private:
    GLFWwindow* window_ = nullptr;
    Shader shader_;
    OrbitCamera camera_;
    int width_;
    int height_;

    uint32_t ground_vao_ = 0;
    uint32_t ground_vbo_ = 0;

    void setup_ground();
    void draw_cylinder(glm::vec3 start, glm::vec3 end,
                       float r_start, float r_end,
                       glm::vec3 color, int segments = 8);
    void draw_leaf(glm::vec3 position, glm::vec3 direction, float size);
};

} // namespace botany
