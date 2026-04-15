#pragma once

#include <string>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include "renderer/shader.h"
#include "renderer/camera.h"
#include "renderer/light_system.h"
#include "engine/light.h"

struct GLFWwindow;

namespace botany {

class Plant;
struct Node;
struct TickSnapshot;

// Function that extracts the chemical value from a node
using ChemicalAccessor = std::function<float(const Node&)>;

class Renderer {
public:
    bool init(int width, int height, const std::string& shader_dir);
    void shutdown();

    GLFWwindow* window() const { return window_; }
    OrbitCamera& camera() { return camera_; }

    void set_color_mode(ChemicalAccessor accessor);
    void set_color_by_type(bool enabled) { color_by_type_ = enabled; }
    void set_color_tint(float tint) { color_tint_ = tint; }  // 0-1 multiplier on all draw colors

    LightSystem& light_system() { return light_system_; }

    void begin_frame();
    void draw_plant(const Plant& plant);
    void draw_snapshot(const TickSnapshot& snapshot);
    void draw_ground();
    void draw_shadow_map(const ShadowMapViz& shadow_map);
    void draw_highlight(const Node& node);
    void end_frame();

private:
    GLFWwindow* window_ = nullptr;
    Shader shader_;
    Shader ground_shader_;  // shadow-aware ground shader (optional — falls back to shader_)
    OrbitCamera camera_;
    LightSystem light_system_;
    int width_;
    int height_;

    ChemicalAccessor chemical_accessor_;
    bool color_by_type_ = false;
    float color_tint_ = 1.0f;

    uint32_t ground_vao_ = 0;
    uint32_t ground_vbo_ = 0;

    void setup_ground();
    void draw_grid();
    glm::vec3 heatmap(float t) const;
    void draw_cylinder(glm::vec3 start, glm::vec3 end,
                       float r_start, float r_end,
                       glm::vec3 color, int segments = 8);
    void draw_leaf(glm::vec3 position, glm::vec3 outward, glm::vec3 facing, float size,
                   glm::vec3 color = glm::vec3(0.2f, 0.6f, 0.15f),
                   const glm::vec3* vertex_colors = nullptr);  // if set: per-corner colors for p0,p1,p2,p3
};

} // namespace botany
