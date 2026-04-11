#include "renderer/renderer.h"
#include <glad/gl.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/geometric.hpp>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <iostream>
#include "engine/plant.h"
#include "engine/node/leaf_node.h"
#include "engine/node/node.h"
#include "serialization/serializer.h"

namespace botany {

bool Renderer::init(int width, int height, const std::string& shader_dir) {
    width_ = width;
    height_ = height;

    if (!glfwInit()) {
        std::cerr << "Failed to init GLFW" << std::endl;
        return false;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    window_ = glfwCreateWindow(width, height, "Botany", nullptr, nullptr);
    if (!window_) {
        std::cerr << "Failed to create window" << std::endl;
        glfwTerminate();
        return false;
    }
    glfwMakeContextCurrent(window_);

    int version = gladLoadGL(glfwGetProcAddress);
    if (!version) {
        std::cerr << "Failed to load OpenGL" << std::endl;
        return false;
    }

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.53f, 0.81f, 0.92f, 1.0f);

    if (!shader_.load(shader_dir + "/plant.vert", shader_dir + "/plant.frag")) {
        return false;
    }

    setup_ground();
    return true;
}

void Renderer::shutdown() {
    if (ground_vao_) { glDeleteVertexArrays(1, &ground_vao_); ground_vao_ = 0; }
    if (ground_vbo_) { glDeleteBuffers(1, &ground_vbo_); ground_vbo_ = 0; }
    if (window_) { glfwDestroyWindow(window_); window_ = nullptr; }
    glfwTerminate();
}

void Renderer::begin_frame() {
    glfwGetFramebufferSize(window_, &width_, &height_);
    glViewport(0, 0, width_, height_);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    shader_.use();
    float aspect = static_cast<float>(width_) / static_cast<float>(height_);
    shader_.set_mat4("uView", camera_.view_matrix());
    shader_.set_mat4("uProjection", camera_.projection_matrix(aspect));
    shader_.set_mat4("uModel", glm::mat4(1.0f));
    shader_.set_vec3("uLightDir", glm::vec3(-0.3f, -1.0f, -0.5f));
    shader_.set_vec3("uAmbient", glm::vec3(0.3f, 0.3f, 0.3f));
}

void Renderer::set_color_mode(ChemicalAccessor accessor) {
    chemical_accessor_ = std::move(accessor);
}

glm::vec3 Renderer::heatmap(float t) const {
    // blue (0) -> cyan (0.25) -> green (0.5) -> yellow (0.75) -> red (1)
    t = glm::clamp(t, 0.0f, 1.0f);
    if (t < 0.25f) {
        float s = t / 0.25f;
        return glm::vec3(0.0f, s, 1.0f);
    } else if (t < 0.5f) {
        float s = (t - 0.25f) / 0.25f;
        return glm::vec3(0.0f, 1.0f, 1.0f - s);
    } else if (t < 0.75f) {
        float s = (t - 0.5f) / 0.25f;
        return glm::vec3(s, 1.0f, 0.0f);
    } else {
        float s = (t - 0.75f) / 0.25f;
        return glm::vec3(1.0f, 1.0f - s, 0.0f);
    }
}

void Renderer::draw_plant(const Plant& plant) {
    // If coloring by chemical, find max value for normalization
    float max_val = 0.0f;
    if (chemical_accessor_) {
        plant.for_each_node([&](const Node& node) {
            float v = chemical_accessor_(node);
            if (v > max_val) max_val = v;
        });
        if (max_val < 1e-6f) max_val = 1.0f; // avoid div-by-zero
    }

    plant.for_each_node([&](const Node& node) {
        if (!node.parent) return;

        if (auto* leaf = node.as_leaf()) {
            glm::vec3 dir = glm::normalize(node.position - node.parent->position);
            glm::vec3 leaf_color;
            if (chemical_accessor_) {
                float v = chemical_accessor_(node);
                leaf_color = heatmap(v / max_val);
            } else {
                // Sunlit leaves are bright green, shaded leaves darken
                glm::vec3 sun_color(0.2f, 0.6f, 0.15f);
                glm::vec3 shade_color(0.08f, 0.25f, 0.06f);
                leaf_color = glm::mix(shade_color, sun_color, leaf->light_exposure);
                if (node.starvation_ticks > 0) {
                    float stress = static_cast<float>(node.starvation_ticks) / 50.0f;
                    stress = glm::clamp(stress, 0.0f, 1.0f);
                    glm::vec3 dead_color(0.4f, 0.3f, 0.1f);
                    leaf_color = glm::mix(leaf_color, dead_color, stress);
                }
            }
            // Senescence: green -> yellow -> brown (overrides other coloring)
            if (leaf->senescence_ticks > 0) {
                float progress = static_cast<float>(leaf->senescence_ticks) / 48.0f;
                progress = glm::clamp(progress, 0.0f, 1.0f);
                glm::vec3 yellow(0.8f, 0.7f, 0.1f);
                glm::vec3 brown(0.4f, 0.25f, 0.05f);
                glm::vec3 senesce_color = glm::mix(yellow, brown, progress);
                leaf_color = glm::mix(leaf_color, senesce_color, progress);
            }
            draw_leaf(node.position, dir, leaf->leaf_size, leaf_color);
            return;
        }

        glm::vec3 color;
        if (chemical_accessor_) {
            float v = chemical_accessor_(node);
            color = heatmap(v / max_val);
        } else if (color_by_type_) {
            if (node.type == NodeType::SHOOT_APICAL) {
                color = glm::vec3(1.0f, 0.2f, 0.2f);   // red
            } else if (node.type == NodeType::SHOOT_AXILLARY) {
                color = glm::vec3(1.0f, 0.6f, 0.0f);   // orange
            } else if (node.type == NodeType::ROOT_APICAL) {
                color = glm::vec3(0.2f, 0.4f, 1.0f);   // blue
            } else if (node.type == NodeType::ROOT_AXILLARY) {
                color = glm::vec3(0.6f, 0.2f, 0.8f);   // purple
            } else if (node.type == NodeType::STEM) {
                glm::vec3 green(0.2f, 0.8f, 0.2f);
                glm::vec3 brown(0.5f, 0.35f, 0.15f);
                float radius_mat = (node.radius - 0.05f) / 0.10f;
                float age_mat = static_cast<float>(node.age) / 360.0f; // 15 days
                float maturity = glm::clamp(std::max(radius_mat, age_mat), 0.0f, 1.0f);
                color = glm::mix(green, brown, maturity);
            } else if (node.type == NodeType::ROOT) {
                color = glm::vec3(0.8f, 0.6f, 0.2f);   // yellow-brown
            } else {
                color = glm::vec3(0.2f, 0.8f, 0.2f);   // green (leaf drawn separately, but fallback)
            }
        } else if (node.type == NodeType::STEM) {
            // Young stems are green, maturing to brown as they thicken or age
            glm::vec3 green(0.25f, 0.55f, 0.15f);
            glm::vec3 brown(0.45f, 0.3f, 0.15f);
            float radius_mat = (node.radius - 0.05f) / 0.10f;
            float age_mat = static_cast<float>(node.age) / 360.0f; // 15 days
            float maturity = glm::clamp(std::max(radius_mat, age_mat), 0.0f, 1.0f);
            color = glm::mix(green, brown, maturity);
        } else {
            color = glm::vec3(0.35f, 0.2f, 0.1f);
        }

        // Starvation visual: interpolate toward grey/brown
        if (node.starvation_ticks > 0) {
            float stress = static_cast<float>(node.starvation_ticks) / 50.0f;
            stress = glm::clamp(stress, 0.0f, 1.0f);
            glm::vec3 dead_color(0.3f, 0.25f, 0.2f);
            color = glm::mix(color, dead_color, stress);
        }

        draw_cylinder(node.parent->position, node.position,
                      node.parent->radius, node.radius, color);
    });
}

void Renderer::draw_snapshot(const TickSnapshot& snapshot) {
    std::unordered_map<uint32_t, size_t> id_to_idx;
    for (size_t i = 0; i < snapshot.nodes.size(); i++) {
        id_to_idx[snapshot.nodes[i].id] = i;
    }

    for (const auto& ns : snapshot.nodes) {
        if (ns.parent_id == UINT32_MAX) continue;

        auto it = id_to_idx.find(ns.parent_id);
        if (it == id_to_idx.end()) continue;
        const auto& parent = snapshot.nodes[it->second];

        if (ns.type == NodeType::LEAF) {
            glm::vec3 dir = glm::normalize(ns.position - parent.position);
            draw_leaf(ns.position, dir, ns.leaf_size);
            continue;
        }

        glm::vec3 color;
        if (ns.type == NodeType::STEM) {
            glm::vec3 green(0.25f, 0.55f, 0.15f);
            glm::vec3 brown(0.45f, 0.3f, 0.15f);
            float maturity = glm::clamp((ns.radius - 0.05f) / 0.10f, 0.0f, 1.0f);
            color = glm::mix(green, brown, maturity);
        } else {
            color = glm::vec3(0.35f, 0.2f, 0.1f);
        }

        draw_cylinder(parent.position, ns.position,
                      parent.radius, ns.radius, color);
    }
}

void Renderer::draw_cylinder(glm::vec3 start, glm::vec3 end,
                              float r_start, float r_end,
                              glm::vec3 color, int segments) {
    glm::vec3 axis = end - start;
    float height = glm::length(axis);
    if (height < 0.0001f) return;
    glm::vec3 axis_norm = axis / height;

    glm::vec3 perp1;
    if (std::abs(axis_norm.y) < 0.9f) {
        perp1 = glm::normalize(glm::cross(axis_norm, glm::vec3(0.0f, 1.0f, 0.0f)));
    } else {
        perp1 = glm::normalize(glm::cross(axis_norm, glm::vec3(1.0f, 0.0f, 0.0f)));
    }
    glm::vec3 perp2 = glm::cross(axis_norm, perp1);

    std::vector<float> vertices;
    for (int i = 0; i < segments; i++) {
        float a0 = (2.0f * 3.14159f * i) / segments;
        float a1 = (2.0f * 3.14159f * (i + 1)) / segments;
        float c0 = std::cos(a0), s0 = std::sin(a0);
        float c1 = std::cos(a1), s1 = std::sin(a1);

        glm::vec3 n0 = perp1 * c0 + perp2 * s0;
        glm::vec3 n1 = perp1 * c1 + perp2 * s1;

        glm::vec3 b0 = start + n0 * r_start;
        glm::vec3 b1 = start + n1 * r_start;
        glm::vec3 t0 = end + n0 * r_end;
        glm::vec3 t1 = end + n1 * r_end;

        auto push = [&](glm::vec3 pos, glm::vec3 norm) {
            vertices.push_back(pos.x); vertices.push_back(pos.y); vertices.push_back(pos.z);
            vertices.push_back(norm.x); vertices.push_back(norm.y); vertices.push_back(norm.z);
            vertices.push_back(color.x); vertices.push_back(color.y); vertices.push_back(color.z);
        };

        push(b0, n0); push(t0, n0); push(b1, n1);
        push(b1, n1); push(t0, n0); push(t1, n1);
    }

    uint32_t vao, vbo;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STREAM_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);

    glDrawArrays(GL_TRIANGLES, 0, static_cast<int>(vertices.size() / 9));

    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
}

void Renderer::draw_leaf(glm::vec3 position, glm::vec3 direction, float size, glm::vec3 color) {
    glm::vec3 perp;
    if (std::abs(direction.y) < 0.9f) {
        perp = glm::normalize(glm::cross(direction, glm::vec3(0.0f, 1.0f, 0.0f)));
    } else {
        perp = glm::normalize(glm::cross(direction, glm::vec3(1.0f, 0.0f, 0.0f)));
    }
    glm::vec3 up = glm::normalize(glm::cross(perp, direction));

    // Leaf attached at position corner, extending outward along direction and perp
    glm::vec3 p0 = position;
    glm::vec3 p1 = position + perp * size;
    glm::vec3 p2 = position + perp * size + up * size;
    glm::vec3 p3 = position + up * size;
    glm::vec3 normal = glm::normalize(glm::cross(p1 - p0, p3 - p0));

    float vertices[] = {
        p0.x, p0.y, p0.z, normal.x, normal.y, normal.z, color.x, color.y, color.z,
        p1.x, p1.y, p1.z, normal.x, normal.y, normal.z, color.x, color.y, color.z,
        p2.x, p2.y, p2.z, normal.x, normal.y, normal.z, color.x, color.y, color.z,
        p0.x, p0.y, p0.z, normal.x, normal.y, normal.z, color.x, color.y, color.z,
        p2.x, p2.y, p2.z, normal.x, normal.y, normal.z, color.x, color.y, color.z,
        p3.x, p3.y, p3.z, normal.x, normal.y, normal.z, color.x, color.y, color.z,
    };

    uint32_t vao, vbo;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STREAM_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);

    glDrawArrays(GL_TRIANGLES, 0, 6);

    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
}

void Renderer::setup_ground() {
    float size = 20.0f;
    glm::vec3 color(0.4f, 0.35f, 0.25f);
    glm::vec3 normal(0.0f, 1.0f, 0.0f);
    float verts[] = {
        -size, 0, -size, normal.x, normal.y, normal.z, color.x, color.y, color.z,
         size, 0, -size, normal.x, normal.y, normal.z, color.x, color.y, color.z,
         size, 0,  size, normal.x, normal.y, normal.z, color.x, color.y, color.z,
        -size, 0, -size, normal.x, normal.y, normal.z, color.x, color.y, color.z,
         size, 0,  size, normal.x, normal.y, normal.z, color.x, color.y, color.z,
        -size, 0,  size, normal.x, normal.y, normal.z, color.x, color.y, color.z,
    };

    glGenVertexArrays(1, &ground_vao_);
    glGenBuffers(1, &ground_vbo_);
    glBindVertexArray(ground_vao_);
    glBindBuffer(GL_ARRAY_BUFFER, ground_vbo_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(verts), verts, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);
}

void Renderer::draw_ground() {
    glBindVertexArray(ground_vao_);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    draw_grid();
}

void Renderer::draw_grid() {
    // Grid on the ground plane (y = 0.001 to avoid z-fighting).
    // Lines every 1 dm (10 cm). Thicker/brighter every 10 dm (1 m).
    float y = 0.001f;
    float extent = 10.0f;  // ±10 dm = ±1 m
    glm::vec3 minor_color(0.35f, 0.30f, 0.20f);
    glm::vec3 major_color(0.30f, 0.25f, 0.15f);
    float minor_r = 0.003f;
    float major_r = 0.006f;

    for (int i = static_cast<int>(-extent); i <= static_cast<int>(extent); i++) {
        float pos = static_cast<float>(i);
        bool major = (i % 10 == 0);
        float r = major ? major_r : minor_r;
        glm::vec3 c = major ? major_color : minor_color;

        // X-parallel line
        draw_cylinder(glm::vec3(-extent, y, pos), glm::vec3(extent, y, pos), r, r, c, 4);
        // Z-parallel line
        draw_cylinder(glm::vec3(pos, y, -extent), glm::vec3(pos, y, extent), r, r, c, 4);
    }
}

void Renderer::end_frame() {
    glfwSwapBuffers(window_);
    glfwPollEvents();
}

} // namespace botany
