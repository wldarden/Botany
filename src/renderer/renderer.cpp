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
#include "engine/node/tissues/apical.h"
#include "engine/node/tissues/root_apical.h"
#include "engine/node/tissues/leaf.h"
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

    // Report actual context version — macOS may give us less than we asked for.
    const char* gl_version = reinterpret_cast<const char*>(glGetString(GL_VERSION));
    const char* gl_renderer = reinterpret_cast<const char*>(glGetString(GL_RENDERER));
    std::cout << "OpenGL: " << (gl_version ? gl_version : "unknown")
              << "  Renderer: " << (gl_renderer ? gl_renderer : "unknown") << std::endl;

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.53f, 0.81f, 0.92f, 1.0f);

    if (!shader_.load(shader_dir + "/plant.vert", shader_dir + "/plant.frag")) {
        return false;
    }
    // Non-fatal — falls back to plain ground rendering without shadows.
    ground_shader_.load(shader_dir + "/ground.vert", shader_dir + "/ground.frag");
    ground_shadow_shader_.load(shader_dir + "/ground_shadow.vert", shader_dir + "/ground_shadow.frag");

    setup_ground();

    if (!light_system_.init(shader_dir)) {
        std::cerr << "Warning: LightSystem GPU shadows unavailable" << std::endl;
    }

    return true;
}

void Renderer::shutdown() {
    light_system_.shutdown();
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
    shader_.set_vec3("uLightDir", light_system_.sun_direction);  // shader does -uLightDir internally; sun_direction already points toward surface
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
    // If coloring by chemical, find per-type max values for normalization.
    // Each tissue category gets its own color scale so leaves don't wash out
    // relative to sugar-rich stems, etc.
    float max_leaf = 0.0f, max_root = 0.0f, max_shoot = 0.0f;
    if (chemical_accessor_) {
        plant.for_each_node([&](const Node& node) {
            float v = chemical_accessor_(node);
            if (node.type == NodeType::LEAF) {
                if (v > max_leaf) max_leaf = v;
            } else if (node.type == NodeType::ROOT || node.type == NodeType::ROOT_APICAL) {
                if (v > max_root) max_root = v;
            } else {
                if (v > max_shoot) max_shoot = v;
            }
        });
        if (max_leaf  < 1e-6f) max_leaf  = 1.0f;
        if (max_root  < 1e-6f) max_root  = 1.0f;
        if (max_shoot < 1e-6f) max_shoot = 1.0f;
    }
    auto max_for = [&](const Node& n) -> float {
        if (n.type == NodeType::LEAF) return max_leaf;
        if (n.type == NodeType::ROOT || n.type == NodeType::ROOT_APICAL) return max_root;
        return max_shoot;
    };

    plant.for_each_node([&](const Node& node) {
        if (!node.parent) return;

        // Skip drawing dormant meristems — there are many and they're invisible
        if (auto* ap = node.as_apical()) { if (!ap->active) return; }
        if (auto* ra = node.as_root_apical()) { if (!ra->active) return; }

        if (auto* leaf = node.as_leaf()) {
            glm::vec3 outward = glm::normalize(node.position - node.parent->position);
            glm::vec3 leaf_color;
            const glm::vec3* vtx_ptr = nullptr;
            glm::vec3 vtx_colors[4];

            if (chemical_accessor_) {
                float v = chemical_accessor_(node);
                leaf_color = heatmap(v / max_for(node));
                // Chemical overlay → single color, no per-vertex shading
            } else {
                glm::vec3 sun_color(0.2f, 0.6f, 0.15f);
                glm::vec3 shade_color(0.08f, 0.25f, 0.06f);
                glm::vec3 dead_color(0.4f, 0.3f, 0.1f);
                float starvation_t = (node.starvation_ticks > 0)
                    ? glm::clamp(static_cast<float>(node.starvation_ticks) / 50.0f, 0.0f, 1.0f)
                    : 0.0f;
                float senesce_t = (leaf->senescence_ticks > 0)
                    ? glm::clamp(static_cast<float>(leaf->senescence_ticks) / 48.0f, 0.0f, 1.0f)
                    : 0.0f;
                glm::vec3 senesce_color = glm::mix(glm::vec3(0.8f, 0.7f, 0.1f),
                                                   glm::vec3(0.4f, 0.25f, 0.05f), senesce_t);

                // Per-corner shading using GPU sample_exposure[0..3] (p0,p1,p2,p3).
                // sample_exposure defaults to 1.0, so this is safe even before first GPU update.
                for (int i = 0; i < 4; i++) {
                    glm::vec3 c = glm::mix(shade_color, sun_color, leaf->sample_exposure[i]);
                    if (starvation_t > 0.0f) c = glm::mix(c, dead_color, starvation_t);
                    if (senesce_t   > 0.0f) c = glm::mix(c, senesce_color, senesce_t);
                    vtx_colors[i] = c * color_tint_;
                }
                vtx_ptr = vtx_colors;
                leaf_color = vtx_colors[0];  // fallback for petiole color
            }

            // Senescence overlay for chemical-accessor case (uniform across leaf)
            if (chemical_accessor_ && leaf->senescence_ticks > 0) {
                float progress = glm::clamp(static_cast<float>(leaf->senescence_ticks) / 48.0f, 0.0f, 1.0f);
                glm::vec3 senesce_color = glm::mix(glm::vec3(0.8f, 0.7f, 0.1f),
                                                   glm::vec3(0.4f, 0.25f, 0.05f), progress);
                leaf_color = glm::mix(leaf_color, senesce_color, progress);
            }

            // Draw petiole (thin stalk from stem to leaf blade)
            float petiole_radius = 0.005f;  // 0.5mm — thin stalk
            glm::vec3 petiole_color = glm::vec3(0.3f, 0.45f, 0.12f) * color_tint_;  // green-brown
            draw_cylinder(node.parent->position, node.position, petiole_radius, petiole_radius, petiole_color);
            draw_leaf(node.position, outward, leaf->facing, leaf->leaf_size,
                      leaf_color * (chemical_accessor_ ? color_tint_ : 1.0f), vtx_ptr);
            return;
        }

        glm::vec3 color;
        if (chemical_accessor_) {
            float v = chemical_accessor_(node);
            color = heatmap(v / max_for(node));
        } else if (color_by_type_) {
            if (node.type == NodeType::APICAL) {
                color = glm::vec3(1.0f, 0.2f, 0.2f);   // red
            } else if (node.type == NodeType::ROOT_APICAL) {
                color = glm::vec3(0.2f, 0.4f, 1.0f);   // blue
            } else if (node.type == NodeType::STEM) {
                glm::vec3 green(0.2f, 0.8f, 0.2f);
                glm::vec3 brown(0.5f, 0.35f, 0.15f);
                float threshold = plant.genome().stem_green_radius_threshold;
                float maturity = glm::clamp(node.radius / threshold, 0.0f, 1.0f);
                color = glm::mix(green, brown, maturity);
            } else if (node.type == NodeType::ROOT) {
                color = glm::vec3(0.8f, 0.6f, 0.2f);   // yellow-brown
            } else {
                color = glm::vec3(0.2f, 0.8f, 0.2f);   // green (leaf drawn separately, but fallback)
            }
        } else if (node.type == NodeType::STEM) {
            // Young stems: green (chlorenchyma visible), mature: brown (bark covered).
            // Transition driven by genome threshold — same value that gates photosynthesis.
            glm::vec3 green(0.2f, 0.5f, 0.15f);
            glm::vec3 brown(0.35f, 0.2f, 0.1f);
            float threshold = plant.genome().stem_green_radius_threshold;
            float maturity = glm::clamp(node.radius / threshold, 0.0f, 1.0f);
            color = glm::mix(green, brown, maturity);
        } else if (node.type == NodeType::ROOT || node.type == NodeType::ROOT_APICAL) {
            // Young roots: white/tan (meristematic zone), mature: darker earth-brown.
            glm::vec3 young(0.90f, 0.82f, 0.68f);   // pale tan — fresh root tip
            glm::vec3 mature(0.55f, 0.38f, 0.22f);  // darker brown — suberized root
            float maturity = glm::clamp(node.radius / 0.05f, 0.0f, 1.0f);
            color = glm::mix(young, mature, maturity);
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

        // Apical meristems: hemisphere facing growth direction instead of cylinder.
        if (node.type == NodeType::APICAL || node.type == NodeType::ROOT_APICAL) {
            glm::vec3 grow_dir = (glm::length(node.offset) > 1e-4f)
                ? glm::normalize(node.offset)
                : (node.type == NodeType::APICAL ? glm::vec3(0.0f, 1.0f, 0.0f)
                                                 : glm::vec3(0.0f, -1.0f, 0.0f));
            draw_hemisphere(node.position, grow_dir, node.radius, color * color_tint_);
            return;
        }

        draw_cylinder(node.parent->position, node.position,
                      node.parent->radius, node.radius, color * color_tint_);
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
            glm::vec3 outward = glm::normalize(ns.position - parent.position);
            draw_leaf(ns.position, outward, ns.facing, ns.leaf_size);
            continue;
        }

        glm::vec3 color;
        if (ns.type == NodeType::STEM) {
            glm::vec3 green(0.2f, 0.5f, 0.15f);
            glm::vec3 brown(0.35f, 0.2f, 0.1f);
            float maturity = glm::clamp(ns.radius / 0.04f, 0.0f, 1.0f);  // default threshold
            color = glm::mix(green, brown, maturity);
        } else if (ns.type == NodeType::ROOT || ns.type == NodeType::ROOT_APICAL) {
            glm::vec3 young(0.90f, 0.82f, 0.68f);
            glm::vec3 mature(0.55f, 0.38f, 0.22f);
            float maturity = glm::clamp(ns.radius / 0.05f, 0.0f, 1.0f);
            color = glm::mix(young, mature, maturity);
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

void Renderer::draw_hemisphere(glm::vec3 center, glm::vec3 dir, float r, glm::vec3 color) {
    const int LON = 8;
    const int LAT = 4;
    const float half_pi = 1.5707963f;
    const float two_pi  = 6.2831853f;

    // Build an orthonormal frame: N = dome axis (outward), T and B fill the disc.
    glm::vec3 N = glm::normalize(dir);
    glm::vec3 T = (std::abs(N.y) < 0.9f)
        ? glm::normalize(glm::cross(N, glm::vec3(0.0f, 1.0f, 0.0f)))
        : glm::normalize(glm::cross(N, glm::vec3(1.0f, 0.0f, 0.0f)));
    glm::vec3 B = glm::cross(N, T);

    // Returns world position + outward normal for grid cell (lat_i, lon_i).
    // lat_i=0 → dome tip (N direction); lat_i=LAT → equator.
    // theta = polar angle from tip (0 → π/2).
    auto make = [&](int lat_i, int lon_i, glm::vec3& out_p, glm::vec3& out_n) {
        float theta = static_cast<float>(lat_i) / LAT * half_pi;
        float phi   = static_cast<float>(lon_i) / LON * two_pi;
        float x = std::sin(theta) * std::cos(phi);
        float y = std::sin(theta) * std::sin(phi);
        float z = std::cos(theta);
        out_n = T * x + B * y + N * z;  // outward unit normal
        out_p = center + out_n * r;
    };

    std::vector<float> verts;
    verts.reserve(LON * LAT * 6 * 9);

    auto push = [&](int lat_i, int lon_i) {
        glm::vec3 p, n;
        make(lat_i, lon_i, p, n);
        verts.push_back(p.x); verts.push_back(p.y); verts.push_back(p.z);
        verts.push_back(n.x); verts.push_back(n.y); verts.push_back(n.z);
        verts.push_back(color.x); verts.push_back(color.y); verts.push_back(color.z);
    };

    // Quad strip from tip (j=0) to equator (j=LAT).
    // Degenerate triangles at j=0 (all tip vertices are the same point) are discarded by GPU.
    for (int j = 0; j < LAT; j++) {
        for (int i = 0; i < LON; i++) {
            int i1 = (i + 1) % LON;
            push(j,   i ); push(j+1, i ); push(j+1, i1);
            push(j,   i ); push(j+1, i1); push(j,   i1);
        }
    }

    uint32_t vao, vbo;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(verts.size() * sizeof(float)),
                 verts.data(), GL_STREAM_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);

    glDrawArrays(GL_TRIANGLES, 0, static_cast<int>(verts.size() / 9));

    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
}

void Renderer::draw_leaf(glm::vec3 position, glm::vec3 outward, glm::vec3 facing, float size,
                          glm::vec3 color, const glm::vec3* vertex_colors) {
    // Project outward onto the leaf blade plane (perpendicular to facing)
    glm::vec3 proj = outward - facing * glm::dot(outward, facing);
    float proj_len = glm::length(proj);
    if (proj_len < 1e-4f) {
        // facing is parallel to outward — pick any direction in the blade plane
        if (std::abs(facing.y) < 0.9f) {
            proj = glm::normalize(glm::cross(facing, glm::vec3(0.0f, 1.0f, 0.0f)));
        } else {
            proj = glm::normalize(glm::cross(facing, glm::vec3(1.0f, 0.0f, 0.0f)));
        }
    } else {
        proj /= proj_len;
    }
    glm::vec3 width = glm::normalize(glm::cross(facing, proj));

    // Diamond: rotated 45 degrees, stem-side corner at node position (petiole attachment)
    float half = size * 0.70711f; // size * sqrt(2)/2
    glm::vec3 base = position;                        // petiole tip = stem-side corner
    glm::vec3 p0 = base;                              // stem-side corner (attachment point)
    glm::vec3 p1 = base + proj * half + width * half; // right corner
    glm::vec3 p2 = base + proj * (half * 2.0f);       // tip corner
    glm::vec3 p3 = base + proj * half - width * half; // left corner

    // Droop side corners for a natural cupped leaf shape
    glm::vec3 down(0.0f, -1.0f, 0.0f);
    float droop = size * 0.15f;
    p1 += down * droop;
    p3 += down * droop;

    // Per-triangle normals (halves are no longer coplanar)
    glm::vec3 n0 = glm::normalize(glm::cross(p1 - p0, p2 - p0));
    glm::vec3 n1 = glm::normalize(glm::cross(p2 - p0, p3 - p0));
    if (glm::dot(n0, facing) < 0.0f) n0 = -n0;
    if (glm::dot(n1, facing) < 0.0f) n1 = -n1;

    // Per-vertex colors: use vertex_colors[0..3] for p0,p1,p2,p3 if provided,
    // otherwise fall back to the single uniform color.
    glm::vec3 c0 = vertex_colors ? vertex_colors[0] : color;
    glm::vec3 c1 = vertex_colors ? vertex_colors[1] : color;
    glm::vec3 c2 = vertex_colors ? vertex_colors[2] : color;
    glm::vec3 c3 = vertex_colors ? vertex_colors[3] : color;

    float vertices[] = {
        p0.x, p0.y, p0.z, n0.x, n0.y, n0.z, c0.x, c0.y, c0.z,
        p1.x, p1.y, p1.z, n0.x, n0.y, n0.z, c1.x, c1.y, c1.z,
        p2.x, p2.y, p2.z, n0.x, n0.y, n0.z, c2.x, c2.y, c2.z,
        p0.x, p0.y, p0.z, n1.x, n1.y, n1.z, c0.x, c0.y, c0.z,
        p2.x, p2.y, p2.z, n1.x, n1.y, n1.z, c2.x, c2.y, c2.z,
        p3.x, p3.y, p3.z, n1.x, n1.y, n1.z, c3.x, c3.y, c3.z,
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

void Renderer::draw_highlight(const Node& node) {
    if (!node.parent) return;

    // Draw a bright thickened cylinder over the hovered segment
    glm::vec3 highlight_color(1.0f, 1.0f, 0.0f); // bright yellow
    float highlight_radius = std::max(node.radius * 2.0f, 0.02f);

    if (auto* leaf = node.as_leaf()) {
        if (leaf->leaf_size > 0.0f) {
            glm::vec3 outward = glm::normalize(node.position - node.parent->position);
            // Use GL_LEQUAL so the highlight passes the depth test when drawn at the same
            // depth as the leaf already in the buffer (GL_LESS would silently discard it).
            glDepthFunc(GL_LEQUAL);
            draw_leaf(node.position, outward, leaf->facing, leaf->leaf_size, highlight_color);
            glDepthFunc(GL_LESS);
            return;
        }
    }

    draw_cylinder(node.parent->position, node.position,
                  highlight_radius, highlight_radius,
                  highlight_color, 12);
}

void Renderer::draw_shadow_map(const ShadowMapViz& shadow_map) {
    if (shadow_map.cells.empty()) return;

    float half = shadow_map.cell_size * 0.5f;
    float y = 0.001f; // slightly above ground to avoid z-fighting
    glm::vec3 normal(0.0f, 1.0f, 0.0f);

    std::vector<float> vertices;
    vertices.reserve(shadow_map.cells.size() * 6 * 9); // 6 verts per quad, 9 floats per vert

    for (const auto& cell : shadow_map.cells) {
        float shade = std::min(cell.coverage, 1.0f);
        // Yellow = full sun, dark purple = fully shaded
        glm::vec3 sun_color(1.0f, 0.95f, 0.4f);
        glm::vec3 shade_color(0.15f, 0.05f, 0.2f);
        glm::vec3 color = glm::mix(sun_color, shade_color, shade);

        float x0 = cell.x - half, x1 = cell.x + half;
        float z0 = cell.z - half, z1 = cell.z + half;

        auto push = [&](float x, float z) {
            vertices.push_back(x); vertices.push_back(y); vertices.push_back(z);
            vertices.push_back(normal.x); vertices.push_back(normal.y); vertices.push_back(normal.z);
            vertices.push_back(color.x); vertices.push_back(color.y); vertices.push_back(color.z);
        };

        push(x0, z0); push(x1, z0); push(x1, z1);
        push(x0, z0); push(x1, z1); push(x0, z1);
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
    if (ground_shader_.is_loaded() && light_system_.is_initialized()) {
        ground_shader_.use();
        float aspect = static_cast<float>(width_) / static_cast<float>(height_);
        ground_shader_.set_mat4("uView",       camera_.view_matrix());
        ground_shader_.set_mat4("uProjection", camera_.projection_matrix(aspect));
        ground_shader_.set_mat4("uModel",      glm::mat4(1.0f));
        ground_shader_.set_vec3("uLightDir",   light_system_.sun_direction);
        ground_shader_.set_vec3("uAmbient",    glm::vec3(0.3f, 0.3f, 0.3f));
        ground_shader_.set_mat4("u_light_pv",  light_system_.light_pv());
        ground_shader_.set_int("u_num_slices", LightSystem::NUM_SLICES);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D_ARRAY, light_system_.slice_tex());
        ground_shader_.set_int("u_slice_array", 0);
    } else {
        shader_.use();
    }

    glBindVertexArray(ground_vao_);
    glDrawArrays(GL_TRIANGLES, 0, 6);

    // Grid always uses the plant shader (thin lines, shadow not worth it).
    shader_.use();
    draw_grid();
}

void Renderer::draw_ground_shadow() {
    if (ground_shadow_shader_.is_loaded() && light_system_.is_initialized()) {
        ground_shadow_shader_.use();
        float aspect = static_cast<float>(width_) / static_cast<float>(height_);
        ground_shadow_shader_.set_mat4("uView",        camera_.view_matrix());
        ground_shadow_shader_.set_mat4("uProjection",  camera_.projection_matrix(aspect));
        ground_shadow_shader_.set_mat4("uModel",       glm::mat4(1.0f));
        ground_shadow_shader_.set_mat4("u_light_view", light_system_.light_view());
        ground_shadow_shader_.set_mat4("u_light_proj", light_system_.light_proj());
        ground_shadow_shader_.set_vec3("u_sun_dir",    glm::normalize(light_system_.sun_direction));
        ground_shadow_shader_.set_float("u_min_y",     light_system_.min_depth());
        ground_shadow_shader_.set_float("u_max_y",     light_system_.max_depth());
        ground_shadow_shader_.set_int("u_num_slices",  LightSystem::NUM_SLICES);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D_ARRAY, light_system_.slice_tex());
        ground_shadow_shader_.set_int("u_slice_array", 0);
    } else {
        // Fallback: plain ground shader when shadow system isn't ready.
        shader_.use();
        float aspect = static_cast<float>(width_) / static_cast<float>(height_);
        shader_.set_mat4("uView",       camera_.view_matrix());
        shader_.set_mat4("uProjection", camera_.projection_matrix(aspect));
        shader_.set_mat4("uModel",      glm::mat4(1.0f));
    }

    glBindVertexArray(ground_vao_);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glBindTexture(GL_TEXTURE_2D_ARRAY, 0);

    // Grid always uses the plant shader.
    shader_.use();
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
