#include "renderer/light_system.h"
#include "engine/plant.h"
#include "engine/node/node.h"
#include "engine/node/tissues/leaf.h"
#include <glad/gl.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/geometric.hpp>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace botany {

// ---------------------------------------------------------------------------
// init / shutdown
// ---------------------------------------------------------------------------

bool LightSystem::init(const std::string& shader_dir) {
    // Shader programs
    if (!shadow_shader_.load(shader_dir + "/shadow_collect.vert",
                             shader_dir + "/shadow_collect.frag")) {
        std::cerr << "[LightSystem] Failed to load shadow_collect shaders\n";
        return false;
    }
    if (!query_shader_.load(shader_dir + "/leaf_query.vert",
                            shader_dir + "/leaf_query.frag")) {
        std::cerr << "[LightSystem] Failed to load leaf_query shaders\n";
        return false;
    }
    if (!debug_shader_.load(shader_dir + "/debug_slice.vert",
                            shader_dir + "/debug_slice.frag")) {
        // Debug shader is optional — log but don't fail init.
        std::cerr << "[LightSystem] Note: debug_slice shaders not found (debug overlay unavailable)\n";
    }

    // --- Slice array texture (512×512×32, GL_R32F) ---
    glGenTextures(1, &slice_array_tex_);
    glBindTexture(GL_TEXTURE_2D_ARRAY, slice_array_tex_);
    glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_R32F,
                 SHADOW_RES, SHADOW_RES, NUM_SLICES,
                 0, GL_RED, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D_ARRAY, 0);

    // --- 4 MRT FBOs, each binding 8 consecutive layers ---
    glGenFramebuffers(4, mrt_fbos_);
    for (int batch = 0; batch < 4; batch++) {
        glBindFramebuffer(GL_FRAMEBUFFER, mrt_fbos_[batch]);
        GLenum draw_buffers[8];
        for (int s = 0; s < 8; s++) {
            int layer = batch * 8 + s;
            glFramebufferTextureLayer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0 + s,
                                      slice_array_tex_, 0, layer);
            draw_buffers[s] = GL_COLOR_ATTACHMENT0 + s;
        }
        glDrawBuffers(8, draw_buffers);
        GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            std::cerr << "[LightSystem] MRT FBO " << batch << " incomplete: " << status << "\n";
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            return false;
        }
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // --- Caster VAO/VBO (leaf quads + stem cylinders from sun POV) ---
    glGenVertexArrays(1, &caster_vao_);
    glGenBuffers(1, &caster_vbo_);
    glBindVertexArray(caster_vao_);
    glBindBuffer(GL_ARRAY_BUFFER, caster_vbo_);
    // layout: vec3 pos (loc 0), float opacity (loc 1)
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(CasterVertex),
                          reinterpret_cast<void*>(offsetof(CasterVertex, x)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(CasterVertex),
                          reinterpret_cast<void*>(offsetof(CasterVertex, opacity)));
    glBindVertexArray(0);

    // --- Query VAO/VBO (5 sample points per leaf, rendered as GL_POINTS) ---
    glGenVertexArrays(1, &query_vao_);
    glGenBuffers(1, &query_vbo_);
    glBindVertexArray(query_vao_);
    glBindBuffer(GL_ARRAY_BUFFER, query_vbo_);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(LeafSamplePoint),
                          reinterpret_cast<void*>(offsetof(LeafSamplePoint, x)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(LeafSamplePoint),
                          reinterpret_cast<void*>(offsetof(LeafSamplePoint, output_index)));
    glBindVertexArray(0);

    // Output texture starts at 0 size; resized in resize_output_texture().
    glGenTextures(1, &output_tex_);
    glGenFramebuffers(1, &query_fbo_);

    // Debug fullscreen triangle
    glGenVertexArrays(1, &debug_vao_);
    glGenBuffers(1, &debug_vbo_);

    initialized_ = true;
    return true;
}

void LightSystem::shutdown() {
    if (!initialized_) return;
    glDeleteTextures(1, &slice_array_tex_);
    glDeleteFramebuffers(4, mrt_fbos_);
    glDeleteVertexArrays(1, &caster_vao_);
    glDeleteBuffers(1, &caster_vbo_);
    glDeleteVertexArrays(1, &query_vao_);
    glDeleteBuffers(1, &query_vbo_);
    glDeleteTextures(1, &output_tex_);
    glDeleteFramebuffers(1, &query_fbo_);
    glDeleteVertexArrays(1, &debug_vao_);
    glDeleteBuffers(1, &debug_vbo_);
    initialized_ = false;
}

// ---------------------------------------------------------------------------
// update — main entry point
// ---------------------------------------------------------------------------

void LightSystem::update(const std::vector<std::unique_ptr<Plant>>& plants) {
    if (!initialized_) return;

    collect_casters(plants);
    build_query_vbo(plants);

    if (n_leaves_ == 0) return;

    update_light_matrices();
    upload_caster_vbo();
    upload_query_vbo();
    resize_output_texture(n_leaves_ * SAMPLES_PER_LEAF);

    shadow_collect_pass();
    leaf_query_pass();
    readback_results();
}

// ---------------------------------------------------------------------------
// collect_casters — build caster_verts_ and track adaptive Y range
// ---------------------------------------------------------------------------

void LightSystem::collect_casters(const std::vector<std::unique_ptr<Plant>>& plants) {
    caster_verts_.clear();
    float lo = +1e9f, hi = -1e9f;

    for (const auto& plant : plants) {
        float opacity = plant->genome().leaf_opacity;
        plant->for_each_node([&](const Node& node) {
            // Skip underground nodes and zero-size nodes
            if (node.position.y < -0.01f) return;
            if (node.radius <= 0.0f && node.type != NodeType::LEAF) return;

            if (auto* leaf = node.as_leaf()) {
                if (leaf->leaf_size <= 0.0f) return;

                // Emit the 6 vertices of the leaf diamond (2 triangles),
                // matching Renderer::draw_leaf() geometry so shadows are accurate.
                glm::vec3 pts[5];
                compute_leaf_corners(node, leaf, pts);

                // Triangle 1: p0, p1, p2 (pts[0], pts[1], pts[2])
                // Triangle 2: p0, p2, p3 (pts[0], pts[2], pts[3])
                // pts[4] = center (not emitted as geometry, only used for query)
                auto push = [&](const glm::vec3& p) {
                    caster_verts_.push_back({p.x, p.y, p.z, opacity});
                    lo = std::min(lo, p.y);
                    hi = std::max(hi, p.y);
                };
                push(pts[0]); push(pts[1]); push(pts[2]);
                push(pts[0]); push(pts[2]); push(pts[3]);
            } else if (node.parent) {
                // Stem/root: emit a ribbon spanning the full segment parent→node
                // with width 2*radius. For a vertical sun the shadow footprint is
                // this ribbon's projection onto the XZ plane.
                float r = node.radius;
                glm::vec3 pa = node.parent->position;
                glm::vec3 pb = node.position;

                // Perpendicular direction in XZ (width axis of the ribbon).
                glm::vec2 seg2d(pb.x - pa.x, pb.z - pa.z);
                float seg_len_2d = glm::length(seg2d);
                glm::vec2 perp2d = (seg_len_2d > 1e-4f)
                    ? glm::vec2(-seg2d.y, seg2d.x) / seg_len_2d
                    : glm::vec2(1.0f, 0.0f);
                glm::vec3 perp(perp2d.x * r, 0.0f, perp2d.y * r);

                glm::vec3 c0 = pa - perp;
                glm::vec3 c1 = pa + perp;
                glm::vec3 c2 = pb + perp;
                glm::vec3 c3 = pb - perp;

                auto push = [&](const glm::vec3& p) {
                    caster_verts_.push_back({p.x, p.y, p.z, 1.0f});
                    lo = std::min(lo, p.y);
                    hi = std::max(hi, p.y);
                };
                push(c0); push(c1); push(c2);
                push(c0); push(c2); push(c3);
            }
        });
    }

    // Adaptive Y range: pad by one-slice worth on each side to include leaf corners
    // that may be slightly above or below the node positions.
    float range = hi - lo;
    float pad = (range > 0.0f) ? range / static_cast<float>(NUM_SLICES) : 1.0f;
    min_y_ = lo - pad;
    max_y_ = hi + pad;
}

// ---------------------------------------------------------------------------
// build_query_vbo — 5 sample points per leaf
// ---------------------------------------------------------------------------

void LightSystem::build_query_vbo(const std::vector<std::unique_ptr<Plant>>& plants) {
    query_verts_.clear();
    leaf_ptrs_.clear();
    n_leaves_ = 0;

    for (const auto& plant : plants) {
        plant->for_each_node_mut([&](Node& node) {
            auto* leaf = node.as_leaf();
            if (!leaf || leaf->leaf_size <= 0.0f || !node.parent) return;

            glm::vec3 pts[5];
            compute_leaf_corners(node, leaf, pts);

            int leaf_id = n_leaves_++;
            leaf_ptrs_.push_back(leaf);

            for (int s = 0; s < SAMPLES_PER_LEAF; s++) {
                query_verts_.push_back({
                    pts[s].x, pts[s].y, pts[s].z,
                    static_cast<float>(leaf_id * SAMPLES_PER_LEAF + s)
                });
            }
        });
    }
}

// ---------------------------------------------------------------------------
// compute_leaf_corners — matches Renderer::draw_leaf() geometry exactly
// pts[0]=p0 (base), pts[1]=p1 (right), pts[2]=p2 (tip), pts[3]=p3 (left), pts[4]=center
// ---------------------------------------------------------------------------

void LightSystem::compute_leaf_corners(const Node& node,
                                       const LeafNode* leaf,
                                       glm::vec3 out_pts[5]) {
    float size = leaf->leaf_size;
    float half = size * 0.70711f;  // size * sqrt(2)/2

    glm::vec3 outward = glm::normalize(node.position - node.parent->position);
    // Project outward onto the blade plane (perpendicular to facing)
    glm::vec3 proj = outward - leaf->facing * glm::dot(outward, leaf->facing);
    float proj_len = glm::length(proj);
    if (proj_len < 1e-4f) {
        proj = (std::abs(leaf->facing.y) < 0.9f)
            ? glm::normalize(glm::cross(leaf->facing, glm::vec3(0.0f, 1.0f, 0.0f)))
            : glm::normalize(glm::cross(leaf->facing, glm::vec3(1.0f, 0.0f, 0.0f)));
    } else {
        proj /= proj_len;
    }
    glm::vec3 width = glm::normalize(glm::cross(leaf->facing, proj));
    glm::vec3 down(0.0f, -1.0f, 0.0f);
    float droop = size * 0.15f;

    glm::vec3 base = node.position;
    out_pts[0] = base;                                              // p0: attachment
    out_pts[1] = base + proj * half + width * half + down * droop; // p1: right (drooped)
    out_pts[2] = base + proj * (half * 2.0f);                      // p2: tip
    out_pts[3] = base + proj * half - width * half + down * droop; // p3: left (drooped)
    out_pts[4] = base + proj * half;                               // center
}

// ---------------------------------------------------------------------------
// update_light_matrices — orthographic projection from directly above
// ---------------------------------------------------------------------------

void LightSystem::update_light_matrices() {
    // Eye just above the highest caster so the ortho near/far planes actually
    // encompass the scene.  With eye at y=100 and near=0/far=6 the objects end
    // up at view-z ≈ -95 to -100, entirely outside the frustum, so everything
    // was silently clipped and nothing was ever written to the slice textures.
    float eye_y      = max_y_ + 1.0f;          // 1 unit above tallest caster
    float near_plane = 1.0f;                    // eye to top of canopy
    float far_plane  = eye_y - min_y_ + 1.0f;  // eye to below deepest root + margin

    glm::vec3 sun_pos(0.0f, eye_y, 0.0f);
    light_view_ = glm::lookAt(sun_pos,
                              glm::vec3(0.0f, 0.0f, 0.0f),
                              glm::vec3(0.0f, 0.0f, -1.0f));
    light_proj_ = glm::ortho(-SCENE_HALF, SCENE_HALF, -SCENE_HALF, SCENE_HALF,
                             near_plane, far_plane);
}

// ---------------------------------------------------------------------------
// upload helpers
// ---------------------------------------------------------------------------

void LightSystem::upload_caster_vbo() {
    glBindBuffer(GL_ARRAY_BUFFER, caster_vbo_);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(caster_verts_.size() * sizeof(CasterVertex)),
                 caster_verts_.data(), GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void LightSystem::upload_query_vbo() {
    glBindBuffer(GL_ARRAY_BUFFER, query_vbo_);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(query_verts_.size() * sizeof(LeafSamplePoint)),
                 query_verts_.data(), GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void LightSystem::resize_output_texture(int n_samples) {
    // Output texture: 1 row × n_samples columns, RGBA32F (R = light_exposure, GBA unused).
    glBindTexture(GL_TEXTURE_2D, output_tex_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, n_samples, 1, 0, GL_RGBA, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, query_fbo_);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, output_tex_, 0);
    GLenum draw_buf = GL_COLOR_ATTACHMENT0;
    glDrawBuffers(1, &draw_buf);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

// ---------------------------------------------------------------------------
// shadow_collect_pass — 4 MRT draw calls, additive blending
// ---------------------------------------------------------------------------

void LightSystem::shadow_collect_pass() {
    glm::mat4 light_pv = light_proj_ * light_view_;

    shadow_shader_.use();
    shadow_shader_.set_mat4("u_light_pv", light_pv);
    shadow_shader_.set_float("u_min_y", min_y_);
    shadow_shader_.set_float("u_max_y", max_y_);

    glViewport(0, 0, SHADOW_RES, SHADOW_RES);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);   // additive accumulation
    glBlendEquation(GL_FUNC_ADD);

    glBindVertexArray(caster_vao_);
    int n_verts = static_cast<int>(caster_verts_.size());

    for (int batch = 0; batch < 4; batch++) {
        glBindFramebuffer(GL_FRAMEBUFFER, mrt_fbos_[batch]);
        // Clear all 8 attachments for this batch to 0 (no coverage yet).
        const float zero = 0.0f;
        for (int s = 0; s < 8; s++) {
            glClearBufferfv(GL_COLOR, s, &zero);
        }
        shadow_shader_.set_int("u_batch", batch);
        glDrawArrays(GL_TRIANGLES, 0, n_verts);
    }

    glBindVertexArray(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
}

// ---------------------------------------------------------------------------
// leaf_query_pass — 5 points per leaf, output to 1D texture
// ---------------------------------------------------------------------------

void LightSystem::leaf_query_pass() {
    int n_samples = n_leaves_ * SAMPLES_PER_LEAF;
    glm::mat4 light_pv = light_proj_ * light_view_;

    query_shader_.use();
    query_shader_.set_mat4("u_light_pv", light_pv);
    query_shader_.set_float("u_min_y", min_y_);
    query_shader_.set_float("u_max_y", max_y_);
    query_shader_.set_float("u_inv_output_width", 1.0f / static_cast<float>(n_samples));

    // Bind the slice array texture on unit 0.
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D_ARRAY, slice_array_tex_);
    query_shader_.set_int("u_slice_array", 0);

    glBindFramebuffer(GL_FRAMEBUFFER, query_fbo_);
    glViewport(0, 0, n_samples, 1);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);

    glEnable(GL_PROGRAM_POINT_SIZE);
    glBindVertexArray(query_vao_);
    glDrawArrays(GL_POINTS, 0, static_cast<int>(query_verts_.size()));
    glBindVertexArray(0);
    glDisable(GL_PROGRAM_POINT_SIZE);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glEnable(GL_DEPTH_TEST);
    glBindTexture(GL_TEXTURE_2D_ARRAY, 0);
}

// ---------------------------------------------------------------------------
// readback_results — read output texture, average 5 samples, write light_exposure
// ---------------------------------------------------------------------------

void LightSystem::readback_results() {
    int n_samples = n_leaves_ * SAMPLES_PER_LEAF;
    std::vector<float> pixels(n_samples * 4);  // RGBA32F

    glBindTexture(GL_TEXTURE_2D, output_tex_);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, pixels.data());
    glBindTexture(GL_TEXTURE_2D, 0);

    for (int i = 0; i < n_leaves_; i++) {
        float sum = 0.0f;
        for (int s = 0; s < SAMPLES_PER_LEAF; s++) {
            float val = pixels[(i * SAMPLES_PER_LEAF + s) * 4];  // R channel
            leaf_ptrs_[i]->sample_exposure[s] = val;
            sum += val;
        }
        leaf_ptrs_[i]->light_exposure = sum / static_cast<float>(SAMPLES_PER_LEAF);
    }
}

// ---------------------------------------------------------------------------
// draw_debug_slice — draw one depth-slice texture as a fullscreen overlay
// ---------------------------------------------------------------------------

void LightSystem::draw_debug_slice(int slice_index) {
    if (!debug_shader_.is_loaded()) return;
    slice_index = std::max(0, std::min(slice_index, NUM_SLICES - 1));

    // Query current viewport to pass correct screen size to shader.
    GLint vp[4];
    glGetIntegerv(GL_VIEWPORT, vp);

    debug_shader_.use();
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D_ARRAY, slice_array_tex_);
    debug_shader_.set_int("u_slice_array", 0);
    debug_shader_.set_int("u_slice_index", slice_index);
    glUniform2f(glGetUniformLocation(debug_shader_.id(), "u_screen_size"),
                static_cast<float>(vp[2]), static_cast<float>(vp[3]));

    glDisable(GL_DEPTH_TEST);
    glBindVertexArray(debug_vao_);
    glDrawArrays(GL_TRIANGLES, 0, 3);  // fullscreen triangle
    glBindVertexArray(0);
    glEnable(GL_DEPTH_TEST);
}

} // namespace botany
