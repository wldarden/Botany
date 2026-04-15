#pragma once

#include <memory>
#include <vector>
#include <glm/vec3.hpp>
#include <glm/mat4x4.hpp>
#include "renderer/shader.h"

namespace botany {

class Plant;
struct LeafNode;

// GPU-accelerated deep shadow map for per-leaf light occlusion.
//
// Algorithm: Opacity Shadow Maps (Kim & Neumann 2001) adapted for foliage.
// - 32 depth slices with power-curve distribution (top-heavy: more slices near canopy top).
// - Slice textures are 512×512 GL_R32F, one per depth slice, stored as a GL_TEXTURE_2D_ARRAY.
// - Shadow collect pass: 4 MRT draw calls (8 slices each) with additive blending.
// - Leaf query pass: 5 sample points per leaf (center + 4 corners), renders to 1D output texture.
// - Transmittance model: T *= max(0, 1 - coverage) per slice, multiplicative (Beer-Lambert).
//
// Coordinate convention:
// - depth_01 = 0: top of canopy (closest to sun, most important for accurate occlusion)
// - depth_01 = 1: bottom of canopy (farthest from sun)
// - Power-curve mapping: slice_index = int(NUM_SLICES * sqrt(depth_01))
//   Concentrates slices near depth_01=0 where each new occluder has the most impact.
//
// Sun direction is currently vertical (0, -1, 0). The sun_direction field is exposed
// for future day/night cycle support — change it and call update(); no shader rewrites needed.
// (At non-vertical angles, depth_01 computation should use dot(pos, -sun_direction) instead
// of world Y. Currently hardcoded to world Y for clarity.)

class LightSystem {
public:
    static constexpr int SHADOW_RES      = 512;   // shadow map width and height in pixels
    static constexpr int NUM_SLICES      = 32;    // depth slices (power-curve distributed)
    static constexpr int SAMPLES_PER_LEAF = 5;    // sample points per leaf (center + 4 corners)
    static constexpr float SCENE_HALF   = 12.0f;  // scene footprint: ±12 dm = 24 dm × 24 dm

    // Sun direction: unit vector pointing TOWARD the scene (downward = overhead sun).
    // Change to any direction for oblique angles / day-night cycle.
    glm::vec3 sun_direction = glm::normalize(glm::vec3(-1.0f, -2.0f, 0.5f));  // angled test: ~60° from horizontal

    bool init(const std::string& shader_dir);
    void shutdown();

    // Main entry point. Collects casters + query points from plants, runs shadow and query
    // passes, writes light_exposure to every LeafNode. Call once per light update interval.
    void update(const std::vector<std::unique_ptr<Plant>>& plants);

    // Debug: draw one depth-slice texture as a full-screen overlay.
    // slice_index in [0, NUM_SLICES). 0 = topmost (closest to sun).
    void draw_debug_slice(int slice_index = 0);

    bool is_initialized() const { return initialized_; }
    uint32_t slice_tex() const { return slice_array_tex_; }
    glm::mat4 light_pv() const { return light_proj_ * light_view_; }

private:
    bool initialized_ = false;

    // --- GL resources ---
    uint32_t slice_array_tex_ = 0;   // GL_TEXTURE_2D_ARRAY, 512×512×32, GL_R32F
    uint32_t mrt_fbos_[4]     = {};  // 4 FBOs, each binding 8 consecutive layers
    uint32_t caster_vao_      = 0;
    uint32_t caster_vbo_      = 0;
    uint32_t query_vao_       = 0;
    uint32_t query_vbo_       = 0;
    uint32_t output_tex_      = 0;   // GL_TEXTURE_2D, 1×(SAMPLES_PER_LEAF*N_leaves), GL_RGBA32F
    uint32_t query_fbo_       = 0;
    uint32_t debug_vao_       = 0;   // fullscreen triangle for debug overlay
    uint32_t debug_vbo_       = 0;

    Shader shadow_shader_;   // shadow_collect.vert + shadow_collect.frag
    Shader query_shader_;    // leaf_query.vert + leaf_query.frag
    Shader debug_shader_;    // debug_slice.vert + debug_slice.frag

    // --- Per-frame state ---
    struct LeafSamplePoint {
        float x, y, z;      // world-space sample position
        float output_index; // flat index into output_tex_
    };

    struct CasterVertex {
        float x, y, z;   // world-space position
        float opacity;    // leaf_opacity from genome (stems use 1.0)
    };

    std::vector<CasterVertex>   caster_verts_;
    std::vector<LeafSamplePoint> query_verts_;
    std::vector<LeafNode*>      leaf_ptrs_;   // parallel to leaf index (leaf_id)
    int n_leaves_ = 0;

    float min_y_ = 0.0f;  // bottom of active Y range (adaptive)
    float max_y_ = 20.0f; // top of active Y range (adaptive)

    glm::mat4 light_view_;
    glm::mat4 light_proj_;

    // --- Internal passes ---
    void collect_casters(const std::vector<std::unique_ptr<Plant>>& plants);
    void build_query_vbo(const std::vector<std::unique_ptr<Plant>>& plants);
    void update_light_matrices();
    void shadow_collect_pass();
    void leaf_query_pass();
    void readback_results();

    // Geometry helpers (same math as Renderer::draw_leaf)
    static void compute_leaf_corners(const struct Node& node,
                                     const LeafNode* leaf,
                                     glm::vec3 out_pts[5]);

    bool check_texture_units() const;
    void upload_caster_vbo();
    void upload_query_vbo();
    void resize_output_texture(int n_samples);
};

} // namespace botany
