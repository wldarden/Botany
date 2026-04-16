#version 410 core

// Ground shadow fragment shader.
// Computes transmittance T by accumulating (1 - coverage) over all depth slices
// above this fragment (same algorithm as leaf_query.frag), then applies:
//   FragColor = ground_color * (0.2 + T * 0.8)
// giving an ambient floor of 0.2 so fully shadowed ground is never black.

in vec3 v_world_pos;
in vec3 v_color;
in vec2 v_shadow_uv;

uniform vec3  u_sun_dir;     // normalized direction light travels (toward scene)
uniform float u_min_y;       // min depth along sun axis over caster set
uniform float u_max_y;       // max depth along sun axis over caster set
uniform sampler2DArray u_slice_array;  // 32-slice shadow array (GL_TEXTURE_2D_ARRAY)
uniform int   u_num_slices;  // NUM_SLICES (32)

out vec4 FragColor;

void main() {
    // Depth of this ground fragment along the sun axis.
    // Matches the convention in shadow_collect.vert / leaf_query.vert.
    float depth    = dot(v_world_pos, -u_sun_dir);
    float range    = max(u_max_y - u_min_y, 0.01);
    float depth_01 = clamp((u_max_y - depth) / range, 0.0, 1.0);

    // Which slice this ground point falls into (power-curve distribution).
    int ground_slice = int(float(u_num_slices) * sqrt(depth_01));
    ground_slice = clamp(ground_slice, 0, u_num_slices - 1);

    // Accumulate transmittance from all slices above this fragment (0 .. ground_slice-1).
    // Slice 0 = top of canopy (closest to sun).  T=1 means no occlusion.
    float T = 1.0;
    for (int s = 0; s < ground_slice; s++) {
        float coverage = texture(u_slice_array, vec3(v_shadow_uv, float(s))).r;
        T *= max(0.0, 1.0 - clamp(coverage, 0.0, 1.0));
    }

    FragColor = vec4(v_color * (0.2 + T * 0.8), 1.0);
}
