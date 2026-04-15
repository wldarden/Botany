#version 410 core

// Leaf query vertex shader.
// Each "vertex" is one of 5 sample points for a leaf (center + 4 corners).
// Placed as a GL_POINT at the correct texel in the 1D output texture.

layout(location = 0) in vec3  a_sample_pos;    // world-space sample position
layout(location = 1) in float a_output_index;  // flat texel index in output texture

uniform mat4  u_light_pv;
uniform float u_min_y;
uniform float u_max_y;
uniform float u_inv_output_width;  // 1.0 / (n_leaves * SAMPLES_PER_LEAF)

out vec2  v_light_uv;     // [0,1]² in shadow map space → which texel column/row to sample
out float v_leaf_depth;   // depth_01 for this sample point

void main() {
    vec4 ls  = u_light_pv * vec4(a_sample_pos, 1.0);
    vec3 ndc = ls.xyz / ls.w;

    v_light_uv = ndc.xy * 0.5 + 0.5;

    float range = max(u_max_y - u_min_y, 0.01);
    v_leaf_depth = clamp((u_max_y - a_sample_pos.y) / range, 0.0, 1.0);

    // Map output_index to a pixel position in the 1D output texture.
    // NDC x: center of texel a_output_index in [0, output_width).
    float x_ndc = (a_output_index + 0.5) * u_inv_output_width * 2.0 - 1.0;
    gl_Position  = vec4(x_ndc, 0.0, 0.0, 1.0);
    gl_PointSize = 1.0;
}
