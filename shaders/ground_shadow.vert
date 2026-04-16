#version 410 core

// Ground shadow vertex shader.
// Passes world position and light-space UV to the fragment shader so it can
// look up per-slice shadow coverage and apply the ambient + transmittance model.

layout(location = 0) in vec3 a_pos;
layout(location = 1) in vec3 a_normal;  // unused — present to match ground VAO layout
layout(location = 2) in vec3 a_color;

uniform mat4 uModel;
uniform mat4 uView;
uniform mat4 uProjection;
uniform mat4 u_light_view;  // light-space view matrix
uniform mat4 u_light_proj;  // light-space projection matrix (orthographic)

out vec3 v_world_pos;
out vec3 v_color;
out vec2 v_shadow_uv;  // [0,1]² in shadow map (slice array) texture space

void main() {
    vec4 world     = uModel * vec4(a_pos, 1.0);
    gl_Position    = uProjection * uView * world;
    v_world_pos    = world.xyz;
    v_color        = a_color;

    // Project world position into light space for slice-array UV lookup.
    // u_light_proj is orthographic so w == 1, but divide anyway for correctness.
    vec4 lc    = u_light_proj * u_light_view * world;
    v_shadow_uv = (lc.xy / lc.w) * 0.5 + 0.5;
}
