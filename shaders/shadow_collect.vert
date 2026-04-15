#version 410 core

// Shadow collect vertex shader.
// Transforms leaf/stem geometry into light (sun) space for depth-slice accumulation.
// depth_01 = 0 means top of canopy (closest to sun), 1 means bottom.

layout(location = 0) in vec3  a_world_pos;
layout(location = 1) in float a_opacity;

uniform mat4  u_light_pv;    // light projection × view (orthographic from sun)
uniform float u_min_y;       // min depth along sun axis (dot(pos, -sun_dir))
uniform float u_max_y;       // max depth along sun axis
uniform vec3  u_sun_dir;     // normalized direction light travels (toward scene)

out float v_depth_01;
out float v_opacity;

void main() {
    gl_Position = u_light_pv * vec4(a_world_pos, 1.0);

    // depth_01: 0 = closest to sun, 1 = farthest. Works for any sun angle.
    float depth = dot(a_world_pos, -u_sun_dir);
    float range = max(u_max_y - u_min_y, 0.01);
    v_depth_01 = clamp((u_max_y - depth) / range, 0.0, 1.0);
    v_opacity  = a_opacity;
}
