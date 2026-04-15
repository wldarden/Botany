#version 410 core

// Shadow collect vertex shader.
// Transforms leaf/stem geometry into light (sun) space for depth-slice accumulation.
// depth_01 = 0 means top of canopy (closest to sun), 1 means bottom.

layout(location = 0) in vec3  a_world_pos;
layout(location = 1) in float a_opacity;

uniform mat4  u_light_pv;    // light projection × view (orthographic from above)
uniform float u_min_y;       // bottom of adaptive Y range
uniform float u_max_y;       // top of adaptive Y range

out float v_depth_01;
out float v_opacity;

void main() {
    gl_Position = u_light_pv * vec4(a_world_pos, 1.0);

    // depth_01: 0 = top of canopy (closest to sun), 1 = bottom.
    // Uses world Y directly (vertical sun). For oblique sun: dot(a_world_pos, -sun_dir).
    float range = max(u_max_y - u_min_y, 0.01);
    v_depth_01 = clamp((u_max_y - a_world_pos.y) / range, 0.0, 1.0);
    v_opacity  = a_opacity;
}
