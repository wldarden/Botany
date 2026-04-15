#version 410 core

layout(location = 0) in vec3 a_pos;
layout(location = 1) in vec3 a_normal;
layout(location = 2) in vec3 a_color;

uniform mat4 uModel;
uniform mat4 uView;
uniform mat4 uProjection;
uniform mat4 u_light_pv;  // light-space projection × view

out vec3 v_color;
out vec3 v_normal;
out vec2 v_shadow_uv;  // [0,1]² in shadow map space

void main() {
    gl_Position = uProjection * uView * uModel * vec4(a_pos, 1.0);
    v_color  = a_color;
    v_normal = a_normal;

    // Project world position into light space to get shadow map UV.
    vec4 lc = u_light_pv * uModel * vec4(a_pos, 1.0);
    v_shadow_uv = lc.xy * 0.5 + 0.5;
}
