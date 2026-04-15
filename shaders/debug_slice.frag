#version 410 core

// Debug overlay: draws one layer of the slice array texture fullscreen.
// Low coverage (no shadow) = yellow, high coverage (deep shadow) = purple.

out vec4 frag_color;

uniform sampler2DArray u_slice_array;
uniform int            u_slice_index;

void main() {
    vec2 uv = gl_FragCoord.xy / vec2(512.0);
    float coverage = texture(u_slice_array, vec3(uv, float(u_slice_index))).r;
    coverage = clamp(coverage, 0.0, 1.0);
    // Yellow (no shadow) → Purple (full shadow)
    vec3 col = mix(vec3(1.0, 1.0, 0.0), vec3(0.5, 0.0, 0.8), coverage);
    frag_color = vec4(col, 0.7);
}
