#version 410 core

// Debug overlay: draws one layer of the slice array texture fullscreen.
// Low coverage (no shadow) = yellow, high coverage (deep shadow) = purple.

out vec4 frag_color;

uniform sampler2DArray u_slice_array;
uniform int            u_slice_index;
uniform vec2           u_screen_size;  // actual window width/height in pixels

void main() {
    // Map fragment to [0,1] UV across the window so the shadow map
    // fills the full screen regardless of window dimensions.
    vec2 uv = gl_FragCoord.xy / u_screen_size;
    float coverage = texture(u_slice_array, vec3(uv, float(u_slice_index))).r;
    coverage = clamp(coverage, 0.0, 1.0);
    // Yellow (no shadow) → Purple (full shadow)
    vec3 col = mix(vec3(1.0, 1.0, 0.0), vec3(0.5, 0.0, 0.8), coverage);
    frag_color = vec4(col, 0.7);
}
