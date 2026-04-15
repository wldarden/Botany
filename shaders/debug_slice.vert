#version 410 core

// Fullscreen triangle for debug slice overlay.
// No input attributes needed — positions are generated from vertex ID.

void main() {
    // Fullscreen triangle vertices in clip space.
    // gl_VertexID: 0 → (-1,-1), 1 → (3,-1), 2 → (-1,3)
    vec2 pos;
    if (gl_VertexID == 0) pos = vec2(-1.0, -1.0);
    else if (gl_VertexID == 1) pos = vec2( 3.0, -1.0);
    else                       pos = vec2(-1.0,  3.0);
    gl_Position = vec4(pos, 0.0, 1.0);
}
