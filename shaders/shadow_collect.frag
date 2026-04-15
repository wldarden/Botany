#version 410 core

// Shadow collect fragment shader.
// Routes each caster fragment into the correct depth slice via MRT.
// 4 batches of 8 slices each. Additive blending accumulates coverage.
//
// Slice index uses power-curve distribution (p=2, via sqrt):
//   slice = int(NUM_SLICES * sqrt(depth_01))
// This concentrates slices near depth_01=0 (top of canopy) where each
// new occluder has the largest impact on transmittance.

in float v_depth_01;
in float v_opacity;

// 8 outputs — one per slice in the current batch.
// GL_R32F per slice; accumulated with additive blending.
layout(location = 0) out float frag0;
layout(location = 1) out float frag1;
layout(location = 2) out float frag2;
layout(location = 3) out float frag3;
layout(location = 4) out float frag4;
layout(location = 5) out float frag5;
layout(location = 6) out float frag6;
layout(location = 7) out float frag7;

uniform int u_batch;  // 0, 1, 2, or 3 — which group of 8 slices this draw covers

void main() {
    // Power-curve slice assignment: more slices near top (small depth_01).
    int slice = int(32.0 * sqrt(v_depth_01));
    slice = clamp(slice, 0, 31);

    int batch = slice / 8;
    if (batch != u_batch) discard;  // not our batch

    int slot = slice % 8;  // 0..7, selects which MRT output

    // Write opacity to the correct slot; 0.0 elsewhere.
    // Additive blending (GL_ONE, GL_ONE) accumulates per pixel.
    float v = v_opacity;
    frag0 = (slot == 0) ? v : 0.0;
    frag1 = (slot == 1) ? v : 0.0;
    frag2 = (slot == 2) ? v : 0.0;
    frag3 = (slot == 3) ? v : 0.0;
    frag4 = (slot == 4) ? v : 0.0;
    frag5 = (slot == 5) ? v : 0.0;
    frag6 = (slot == 6) ? v : 0.0;
    frag7 = (slot == 7) ? v : 0.0;
}
