#version 410 core

// Leaf query fragment shader.
// For each sample point, computes transmittance T by multiplying (1 - coverage) across
// all depth slices above this point's depth. T is the fraction of sunlight that reaches
// this leaf sample. Writes T to the R channel of the 1D output texture.

in vec2  v_light_uv;
in float v_leaf_depth;

// All 32 depth slices as a 2D array texture (single sampler, no texture-unit limit issue).
// Texture coordinates: (u, v, layer) where layer = slice index.
uniform sampler2DArray u_slice_array;

out vec4 frag_output;  // R = light_exposure (T), GBA unused

void main() {
    // Power-curve slice for this sample: same formula as shadow_collect.frag.
    int leaf_slice = int(32.0 * sqrt(v_leaf_depth));
    leaf_slice = clamp(leaf_slice, 0, 31);

    // Accumulate transmittance from all slices above this sample (slices 0..leaf_slice-1).
    // Slice 0 = top of canopy (closest to sun). T starts at 1.0 (full sun).
    float T = 1.0;
    for (int s = 0; s < leaf_slice; s++) {
        float coverage = texture(u_slice_array, vec3(v_light_uv, float(s))).r;
        T *= max(0.0, 1.0 - clamp(coverage, 0.0, 1.0));
    }

    frag_output = vec4(T, 0.0, 0.0, 1.0);
}
