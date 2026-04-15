#version 410 core

in vec3 v_color;
in vec3 v_normal;
in vec2 v_shadow_uv;

uniform vec3 uLightDir;
uniform vec3 uAmbient;
uniform sampler2DArray u_slice_array;  // 32-slice shadow array
uniform int  u_num_slices;             // NUM_SLICES (32)

out vec4 FragColor;

void main() {
    // Accumulate transmittance through all slices above ground (all of them).
    // T = product of max(0, 1 - coverage_s) — Beer-Lambert, same as leaf query.
    float T = 1.0;
    for (int s = 0; s < u_num_slices; s++) {
        float coverage = texture(u_slice_array, vec3(v_shadow_uv, float(s))).r;
        T *= max(0.0, 1.0 - clamp(coverage, 0.0, 1.0));
    }

    // Standard diffuse lighting, with direct component attenuated by transmittance.
    vec3  norm     = normalize(v_normal);
    vec3  lightDir = normalize(-uLightDir);
    float diff     = max(dot(norm, lightDir), 0.0);
    vec3  result   = (uAmbient + diff * T * vec3(1.0)) * v_color;
    FragColor = vec4(result, 1.0);
}
