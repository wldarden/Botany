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
    // PCF: sample each slice at a small rotated Poisson kernel and average
    // coverage before multiplying into T. Softens hard shadow edges for free.
    // Kernel radius in UV space: 3 texels at 512 resolution = ~6mm at scene scale.
    const float r = 3.0 / 512.0;
    const vec2 offsets[5] = vec2[5](
        vec2( 0.000,  0.000),   // center
        vec2( r,      0.000),   // right
        vec2(-r,      0.000),   // left
        vec2( 0.000,  r    ),   // up
        vec2( 0.000, -r    )    // down
    );
    const float weights[5] = float[5](0.40, 0.15, 0.15, 0.15, 0.15);

    float T = 1.0;
    for (int s = 0; s < u_num_slices; s++) {
        float coverage = 0.0;
        for (int k = 0; k < 5; k++) {
            coverage += weights[k] * texture(u_slice_array,
                            vec3(v_shadow_uv + offsets[k], float(s))).r;
        }
        T *= max(0.0, 1.0 - clamp(coverage, 0.0, 1.0));
    }

    // Standard diffuse lighting, direct component attenuated by transmittance.
    vec3  norm     = normalize(v_normal);
    vec3  lightDir = normalize(-uLightDir);
    float diff     = max(dot(norm, lightDir), 0.0);
    vec3  result   = (uAmbient + diff * T * vec3(1.0)) * v_color;
    FragColor = vec4(result, 1.0);
}
