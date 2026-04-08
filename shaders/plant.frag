#version 410 core

in vec3 vNormal;
in vec3 vFragPos;
in vec3 vColor;

out vec4 FragColor;

uniform vec3 uLightDir;
uniform vec3 uAmbient;

void main() {
    vec3 norm = normalize(vNormal);
    vec3 lightDir = normalize(-uLightDir);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 result = (uAmbient + diff * vec3(1.0)) * vColor;
    FragColor = vec4(result, 1.0);
}
