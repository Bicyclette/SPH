#version 460 core

out vec4 color;
in vec3 normal;
in vec3 fPos;

uniform vec3 camPos;
uniform bool draw_domain;

struct DirectionalLight
{
    vec3 direction;
    vec3 color;
    float intensity;
};

DirectionalLight sun;

vec3 tone_mapping(vec3 c)
{
    return c / (vec3(1.0f) + c);
}

vec3 gamma_correction(vec3 c)
{
    return pow(c, vec3(1.0f/2.2f));
}

void main()
{
    sun.direction = normalize(vec3(-0.25f, -1.0f, 0.25f));
    sun.color = vec3(1.0f, 0.89f, 0.51f);
    sun.intensity = 1.0f;

    if(draw_domain)
    {
        color = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    }
    else
    {
        vec3 blue = vec3(0.0f, 0.2f, 0.75f);
        vec3 viewDir = normalize(camPos - fPos);
        vec3 n = normalize(normal);
        vec3 h = normalize(viewDir -sun.direction);
        float cosTheta = max(0.0f, dot(-sun.direction, n));
        float spec = pow(max(dot(n, h), 0.0f), 256.0f) * 0.5f;
        vec3 specular = sun.color * spec;
        vec3 colorResponse = cosTheta * sun.color * sun.intensity * blue + blue * 0.2f + spec;
        colorResponse = tone_mapping(colorResponse);
        colorResponse = gamma_correction(colorResponse);
        color = vec4(colorResponse, 1.0f);
    }
}