#version 460 core

out vec4 color;

uniform bool draw_domain;

void main()
{
    if(draw_domain)
    {
        color = vec4(vec3(0.0f), 1.0f);
    }
    else
    {
        color = vec4(0.05f, 0.05f, 0.85f, 1.0f);
    }
}