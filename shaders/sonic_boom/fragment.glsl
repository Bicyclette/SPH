#version 460 core

out vec4 color;

in VS_OUT
{
	vec3 fPos;
	flat bool draw_domain;
} fs_in;

void main()
{
    if(fs_in.draw_domain)
    {
        color = vec4(vec3(0.0f), 1.0f);
    }
    else
    {
        color = vec4(0.05f, 0.05f, 0.85f, 1.0f);
    }
}