#version 460 core

layout (location = 0) in vec3 position;

uniform mat4 view;
uniform mat4 proj;

out VS_OUT
{
    mat4 view;
    mat4 proj;
} vs_out;

void main()
{
	gl_Position = proj * view * vec4(position, 1.0f);
	vs_out.view = view;
	vs_out.proj = proj;
}