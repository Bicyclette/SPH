#version 460 core

layout (location = 0) in vec3 domainPosition;

layout (binding = 0, std430) buffer position
{
	float _pos[];
};

uniform mat4 view;
uniform mat4 proj;

out VS_OUT
{
	vec3 pos;
	mat4 view;
	mat4 proj;
} vs_out;

void main()
{
	float x = _pos[gl_InstanceID * 4];
	float y = _pos[gl_InstanceID * 4 + 1];
	float z = _pos[gl_InstanceID * 4 + 2];
	vs_out.pos = vec3(x, y, z);

	gl_Position = proj * view * vec4(vs_out.pos, 1.0f);
	vs_out.view = view;
	vs_out.proj = proj;
}