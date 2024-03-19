#version 460 core

layout (location = 0) in vec3 domainPosition;

layout (binding = 0, std430) buffer position
{
	float _pos[];
};

uniform bool draw_domain;

uniform mat4 view;
uniform mat4 proj;

out VS_OUT
{
	vec3 fPos;
	bool draw_domain;
} vs_out;

void main()
{
	float x = _pos[gl_InstanceID * 4];
	float y = _pos[gl_InstanceID * 4 + 1];
	float z = _pos[gl_InstanceID * 4 + 2];

	if(draw_domain)
	{
		vs_out.fPos = domainPosition;
	}
	else
	{
		vs_out.fPos = vec3(x, y, z);
	}

	gl_Position = proj * view * vec4(vs_out.fPos, 1.0f);
	vs_out.draw_domain = draw_domain;
}