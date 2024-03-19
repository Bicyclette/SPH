#version 460 core

layout (location = 0) in vec3 pos;
layout (location = 1) in vec3 normal;

layout (binding = 0, std430) buffer position
{
	float _pos[];
};

out VS_OUT
{
	vec3 fPos;
	vec3 fNormal;
} vs_out;

mat4 model_domain = mat4(1.0f);
uniform bool draw_domain;

uniform mat4 view;
uniform mat4 proj;
uniform float particleScale;

void main()
{
	float x = _pos[gl_InstanceID * 4];
	float y = _pos[gl_InstanceID * 4 + 1];
	float z = _pos[gl_InstanceID * 4 + 2];

	vs_out.fPos = pos + vec3(x, y, z);
	if(draw_domain)
	{
		gl_Position = proj * view * model_domain * vec4(pos, 1.0f);
		vs_out.fNormal = normalize(vec3(transpose(inverse(model_domain)) * vec4(normal, 1.0f)));
	}
	else
	{
		vec4 c1 = vec4(particleScale, 0.0f, 0.0f, 0.0f);
		vec4 c2 = vec4(0.0f, particleScale, 0.0f, 0.0f);
		vec4 c3 = vec4(0.0f, 0.0f, particleScale, 0.0f);
		vec4 c4 = vec4(x, y, z, 1.0f);
		mat4 model = mat4(c1, c2, c3, c4);

		gl_Position = proj * view * model * vec4(pos, 1.0f);
		vs_out.fNormal = normalize(vec3(transpose(inverse(model)) * vec4(normal, 1.0f)));
	}
}