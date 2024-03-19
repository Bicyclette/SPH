#version 460 core

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;
layout (location = 2) in mat4 instanceModel;

uniform bool instanced_rendering;
uniform mat4 model;
uniform mat4 view;
uniform mat4 proj;

out VS_OUT
{
	vec3 fPos;
	vec3 fNormal;
} vs_out;

void main()
{
	mat4 modelMat;
	if(instanced_rendering)
	{
		modelMat = instanceModel;
	}
	else
	{
		modelMat = model;
	}

	gl_Position = proj * view * modelMat * vec4(position, 1.0f);
	vs_out.fPos = (modelMat * vec4(position, 1.0f)).xyz;
	vs_out.fNormal = normalize(vec3(transpose(inverse(modelMat)) * vec4(normal, 1.0f)));
}