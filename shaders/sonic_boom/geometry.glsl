#version 460 core

layout (points) in;
layout (triangle_strip, max_vertices=60) out;

uniform float particleRadius;

in VS_OUT
{
	vec3 pos;
    mat4 view;
    mat4 proj;
} gs_in[];

out vec3 normal;
out vec3 fPos;

void main()
{
	float scale = particleRadius * 2.0f;
	mat4 viewProj = gs_in[0].proj * gs_in[0].view;
	vec4 pos = vec4(gs_in[0].pos, 1.0f);
	fPos = gs_in[0].pos;
	
	gl_Position = viewProj * (pos + vec4(0.0f, -1.0f, 0.0f, 0.0f) * scale);
	normal = vec3(0.0f, -1.0f, 0.0f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.723600f, -0.447215f, 0.525720f, 0.0f) * scale);
	normal = vec3(0.723600f, -0.447215f, 0.525720f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.276385f, -0.447215f, 0.850640f, 0.0f) * scale);
	normal = vec3(-0.276385f, -0.447215f, 0.850640f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(0.723600f, -0.447215f, 0.525720f, 0.0f) * scale);
	normal = vec3(0.723600f, -0.447215f, 0.525720f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.0f, -1.0f, 0.0f, 0.0f) * scale);
	normal = vec3(0.0f, -1.0f, 0.0f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.723600f, -0.447215f, -0.525720f, 0.0f) * scale);
	normal = vec3(0.723600f, -0.447215f, -0.525720f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(0.0f, -1.0f, 0.0f, 0.0f) * scale);
	normal = vec3(0.0f, -1.0f, 0.0f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.276385f, -0.447215f, 0.850640f, 0.0f) * scale);
	normal = vec3(-0.276385f, -0.447215f, 0.850640f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.894425f, -0.447215f, 0.0f, 0.0f) * scale);
	normal = vec3(-0.894425f, -0.447215f, 0.0f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(0.0f, -1.0f, 0.0f, 0.0f) * scale);
	normal = vec3(0.0f, -1.0f, 0.0f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.894425f, -0.447215f, 0.0f, 0.0f) * scale);
	normal = vec3(-0.894425f, -0.447215f, 0.0f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.276385f, -0.447215f, -0.850640f, 0.0f) * scale);
	normal = vec3(-0.276385f, -0.447215f, -0.850640f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(0.0f, -1.0f, 0.0f, 0.0f) * scale);
	normal = vec3(0.0f, -1.0f, 0.0f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.276385f, -0.447215f, -0.850640f, 0.0f) * scale);
	normal = vec3(-0.276385f, -0.447215f, -0.850640f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.723600f, -0.447215f, -0.525720f, 0.0f) * scale);
	normal = vec3(0.723600f, -0.447215f, -0.525720f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(0.723600f, -0.447215f, 0.525720f, 0.0f) * scale);
	normal = vec3(0.723600f, -0.447215f, 0.525720f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.723600f, -0.447215f, -0.525720f, 0.0f) * scale);
	normal = vec3(0.723600f, -0.447215f, -0.525720f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.894425f, 0.447215f, 0.0f, 0.0f) * scale);
	normal = vec3(0.894425f, 0.447215f, 0.0f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(-0.276385f, -0.447215f, 0.850640f, 0.0f) * scale);
	normal = vec3(-0.276385f, -0.447215f, 0.850640f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.723600f, -0.447215f, 0.525720f, 0.0f) * scale);
	normal = vec3(0.723600f, -0.447215f, 0.525720f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.276385f, 0.447215f, 0.850640f, 0.0f) * scale);
	normal = vec3(0.276385f, 0.447215f, 0.850640f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(-0.894425f, -0.447215f, 0.0f, 0.0f) * scale);
	normal = vec3(-0.894425f, -0.447215f, 0.0f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.276385f, -0.447215f, 0.850640f, 0.0f) * scale);
	normal = vec3(-0.276385f, -0.447215f, 0.850640f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.723600f, 0.447215f, 0.525720f, 0.0f) * scale);
	normal = vec3(-0.723600f, 0.447215f, 0.525720f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(-0.276385f, -0.447215f, -0.850640f, 0.0f) * scale);
	normal = vec3(-0.276385f, -0.447215f, -0.850640f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.894425f, -0.447215f, 0.0f, 0.0f) * scale);
	normal = vec3(-0.894425f, -0.447215f, 0.0f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.723600f, 0.447215f, -0.525720f, 0.0f) * scale);
	normal = vec3(-0.723600f, 0.447215f, -0.525720f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(0.723600f, -0.447215f, -0.525720f, 0.0f) * scale);
	normal = vec3(0.723600f, -0.447215f, -0.525720f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.276385f, -0.447215f, -0.850640f, 0.0f) * scale);
	normal = vec3(-0.276385f, -0.447215f, -0.850640f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.276385f, 0.447215f, -0.850640f, 0.0f) * scale);
	normal = vec3(0.276385f, 0.447215f, -0.850640f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(0.723600f, -0.447215f, 0.525720f, 0.0f) * scale);
	normal = vec3(0.723600f, -0.447215f, 0.525720f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.894425f, 0.447215f, 0.0f, 0.0f) * scale);
	normal = vec3(0.894425f, 0.447215f, 0.0f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.276385f, 0.447215f, 0.850640f, 0.0f) * scale);
	normal = vec3(0.276385f, 0.447215f, 0.850640f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(-0.276385f, -0.447215f, 0.850640f, 0.0f) * scale);
	normal = vec3(-0.276385f, -0.447215f, 0.850640f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.276385f, 0.447215f, 0.850640f, 0.0f) * scale);
	normal = vec3(0.276385f, 0.447215f, 0.850640f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.723600f, 0.447215f, 0.525720f, 0.0f) * scale);
	normal = vec3(-0.723600f, 0.447215f, 0.525720f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(-0.894425f, -0.447215f, 0.0f, 0.0f) * scale);
	normal = vec3(-0.894425f, -0.447215f, 0.0f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.723600f, 0.447215f, 0.525720f, 0.0f) * scale);
	normal = vec3(-0.723600f, 0.447215f, 0.525720f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.723600f, 0.447215f, -0.525720f, 0.0f) * scale);
	normal = vec3(-0.723600f, 0.447215f, -0.525720f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(-0.276385f, -0.447215f, -0.850640f, 0.0f) * scale);
	normal = vec3(-0.276385f, -0.447215f, -0.850640f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.723600f, 0.447215f, -0.525720f, 0.0f) * scale);
	normal = vec3(-0.723600f, 0.447215f, -0.525720f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.276385f, 0.447215f, -0.850640f, 0.0f) * scale);
	normal = vec3(0.276385f, 0.447215f, -0.850640f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(0.723600f, -0.447215f, -0.525720f, 0.0f) * scale);
	normal = vec3(0.723600f, -0.447215f, -0.525720f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.276385f, 0.447215f, -0.850640f, 0.0f) * scale);
	normal = vec3(0.276385f, 0.447215f, -0.850640f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.894425f, 0.447215f, 0.0f, 0.0f) * scale);
	normal = vec3(0.894425f, 0.447215f, 0.0f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(0.276385f, 0.447215f, 0.850640f, 0.0f) * scale);
	normal = vec3(0.276385f, 0.447215f, 0.850640f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.894425f, 0.447215f, 0.0f, 0.0f) * scale);
	normal = vec3(0.894425f, 0.447215f, 0.0f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.0f, 1.0f, 0.0f, 0.0f) * scale);
	normal = vec3(0.0f, 1.0f, 0.0f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(-0.723600f, 0.447215f, 0.525720f, 0.0f) * scale);
	normal = vec3(-0.723600f, 0.447215f, 0.525720f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.276385f, 0.447215f, 0.850640f, 0.0f) * scale);
	normal = vec3(0.276385f, 0.447215f, 0.850640f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.0f, 1.0f, 0.0f, 0.0f) * scale);
	normal = vec3(0.0f, 1.0f, 0.0f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(-0.723600f, 0.447215f, -0.525720f, 0.0f) * scale);
	normal = vec3(-0.723600f, 0.447215f, -0.525720f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.723600f, 0.447215f, 0.525720f, 0.0f) * scale);
	normal = vec3(-0.723600f, 0.447215f, 0.525720f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.0f, 1.0f, 0.0f, 0.0f) * scale);
	normal = vec3(0.0f, 1.0f, 0.0f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(0.276385f, 0.447215f, -0.850640f, 0.0f) * scale);
	normal = vec3(0.276385f, 0.447215f, -0.850640f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(-0.723600f, 0.447215f, -0.525720f, 0.0f) * scale);
	normal = vec3(-0.723600f, 0.447215f, -0.525720f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.0f, 1.0f, 0.0f, 0.0f) * scale);
	normal = vec3(0.0f, 1.0f, 0.0f);
	EmitVertex();
	EndPrimitive();

	gl_Position = viewProj * (pos + vec4(0.894425f, 0.447215f, 0.0f, 0.0f) * scale);
	normal = vec3(0.894425f, 0.447215f, 0.0f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.276385f, 0.447215f, -0.850640f, 0.0f) * scale);
	normal = vec3(0.276385f, 0.447215f, -0.850640f);
	EmitVertex();
	gl_Position = viewProj * (pos + vec4(0.0f, 1.0f, 0.0f, 0.0f) * scale);
	normal = vec3(0.0f, 1.0f, 0.0f);
	EmitVertex();
	EndPrimitive();
}