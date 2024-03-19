#version 460 core

#define M_PI 3.1415926535897932384626433832795

struct DirectionalLight
{
	vec3 direction;
	vec3 color;
	float intensity;
};

struct PointLight
{
	vec3 position;
	vec3 color;
	float intensity;
	float ac;
	float al;
	float aq;
};

struct Material
{
	vec3 albedo;
	float roughness;
	float metallic;
};

in VS_OUT
{
	vec3 fPos;
	vec3 fNormal;
} fs_in;

out vec4 color_response;

uniform DirectionalLight directionalLight[10];
uniform int numDirectionalLight;
uniform PointLight pointLight[10];
uniform int numPointLight;

uniform vec3 camPos;
uniform Material material;

float normal_distribution(vec3 wi, vec3 wh, vec3 n, float roughness) // GGX distribution
{
	float r2 = roughness * roughness;
	float n_dot_wh = max(dot(n, wh), 0.0f);
	float n_dot_wh2 = n_dot_wh * n_dot_wh;
	float denom = 1.0f + (r2 - 1.0f) * n_dot_wh2;
	return r2 / (M_PI * denom * denom);
}

vec3 fresnel_term(vec3 wi, vec3 wh, vec3 F0)
{
	float wi_dot_wh = max(dot(wi, wh), 0.0f);
	return F0 + (1.0f - F0) * pow((1.0f - max(0.0f, wi_dot_wh)), 5);
}

float GGX_schlick(vec3 w, vec3 n, float roughness)
{
	float k = roughness * sqrt(2.0f / M_PI);
	float n_dot_w = max(dot(n, w), 0.0f);
	return n_dot_w / (n_dot_w * (1.0f - k) + k);
}

float geometric_term(vec3 wi, vec3 wo, vec3 n, float roughness)
{
	return GGX_schlick(wi, n, roughness) * GGX_schlick(wo, n, roughness);
}

vec3 BRDF(vec3 wi, vec3 wo, vec3 wh, vec3 n, vec3 F0, vec3 albedo, float metallic, float roughness)
{
	float n_dot_wi = max(dot(n, wi), 0.0f);
	float n_dot_wo = max(dot(n, wo), 0.0f);

	// specular part
	float D = normal_distribution(wi, wh, n, roughness);
	vec3 F = fresnel_term(wi, wh, F0);
	float G = geometric_term(wi, wo, n, roughness);
	vec3 specular = (D * F * G) / (4 * n_dot_wi * n_dot_wo + 0.001);

	// diffuse part
	vec3 Kd = (vec3(1.0f) - F);
	Kd *= (1.0f - metallic);
	vec3 diffuse = Kd / M_PI;

	// result
	return diffuse * albedo + specular;
}

void main()
{
	vec3 albedo = material.albedo;
	float metallic = material.metallic;
	float roughness = material.roughness;

	vec3 Lo = vec3(0.0f);
	vec3 n = fs_in.fNormal;
	vec3 wo = normalize(camPos - fs_in.fPos);
	vec3 F0 = mix(vec3(0.04f), albedo, metallic);

	for(int i = 0; i < numDirectionalLight; ++i)
	{
		DirectionalLight light = directionalLight[i];

		// radiance
		vec3 wi = normalize(-light.direction);
		vec3 wh = normalize(wi + wo);
		float attenuation = 1.0f;
		vec3 radiance = light.color * light.intensity * attenuation;
		float n_dot_wi = max(dot(n, wi), 0.0f);

		// BRDF
		vec3 brdf = BRDF(wi, wo, wh, n, F0, albedo, metallic, roughness);

		// add diffuse and specular parts
		Lo += radiance * brdf * n_dot_wi;
	}

	for(int i = 0; i < numPointLight; ++i)
	{
		PointLight light = pointLight[i];

		// radiance
		vec3 wi = normalize(light.position - fs_in.fPos);
		vec3 wh = normalize(wi + wo);
		float d = length(light.position - fs_in.fPos); // distance of fragment to light
		float attenuation = light.intensity / ( light.ac + (light.al * d) + (light.aq * d * d) );		
		vec3 radiance = light.color * attenuation;
		float n_dot_wi = max(dot(n, wi), 0.0f);

		// BRDF
		vec3 brdf = BRDF(wi, wo, wh, n, F0, albedo, metallic, roughness);

		// add diffuse and specular parts
		Lo += radiance * brdf * n_dot_wi;
	}

	// ambient color
	vec3 ambient = albedo * vec3(0.03f);
	vec3 color = ambient + Lo;

	// to LDR
	color = color / (color + vec3(1.0f));

	// gamma correction
	color = pow(color, vec3(1.0f/2.2f));

	color_response = vec4(color, 1.0f);
}