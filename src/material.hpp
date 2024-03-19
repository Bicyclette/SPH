#pragma once

#include "utils.h"
#include "GLCommon.h"

enum class TextureType
{
	ALBEDO,
	METALLIC,
	ROUGHNESS,
	NORMAL
};

struct Texture2D
{
	TextureType m_type;
	GLuint m_tex;
};

class Material
{
public:
	Material()
	{
		m_albedo = glm::vec3(1.0f);
		m_metallic = 0.0f;
		m_roughness = 0.05f;
	}

	Material(glm::vec3 const& iAlbedo, float const & iMetallic, float const & iRoughness)
	{
		m_albedo = iAlbedo;
		m_metallic = iMetallic;
		m_roughness = iRoughness;
	}

	glm::vec3 m_albedo;
	float m_metallic;
	float m_roughness;
};