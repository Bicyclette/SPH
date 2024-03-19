#pragma once

#include "GLCommon.h"

class DirectionalLight
{
public:
	DirectionalLight(glm::vec3 const& iDirection, glm::vec3 const& iColor, float const& iIntensity)
	{
		m_direction = iDirection;
		m_color = iColor;
		m_intensity = iIntensity;
	}

	glm::vec3 m_direction;
	glm::vec3 m_color;
	float m_intensity;
};


class PointLight
{
public:
	static const glm::vec3 attenuation[12];

public:
	PointLight(glm::vec3 const& iPosition, glm::vec3 const& iColor, float iIntensity, glm::vec3 const& iAttenuation)
	{
		m_position = iPosition;
		m_color = iColor;
		m_intensity = iIntensity;
		m_ac = iAttenuation[0];
		m_al = iAttenuation[1];
		m_aq = iAttenuation[2];
	}

	glm::vec3 m_position;
	glm::vec3 m_color;
	float m_intensity;
	float m_ac;			// constant attenuation coefficient
	float m_al;			// linear attenuation coefficient
	float m_aq;			// quadratic attenuation coefficient
};