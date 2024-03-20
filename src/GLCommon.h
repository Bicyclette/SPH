#pragma once

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/hash.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/epsilon.hpp>

#define VIEWPORT_WIDTH 1280
#define VIEWPORT_HEIGHT 720

constexpr uint32_t c_max_particle_count = 200'000;
constexpr uint32_t c_max_particle_neighbors = 100;

struct Mouse
{
	Mouse()
	{
		m_prevX = 0.0;
		m_prevY = 0.0;
		m_currX = 0.0;
		m_currY = 0.0;
	}

	double m_prevX;
	double m_prevY;
	double m_currX;
	double m_currY;
};

enum NEIGHBORS_COMPUTE_METHOD
{
	EGGMAN,
	SONIC,
	SONIC_BOOM
};

struct UI
{
	int m_draw_mode; // 0 = wireframe, 1 = shaded
	float m_particle_radius;
	bool m_record;
};