#pragma once
#include "GLCommon.h"
#include "utils.h"
#include "material.hpp"
#include "shader.hpp"

struct AABB
{
	float min_x;
	float max_x;
	float min_y;
	float max_y;
	float min_z;
	float max_z;
};

struct Vertex
{
	void set_position(float const & x, float const & y, float const & z)
	{
		m_position[0] = x;
		m_position[1] = y;
		m_position[2] = z;
	}

	void set_normal(float const& nx, float const& ny, float const& nz)
	{
		m_normal[0] = nx;
		m_normal[1] = ny;
		m_normal[2] = nz;
	}

	bool operator==(const Vertex & other) const
	{
		float eps = 0.001f;
		bool cond_position = glm::all(glm::epsilonEqual(other.m_position, m_position, eps));
		bool cond_normal = glm::all(glm::epsilonEqual(other.m_normal, m_normal, eps));
		return cond_position && cond_normal;
	}

	glm::vec3 m_position;
	glm::vec3 m_normal;
};

struct CPUGeometry
{
	std::vector<Vertex> m_vertices;
	std::vector<uint32_t> m_triangleIndices;
};

struct GPUGeometry
{
	GLuint m_vao;
	GLuint m_pos_vbo;
	GLuint m_normal_vbo;
	GLuint m_instance_vbo;
	GLuint m_ibo;

	~GPUGeometry()
	{
		glDeleteVertexArrays(1, &m_vao);
		glDeleteBuffers(1, &m_pos_vbo);
		glDeleteBuffers(1, &m_normal_vbo);
		glDeleteBuffers(1, &m_instance_vbo);
		glDeleteBuffers(1, &m_ibo);
	}
};

class Mesh
{
public:
	Mesh(std::string const& iFilePath);
	void create_GPU_objects();
	void draw(std::shared_ptr<Shader> iShader, bool use_sonic_boom = false);
	void deactivate_instance_rendering();
	void set_instance_rendering(std::vector<glm::vec3> const & iPositions, float const& iParticleRadius);
	AABB compute_axis_aligned_bounding_box();
	void link_to_pos_SSBO(GLuint ssbo);

	Material m_material;
	CPUGeometry m_cpu_geometry;
	GPUGeometry m_gpu_geometry;
	bool m_instance_rendering;
	uint32_t m_instance_count;
	glm::mat4 m_model;
};