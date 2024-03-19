#include "mesh.hpp"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

namespace std
{
	template<> struct hash<Vertex>
	{
		size_t operator()(Vertex const& vertex) const
		{
			return ((hash<glm::vec3>()(vertex.m_position) ^ (hash<glm::vec3>()(vertex.m_normal) << 1)) >> 1);
		}
	};
}

Mesh::Mesh(std::string const& iFilePath)
{
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string warn;
	std::string err;

	if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, iFilePath.c_str()))
	{
		std::cerr << "Error loading obj file \"" << iFilePath << "\" : " << err << std::endl;
		std::exit(-1);
	}

	std::unordered_map<Vertex, uint32_t> unique_vertices;
	for (const auto& shape : shapes)
	{
		for (const auto& index : shape.mesh.indices)
		{
			// position
			float px = attrib.vertices[3 * index.vertex_index];
			float py = attrib.vertices[3 * index.vertex_index + 1];
			float pz = attrib.vertices[3 * index.vertex_index + 2];
			// normal
			float nx = attrib.normals[3 * index.normal_index];
			float ny = attrib.normals[3 * index.normal_index + 1];
			float nz = attrib.normals[3 * index.normal_index + 2];
			
			Vertex vertex;
			vertex.set_position(px, py, pz);
			vertex.set_normal(nx, ny, nz);

			if (unique_vertices.count(vertex) == 0)
			{
				unique_vertices[vertex] = static_cast<uint32_t>(m_cpu_geometry.m_vertices.size());
				m_cpu_geometry.m_vertices.push_back(vertex);
			}
			m_cpu_geometry.m_triangleIndices.push_back(unique_vertices[vertex]);
		}
	}

	// other data
	m_model = glm::mat4(1.0f);
	m_instance_rendering = false;
	m_instance_count = 0;
	create_GPU_objects();
}

void Mesh::create_GPU_objects()
{
	// CREATE DATA ARRAYS
	std::vector<float> positionBuffer;
	std::vector<float> normalBuffer;

	for (Vertex const& v : m_cpu_geometry.m_vertices)
	{
		positionBuffer.push_back(v.m_position[0]);
		positionBuffer.push_back(v.m_position[1]);
		positionBuffer.push_back(v.m_position[2]);

		normalBuffer.push_back(v.m_normal[0]);
		normalBuffer.push_back(v.m_normal[1]);
		normalBuffer.push_back(v.m_normal[2]);
	}

	// VAO
	glGenVertexArrays(1, &m_gpu_geometry.m_vao);
	glBindVertexArray(m_gpu_geometry.m_vao);

	// POSITION VBO
	glGenBuffers(1, &m_gpu_geometry.m_pos_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, m_gpu_geometry.m_pos_vbo);
	glBufferData(GL_ARRAY_BUFFER, positionBuffer.size() * sizeof(float), positionBuffer.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);

	// NORMAL VBO
	glGenBuffers(1, &m_gpu_geometry.m_normal_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, m_gpu_geometry.m_normal_vbo);
	glBufferData(GL_ARRAY_BUFFER, normalBuffer.size() * sizeof(float), normalBuffer.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);

	// INSTANCE VBO
	glGenBuffers(1, &m_gpu_geometry.m_instance_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, m_gpu_geometry.m_instance_vbo);
	glBufferData(GL_ARRAY_BUFFER, c_max_particle_count * sizeof(glm::mat4), nullptr, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)0);
	glVertexAttribDivisor(2, 1);
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(sizeof(glm::vec4)));
	glVertexAttribDivisor(3, 1);
	glEnableVertexAttribArray(4);
	glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(2 * sizeof(glm::vec4)));
	glVertexAttribDivisor(4, 1);
	glEnableVertexAttribArray(5);
	glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(glm::vec4), (void*)(3 * sizeof(glm::vec4)));
	glVertexAttribDivisor(5, 1);

	// ELEMENT BUFFER
	glGenBuffers(1, &m_gpu_geometry.m_ibo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_gpu_geometry.m_ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_cpu_geometry.m_triangleIndices.size() * sizeof(unsigned int), m_cpu_geometry.m_triangleIndices.data(), GL_STATIC_DRAW);

	glBindVertexArray(0);
}

void Mesh::draw(std::shared_ptr<Shader> iShader, bool use_sonic_boom)
{
	iShader->use();

	glBindVertexArray(m_gpu_geometry.m_vao);
	if (!m_instance_rendering)
	{
		glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(m_cpu_geometry.m_triangleIndices.size()), GL_UNSIGNED_INT, 0);
	}
	else
	{
		glDrawElementsInstanced(GL_TRIANGLES, static_cast<GLsizei>(m_cpu_geometry.m_triangleIndices.size()), GL_UNSIGNED_INT, 0, m_instance_count);
	}
}

AABB Mesh::compute_axis_aligned_bounding_box()
{
	std::vector<Vertex> const& vertices = m_cpu_geometry.m_vertices;
	
	float x = vertices[0].m_position.x;
	float y = vertices[0].m_position.y;
	float z = vertices[0].m_position.z;
	
	AABB res;
	res.min_x = x;
	res.max_x = x;
	res.min_y = y;
	res.max_y = y;
	res.min_z = z;
	res.max_z = z;

	for (size_t i = 1; i < m_cpu_geometry.m_vertices.size(); ++i)
	{
		glm::vec3 vertex = vertices[i].m_position;

		if (res.min_x > vertex.x) { res.min_x = vertex.x; }
		else if (res.max_x < vertex.x) { res.max_x = vertex.x; }

		if (res.min_y > vertex.y) { res.min_y = vertex.y; }
		else if (res.max_y < vertex.y) { res.max_y = vertex.y; }

		if (res.min_z > vertex.z) { res.min_z = vertex.z; }
		else if (res.max_z < vertex.z) { res.max_z = vertex.z; }
	}

	return res;
}

void Mesh::deactivate_instance_rendering()
{
	m_instance_rendering = false;
	m_instance_count = 0;
}

void Mesh::set_instance_rendering(std::vector<glm::vec3> const& iPositions, float const & iParticleRadius)
{
	m_instance_rendering = true;
	m_instance_count = iPositions.size();

	std::vector<glm::mat4> instance_model;
	instance_model.resize(iPositions.size());

	glm::vec3 const scale(iParticleRadius);
	
	for (size_t i = 0; i < iPositions.size(); ++i)
	{
		glm::mat4 model(1.0f);
		model = glm::translate(model, iPositions[i]);
		model = glm::scale(model, scale);
		instance_model[i] = model;
	}

	// INSTANCE VBO
	glBindBuffer(GL_ARRAY_BUFFER, m_gpu_geometry.m_instance_vbo);
	glBufferSubData(GL_ARRAY_BUFFER, 0, instance_model.size() * sizeof(glm::mat4), instance_model.data());
}

void Mesh::link_to_pos_SSBO(GLuint ssbo)
{
	glBindVertexArray(m_gpu_geometry.m_vao);
	glBindBuffer(GL_ARRAY_BUFFER, ssbo);
}