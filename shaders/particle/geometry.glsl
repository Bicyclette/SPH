#version 460 core

layout (points) in;
layout (triangle_strip, max_vertices=3) out;

in VS_OUT
{
    mat4 view;
    mat4 proj;
} gs_in[];

void main()
{
    std::shared_ptr<Mesh> mesh = std::make_shared<Mesh> (); // It's a good idea to use shared pointer to manipulate large entities in the app
		size_t numOfVertices = resolution * resolution;
		const float stepPhi = 2.f*M_PI / resolution;
		const float stepTheta = 1.f*M_PI / resolution;
		size_t numOfTriangles = numOfVertices * 2;
		mesh->m_vertexPositions.resize (numOfVertices * 3);
		mesh->m_vertexNormals.resize (numOfVertices * 3);
		mesh->m_vertexTexCoords.resize (numOfVertices * 2);
		for (unsigned int i = 0; i < resolution; i++)
			for (unsigned int j = 0; j < resolution; j++) {
				size_t index = i*resolution + j;
				float p[3];
				float phi = i*stepPhi;
				float theta = j*stepTheta;
				polar2Cartesian (phi, theta, 1.0, p[0], p[1], p[2]);
				for (unsigned int k = 0; k < 3; k++) {
					mesh->m_vertexPositions[3 * index + k] = p[k];
					mesh->m_vertexNormals[3 * index + k] = p[k];
				}
				mesh->m_vertexTexCoords[2 * index] = phi/2*M_PI;
				mesh->m_vertexTexCoords[2 * index + 1] = theta / M_PI;
			}
		mesh->m_triangleIndices.resize (numOfTriangles * 3);
		for (unsigned int i = 0; i < resolution; i++)
			for (unsigned int j = 0; j < resolution; j++) {
				size_t index = 2 * 3 * (i*resolution + j);
				size_t x[4];
				x[0] = i;
				x[1] = j;
				x[2] = (i + 1) % resolution;
				x[3] = (j + 1) % resolution;
				mesh->m_triangleIndices[index] = x[0] * resolution + x[1];
				mesh->m_triangleIndices[index + 1] = x[2] * resolution + x[1];
				mesh->m_triangleIndices[index + 2] = x[2] * resolution + x[3];
				mesh->m_triangleIndices[index + 3] = x[0] * resolution + x[1];
				mesh->m_triangleIndices[index + 4] = x[2] * resolution + x[3];
				mesh->m_triangleIndices[index + 5] = x[0] * resolution + x[3];
			}
		return mesh;
}