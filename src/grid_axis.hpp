#ifndef GRID_AXIS_HPP
#define GRID_AXIS_HPP

#include "GLCommon.h"
#include "shader.hpp"

class GridAxis
{
	public:

		GridAxis(int gridDim = 30);
		~GridAxis();
		void draw(glm::mat4 view, glm::mat4 projection);

	private:

		// Grid data
		GLuint vaoG;
		GLuint vboG;
		GLuint eboG;
		Shader gridShader;
		float* grid;
		int* indices;
		int dim;
		int nbPoints;
		int nbIndices;

		// Axis Data
		GLuint vaoA;
		GLuint vboA;
		Shader axisShader;
		float* axis;
};

#endif
