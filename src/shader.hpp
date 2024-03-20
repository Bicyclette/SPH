#pragma once

#include "GLCommon.h"
#include "utils.h"

class Shader
{
public:
	Shader() = delete;
	Shader(std::string const& iCompute);
	Shader(std::string const & iVertex, std::string const & iFragment);
	Shader(std::string const & iVertex, std::string const & iFragment, std::string const & iGeometry);
	~Shader();
	void checkCompileError(GLuint const & iStage, GLenum iType, std::string const & path = "");
	bool checkLinkError();
	void set(std::string const & iUniform, bool const & iVal);
	void set(std::string const & iUniform, int const & iVal);
	void set(std::string const & iUniform, float const & iVal);
	void set(std::string const & iUniform, glm::vec2 const & iVec);
	void set(std::string const & iUniform, glm::vec3 const & iVec);
	void set(std::string const & iUniform, glm::ivec3 const & iVec);
	void set(std::string const & iUniform, glm::mat3 const & iMat);
	void set(std::string const & iUniform, glm::mat4 const & iMat);
	void use();

private:
	GLuint m_program;
};