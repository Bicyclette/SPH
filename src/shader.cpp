#include "shader.hpp"

Shader::Shader(std::string const& iCompute)
{
	GLuint cShader = glCreateShader(GL_COMPUTE_SHADER);
	std::string cShaderString = file2String(iCompute);
	const GLchar* cShaderSource = reinterpret_cast<const GLchar*>(cShaderString.c_str());
	glShaderSource(cShader, 1, &cShaderSource, NULL);
	glCompileShader(cShader);

	// Check for errors
	checkCompileError(cShader, GL_COMPUTE_SHADER, iCompute);

	// create program
	m_program = glCreateProgram();
	glAttachShader(m_program, cShader);
	glLinkProgram(m_program);
	if (!checkLinkError())
	{
		glDeleteShader(cShader);
	}
	glDetachShader(m_program, cShader);
}

Shader::Shader(std::string const & iVertex, std::string const & iFragment)
{
	GLuint vShader = glCreateShader(GL_VERTEX_SHADER);
	std::string vShaderString = file2String(iVertex);
	const GLchar* vShaderSource = reinterpret_cast<const GLchar*>(vShaderString.c_str());
	glShaderSource(vShader, 1, &vShaderSource, NULL);
	glCompileShader(vShader);

	GLuint fShader = glCreateShader(GL_FRAGMENT_SHADER);
	std::string fShaderString = file2String(iFragment);
	const GLchar* fShaderSource = reinterpret_cast<const GLchar*>(fShaderString.c_str());
	glShaderSource(fShader, 1, &fShaderSource, NULL);
	glCompileShader(fShader);

	// Check for errors
	checkCompileError(vShader, GL_VERTEX_SHADER, iVertex);
	checkCompileError(fShader, GL_FRAGMENT_SHADER, iFragment);

	// create program
	m_program = glCreateProgram();
	glAttachShader(m_program, vShader);
	glAttachShader(m_program, fShader);
	glLinkProgram(m_program);
	if (!checkLinkError())
	{
		glDeleteShader(vShader);
		glDeleteShader(fShader);
	}
	glDetachShader(m_program, vShader);
	glDetachShader(m_program, fShader);
}

Shader::Shader(std::string const & iVertex, std::string const & iFragment, std::string const & iGeometry)
{
	GLuint vShader = glCreateShader(GL_VERTEX_SHADER);
	std::string vShaderString = file2String(iVertex);
	const GLchar* vShaderSource = reinterpret_cast<const GLchar*>(vShaderString.c_str());
	glShaderSource(vShader, 1, &vShaderSource, NULL);
	glCompileShader(vShader);

	GLuint gShader = glCreateShader(GL_GEOMETRY_SHADER);
	std::string gShaderString = file2String(iGeometry);
	const GLchar* gShaderSource = reinterpret_cast<const GLchar*>(gShaderString.c_str());
	glShaderSource(gShader, 1, &gShaderSource, NULL);
	glCompileShader(gShader);

	GLuint fShader = glCreateShader(GL_FRAGMENT_SHADER);
	std::string fShaderString = file2String(iFragment);
	const GLchar* fShaderSource = reinterpret_cast<const GLchar*>(fShaderString.c_str());
	glShaderSource(fShader, 1, &fShaderSource, NULL);
	glCompileShader(fShader);
	
	// Check for errors
	checkCompileError(vShader, GL_VERTEX_SHADER, iVertex);
	checkCompileError(fShader, GL_GEOMETRY_SHADER, iFragment);
	checkCompileError(fShader, GL_FRAGMENT_SHADER, iGeometry);

	// create program
	m_program = glCreateProgram();
	glAttachShader(m_program, vShader);
	glAttachShader(m_program, gShader);
	glAttachShader(m_program, fShader);
	glLinkProgram(m_program);
	if (!checkLinkError())
	{
		glDeleteShader(vShader);
		glDeleteShader(gShader);
		glDeleteShader(fShader);
	}
	glDetachShader(m_program, vShader);
	glDetachShader(m_program, gShader);
	glDetachShader(m_program, fShader);
}

Shader::~Shader()
{
	glDeleteProgram(m_program);
}

void Shader::checkCompileError(GLuint const & iStage, GLenum iType, std::string const & path)
{
	int success;
	int logLength;
	std::unique_ptr<char> log;

	glGetShaderiv(iStage, GL_COMPILE_STATUS, &success);
	if (success == GL_FALSE)
	{
		glGetShaderiv(iStage, GL_INFO_LOG_LENGTH, &logLength);
		log = std::make_unique<char>(logLength);
		glGetShaderInfoLog(iStage, logLength, nullptr, log.get());
		if(iType == GL_VERTEX_SHADER)
		{
			std::cerr << "Error while compiling the vertex shader " << path << ": " << log.get() << std::endl;
		}
		else if(iType == GL_FRAGMENT_SHADER)
		{
			std::cerr << "Error while compiling the fragment shader " << path << ": " << log.get() << std::endl;
		}
		else if(iType == GL_GEOMETRY_SHADER)
		{
			std::cerr << "Error while compiling the geometry shader " << path << ": " << log.get() << std::endl;
		}
		else if (iType == GL_COMPUTE_SHADER)
		{
			std::cerr << "Error while compiling the compute shader " << path << ": " << log.get() << std::endl;
		}
		glDeleteShader(iStage);
	}
}

bool Shader::checkLinkError()
{
	int success;
	int logLength;
	std::unique_ptr<char> log;

	glGetProgramiv(m_program, GL_LINK_STATUS, &success);
	if (success == GL_FALSE)
	{
		glGetProgramiv(m_program, GL_INFO_LOG_LENGTH, &logLength);
		log = std::make_unique<char>(logLength);
		glGetProgramInfoLog(m_program, logLength, nullptr, log.get());
		std::cerr << "Error while linking shaders into a program : " << log.get() << std::endl;
		return false;
	}
	return true;
}

void Shader::set(std::string const & iUniform, bool const & iVal)
{
	GLint location = glGetUniformLocation(m_program, iUniform.c_str());
	glUniform1i(location, iVal);
}

void Shader::set(std::string const & iUniform, int const & iVal)
{
	GLint location = glGetUniformLocation(m_program, iUniform.c_str());
	glUniform1i(location, iVal);
}

void Shader::set(std::string const & iUniform, float const & iVal)
{
	GLint location = glGetUniformLocation(m_program, iUniform.c_str());
	glUniform1f(location, iVal);
}

void Shader::set(std::string const & iUniform, glm::vec2 const & iVec)
{
	GLint location = glGetUniformLocation(m_program, iUniform.c_str());
	glUniform2fv(location, 1, glm::value_ptr(iVec));
}

void Shader::set(std::string const & iUniform, glm::vec3 const & iVec)
{
	GLint location = glGetUniformLocation(m_program, iUniform.c_str());
	glUniform3fv(location, 1, glm::value_ptr(iVec));
}

void Shader::set(std::string const & iUniform, glm::ivec3 const& iVec)
{
	GLint location = glGetUniformLocation(m_program, iUniform.c_str());
	glUniform3iv(location, 1, glm::value_ptr(iVec));
}

void Shader::set(std::string const & iUniform, glm::mat3 const & iMat)
{
	GLint location = glGetUniformLocation(m_program, iUniform.c_str());
	glUniformMatrix3fv(location, 1, GL_FALSE, glm::value_ptr(iMat));
}

void Shader::set(std::string const & iUniform, glm::mat4 const & iMat)
{
	GLint location = glGetUniformLocation(m_program, iUniform.c_str());
	glUniformMatrix4fv(location, 1, GL_FALSE, glm::value_ptr(iMat));
}

void Shader::use()
{
	glUseProgram(m_program);
}
