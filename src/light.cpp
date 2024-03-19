#include "light.hpp"

const glm::vec3 PointLight::attenuation[12] =
{
		glm::vec3(1.0f, 0.7f, 1.8f),
		glm::vec3(1.0f, 0.35f, 0.44f),
		glm::vec3(1.0f, 0.22f, 0.20f),
		glm::vec3(1.0f, 0.14f, 0.07f),
		glm::vec3(1.0f, 0.09f, 0.032f),
		glm::vec3(1.0f, 0.07f, 0.017f),
		glm::vec3(1.0f, 0.045f, 0.0075f),
		glm::vec3(1.0f, 0.027f, 0.0028f),
		glm::vec3(1.0f, 0.022f, 0.0019f),
		glm::vec3(1.0f, 0.014f, 0.0007f),
		glm::vec3(1.0f, 0.007f, 0.0002f),
		glm::vec3(1.0f, 0.0014f, 0.00007f)
};