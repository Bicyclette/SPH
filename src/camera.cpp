#include "camera.hpp"

const float Camera::m_move_dt = 0.15f;

Camera::Camera(glm::vec3 iPos, glm::vec3 iLookAt, glm::vec3 iUp, int iViewportWidth, int iViewportHeight)
{
	m_fov = 45.0f;
	m_near = 0.1f;
	m_far = 1000.0f;
	m_position = iPos;
	m_lookAt = iLookAt;
	m_up = iUp;
	m_viewport_width = iViewportWidth;
	m_viewport_height = iViewportHeight;
	m_aspectRatio = static_cast<float>(m_viewport_width) / static_cast<float>(m_viewport_height);
	m_view = glm::lookAt(m_position, m_lookAt, m_up);
	m_proj = glm::perspective(m_fov, m_aspectRatio, m_near, m_far);
}

void Camera::rotateView(struct Mouse const & iMouse)
{
	glm::vec4 eye = glm::vec4(m_position - m_lookAt, 1.0f);
	glm::vec4 center = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
	glm::vec4 right = glm::vec4(getRightVector(), 1.0f);

	float dtx = static_cast<float>((2 * M_PI) / m_viewport_width);
	float dty = static_cast<float>(M_PI / m_viewport_height);
	float azimuth = static_cast<float>((iMouse.m_prevX - iMouse.m_currX) * dtx);
	float zenith = static_cast<float>((iMouse.m_prevY - iMouse.m_currY) * dty);

	glm::mat4 R_horizontal = glm::rotate(glm::mat4(1.0f), azimuth, m_up);
	glm::mat4 R_vertical = glm::rotate(glm::mat4(1.0f), zenith, getRightVector());
	
	eye = R_horizontal * eye;
	m_position = R_vertical * eye;

	// update up vector
	right = R_horizontal * right;
	right = R_vertical * right;
	m_up = glm::normalize(glm::cross(glm::vec3(right), -m_position));

	// translate to actual position
	m_position += m_lookAt;

	updateViewMatrix(m_position, m_lookAt);
}

void Camera::zoomView(float direction)
{
	glm::vec3 new_position = m_position + (direction * m_move_dt * 5.0f) * glm::normalize(getViewDirection());
	float d = glm::length(m_lookAt - new_position);
	if (d > (m_move_dt * 2.0f))
	{
		m_position = new_position;
		m_view = glm::lookAt(m_position, m_lookAt, m_up);
	}
}

void Camera::panView(struct Mouse const& iMouse)
{
	float dtx = static_cast<float>((iMouse.m_prevX - iMouse.m_currX) * m_move_dt);
	float dty = static_cast<float>(-(iMouse.m_prevY - iMouse.m_currY) * m_move_dt);
	glm::vec3 right = glm::normalize(getRightVector());

	m_position += right * dtx;
	m_position += m_up * dty;
	m_lookAt += right * dtx;
	m_lookAt += m_up * dty;

	m_view = glm::lookAt(m_position, m_lookAt, m_up);
}

void Camera::updateViewMatrix(glm::vec3 iPos, glm::vec3 iLookAt)
{
	m_position = iPos;
	m_lookAt = iLookAt;
	m_view = glm::lookAt(m_position, m_lookAt, m_up);
}

void Camera::updateProjectionMatrix(int iViewportWidth, int iViewportHeight)
{
	m_viewport_width = iViewportWidth;
	m_viewport_height = iViewportHeight;
	m_aspectRatio = static_cast<float>(m_viewport_width) / static_cast<float>(m_viewport_height);
	m_proj = glm::perspective(m_fov, m_aspectRatio, m_near, m_far);
}

glm::vec3 Camera::getViewDirection()
{
	return -glm::transpose(m_view)[2];
}

glm::vec3 Camera::getRightVector()
{
	return glm::transpose(m_view)[0];
}

glm::mat4 Camera::getProjectionMatrix()
{
	return m_proj;
}

glm::mat4 Camera::getViewMatrix()
{
	return m_view;
}

glm::vec3 Camera::getPosition()
{
	return m_position;
}

float Camera::getNear()
{
	return m_near;
}

float Camera::getFar()
{
	return m_far;
}

glm::ivec2 Camera::getViewportDimensions()
{
	return glm::ivec2(m_viewport_width, m_viewport_height);
}