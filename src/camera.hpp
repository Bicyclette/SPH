#pragma once

#include "utils.h"
#include "GLCommon.h"

class Camera
{
public:
	static const float m_move_dt;

public:
	Camera(glm::vec3 iPos, glm::vec3 iLookAt, glm::vec3 iUp, int iViewportWidth, int iViewportHeight);
	void rotateView(struct Mouse const & iMouse);
	void zoomView(float direction);
	void panView(struct Mouse const& iMouse);
	void updateViewMatrix(glm::vec3 iPos, glm::vec3 iLookAt);
	void updateProjectionMatrix(int iViewportWidth, int iViewportHeight);
	glm::vec3 getViewDirection();
	glm::vec3 getRightVector();
	glm::mat4 getProjectionMatrix();
	glm::mat4 getViewMatrix();
	glm::vec3 getPosition();
	float getNear();
	float getFar();

private:
	float m_fov;
	float m_near;
	float m_far;
	float m_aspectRatio;
	int m_viewport_width;
	int m_viewport_height;
	glm::vec3 m_position;
	glm::vec3 m_lookAt;
	glm::vec3 m_up;
	glm::mat4 m_proj;
	glm::mat4 m_view;
};