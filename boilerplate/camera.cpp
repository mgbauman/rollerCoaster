#include "camera.h"

using namespace glm;

mat4 rotateAbout(vec3 axis, float radians)
{
	mat4 matrix;

	matrix[0][0] = cos(radians) + axis.x*axis.x*(1-cos(radians));
	matrix[1][0] = axis.x*axis.y*(1-cos(radians))-axis.z*sin(radians);
	matrix[2][0] = axis.x*axis.z*(1-cos(radians)) + axis.y*sin(radians);

	matrix[0][1] = axis.y*axis.x*(1-cos(radians)) + axis.z*sin(radians);
	matrix[1][1] = cos(radians) + axis.y*axis.y*(1-cos(radians));
	matrix[2][1] = axis.y*axis.z*(1-cos(radians)) - axis.x*sin(radians);

	matrix[0][2] = axis.z*axis.x*(1-cos(radians)) - axis.y*sin(radians);
	matrix[1][2] = axis.z*axis.y*(1-cos(radians)) + axis.x*sin(radians);
	matrix[2][2] = cos(radians) + axis.z*axis.z*(1-cos(radians));

	return matrix;
}

void Camera::trackballUp(float radians)
{
	mat4 rotation = rotateAbout(right, -radians);

	vec4 newPos = rotation*vec4(pos, 1);
	pos = vec3(newPos.x, newPos.y, newPos.z);

	vec4 newUp = rotation*vec4(up, 1);
	up = normalize(vec3(newUp.x, newUp.y, newUp.z));

	vec4 newDir = rotation*vec4(dir, 1);
	dir = normalize(vec3(newDir.x, newDir.y, newDir.z));
}

void Camera::trackballRight(float radians)
{
	mat4 rotation = rotateAbout(up, radians);

	vec4 newPos = rotation*vec4(pos, 1);
	pos = vec3(newPos.x, newPos.y, newPos.z);

	vec4 newRight = rotation*vec4(right, 1);
	right = normalize(vec3(newRight.x, newRight.y, newRight.z));

	vec4 newDir = rotation*vec4(dir, 1);
	dir = normalize(vec3(newDir.x, newDir.y, newDir.z));
}

void Camera::zoom(float factor)
{
	pos = -dir*length(pos)*factor;
}

mat4 Camera::getMatrix()
{
	mat4 cameraRotation = mat4(
			vec4(right, 0),
			vec4(up, 0),
			vec4(-dir, 0),
			vec4(0, 0, 0, 1));

	mat4 translation = mat4 (
			vec4(1, 0, 0, 0),
			vec4(0, 1, 0, 0),
			vec4(0, 0, 1, 0),
			vec4(-pos, 1));

	return transpose(cameraRotation)*translation;
}