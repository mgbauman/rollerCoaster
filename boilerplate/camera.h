#ifndef CAMERA_H
#define CAMERA_H


#include "glm/glm.hpp"
#include <cstdio>

//using namespace glm;

class Camera{
public:
	glm::vec3 dir;
	glm::vec3 right;
	glm::vec3 up;

	glm::vec3 pos;

	Camera():	dir(glm::vec3(0, 0, -1)), right(glm::vec3(1, 0, 0)), up(glm::vec3(0, 1, 0)),
				pos(glm::vec3(0, 0, 1.f)){}

	Camera(glm::vec3 _dir, glm::vec3 _pos):dir(_dir), pos(_pos){
		right = normalize(cross(dir, glm::vec3(0, 1, 0)));
		up = normalize(cross(right, dir));
	}

	void trackballUp(float radians);
	void trackballRight(float radians);
	void zoom(float factor);
	glm::mat4 getMatrix();
};

#endif