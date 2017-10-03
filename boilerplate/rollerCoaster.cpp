// ==========================================================================
// Barebones OpenGL Core Profile Boilerplate
//    using the GLFW windowing system (http://www.glfw.org)
//
// Loosely based on
//  - Chris Wellons' example (https://github.com/skeeto/opengl-demo) and
//  - Camilla Berglund's example (http://www.glfw.org/docs/latest/quick.html)
//
// Author:  Sonny Chan, University of Calgary
// Date:    December 2015
// ==========================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <algorithm>
#include <vector>
#include <cstdlib>

#include "glm/glm.hpp"
#include "glm/mat4x4.hpp"
#include "glm/gtc/matrix_transform.hpp"

#include "Vec3f.h"
#include "Vec3f_FileIO.h"

// specify that we want the OpenGL core profile before including GLFW headers
#include "glad/glad.h"
#include <GLFW/glfw3.h>

#include "camera.h"
#define PI 3.14159265359

using namespace std;
using namespace glm;

//Forward definitions
bool CheckGLErrors(string location);
void QueryGLVersion();
string LoadSource(const string &filename);
GLuint CompileShader(GLenum shaderType, const string &source);
GLuint LinkProgram(GLuint vertexShader, GLuint fragmentShader);
vector<vec3> loadInVecsFromFile();

vec2 mousePos;
bool leftmousePressed = false;
bool rightmousePressed = false;
//int t = 0;
float g = 0.000981*2;
Camera* activeCamera;

Camera cam = Camera(vec3(0, 0, -1), vec3(0, 0, 7));
// Remember to start this at false for final submission and demo
bool animate = true;
bool drawTrack = false;
bool firstPerson = false;
GLFWwindow* window = 0;

mat4 winRatio = mat4(1.f);

// --------------------------------------------------------------------------
// GLFW callback functions

// reports GLFW errors
void ErrorCallback(int error, const char* description)
{
    cout << "GLFW ERROR " << error << ":" << endl;
    cout << description << endl;
}

// handles keyboard input events
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
        
    if (key == GLFW_KEY_SPACE && action == GLFW_PRESS){
		animate = !animate;
	}
	if (key == GLFW_KEY_T && action == GLFW_PRESS){
		drawTrack = !drawTrack;
	}
	if (key == GLFW_KEY_F && action == GLFW_PRESS){
		firstPerson = !firstPerson;
		if ( !firstPerson){
			activeCamera->pos = vec3(0,0,7);
			activeCamera->dir = vec3(0,0,-1);
			activeCamera->up = vec3(0,1,0);
			activeCamera->right = vec3(1,0,0);
		} 
	}
}

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
	if( (action == GLFW_PRESS) || (action == GLFW_RELEASE) ){
		if(button == GLFW_MOUSE_BUTTON_LEFT)
			leftmousePressed = !leftmousePressed;
		else if(button == GLFW_MOUSE_BUTTON_RIGHT)
			rightmousePressed = !rightmousePressed;
	}
}

void mousePosCallback(GLFWwindow* window, double xpos, double ypos)
{
	int vp[4];
	glGetIntegerv(GL_VIEWPORT, vp);

	vec2 newPos = vec2(xpos/(double)vp[2], -ypos/(double)vp[3])*2.f - vec2(1.f);

	vec2 diff = newPos - mousePos;
	if(leftmousePressed){
		activeCamera->trackballRight(-diff.x);
		activeCamera->trackballUp(-diff.y);
	}
	else if(rightmousePressed){
		float zoomBase = (diff.y > 0) ? 1.f/2.f : 2.f;

		activeCamera->zoom(pow(zoomBase, abs(diff.y)));
	}

	mousePos = newPos;
}

void resizeCallback(GLFWwindow* window, int width, int height)
{
	int vp[4];
	glGetIntegerv(GL_VIEWPORT, vp);

	glViewport(0, 0, width, height);

	float minDim = float(min(width, height));

	winRatio[0][0] = minDim/float(width);
	winRatio[1][1] = minDim/float(height);
}

//==========================================================================
// TUTORIAL STUFF


//vec2 and vec3 are part of the glm math library. 
//Include in your own project by putting the glm directory in your project, 
//and including glm/glm.hpp as I have at the top of the file.
//"using namespace glm;" will allow you to avoid writing everyting as glm::vec2

struct VertexBuffers{
	enum{ VERTICES=0, NORMALS, INDICES, COUNT};

	GLuint id[COUNT];
};

//Describe the setup of the Vertex Array Object
bool initVAO(GLuint vao, const VertexBuffers& vbo)
{
	glBindVertexArray(vao);		//Set the active Vertex Array

	glEnableVertexAttribArray(0);		//Tell opengl you're using layout attribute 0 (For shader input)
	glBindBuffer( GL_ARRAY_BUFFER, vbo.id[VertexBuffers::VERTICES] );		//Set the active Vertex Buffer
	glVertexAttribPointer(
		0,				//Attribute
		3,				//Size # Components
		GL_FLOAT,	//Type
		GL_FALSE, 	//Normalized?
		sizeof(vec3),	//Stride
		(void*)0			//Offset
		);

	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, vbo.id[VertexBuffers::NORMALS]);
	glVertexAttribPointer(
		1,				//Attribute
		3,				//Size # Components
		GL_FLOAT,	//Type
		GL_FALSE, 	//Normalized?
		sizeof(vec3),	//Stride
		(void*)0			//Offset
		);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo.id[VertexBuffers::INDICES]);

	return !CheckGLErrors("initVAO");		//Check for errors in initialize
}


//Loads buffers with data
bool loadBuffer(const VertexBuffers& vbo, 
				const vector<vec3>& points, 
				const vector<vec3> normals, 
				const vector<unsigned int>& indices)
{
	glBindBuffer(GL_ARRAY_BUFFER, vbo.id[VertexBuffers::VERTICES]);
	glBufferData(
		GL_ARRAY_BUFFER,				//Which buffer you're loading too
		sizeof(vec3)*points.size(),		//Size of data in array (in bytes)
		&points[0],						//Start of array (&points[0] will give you pointer to start of vector)
		GL_STATIC_DRAW					//GL_DYNAMIC_DRAW if you're changing the data often
										//GL_STATIC_DRAW if you're changing seldomly
		);

	glBindBuffer(GL_ARRAY_BUFFER, vbo.id[VertexBuffers::NORMALS]);
	glBufferData(
		GL_ARRAY_BUFFER,				//Which buffer you're loading too
		sizeof(vec3)*normals.size(),	//Size of data in array (in bytes)
		&normals[0],					//Start of array (&points[0] will give you pointer to start of vector)
		GL_STATIC_DRAW					//GL_DYNAMIC_DRAW if you're changing the data often
										//GL_STATIC_DRAW if you're changing seldomly
		);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo.id[VertexBuffers::INDICES]);
	glBufferData(
		GL_ELEMENT_ARRAY_BUFFER,
		sizeof(unsigned int)*indices.size(),
		&indices[0],
		GL_STATIC_DRAW
		);

	return !CheckGLErrors("loadBuffer");	
}

//Compile and link shaders, storing the program ID in shader array
GLuint initShader(string vertexName, string fragmentName)
{	
	string vertexSource = LoadSource(vertexName);		//Put vertex file text into string
	string fragmentSource = LoadSource(fragmentName);		//Put fragment file text into string

	GLuint vertexID = CompileShader(GL_VERTEX_SHADER, vertexSource);
	GLuint fragmentID = CompileShader(GL_FRAGMENT_SHADER, fragmentSource);
	
	return LinkProgram(vertexID, fragmentID);	//Link and store program ID in shader array
}

//Initialization
void initGL()
{
	glEnable(GL_DEPTH_TEST);
	glPointSize(25);
	glDepthFunc(GL_LEQUAL);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glClearColor(0.f, 0.f, 0.f, 0.f);		//Color to clear the screen with (R, G, B, Alpha)
}

bool loadUniforms(GLuint program, mat4 perspective, mat4 modelview)
{
	glUseProgram(program);

	glUniformMatrix4fv(glGetUniformLocation(program, "modelviewMatrix"),
						1,
						false,
						&modelview[0][0]);

	glUniformMatrix4fv(glGetUniformLocation(program, "perspectiveMatrix"),
						1,
						false,
						&perspective[0][0]);

	return !CheckGLErrors("loadUniforms");
}

//Draws buffers to screen
void render(GLuint vao, int startElement, int numElements)
{
	glBindVertexArray(vao);		//Use the LINES vertex array

	glDrawElements(
			GL_TRIANGLES,		//What shape we're drawing	- GL_TRIANGLES, GL_LINES, GL_POINTS, GL_QUADS, GL_TRIANGLE_STRIP
			numElements,		//How many indices
			GL_UNSIGNED_INT,	//Type
			(void*)0			//Offset
			);

	CheckGLErrors("render");
}


//Draws buffers to screen
void renderTrack(GLuint vao, int startElement, int numElements)
{
	glBindVertexArray(vao);

	glDrawElements(
			GL_LINE_LOOP,		//What shape we're drawing	- GL_TRIANGLES, GL_LINES, GL_POINTS, GL_QUADS, GL_TRIANGLE_STRIP
			numElements,		//How many indices
			GL_UNSIGNED_INT,	//Type
			(void*)0			//Offset
			);

	CheckGLErrors("renderTrack");
}

//Draws buffers to screen
void renderSquare(GLuint vao, int startElement, int numElements)
{
	glBindVertexArray(vao);

	glDrawElements(
			GL_TRIANGLES,		//What shape we're drawing	- GL_TRIANGLES, GL_LINES, GL_POINTS, GL_QUADS, GL_TRIANGLE_STRIP
			numElements,		//How many indices
			GL_UNSIGNED_INT,	//Type
			(void*)0			//Offset
			);

	CheckGLErrors("renderSquare");
}



//Draws buffers to screen
void renderBox(GLuint vao, int startElement, int numElements)
{
	glBindVertexArray(vao);

	glDrawElements(
			GL_TRIANGLES,		//What shape we're drawing	- GL_TRIANGLES, GL_LINES, GL_POINTS, GL_QUADS, GL_TRIANGLE_STRIP
			numElements,		//How many indices
			GL_UNSIGNED_INT,	//Type
			(void*)0			//Offset
			);

	CheckGLErrors("renderBox");
}


void generateSquare(vector<vec3>& vertices, vector<vec3>& normals, 
					vector<unsigned int>& indices, float width)
{
	vertices.push_back(vec3(-width*0.5f,-0.06, -width*0.5f));
	vertices.push_back(vec3(width*0.5f,-0.06, -width*0.5f));
	vertices.push_back(vec3(width*0.5f,-0.06, width*0.5f));
	vertices.push_back(vec3(-width*0.5f,-0.06, width*0.5f));

	normals.push_back(vec3(0.8f, 0.8f, 0.8f));
	normals.push_back(vec3(0.8f, 0.8f, 0.8f));
	normals.push_back(vec3(0.8f, 0.8f, 0.8f));
	normals.push_back(vec3(0.8f, 0.8f, 0.8f));

	//First triangle
	indices.push_back(0);
	indices.push_back(1);
	indices.push_back(2);
	//Second triangle
	indices.push_back(2);
	indices.push_back(3);
	indices.push_back(0);
}

GLFWwindow* createGLFWWindow()
{
	// initialize the GLFW windowing system
    if (!glfwInit()) {
        cout << "ERROR: GLFW failed to initialize, TERMINATING" << endl;
        return NULL;
    }
    glfwSetErrorCallback(ErrorCallback);

    // attempt to create a window with an OpenGL 4.1 core profile context
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    window = glfwCreateWindow(512, 512, "OpenGL Example", 0, 0);
    if (!window) {
        cout << "Program failed to create GLFW window, TERMINATING" << endl;
        glfwTerminate();
        return NULL;
    }

    // set keyboard callback function and make our context current (active)
    glfwSetKeyCallback(window, keyCallback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetCursorPosCallback(window, mousePosCallback);
    glfwSetWindowSizeCallback(window, resizeCallback);
    glfwMakeContextCurrent(window);

    return window;
}

void generateSphere(vector<vec3>* positions, vector<vec3>* normals, vector<unsigned int>* indices,
					float r, vec3 center, int uDivisions, int vDivisions)
{
	// udivisions will be theta
	// vdivisions will be phi	
	float uStep = 1.f/(float)(uDivisions-1);
	float vStep = 1.f/(float)(vDivisions-1);

	float u = 0.f;
	
    // Iterate through phi and theta
    for (double phi = 0.; phi < uDivisions; phi ++) // Azimuth [0, 2PI]
    {
		float v = 0.f;
        for (double theta = 0.; theta < vDivisions; theta++) // Elevation [0, PI]
        {
            vec3 point;
            point.x = r * cos(v*2*PI) * sin(u*PI) + center.x;
            point.y = r * sin(v*2*PI) * sin(u*PI) + center.y;
            point.z = r               * cos(u*PI) + center.z;
            
            vec3 normal = normalize(point - center);
            
            positions->push_back(point);
            normals->push_back(normal);
            
            v+=vStep;
        }
        u+=uStep;
    }
   
    for(int i=0; i<uDivisions-1; i++)
	{
		for(int j=0; j<vDivisions -1; j++)
		{
			unsigned int p00 = i*vDivisions+j;
			unsigned int p01 = i*vDivisions+j+1;
			unsigned int p10 = (i+1)*vDivisions + j;
			unsigned int p11 = (i+1)*vDivisions + j + 1;

			indices->push_back(p00);
			indices->push_back(p10);
			indices->push_back(p01);

			indices->push_back(p01);
			indices->push_back(p10);
			indices->push_back(p11);
		}
	}
    
}

void generateCylinder(vector<vec3>* positions, vector<vec3>* normals, vector<unsigned int>* indices,
					float r, vec3 center, int uDivisions, int vDivisions, float h)
{
	// udivisions will be theta
	// vdivisions will be phi	
	float uStep = 1.f/(float)(uDivisions-1);
	float vStep = 1.f/(float)(vDivisions-1);

	float u = 0.f;
	
    // Iterate through phi and theta
    for (double phi = 0.; phi < uDivisions; phi ++) // Azimuth [0, 2PI]
    {
		float v = 0.f;
        for (double theta = 0.; theta < vDivisions; theta++) // Elevation [0, PI]
        {
            vec3 point;
            point.x = r * 2 * cos(v*2*PI) - 0.05;
            point.y = r * 2 * sin(v*2*PI) - 4*r ;
            point.z = -(u*h) + (h/2);
            
            vec3 normal = normalize(point - center);
            normal = vec3 ( 1, 0 ,0 );
            positions->push_back(point);
            normals->push_back(normal);
            
            v+=vStep;
        }
        u+=uStep;
    }
	u = 0.f;
    // Iterate through phi and theta
    for (double phi = 0.; phi < uDivisions; phi ++) // Azimuth [0, 2PI]
    {
		float v = 0.f;
        for (double theta = 0.; theta < vDivisions; theta++) // Elevation [0, PI]
        {
            vec3 point;
            point.x = r * 2 * cos(v*2*PI) + 0.05;
            point.y = r * 2 * sin(v*2*PI) - 4*r ;
            point.z = -(u*h) + (h/2);
            
            vec3 normal = normalize(point - center);
            normal = vec3 ( 1, 0 ,0 );
            positions->push_back(point);
            normals->push_back(normal);
            
            v+=vStep;
        }
        u+=uStep;
    }
   
    for(int i=0; i<(uDivisions*2)-2; i++)
	{
		for(int j=0; j<(vDivisions*2) -2; j++)
		{
			unsigned int p00 = i*vDivisions+j;
			unsigned int p01 = i*vDivisions+j+1;
			unsigned int p10 = (i+1)*vDivisions + j;
			unsigned int p11 = (i+1)*vDivisions + j + 1;

			indices->push_back(p00);
			indices->push_back(p10);
			indices->push_back(p01);

			indices->push_back(p01);
			indices->push_back(p10);
			indices->push_back(p11);
		}
	}
}


void generateBox(vector<vec3>* vertices, vector<vec3>* normals, vector<unsigned int>* indices, float sideLength){
	vec3 nxnyz = vec3(-sideLength, 0, sideLength);
	vec3 nxnynz = vec3(-sideLength, 0, -sideLength);
	vec3 nxynz = vec3(-sideLength, sideLength, -sideLength);
	vec3 nxyz = vec3(-sideLength, sideLength, sideLength);
	
	vec3 xnyz = vec3(sideLength, 0, sideLength);
	vec3 xnynz = vec3(sideLength, 0, -sideLength);
	vec3 xynz = vec3(sideLength, sideLength, -sideLength);
	vec3 xyz = vec3(sideLength, sideLength, sideLength);
	
	vertices->push_back(nxnyz); 
	vertices->push_back(nxnynz);
	vertices->push_back(nxynz);
	vertices->push_back(nxyz);
	
	vertices->push_back(xnyz);
	vertices->push_back(xnynz);
	vertices->push_back(xynz);
	vertices->push_back(xyz);
	
	
	
	for( int i =0; i<1; i++){
		normals->push_back(vec3(1.f,1.f,0.f));
	}
	for( int i =0; i<2; i++){
		normals->push_back(vec3(0,0,1));
	}
	for( int i =0; i<2; i++){
		normals->push_back(vec3(1.f,1.f,0.f));
	}
	for( int i =0; i<2; i++){
		normals->push_back(vec3(0,0,1));
	}
	for( int i =0; i<1; i++){
		normals->push_back(vec3(1.f,1.f,0.f));
	}
	
	
	//TOP
	indices->push_back(3);indices->push_back(2);indices->push_back(6);
	indices->push_back(6);indices->push_back(7);indices->push_back(3);
	
	//RIGHT
	indices->push_back(5);indices->push_back(6);indices->push_back(7);
	indices->push_back(4);indices->push_back(5);indices->push_back(7);
	
	//BACK
	indices->push_back(1);indices->push_back(2);indices->push_back(6);
	indices->push_back(1);indices->push_back(5);indices->push_back(6);

	//LEFT
	indices->push_back(1);indices->push_back(2);indices->push_back(3);
	indices->push_back(0);indices->push_back(1);indices->push_back(3);

	//FRONT
	indices->push_back(0);indices->push_back(4);indices->push_back(7);
	indices->push_back(0);indices->push_back(7);indices->push_back(3);

	//BOTTOM
	indices->push_back(0);indices->push_back(1);indices->push_back(5);
	indices->push_back(0);indices->push_back(5);indices->push_back(4);
}

// Subdivide the track @subdivison times using Chaikin Subdivision
std::vector<vec3> subdivideTrack(std::vector<vec3> controlPoints, vector<vec3>* normals, vector<unsigned int>* indices, int subdivisions){
	for( int i =0; i < subdivisions; i++){
		
		std::vector<vec3> splitPoints;
		std::vector<vec3> repoPoints;
	
		//split
		for( unsigned int j =0; j< controlPoints.size(); j++){
			vec3 midPoint = (controlPoints[(j+1)%controlPoints.size()] + controlPoints[j])/2.0f;
			splitPoints.push_back(controlPoints[j]);
			splitPoints.push_back(midPoint);
		}
		
		//reposition
		for( unsigned int j =0; j< splitPoints.size(); j++){
			vec3 newPoint = (splitPoints[(j+1)%splitPoints.size()] + splitPoints[j])/2.0f;
			repoPoints.push_back(newPoint);
		}
		
		controlPoints = repoPoints;
	}
	
	for (unsigned int i = 0; i < controlPoints.size(); i++){
		indices->push_back(i);
		normals->push_back(vec3(0,0,1.0));
	}

	return controlPoints;
}

float magnitude(vec3 v){
	return sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
}

vec3 arcLengthParamaterize(vec3 center, int &i, vector<vec3> p, float s){
	// This will ensure that the array stays in bounds and loops once it reaches the upper bound
	int ip = (i+1)%p.size();
	if( magnitude(p[ip] - center) > s){
		vec3 next = center + s * ((p[ip] - center) / (magnitude(p[ip] - center)));
		return next;
	} else {
		float sNew = magnitude(p[ip] - center);
		i= (i+1)%p.size();
		ip= (i+1)%p.size();
		while(sNew + magnitude(p[ip] - p[i]) < s){
			sNew = sNew + magnitude(p[ip] - p[i]);
			i= (i+1)%p.size();
			ip= (i+1)%p.size();
		}
		
		return p[i]+ (s-sNew) * ((p[ip] - p[i]) / (magnitude(p[ip] - p[i])));
	}
}

vec3 getMaxHeight(vector<vec3> positions, int &i){
	vec3 max = positions[0];
	vec3 temp;
	for( unsigned int j = 0; j < positions.size(); j++){
		if ( positions[j].y >= max.y ){
			max = positions[j];
			i = j; 
		}
	}
	return max;
}

vec3 getMinHeight(vector<vec3> positions, int &i){
	vec3 min = positions[0];
	vec3 temp;
	for( unsigned int j = 0; j < positions.size(); j++){
		if ( positions[j].y < min.y ){
			min = positions[j];
			i =j;
		}
	}
	return min;
}

// Play with gravity here to make the movement appear more correct
void updateVelocity(vec3 currentPosition, float &v, vec3 maxHeight){
	float h = currentPosition.y;
	float H = maxHeight.y;
	
	v = sqrt( 2 * g * (H-h));
}

float circumference(float r){
	return 2 * PI * r;
}

float calculateCurvature(float t, float dt){
	return ((circumference(t+dt) - (2* circumference(t)) + circumference(t - dt))/ ( dt * dt));
}


vec3 calculateCurvature2(vec3 pi1, vec3 pi, vec3 pm1){
	// divide by multplication of the two differences 
	return (pi1 - (2.f * pi) + pm1) / ((pi1-(pi) * (pi-pm1)));
} 

vec3 getNormal(float v, float a){
	return ((v*v)*a)*vec3(1.f) + vec3(0,g,0);
	//return (a -g)*vec3(1.f);
}


vec3 getNormal2(float v, vec3 a){
	return ((v*v)*a) + vec3(0,g,0);
	//return (a -g)*vec3(1.f);
}

float getTangentCo(float t, float dt){
	return (circumference(t+dt) - circumference(t))/ dt;
}
vec3 getTangentCo2(vec3 pi1, vec3 pi){
	return pi1-pi;
}

vec3 calculateCentripedalAcceleration(vec3 pi1, vec3 pi, vec3 pm1){
	vec3 n = (pi1- (2.f * pi) + pm1) / magnitude((pi1- (2.f * pi) + pm1));
		
	float x = magnitude(pi1- (2.f * pi) + pm1) /2;
	float c = magnitude(pi1 - pm1) /2;
	
	float k = (2*x) / ((x*x)+(c*c));
	
	return k*n;
}

float getVelocity(vec3 currentPosition, vec3 maxHeight){
	float h = currentPosition.y;
	float H = maxHeight.y;
	
	return sqrt( 2 * g * (H-h));
}

float lacc(vec3 pi, vec3 pi1, vec3 maxHeight){
	float num = getVelocity(pi1, maxHeight) - getVelocity(pi, maxHeight);
	float den = magnitude(pi1 - pi);
	
	return num/den;
}


// ==========================================================================
// PROGRAM ENTRY POINT

int main(int argc, char *argv[])
{   
    window = createGLFWWindow();
    if(window == NULL)
    	return -1;

    //Initialize glad
    if (!gladLoadGL())
	{
		cout << "GLAD init failed" << endl;
		return -1;
	}

    // query and print out information about our OpenGL environment
    QueryGLVersion();

	initGL();

	//Initialize shader
	GLuint program = initShader("vertex.glsl", "fragment.glsl");
	GLuint vao, vaoSphere, vaoSquare, vaoBox;
	VertexBuffers vbo, vboSphere, vboSquare, vboBox;

	//Generate object ids
	glGenVertexArrays(1, &vao);
	glGenBuffers(VertexBuffers::COUNT, vbo.id);

	//Generate object ids
	glGenVertexArrays(1, &vaoSphere);
	glGenBuffers(VertexBuffers::COUNT, vboSphere.id);

	//Generate object ids for the "floor"
	glGenVertexArrays(1, &vaoSquare);
	glGenBuffers(VertexBuffers::COUNT, vboSquare.id);
	
	//Generate object ids for the cart
	glGenVertexArrays(1, &vaoBox);
	glGenBuffers(VertexBuffers::COUNT, vboBox.id);

	initVAO(vao, vbo);

	initVAO(vaoSphere, vboSphere);

	initVAO(vaoSquare, vboSquare);
	
	initVAO(vaoBox, vboBox);
	
	//Geometry information
	vector<vec3> points, normals;
	vector<unsigned int> indices;
	
	//Sphere specifically
	vector<vec3> pointsSphere, normalsSphere;
	vector<unsigned int> indicesSphere;

	//Square Specifically
	vector<vec3> pointsSquare, normalsSquare;
	vector<unsigned int> indicesSquare;
	
	//Cube Specifically
	vector<vec3> pointsBox, normalsBox;
	vector<unsigned int> indicesBox;
	
	//READ IN TRACK FROM FILES
	vector<vec3> controlPoints = loadInVecsFromFile();
	vector<vec3> temp;
	for (vec3 c : controlPoints){
		temp.push_back(c*10.f);
	}
	controlPoints = temp;
	
	int subdivisions = 6;

	points = subdivideTrack(controlPoints, &normals, &indices,  subdivisions);
	loadBuffer(vbo, points, normals, indices);
	
	generateCylinder(&pointsSphere, &normalsSphere, &indicesSphere, 0.01, vec3(0,0,0), 100,100, 0.2);
	loadBuffer(vboSphere, pointsSphere, normalsSphere, indicesSphere);

	generateSquare(pointsSquare, normalsSquare, indicesSquare,0.5f);
	loadBuffer(vboSquare, pointsSquare, normalsSquare, indicesSquare);
	
	generateBox(&pointsBox, &normalsBox, &indicesBox, 0.05);
	loadBuffer(vboBox, pointsBox, normalsBox, indicesBox);
	
	
	Camera cam = Camera(vec3(0, 0, -1), vec3(0, 0, 7));
	
	activeCamera = &cam;
	
	//float fovy, float aspect, float zNear, float zFar
	mat4 perspectiveMatrix = perspective(radians(80.f), 1.f, 0.1f, 20.f);
	
	int i =0;
	
	// Only calculate one of highPoint and lowPoint at a time
	vec3 highPoint = getMaxHeight(points, i);
	//vec3 lowPoint = getMinHeight(points , i);
	
	vec3 start = highPoint;
	
	mat4 M2 = mat4(1.f);
	vector<mat4> trackRotations;
	
	float v = 0.1;

	float tempT = 0;
	float tempDT =0.01;
	for ( int trackIndex = 0; trackIndex < (int)indices.size(); trackIndex++){
		tempT += tempDT;
		vec3 kn = calculateCentripedalAcceleration(points[(trackIndex+3)%indices.size()], points[trackIndex], points[((trackIndex-3)<0)? points.size()-3:(trackIndex-3)]);
		updateVelocity(points[trackIndex], v, highPoint);
		
		vec3 currentTangent = normalize(getTangentCo2(points[(trackIndex+1)%points.size()], points[((trackIndex-1)<0)?points.size()-1:(trackIndex-1)]));	
		vec3 norm = normalize((v * v) * kn + vec3(0,g,0) + lacc(points[trackIndex], points[(trackIndex+1)%points.size()], highPoint) * currentTangent);
	
		vec3 B = normalize(cross(currentTangent, norm));	
				
		vec3 finalT = currentTangent;//normalize(cross(norm, B));
		vec3 finalN = normalize(cross( B, currentTangent));
	
		mat4 face = mat4(
				vec4(B,0), 
				vec4(finalN,0),
				vec4(finalT,0),
				vec4(points[indices[trackIndex]],1));
		//face = scale(face, vec3(2.f));
		trackRotations.push_back(face);
	}

	
	v = 0.01;
	float t = 0;
	float dt = 0.00001;

	mat4 fM = mat4(1.f);
	
    // run an event-triggered main loop
    while (!glfwWindowShouldClose(window))
    {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		//Clear color and depth buffers (Haven't covered yet)
		
		loadUniforms(program, winRatio*perspectiveMatrix*cam.getMatrix(), M2 * scale(M2, vec3(50,1,50)));
		
		renderSquare(vaoSquare, 0, indicesSquare.size());
		loadUniforms(program, winRatio*perspectiveMatrix*cam.getMatrix(), M2);
		if(drawTrack){
			renderTrack(vao, 0, indices.size());
		}
		for( int i = 0; i < (int)indices.size(); i++){
			mat4 face = trackRotations[indices[i]];
			loadUniforms(program, winRatio*perspectiveMatrix*cam.getMatrix(), face);
			render(vaoSphere, 0 , indicesSphere.size());
		}
		
		if (animate){
			t += dt;
			vec3 a = calculateCentripedalAcceleration(points[(i+2)%points.size()], points[i], points[((i-2)<0)?points.size()-2:(i-2)]);
			updateVelocity(start, v, highPoint);
			start = arcLengthParamaterize(start, i, points, v+0.01);
			vec3 currentTangent = normalize(getTangentCo2(points[(i+1)%points.size()], points[i]));
			vec3 norm = normalize((v * v) * a + vec3(0,g,0) + lacc(points[i], points[(i+1)%points.size()], highPoint) * currentTangent);
			
			vec3 B = normalize(cross(currentTangent, norm));
			
			vec3 finalN = normalize(cross( B, currentTangent));
			fM = mat4(
					vec4(B,0), 
					vec4(finalN,0),
					vec4(currentTangent,0),
					vec4(start,1));
			if ( firstPerson){
				cam.pos = start + vec3(0,0.09,0);
				cam.dir = currentTangent;
				cam.up = finalN;
				cam.right = B;
			} 
		}
		
		loadUniforms(program, winRatio*perspectiveMatrix*cam.getMatrix(), fM);
		renderBox(vaoBox, 0, indicesBox.size());
		
        // scene is rendered to the back buffer, so swap to front for display
        glfwSwapBuffers(window);

        // sleep until next event before drawing again
        glfwPollEvents();
	}

	// clean up allocated resources before exit
	glDeleteVertexArrays(1, &vao);
	glDeleteBuffers(VertexBuffers::COUNT, vbo.id);
	glDeleteProgram(program);


	glfwDestroyWindow(window);
   glfwTerminate();

   return 0;
}

// ==========================================================================
// SUPPORT FUNCTION DEFINITIONS

// --------------------------------------------------------------------------
// OpenGL utility functions

void QueryGLVersion()
{
    // query opengl version and renderer information
    string version  = reinterpret_cast<const char *>(glGetString(GL_VERSION));
    string glslver  = reinterpret_cast<const char *>(glGetString(GL_SHADING_LANGUAGE_VERSION));
    string renderer = reinterpret_cast<const char *>(glGetString(GL_RENDERER));

    cout << "OpenGL [ " << version << " ] "
         << "with GLSL [ " << glslver << " ] "
         << "on renderer [ " << renderer << " ]" << endl;
}

bool CheckGLErrors(string location)
{
    bool error = false;
    for (GLenum flag = glGetError(); flag != GL_NO_ERROR; flag = glGetError())
    {
        cout << "OpenGL ERROR:  ";
        switch (flag) {
        case GL_INVALID_ENUM:
            cout << location << ": " << "GL_INVALID_ENUM" << endl; break;
        case GL_INVALID_VALUE:
            cout << location << ": " << "GL_INVALID_VALUE" << endl; break;
        case GL_INVALID_OPERATION:
            cout << location << ": " << "GL_INVALID_OPERATION" << endl; break;
        case GL_INVALID_FRAMEBUFFER_OPERATION:
            cout << location << ": " << "GL_INVALID_FRAMEBUFFER_OPERATION" << endl; break;
        case GL_OUT_OF_MEMORY:
            cout << location << ": " << "GL_OUT_OF_MEMORY" << endl; break;
        default:
            cout << "[unknown error code]" << endl;
        }
        error = true;
    }
    return error;
}

// --------------------------------------------------------------------------
// OpenGL shader support functions

// reads a text file with the given name into a string
string LoadSource(const string &filename)
{
    string source;

    ifstream input(filename.c_str());
    if (input) {
        copy(istreambuf_iterator<char>(input),
             istreambuf_iterator<char>(),
             back_inserter(source));
        input.close();
    }
    else {
        cout << "ERROR: Could not load shader source from file "
             << filename << endl;
    }

    return source;
}

// creates and returns a shader object compiled from the given source
GLuint CompileShader(GLenum shaderType, const string &source)
{
    // allocate shader object name
    GLuint shaderObject = glCreateShader(shaderType);

    // try compiling the source as a shader of the given type
    const GLchar *source_ptr = source.c_str();
    glShaderSource(shaderObject, 1, &source_ptr, 0);
    glCompileShader(shaderObject);

    // retrieve compile status
    GLint status;
    glGetShaderiv(shaderObject, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE)
    {
        GLint length;
        glGetShaderiv(shaderObject, GL_INFO_LOG_LENGTH, &length);
        string info(length, ' ');
        glGetShaderInfoLog(shaderObject, info.length(), &length, &info[0]);
        cout << "ERROR compiling shader:" << endl << endl;
        cout << source << endl;
        cout << info << endl;
    }

    return shaderObject;
}

// creates and returns a program object linked from vertex and fragment shaders
GLuint LinkProgram(GLuint vertexShader, GLuint fragmentShader)
{
    // allocate program object name
    GLuint programObject = glCreateProgram();

    // attach provided shader objects to this program
    if (vertexShader)   glAttachShader(programObject, vertexShader);
    if (fragmentShader) glAttachShader(programObject, fragmentShader);

    // try linking the program with given attachments
    glLinkProgram(programObject);

    // retrieve link status
    GLint status;
    glGetProgramiv(programObject, GL_LINK_STATUS, &status);
    if (status == GL_FALSE)
    {
        GLint length;
        glGetProgramiv(programObject, GL_INFO_LOG_LENGTH, &length);
        string info(length, ' ');
        glGetProgramInfoLog(programObject, info.length(), &length, &info[0]);
        cout << "ERROR linking shader program:" << endl;
        cout << info << endl;
    }

    return programObject;
}

vector<vec3> loadInVecsFromFile(){
	std::string file( "./track3.con" );
	//std::string file( "./infin.con" );
	VectorContainerVec3f vecs;
	loadVec3fFromFile( vecs, file );
	
	vector<vec3> controlPoints;

	for( Vec3f &v: vecs )
	{
		controlPoints.push_back(vec3(v.m_x,v.m_y,v.m_z));
	}
	
	return controlPoints;
}

// ==========================================================================
