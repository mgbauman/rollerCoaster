// ==========================================================================
// Vertex program for barebones GLFW boilerplate
//
// Author:  Sonny Chan, University of Calgary
// Date:    December 2015
// ==========================================================================
#version 410

// interpolated colour received from vertex stage

// first output is mapped to the framebuffer's colour index by default
out vec4 FragmentColour;

in vec3 FragNormal;

void main(void)
{
	FragmentColour = vec4(FragNormal, 1.0);
	//FragmentColour = vec4(0,0,1,1);
}
