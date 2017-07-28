#include <GL\glew.h>
#include <GL\glut.h>
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

struct vec2
{
	float x;
	float y;
};
struct vec4
{
	float r;
	float g;
	float b;
	float a;
};

class Particle
{
public:
	Particle()
	{
		Particle(vec2());
	}
	vec2 pos;
	vec2 vel;
	vec2 force;
	float mass;
	float den;
	float pres;

	vec2 pos_p;
	vec2 vel_p;
	float den_p;

	Particle(vec2 position)
	{
		pos = position;
		vel = vec2();
		force = vec2();
		mass = 0;
		den = 0;
		pres = 0;

		pos_p = vec2();
		vel_p = vec2();
		den_p = 0;
	}
};
