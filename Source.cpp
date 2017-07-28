/*Few computations for sph and opengl apis referred from online sources*/
/*Image Loader is taken from Opengl tutorials*/
#include <iostream>
#include <conio.h>
#include "sph.h"
#include "imageloader.h"
using namespace std;
GLuint _textureId;
const float step1 = 0.0002f;
const float step2 = 0.00015f;
float step_;
int FLAG = -100;
sph s1 = sph();
int algo = 0;
static bool paused = false;
int fact = 0;
void drawBox()
{
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glBegin(GL_LINES);
	glLineWidth(5.0);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex2f(-3, -3);
	glVertex2f(203, -3);
	glVertex2f(203, -3);
	glVertex2f(203, 150);
	glVertex2f(-3, -3);
	glVertex2f(-3, 150);
	glEnd();
	glDisable(GL_BLEND);
}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	drawBox();

	if(algo == 0)
		s1.updateSPHSystem(step_);
	else if (algo == 1)
		s1.updatePCISPHSystem(step_);

	s1.renderParticles(FLAG);
	glutSwapBuffers();
}

void initWithoutTexture(float width, float height)
{
	FLAG = -100;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0.8f, 0.55f, 0.0f, 0.0f);
	glEnable(GL_DEPTH_TEST);
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(90.0, (float)width / height, 0.001, 1000.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(-100.0f, -100.0f, -115.0f);
}

void initWithTexture(float width, float height)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0.8f, 0.55f, 0.0f, 0.0f);
	glEnable(GL_DEPTH_TEST);
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(90.0, (float)width / height, 0.001, 1000.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(-100.0f, -100.0f, -115.0f);
	glEnable(GL_COLOR_MATERIAL);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	Image* image = loadBMP("water.bmp");
	Image* alphaImage = loadBMP("circlealpha.bmp");
	_textureId = loadAlphaTexture(image, alphaImage);
	FLAG = _textureId;
	delete image;
	delete alphaImage;
}

void keyboard_func(unsigned char key, int x, int y)
{
	/*SPH*/
	if (key == ' ')
	{
		s1 = sph(25, 3.5, step1);
		step_ = step1;
		algo = 0;
		initWithTexture(600, 600);
	}
	if (key == '1')
	{
		s1 = sph(20, 5, step1);
		step_ = step1;
		algo = 0;
	}
	if (key == '2')
	{
		s1 = sph(25, 3.5, step1);
		step_ = step1;
		algo = 0;
	}
	if (key == '3')
	{
		s1 = sph(30, 3, step1);
		step_ = step1;
		algo = 0;
	}
	if (key == '4')
	{
		s1 = sph(40, 2.0, step2);
		step_ = step1;
		algo = 0;
	}
	/*if (key == '5')
	{
		s1 = sph(46, 1, step2);
		step_ = step2;
		algo = 0;
	}*/

	/*PCISPH*/
	if (key == 'p' || key == 'P')
	{
		s1 = sph(25, 3.5, step1);
		step_ = step1;
		algo = 1;
		initWithTexture(600, 600);
	}
	/*if (key == '0')
	{
		s1 = sph(10, 8, step1);
		step_ = step1;
		algo = 1;
	}*/
	if (key == '9')
	{
		s1 = sph(20, 3.5, step1);
		step_ = step1;
		algo = 1;
	}
	if (key == '8')
	{
		s1 = sph(30, 3.3, step1);
		step_ = step1;
		algo = 1;
	}
	/*if (key == '7')
	{
		s1 = sph(40, 2.0, step1);
		step_ = step1;
		algo = 1;
	}*/

	if (key == 'z' || key == 'Z')
	{
		initWithTexture(600, 600);
	}

	if (key == 'x' || key == 'X')
	{
		initWithoutTexture(600, 600);
	}

	if (key == 'q' || key == 'Q')
	{
		exit(0);
	}

	glutPostRedisplay();
}

void idle_func()
{
	glutPostRedisplay();
}

int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(600, 600);
	glutCreateWindow("Demo");
	s1 = sph(25, 3.5, step1);
	step_ = step1;
	algo = 0;
	initWithTexture(600, 600);
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard_func);
	glutIdleFunc(idle_func);
	glutMainLoop();
	return 0;
}