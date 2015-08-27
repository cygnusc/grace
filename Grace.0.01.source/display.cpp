/*

display.cpp - Grace Micromagnetic Simulator written by Ru Zhu, zhu@sting.graceland.edu

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "display.h"
#include "Configuration.h"
//#include "utils.h"
//#include "magnetic.h"

//#define INDEX(i,j,k) (k * myConfig->disDimX * myConfig->disDimY + j * myConfig->disDimX + i)
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#endif

#include <cstdlib>
#include <cmath>

#define COLORVEC 0

extern Configuration * myConfig;

float shiftZ = 2.f;
float eyex = 0., eyey = 0., eyez = shiftZ; // camera position
float centerx = 0., centery = 0., centerz = 0.; // actual vector representing the camera's direction
float deltaAnglex = 0.0f;
float deltaAngley = 0.0f;
int xOrigin = -1;
int yOrigin = -1;
// angle of rotation for the camera direction
float anglex = 0.0f;
float angley = 0.0f;
static bool paused = false;

void GLInit(int argc, char **argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(640, 640);
	glutCreateWindow("Magnetic2");
}

void idleDisplay() {
	
	if (myConfig->frame == 0) paused = true;
	if (!paused) {
		myConfig->readNext();
		static char buffer[10];
		sprintf_s(buffer, "%3.0f %%", 100.f - float (myConfig->frame) / myConfig->totalFrames * 100.f); 
		glutSetWindowTitle(buffer);
	
		glutPostRedisplay();
	} else { 
		Sleep (10);
	}
	
}

void reShape(int w, int h) {

	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if (h == 0)
		h = 1;
	float ratio =  w * 1.0f / h;

	// Use the Projection Matrix
	glMatrixMode(GL_PROJECTION);

	// Reset Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set the correct perspective.
	gluPerspective(45.0f, ratio, 0.1f, 100.0f);

	// Get Back to the Modelview
	glMatrixMode(GL_MODELVIEW);
}

void keyboardFunc(unsigned char key, int x, int y) {
	if (key == 27 || key == 113) // ESC or 'q'
		exit(0);
	if (key == 114) { // ASCII 'r'
		eyex = eyey = 0.;
		eyez = shiftZ; // reset camera position
		centerx = centery = centerz = 0.; 
		glutPostRedisplay();
	}
	if (key == 'p' || key == ' ') { // p or SPACE to pause
		paused = !paused;
	}
}

void displayVec(float X, float Y, float Z, float coordx, float coordy, float coordz) {
#if COLORVEC
	glColor3f(X, Y, Z);
#else
	glColor3d(1,0,0);
#endif
	float mag = sqrt (X * X + Y * Y + Z * Z);
	X /= mag;
	Y /= mag;
	Z /= mag;
	
	float PHI = 0.f, THETA = 0.f;
	if (X == 0.f) {
		if (Y >= 0.f) {
			PHI = 90.f;
		} else {
			PHI = -90.f;
		}
	} else {
			PHI = atan (Y / X) * 180.f / PI;
			if (X < 0.f) {
				PHI += 180.f;
			}
	}
	THETA = acos (Z / sqrt(X * X + Y * Y + Z * Z)) * 180.f / PI;
	glPushMatrix();
		//glTranslated(coordx, coordy, 3 * shiftZ + coordz);
	glTranslated(coordx, coordy, coordz);
		glRotated(PHI, 0., 0., 1.);
		//glRotated(180.f - THETA, 0., 1., 0.);
		glRotated(THETA, 0., 1., 0.);
		glutSolidCone(.01, .05,3,3);
    glPopMatrix(); 
}

/*void renderBitmapString(float x, float y, void *font,const char *string){     
	const char *c;     
	glRasterPos2f(x, y);     
	for (c=string; *c != '\0'; c++) { 
		glutBitmapCharacter(font, *c);    
	}
}
*/
/* code taken from 
http://programming-technique.blogspot.com/2012/05/font-rendering-in-glut-using-bitmap.html*/

/*void setOrthographicProjection() { 
	glMatrixMode(GL_PROJECTION); 
	glPushMatrix();   
	glLoadIdentity();     
	gluOrtho2D(0, w, 0, h);    
	glScalef(1, -1, 1);    
	glTranslatef(0, -h, 0);  
	glMatrixMode(GL_MODELVIEW); 
}  */



void displayFunc() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glColor3d(1,0,0); 

	
	glLoadIdentity();

	/*if (angley + deltaAngley < -PI || angley + deltaAngley > PI)
		gluLookAt(eyex, eyey, eyez, centerx, centery, centerz, 0., 1., 0.);
	else*/
	gluLookAt(eyex, eyey, eyez, centerx, centery, centerz, 0., 1., 0.);
#if DISPLAY_TEST
	displayVec(1., 0., 0., 0., 0., 0.); // display test
#else

	for (int i = 0; i < myConfig->disDimX; i+=1) {
	for (int j = 0; j < myConfig->disDimY; j+=1) {
	for (int k = 0; k < myConfig->disDimZ; k+=1) {
		displayVec(
				   myConfig->vecX[myConfig->INDEX(i,j,k)], myConfig->vecY[myConfig->INDEX(i,j,k)], myConfig->vecZ[myConfig->INDEX(i,j,k)],
				   myConfig->X[myConfig->INDEX(i,j,k)] / (myConfig->maxDimX) - 0.5f, 
				   myConfig->Y[myConfig->INDEX(i,j,k)] / (myConfig->maxDimX) - 0.5f * myConfig->maxDimY / myConfig->maxDimX,
				   myConfig->Z[myConfig->INDEX(i,j,k)] / (myConfig->maxDimX) - 0.5f * myConfig->maxDimZ / myConfig->maxDimX
                   /*
				   myConfig->X[myConfig->INDEX(i,j,k)] / (myConfig->maxDimX) - 0.5, 
				   myConfig->Y[myConfig->INDEX(i,j,k)] / (myConfig->maxDimY) - 0.5, 
				   myConfig->Z[myConfig->INDEX(i,j,k)] / (myConfig->maxDimZ) - 0.5
                   */
		);
	}
	}
	}
	
#endif
	//Sleep(10);
	//Sleep (100);
    glutSwapBuffers();
}

void mouseButton(int button, int state, int x, int y) {

	// only start motion if the left button is pressed
	if (button == GLUT_LEFT_BUTTON) {

		// when the button is released
		if (state == GLUT_UP) {
			anglex += deltaAnglex;
			angley += deltaAngley;
			xOrigin = -1;
			yOrigin = -1;
		}
		else  {// state = GLUT_DOWN
			xOrigin = x;
			yOrigin = y;
		}
	}
}

void mouseMove(int x, int y) {

	// this will only be true when the left button is down
	if (xOrigin >= 0) {

		// update deltaAngle
		deltaAnglex = (x - xOrigin) * 0.01f;
		deltaAngley = (y - yOrigin) * 0.01f;

		// update camera's direction
		/*centerx = sin(anglex + deltaAnglex);
		centery = sin(angley + deltaAngley);
		centerz = shiftZ + 1. - cos(anglex + deltaAnglex);*/
		eyex = shiftZ * cos(angley + deltaAngley) * -sin(anglex + deltaAnglex);
		eyey = shiftZ * sin(angley + deltaAngley);
		eyez = shiftZ * cos(angley + deltaAngley) * cos(anglex + deltaAnglex);
	}
}



void Display::Run(int argc, char** argv) {
	GLInit(argc, argv); // glut init. functions
	glutKeyboardFunc(keyboardFunc); // press esc or q to quit, r to redraw display
//	if (!paused) {
	glutDisplayFunc(displayFunc); // set vector positions and magnitudes
//	glutTimerFunc(10, idle, 0);
//	}
	glutIdleFunc(idleDisplay); // update system parameters
	glutReshapeFunc(reShape); // handles display resizing
	glutMouseFunc(mouseButton); // handles mouse button clicking
	glutMotionFunc(mouseMove); // handles mouse moving

	// OpenGL init
	glEnable(GL_DEPTH_TEST);

	// enter GLUT event processing cycle
	glutMainLoop();
}

