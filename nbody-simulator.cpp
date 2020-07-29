#include <stdio.h>
#include <windows.h>
#include <GL/glut.h>
#include <ctime>
#include "barnes_hut.h"

double dims = 1000.0;
int bodyCnt = 10;

void init(int N_, double G){
	bodyCnt = N_;
	setSimulationParams(bodyCnt, G, 2.0, 1000, dims);
	initBodies(); //initialize all bodies for grav simulation	
	glClearColor(0.0, 0.0, 0.0, 1.0);
}

void display(){
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();
	
	glPointSize(2.0);
	
	double *points = getBodyPoints();
	
	glBegin(GL_POINTS);
	
		for( int i = 0; i < bodyCnt; i++ ){
			glVertex2f(points[i*2], points[i*2+1]);
		}
	
	glEnd();
	evolveBodiesStep(1.0);
	glutPostRedisplay();
	glFlush();
}

void reshape(int wid, int hig){
	glViewport(0, 0, wid, hig);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//left, right, bottom ,top
	gluOrtho2D(-dims, dims, -dims, dims);
	glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB);
	
	double G;
	int N_;
	printf("Enter number of bodies : ");
	scanf("%d", &N_);
	printf("Enter G value as a double: ");
	scanf("%lf", &G);
	
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(700, 700);
	
	glutCreateWindow("N Body Simulator");
	
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	init(N_, G);
	
	glutMainLoop();
}