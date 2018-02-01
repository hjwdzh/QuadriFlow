#include "gldraw.hpp"

#include <stdio.h> // printf
#include <stdlib.h> // exit
#include <string.h> // strlen

#if defined __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

static int window_width = 800;
static int window_height = 800;
static int cursor_x = window_width / 2;
static int cursor_y = window_height / 2;

static GLint glyphs_display_list = 0;

static void glutInitGlyphs()
{
	//just doing 7-bit ascii
	glyphs_display_list = glGenLists(256);

	for (int i = 0; i < 256; i++) {
		glNewList(glyphs_display_list + i, GL_COMPILE);
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, i);
		glEndList();
	}
}

static void glutDrawString(int x, int y, const char * text)
{
	glListBase(glyphs_display_list);
	glRasterPos2i(x + 1, y + 4);
	glCallLists((GLsizei)strlen(text), GL_UNSIGNED_BYTE, text);
}

static void glutKeyboardUpCallback(unsigned char key, int x, int y)
{
	//printf( "keyup=%i\n", key );
}


static void glutDisplayCallback(void)
{
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glBegin(GL_LINES);
	glVertex2f(cursor_x, cursor_y);
	glVertex2f(window_width / 2, window_height / 2);
	glEnd();

	glutDrawString(16, 16, "Hello GLUT");

	glutSwapBuffers();
}

static void glutReshapeCallback(int width, int height)
{
	window_width = width;
	window_height = height;

	glViewport(0, 0, window_width, window_height);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, /*near=*/-1, /*far=*/1);

	glutPostRedisplay();
}

/*static void glutIdleCallback(void)
{
glutPostRedisplay();
}*/

/*static void glutTimerFunc(int value)
{
glutPostRedisplay();
}*/



void gldraw(void(*mouse_callback)(int, int, int, int),
	void(*render_callback)(void),
	void(*motion_callback)(int, int),
	void(*keyboard_callback)(unsigned char, int, int))
{
	int argc = 1;
	char dummy[6] = "dummy";
	char* argv[1];
	argv[0] = dummy;
	glutInit(&argc, argv);

	//glutInitDisplayString("rgb depth double samples=4");
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(window_width, window_height);
	glutCreateWindow("glut-template");
	GLfloat light_ambient[] =
	{ 0.2, 0.2, 0.2, 1.0 };
	GLfloat light_diffuse[] =
	{ 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular[] =
	{ 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] =
	{ 1.0, 1.0, 1.0, 0.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, light_position);

	glEnable(GL_LIGHT0);
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);

	glutKeyboardFunc(keyboard_callback);
	glutKeyboardUpFunc(glutKeyboardUpCallback);
	glutMouseFunc(mouse_callback);
	glutMotionFunc(motion_callback);

	//glutPassiveMotionFunc(glutMotionCallback);
	glutDisplayFunc(render_callback);
	glutReshapeFunc(glutReshapeCallback);
	//glutIdleFunc( glutIdleCallback );
	//glutTimerFunc( 100, glutTimerFunc, 0 );

	glutInitGlyphs();

	glutPostRedisplay();

	glutMainLoop();

	return;
}
