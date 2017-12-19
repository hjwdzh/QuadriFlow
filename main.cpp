#include "Parametrizer.h"
#include "Optimizer.h"
#include <GL/glut.h>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include "glDraw.h"

glm::dvec3 render_translation;
glm::dmat4 render_rotation = glm::dmat4(1.0);
double render_scale = 0.0;
int render_wireframe = 0;
int mouse_state = 0;
int mouse_x = 0, mouse_y = 0;
int level = 0;
int show_field = 0;
int show_mesh = 1;
int show_quad = 0;
int show_hierarchy = 0;
Parametrizer field;

static void render_callback(void)
{
	glClearColor(0.0, 191.0 / 255.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	GLfloat light_position[] = { 1, 1, 1, 0 };
	glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, light_position);
	// now you can setup view matrix (gluLookAt())

	glPushMatrix();
	glTranslated(render_translation.x, render_translation.y, render_translation.z);
	glMultMatrixd((double*)&render_rotation);
	double model_scale = exp(render_scale);
	glScaled(model_scale, model_scale, model_scale);
	if (render_wireframe) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

	auto& mF = field.hierarchy.mF;
	auto& mV = field.hierarchy.mV[0];
	auto& mN = field.hierarchy.mN[0];
	auto& mVq = field.mV_extracted;
	if (show_mesh) {
		static GLfloat white[4] =
		{ 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
		glBegin(GL_TRIANGLES);
		for (int i = 0; i < mF.cols(); ++i) {
			for (int j = 0; j < 3; ++j) {
				glNormal3fv(&mN(0, mF(j, i)));
				glVertex3fv(&mV(0, mF(j, i)));
			}
		}
		glEnd();

	}
	glDisable(GL_LIGHTING);
	glPointSize(5.0f);
	if (show_quad) {
		
		glColor3f(0.0f, 1.0f, 0.0f);
		glBegin(GL_POINTS);
		for (auto& p : field.vertex_singularities) {
			if (p.second == 1) {
				Vector3f v = (mVq.col(p.first));
				glVertex3f(v.x(), v.y(), v.z());
			}
		}
		glEnd();
		glColor3f(0.0f, 0.0f, 1.0f);
		glBegin(GL_POINTS);
		for (auto& p : field.vertex_singularities) {
			if (p.second == 3) {
				Vector3f v = (mVq.col(p.first));
				glVertex3f(v.x(), v.y(), v.z());
			}
		}
		glEnd();
	}
	else {
		glColor3f(0.0f, 1.0f, 0.0f);
		glBegin(GL_POINTS);
		for (auto& p : field.singularities) {
			if (p.second == 1) {
				Vector3f v = (mV.col(mF(0, p.first))
					+ mV.col(mF(1, p.first))
					+ mV.col(mF(2, p.first))) / 3.0f;
				glVertex3f(v.x(), v.y(), v.z());
			}
		}
		glEnd();
		glColor3f(0.0f, 0.0f, 1.0f);
		glBegin(GL_POINTS);
		for (auto& p : field.singularities) {
			if (p.second == 3) {
				Vector3f v = (mV.col(mF(0, p.first))
					+ mV.col(mF(1, p.first))
					+ mV.col(mF(2, p.first))) / 3.0f;
				glVertex3f(v.x(), v.y(), v.z());
			}
		}
		glEnd();
	}
	if (show_hierarchy) {
		glColor3f(0, 0, 1);
		glBegin(GL_LINES);
		auto& mV = field.hierarchy.mV[level];
		auto& mN = field.hierarchy.mN[level];
		auto& mQ = field.hierarchy.mQ[level];
		auto& adj = field.hierarchy.mAdj[level];
		for (int i = 0; i < adj.size(); ++i) {
			for (auto& l : adj[i]) {
				int j = l.id;
				glVertex3fv(&mV(0, i));
				glVertex3fv(&mV(0, j));
			}
		}
		glEnd();
	}
	
	if (show_quad) {
		glColor3f(1, 0, 0);
		glBegin(GL_LINES);
		for (auto& e : field.edge_idmap) {
			glVertex3fv(&field.mV_extracted(0, e.first.first));
			glVertex3fv(&field.mV_extracted(0, e.first.second));
		}
		/*
		for (int i = 0; i < field.adj_extracted.size(); ++i) {
			for (auto& l : field.adj_extracted[i]) {
				int j = l.id;
				glVertex3fv(&field.mV_extracted(0, i));
				glVertex3fv(&field.mV_extracted(0, j));
			}
		}
		*/
		glEnd();	
	}

	if (show_field) {
		int l = (show_hierarchy) ? level : 0;
		auto& mV = field.hierarchy.mV[level];
		auto& mN = field.hierarchy.mN[level];
		auto& mQ = field.hierarchy.mQ[level];
		auto& adj = field.hierarchy.mAdj[level];
		glColor3f(1, 0, 0);
		float len = field.scale * 0.2;
		glBegin(GL_LINES);
		for (int i = 0; i < mQ.cols(); ++i) {
			glm::vec3 p(mV(0, i), mV(1, i), mV(2, i));
			glm::vec3 n(mN(0, i), mN(1, i), mN(2, i));
			glm::vec3 tangent1 = glm::vec3(mQ(0, i), mQ(1, i), mQ(2, i)) * len;
			glm::vec3 tangent2 = glm::normalize(glm::cross(n, tangent1)) * len;
			auto l = p - tangent1;
			auto r = p + tangent1;
			auto u = p - tangent2;
			auto d = p + tangent2;
			glVertex3f(l.x, l.y, l.z);
			glVertex3f(r.x, r.y, r.z);
			glVertex3f(u.x, u.y, u.z);
			glVertex3f(d.x, d.y, d.z);
		}
		glEnd();
	}
	glEnable(GL_LIGHTING);
	glPopMatrix();
	glutSwapBuffers();
}

static void motion_callback(int x, int y)
{
	if (mouse_state == 1) {
		render_rotation = glm::rotate((double)(x - mouse_x) / 200, glm::dvec3(0, 1, 0)) * render_rotation;
		render_rotation = glm::rotate((double)(y - mouse_y) / 200, glm::dvec3(1, 0, 0)) * render_rotation;
		mouse_x = x;
		mouse_y = y;
	}
	if (mouse_state == 2) {
		glm::dmat4 rot = glm::inverse(render_rotation);
		render_translation += (double)(x - mouse_x) / 400 * glm::dvec3(1, 0, 0);
		render_translation -= (double)(y - mouse_y) / 400 * glm::dvec3(0, 1, 0);
		mouse_x = x;
		mouse_y = y;
	}
	else {
		render_scale += (double)(y - mouse_y) / 400;
		mouse_x = x;
		mouse_y = y;
	}
	glutPostRedisplay();
}

static void mouse_callback(int button, int state, int x, int y)
{
	if (state == 1) {
		mouse_state = 0;
		return;
	}
	if (mouse_state == 0) {
		mouse_x = x;
		mouse_y = y;
	}

	if (button == GLUT_LEFT_BUTTON) {
		mouse_state = 1;
	}
	else if (button == GLUT_RIGHT_BUTTON) {
		mouse_state = 2;
	}
	else if (button == GLUT_MIDDLE_BUTTON) {
		mouse_state = 3;
	}

	glutPostRedisplay();
}

static void keyboard_callback(unsigned char key, int x, int y)
{
	int modifiers = glutGetModifiers();

	if (key == 'z') {
		render_wireframe = 1 - render_wireframe;
		glutPostRedisplay();
		return;
	}

	if (key == 'u') {
		level += 1;
		if (level >= field.hierarchy.mV.size())
			level = field.hierarchy.mV.size() - 1;
		printf("vertices = %d\n", field.hierarchy.mV[level].size());
		glutPostRedisplay();
	}
	if (key == 'd') {
		level -= 1;
		if (level < 0)
			level = 0;
		glutPostRedisplay();
	}

	if (key == 'f') {
		show_field = 1 - show_field;
		glutPostRedisplay();
	}

	if (key == 'm') {
		show_mesh = 1 - show_mesh;
		glutPostRedisplay();
	}

	if (key == 'q') {
		show_quad = 1 - show_quad;
		glutPostRedisplay();
	}

	if (key == 'h') {
		show_hierarchy = 1 - show_hierarchy;
		glutPostRedisplay();
	}

	if (modifiers & GLUT_ACTIVE_SHIFT) {

	}
	if (modifiers & GLUT_ACTIVE_CTRL) {

	}
	if (modifiers & GLUT_ACTIVE_ALT) {

	}

	if (key == 27) {
		exit(0);
	}
}

int main(int argc, char** argv)
{
	field.Load(argv[1]);
	field.Initialize();
	Optimizer::optimize_orientations(field.hierarchy);
	field.ComputeOrientationSingularities();

	Optimizer::optimize_positions(field.hierarchy);

	field.ExtractMesh();
	gldraw(mouse_callback, render_callback, motion_callback, keyboard_callback);
	return 0;
}