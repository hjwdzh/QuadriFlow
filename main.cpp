#include "Parametrizer.h"
#include "Optimizer.h"
#include <GL/glut.h>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include "glDraw.h"
#include "field_math.h"
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
int show_loop = 0;
int show_singularity = 0;
int show_color = 0;
std::vector<Vector3d> color;

Parametrizer field;

Vector3d Gray2HSV(double gray)
{
	Vector3d res;
	double* r = &res[0];
	double* g = &res[1];
	double* b = &res[2];
	int i;
	int h = gray * 360.0;
	unsigned char s = 100, v = 100;
	double RGB_min, RGB_max;
	RGB_max = v*2.55f;
	RGB_min = RGB_max*(100 - s) / 100.0f;

	i = h / 60;
	int difs  = h  % 60; // factorial part of h  

	// RGB adjustment amount by hue   
	double RGB_Adj  = (RGB_max  - RGB_min)*difs  / 60.0f;

	switch (i) {
	case 0:
		*r  = RGB_max;
		*g  = RGB_min  + RGB_Adj;
		*b  = RGB_min;
		break;
	case 1:
		*r  = RGB_max  - RGB_Adj;
		*g  = RGB_max;
		*b  = RGB_min;
		break;
	case 2:
		*r  = RGB_min;
		*g  = RGB_max;
		*b  = RGB_min  + RGB_Adj;
		break;
	case 3:
		*r  = RGB_min;
		*g  = RGB_max  - RGB_Adj;
		*b  = RGB_max;
		break;
	case 4:
		*r  = RGB_min  + RGB_Adj;
		*g  = RGB_min;
		*b  = RGB_max;
		break;
    default:        // case 5:  
		*r  = RGB_max;
		*g  = RGB_min;
		*b  = RGB_max  - RGB_Adj;
		break;
	}
	return res / 255.0f;
}

void render_test_travel(int f)
{
	auto& mF = field.hierarchy.mF;
	auto& mV = field.hierarchy.mV[0];
	auto& mN = field.hierarchy.mN[0];
	auto& mVq = field.mV_extracted;

	Vector3d p = mV.col(mF(0, f)) + mV.col(mF(1, f)) + mV.col(mF(2, f));
	p *= 1.0f / 3;
	double len = field.hierarchy.mScale * 10;
	int f1 = f;
	Vector3d q = Travel(p, field.hierarchy.mQ[0].col(mF(0, f)), len, f1, field.hierarchy.mE2E, mV, mF, field.Nf, field.triangle_space);

	glPointSize(10.0f);
	glBegin(GL_POINTS);
	glColor3f(1, 0, 0);
	glVertex3f(p.x(), p.y(), p.z());
	glColor3f(0, 1, 0);
	glVertex3f(q.x(), q.y(), q.z());
	glEnd();
}

static void render_mesh()
{
	if (show_color == 0) {
		auto& mF = field.hierarchy.mF;
		auto& mV = field.hierarchy.mV[0];
		auto& mN = field.hierarchy.mN[0];
		auto& mVq = field.mV_extracted;

		glEnable(GL_LIGHTING);
		if (show_mesh) {
			static GLfloat white[4] =
			{ 1.0, 1.0, 1.0, 1.0 };
			glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
			glBegin(GL_TRIANGLES);
			for (int i = 0; i < mF.cols(); ++i) {
				for (int j = 0; j < 3; ++j) {
					glNormal3dv(&mN(0, mF(j, i)));
					glVertex3dv(&mV(0, mF(j, i)));
				}
			}
			glEnd();
		}
		glDisable(GL_LIGHTING);
	}
	else {
		glDisable(GL_LIGHTING);
		auto& mF = field.hierarchy.mF;
		auto& mV = field.hierarchy.mV[0];
		auto& mN = field.hierarchy.mN[0];
		auto& mVq = field.mV_extracted;
		if (show_mesh) {
			glBegin(GL_TRIANGLES);
			for (int i = 0; i < mF.cols(); ++i) {
				for (int j = 0; j < 3; ++j) {
					glColor3d(color[mF(j, i)].x(), color[mF(j, i)].y(), color[mF(j, i)].z());
					glNormal3dv(&mN(0, mF(j, i)));
					glVertex3dv(&mV(0, mF(j, i)));
				}
			}
			glEnd();
		}
	}
}

static void render_singularities()
{
	if (show_singularity == 0)
		return;
	auto& mF = field.hierarchy.mF;
	auto& mV = field.hierarchy.mV[0];
	auto& mN = field.hierarchy.mN[0];
	auto& mVq = field.mV_extracted;
	glPointSize(5.0f);
	if (show_quad) {

		glColor3f(0.0f, 1.0f, 0.0f);
		glBegin(GL_POINTS);
		for (auto& p : field.vertex_singularities) {
			if (p.second == 1) {
				Vector3d v = (mVq.col(p.first));
				glVertex3d(v.x(), v.y(), v.z());
			}
		}
		glEnd();
		glColor3f(0.0f, 0.0f, 1.0f);
		glBegin(GL_POINTS);
		for (auto& p : field.vertex_singularities) {
			if (p.second == 3) {
				Vector3d v = (mVq.col(p.first));
				glVertex3d(v.x(), v.y(), v.z());
			}
		}
		glEnd();
	}
	else {
		glColor3f(0.0f, 1.0f, 0.0f);
		glBegin(GL_POINTS);
		for (auto& p : field.singularities) {
			if (p.second == 1) {
				Vector3d v = (mV.col(mF(0, p.first))
					+ mV.col(mF(1, p.first))
					+ mV.col(mF(2, p.first))) / 3.0f;
				glVertex3d(v.x(), v.y(), v.z());
			}
		}
		glEnd();
		glColor3f(0.0f, 0.0f, 1.0f);
		glBegin(GL_POINTS);
		for (auto& p : field.singularities) {
			if (p.second == 3) {
				Vector3d v = (mV.col(mF(0, p.first))
					+ mV.col(mF(1, p.first))
					+ mV.col(mF(2, p.first))) / 3.0f;
				glVertex3d(v.x(), v.y(), v.z());
			}
		}
		glEnd();
	}
}

void render_hierarchy()
{
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
				glVertex3dv(&mV(0, i));
				glVertex3dv(&mV(0, j));
			}
		}
		glEnd();
	}
}

static void render_quadmesh()
{
	if (show_quad) {
		glColor3f(1, 0, 0);
		glBegin(GL_LINES);
		for (auto& e : field.edge_idmap) {
			glVertex3dv(&field.mV_extracted(0, e.first.first));
			glVertex3dv(&field.mV_extracted(0, e.first.second));
		}
		glEnd();
	}
}

static void render_crossfield()
{
	if (show_field) {
		int l = (show_hierarchy) ? level : 0;
		auto& mV = field.hierarchy.mV[level];
		auto& mN = field.hierarchy.mN[level];
		auto& mQ = field.hierarchy.mQ[level];
		auto& adj = field.hierarchy.mAdj[level];
		glColor3f(1, 0, 0);
		double len = field.scale * 0.2;
		glBegin(GL_LINES);
		for (int i = 0; i < mQ.cols(); ++i) {
			glm::dvec3 p(mV(0, i), mV(1, i), mV(2, i));
			glm::dvec3 n(mN(0, i), mN(1, i), mN(2, i));
			glm::dvec3 tangent1 = glm::dvec3(mQ(0, i), mQ(1, i), mQ(2, i)) * len;
			glm::dvec3 tangent2 = glm::normalize(glm::cross(n, tangent1)) * len;
			auto l = p - tangent1;
			auto r = p + tangent1;
			auto u = p - tangent2;
			auto d = p + tangent2;
			glVertex3d(l.x, l.y, l.z);
			glVertex3d(r.x, r.y, r.z);
			glVertex3d(u.x, u.y, u.z);
			glVertex3d(d.x, d.y, d.z);
		}
		glEnd();
	}
}

static void render_loop() {
	if (show_loop) {
		glLineWidth(1.0f);
		glColor3f(1, 0, 0);
		glBegin(GL_LINES);
		for (auto& strip : field.edge_strips) {
			for (auto& e : strip) {
				Vector3d p = field.mV_extracted.col(field.qE[e].first);
				glVertex3d(p.x(), p.y(), p.z());
				p = field.mV_extracted.col(field.qE[e].second);
				glVertex3d(p.x(), p.y(), p.z());
			}
		}
		glEnd();
		
		glPointSize(3.0f);
		glColor3f(0, 0, 1);
		glBegin(GL_POINTS);
		for (int i = 0; i < field.sin_graph.size(); ++i) {
			if (field.vertex_singularities.count(i) == 0 && field.sin_graph[i].size() > 0) {
				Vector3d p = field.mV_extracted.col(i);
				glVertex3d(p.x(), p.y(), p.z());
			}
		}
		glEnd();
		
		glPointSize(10.0f);
		glColor3f(0, 1, 0);
		glBegin(GL_POINTS);
		for (auto& v : field.vertex_singularities) {
			Vector3d p = field.mV_extracted.col(v.first);
			glVertex3d(p.x(), p.y(), p.z());
		}
		glEnd();
	}
}

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

	render_mesh();
	render_singularities();
	render_hierarchy();
	render_quadmesh();
	render_crossfield();
	render_loop();

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

	if (key == 's') {
		show_singularity = 1 - show_singularity;
		glutPostRedisplay();
	}

	if (key == 'l') {
		show_loop = 1 - show_loop;
		glutPostRedisplay();
	}

	if (key == 'r') {
		for (int i = 0; i < 10; ++i) {
			field.LoopFace(1);
		}
		glutPostRedisplay();
	}

	if (key == 'c') {
		show_color = (show_color + 1) % 3;
		if (color.size() != field.hierarchy.mV[0].cols()) {
			color.resize(field.hierarchy.mV[0].cols());
		}
		if (show_color >= 1) {
			for (int i = 0; i < color.size(); ++i) {
				color[i] = Gray2HSV(field.hierarchy.mS[0](show_color - 1, i));
			}
		}
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

	Optimizer::optimize_scale(field.hierarchy);
	Optimizer::optimize_positions(field.hierarchy);
	field.ExtractMesh();
	/*
	printf("save\n");
	FILE* fp_w = fopen("result.txt", "wb");
	field.SaveToFile(fp_w);
	fclose(fp_w);
	printf("save finish\n");

	FILE* fp = fopen("result.txt", "rb");
	field.LoadFromFile(fp);
	fclose(fp);
	*/
	//	field.LoopFace(2);
	gldraw(mouse_callback, render_callback, motion_callback, keyboard_callback);
	return 0;
}