#include "gldraw.hpp"

#include <stdio.h> // printf
#include <stdlib.h> // exit
#include <string.h> // strlen
#include <vector>
#include "parametrizer.hpp"
#ifdef WITH_OPENGL
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <Eigen/Core>
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
int show_index = 0;
int show_v = -1;
int show_cuts;
int select_mode = 0;
int show_f;
std::vector<Vector3d> color;

extern Parametrizer field;

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

void render_cuts() {
    if (show_cuts == 0)
        return;
    auto& V = field.hierarchy.mV[0];
    glColor3f(1, 0, 0);
    glBegin(GL_LINES);
    for (auto& p : field.cuts) {
        glVertex3dv(&V(0, p.x));
        glVertex3dv(&V(0, p.y));
    }
    glEnd();
}

void render_test_travel(int f)
{
    auto& mF = field.hierarchy.mF;
    auto& mV = field.hierarchy.mV[0];
    auto& mN = field.hierarchy.mN[0];
    
    Vector3d p = mV.col(mF(0, f)) + mV.col(mF(1, f)) + mV.col(mF(2, f));
    p *= 1.0f / 3;
    
    const Vector3d &n = field.Nf.col(f);
    Vector3d q_x = field.FQ.col(f), q_y = n.cross(q_x);
    Vector3d q_xl = -q_x, q_xr = q_x;
    Vector3d q_yl_unfold = q_y, q_yr_unfold = q_y, q_xl_unfold = q_x, q_xr_unfold = q_x;
    double tx, ty, len;
    double step = field.scale;
    int i;
    i = f; len = step;
    Vector3d q1 = TravelField(p, q_xl, len, i, field.hierarchy.mE2E, mV, mF, field.Nf, field.FQ, field.hierarchy.mQ[0], mN, field.triangle_space, &tx, &ty, &q_yl_unfold);
    
    i = f; len = step;
    Vector3d q2 = TravelField(p, q_xr, len, i, field.hierarchy.mE2E, mV, mF, field.Nf, field.FQ, field.hierarchy.mQ[0], mN, field.triangle_space, &tx, &ty, &q_yr_unfold);
    
    printf("scale %lf\n", (q_yr_unfold - q_yl_unfold).dot(q_x) / (2 * step));
    glPointSize(2.0f);
    glBegin(GL_POINTS);
    glColor3f(1, 0, 0);
    glVertex3f(p.x(), p.y(), p.z());
    glColor3f(0, 1, 0);
    glVertex3f(q1.x(), q1.y(), q1.z());
    glVertex3f(q2.x(), q2.y(), q2.z());
    glEnd();
    
    Vector3d pl = p - q_x * step * 5;
    Vector3d pr = p + q_x * step * 5;
    Vector3d pu = p - q_y * step * 5;
    Vector3d pd = p + q_y * step * 5;
    
    Vector3d pul = p - q_yl_unfold * step * 5;
    Vector3d pdl = p + q_yl_unfold * step *  5;
    Vector3d pur = p - q_yr_unfold * step * 5;
    Vector3d pdr = p + q_yr_unfold * step * 0.5;
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    glColor3f(0, 0, 1);
    glVertex3dv(&pl.x());
    glVertex3dv(&pr.x());
    glColor3f(1, 0, 0);
    glVertex3dv(&pu.x());
    glVertex3dv(&pd.x());
    glColor3f(0, 1, 0);
    glVertex3dv(&pul.x());
    glVertex3dv(&pdl.x());
    glColor3f(1, 1, 0);
    glVertex3dv(&pur.x());
    glVertex3dv(&pdr.x());
    glEnd();
    
}

static void render_mesh()
{
    //    render_test_travel(348394);
    if (show_color == 0 && show_index == 0) {
        auto& mF = field.hierarchy.mF;
        auto& mV = field.hierarchy.mV[0];
        auto& mN = field.hierarchy.mN[0];
        
        glDisable(GL_LIGHTING);
        if (show_mesh) {
            //            static GLfloat white[4] =
            //            { 1.0, 1.0, 1.0, 1.0 };
            //            glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
            glBegin(GL_TRIANGLES);
            for (int i = 0; i < mF.cols(); ++i) {
                for (int j = 0; j < 3; ++j) {
                    int t1 = field.sing.count(mF(j, i));
                    int t2 = field.flipped.count(mF(j, i));
                    Vector3d c(0, 0, 0);
                    if (t1 == 1 && t2 == 0)
                        c[0] = 1;
                    else if (t1 == 1 && t2 == 1)
                        c[1] = 1;
                    else if (t1 == 0 && t2 == 1)
                        c[2] = 1;
                    glColor3f(c[0], c[1], c[2]);
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
    glPointSize(5.0f);
    if (show_quad) {
        glBegin(GL_POINTS);
        
        for (auto& p : field.singularities) {
            if (p.second == 1)
                glColor3f(0, 1, 0);
            else
                glColor3f(0, 0, 1);
            Vector3d v = (mV.col(mF(0, p.first))
                          + mV.col(mF(1, p.first))
                          + mV.col(mF(2, p.first))) / 3.0f;
            glVertex3d(v.x(), v.y(), v.z());
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
        auto& O = field.O_compact;
        auto& ED = field.edge_diff;
        auto& EV = field.edge_values;
        
        glPointSize(3.0);
        glColor3f(1, 0, 0);
        glBegin(GL_POINTS);
        for (int i = 0; i < O.size(); ++i) {
            glVertex3dv(&O[i][0]);
        }
        glEnd();
        glLineWidth(1.0);
        glBegin(GL_LINES);
        for (int i = 0; i < ED.size(); ++i) {
            if (abs(ED[i][0]) > 1 || abs(ED[i][1]) > 1) {
                printf("Long edge!\n");
                exit(0);
            }
            int l = abs(ED[i][0]) + abs(ED[i][1]);
            if (l > 1) {
                glColor3f(0, 0, 1);
            }
            else if (l == 1) {
                glColor3f(1, 0, 0);
            }
            else
                continue;
            if (field.fixed_cuts.count(
                                       DEdge(field.disajoint_tree.Index(EV[i].x),
                                             field.disajoint_tree.Index(EV[i].y)))) {
                                           //                glColor3f(0, 1, 0);
                                       }
            /*            if (field.fixed[field.disajoint_tree.Index(EV[i].x)] &&
             field.fixed[field.disajoint_tree.Index(EV[i].y)]) {
             glColor3f(0, 1, 0);
             }*/
            glVertex3dv(&O[field.disajoint_tree.Index(EV[i].x)][0]);
            glVertex3dv(&O[field.disajoint_tree.Index(EV[i].y)][0]);
        }
        glEnd();
        
        glColor3f(0.5, 0.5, 0.5);
/*        glBegin(GL_TRIANGLES);
        for (auto& p : field.flipped) {
            for (int j = 0; j < 3; ++j) {
                glVertex3dv(&O[p[j]][0]);
            }
        }
        glEnd();*/
    }
}

static void render_crossfield()
{
    if (show_field) {
        auto& mV = field.hierarchy.mV[level];
        auto& mN = field.hierarchy.mN[level];
        auto& mQ = field.hierarchy.mQ[level];
        
        double len = field.scale * 0.2f;
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
            glColor3f(1, 0, 0);
            glVertex3d(l.x, l.y, l.z);
            glVertex3d(r.x, r.y, r.z);
            glColor3f(0, 0, 1);
            glVertex3d(u.x, u.y, u.z);
            glVertex3d(d.x, d.y, d.z);
        }
        /*
         for (int i = 0; i < field.FQ.cols(); ++i) {
         Vector3d p = (mV.col(mF(0, i)) + mV.col(mF(1, i)) + mV.col(mF(2, i))) * (1.0 / 3.0);
         Vector3d t1 = field.FQ.col(i);
         Vector3d n = field.Nf.col(i);
         Vector3d t2 = n.cross(t1);
         auto l = p - t1 * len;
         auto r = p + t1 * len;
         auto u = p - t2 * len;
         auto d = p + t2 * len;
         glVertex3d(l.x(), l.y(), l.z());
         glVertex3d(r.x(), r.y(), r.z());
         glVertex3d(u.x(), u.y(), u.z());
         glVertex3d(d.x(), d.y(), d.z());
         }
         */
        glEnd();
    }
}

static void render_callback(void)
{
    glClearColor(0.0, 191.0 / 255.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1, 1, -1, 1, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    GLfloat light_position[] = { 1, 1, 1, 0 };
    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, light_position);
    // now you can setup view matrix (gluLookAt())
    
    glPushMatrix();
    glTranslated(render_translation.x, render_translation.y, render_translation.z);
    double model_scale = exp(render_scale);
    glScaled(model_scale, model_scale, 0.9);
    glMultMatrixd((double*)&render_rotation);
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
    render_cuts();
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
    if (key == 'v' || key == 'b') {
        show_v = show_v + (key == 'b' ? 1 : -1);
        printf("show_v: %d\n", show_v);
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
        show_f = 1 - show_f;
        glutPostRedisplay();
    }
    
    if (key == 'm') {
        show_mesh = 1 - show_mesh;
        glutPostRedisplay();
    }
    
    if (key == 'u') {
        show_cuts = 1 - show_cuts;
        glutPostRedisplay();
    }
    
    if (key == 'c') {
        if (show_v != -1) {
            select_mode = 0;
            show_v = -1;
        }
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
    
    if (key == 'a') {
        select_mode = 1 - select_mode;
        glutPostRedisplay();
    }
    
    if (key == 'l') {
        show_loop = 1 - show_loop;
        glutPostRedisplay();
    }
    
    if (key == 'c') {
        show_color = (show_color + 1) % 3;
        if (color.size() != field.hierarchy.mV[0].cols()) {
            color.resize(field.hierarchy.mV[0].cols());
        }
        if (show_color >= 1) {
            double minX = 1e30, maxX = -1e30;
            for (int i = 0; i < color.size(); ++i) {
                minX = std::min(minX, field.hierarchy.mS[0](show_color - 1, i));
                maxX = std::max(maxX, field.hierarchy.mS[0](show_color - 1, i));
            }
            for (int i = 0; i < color.size(); ++i) {
                double t = field.hierarchy.mS[0](show_color - 1, i);
                t = (t - minX) / (maxX - minX);
                color[i] = Vector3d(t, t, t);//Gray2HSV(t);
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

void gldraw()
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

	glutPostRedisplay();

	glutMainLoop();

	return;
}

#endif
