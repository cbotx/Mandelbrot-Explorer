#include <GL/glew.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <GL/freeglut.h>

#include <array>

#include "glut_canvas.h"

static int timer_interval;

void glut_timer_func(int v) {
    glutPostRedisplay();
    glutTimerFunc(timer_interval, glut_timer_func, v);
}

GlutCanvas::GlutCanvas() {
}

void GlutCanvas::GlutInit(int argc, char** argv) {
    glutInit(&argc, argv);
}

int GlutCanvas::GlutCreateWindow(int width, int height, int left, int top, const char* name) {
    _w = width;
    _h = height;
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(width, height);
    glutInitWindowPosition(left, top);
    return glutCreateWindow(name);
}

void GlutCanvas::GlutBindDisplayFunc(void (*func)(void)) {
    glutDisplayFunc(func);
}

void GlutCanvas::GlutBindReshapeFunc(void (*func)(int, int)) {
    glutReshapeFunc(func);
}

void GlutCanvas::GlutBindMouseFunc(void (*func)(int, int, int, int)) {
    glutMouseFunc(func);
}

void GlutCanvas::GlutBindKeyboardFunc(void (*func)(unsigned char, int, int)) {
    glutKeyboardFunc(func);
}

void GlutCanvas::GlutBindMotionFunc(void (*func)(int, int)) {
    glutMotionFunc(func);
}

void GlutCanvas::GlutSetRenderInterval(int interval) {
    timer_interval = interval;
    glutTimerFunc(timer_interval, glut_timer_func, 0);
}

void GlutCanvas::GlutResetWindowSize() {
    glutReshapeWindow(_w, _h);
}

void GlutCanvas::GlutRun() {
    glutMainLoop();
}

void GlutCanvas::BeginDraw() {
    glClearColor(1, 1, 1, 1);
    glClear(GL_COLOR_BUFFER_BIT);
}

void GlutCanvas::EndDraw() {
    glFlush();
    glutSwapBuffers();
}

GlutCanvas::Texture GlutCanvas::CreateTexture(bool interpolation) {
    Texture texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    
    const float borderColor[] = { 0.f, 0.f, 0.f, 1.f };
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

    if (interpolation) {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    } else {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    }
    glBindTexture(GL_TEXTURE_2D, 0);
    return texture;
}

void GlutCanvas::BitmapToTexture(const Texture &texture, const uint8_t* bitmap) {
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, _w, _h, 0, GL_RGB, GL_UNSIGNED_BYTE, bitmap);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void GlutCanvas::BitmapToTexture(const Texture &texture, const uint8_t* bitmap, const int width, const int height) {
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, bitmap);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void GlutCanvas::RenderTexture(const Texture &texture, std::array<std::array<double, 2>, 4> src_bound, std::array<std::array<double, 2>, 4> tgt_bound) {
    glBindTexture(GL_TEXTURE_2D, texture);
    glEnable(GL_TEXTURE_2D);
    glBegin(GL_QUADS);
    for (int i = 0; i < 4; ++i) {
        glTexCoord2f(src_bound[i][0], src_bound[i][1]);
        glVertex2f(tgt_bound[i][0], tgt_bound[i][1]);
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void GlutCanvas::RenderTextureP(const Texture &texture, std::array<std::array<double, 2>, 4> src_bound, std::array<std::array<int, 2>, 4> tgt_bound) {
    glColor3f(1, 1, 1);
    glBindTexture(GL_TEXTURE_2D, texture);
    glEnable(GL_TEXTURE_2D);
    glBegin(GL_QUADS);
    for (int i = 0; i < 4; ++i) {
        glTexCoord2f(src_bound[i][0], src_bound[i][1]);
        glVertex2f(2.0 * tgt_bound[i][0] / _w - 1.0, 2.0 * tgt_bound[i][1] / _h - 1.0);
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);
}


void GlutCanvas::DrawPolygon(std::vector<std::array<float, 2>>& bound, const float r, const float g, const float b, bool fill) {
    glColor3f(r, g, b);
    glBegin(fill ? GL_POLYGON : GL_LINE_LOOP);
    for (auto& p : bound) {
        glVertex2f(2.0 * p[0] / _w - 1.0, 2.0 * p[1] / _h - 1.0);
    }
    glEnd();
}

void GlutCanvas::RenderText(const int x, const int y, const float r, const float g, const float b, const unsigned char* text) {
    glColor3f(r, g, b);
    glRasterPos2f(2.0 * x / _w - 1.0, 2.0 * y / _h - 1.0);
    glutBitmapString(GLUT_BITMAP_HELVETICA_18, text);
}

void GlutCanvas::ScreenToBitmap(uint8_t* bitmap) {
    glReadPixels(0, 0, _w, _h, GL_RGB, GL_UNSIGNED_BYTE, bitmap);
}

void GlutCanvas::Redisplay() {
    glutPostRedisplay();
}

void GlutCanvas::Redisplay(int window) {
    glutSetWindow(window);
    glutPostRedisplay();
}