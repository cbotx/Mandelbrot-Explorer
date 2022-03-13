#ifndef __GLUT_CANVAS_H__
#define __GLUT_CANVAS_H__


#include <GL/glew.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <array>
#include <vector>

class GlutCanvas {
private:
    int _w, _h;
public:
    typedef GLuint Texture;

    GlutCanvas();

    static void GlutInit(int argc, char** argv);

    int GlutCreateWindow(int width, int height, int left, int top, const char* name);

    void GlutBindDisplayFunc(void (*func)(void));

    void GlutBindReshapeFunc(void (*func)(int, int));

    void GlutBindMouseFunc(void (*func)(int, int, int, int));

    void GlutBindMotionFunc(void (*func)(int, int));

    void GlutBindKeyboardFunc(void (*func)(unsigned char, int, int));

    void GlutSetRenderInterval(int interval);

    void GlutResetWindowSize();

    void BeginDraw();

    void EndDraw();

    void GlutRun();

    Texture CreateTexture(bool interpolation=true);

    void BitmapToTexture(const Texture &texture, const uint8_t* bitmap);

    void BitmapToTexture(const Texture &texture, const uint8_t* bitmap, const int width, const int height);

    void RenderTexture(const Texture &texture, std::array<std::array<double, 2>, 4> src_bound, std::array<std::array<double, 2>, 4> tgt_bound);

    void RenderTextureP(const Texture &texture, std::array<std::array<double, 2>, 4> src_bound, std::array<std::array<int, 2>, 4> tgt_bound);

    void DrawPolygon(std::vector<std::array<float, 2>>& bound, const float r, const float g, const float b, bool fill = false);

    void RenderText(const int x, const int y, const float r, const float g, const float b, const unsigned char* text);

    void ScreenToBitmap(uint8_t* bitmap);

    void Redisplay();

    void Redisplay(int window);
    
};


#endif