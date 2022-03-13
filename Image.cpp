#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>

#include "image.h"
#include "interpolate.h"
#include "mandel_perturbation.h"

#define ADDBORDER 1


float color_map[3][colP];
float color_density = 60;

static inline float color_func(float it, int c_method=0) {
    if (c_method & ColoringMethod::EXTERIOR_DIST_EST)
        return tanh(it / color_density * 5);
    return pow(log(it + 2), log(it + 2) * log(it + 2) / color_density);
}

void getColor(float iteration, uint8_t& r, uint8_t& g, uint8_t& b, int c_method) {
    float fr, fg, fb;
    getColor(iteration, fr, fg, fb, c_method);
    r = (int)fr;
    g = (int)fg;
    b = (int)fb;
}

void getColor(float iteration, float& fr, float& fg, float& fb, int c_method) {
    if (iteration == -3) {
        fr = 255.f;
        fg = 0.f;
        fb = 0.f;
    } else if (iteration == -2) {
        fr = 0.f;
        fg = 0.f;
        fb = 0.f;
    } else {
        float f_it = color_func(iteration, c_method);
        int x = (int)(f_it * colP) % colP;
        float rad = f_it / 100;
        fr = color_map[0][x];
        fg = color_map[1][x];
        fb = color_map[2][x];
        // rgbRotate(fr, fg, fb, rad);
    }
}

void rgbRotate(float& fr, float& fg, float& fb, float rad) {
    float c = cos(rad);
    float s = sin(rad);
    float th = 1.f / 3.f;
    float sqth = sqrt(th);
    float m1 = c + th * (1.f - c);
    float m2 = th * (1.f - c) - sqth * s;
    float m3 = th * (1.f - c) + sqth * s;
    float rx = fr * m1 + fg * m2 + fb * m3;
    float gx = fr * m3 + fg * m1 + fb * m2;
    float bx = fr * m2 + fg * m3 + fb * m1;
    fr = std::max(0.f, std::min(255.f, rx));
    fg = std::max(0.f, std::min(255.f, gx));
    fb = std::max(0.f, std::min(255.f, bx));
}

void colorMapInitialize() {
    static const int N = 6;
    static const float pos[N] = { 0.0, 0.16, 0.42, 0.6425, 0.8575, 1.0 };
    static const float col[3][N] = { {0, 32, 237, 255, 0, 0}, {70, 107, 255, 170, 2, 70}, {100, 203, 255, 0, 0, 100} };
    for (int i = 0; i < 3; ++i) {
        mono_cubic_interpolate(pos, col[i], N, color_map[i], colP);
    }
}

void writePNGImage(uint8_t* data, int width, int height, const char* filename) {
    png::image<png::rgb_pixel> image(width, height);
    for (png::uint_32 y = 0; y < image.get_height(); ++y) {
        for (png::uint_32 x = 0; x < image.get_width(); ++x) {
            int idx = (y * width + x) * 3;
            image[y][x] = png::rgb_pixel(data[idx + 0], data[idx + 1], data[idx + 2]);
        }
    }
    image.write(filename);
}