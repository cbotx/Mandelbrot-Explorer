#include <array>
#include <math.h>

#include "image.h"
#include "renderer.h"

Renderer::Renderer() {
    
}

void Renderer::outputImage(int w, int h, int sub, float* iter, char fname[]) {
    uint8_t* img = new uint8_t[w * h * sub * sub * 3];
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            if (sub == 1) {
                setColor(img, i * w + j, iter[i * w + j]);
            } else {
                if (iter[getIndex({i, j, -1, -1}, w, sub)] == -3) {
                    setColor(img, i * w + j, iter[getIndex({i, j, 0, 0}, w, sub)]);
                } else {
                    mixColor(img, iter, {i, j}, w);
                }
            }
        }
    }
    writePNGImage(img, w, h, fname);
}

void Renderer::mixColor(uint8_t* img, float* iter, std::array<int, 2> p, int w) const {
    float rs, gs, bs;
    float ws = 0;

    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            std::array<int, 4> arr = { p[0], p[1], i, j };
            float r, g, b;
            getColor(iter[getIndex(arr, w, 3)], r, g, b);

            rs += r * r;
            gs += g * g;
            bs += b * b;
            ws += 1;
        }
    }
    rs /= ws;
    gs /= ws;
    bs /= ws;
    setColor(img, p[0] * w + p[1], sqrt(rs), sqrt(gs), sqrt(bs));
}


void Renderer::setColor(uint8_t* img, int pos, float iteration) const {
    getColor(iteration, img[pos * 3 + 0], img[pos * 3 + 1], img[pos * 3 + 2]);
}

void Renderer::setColor(uint8_t* img, int pos, uint8_t r, uint8_t g, uint8_t b) const {
    img[pos * 3 + 0] = r;
    img[pos * 3 + 1] = g;
    img[pos * 3 + 2] = b;
}

inline int Renderer::getIndex(std::array<int, 4> arr, int w, int sub) const {
    return (arr[0] * sub + sub / 2 + arr[2]) * w * sub + (arr[1] * sub + sub / 2 + arr[3]);
}