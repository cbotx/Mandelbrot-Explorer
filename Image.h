#ifndef __IMAGE_H__
#define __IMAGE_H__

#include "png++/png.hpp"

static const int colP = 2048;
extern float color_map[3][colP];
extern float color_density;

void writePPMImage(int* data, int width, int height, const char* filename, int maxIterations);

void writePNGImage(uint8_t* img, int width, int height, const char* filename);

void colorMapInitialize();

void getColor(float iteration, uint8_t& r, uint8_t& g, uint8_t& b, int c_method=0);

void getColor(float iteration, float& fr, float& fg, float& fb, int c_method=0);

void rgbRotate(float& fr, float& fg, float& fb, float rad);

#endif