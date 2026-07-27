#ifndef __IMAGE_H__
#define __IMAGE_H__

#include <cstdint>

static const int colP = 2048;
extern float color_map[3][colP];
extern float color_density;
// Palette-cycle phase in palette-index units [0, colP). Added to every colour
// lookup so the gradient can be rotated (animated) without recomputing the
// fractal. Wraps modulo colP.
extern float color_phase;

void writePPMImage(int* data, int width, int height, const char* filename, int maxIterations);

void writePNGImage(uint8_t* img, int width, int height, const char* filename);

void colorMapInitialize();

void getColor(float iteration, uint8_t& r, uint8_t& g, uint8_t& b, int c_method=0);

void getColor(float iteration, float& fr, float& fg, float& fb, int c_method=0);

// Analytic anti-aliasing: rebuild the palette box-filter integral (call once per
// frame before AA coloring), then colour a pixel as the average colour over the
// palette range its neighbourhood spans -- approximating infinite supersampling
// at one sample/pixel. vL/vR/vU/vD are the neighbouring smooth values (pass
// EMPTYPIXEL / a negative value for missing/interior neighbours).
void prepareColorFilter();
void getColorAA(float v, float vL, float vR, float vU, float vD,
                uint8_t& r, uint8_t& g, uint8_t& b, int c_method=0);

// Split of the analytic-AA colouring into a phase-INDEPENDENT analysis and a
// cheap phase-DEPENDENT shading, so an animated palette phase can re-shade the
// cached per-pixel (baseU, width) each frame without re-running colorFunction or
// the neighbourhood gradient. getColorAA == colorShadeAA(colorAnalyzeAA(...), color_phase).
//   baseU  : palette-index centre (colorFunction*colP), or <0 for interior/empty.
//   width  : palette entries spanned by the pixel footprint (0 => point sample).
void colorAnalyzeAA(float v, float vL, float vR, float vU, float vD, int c_method,
                    float& baseU, float& width);
void colorShadeAA(float baseU, float width, float phase, uint8_t& r, uint8_t& g, uint8_t& b);
// Phase-independent base palette index for a point sample (Feather / non-AA);
// returns <0 for interior/empty pixels.
float colorBaseIndex(float iteration, int c_method);

// Gamma-correct (sRGB) averaging: colours are averaged in LINEAR light so a
// downsampled pixel matches what the eye sees "from far away" (radiance averages
// linearly). g_srgb2lin maps an 8-bit sRGB value to linear [0,1]; srgbEncode maps
// linear [0,1] back to an sRGB [0,255] value.
extern float g_srgb2lin[256];
float srgbEncode(double lin);

void rgbRotate(float& fr, float& fg, float& fb, float rad);

#endif