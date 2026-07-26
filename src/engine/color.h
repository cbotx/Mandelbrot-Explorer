#ifndef __COLOR_H__
#define __COLOR_H__


void lab2rgb( float l_s, float a_s, float b_s, float& R, float& G, float& B );

void rgb2lab( float R, float G, float B, float & l_s, float &a_s, float &b_s);

void hsv2rgb(float H, float S, float V, float& R, float& G, float& B);

void rgb2hsv(float R, float G, float B, float& H, float& S, float& V);

#endif
