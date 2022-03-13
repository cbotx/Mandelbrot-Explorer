#ifndef __RENDERER_H__
#define __RENDERER_H__

class Renderer {
public:
    Renderer();
    void outputImage(int w, int h, int sub, float* iter, char fname[]);
    void mixColor(uint8_t* img, float* iter, std::array<int, 2> p, int w) const;
    void setColor(uint8_t* img, int pos, float iteration) const;
    void setColor(uint8_t* img, int pos, uint8_t r, uint8_t g, uint8_t b) const;
    inline int getIndex(std::array<int, 4> arr, int w, int sub) const;
};


#endif