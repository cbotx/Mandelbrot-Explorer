#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include <mpir.h>
#include <chrono>
#include <bitset>

#include "image.h"
#include "mandel_perturbation.h"
#include "test_cases.h"
#include "renderer.h"


int main() {
    // std::freopen("accurate.txt", "w", stdout);
    srand(0);
    colorMapInitialize();
    static const bool SHOW_DEEPEST = false;
    static const bool SHOW_CENTER = false;
    int width = 300;
    int height = 200;
    int max_iteration = 1000000;
    const char* const& xc = "0";
    const char* const& yc = "0";
    std::string scale_s = "1";
    // const char* const& xc = testcases::deep1_x;
    // const char* const& yc = testcases::deep1_y;
    // std::string scale_s = testcases::deep1_z;
    // scale_s = scale_s.substr(0, 860);
    int sub = 1;
    float* iter = new float[width * height * sub * sub];
    Mandel mandel(width, height, max_iteration, sub, iter);
    Renderer renderer;

#ifdef AUTO_ZOOM
    while(true) {
#endif
        auto t_begin = std::chrono::high_resolution_clock::now();
        int precision = static_cast<int>(scale_s.size() * log(10) / log(2)) + 30;
        std::cout << "Precision: " << precision << std::endl;
        
        // HPComp c{HPFloat(xc, precision), HPFloat(yc, precision)};
        mpf_t c_re, c_im, scale;
        mpf_set_default_prec(precision);
        mandel.setPrecision(precision);
        mpf_init_set_str(c_re, xc, 10);
        mpf_init_set_str(c_im, yc, 10);
        // HPComp p = c;
        // HPFloat scale{ scale_s, precision };
        mpf_init_set_str(scale, scale_s.c_str(), 10);
        mandel.Compute(c_re, c_im, scale, max_iteration);
        auto t_end = std::chrono::high_resolution_clock::now();
        auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_begin);
        std::cout << std::setprecision(5) << "Time: " << time_span.count() / 1000.0 << "s\n";
        
        char fname[256] = "test3.png";
        renderer.outputImage(width, height, sub, iter, fname);
 #ifdef AUTO_ZOOM
        scale_s.push_back('0');
        c = p;
        // std::cout << std::setprecision(1300) << c.get_real_imag().first << '\n' << c.get_real_imag().second << "i\n" << scale << "\n";
    }
#endif
}