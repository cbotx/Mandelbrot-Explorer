#ifndef __MANDEL_PERTURBATION_H__
#define __MANDEL_PERTURBATION_H__

#include <complex>
#include <mpir.h>

#include <set>
#include <algorithm>
#include <array>

constexpr float EMPTYPIXEL = -10;

struct HPComp {
    mpf_t re;
    mpf_t im;
};

enum ColoringMethod {
    SUPER_SAMPLING = 1,
    EXTERIOR_DIST_EST = 2,
    INTERIOR_DIST_EST = 4
};

class Mandel {
public:
    typedef long double Float;
    typedef std::complex<Float> Comp;
public:
    Mandel(int width, int height, int max_iteration, int sub, float* iter);
    virtual ~Mandel();

    void Compute(mpf_t c_re, mpf_t c_im, mpf_t scale, int mxit, int c_method = 0);
    void Output(char fname[]) const;
    void setPrecision(int precision);
    HPComp getHighIterationPoint() const;
    HPComp getSymmetryCenter() const;
    HPComp getRandomZoomPoint() const;
    uint8_t* getImage() const;

    void SetHalt(bool flag);

private:
    inline bool escape(Comp& z) const;
    inline bool escape(mpf_t& z_re, mpf_t& z_im) const;
    inline float getEscapeTime(Comp& z, int i) const;
    inline float getEscapeTime(mpf_t& z_re, mpf_t& z_im, int i);

    inline int getIndex(std::array<int, 4>& arr) const;
    inline int getIndex(int i, int j, int u, int v) const;

    void setPixel(std::array<int, 4> p, float iteration) const;
    void stepParallel(std::set<std::array<int, 4>>& s, int mx_ref_it, int mxit, int c_method = 0);
    int createRef(std::set<std::array<int, 4>>& s, int pr_it, int mxit, bool random, int c_method = 0);
    bool calCoefficient(int i, int pr_it, int c_method = 0);
    int SACheckMagnitude() const;
    float accuratePointCompute(mpf_t c_re, mpf_t c_im, int mxit, int c_method = 0) const;
    double floatPointCompute(Float c_re, Float c_im, int mxit, int c_method = 0) const;
    bool attractor(double z_in_re, double z_in_im, const double c_re, const double c_im, int period) const;

private:
    volatile bool _flag_halt = false;
    volatile bool _sub_flag;
    const float _ESCAPE_RADIUS = 100;
    const double _TOL = 1e12;
    static const int _SA_N = 30;

    int _sub;
    float* _iter;
    const int _w, _h;
    const int _mxit;
    mpf_t* _z_re, * _z_im;  // reference orbit
    Comp* _zf, * _delta_0;

    // Exterior DE
    mpf_t* _d_re, * _d_im;
    Comp* _df;
    Float* _dfr, * _dfi;

    // Series approximation context
    Comp _Adf_old[_SA_N], _Adf_new[_SA_N];
    Comp _Bdf_old[_SA_N], _Bdf_new[_SA_N];
    Comp _SA_delta;
    bool _SA_flag;
    int _SA_order;
    int _SA_it;

    Float* _zfr, * _zfi;
    mpf_t _c0_re, _c0_im;  // bottom left point coordinate
    mpf_t _dx;
    mpf_t _dy;
    mpf_t _scale;
    std::array<int, 4> _ref;  // reference pixel index
    std::array<int, 4> _new_ref;  // new reference pixel index
    mpf_t _ref_z_re, _ref_z_im;  // reference point coordinate
    Comp _ref_z_f;
    int _ref_cnt;
    bool* _done;
    Float* _z_m3;
    int _mx_coef;

    // temporary vavriables
    mpf_t _t1, _t2;
};

#endif