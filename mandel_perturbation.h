#ifndef __MANDEL_PERTURBATION_H__
#define __MANDEL_PERTURBATION_H__

#include <complex>
#include <gmp.h>

#include "floatexp.h"

#include <set>
#include <algorithm>
#include <array>
#include <vector>

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
    typedef double Float;
    typedef std::complex<Float> Comp;
public:
    Mandel(int width, int height, int max_iteration, int sub, float* iter);
    virtual ~Mandel();

    void Compute(mpf_t c_re, mpf_t c_im, mpf_t scale, int mxit, int c_method = 0);
    // Verification helper: brute-force high-precision escape time into out[],
    // reusing the grid (_c0/_dx/_dy) set by the most recent Compute() call.
    // step>1 samples only pixels on a (step x step) grid (others left untouched)
    // so the O(mxit) full-precision oracle stays feasible at extreme depth.
    // Requires sub == 1.
    void ComputeDirect(int mxit, float* out, int step = 1, int c_method = 0);
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
    // AVX2 kernel: iterate a group of up to 4 pixels (non-EDE path). Mirrors the
    // scalar delta loop op-for-op so results are bit-identical. lanes<=4.
    void solveSimd4(const std::array<int, 4>* v, int g, int lanes,
                    const double* dcr, const double* dci,
                    const double* dzr, const double* dzi,
                    int mx_ref_it, int mxit, Float* glitch_p);

    // AVX2 shallow (non-perturbation) Mandelbrot with lane refilling: computes a
    // LIST of `count` pixels (arbitrary c = cre[k] + i*cim[k]), writing
    // floatPointCompute-equivalent values into out[0..count). Used for the base row
    // and the adaptive supersample sub-pixels.
    void solveShallowSimdList(const double* cre, const double* cim, int count,
                              float* out, int mxit, int c_method) const;

    // Bivariate Linear Approximation (Zhuoran / Fraktaler-3). A BLA skips l
    // reference iterations via dz -> A*dz + B*dc when |dz| < R. The hot entry
    // (A, B, r2, l) is what every tryBLA/tryBLAfe lookup reads; splitting it from
    // the EDE-only derivative couplings (C, E in a parallel _blaD) shrinks the
    // per-entry footprint ~40% (80->48 B) so the large table streams from cache
    // better in the memory-bound delta loop.
    struct BLAEntry { double ar, ai, br, bi, r2; int l; };
    // Derivative deltas carried through a skip under EDE: dd -> C*dz + A*dd + E*dc.
    struct BLADeriv { double cr, ci, er, ei; };
    std::vector<std::vector<BLAEntry>> _bla;
    std::vector<std::vector<BLADeriv>> _blaD;   // parallel to _bla; only built under EDE
    double _bla_eps = 0.0;
    double _bla_rmax2 = 0.0;      // largest BLA validity radius^2 (gate for tryBLA)
    bool _use_bla = false;
    int _bla_minlevel = 0;        // only apply BLAs with skip >= 2^_bla_minlevel
    bool _ref_bounded = false;    // reference orbit never escapes (minibrot center)
    bool _use_interior = false;   // periodicity-based interior detection
    double _interior_eps2 = 0.0;  // |z_n - z_saved|^2 threshold for a detected cycle
    int _interior_confirm = 30;   // consecutive shrinking periods required to confirm
    void buildBLA(int reflen, bool ede);
    // Try the largest valid BLA starting at reference index s; on success updates
    // dz (and, under ede, dd) and returns the skip (>0), else returns 0. Rejects
    // skips that would land escaped or past mx_ref_it.
    int tryBLA(int s, double& dzr, double& dzi, double& ddr, double& ddi,
               double dcr, double dci, bool ede, double ESC2, int mx_ref_it) const;
    // Floatexp BLA for the rescaled path: dz is stored as dz = S*(wr,wi). On a
    // valid skip at reference index s it advances dz by the entry's length
    // (updating S/wr/wi) and returns the skip (>0), else 0. dz-only (the rescaled
    // kernel returns smooth iteration, so no EDE derivative is carried).
    int tryBLAfe(int s, FloatExp& S, FloatExp S2, double& wr, double& wi,
                 FloatExp dcr, FloatExp dci, double ESC2, int mx_ref_it) const;
    int createRef(std::set<std::array<int, 4>>& s, int pr_it, int mxit, bool random,
                  int c_method = 0, bool view_center = false);
    bool calCoefficient(int i, int pr_it, int c_method = 0);
    // Deep-zoom rescaled perturbation (Zhuoran z = S w): the delta w stays an
    // O(1) double while the floatexp scale S carries the deep exponent, so the
    // inner loop is native-double yet correct far past double's ~1e320 underflow.
    // Used when _use_floatexp; returns the pixel escape value (interior -> -2).
    float pixelRescaled(FloatExp dcr, FloatExp dci, int mx_ref_it, int mxit, int c_method) const;
    // AVX2 4-wide version of pixelRescaled: iterates up to 4 deep-zoom pixels in
    // lockstep, vectorising the rescaled double step while keeping BLA / rescale /
    // rebase as cheap per-lane scalar events. Writes smooth-iteration values into
    // out[g..g+lanes). Mirrors pixelRescaled's math; no EDE derivative (dz-only).
    void solveRescaledSimd4(const FloatExp* Dcr, const FloatExp* Dci, int g, int lanes,
                            int mx_ref_it, int mxit, int c_method, float* out) const;
    int SACheckMagnitude() const;
    float accuratePointCompute(mpf_t c_re, mpf_t c_im, int mxit, int c_method = 0) const;
    double floatPointCompute(Float c_re, Float c_im, int mxit, int c_method = 0) const;
    bool attractor(double z_in_re, double z_in_im, const double c_re, const double c_im, int period) const;

private:
    volatile bool _flag_halt = false;
    volatile bool _sub_flag;
    float _ESCAPE_RADIUS = 1e8f;
    const double _TOL = 1e12;
    static const int _SA_N = 30;

    int _sub;
    float* _iter;
    const int _w, _h;
    const int _mxit;
    mpf_t* _z_re, * _z_im;  // reference orbit
    Comp* _zf;

    // Exterior DE
    mpf_t* _d_re, * _d_im;
    Comp* _df;
    Float* _dfr, * _dfi;
    // Reference derivative D = dZ/dc carried in the reference build. It only feeds
    // the double shadow _dfr/_dfi (the delta loop + BLA never read the full-
    // precision mpf derivative), so it is iterated directly here instead of in
    // mpf: floatexp for the deep path (Z underflows the double shadow at a
    // minibrot's near-zero passes) and plain double otherwise. Saves 3 mpf muls
    // per reference iteration.
    FloatExp _dfe_r{ 1.0, 0 }, _dfe_i{ 0.0, 0 };

    // Series approximation context
    Comp _Adf_old[_SA_N], _Adf_new[_SA_N];
    Comp _Bdf_old[_SA_N], _Bdf_new[_SA_N];
    Comp _SA_delta;
    bool _SA_flag;
    int _SA_order;
    int _SA_it;

    Float* _zfr, * _zfi;
    // floatexp reference orbit (filled only when _use_floatexp): accurate where
    // |Z| underflows the double shadow _zfr/_zfi, needed for deep-zoom rebasing.
    FloatExp* _zfr_fe = nullptr, * _zfi_fe = nullptr;
    bool _use_floatexp = false;
    bool _fe_cutoff_sensitive = false;
    FloatExp _dxfe{ 0.0, 0 }, _dyfe{ 0.0, 0 };   // pixel spacing in floatexp
    mpf_t _c0_re, _c0_im;  // bottom left point coordinate
    mpf_t _dx;
    mpf_t _dy;
    mpf_t _scale;
    std::array<int, 4> _ref;  // reference pixel index
    double _ref_x = 0.0, _ref_y = 0.0;  // reference position in pixel coordinates
    bool _ref_virtual = false;           // exact view center, not an actual pixel
    std::array<int, 4> _new_ref;  // new reference pixel index
    mpf_t _ref_z_re, _ref_z_im;  // reference point coordinate
    Comp _ref_z_f;
    int _ref_cnt;
    bool* _done;
    Float* _z_m3;
    int _mx_coef;

    // temporary vavriables
    mpf_t _t1, _t2, _t3, _t4;
};

#endif