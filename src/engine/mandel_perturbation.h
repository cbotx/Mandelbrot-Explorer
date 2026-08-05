#ifndef __MANDEL_PERTURBATION_H__
#define __MANDEL_PERTURBATION_H__

#include <complex>
#include <gmp.h>

#include "floatexp.h"
#include "bigfixed.h"
#include "formula_spec.h"
#include "formula_expression.h"

#include <set>
#include <algorithm>
#include <atomic>
#include <array>
#include <vector>

namespace formula { class ExpressionJit4; }

constexpr float EMPTYPIXEL = -10;

struct HPComp {
    mpf_t re;
    mpf_t im;
};

enum ColoringMethod {
    SUPER_SAMPLING = 1,
    EXTERIOR_DIST_EST = 2,
    INTERIOR_DIST_EST = 4,
    STRIPE_AVERAGE = 8,
    // Analytic normal-map lighting: outputs the smooth escape value into `_iter`
    // (as method 0) AND the surface-normal angle arg(z) - arg(dz/dc) into a second
    // buffer `_normal`, from the same derivative EDE tracks. The GUI shades the
    // base colour by a Lambert light using that normal (light stays animatable).
    NORMAL_MAP = 16,
    // Orbit-trap colouring: tracks the orbit's closest approach to a composite
    // trap (point + Pickover cross + ring) and outputs a palette coordinate. Needs
    // the full orbit (BLA is disabled), so it is slow at deep zoom.
    ORBIT_TRAP = 32,
    // Smooth colour + exterior distance estimate overlay: outputs the smooth value
    // into `_iter` (base colour) AND the pixel-normalised distance estimate into
    // `_normal`, so the GUI can draw the B&W filament layer on top of the smooth
    // colouring (Wikipedia's "DE + smooth"). Reuses the EDE derivative.
    DE_OVERLAY = 64
};

class Mandel {
public:
    typedef double Float;
    typedef std::complex<Float> Comp;
public:
    Mandel(int width, int height, int max_iteration, int sub, float* iter);
    virtual ~Mandel();

    void Compute(mpf_t c_re, mpf_t c_im, mpf_t scale, int mxit, int c_method = 0,
                 int full_h = 0, int row_base = 0);
    // Quadratic Julia z0-plane: each pixel supplies z0 while fixedC is constant.
    // Phase 1 is a dedicated direct AVX2 kernel (sub==1, Smooth/EDE); advanced
    // perturbation/SA/BLA capabilities remain disabled in quadraticJulia().
    void ComputeJulia(mpf_t z0_re, mpf_t z0_im, mpf_t scale,
                      mpf_t fixed_c_re, mpf_t fixed_c_im,
                      int mxit, int c_method = 0);
    // Generic direct-expression backend. Each pixel can bind either c or z0;
    // other context values/parameters stay fixed. Returns false for invalid
    // bytecode/bindings. Output is raw escape iteration (no formula-specific
    // smoothing/DE/perturbation claims).
    bool ComputeExpression(mpf_t center_re, mpf_t center_im, mpf_t scale,
                           const formula::ExpressionProgram& program,
                           const formula::ExpressionContext& fixed,
                           FormulaParameter pixelParameter,
                           int mxit, double bailout = 4.0,
                           formula::ExpressionColoring coloring =
                               formula::ExpressionColoring::Raw,
                           const formula::ExpressionJit4* jit = nullptr);
    // Exact residual perturbation prototype:
    // dz' = F(Z + dz, pixel parameters) - Z'. Falls back to the direct
    // expression renderer when the center reference escapes before mxit.
    bool ComputeExpressionResidual(mpf_t center_re, mpf_t center_im, mpf_t scale,
                                   const formula::ExpressionProgram& program,
                                   const formula::ExpressionContext& fixed,
                                   FormulaParameter pixelParameter,
                                   int mxit, double bailout,
                                   bool* usedPerturbation = nullptr,
                                   formula::ExpressionColoring coloring =
                                       formula::ExpressionColoring::Raw);
    // Exact shallow scalar used to resolve the small set of numerically
    // sensitive GPU pixels without re-running the full frame on the CPU.
    float ComputeShallowPoint(double cRe, double cIm, int mxit) const;
    // Batched equivalent used by the GPU tail: AVX2 lane refill preserves the
    // shallow CPU result while avoiding one scalar call per unresolved pixel.
    void ComputeShallowPoints(const double* cRe, const double* cIm,
                              int count, float* output, int mxit) const;
    // Verification helper: brute-force high-precision escape time into out[],
    // reusing the grid (_c0/_dx/_dy) set by the most recent Compute() call.
    // step>1 samples only pixels on a (step x step) grid (others left untouched)
    // so the O(mxit) full-precision oracle stays feasible at extreme depth.
    // Requires sub == 1.
    void ComputeDirect(int mxit, float* out, int step = 1, int c_method = 0);
    void ComputeJuliaDirect(mpf_t fixed_c_re, mpf_t fixed_c_im, int mxit,
                            float* out, int step = 1, int c_method = 0);
    void Output(char fname[]) const;
    void setPrecision(int precision);
    HPComp getHighIterationPoint() const;
    HPComp getSymmetryCenter() const;
    HPComp getRandomZoomPoint() const;
    uint8_t* getImage() const;

    void SetHalt(bool flag);
    void SetProgress(std::atomic<float>* progress, float offset = 0.0f, float scale = 1.0f);
    double escapeRadius() const { return _ESCAPE_RADIUS; }
    double escapeRadiusSquared() const {
        return static_cast<double>(_ESCAPE_RADIUS * _ESCAPE_RADIUS);
    }

    // Palette density, used only by the SAC/Feather adaptive-supersampling
    // detector to flag pixels by their actual colour change (colour-index change
    // per pixel ~= stripe_diff * density/20 * colP). The GUI sets it from the live
    // colour density before each frame; headless tools leave the default.
    float _ss_density = 60.0f;
    void setDensity(float d) { _ss_density = d; }
    // Optional second output buffer for NORMAL_MAP: per-pixel surface-normal angle
    // (same layout as the iter buffer). nullptr => no normal output.
    void setNormalBuffer(float* n) { _normal = n; }

private:
    inline bool escape(Comp& z) const;
    inline bool escape(mpf_t& z_re, mpf_t& z_im) const;
    inline float getEscapeTime(Comp& z, int i) const;
    inline float getEscapeTime(mpf_t& z_re, mpf_t& z_im, int i);

    inline int getIndex(std::array<int, 4>& arr) const;
    inline int getIndex(int i, int j, int u, int v) const;

    void setPixel(std::array<int, 4> p, float iteration) const;
    void progressSet(double local);
    void progressBegin(int total, double begin, double span);
    void progressAdvance();
    inline void markDone(int i);
    void stepParallel(std::set<std::array<int, 4>>& s, int mx_ref_it, int mxit,
                      int c_method = 0,
                      const std::vector<std::array<int, 4>>* contiguous = nullptr);
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
    void solveJuliaShallowSimdList(const double* z0re, const double* z0im,
                                   int count, double cre, double cim,
                                   float* out, int mxit, bool ede) const;
    float accurateJuliaPoint(mpf_t z0_re, mpf_t z0_im, mpf_t c_re, mpf_t c_im,
                             int mxit, bool ede) const;

    // Bivariate Linear Approximation (Zhuoran / Fraktaler-3). A BLA skips l
    // reference iterations via dz -> A*dz + B*dc when |dz| < R. Searches touch
    // only the compact validity metadata; A/B coefficients are fetched after a
    // candidate passes, avoiding 32 B of cache traffic per rejected level.
    struct BLAEntry { double r2; int l; };
    struct BLACoeff { double ar, ai, br, bi; };
    // Derivative deltas carried through a skip under EDE: dd -> C*dz + A*dd + E*dc.
    struct BLADeriv { double cr, ci, er, ei; };
    std::vector<std::vector<BLAEntry>> _bla;
    std::vector<std::vector<BLACoeff>> _blaCoeff;
    std::vector<std::vector<BLADeriv>> _blaD;   // parallel to _bla; only built under EDE
    double _bla_eps = 0.0;
    double _bla_rmax2 = 0.0;      // largest BLA validity radius^2 (gate for tryBLA)
    double _bla_dcmax2 = 0.0;     // max |dc|^2 over the frame (BLA-effectiveness gate)
    bool _simd_bla_idle = false;  // measured: BLA barely skips -> compute-bound -> SIMD ok
    bool _simd_measured = false;  // whether _simd_bla_idle has been measured this frame
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
                 double dr, double di, double ESC2, int mx_ref_it,
                 double* outAB = nullptr) const;
    int createRef(std::set<std::array<int, 4>>& s, int pr_it, int mxit, bool random,
                  int c_method = 0, bool view_center = false);
    // Periodic-reference experiment (opt-in via MANDEL_PERIODIC, floatexp non-EDE
    // path). Detect a minibrot nucleus near the view centre; returns its period
    // (0 if none in view) and fills the nucleus coordinate.
    int findNucleus(mpf_t nuc_re, mpf_t nuc_im, int maxp);
    // Build the reference from ONE period at the nucleus (GMP), then replicate it
    // over [0, mxit] so the existing delta loop / BLA are unchanged. Returns the
    // effective reference length (mxit). Only the Z orbit is periodic (no EDE
    // derivative), so this is non-EDE only.
    int createPeriodicRef(int period, int mxit, int c_method);
    // Minibrot linear size estimate (Heiland-Allen / Kalles-Fraktaler) for the
    // just-built period-_ref_period nucleus orbit, read from the floatexp shadow
    // _z*_fe. Returns |size| in c-space; the island is "sub-pixel" (unresolvable,
    // glitch-prone under a bounded reference) when this is below ~one pixel. Pure
    // FloatExp over the existing orbit, so essentially free.
    FloatExp nucleusSizeMag() const;
    bool calCoefficient(int i, int pr_it, int c_method = 0);
    // Deep-zoom rescaled perturbation (Zhuoran z = S w): the delta w stays an
    // O(1) double while the floatexp scale S carries the deep exponent, so the
    // inner loop is native-double yet correct far past double's ~1e320 underflow.
    // Used when _use_floatexp; returns the pixel escape value (interior -> -2).
    float pixelRescaled(FloatExp dcr, FloatExp dci, int mx_ref_it, int mxit, int c_method, float* normalOut = nullptr) const;
    // AVX2 4-wide version of pixelRescaled: iterates up to 4 deep-zoom pixels in
    // lockstep, vectorising the rescaled double step while keeping BLA / rescale /
    // rebase as cheap per-lane scalar events. Writes smooth-iteration values into
    // out[g..g+lanes). Mirrors pixelRescaled's math; no EDE derivative (dz-only).
    void solveRescaledSimd4(const FloatExp* Dcr, const FloatExp* Dci, int g, int lanes,
                            int mx_ref_it, int mxit, int c_method, float* out) const;
    int SACheckMagnitude() const;
    float accuratePointCompute(mpf_t c_re, mpf_t c_im, int mxit, int c_method = 0) const;
    bool refTailEscapes(int last, int extra) const;
    // Builds the deep-zoom reference orbit (periodic-nucleus or full), applies the
    // content-aware floatexp escalation and the mxit-boundary sensitivity rebuild,
    // and returns its length. Extracted from Compute; sets _ref_period/_ref_bounded/
    // _use_floatexp, seeds the pixel set `s`, and accumulates the timer `pf_ref`.
    int buildReferenceOrbit(std::set<std::array<int, 4>>& s, mpf_t scale, int mxit,
                            int c_method, bool profile, double& pf_ref);
    double floatPointCompute(Float c_re, Float c_im, int mxit, int c_method = 0, double* normalOut = nullptr) const;
    bool attractor(double z_in_re, double z_in_im, const double c_re, const double c_im, int period) const;

private:
    volatile bool _flag_halt = false;
    std::atomic<float>* _progress = nullptr;
    std::atomic<int> _progress_done{0};
    int _progress_total = 0, _progress_report_step = 1;
    double _progress_offset = 0.0, _progress_scale = 1.0;
    double _progress_begin = 0.0, _progress_span = 1.0;
    volatile bool _sub_flag;
    float _ESCAPE_RADIUS = 1e8f;
    const double _TOL = 1e12;
    static const int _SA_N = 30;

    int _sub;
    float* _iter;
    float* _normal = nullptr;   // optional NORMAL_MAP output (surface-normal angle)
    const int _w, _h;
    const int _mxit;
    mpf_t* _z_re, * _z_im;  // reference orbit
    Comp* _zf;
    // BigFixed reference orbit (env MANDEL_BIGFIXED): a Mandelbrot-specific fixed-
    // point bignum that replaces the generic mpf_t orbit recurrence (~1.5x faster
    // multiply via high-half short product). Auto-enabled on the deep floatexp path
    // (where it wins); env forces on(>=1)/off(0). Set per reference build in createRef.
    bool _use_bigfixed = false;
    // Rotating 2-buffer like _z_re/_z_im; _bfL limbs (with guard for near-zero).
    BigFixed _bz_re[2], _bz_im[2], _bc_re, _bc_im, _bt1, _bt2, _bab, _bre, _bim;
    std::vector<uint64_t> _bftmp;   // 2*_bfL scratch for bf_mul
    int _bfL = 0;

    // Exterior DE
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
    // Prefix sum of the reference orbit's stripe values (freq-7 SAC), P[m] =
    // sum_{k=1..m} stripe(X_k). Lets pixelRescaled add a BLA skip's omitted
    // stripe contributions back (z ~ X during a valid skip), so Feather stays
    // reference-independent (no pan-dependent radial halo). Built only for the
    // floatexp+BLA+SAC path.
    std::vector<double> _sacRefPre;
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
    // Periodic-reference experiment: when >0 the reference is a period-p minibrot
    // nucleus and the floatexp delta loop indexes it modulo p (only p entries are
    // built, not mxit). 0 = normal (full) reference.
    int _ref_period = 0;
    // Set when a periodic-nucleus reference is used but the nucleus's minibrot body is
    // smaller than a pixel. The BLA skip then damps the near-nucleus ("atom domain")
    // pixels' divergence and paints them as a false-interior body, so BLA is disabled
    // for such references (Compute). Detected from the nucleus-size estimate.
    bool _ref_subpixel = false;
    // A full escaping reference whose orbit passes below double's exponent range.
    // Its double BLA shadow cannot safely represent the deep-zero crossing; use the
    // exact rescaled AVX2 path without BLA instead.
    bool _ref_deep_zero = false;
    // Sub-pixel offset of the nucleus from its nearest pixel (_ref), so the double
    // path's integer-pixel dc can be corrected to dc = pixel_c - nucleus.
    Float _ref_frac_re = 0.0, _ref_frac_im = 0.0;
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