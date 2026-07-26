// Diff-image harness: for each test case, render the perturbation engine three
// ways and compare against a full-resolution brute-force ground-truth oracle:
//   * engine, BLA OFF (EDE)   -> diff vs GT
//   * engine, BLA ON  (EDE)   -> diff vs GT
// Emits, per case, into build/diffs/:
//   <case>_gt.bmp        ground-truth (brute-force full-precision) render
//   <case>_blaoff.bmp    engine render, BLA off
//   <case>_blaon.bmp     engine render, BLA on
//   <case>_diff_off.bmp  |blaoff - gt| colour diff (amplified), class-flip in red
//   <case>_diff_on.bmp   |blaon  - gt| colour diff (amplified), class-flip in red
//
// The oracle (the expensive part) is computed ONCE per case; the two engine
// passes toggle BLA via _putenv_s between Compute() calls. EDE is on by default.
//
// Usage: diff_img [case] [W H] [mxit]
//   case = shallow | deep51 | ticktock | flake | all   (default: all)

#include <gmp.h>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>

#include "mandel_perturbation.h"
#include "interpolate.h"
#include "test_cases.h"

static const int colP = 2048;
static float color_map[3][colP];
static const float color_density = 60.0f;

static void colorMapInit() {
    static const int N = 6;
    static const float pos[N] = { 0.0f, 0.16f, 0.42f, 0.6425f, 0.8575f, 1.0f };
    static const float col[3][N] = { {0,32,237,255,0,0}, {70,107,255,170,2,70}, {100,203,255,0,0,100} };
    for (int i = 0; i < 3; ++i)
        mono_cubic_interpolate(pos, col[i], N, color_map[i], colP);
}

static int g_method = 0;   // set per run; selects the palette mapping

// Mirrors win32_main.cpp::colorFunction: EDE distance estimates are tone-mapped
// with tanh, iteration counts with the log-power curve.
static inline float color_func(float it) {
    if (g_method & ColoringMethod::EXTERIOR_DIST_EST)
        return tanhf(it / color_density * 5.0f);
    float l = logf(it + 2.0f);
    return powf(l, l * l / color_density);
}

static void getColor(float it, uint8_t& r, uint8_t& g, uint8_t& b) {
    if (it < 0) { r = g = b = 0; return; }   // interior -> black
    float f = color_func(it);
    int x = (int)(f * colP) % colP;
    if (x < 0) x += colP;
    r = (uint8_t)(color_map[0][x] + 0.5f);
    g = (uint8_t)(color_map[1][x] + 0.5f);
    b = (uint8_t)(color_map[2][x] + 0.5f);
}

static void writeBMP(const char* path, const std::vector<uint8_t>& rgb, int W, int H) {
    int row = (W * 3 + 3) & ~3;
    int dataSize = row * H;
    uint8_t fh[14] = { 'B','M' };
    uint32_t fsize = 14 + 40 + dataSize;
    memcpy(fh + 2, &fsize, 4);
    uint32_t off = 54; memcpy(fh + 10, &off, 4);
    uint8_t ih[40] = { 0 }; uint32_t v;
    v = 40; memcpy(ih + 0, &v, 4);
    v = W; memcpy(ih + 4, &v, 4);
    v = H; memcpy(ih + 8, &v, 4);
    uint16_t planes = 1, bpp = 24; memcpy(ih + 12, &planes, 2); memcpy(ih + 14, &bpp, 2);
    v = dataSize; memcpy(ih + 20, &v, 4);
    FILE* f = fopen(path, "wb");
    fwrite(fh, 1, 14, f); fwrite(ih, 1, 40, f);
    std::vector<uint8_t> line(row, 0);
    for (int y = H - 1; y >= 0; --y) {
        for (int x = 0; x < W; ++x) {
            const uint8_t* p = &rgb[(y * W + x) * 3];
            line[x * 3 + 0] = p[2]; line[x * 3 + 1] = p[1]; line[x * 3 + 2] = p[0];
        }
        fwrite(line.data(), 1, row, f);
    }
    fclose(f);
}

static void colorImage(const float* buf, int n, std::vector<uint8_t>& img) {
    img.resize(n * 3);
    for (int i = 0; i < n; ++i)
        getColor(buf[i], img[i * 3], img[i * 3 + 1], img[i * 3 + 2]);
}

struct DiffStats { long class_mismatch; double max_cdiff; double mean_cdiff; long ext_both; };

// Amplified colour-diff heatmap. Class flips (interior<->exterior) painted pure
// red so they can't be missed; other pixels show |ΔRGB| amplified on black.
static DiffStats diffImage(const float* pert, const float* gt,
                           const std::vector<uint8_t>& pc, const std::vector<uint8_t>& gc,
                           int n, int amp, std::vector<uint8_t>& out) {
    out.assign(n * 3, 0);
    DiffStats s{ 0, 0, 0, 0 };
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        bool pi = pert[i] < 0, gi = gt[i] < 0;
        if (pi != gi) {
            ++s.class_mismatch;
            out[i * 3] = 255; out[i * 3 + 1] = 0; out[i * 3 + 2] = 0;   // red
            continue;
        }
        int dr = std::abs((int)pc[i * 3] - (int)gc[i * 3]);
        int dg = std::abs((int)pc[i * 3 + 1] - (int)gc[i * 3 + 1]);
        int db = std::abs((int)pc[i * 3 + 2] - (int)gc[i * 3 + 2]);
        double cd = (dr + dg + db) / 3.0;
        if (!pi) { s.mean_cdiff += cd; sum = s.mean_cdiff; ++s.ext_both; }
        if (cd > s.max_cdiff) s.max_cdiff = cd;
        int a = dr * amp; if (a > 255) a = 255;
        int bb = dg * amp; if (bb > 255) bb = 255;
        int c = db * amp; if (c > 255) c = 255;
        out[i * 3] = (uint8_t)a; out[i * 3 + 1] = (uint8_t)bb; out[i * 3 + 2] = (uint8_t)c;
    }
    s.mean_cdiff = s.ext_both ? sum / s.ext_both : 0.0;
    return s;
}

struct Case { const char* name; std::string cx, cy, scale; int mxit, W, H; };

static std::string pow10(int n) { std::string s = "1"; s.append(n, '0'); return s; }

int main(int argc, char** argv) {
    std::string which = (argc > 1) ? argv[1] : "all";
    int Warg = (argc > 3) ? atoi(argv[2]) : 0;
    int Harg = (argc > 3) ? atoi(argv[3]) : 0;
    int mxarg = (argc > 4) ? atoi(argv[4]) : 0;

    const int amp = 8;   // colour-diff amplification for visibility
    int c_method = (getenv("MANDEL_NOEDE") != nullptr) ? 0 : ColoringMethod::EXTERIOR_DIST_EST;
    g_method = c_method;

    std::string deep51_scale = "3831277"; deep51_scale.append(45, '0');   // 3.831277e51
    std::vector<Case> cases = {
        { "shallow",  "-0.5", "0",                 pow10(0),      5000,   256, 192 },
        { "deep51",   testcases::deep51_x, testcases::deep51_y, deep51_scale, 2000000, 128, 96 },
        { "ticktock", testcases::ticktock_x, testcases::ticktock_y, pow10(141), 200000, 160, 120 },
        { "flake",    testcases::flake_x, testcases::flake_y,   pow10(157),   200000, 160, 120 },
    };

    colorMapInit();
    system("if not exist diffs mkdir diffs");

    // Custom case from env: MANDEL_CX / MANDEL_CY / MANDEL_SCALE (e.g. 3.131699e286).
    // Expands e-notation/decimal scale to an integer digit string like render.cpp.
    if (const char* cxs = getenv("MANDEL_CX")) {
        std::string sa = getenv("MANDEL_SCALE") ? getenv("MANDEL_SCALE") : "1", scl;
        if (sa.find_first_of("eE.") != std::string::npos) {
            size_t ep = sa.find_first_of("eE");
            std::string mant = (ep == std::string::npos) ? sa : sa.substr(0, ep);
            int e10 = (ep == std::string::npos) ? 0 : atoi(sa.substr(ep + 1).c_str());
            std::string sign;
            if (!mant.empty() && (mant[0] == '-' || mant[0] == '+')) { sign = (mant[0] == '-') ? "-" : ""; mant.erase(0, 1); }
            size_t dp = mant.find('.');
            int frac = (dp == std::string::npos) ? 0 : (int)(mant.size() - dp - 1);
            if (dp != std::string::npos) mant.erase(dp, 1);
            int zeros = e10 - frac; scl = sign + mant; if (zeros > 0) scl.append(zeros, '0');
        } else { int se = atoi(sa.c_str()); scl = "1"; for (int i = 0; i < se; ++i) scl += "0"; }
        cases.clear();
        cases.push_back({ "custom", cxs, getenv("MANDEL_CY"), scl,
                          mxarg ? mxarg : 100000, Warg ? Warg : 160, Harg ? Harg : 120 });
        which = "custom";
    }

    printf("mode: %s\n", c_method ? "EDE" : "raw-iter");
    printf("%-10s %6s %8s  %-18s %-18s\n", "case", "WxH", "mxit",
           "BLA-off vs GT", "BLA-on vs GT");
    printf("%-10s %6s %8s  %-18s %-18s\n", "----", "---", "----",
           "cls / maxΔ / meanΔ", "cls / maxΔ / meanΔ");

    for (auto& tc : cases) {
        if (which != "all" && which != tc.name) continue;
        if (Warg) { tc.W = Warg; tc.H = Harg; }
        if (mxarg) tc.mxit = mxarg;

        int precision = (int)(tc.scale.size() * log(10) / log(2)) + 64;
        mpf_set_default_prec(precision);
        mpf_t cre, cim, scale;
        mpf_init_set_str(cre, tc.cx.c_str(), 10);
        mpf_init_set_str(cim, tc.cy.c_str(), 10);
        mpf_init_set_str(scale, tc.scale.c_str(), 10);

        int n = tc.W * tc.H;
        std::vector<float> itp(n), itd(n), off(n), on(n);

        Mandel mandel(tc.W, tc.H, tc.mxit, 1, itp.data());
        mandel.setPrecision(precision);

        const bool skipOracle = (getenv("MANDEL_NOORACLE") != nullptr);

        // Pass 1 + oracle + pass 2. Order can be reversed (MANDEL_DIFF_REVERSE)
        // to prove the per-BLA result is not an artefact of shared instance state.
        const char* rev = getenv("MANDEL_DIFF_REVERSE");
        const char* first = (rev && atoi(rev)) ? "1" : "0";
        const char* second = (rev && atoi(rev)) ? "0" : "1";

        _putenv_s("MANDEL_BLA", first);
        mandel.Compute(cre, cim, scale, tc.mxit, c_method);
        std::vector<float> firstBuf = itp;

        // Ground truth: full-resolution brute-force oracle (reuses grid).
        for (int i = 0; i < n; ++i) itd[i] = EMPTYPIXEL;
        if (!skipOracle) mandel.ComputeDirect(tc.mxit, itd.data(), 1, c_method);

        _putenv_s("MANDEL_BLA", second);
        mandel.Compute(cre, cim, scale, tc.mxit, c_method);
        std::vector<float> secondBuf = itp;
        _putenv_s("MANDEL_BLA", "0");

        // Map first/second back to off/on regardless of run order.
        if (rev && atoi(rev)) { on = firstBuf; off = secondBuf; }
        else                  { off = firstBuf; on = secondBuf; }

        // Report the mean exterior value (raw DE in EDE mode) so a bailout change
        // can be checked for a colour-density shift.
        { double s = 0; long c = 0; for (int i = 0; i < n; ++i) if (on[i] >= 0) { s += on[i]; ++c; }
          printf("  [%s] mean exterior value (BLA-on) = %.6g over %ld px\n",
                 tc.name, c ? s / c : 0.0, c); fflush(stdout); }

        if (skipOracle) {
            std::vector<uint8_t> onc2; colorImage(on.data(), n, onc2);
            std::string d2 = "diffs\\";
            std::string t2 = std::string(tc.name) + (c_method ? "_ede" : "_iter");
            writeBMP((d2 + t2 + "_blaon.bmp").c_str(), onc2, tc.W, tc.H);
            mpf_clear(cre); mpf_clear(cim); mpf_clear(scale);
            continue;
        }

        std::vector<uint8_t> gc, offc, onc, doff, don;
        colorImage(itd.data(), n, gc);
        colorImage(off.data(), n, offc);
        colorImage(on.data(), n, onc);
        DiffStats soff = diffImage(off.data(), itd.data(), offc, gc, n, amp, doff);
        DiffStats son  = diffImage(on.data(),  itd.data(), onc,  gc, n, amp, don);

        std::string dir = "diffs\\";
        std::string tag = std::string(tc.name) + (c_method ? "_ede" : "_iter");
        writeBMP((dir + tag + "_gt.bmp").c_str(),       gc,   tc.W, tc.H);
        writeBMP((dir + tag + "_blaoff.bmp").c_str(),   offc, tc.W, tc.H);
        writeBMP((dir + tag + "_blaon.bmp").c_str(),    onc,  tc.W, tc.H);
        writeBMP((dir + tag + "_diff_off.bmp").c_str(), doff, tc.W, tc.H);
        writeBMP((dir + tag + "_diff_on.bmp").c_str(),  don,  tc.W, tc.H);

        printf("%-10s %3dx%-3d %8d  %3ld /%5.1f /%5.2f   %3ld /%5.1f /%5.2f\n",
               tc.name, tc.W, tc.H, tc.mxit,
               soff.class_mismatch, soff.max_cdiff, soff.mean_cdiff,
               son.class_mismatch,  son.max_cdiff,  son.mean_cdiff);
        fflush(stdout);

        mpf_clear(cre); mpf_clear(cim); mpf_clear(scale);
    }
    printf("(diff amplification x%d; class flips shown red)\n", amp);
    return 0;
}
