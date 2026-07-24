// Offline validation for analytic AA + adaptive supersampling. All variants are
// rendered on the engine's NATIVE per-pixel sub-sampling grid so samples align:
//   gt.bmp       - engine SS at high sub, flag-all  (ground truth ~ infinite SS)
//   point.bmp    - 1 sample/pixel, no AA            (noisy baseline)
//   aa.bmp       - 1 sample/pixel + analytic AA
//   full3x.bmp   - engine 3x SS on every pixel      (the "3x supersampling" ref)
//   adaptive.bmp - adaptive SS (sub=MANDEL_ADAPT_SUB) + AA base  (the shipped path)
// Prints per-channel mean |err| (0..255) of each vs GT, and % pixels supersampled.
//
// Usage: aa_bench <cx> <cy> <scale> <mxit> [W H GTsub]
//   env: MANDEL_ADAPT_SUB (adaptive sub, default 3), MANDEL_SS_K (flag threshold).

#include <gmp.h>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

#include "mandel_perturbation.h"
#include "interpolate.h"

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

static inline float colorFunction(float it) {
    return powf(logf(it + 2), logf(it + 2) * logf(it + 2) / color_density);
}

static void getColor(float it, float& r, float& g, float& b) {
    if (it < 0) { r = g = b = 0; return; }
    float f = colorFunction(it);
    int x = (int)(f * colP) % colP; if (x < 0) x += colP;
    r = color_map[0][x]; g = color_map[1][x]; b = color_map[2][x];
}

// ---- sRGB linear-light averaging (matches win32_main.cpp) ------------------
static inline double srgb2linF(double v255) {
    double c = v255 / 255.0; if (c < 0) c = 0; else if (c > 1) c = 1;
    return c <= 0.04045 ? c / 12.92 : std::pow((c + 0.055) / 1.055, 2.4);
}
static inline double srgbEncode(double L) {
    if (L <= 0.0) return 0.0; if (L >= 1.0) return 1.0;
    return L <= 0.0031308 ? 12.92 * L : 1.055 * std::pow(L, 1.0 / 2.4) - 0.055;
}
static inline double enc255(double L) { return srgbEncode(L) * 255.0; }

// ---- analytic AA (must mirror win32_main.cpp exactly) ----------------------
static double g_palInt[3][colP + 1];
static double g_palMean[3];
static void prepareColorFilter() {
    for (int c = 0; c < 3; ++c) {
        double s = 0; g_palInt[c][0] = 0;
        for (int i = 0; i < colP; ++i) { s += srgb2linF(color_map[c][i]); g_palInt[c][i + 1] = s; }
        g_palMean[c] = enc255(s / colP);
    }
}
static double palPrefix(int c, double x) {
    double q = std::floor(x / colP);
    double r = x - q * colP;
    int i = (int)r; if (i >= colP) i = colP - 1;
    double t = r - i;
    double part = g_palInt[c][i] * (1.0 - t) + g_palInt[c][i + 1] * t;
    return q * g_palInt[c][colP] + part;
}
static void getColorAA(float v, float vL, float vR, float vU, float vD, float& r, float& g, float& b) {
    if (v < 0) { r = g = b = 0; return; }
    auto ext = [](float x) { return x != EMPTYPIXEL && x >= 0.0f; };
    float f = colorFunction(v);
    auto cf = [&](float nv) { return colorFunction(nv); };
    float fx = 0, fy = 0;
    if (ext(vL) && ext(vR)) fx = 0.5f * (cf(vR) - cf(vL));
    else if (ext(vR)) fx = cf(vR) - f; else if (ext(vL)) fx = f - cf(vL);
    if (ext(vU) && ext(vD)) fy = 0.5f * (cf(vD) - cf(vU));
    else if (ext(vD)) fy = cf(vD) - f; else if (ext(vU)) fy = f - cf(vU);
    float gradf = std::sqrt(fx * fx + fy * fy);
    double width = (double)gradf * colP;
    double centerU = (double)f * colP;
    if (width < 1.0) {
        int x = (int)centerU % colP; if (x < 0) x += colP;
        r = color_map[0][x]; g = color_map[1][x]; b = color_map[2][x]; return;
    }
    if (width >= colP) { r = g_palMean[0]; g = g_palMean[1]; b = g_palMean[2]; return; }
    double a = centerU - 0.5 * width, e = centerU + 0.5 * width;
    auto avg = [&](int c) {
        double m = (palPrefix(c, e) - palPrefix(c, a)) / width;
        if (m < 0) m = 0;
        double s = enc255(m);
        return s > 255.0 ? 255.0f : (float)s;
    };
    r = avg(0); g = avg(1); b = avg(2);
}

static void writeBMP(const char* path, const std::vector<uint8_t>& rgb, int W, int H) {
    int row = (W * 3 + 3) & ~3, dataSize = row * H;
    uint8_t fh[14] = { 'B','M' }; uint32_t fsize = 14 + 40 + dataSize;
    memcpy(fh + 2, &fsize, 4); uint32_t off = 54; memcpy(fh + 10, &off, 4);
    uint8_t ih[40] = { 0 }; uint32_t v;
    v = 40; memcpy(ih, &v, 4); v = W; memcpy(ih + 4, &v, 4); v = H; memcpy(ih + 8, &v, 4);
    uint16_t planes = 1, bpp = 24; memcpy(ih + 12, &planes, 2); memcpy(ih + 14, &bpp, 2);
    v = dataSize; memcpy(ih + 20, &v, 4);
    FILE* f = fopen(path, "wb"); fwrite(fh, 1, 14, f); fwrite(ih, 1, 40, f);
    std::vector<uint8_t> line(row, 0);
    for (int y = H - 1; y >= 0; --y) {
        for (int x = 0; x < W; ++x) {
            const uint8_t* p = &rgb[(y * W + x) * 3];
            line[x * 3] = p[2]; line[x * 3 + 1] = p[1]; line[x * 3 + 2] = p[0];
        }
        fwrite(line.data(), 1, row, f);
    }
    fclose(f);
}

// Render an engine supersampled grid at the given sub with a flag threshold K
// (passed via MANDEL_SS_K). Caller owns the returned buffer.
static float* renderSS(int W, int H, int mxit, int sub, int precision,
                       mpf_t mcx, mpf_t mcy, mpf_t msc, const char* kval) {
    _putenv_s("MANDEL_SS_K", kval);
    size_t N = (size_t)W * H * sub * sub;
    float* it = new float[N];
    for (size_t i = 0; i < N; ++i) it[i] = EMPTYPIXEL;
    Mandel m(W, H, mxit, sub, it);
    m.setPrecision(precision);
    m.Compute(mcx, mcy, msc, mxit, ColoringMethod::SUPER_SAMPLING);
    return it;
}

// Colour one pixel from an engine SS grid: RMS-average the sub-block for a
// flagged pixel, else single centre sample. Returns true if the pixel was
// supersampled (flagged).
static bool colorBlock(const float* it, int W, int sub, int i, int j,
                       float& r, float& g, float& b) {
    int c = sub / 2, stride = W * sub;
    int idx = (i * sub + c) * stride + (j * sub + c);
    if (it[idx] == EMPTYPIXEL) { r = g = b = 0; return false; }
    if (it[idx - (stride + 1) * c] != EMPTYPIXEL) {          // flagged
        double rl = 0, gl = 0, bl = 0; int n = 0;
        for (int a = -c; a <= c; ++a)
            for (int bb = -c; bb <= c; ++bb) {
                float rr, gg, bbb; getColor(it[idx + a * stride + bb], rr, gg, bbb);
                rl += srgb2linF(rr); gl += srgb2linF(gg); bl += srgb2linF(bbb); n++;
            }
        r = (float)enc255(rl / n); g = (float)enc255(gl / n); b = (float)enc255(bl / n);
        return true;
    }
    getColor(it[idx], r, g, b);                              // uniform: centre
    return false;
}

// RMS (gamma-2) variant of colorBlock, i.e. the OLD averaging, for A/B comparison.
static void colorBlockRMS(const float* it, int W, int sub, int i, int j,
                          float& r, float& g, float& b) {
    int c = sub / 2, stride = W * sub;
    int idx = (i * sub + c) * stride + (j * sub + c);
    if (it[idx] == EMPTYPIXEL) { r = g = b = 0; return; }
    if (it[idx - (stride + 1) * c] != EMPTYPIXEL) {
        double r2 = 0, g2 = 0, b2 = 0; int n = 0;
        for (int a = -c; a <= c; ++a)
            for (int bb = -c; bb <= c; ++bb) {
                float rr, gg, bbb; getColor(it[idx + a * stride + bb], rr, gg, bbb);
                r2 += rr * rr; g2 += gg * gg; b2 += bbb * bbb; n++;
            }
        r = (float)std::sqrt(r2 / n); g = (float)std::sqrt(g2 / n); b = (float)std::sqrt(b2 / n);
        return;
    }
    getColor(it[idx], r, g, b);
}

// Like colorBlock, but each sub-sample is itself analytic-AA filtered using its
// neighbouring sub-samples (spacing = one sub-cell). This resolves colour cycling
// analytically per sub-cell, so a low sub reaches a much higher effective quality
// for the same number of fractal evaluations. Returns true if flagged.
static bool colorBlockAA(const float* it, int W, int H, int sub, int i, int j,
                         float& r, float& g, float& b) {
    int c = sub / 2, stride = W * sub, rowmax = H * sub, colmax = W * sub;
    int idx = (i * sub + c) * stride + (j * sub + c);
    if (it[idx] == EMPTYPIXEL) { r = g = b = 0; return false; }
    if (it[idx - (stride + 1) * c] == EMPTYPIXEL) { getColor(it[idx], r, g, b); return false; }
    double r2 = 0, g2 = 0, b2 = 0; int n = 0;
    for (int a = -c; a <= c; ++a)
        for (int bb = -c; bb <= c; ++bb) {
            int si = i * sub + c + a, sj = j * sub + c + bb;
            int sidx = si * stride + sj;
            float v  = it[sidx];
            float vL = sj > 0          ? it[sidx - 1]      : EMPTYPIXEL;
            float vR = sj < colmax - 1 ? it[sidx + 1]      : EMPTYPIXEL;
            float vU = si > 0          ? it[sidx - stride] : EMPTYPIXEL;
            float vD = si < rowmax - 1 ? it[sidx + stride] : EMPTYPIXEL;
            float rr, gg, bbb; getColorAA(v, vL, vR, vU, vD, rr, gg, bbb);
            r2 += srgb2linF(rr); g2 += srgb2linF(gg); b2 += srgb2linF(bbb); n++;
        }
    r = (float)enc255(r2 / n); g = (float)enc255(g2 / n); b = (float)enc255(b2 / n);
    return true;
}

// Sample an aligned s x s grid for pixel (i,j) from the GT super-grid; return the
// RMS colour and the max-channel colour standard deviation across the s*s samples.
static void sampleGrid(const float* it, int W, int GTsub, int s, int i, int j,
                       float& r, float& g, float& b, double& sd) {
    int stride = W * GTsub;
    double rl=0,gl=0,bl=0, s2=0,sm=0; int n=0;
    for (int k = 0; k < s; ++k)
        for (int l = 0; l < s; ++l) {
            int rr = (int)((k + 0.5) * GTsub / s), cc = (int)((l + 0.5) * GTsub / s);
            int sidx = (i * GTsub + rr) * stride + (j * GTsub + cc);
            float R,G,B; getColor(it[sidx], R,G,B);
            rl+=srgb2linF(R); gl+=srgb2linF(G); bl+=srgb2linF(B);
            double lum = 0.299*R + 0.587*G + 0.114*B; s2 += lum*lum; sm += lum; n++;
        }
    r=(float)enc255(rl/n); g=(float)enc255(gl/n); b=(float)enc255(bl/n);
    double vv = s2/n - (sm/n)*(sm/n); sd = std::sqrt(vv > 0 ? vv : 0);
}

int main(int argc, char** argv) {
    if (argc < 5) { fprintf(stderr, "usage: aa_bench cx cy scale mxit [W H GTsub]\n"); return 1; }
    const char* cx = argv[1]; const char* cy = argv[2];
    const char* scaleStr = argv[3]; int mxit = atoi(argv[4]);
    int W = argc > 5 ? atoi(argv[5]) : 480;
    int H = argc > 6 ? atoi(argv[6]) : 270;
    int GTsub = argc > 7 ? atoi(argv[7]) : 15;
    if (GTsub % 2 == 0) GTsub++;
    int subA = 3; { const char* e = getenv("MANDEL_ADAPT_SUB"); if (e) subA = atoi(e); }
    if (subA % 2 == 0) subA++;
    const char* Kstr = getenv("MANDEL_SS_K"); if (!Kstr) Kstr = "8";

    double sd = atof(scaleStr);
    int scaleExp = (sd > 1.0) ? (int)(log10(sd) + 0.5) : 0;
    int precision = (int)((scaleExp + 20) * log(10) / log(2)) + 64;
    mpf_set_default_prec(precision);
    colorMapInit(); prepareColorFilter();

    mpf_t mcx, mcy, msc;
    mpf_init_set_str(mcx, cx, 10); mpf_init_set_str(mcy, cy, 10); mpf_init_set_str(msc, scaleStr, 10);
    printf("Location scale=%s mxit=%d prec=%d  W=%d H=%d  GTsub=%d  adaptSub=%d  K=%s\n",
           scaleStr, mxit, precision, W, H, GTsub, subA, Kstr);

    // Ground truth: engine SS at high sub, flag everything (K huge).
    float* gtIt = renderSS(W, H, mxit, GTsub, precision, mcx, mcy, msc, "100000");
    // "Full 3x SS": engine 3x SS on every pixel (the user's reference).
    float* f3It = renderSS(W, H, mxit, 3, precision, mcx, mcy, msc, "100000");
    // Adaptive: chosen sub, chosen K.
    float* adIt = renderSS(W, H, mxit, subA, precision, mcx, mcy, msc, Kstr);

    std::vector<uint8_t> imgG(W*H*3), imgP(W*H*3), imgA(W*H*3), imgF(W*H*3), imgD(W*H*3), imgE(W*H*3);
    std::vector<float> gtF(W*H*3);
    std::vector<char> gtFlag(W*H, 0);
    double sP=0,sA=0,sF=0,sD=0,sE=0, mP=0,mA=0,mF=0,mD=0,mE=0; long cnt=0, ssPix=0, gtSS=0;
    double sRms=0, mRms=0;   // OLD RMS/gamma-2 GT vs correct sRGB-linear GT

    const int cG = GTsub/2, strideG = W*GTsub;
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            // --- ground truth
            float Gr,Gg,Gb; bool gflag = colorBlock(gtIt, W, GTsub, i, j, Gr,Gg,Gb); if (gflag) gtSS++;
            gtFlag[i*W+j] = gflag ? 1 : 0;
            float* gf = &gtF[(i*W+j)*3]; gf[0]=Gr; gf[1]=Gg; gf[2]=Gb;
            uint8_t* pG=&imgG[(i*W+j)*3]; pG[0]=(uint8_t)(Gr+.5f);pG[1]=(uint8_t)(Gg+.5f);pG[2]=(uint8_t)(Gb+.5f);
            // OLD RMS/gamma-2 downsample of the same high-SS grid, for A/B accuracy.
            float Rr,Rg,Rb; colorBlockRMS(gtIt, W, GTsub, i, j, Rr,Rg,Rb);
            sRms += (std::fabs(Rr-Gr)+std::fabs(Rg-Gg)+std::fabs(Rb-Gb))/3.0;
            if (std::fabs(Rr-Gr)>mRms) mRms=std::fabs(Rr-Gr);

            // --- point + AA, from the GT grid centre subpixel (same centring)
            int idx = (i*GTsub+cG)*strideG + (j*GTsub+cG);
            float vC = gtIt[idx];
            float vL = j>0     ? gtIt[idx-GTsub]          : EMPTYPIXEL;
            float vR = j<W-1   ? gtIt[idx+GTsub]          : EMPTYPIXEL;
            float vU = i>0     ? gtIt[idx-strideG*GTsub]  : EMPTYPIXEL;
            float vD = i<H-1   ? gtIt[idx+strideG*GTsub]  : EMPTYPIXEL;
            float pr,pg,pb; getColor(vC,pr,pg,pb);
            float ar,ag,ab; getColorAA(vC,vL,vR,vU,vD,ar,ag,ab);
            uint8_t* pP=&imgP[(i*W+j)*3]; pP[0]=(uint8_t)(pr+.5f);pP[1]=(uint8_t)(pg+.5f);pP[2]=(uint8_t)(pb+.5f);
            uint8_t* pA=&imgA[(i*W+j)*3]; pA[0]=(uint8_t)(ar+.5f);pA[1]=(uint8_t)(ag+.5f);pA[2]=(uint8_t)(ab+.5f);

            // --- full 3x SS
            float Fr,Fg,Fb; colorBlock(f3It, W, 3, i, j, Fr,Fg,Fb);
            uint8_t* pF=&imgF[(i*W+j)*3]; pF[0]=(uint8_t)(Fr+.5f);pF[1]=(uint8_t)(Fg+.5f);pF[2]=(uint8_t)(Fb+.5f);

            // --- adaptive: block RMS if flagged, else analytic AA base
            float Dr,Dg,Db;
            {
                int cA=subA/2, sA_stride=W*subA;
                int ia=(i*subA+cA)*sA_stride+(j*subA+cA);
                if (adIt[ia]==EMPTYPIXEL) { Dr=Dg=Db=0; }
                else if (adIt[ia-(sA_stride+1)*cA]!=EMPTYPIXEL) { colorBlock(adIt,W,subA,i,j,Dr,Dg,Db); ssPix++; }
                else {
                    float aL=j>0?adIt[ia-subA]:EMPTYPIXEL, aR=j<W-1?adIt[ia+subA]:EMPTYPIXEL;
                    float aU=i>0?adIt[ia-sA_stride*subA]:EMPTYPIXEL, aD=i<H-1?adIt[ia+sA_stride*subA]:EMPTYPIXEL;
                    getColorAA(adIt[ia],aL,aR,aU,aD,Dr,Dg,Db);
                }
            }
            uint8_t* pD=&imgD[(i*W+j)*3]; pD[0]=(uint8_t)(Dr+.5f);pD[1]=(uint8_t)(Dg+.5f);pD[2]=(uint8_t)(Db+.5f);

            // --- adaptive with per-subsample AA (same fractal evals, better colour)
            float Er,Eg,Eb;
            {
                int cA=subA/2, sA_stride=W*subA;
                int ia=(i*subA+cA)*sA_stride+(j*subA+cA);
                if (adIt[ia]==EMPTYPIXEL) { Er=Eg=Eb=0; }
                else if (!colorBlockAA(adIt,W,H,subA,i,j,Er,Eg,Eb)) {   // unflagged -> pixel AA
                    float aL=j>0?adIt[ia-subA]:EMPTYPIXEL, aR=j<W-1?adIt[ia+subA]:EMPTYPIXEL;
                    float aU=i>0?adIt[ia-sA_stride*subA]:EMPTYPIXEL, aD=i<H-1?adIt[ia+sA_stride*subA]:EMPTYPIXEL;
                    getColorAA(adIt[ia],aL,aR,aU,aD,Er,Eg,Eb);
                }
            }
            uint8_t* pE=&imgE[(i*W+j)*3]; pE[0]=(uint8_t)(Er+.5f);pE[1]=(uint8_t)(Eg+.5f);pE[2]=(uint8_t)(Eb+.5f);

            auto err=[&](float r,float g,float b){ return (std::fabs(r-Gr)+std::fabs(g-Gg)+std::fabs(b-Gb))/3.0; };
            double eP=err(pr,pg,pb), eA=err(ar,ag,ab), eF=err(Fr,Fg,Fb), eD=err(Dr,Dg,Db), eE=err(Er,Eg,Eb);
            sP+=eP; sA+=eA; sF+=eF; sD+=eD; sE+=eE;
            if(eP>mP)mP=eP; if(eA>mA)mA=eA; if(eF>mF)mF=eF; if(eD>mD)mD=eD; if(eE>mE)mE=eE; cnt++;
        }

    writeBMP("gt.bmp", imgG, W, H);
    writeBMP("point.bmp", imgP, W, H);
    writeBMP("aa.bmp", imgA, W, H);
    writeBMP("full3x.bmp", imgF, W, H);
    writeBMP("adaptive.bmp", imgD, W, H);
    writeBMP("adaptive_saa.bmp", imgE, W, H);

    printf("GT boundary pixels supersampled: %.1f%%\n", 100.0*gtSS/cnt);
    printf("OLD RMS/gamma-2 avg vs correct sRGB-linear avg (both %dx SS): mean=%.2f max=%.1f  <- color-accuracy gain\n",
           GTsub, sRms/cnt, mRms);
    printf("mean |err| vs GT (0..255):\n");
    printf("   point                 = %6.2f   (max %.1f)\n", sP/cnt, mP);
    printf("   analytic AA           = %6.2f   (max %.1f)\n", sA/cnt, mA);
    printf("   full 3x SS            = %6.2f   (max %.1f)\n", sF/cnt, mF);
    printf("   adaptive+AA(s%d)       = %6.2f   (max %.1f)   supersampled=%.1f%%\n",
           subA, sD/cnt, mD, 100.0*ssPix/cnt);
    printf("   adaptive+per-subAA(s%d)= %6.2f   (max %.1f)   [same %d evals/flagged px]\n",
           subA, sE/cnt, mE, subA*subA);

    // --- progressive within-pixel refinement: flagged pixels start at 3x3; escalate
    // to 5x5 only when the 3x3 colour std-dev exceeds T (nested design reuses the 9).
    // Simulated on the aligned GT grid so all sampling is apples-to-apples.
    {
        const int NT = 5; const double Th[NT] = { 6, 10, 16, 24, 1e9 };
        double prErr[NT] = {0}, prEv[NT] = {0};
        double u3=0,u5=0, u3ev=0,u5ev=0;
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) {
                float* gf = &gtF[(i*W+j)*3];
                auto er=[&](float r,float g,float b){ return (std::fabs(r-gf[0])+std::fabs(g-gf[1])+std::fabs(b-gf[2]))/3.0; };
                if (!gtFlag[i*W+j]) {                       // uniform: 1 eval everywhere
                    float r,g,b; double sd0; sampleGrid(gtIt,W,GTsub,1,i,j,r,g,b,sd0);
                    double e=er(r,g,b);
                    u3+=e; u5+=e; u3ev+=1; u5ev+=1;
                    for (int t=0;t<NT;++t){ prErr[t]+=e; prEv[t]+=1; }
                    continue;
                }
                float r3,g3,b3; double sd3; sampleGrid(gtIt,W,GTsub,3,i,j,r3,g3,b3,sd3);
                float r5,g5,b5; double sd5; sampleGrid(gtIt,W,GTsub,5,i,j,r5,g5,b5,sd5);
                u3+=er(r3,g3,b3); u3ev+=9; u5+=er(r5,g5,b5); u5ev+=25;
                for (int t=0;t<NT;++t) {
                    if (sd3 <= Th[t]) { prErr[t]+=er(r3,g3,b3); prEv[t]+=9; }
                    else              { prErr[t]+=er(r5,g5,b5); prEv[t]+=25; }
                }
            }
        printf("\ncompute/quality frontier (avg fractal evals per pixel vs err):\n");
        printf("   uniform sub=3     : %.1f evals/px   err=%.2f\n", u3ev/cnt, u3/cnt);
        printf("   uniform sub=5     : %.1f evals/px   err=%.2f\n", u5ev/cnt, u5/cnt);
        for (int t=0;t<NT;++t)
            printf("   progressive T=%-4.0f: %.1f evals/px   err=%.2f\n", Th[t], prEv[t]/cnt, prErr[t]/cnt);
    }

    printf("Wrote gt/point/aa/full3x/adaptive/adaptive_saa .bmp\n");
    delete[] gtIt; delete[] f3It; delete[] adIt;
    return 0;
}
