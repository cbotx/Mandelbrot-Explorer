// Coloring-method demo (brute-force double, shallow views only) to compare
// visual styles: smooth-iteration vs Stripe Average Coloring (SAC). Reuses the
// GUI's LAB palette rebuild. Not for deep zoom -- just to eyeball the styles.
//
// Usage: coloring_demo out.bmp W H cx cy scale mxit palIdx mode [SS]
//   mode: 0 = smooth iteration, 1 = stripe average coloring

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include "interpolate.h"
#include "color.h"

static const int CP = 2048;
static float color_map[3][CP];
struct Stop { float pos, r, g, b; };
static std::vector<std::vector<Stop>> pals = {
    { {0,0,70,100},{0.16f,32,107,203},{0.42f,237,255,255},{0.6425f,255,170,0},{0.8575f,0,2,0} },                                   // sunrise
    { {0.05f,30,90,200},{0.17f,242,248,255},{0.30f,25,155,95},{0.42f,10,12,20},{0.55f,125,70,185},{0.67f,250,244,252},{0.80f,240,180,55},{0.92f,18,12,8} }, // gems
    { {0,4,12,38},{0.22f,10,64,120},{0.46f,36,150,190},{0.66f,150,220,230},{0.85f,244,252,255} },                                  // ocean
};
static void rebuild(const std::vector<Stop>& st) {
    int n = (int)st.size();
    std::vector<float> xs(n+2), ys(n+2), out(CP);
    std::vector<std::array<float,3>> lab(n);
    for (int i=0;i<n;++i) rgb2lab(st[i].r, st[i].g, st[i].b, lab[i][0], lab[i][1], lab[i][2]);
    xs[0]=st.back().pos-1; xs[n+1]=st.front().pos+1;
    for (int i=0;i<n;++i) xs[i+1]=st[i].pos;
    for (int c=0;c<3;++c){ ys[0]=lab[n-1][c]; ys[n+1]=lab[0][c]; for(int i=0;i<n;++i) ys[i+1]=lab[i][c];
        mono_cubic_interpolate(xs.data(), ys.data(), n+2, out.data(), CP); for(int i=0;i<CP;++i) color_map[c][i]=out[i]; }
    for (int i=0;i<CP;++i){ float r,g,b; lab2rgb(color_map[0][i],color_map[1][i],color_map[2][i],r,g,b); color_map[0][i]=r;color_map[1][i]=g;color_map[2][i]=b; }
}
static void palColor(float t, float& r, float& g, float& b) {   // t in [0,1)
    int x = (int)(t * CP) % CP; if (x < 0) x += CP;
    r = color_map[0][x]; g = color_map[1][x]; b = color_map[2][x];
}
static inline float toLin(float c){ return powf(c/255.f, 2.2f); }
static inline float toGam(float l){ return 255.f*powf(l<0?0:l, 1/2.2f); }

// mode 0: smooth iteration -> palette via log-power. mode 1: stripe average.
static void pixel(double cr, double ci, int mxit, int mode, float& R, float& G, float& B) {
    double zr=0, zi=0, sum=0, lastAdd=0; int cnt=0;
    const double freq = 7.0;
    const double R2 = 1e20, logR = log(sqrt(R2));
    for (int i=0;i<mxit;++i){
        double zr2=zr*zr, zi2=zi*zi;
        double nzr=zr2-zi2+cr, nzi=2*zr*zi+ci; zr=nzr; zi=nzi;
        double add = 0.5 + 0.5*sin(freq*atan2(zi,zr)); lastAdd=add; sum+=add; ++cnt;
        double m2=zr*zr+zi*zi;
        if (m2 > R2){
            if (mode==1){
                double aw=sum/cnt, awo=cnt>1?(sum-lastAdd)/(cnt-1):aw;
                double frac = 1.0 - log(log(sqrt(m2))/logR)/log(2.0);
                if (frac<0)frac=0; if(frac>1)frac=1;
                double avg = awo + (aw-awo)*frac;
                palColor((float)fmod(avg,1.0), R,G,B); return;
            } else {
                double it = i + 1 - log(log(m2)/2/log(2))/log(2);
                float l=logf((float)it+2); float t=powf(l, l*l/60.f);
                palColor(t - floorf(t), R,G,B); return;
            }
        }
    }
    R=G=B=0;
}

static void writeBMP(const char* p, const std::vector<uint8_t>& rgb, int W, int H){
    int row=(W*3+3)&~3, ds=row*H; uint8_t fh[14]={'B','M'}; uint32_t fs=14+40+ds; memcpy(fh+2,&fs,4); uint32_t off=54; memcpy(fh+10,&off,4);
    uint8_t ih[40]={0}; uint32_t v; v=40;memcpy(ih,&v,4); v=W;memcpy(ih+4,&v,4); v=H;memcpy(ih+8,&v,4); uint16_t pl=1,bp=24; memcpy(ih+12,&pl,2); memcpy(ih+14,&bp,2); v=ds;memcpy(ih+20,&v,4);
    FILE* f=fopen(p,"wb"); fwrite(fh,1,14,f); fwrite(ih,1,40,f); std::vector<uint8_t> ln(row,0);
    for(int y=H-1;y>=0;--y){ for(int x=0;x<W;++x){ const uint8_t* q=&rgb[(y*W+x)*3]; ln[x*3]=q[2]; ln[x*3+1]=q[1]; ln[x*3+2]=q[0]; } fwrite(ln.data(),1,row,f); } fclose(f);
}

int main(int argc, char** argv){
    if (argc < 10){ fprintf(stderr,"usage: coloring_demo out W H cx cy scale mxit palIdx mode [SS]\n"); return 1; }
    const char* out=argv[1]; int W=atoi(argv[2]), H=atoi(argv[3]);
    double cx=atof(argv[4]), cy=atof(argv[5]), scale=atof(argv[6]); int mxit=atoi(argv[7]);
    int palIdx=atoi(argv[8]), mode=atoi(argv[9]); int SS=argc>10?atoi(argv[10]):3;
    if (palIdx<0||palIdx>=(int)pals.size()) palIdx=0;
    rebuild(pals[palIdx]);
    int Ws=W*SS, Hs=H*SS;
    double half = 2.0/scale;          // view half-width in complex units (approx)
    std::vector<uint8_t> img(W*H*3);
    #pragma omp parallel for schedule(dynamic,1)
    for (int i=0;i<H;++i) for (int j=0;j<W;++j){
        float lr=0,lg=0,lb=0;
        for (int a=0;a<SS;++a) for (int b=0;b<SS;++b){
            double px = (j*SS+b+0.5)/Ws, py=(i*SS+a+0.5)/Hs;
            double cr = cx + (px*2-1)*half;
            double ci = cy - (py*2-1)*half*H/W;
            float R,G,B; pixel(cr,ci,mxit,mode,R,G,B); lr+=toLin(R); lg+=toLin(G); lb+=toLin(B);
        }
        int n=SS*SS; uint8_t* p=&img[(i*W+j)*3];
        p[0]=(uint8_t)(toGam(lr/n)+0.5f); p[1]=(uint8_t)(toGam(lg/n)+0.5f); p[2]=(uint8_t)(toGam(lb/n)+0.5f);
    }
    writeBMP(out, img, W, H);
    printf("wrote %s mode=%s\n", out, mode?"stripe-average":"smooth-iter");
    return 0;
}
