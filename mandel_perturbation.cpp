#include <iostream>
#include <iomanip>
#include <complex>
#include <mpir.h>
#include <assert.h>
#include <chrono>

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <array>
#include <utility>
#include <cstring>

#include "mandel_perturbation.h"
#include "float_math.h"
#include "dll_interface.h"


Mandel::Mandel(int width, int height, int max_iteration, int sub, float* iter) : _w(width), _h(height), _mxit(max_iteration), _sub(sub), _iter(iter) {
    assert(width > 0);
    assert(height > 0);
    assert(max_iteration > 0);
    assert(sub % 2);
    _z_re = new mpf_t[_mxit + 1];
    _z_im = new mpf_t[_mxit + 1];
    _zf = new Comp[_mxit + 1];
    _zfr = new Float[_mxit + 1];
    _zfi = new Float[_mxit + 1];

    _d_re = new mpf_t[_mxit + 1];
    _d_im = new mpf_t[_mxit + 1];
    _df = new Comp[_mxit + 1];
    _dfr = new Float[_mxit + 1];
    _dfi = new Float[_mxit + 1];

    _delta_0 = new Comp[_mxit + 1];
    _done = new bool[width * height * sub * sub];
    _z_m3 = new Float[_mxit + 1];

    mpf_init(_c0_re);
    mpf_init(_c0_im);
    mpf_init(_dx);
    mpf_init(_dy);
    mpf_init(_scale);
    mpf_init(_ref_z_re);
    mpf_init(_ref_z_im);

    mpf_init(_t1);
    mpf_init(_t2);
    for (int i = 0; i <= _mxit; ++i) {
        mpf_init(_z_re[i]);
        mpf_init(_z_im[i]);
        mpf_init(_d_re[i]);
        mpf_init(_d_im[i]);
    }
}

Mandel::~Mandel() {
    delete[] _delta_0;
    delete[] _done;
    delete[] _z_m3;

    mpf_clear(_c0_re);
    mpf_clear(_c0_im);
    mpf_clear(_dx);
    mpf_clear(_dy);
    mpf_clear(_scale);
    mpf_clear(_ref_z_re);
    mpf_clear(_ref_z_im);

    mpf_clear(_t1);
    mpf_clear(_t2);
    for (int i = 0; i <= _mxit; ++i) {
        mpf_clear(_z_re[i]);
        mpf_clear(_z_im[i]);
        mpf_clear(_d_re[i]);
        mpf_clear(_d_im[i]);
    }
    delete[] _zf;
    delete[] _z_re;
    delete[] _z_im;

    delete[] _df;
    delete[] _d_re;
    delete[] _d_im;
}

bool Mandel::attractor(double zz_re, double zz_im, const double c_re, const double c_im, int period) const {
    const double epsilon = 1e-20;
    double z_re, z_im;
    double dz_re, dz_im;
    double zz1_re, zz1_im;
    double tmp, t2;
    for (int j = 0; j < 64; ++j) {
        z_re = zz_re;
        z_im = zz_im;
        dz_re = 1;
        dz_im = 0;
        for (int i = 0; i < period; ++i) {
            tmp = (z_re * dz_re - z_im * dz_im) * 2;
            dz_im = (z_re * dz_im + z_im * dz_re) * 2;
            dz_re = tmp;

            tmp = (z_re * z_re - z_im * z_im) + c_re;
            z_im = (z_re * z_im) * 2 + c_im;
            z_re = tmp;
        }
        z_re -= zz_re;
        z_im -= zz_im;
        tmp = dz_re - 1;
        t2 = tmp * tmp + dz_im * dz_im;
        zz1_re = (z_re * tmp + z_im * dz_im) / t2;
        zz1_im = (z_im * tmp - z_re * dz_im) / t2;
        if (zz1_re * zz1_re + zz1_im * zz1_im < epsilon) {
            return (dz_re * dz_re + dz_im * dz_im <= 1);
        }
        zz_re -= zz1_re;
        zz_im -= zz1_im;
    }
    return false;
}


double Mandel::floatPointCompute(Float c_re, Float c_im, int mxit, int c_method) const {
    Float z_re = c_re;
    Float z_im = c_im;
    Float d_re = 2;
    Float d_im = 0;
    Float dc_re = 1;
    Float dc_im = 0;
    Float tmp;
    int i = 1;
    while (i < mxit) {
        tmp = 2.0 * (d_re * z_re - d_im * z_im);
        d_im = 2.0 * (d_re * z_im + d_im * z_re);
        d_re = tmp;

        if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
            tmp = 2.0 * (dc_re * z_re - dc_im * z_im) + 1;
            dc_im = 2.0 * (dc_re * z_im + dc_im * z_re);
            dc_re = tmp;
        }
        
        tmp = z_re * z_re - z_im * z_im + c_re;
        z_im = 2.0 * z_re * z_im + c_im;
        z_re = tmp;

        tmp = z_re * z_re + z_im * z_im;
        if (tmp > _ESCAPE_RADIUS * _ESCAPE_RADIUS) {
            if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
                return sqrt(tmp) * log(tmp) / sqrt(dc_re * dc_re + dc_im * dc_im);
            } else {
                return (i + 1 - log(log(tmp) / 2 / log(2)) / log(2));
            }
        }

        tmp = d_re * d_re + d_im * d_im;
        if (tmp < 0.000000001) return -2;
        ++i;
    }
    return -2;

}

float Mandel::accuratePointCompute(mpf_t c_re, mpf_t c_im, int mxit, int c_method) const {
    mpf_t t1, t2;
    mpf_init(t1);
    mpf_init(t2);

    // z = c;
    mpf_t z_re, z_im;
    mpf_init_set(z_re, c_re);
    mpf_init_set(z_im, c_im);
    
    // d = 2.0;
    mpf_t d_re, d_im;
    mpf_init(d_re);
    mpf_init(d_im);
    mpf_set_str(d_re, "2", 10);
    int i = 1;
    float res = -2;
    while (i < mxit) {
        // d = 2.0 * d * z;
        mpf_mul(t1, d_re, z_im);
        mpf_mul(d_re, d_re, z_re);
        mpf_mul(t2, d_im, z_im);
        mpf_mul(d_im, d_im, z_re);
        mpf_sub(d_re, d_re, t2);
        mpf_mul_ui(d_re, d_re, 2);
        mpf_add(d_im, d_im, t1);
        mpf_mul_ui(d_im, d_im, 2);

        // z = z * z + c;
        mpf_mul(t1, z_re, z_im);
        mpf_mul(z_re, z_re, z_re);
        mpf_mul(z_im, z_im, z_im);
        mpf_sub(z_re, z_re, z_im);
        mpf_add(z_re, z_re, c_re);
        mpf_mul_ui(z_im, t1, 2);
        mpf_add(z_im, z_im, c_im);

        // auto re_im = z.get_real_imag();
        // if (re_im.first * re_im.first + re_im.second * re_im.second > _ESCAPE_RADIUS) return (i + 1 - log(log(static_cast<float>(re_im.first * re_im.first + re_im.second * re_im.second))) / log(2));
        mpf_mul(t1, z_re, z_re);
        mpf_mul(t2, z_im, z_im);
        mpf_add(t1, t1, t2);
            
        
        double rad = mpf_get_d(t1);

        mpf_sub(t1, z_re, _z_re[i]);
        mpf_sub(t2, z_im, _z_im[i]);

        if (rad > _ESCAPE_RADIUS) {
            res = (i + 1 - log(log(rad) / 2 / log(2)) / log(2));
            break;
        }

        // re_im = d.get_real_imag();
        // if (re_im.first * re_im.first + re_im.second * re_im.second < 0.000000001) return -2;
        mpf_mul(t1, d_re, d_re);
        mpf_mul(t2, d_im, d_im);
        mpf_add(t1, t1, t2);
        rad = mpf_get_d(t1);
        if (rad < 0.000000001) break;
        
        ++i;
    }
    mpf_clear(d_re);
    mpf_clear(d_im);
    mpf_clear(z_re);
    mpf_clear(z_im);
    mpf_clear(t1);
    mpf_clear(t2);
    return res;
}


void Mandel::Compute(mpf_t c_re, mpf_t c_im, mpf_t scale, int mxit, int c_method) {
    std::cout << "mxit: " << mxit << '\n';
    _mx_coef = -1;
    _ref_cnt = 0;
    int iteration = 0;

    // _scale = scale;
    mpf_set(_scale, scale);

    mpf_t dw, dh;
    mpf_init_set_ui(dw, 2);
    mpf_div(dw, dw, scale);
    mpf_init_set(dh, dw);
    mpf_div_ui(dh, dh, _w);
    mpf_mul_ui(dh, dh, _h);
    mpf_sub(_c0_re, c_re, dw);
    mpf_sub(_c0_im, c_im, dh);

    mpf_div_ui(_dx, dw, _w - 1);
    mpf_mul_ui(_dx, _dx, 2);
    mpf_div_ui(_dy, dh, _h - 1);
    mpf_mul_ui(_dy, _dy, 2);

    mpf_clear(dw);
    mpf_clear(dh);
    
    std::set<std::array<int, 4>> s;

    int method = 0;

    // Shallow zoom: method = 0, basic algorithm
    // Deep zoom: method = 1, Perturbation + Series Approximation + Rebase
    if (mpf_cmp_d(scale, 1e6) > 0) method = 1;
    
    if (_flag_halt) return;
    int ref_it = 0, pr_it = 0;
    if (method == 0) {
        Float c0_re_f = mpf_get_ld(_c0_re);
        Float c0_im_f = mpf_get_ld(_c0_im);
        Float dx_f = mpf_get_ld(_dx);
        Float dy_f = mpf_get_ld(_dy);
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < _h; ++i) {
            if (_flag_halt) continue;
            for (int j = 0; j < _w; ++j) {
                if (_flag_halt) break;
                int idx = getIndex(i, j, 0, 0);
                _iter[idx] = floatPointCompute(c0_re_f + dx_f * j, c0_im_f + dy_f * i, mxit, c_method);
                if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
                    if (_iter[idx] >= 0) _iter[idx] /= dx_f;
                }
            }
        }
    }
    else if (method == 1) {
        Float c0_re_f = mpf_get_ld(_c0_re);
        Float c0_im_f = mpf_get_ld(_c0_im);
        floatPointCompute(c0_re_f, c0_im_f, mxit, c_method);
        s.insert({ _h / 2, _w / 2, 0, 0 });
        ref_it = createRef(s, mxit, mxit, true, c_method);
        for (int i = 0; i < _h; ++i) {
            for (int j = 0; j < _w; ++j) {
                s.insert({ i, j, 0, 0 });
            }
        }
        s.erase({ _w / 2, _h / 2, 0, 0 });
        
        while (!s.empty()) {
            // std::cout << ref_it << " === " << _SA_it << ' ' << _SA_order << ' ';
            if (_flag_halt) return;
            stepParallel(s, ref_it, mxit, c_method);
            if (_flag_halt) return;
            if (s.empty()) break;
            ref_it = createRef(s, mxit, mxit, false, c_method);
        }
    }
    bool is_super_sampling = (_sub > 1) && (c_method & ColoringMethod::SUPER_SAMPLING);
    if (!is_super_sampling) return;
    if (_flag_halt) return;
    
    // Oversampling pixels that differs from neighbours with (_sub x _sub) pixels
    std::vector<std::array<int, 2>> v;
    int mix_cnt = 0;
    double mx_diff = log(mxit) / 8;
    if (c_method & ColoringMethod::SUPER_SAMPLING) mx_diff = 1;
    for (int i = 0; i < _h; ++i) {
        for (int j = 0; j < _w; ++j) {
            bool need_sample = false;
            if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
                for (int xi = -1; xi <= 1; ++xi) {
                    for (int yi = -1; yi <= 1; ++yi) {
                        if (xi == 0 && yi == 0) continue;
                        int ny = i + yi;
                        int nx = j + xi;
                        if (nx >= 0 && nx < _w && ny >= 0 && ny < _h) {
                            if (_iter[getIndex(i, j, 0, 0)] * _iter[getIndex(ny, nx, 0, 0)] < 0) {
                                need_sample = true;
                                break;
                            }
                        }
                    }
                    if (need_sample) break;
                }
            } else {
                float diff = 0;
                for (int xi = -1; xi <= 1; ++xi) {
                    for (int yi = -1; yi <= 1; ++yi) {
                        if (xi == 0 && yi == 0) continue;
                        int ny = i + yi;
                        int nx = j + xi;
                        if (nx >= 0 && nx < _w && ny >= 0 && ny < _h) {
                            diff = std::max(diff, std::abs(_iter[getIndex(i, j, 0, 0)] - _iter[getIndex(ny, nx, 0, 0)]));
                        }
                    }
                }
                need_sample = (log(diff) > log(mxit) / 8);
            }
            if (need_sample) {
                ++mix_cnt;
                v.push_back({ i, j });
                for (int xi = -_sub / 2; xi <= _sub / 2; ++xi) {
                    for (int yi = -_sub / 2; yi <= _sub / 2; ++yi) {
                        if (xi == 0 && yi == 0) continue;
                        s.insert({ i, j, yi, xi });
                    }
                }
            }
        }
    }
    // printf("Mixing %0.2f%% pixels\n", 1.f * mix_cnt / _w / _h * 100);
    
    if (_flag_halt) return;
    if (method == 0) {
        Float c0_re_f = mpf_get_ld(_c0_re);
        Float c0_im_f = mpf_get_ld(_c0_im);
        Float dx_f = mpf_get_ld(_dx);
        Float dy_f = mpf_get_ld(_dy);
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < v.size(); ++i) {
            if (_flag_halt) continue;
            for (int xi = -_sub / 2; xi <= _sub / 2; ++xi) {
                for (int yi = -_sub / 2; yi <= _sub / 2; ++yi) {
                    if (xi == 0 && yi == 0) continue;
                    std::array<int, 4> arr = { v[i][0], v[i][1], yi, xi };
                    double iteration = floatPointCompute(c0_re_f + dx_f * v[i][1] + dx_f * xi / _sub, c0_im_f + dy_f * v[i][0] + dy_f * yi / _sub, mxit, c_method);
                    if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
                        if (iteration >= 0) iteration /= dx_f;
                    }
                    setPixel(arr, iteration);
                }
            }
        }
    }
    else if (method == 1) {
        while (!s.empty()) {
            // std::cout << ref_it << " === " << _SA_it << ' ' << _SA_order << ' ';
            if (_flag_halt) return;
            stepParallel(s, ref_it, mxit, c_method);
            if (_flag_halt) return;
            if (s.empty()) break;
            ref_it = createRef(s, mxit, mxit, false, c_method);
        }
    }

}


void Mandel::setPrecision(int precision) {
    mpf_set_prec(_c0_re, precision);
    mpf_set_prec(_c0_im, precision);
    mpf_set_prec(_dx, precision);
    mpf_set_prec(_dy, precision);
    mpf_set_prec(_scale, precision);
    mpf_set_prec(_ref_z_re, precision);
    mpf_set_prec(_ref_z_im, precision);

    mpf_set_prec(_t1, precision);
    mpf_set_prec(_t2, precision);
    for (int i = 0; i <= _mxit; ++i) {
        mpf_set_prec(_z_re[i], precision);
        mpf_set_prec(_z_im[i], precision);
    }
}

inline bool Mandel::escape(Comp& z) const { return (std::abs(z) > _ESCAPE_RADIUS); }
inline bool Mandel::escape(mpf_t& z_re, mpf_t& z_im) const {
    // auto ri = z.get_real_imag();
    mpf_t t1, t2;
    mpf_init(t1);
    mpf_init(t2);
    mpf_mul(t1, z_re, z_re);
    mpf_mul(t2, z_im, z_im);
    mpf_add(t1, t1, t2);
    double rad = mpf_get_d(t1);
    mpf_clear(t1);
    mpf_clear(t2);
    return (rad > _ESCAPE_RADIUS);
}

inline float Mandel::getEscapeTime(Comp& z, int i) const {
    double z_re = (double)z.real();
    double z_im = (double)z.imag();
    double rad = z_re * z_re + z_im * z_im;
    return (i + 1 - log(log(rad) / 2 / log(2)) / log(2));
}

inline float Mandel::getEscapeTime(mpf_t& z_re, mpf_t& z_im, int i) {
    // auto ri = z.get_real_imag();
    mpf_mul(_t1, z_re, z_re);
    mpf_mul(_t2, z_im, z_im);
    mpf_add(_t1, _t1, _t2);
    double rad = mpf_get_d(_t1);
    return (i + 1 - log(log(rad) / 2 / log(2)) / log(2));
}

void Mandel::setPixel(std::array<int, 4> p, float iteration) const {
    _iter[getIndex(p)] = iteration;
}

void Mandel::stepParallel(std::set<std::array<int, 4>>& s, int mx_ref_it, int mxit, int c_method) {
    std::vector<std::array<int, 4>> v;
    int n = s.size();
    for (auto& p : s) v.push_back(p);
    memset(_done, 0, n * sizeof(bool));
    auto glitch_p = std::unique_ptr<Float[]>(new Float[n]);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < v.size(); ++i) {
        if (_flag_halt) continue;
        auto arr = v[i];
        mpf_t t1, t2;
        mpf_init(t1);
        mpf_init(t2);
        glitch_p[i] = -1;
        // Comp dc = static_cast<Comp>(_dx * (v[i][1] - _ref[1]) + _dx / _sub * (v[i][3] - _ref[3]) + _dy * (v[i][0] - _ref[0]) + _dy / _sub * (v[i][2] - _ref[2]));
        int xpix = _sub * (arr[1] - _ref[1]) + (arr[3] - _ref[3]);
        int ypix = _sub * (arr[0] - _ref[0]) + (arr[2] - _ref[2]);
        mpf_mul_ui(t1, _dx, abs(xpix));
        mpf_mul_ui(t2, _dy, abs(ypix));
        if (xpix < 0) mpf_neg(t1, t1);
        if (ypix < 0) mpf_neg(t2, t2);
        if (_sub > 1) {
            mpf_div_ui(t1, t1, _sub);
            mpf_div_ui(t2, t2, _sub);
        }
        Comp dc{ mpf_get_ld(t1), mpf_get_ld(t2) };
        mpf_clear(t1);
        mpf_clear(t2);
        Float dxf = mpf_get_ld(_dx);
        Comp dz = { 0 };
        Comp dd = { 0 };
        for (int x = _SA_order; x >= 0; --x) {
            dz += _Adf_old[x];
            dz *= dc;
            dz /= _SA_delta;
        }
        if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
            for (int x = _SA_order; x >= 0; --x) {
                dd += _Bdf_old[x];
                dd *= dc;
                dd /= _SA_delta;
            }
        }
        int j = _SA_it + 1;

        Float dzr = dz.real();
        Float dzi = dz.imag();
        Float dcr = dc.real();
        Float dci = dc.imag();
        Float ddr = dd.real();
        Float ddi = dd.imag();
        Float tmp;
        Float zr, zi, dr, di;
        Float dzr2, dzi2, zrad;
        

        int k = j;
        Float m = 1e100;
        while (j < mxit) {
            if (_flag_halt) break;
            // dd = 2 * (dd*z + dz*(d+dd))
            if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
                if (k == 0) {
                    tmp = (dzr * (1 + ddr) - dzi * (1 + ddi)) * 2;
                    ddi = (dzr * (1 + ddi) + dzi * (1 + ddr)) * 2;
                } else {
                    tmp = (ddr * _zfr[k - 1] - ddi * _zfi[k - 1] + dzr * (_dfr[k - 1] + ddr) - dzi * (_dfi[k - 1] + ddi)) * 2;
                    ddi = (ddr * _zfi[k - 1] + ddi * _zfr[k - 1] + dzr * (_dfi[k - 1] + ddi) + dzi * (_dfr[k - 1] + ddr)) * 2;
                }
                ddr = tmp;
                dr = ddr + _dfr[k];
                di = ddi + _dfi[k];
            }
            // dz = dz^2 + 2*dz*z + dc
            dzr2 = dzr * dzr;
            dzi2 = dzi * dzi;
            if (k == 0) {
                tmp = dzr2 - dzi2 + dcr;
                dzi = dzr * dzi * 2 + dci;
            }
            else {
                tmp = (dzr * _zfr[k - 1] - dzi * _zfi[k - 1]) * 2 + dzr2 - dzi2 + dcr;
                dzi = (dzr * _zfi[k - 1] + dzi * _zfr[k - 1]) * 2 + dzr * dzi * 2 + dci;
            }
            dzr = tmp;
            zr = dzr + _zfr[k];
            zi = dzi + _zfi[k];
            zrad = zr * zr + zi * zi;


            ++k;
            if (zrad > _ESCAPE_RADIUS * _ESCAPE_RADIUS) {
                
                if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
                    setPixel(arr, sqrt(zrad) / dxf * log(zrad) / sqrt(dr * dr + di * di));
                } else {
                    setPixel(arr, j + 1 - log(log((double)zrad) / 2 / log(2)) / log(2));
                }
                _done[i] = true;
                break;
            }
            // ** Zhuoran method: rebase to original reference with index = 0
            if (zrad < dzr * dzr + dzi * dzi || k == mx_ref_it) {
                if ((dzr * dzr + dzi * dzi) / zrad > 10000000) {
                    // Significant magnitude change causes precision loss.
                    glitch_p[i] = zrad / _z_m3[j];
                    break;
                }
                dzr = zr;
                dzi = zi;
                if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
                    ddr = dr - 1;
                    ddi = di;
                }
                k = 0;
            }
            if (k == mx_ref_it) {
                glitch_p[i] = (double)zrad / _z_m3[j];
                break;
            }
            ++j;
        }
        if (j >= mxit) {
            setPixel(arr, -2);
            _done[i] = true;
        }
    }
    if (_flag_halt) return;
    // get new reference with minimum |z + dz|
    Float min_orbit = 1e100;
    for (int i = 0; i < n; ++i) {
        if (glitch_p[i] >= 0 && glitch_p[i] < min_orbit) {
            min_orbit = glitch_p[i];
            _new_ref = v[i];
        }
    }

    int i = 0;
    int remove_cnt = 0;
    for (auto iter = s.begin(); iter != s.end(); ++i) {
        if (_done[i]) {
            iter = s.erase(iter);
            ++remove_cnt;
        }
        else {
            ++iter;
        }
    }

    // printf("%.2f%% %d\n", 100.0 * remove_cnt / (_w * _h), remove_cnt);
}

int Mandel::createRef(std::set<std::array<int, 4>>& s, int pr_it, int mxit, bool random, int c_method) {
    random = true;
    if (s.empty()) return false;
    ++_ref_cnt;
    std::array<int, 4> p;
    if (random) {
        auto r = rand() % s.size();
        auto iter = std::begin(s);
        std::advance(iter, r);
        p = *iter;
        s.erase(iter);
    }
    else {
        s.erase(_new_ref);
        p = _new_ref;
    }

    // _ref_z = _c0 + p[0] * _dy + p[2] * _dy / _sub + p[1] * _dx + p[3] * _dx / _sub;
    mpf_set(_ref_z_re, _c0_re);
    mpf_set(_ref_z_im, _c0_im);
    mpf_mul_ui(_t1, _dx, p[1]);
    mpf_add(_ref_z_re, _ref_z_re, _t1);
    mpf_mul_ui(_t1, _dx, abs(p[3]));
    if (p[3] < 0) mpf_neg(_t1, _t1);
    mpf_div_ui(_t1, _t1, _sub);
    mpf_add(_ref_z_re, _ref_z_re, _t1);
    mpf_mul_ui(_t1, _dy, p[0]);
    mpf_add(_ref_z_im, _ref_z_im, _t1);
    mpf_mul_ui(_t1, _dy, abs(p[2]));
    if (p[2] < 0) mpf_neg(_t1, _t1);
    mpf_div_ui(_t1, _t1, _sub);
    mpf_add(_ref_z_im, _ref_z_im, _t1);

    // _ref_z_f = static_cast<Comp>(_ref_z);
    _ref_z_f = Comp{ mpf_get_ld(_ref_z_re), mpf_get_ld(_ref_z_im) };
    
    _ref = p;
    
    // _z[0] = _ref_z;
    mpf_set(_z_re[0], _ref_z_re);
    mpf_set(_z_im[0], _ref_z_im);

    _zf[0] = _ref_z_f;
    _zfr[0] = _ref_z_f.real();
    _zfi[0] = _ref_z_f.imag();

    mpf_set_ui(_d_re[0], 1);
    mpf_set_ui(_d_im[0], 0);
    _df[0] = Comp{ 1 };

    // _SA_delta = static_cast<Comp>(_dx * _w / 2 + _dy * _h / 2);
    mpf_mul_ui(_t1, _dx, _w);
    mpf_div_ui(_t1, _t1, 2);
    mpf_mul_ui(_t2, _dy, _h);
    mpf_div_ui(_t2, _t2, 2);
    _SA_delta = { mpf_get_ld(_t1), mpf_get_ld(_t2) };
    for (int i = 0; i < _SA_N; ++i) _Adf_old[i] = _Bdf_old[i] = 0;
    _Adf_old[0] = _SA_delta;

    _SA_flag = true;
    _SA_it = 0;

    for (int i = 1; i <= mxit; ++i) {
        if (_flag_halt) break;
        if (!calCoefficient(i, pr_it, c_method)) {
            setPixel(_ref, getEscapeTime(_z_re[i], _z_im[i], i));
            return i;
        }
    }
    return mxit;
}

bool Mandel::calCoefficient(int i, int pr_it, int c_method) {
    // _z[i] = _z[i - 1] * _z[i - 1] + _ref_z;
    mpf_mul(_z_re[i], _z_re[i - 1], _z_re[i - 1]);
    mpf_mul(_t1, _z_im[i - 1], _z_im[i - 1]);
    mpf_sub(_z_re[i], _z_re[i], _t1);
    mpf_mul(_z_im[i], _z_re[i - 1], _z_im[i - 1]);
    mpf_mul_ui(_z_im[i], _z_im[i], 2);
    mpf_add(_z_re[i], _z_re[i], _ref_z_re);
    mpf_add(_z_im[i], _z_im[i], _ref_z_im);

    // _zf[i] = static_cast<Comp>(_z[i]);
    _zf[i] = Comp{ mpf_get_ld(_z_re[i]), mpf_get_ld(_z_im[i]) };
    _zfr[i] = _zf[i].real();
    _zfi[i] = _zf[i].imag();

    // _z_m3[i] = tmp_z.abs().get_real_imag().first / 1000; // for Pauldelbrot condition
    _z_m3[i] = { (_zfr[i] * _zfr[i] + _zfi[i] * _zfi[i]) / 1000000 };

    if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
        // _d[i] = 2 * _d[i - 1] * _z[i - 1] + 1;
        mpf_mul(_d_re[i], _d_re[i - 1], _z_re[i - 1]);
        mpf_mul(_t1, _d_im[i - 1], _z_im[i - 1]);
        mpf_sub(_d_re[i], _d_re[i], _t1);
        mpf_mul_ui(_d_re[i], _d_re[i], 2);
        mpf_add_ui(_d_re[i], _d_re[i], 1);
        mpf_mul(_d_im[i], _d_re[i - 1], _z_im[i - 1]);
        mpf_mul(_t1, _d_im[i - 1], _z_re[i - 1]);
        mpf_add(_d_im[i], _d_im[i], _t1);
        mpf_mul_ui(_d_im[i], _d_im[i], 2);

        
        _df[i] = Comp{ mpf_get_ld(_d_re[i]), mpf_get_ld(_d_im[i]) };
        _dfr[i] = _df[i].real();
        _dfi[i] = _df[i].imag();
    }

    if (escape(_zf[i])) return false;
    if (i <= pr_it) {
        if (_SA_flag) {
            for (int j = 0; j < _SA_N; ++j) {
                _Adf_new[j] = Comp{ 2 } * _zf[i - 1] * _Adf_old[j];
                for (int k = 0; k < j / 2; ++k) {
                    _Adf_new[j] += Comp{ 2 } * _Adf_old[k] * _Adf_old[j - k - 1];
                }
                if (j % 2) _Adf_new[j] += _Adf_old[j / 2] * _Adf_old[j / 2];
            }
            
            if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
                for (int j = 0; j < _SA_N; ++j) {
                    _Bdf_new[j] = _Bdf_old[j] * _zf[i - 1] + _Adf_old[j] * _df[i - 1];
                    for (int k = 0; k < j; ++k) {
                        _Bdf_new[j] += _Adf_old[k] * _Bdf_old[j - k - 1];
                    }
                    _Bdf_new[j] *= Comp{ 2 };
                }
            }
            _Adf_new[0] += _SA_delta;
            int order = SACheckMagnitude();
            
            if (order < 0) {
                _SA_flag = false;
            }
            else {
                _SA_it = i;
                _SA_order = order;
                for (int j = 0; j < _SA_N; ++j) _Adf_old[j] = _Adf_new[j];
                if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
                    for (int j = 0; j < _SA_N; ++j) _Bdf_old[j] = _Bdf_new[j];
                }
            }
        }
        _mx_coef = i;
    }
    return true;
}

int Mandel::SACheckMagnitude() const {
    int pre_mn[_SA_N], suf_mx[_SA_N];
    int re, im;
    re = get_exp(_Adf_new[0].real());
    im = get_exp(_Adf_new[0].imag());
    pre_mn[0] = std::min(re, im);
    
    re = get_exp(_Adf_new[_SA_N - 1].real());
    im = get_exp(_Adf_new[_SA_N - 1].imag());
    suf_mx[_SA_N - 1] = std::max(re, im);
    for (int i = 1; i < _SA_N; ++i) {
        re = get_exp(_Adf_new[i].real());
        im = get_exp(_Adf_new[i].imag());
        pre_mn[i] = std::max(pre_mn[i - 1], std::min(re, im));
    }
    for (int i = _SA_N - 2; i >= 0; --i) {
        re = get_exp(_Adf_new[i].real());
        im = get_exp(_Adf_new[i].imag());
        suf_mx[i] = std::max(suf_mx[i + 1], std::max(re, im));
    }
    
    for (int i = 0; i < _SA_N - 10; ++i) {
        if (pre_mn[i] - suf_mx[i + 1] >= 80) return std::max(i, 5);
    }
    return -1;
}

inline int Mandel::getIndex(std::array<int, 4>& arr) const {
    return (arr[0] * _sub + _sub / 2 + arr[2]) * _w * _sub + (arr[1] * _sub + _sub / 2 + arr[3]);
}

inline int Mandel::getIndex(int i, int j, int u, int v) const {
    return (i * _sub + _sub / 2 + u) * _w * _sub + (j * _sub + _sub / 2 + v);
}

void Mandel::SetHalt(bool flag) {
    _flag_halt = flag;
}
