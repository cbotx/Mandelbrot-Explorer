// Isolated experiment for finding a minibrot nucleus and measuring the GMP
// reference-orbit work removed by storing one period.  This does not alter or
// exercise the production renderer.
//
//   periodic_ref_bench find  cx cy scaleExp [maxPeriod] [precisionBits]
//   periodic_ref_bench bench cx cy scaleExp mxit [maxPeriod] [reps] [precisionBits]
//
// scaleExp is the decimal zoom exponent.  As in the original period tool, the
// searched parameter disk has radius 2 * 10^-scaleExp.

#include <gmp.h>

#include "floatexp.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#include <vector>

static mp_bitcnt_t g_precision_bits = 256;
static volatile uint64_t g_replay_sink = 0;

struct Mpf {
    mpf_t v;
    Mpf() { mpf_init2(v, g_precision_bits); }
    ~Mpf() { mpf_clear(v); }
    Mpf(const Mpf&) = delete;
    Mpf& operator=(const Mpf&) = delete;
};

struct Cmpf {
    mpf_t re, im;
    Cmpf() {
        mpf_init2(re, g_precision_bits);
        mpf_init2(im, g_precision_bits);
    }
    ~Cmpf() {
        mpf_clear(re);
        mpf_clear(im);
    }
    Cmpf(const Cmpf&) = delete;
    Cmpf& operator=(const Cmpf&) = delete;
};

struct Scratch {
    Mpf a, b, c, d;
};

static double now_seconds() {
    return std::chrono::duration_cast<std::chrono::duration<double>>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}

static void complex_zero(Cmpf& z) {
    mpf_set_ui(z.re, 0);
    mpf_set_ui(z.im, 0);
}

static void complex_set(Cmpf& dst, const Cmpf& src) {
    mpf_set(dst.re, src.re);
    mpf_set(dst.im, src.im);
}

// In-place z = z^2 + c, using the same two-multiply Karatsuba square as the
// production reference builder.
static void complex_square_add(Cmpf& z, const Cmpf& c, Scratch& s) {
    mpf_add(s.a.v, z.re, z.im);
    mpf_sub(s.b.v, z.re, z.im);
    mpf_mul(s.c.v, z.re, z.im);
    mpf_mul(s.d.v, s.a.v, s.b.v);
    mpf_mul_ui(s.c.v, s.c.v, 2);
    mpf_add(z.re, s.d.v, c.re);
    mpf_add(z.im, s.c.v, c.im);
}

static void complex_square_add_to(Cmpf& out, const Cmpf& z, const Cmpf& c,
                                  Scratch& s) {
    mpf_add(s.a.v, z.re, z.im);
    mpf_sub(s.b.v, z.re, z.im);
    mpf_mul(s.c.v, z.re, z.im);
    mpf_mul(s.d.v, s.a.v, s.b.v);
    mpf_mul_ui(s.c.v, s.c.v, 2);
    mpf_add(out.re, s.d.v, c.re);
    mpf_add(out.im, s.c.v, c.im);
}

static void complex_abs(mpf_ptr out, const Cmpf& z, Scratch& s) {
    mpf_mul(s.a.v, z.re, z.re);
    mpf_mul(s.b.v, z.im, z.im);
    mpf_add(s.a.v, s.a.v, s.b.v);
    mpf_sqrt(out, s.a.v);
}

static void complex_distance(mpf_ptr out, const Cmpf& a, const Cmpf& b, Scratch& s) {
    mpf_sub(s.a.v, a.re, b.re);
    mpf_sub(s.b.v, a.im, b.im);
    mpf_mul(s.a.v, s.a.v, s.a.v);
    mpf_mul(s.b.v, s.b.v, s.b.v);
    mpf_add(s.a.v, s.a.v, s.b.v);
    mpf_sqrt(out, s.a.v);
}

// Evaluate z_p(c) and dz_p/dc from z_0=0.
static void evaluate_period(const Cmpf& c, uint64_t period, Cmpf& z, Cmpf& dz,
                            Scratch& s) {
    complex_zero(z);
    complex_zero(dz);
    for (uint64_t i = 0; i < period; ++i) {
        // dz' = 2*z*dz + 1.  Calculate it before overwriting z.
        mpf_mul(s.a.v, z.re, dz.re);
        mpf_mul(s.b.v, z.im, dz.im);
        mpf_sub(s.c.v, s.a.v, s.b.v);
        mpf_mul_ui(s.c.v, s.c.v, 2);
        mpf_add_ui(s.c.v, s.c.v, 1);

        mpf_mul(s.a.v, z.re, dz.im);
        mpf_mul(s.b.v, z.im, dz.re);
        mpf_add(s.d.v, s.a.v, s.b.v);
        mpf_mul_ui(s.d.v, s.d.v, 2);

        mpf_set(dz.re, s.c.v);
        mpf_set(dz.im, s.d.v);
        complex_square_add(z, c, s);
    }
}

enum class SearchStop {
    Found,
    EscapedDisk,
    MaxPeriod
};

struct SearchResult {
    SearchStop stop = SearchStop::MaxPeriod;
    uint64_t period = 0;
    uint64_t iterations = 0;
};

// Propagate the parameter disk through z -> z^2+c:
//   R' = (2|z| + R)R + r.
// All magnitudes and radii remain in GMP.  A small precision-scaled local
// roundoff allowance is propagated with the ball; mpf is not directed-rounding,
// so this remains a detector rather than a formal interval proof.
static SearchResult find_period(const Cmpf& center, mpf_srcptr radius,
                                uint64_t max_period) {
    Cmpf z;
    Scratch s;
    Mpf orbit_abs, next_radius, center_abs, escape_bound, round_unit, round_error;
    complex_zero(z);
    complex_abs(center_abs.v, center, s);

    Mpf ball_radius;
    mpf_set_ui(ball_radius.v, 0);
    mpf_set_ui(orbit_abs.v, 0);
    mpf_set_ui(round_unit.v, 1);
    const mp_bitcnt_t guard = g_precision_bits > 64 ? g_precision_bits - 32
                                                    : g_precision_bits / 2;
    mpf_div_2exp(round_unit.v, round_unit.v, guard);

    SearchResult result;
    for (uint64_t p = 1; p <= max_period; ++p) {
        // Exact disk image bound, plus an allowance for this GMP step's rounding.
        mpf_mul_ui(s.a.v, orbit_abs.v, 2);
        mpf_add(s.a.v, s.a.v, ball_radius.v);
        mpf_mul(next_radius.v, s.a.v, ball_radius.v);
        mpf_add(next_radius.v, next_radius.v, radius);

        mpf_mul(round_error.v, orbit_abs.v, orbit_abs.v);
        mpf_add(round_error.v, round_error.v, center_abs.v);
        mpf_add_ui(round_error.v, round_error.v, 1);
        mpf_mul(round_error.v, round_error.v, round_unit.v);
        mpf_add(next_radius.v, next_radius.v, round_error.v);
        mpf_set(ball_radius.v, next_radius.v);

        complex_square_add(z, center, s);
        complex_abs(orbit_abs.v, z, s);
        result.iterations = p;

        if (mpf_cmp(orbit_abs.v, ball_radius.v) <= 0) {
            result.stop = SearchStop::Found;
            result.period = p;
            return result;
        }

        // The complete propagated z disk is beyond a conservative escape bound:
        // |z|-R > |c_center|+r+2.
        mpf_add(escape_bound.v, ball_radius.v, center_abs.v);
        mpf_add(escape_bound.v, escape_bound.v, radius);
        mpf_add_ui(escape_bound.v, escape_bound.v, 2);
        if (mpf_cmp(orbit_abs.v, escape_bound.v) > 0) {
            result.stop = SearchStop::EscapedDisk;
            return result;
        }
    }
    return result;
}

struct NewtonResult {
    bool converged = false;
    bool singular = false;
    uint64_t iterations = 0;
    mp_bitcnt_t convergence_bits = 0;
};

static NewtonResult refine_nucleus(Cmpf& c, uint64_t period, uint64_t max_iterations,
                                   mpf_ptr residual) {
    Cmpf z, dz, delta, previous;
    Scratch s;
    Mpf denominator, step2, tolerance, tolerance2;

    NewtonResult result;
    result.convergence_bits = g_precision_bits > 64 ? g_precision_bits - 64
                                                    : g_precision_bits / 2;
    mpf_set_ui(tolerance.v, 1);
    mpf_div_2exp(tolerance.v, tolerance.v, result.convergence_bits);
    mpf_mul(tolerance2.v, tolerance.v, tolerance.v);

    for (uint64_t it = 0; it < max_iterations; ++it) {
        evaluate_period(c, period, z, dz, s);

        mpf_mul(s.a.v, dz.re, dz.re);
        mpf_mul(s.b.v, dz.im, dz.im);
        mpf_add(denominator.v, s.a.v, s.b.v);
        if (mpf_sgn(denominator.v) == 0) {
            result.singular = true;
            break;
        }

        // delta = z / dz = z*conj(dz)/|dz|^2.
        mpf_mul(s.a.v, z.re, dz.re);
        mpf_mul(s.b.v, z.im, dz.im);
        mpf_add(s.c.v, s.a.v, s.b.v);
        mpf_div(delta.re, s.c.v, denominator.v);

        mpf_mul(s.a.v, z.im, dz.re);
        mpf_mul(s.b.v, z.re, dz.im);
        mpf_sub(s.c.v, s.a.v, s.b.v);
        mpf_div(delta.im, s.c.v, denominator.v);

        complex_set(previous, c);
        mpf_sub(c.re, c.re, delta.re);
        mpf_sub(c.im, c.im, delta.im);
        result.iterations = it + 1;

        mpf_mul(s.a.v, delta.re, delta.re);
        mpf_mul(s.b.v, delta.im, delta.im);
        mpf_add(step2.v, s.a.v, s.b.v);
        const bool unchanged =
            mpf_cmp(c.re, previous.re) == 0 && mpf_cmp(c.im, previous.im) == 0;
        if (unchanged || mpf_cmp(step2.v, tolerance2.v) <= 0) {
            result.converged = true;
            break;
        }
    }

    evaluate_period(c, period, z, dz, s);
    complex_abs(residual, z, s);
    return result;
}

struct FeShadow {
    FloatExp re, im;
};

static FloatExp mpf_to_fe(mpf_srcptr value) {
    if (mpf_sgn(value) == 0) return FloatExp{ 0.0, 0 };
    long exponent;
    return FloatExp{ mpf_get_d_2exp(&exponent, value), static_cast<int64_t>(exponent) };
}

class ShadowStorage {
public:
    ShadowStorage(size_t count, bool deep)
        : doubles_(count), floatexp_(deep ? count : 0), deep_(deep) {}

    void store(size_t i, const Cmpf& z) {
        doubles_[i] = std::complex<double>(mpf_get_d(z.re), mpf_get_d(z.im));
        if (deep_) floatexp_[i] = FeShadow{ mpf_to_fe(z.re), mpf_to_fe(z.im) };
    }

    size_t size() const { return doubles_.size(); }
    bool deep() const { return deep_; }
    const std::complex<double>& value(size_t i) const { return doubles_[i]; }
    const FeShadow& fe_value(size_t i) const { return floatexp_[i]; }

private:
    std::vector<std::complex<double>> doubles_;
    std::vector<FeShadow> floatexp_;
    bool deep_;
};

// Exactly two rotating GMP complex values are retained.  Per-iteration storage
// is only the renderer-facing double/floatexp shadow.
static uint64_t build_reference(size_t count, const Cmpf& c, Cmpf (&z)[2], Scratch& s,
                                ShadowStorage* shadows) {
    complex_zero(z[0]);
    complex_zero(z[1]);
    uint64_t checksum = 0;
    for (size_t i = 0; i < count; ++i) {
        const int previous = static_cast<int>(i & 1);
        const int current = previous ^ 1;
        complex_square_add_to(z[current], z[previous], c, s);
        if (shadows != nullptr) shadows->store(i, z[current]);
        checksum ^= static_cast<uint64_t>(mpf_sgn(z[current].re) + 2) +
                    (static_cast<uint64_t>(mpf_sgn(z[current].im) + 2) << 3) + i;
    }
    g_replay_sink ^= checksum;
    return checksum;
}

static uint64_t replay_period(const ShadowStorage& period, size_t mxit) {
    uint64_t checksum = 0x9e3779b97f4a7c15ULL;
    size_t index = 0;
    for (size_t i = 0; i < mxit; ++i) {
        const std::complex<double>& z = period.value(index);
        checksum ^= static_cast<uint64_t>(std::signbit(z.real())) +
                    (static_cast<uint64_t>(std::signbit(z.imag())) << 1) + i;
        if (period.deep()) {
            const FeShadow& fe = period.fe_value(index);
            checksum ^= static_cast<uint64_t>(fe.re.e) +
                        (static_cast<uint64_t>(fe.im.e) << 1) +
                        (static_cast<uint64_t>(std::signbit(fe.re.m)) << 2) +
                        (static_cast<uint64_t>(std::signbit(fe.im.m)) << 3);
        }
        checksum = (checksum << 7) | (checksum >> (64 - 7));
        if (++index == period.size()) index = 0; // index == i mod p
    }
    g_replay_sink ^= checksum;
    return checksum;
}

static long double estimated_gmp_orbit_bytes() {
    const uint64_t limbs =
        (static_cast<uint64_t>(g_precision_bits) + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS + 1;
    const long double mpf_bytes = static_cast<long double>(sizeof(__mpf_struct)) +
        limbs * static_cast<long double>(sizeof(mp_limb_t));
    return 4.0L * mpf_bytes; // two complex buffers, two mpf_t values each
}

struct BenchmarkTimes {
    double full_shadow_allocate = std::numeric_limits<double>::infinity();
    double full_gmp_only = std::numeric_limits<double>::infinity();
    double full_build_store = std::numeric_limits<double>::infinity();
    double period_shadow_allocate = std::numeric_limits<double>::infinity();
    double period_gmp_only = std::numeric_limits<double>::infinity();
    double period_build_store = std::numeric_limits<double>::infinity();
    double replay = std::numeric_limits<double>::infinity();
    uint64_t checksum = 0;
};

static BenchmarkTimes benchmark_references(const Cmpf& view_center, const Cmpf& nucleus,
                                           uint64_t period, size_t mxit,
                                           int repetitions, bool deep) {
    BenchmarkTimes best;
    Cmpf z[2];
    Scratch s;

    for (int rep = 0; rep < repetitions; ++rep) {
        {
            const double start_allocate = now_seconds();
            ShadowStorage full(mxit, deep);
            const double end_allocate = now_seconds();
            best.checksum ^= build_reference(mxit, view_center, z, s, nullptr);
            const double end_gmp_only = now_seconds();
            best.checksum ^= build_reference(mxit, view_center, z, s, &full);
            const double end_build_store = now_seconds();
            best.full_shadow_allocate =
                std::min(best.full_shadow_allocate, end_allocate - start_allocate);
            best.full_gmp_only =
                std::min(best.full_gmp_only, end_gmp_only - end_allocate);
            best.full_build_store =
                std::min(best.full_build_store, end_build_store - end_gmp_only);
        }

        {
            const double start_allocate = now_seconds();
            ShadowStorage one_period(static_cast<size_t>(period), deep);
            const double end_allocate = now_seconds();
            best.checksum ^= build_reference(
                static_cast<size_t>(period), nucleus, z, s, nullptr);
            const double end_gmp_only = now_seconds();
            best.checksum ^= build_reference(
                static_cast<size_t>(period), nucleus, z, s, &one_period);
            const double end_build_store = now_seconds();
            const uint64_t checksum = replay_period(one_period, mxit);
            const double end_replay = now_seconds();
            best.period_shadow_allocate =
                std::min(best.period_shadow_allocate, end_allocate - start_allocate);
            best.period_gmp_only =
                std::min(best.period_gmp_only, end_gmp_only - end_allocate);
            best.period_build_store =
                std::min(best.period_build_store, end_build_store - end_gmp_only);
            best.replay = std::min(best.replay, end_replay - end_build_store);
            best.checksum ^= checksum;
        }
    }
    return best;
}

static bool parse_u64(const char* text, uint64_t& value) {
    if (text == nullptr || *text == '\0' || *text == '-') return false;
    char* end = nullptr;
    const unsigned long long parsed = std::strtoull(text, &end, 10);
    if (*end != '\0') return false;
    value = static_cast<uint64_t>(parsed);
    return true;
}

static mp_bitcnt_t default_precision(uint64_t scale_exp) {
    const long double bits =
        2.0L * static_cast<long double>(scale_exp) * std::log2(10.0L) + 256.0L;
    const long double limit =
        static_cast<long double>(std::numeric_limits<mp_bitcnt_t>::max());
    if (bits > limit) return std::numeric_limits<mp_bitcnt_t>::max();
    return static_cast<mp_bitcnt_t>(std::ceil(bits));
}

static void print_usage(const char* exe) {
    std::fprintf(stderr,
        "usage:\n"
        "  %s find  cx cy scaleExp [maxPeriod] [precisionBits]\n"
        "  %s bench cx cy scaleExp mxit [maxPeriod] [reps] [precisionBits]\n",
        exe, exe);
}

static void print_scientific(const char* label, mpf_srcptr value) {
    std::printf("%s", label);
    gmp_printf("%.16Fe\n", value);
}

int main(int argc, char** argv) {
    if (argc < 5 || (std::strcmp(argv[1], "find") != 0 &&
                     std::strcmp(argv[1], "bench") != 0)) {
        print_usage(argv[0]);
        return 1;
    }

    const bool bench_mode = std::strcmp(argv[1], "bench") == 0;
    if (bench_mode && argc < 6) {
        print_usage(argv[0]);
        return 1;
    }

    uint64_t scale_exp = 0;
    if (!parse_u64(argv[4], scale_exp)) {
        std::fprintf(stderr, "invalid non-negative scaleExp: %s\n", argv[4]);
        return 1;
    }

    size_t mxit = 0;
    uint64_t max_period = 4000000;
    uint64_t repetitions_u64 = 1;
    uint64_t precision_u64 = default_precision(scale_exp);
    if (bench_mode) {
        uint64_t mxit_u64 = 0;
        if (!parse_u64(argv[5], mxit_u64) || mxit_u64 == 0 ||
            mxit_u64 > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
            std::fprintf(stderr, "invalid mxit: %s\n", argv[5]);
            return 1;
        }
        mxit = static_cast<size_t>(mxit_u64);
        if (argc > 6 && !parse_u64(argv[6], max_period)) return 1;
        if (argc > 7 && !parse_u64(argv[7], repetitions_u64)) return 1;
        if (argc > 8 && !parse_u64(argv[8], precision_u64)) return 1;
    } else {
        if (argc > 5 && !parse_u64(argv[5], max_period)) return 1;
        if (argc > 6 && !parse_u64(argv[6], precision_u64)) return 1;
    }

    if (max_period == 0 || repetitions_u64 == 0 || repetitions_u64 > 100 ||
        precision_u64 < 128 ||
        precision_u64 > static_cast<uint64_t>(std::numeric_limits<mp_bitcnt_t>::max())) {
        std::fprintf(stderr,
            "maxPeriod must be positive, reps must be 1..100, and precisionBits >= 128\n");
        return 1;
    }
    g_precision_bits = static_cast<mp_bitcnt_t>(precision_u64);

    try {
        Cmpf center, nucleus;
        if (mpf_set_str(center.re, argv[2], 10) != 0 ||
            mpf_set_str(center.im, argv[3], 10) != 0) {
            std::fprintf(stderr, "invalid GMP decimal coordinate\n");
            return 1;
        }
        complex_set(nucleus, center);

        Mpf radius, scale;
        if (scale_exp > static_cast<uint64_t>(std::numeric_limits<unsigned long>::max())) {
            std::fprintf(stderr, "scaleExp is too large for this GMP build\n");
            return 1;
        }
        mpf_set_ui(scale.v, 10);
        mpf_pow_ui(scale.v, scale.v, static_cast<unsigned long>(scale_exp));
        mpf_set_ui(radius.v, 2);
        mpf_div(radius.v, radius.v, scale.v);

        std::printf("mode=%s  zoom=1e%llu  precision=%llu bits  maxPeriod=%llu\n",
            bench_mode ? "bench" : "find",
            static_cast<unsigned long long>(scale_exp),
            static_cast<unsigned long long>(g_precision_bits),
            static_cast<unsigned long long>(max_period));
        print_scientific("view radius = ", radius.v);

        const double search_start = now_seconds();
        const SearchResult search = find_period(center, radius.v, max_period);
        const double search_seconds = now_seconds() - search_start;
        if (search.stop != SearchStop::Found) {
            if (search.stop == SearchStop::EscapedDisk) {
                std::printf(
                    "period search: no nucleus detected; the complete propagated "
                    "view disk escaped at iteration %llu (%.6f s)\n",
                    static_cast<unsigned long long>(search.iterations), search_seconds);
            } else {
                std::printf(
                    "period search: no nucleus detected through maxPeriod=%llu "
                    "(%.6f s); a higher-period nucleus is not ruled out\n",
                    static_cast<unsigned long long>(max_period), search_seconds);
            }
            if (bench_mode)
                std::printf("benchmark skipped: periodic storage requires a nucleus in view\n");
            return bench_mode ? 2 : 0;
        }

        std::printf("period search: period=%llu  iterations=%llu  time=%.6f s\n",
            static_cast<unsigned long long>(search.period),
            static_cast<unsigned long long>(search.iterations), search_seconds);

        Mpf residual;
        const double newton_start = now_seconds();
        const NewtonResult newton =
            refine_nucleus(nucleus, search.period, 200, residual.v);
        const double newton_seconds = now_seconds() - newton_start;
        std::printf(
            "Newton: iterations=%llu  converged=%s  target=%llu bits  time=%.6f s\n",
            static_cast<unsigned long long>(newton.iterations),
            newton.converged ? "yes" : "no",
            static_cast<unsigned long long>(newton.convergence_bits), newton_seconds);
        print_scientific("Newton residual |z_p(c)| = ", residual.v);

        Scratch distance_scratch;
        Mpf nucleus_distance;
        complex_distance(nucleus_distance.v, nucleus, center, distance_scratch);
        print_scientific("nucleus distance from view center = ", nucleus_distance.v);
        const bool inside_view = mpf_cmp(nucleus_distance.v, radius.v) <= 0;
        std::printf("nucleus inside searched disk = %s\n", inside_view ? "yes" : "no");
        gmp_printf("nucleus re = %.40Fe\n", nucleus.re);
        gmp_printf("nucleus im = %.40Fe\n", nucleus.im);

        if (!bench_mode) return newton.converged && inside_view ? 0 : 2;
        if (!newton.converged || newton.singular || !inside_view) {
            std::printf("benchmark skipped: Newton did not produce a usable in-view nucleus\n");
            return 2;
        }
        if (search.period > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
            std::fprintf(stderr, "period does not fit this platform's size_t\n");
            return 2;
        }

        const bool deep_shadows = scale_exp > 280;
        const long double shallow_sample_bytes = sizeof(std::complex<double>);
        const long double deep_sample_bytes =
            sizeof(std::complex<double>) + sizeof(FeShadow);
        const long double full_shallow_bytes =
            shallow_sample_bytes * static_cast<long double>(mxit);
        const long double period_shallow_bytes =
            shallow_sample_bytes * static_cast<long double>(search.period);
        const long double full_deep_bytes =
            deep_sample_bytes * static_cast<long double>(mxit);
        const long double period_deep_bytes =
            deep_sample_bytes * static_cast<long double>(search.period);
        const long double mib = 1024.0L * 1024.0L;
        std::printf("\nreference-orbit microbenchmark: mxit=%llu  period=%llu  reps=%llu\n",
            static_cast<unsigned long long>(mxit),
            static_cast<unsigned long long>(search.period),
            static_cast<unsigned long long>(repetitions_u64));
        std::printf(
            "constant GMP orbit state: two rotating complex buffers (4 mpf_t), "
            "estimated %.3Lf KiB\n"
            "shadow payload estimates (vector metadata excluded):\n"
            "  shallow complex<double> (%zu B/sample): full %.6Lf MiB, period %.6Lf MiB\n"
            "  deep + two FloatExp      (%zu B/sample): full %.6Lf MiB, period %.6Lf MiB\n"
            "  active shadow mode: %s\n"
            "  shadow storage ratio: %.2Lfx\n",
            estimated_gmp_orbit_bytes() / 1024.0L,
            sizeof(std::complex<double>), full_shallow_bytes / mib,
            period_shallow_bytes / mib,
            sizeof(std::complex<double>) + sizeof(FeShadow), full_deep_bytes / mib,
            period_deep_bytes / mib, deep_shadows ? "deep" : "shallow",
            static_cast<long double>(mxit) / static_cast<long double>(search.period));

        const BenchmarkTimes timing = benchmark_references(
            center, nucleus, search.period, mxit, static_cast<int>(repetitions_u64),
            deep_shadows);
        std::printf(
            "best timings:\n"
            "  full shadow allocation/init : %.6f s\n"
            "  full GMP arithmetic only    : %.6f s\n"
            "  full GMP build+shadow write : %.6f s\n"
            "  period shadow allocation    : %.6f s\n"
            "  period GMP arithmetic only  : %.6f s\n"
            "  period GMP build+shadow write: %.6f s\n"
            "  shadow modulo replay        : %.6f s\n"
            "  GMP arithmetic ratio        : %.2fx\n"
            "  build+shadow-write ratio    : %.2fx\n",
            timing.full_shadow_allocate, timing.full_gmp_only, timing.full_build_store,
            timing.period_shadow_allocate, timing.period_gmp_only,
            timing.period_build_store, timing.replay,
            timing.period_gmp_only == 0.0 ? 0.0 :
                timing.full_gmp_only / timing.period_gmp_only,
            timing.period_build_store == 0.0 ? 0.0 :
                timing.full_build_store / timing.period_build_store);
        std::printf("  checksum            : %llu\n",
            static_cast<unsigned long long>(timing.checksum));

        std::printf(
            "\nCAVEAT: this compares a forced-mxit full reference at the view center "
            "with one GMP-precision period at the refined nucleus. The residual above "
            "is the finite-precision closure error, not an algebraic certificate. "
            "Allocation, GMP arithmetic, shadow population, and replay are deliberately "
            "reported separately; no aggregate speedup is claimed. The timing excludes "
            "period search, Newton, perturbation, series approximation, BLA, rebasing, "
            "pixel work, and renderer correctness. Modulo replay alone is not an "
            "end-to-end renderer speedup.\n");
    } catch (const std::bad_alloc&) {
        std::fprintf(stderr, "allocation failed; lower mxit/period or precision\n");
        return 2;
    }
    return 0;
}
