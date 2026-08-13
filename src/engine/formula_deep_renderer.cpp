#include "formula_deep_renderer.h"
#include "bigfixed.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <mutex>
#include <new>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace formula {

namespace {

using Clock = std::chrono::steady_clock;
constexpr mpfr_rnd_t RND = MPFR_RNDN;
constexpr uint8_t NoFallbackReason = static_cast<uint8_t>(ExpressionDeepFallbackReason::Count);

double secondsSince(const Clock::time_point& start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

bool finiteComplex(Complex value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

bool sameMpfrStateComponent(mpfr_srcptr left, mpfr_srcptr right) {
    return mpfr_equal_p(left, right) && (!mpfr_zero_p(left) || mpfr_signbit(left) == mpfr_signbit(right));
}

bool decimalSignificandIsNonzero(const std::string& text) {
    size_t end = text.find_first_of("eE");
    if (end == std::string::npos) end = text.size();
    for (size_t i = 0; i < end; ++i)
        if (text[i] >= '1' && text[i] <= '9') return true;
    return false;
}

size_t decimalDigits(const std::string& text) {
    size_t digits = 0;
    size_t end = text.find_first_of("eE");
    if (end == std::string::npos) end = text.size();
    for (size_t i = 0; i < end; ++i)
        if (text[i] >= '0' && text[i] <= '9') ++digits;
    return digits;
}

uint64_t decimalBits(const std::string& text) {
    return static_cast<uint64_t>(std::ceil(decimalDigits(text) * 3.32192809488736234787)) + 2;
}

uint64_t mpfBits(mpf_srcptr value) {
    if (!value) return 0;
    const mp_size_t limbCount = mpf_size(value);
    if (limbCount == 0) return 0;
    const size_t highBit = mpn_sizeinbase(value->_mp_d, limbCount, 2);
    const mp_bitcnt_t lowBit = mpn_scan1(value->_mp_d, 0);
    return highBit > lowBit ? static_cast<uint64_t>(highBit - lowBit) : 1;
}

bool validCenterRepresentation(const ExpressionReferenceExactInput& center) {
    if (center.usesMpf()) return center.realMpf && center.imaginaryMpf && !center.usesDecimal();
    return !center.realDecimal.empty() && !center.usesMpf();
}

bool validScaleRepresentation(const ExpressionDeepExactReal& scale) {
    return scale.usesMpf() ? !scale.usesDecimal() : scale.usesDecimal();
}

bool setDecimal(mpfr_ptr output, const std::string& text) {
    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);
    const int status = mpfr_set_str(output, text.c_str(), 10, RND);
    const mpfr_flags_t rangeFlags = mpfr_flags_test(MPFR_FLAGS_UNDERFLOW | MPFR_FLAGS_OVERFLOW | MPFR_FLAGS_ERANGE);
    const bool okay = status == 0 && rangeFlags == 0 && mpfr_number_p(output) && !(mpfr_zero_p(output) && decimalSignificandIsNonzero(text));
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    return okay;
}

bool setExactMpf(mpfr_ptr output, mpf_srcptr input) {
    if (!input) return false;
    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);
    mpfr_set_f(output, input, RND);
    const mpfr_flags_t rangeFlags = mpfr_flags_test(MPFR_FLAGS_UNDERFLOW | MPFR_FLAGS_OVERFLOW | MPFR_FLAGS_ERANGE);
    const bool okay = rangeFlags == 0 && mpfr_number_p(output) && !(mpf_sgn(input) != 0 && mpfr_zero_p(output));
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    return okay;
}

bool setExactReal(mpfr_ptr output, const ExpressionDeepExactReal& input) {
    if (input.usesMpf()) {
        if (input.usesDecimal()) return false;
        return setExactMpf(output, input.mpf);
    }
    return input.usesDecimal() && setDecimal(output, input.decimal);
}

bool setExactCenter(MpfrComplex& output, const ExpressionReferenceExactInput& input) {
    if (input.usesMpf()) {
        if (!input.realMpf || !input.imaginaryMpf || input.usesDecimal()) return false;
        if (!setExactMpf(output.re, input.realMpf) || !setExactMpf(output.im, input.imaginaryMpf)) return false;
    } else {
        if (input.realDecimal.empty()) return false;
        const std::string imaginary = input.imaginaryDecimal.empty() ? "0" : input.imaginaryDecimal;
        if (!setDecimal(output.re, input.realDecimal) || !setDecimal(output.im, imaginary)) return false;
    }
    return mpfr_number_p(output.re) && mpfr_number_p(output.im);
}

class ExactGeometry {
  public:
    explicit ExactGeometry(mpfr_prec_t precision) : center(precision) { mpfr_inits2(precision, scale, dxHalf, dyHalf, temporary, (mpfr_ptr)0); }

    ~ExactGeometry() { mpfr_clears(scale, dxHalf, dyHalf, temporary, (mpfr_ptr)0); }

    bool initialize(const ExpressionDeepRenderRequest& request) {
        if (!setExactCenter(center, request.center) || !setExactReal(scale, request.scale) || mpfr_sgn(scale) <= 0) return false;

        // dx/2 = 2/(scale*(width-1)).
        mpfr_mul_ui(temporary, scale, static_cast<unsigned long>(request.width - 1), RND);
        mpfr_ui_div(dxHalf, 2, temporary, RND);

        // dy/2 = 2*height/(scale*width*(height-1)).
        mpfr_mul_ui(temporary, scale, static_cast<unsigned long>(request.width), RND);
        mpfr_mul_ui(temporary, temporary, static_cast<unsigned long>(request.height - 1), RND);
        mpfr_ui_div(dyHalf, static_cast<unsigned long>(request.height), temporary, RND);
        mpfr_mul_ui(dyHalf, dyHalf, 2, RND);
        return mpfr_number_p(dxHalf) && mpfr_number_p(dyHalf) && !mpfr_zero_p(dxHalf) && !mpfr_zero_p(dyHalf);
    }

    bool coordinate(int x, int y, const ExpressionDeepRenderRequest& request, MpfrComplex& output) {
        const long centeredX = static_cast<long>(2LL * x - (request.width - 1LL));
        const long centeredY = static_cast<long>(2LL * y - (request.height - 1LL));
        if (centeredX == 0) {
            mpfr_set(output.re, center.re, RND);
        } else {
            mpfr_mul_si(output.re, dxHalf, centeredX, RND);
            mpfr_add(output.re, output.re, center.re, RND);
        }
        if (centeredY == 0) {
            mpfr_set(output.im, center.im, RND);
        } else {
            mpfr_mul_si(output.im, dyHalf, centeredY, RND);
            mpfr_add(output.im, output.im, center.im, RND);
        }
        return mpfr_number_p(output.re) && mpfr_number_p(output.im);
    }

    MpfrComplex center;
    mpfr_t scale;
    mpfr_t dxHalf;
    mpfr_t dyHalf;
    mpfr_t temporary;
};

struct FallbackWorkerWorkspace {
    FallbackWorkerWorkspace(mpfr_prec_t precision, const ExpressionDeepRenderRequest& request) : geometry(precision), geometryReady(geometry.initialize(request)), context(precision), pixel(precision), next(precision), periodicState(precision), piecewiseSquares(precision), magnitudeStorage(precision) {
        for (mpfr_ptr value : piecewiseScratch) mpfr_init2(value, precision);
        mpfr_init2(piecewiseBailoutSquared, precision);
        mpfr_set_d(piecewiseBailoutSquared, request.bailout, RND);
        mpfr_sqr(piecewiseBailoutSquared, piecewiseBailoutSquared, MPFR_RNDD);
        piecewiseInsideSquareExponent = mpfr_get_exp(piecewiseBailoutSquared) - 2;
    }

    ~FallbackWorkerWorkspace() {
        for (mpfr_ptr value : piecewiseScratch) mpfr_clear(value);
        mpfr_clear(piecewiseBailoutSquared);
    }

    ExactGeometry geometry;
    bool geometryReady = false;
    ExpressionOracleContext context;
    MpfrComplex pixel;
    MpfrComplex next;
    MpfrComplex periodicState;
    MpfrComplex piecewiseSquares;
    MpfrComplex magnitudeStorage;
    mpfr_t piecewiseScratch[5];
    mpfr_t piecewiseBailoutSquared;
    mpfr_exp_t piecewiseInsideSquareExponent = 0;
};

struct BigFixedPixelWorkspace {
    BigFixedPixelWorkspace(int limbCount, double bailout, mpfr_prec_t fallbackPrecision)
        : limbs(limbCount), mpfrPrecision(fallbackPrecision), fractionalBits(64LL * (limbs - 1)), bailoutSquared(bailout * bailout), real(limbs), imaginary(limbs), constantReal(limbs), constantImaginary(limbs), realSquared(limbs), imaginarySquared(limbs), product(limbs), realPart(limbs), imaginaryPart(limbs), nextReal(limbs), nextImaginary(limbs), scratch(static_cast<size_t>(2 * limbs)), limbScales(static_cast<size_t>(limbs)) {
        mpf_init2(bridge, static_cast<mp_bitcnt_t>(64 * (limbs + 1)));
        for (int limb = 0; limb < limbs; ++limb) limbScales[static_cast<size_t>(limb)] = std::ldexp(1.0, 64 * (limb - (limbs - 1)));
    }

    ~BigFixedPixelWorkspace() { mpf_clear(bridge); }

    bool set(BigFixed& output, mpfr_srcptr value) {
        if (!mpfr_number_p(value) || mpfr_cmpabs_ui(value, 1UL << 30) > 0) return false;
        mpfr_get_f(bridge, value, MPFR_RNDN);
        bf_from_mpf(output, bridge, limbs);
        if (!mpfr_zero_p(value) && std::none_of(output.m.begin(), output.m.end(), [](uint64_t limb) { return limb != 0; })) return false;
        return true;
    }

    int limbs;
    mpfr_prec_t mpfrPrecision;
    int64_t fractionalBits;
    double bailoutSquared;
    BigFixed real, imaginary, constantReal, constantImaginary;
    BigFixed realSquared, imaginarySquared, product, realPart, imaginaryPart;
    BigFixed nextReal, nextImaginary;
    std::vector<uint64_t> scratch;
    std::vector<double> limbScales;
    mpf_t bridge;
    double fixedUnit = 0.0;
    double orbitError = 0.0;
    double realMagnitude = 0.0;
    double imaginaryMagnitude = 0.0;
    double constantRealMagnitude = 0.0;
    double constantImaginaryMagnitude = 0.0;
};

enum class BigFixedPixelStatus : uint8_t {
    Success,
    Cancelled,
    Unavailable
};

enum class BigFixedBailoutDecision : uint8_t {
    Inside,
    Outside,
    Ambiguous
};

double inflatePositiveDouble(double value, uint64_t ulps = 1) {
    if (value == 0.0 || !std::isfinite(value)) return value;
    uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    constexpr uint64_t InfinityBits = 0x7ff0000000000000ULL;
    bits = bits >= InfinityBits - ulps ? InfinityBits : bits + ulps;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

double deflatePositiveDouble(double value, uint64_t ulps = 1) {
    if (value == 0.0) return 0.0;
    if (!std::isfinite(value)) return value;
    uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    bits = bits > ulps ? bits - ulps : 0;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

double positiveBigFixedUpper(const BigFixedPixelWorkspace& workspace, const BigFixed& value) {
    int top = workspace.limbs - 1;
    while (top >= 0 && value.m[static_cast<size_t>(top)] == 0) --top;
    if (top < 0) return 0.0;
    double magnitude = static_cast<double>(value.m[static_cast<size_t>(top)]) * workspace.limbScales[static_cast<size_t>(top)];
    if (top > 0) magnitude += static_cast<double>(value.m[static_cast<size_t>(top - 1)]) * workspace.limbScales[static_cast<size_t>(top - 1)];
    return inflatePositiveDouble(magnitude, 8);
}

void updateBigFixedMagnitudes(BigFixedPixelWorkspace& workspace) {
    workspace.realMagnitude = positiveBigFixedUpper(workspace, workspace.real);
    workspace.imaginaryMagnitude = positiveBigFixedUpper(workspace, workspace.imaginary);
}

void initializeBigFixedErrorBounds(BigFixedPixelWorkspace& workspace) {
    workspace.fixedUnit = std::ldexp(1.0, -static_cast<int>(workspace.fractionalBits));
    workspace.orbitError = workspace.fixedUnit;
    workspace.constantRealMagnitude = inflatePositiveDouble(positiveBigFixedUpper(workspace, workspace.constantReal) + workspace.fixedUnit, 8);
    workspace.constantImaginaryMagnitude = inflatePositiveDouble(positiveBigFixedUpper(workspace, workspace.constantImaginary) + workspace.fixedUnit, 8);
    updateBigFixedMagnitudes(workspace);
}

void advanceBigFixedErrorBounds(BigFixedPixelWorkspace& workspace) {
    const double realWithError = workspace.realMagnitude + workspace.orbitError;
    const double imaginaryWithError = workspace.imaginaryMagnitude + workspace.orbitError;
    const double fixedLocalError = 2.0 * workspace.fixedUnit;
    const double realOperationMagnitude = realWithError * realWithError + imaginaryWithError * imaginaryWithError + workspace.constantRealMagnitude;
    const double imaginaryOperationMagnitude = 2.0 * realWithError * imaginaryWithError + workspace.constantImaginaryMagnitude;
    const double propagation = 2.0 * (workspace.realMagnitude + workspace.imaginaryMagnitude) * workspace.orbitError + 2.0 * workspace.orbitError * workspace.orbitError;
    const double operationMagnitude = std::max(realOperationMagnitude, imaginaryOperationMagnitude);
    const double nextError = propagation + fixedLocalError + std::ldexp(operationMagnitude, 4 - static_cast<int>(workspace.mpfrPrecision));

    // All terms are nonnegative. Sixty-four ulps cover the rounding of the
    // short scalar expressions above while remaining negligible beside the
    // deliberately conservative MPFR operation allowance.
    workspace.orbitError = inflatePositiveDouble(nextError, 64);
}

bool bigFixedComponentInterval(double magnitude, double error, double& lower, double& upper) {
    if (!std::isfinite(magnitude)) return false;
    const double upwardSpacing = std::nextafter(magnitude, std::numeric_limits<double>::infinity()) - magnitude;
    const double downwardSpacing = magnitude - std::nextafter(magnitude, -std::numeric_limits<double>::infinity());
    const double conversionError = 16.0 * std::max(upwardSpacing, downwardSpacing);
    const double totalError = inflatePositiveDouble(error + conversionError, 8);
    if (!std::isfinite(totalError)) return false;
    lower = deflatePositiveDouble(std::max(0.0, magnitude - totalError), 8);
    upper = inflatePositiveDouble(magnitude + totalError, 8);
    return std::isfinite(upper);
}

BigFixedBailoutDecision decideBigFixedBailout(const BigFixedPixelWorkspace& workspace) {
    double realLower = 0.0;
    double realUpper = 0.0;
    double imaginaryLower = 0.0;
    double imaginaryUpper = 0.0;
    if (!bigFixedComponentInterval(workspace.realMagnitude, workspace.orbitError, realLower, realUpper) || !bigFixedComponentInterval(workspace.imaginaryMagnitude, workspace.orbitError, imaginaryLower, imaginaryUpper)) return BigFixedBailoutDecision::Ambiguous;

    const double lowerRealSquared = deflatePositiveDouble(realLower * realLower, 8);
    const double lowerImaginarySquared = deflatePositiveDouble(imaginaryLower * imaginaryLower, 8);
    const double lowerMagnitudeSquared = deflatePositiveDouble(lowerRealSquared + lowerImaginarySquared, 8);
    if (lowerMagnitudeSquared > workspace.bailoutSquared) return BigFixedBailoutDecision::Outside;

    const double upperRealSquared = inflatePositiveDouble(realUpper * realUpper, 8);
    const double upperImaginarySquared = inflatePositiveDouble(imaginaryUpper * imaginaryUpper, 8);
    const double upperMagnitudeSquared = inflatePositiveDouble(upperRealSquared + upperImaginarySquared, 8);
    if (upperMagnitudeSquared <= workspace.bailoutSquared) return BigFixedBailoutDecision::Inside;
    return BigFixedBailoutDecision::Ambiguous;
}

template <typename Cancel>
BigFixedPixelStatus evaluateBurningShipBigFixedPixel(const ExpressionDeepRenderRequest& request, const ExpressionOracleContext& context, BigFixedPixelWorkspace& workspace, float& output, uint64_t& iterations, Cancel&& shouldCancel) {
    if (!workspace.set(workspace.real, context.z.re) || !workspace.set(workspace.imaginary, context.z.im) || !workspace.set(workspace.constantReal, context.c.re) || !workspace.set(workspace.constantImaginary, context.c.im)) return BigFixedPixelStatus::Unavailable;
    initializeBigFixedErrorBounds(workspace);
    bf_sqr(workspace.realSquared, workspace.real, workspace.scratch.data());
    bf_sqr(workspace.imaginarySquared, workspace.imaginary, workspace.scratch.data());
    const BigFixedBailoutDecision initialDecision = decideBigFixedBailout(workspace);
    if (initialDecision == BigFixedBailoutDecision::Ambiguous) return BigFixedPixelStatus::Unavailable;
    if (initialDecision == BigFixedBailoutDecision::Outside) {
        output = 0.0f;
        return BigFixedPixelStatus::Success;
    }
    output = ExpressionDeepInteriorPixel;
    for (int iteration = 0; iteration < request.maxIterations; ++iteration) {
        if ((iteration & 15) == 0 && shouldCancel()) return BigFixedPixelStatus::Cancelled;
        advanceBigFixedErrorBounds(workspace);
        bf_sub(workspace.realPart, workspace.realSquared, workspace.imaginarySquared);
        bf_mul(workspace.product, workspace.real, workspace.imaginary, workspace.scratch.data());
        if (!workspace.product.isZero()) workspace.product.sign = 1;
        bf_add(workspace.imaginaryPart, workspace.product, workspace.product);
        bf_add(workspace.nextReal, workspace.realPart, workspace.constantReal);
        bf_add(workspace.nextImaginary, workspace.imaginaryPart, workspace.constantImaginary);
        bf_sqr(workspace.realSquared, workspace.nextReal, workspace.scratch.data());
        bf_sqr(workspace.imaginarySquared, workspace.nextImaginary, workspace.scratch.data());
        workspace.real.m.swap(workspace.nextReal.m);
        workspace.real.sign = workspace.nextReal.sign;
        workspace.imaginary.m.swap(workspace.nextImaginary.m);
        workspace.imaginary.sign = workspace.nextImaginary.sign;
        updateBigFixedMagnitudes(workspace);
        ++iterations;
        const BigFixedBailoutDecision decision = decideBigFixedBailout(workspace);
        if (decision == BigFixedBailoutDecision::Ambiguous) return BigFixedPixelStatus::Unavailable;
        if (decision == BigFixedBailoutDecision::Outside) {
            output = static_cast<float>(iteration + 1);
            return BigFixedPixelStatus::Success;
        }
    }
    return BigFixedPixelStatus::Success;
}

bool evaluatePiecewiseQuadraticMpfr(ExpressionPiecewiseQuadraticKind kind, const ExpressionOracleContext& context, const MpfrComplex& componentSquares, MpfrComplex& output, mpfr_t (&scratch)[5]) {
    if (kind != ExpressionPiecewiseQuadraticKind::BurningShip) return false;
    // The generic oracle computes hypot(component, +0) after real()/imag().
    // Its exact mathematical result is representable, so abs is identical and
    // avoids two general square-root calls.
    mpfr_sub(output.re, componentSquares.re, componentSquares.im, RND);
    mpfr_mul(scratch[3], context.z.re, context.z.im, RND);
    mpfr_abs(scratch[3], scratch[3], RND);
    // The generic complex square adds two identical rounded products.
    mpfr_mul_2ui(output.im, scratch[3], 1, RND);
    mpfr_add(output.re, output.re, context.c.re, RND);
    mpfr_add(output.im, output.im, context.c.im, RND);
    return true;
}

bool piecewiseOutsideBailout(const MpfrComplex& value, mpfr_srcptr bailoutSquared, double bailout, mpfr_ptr magnitude, MpfrComplex& componentSquares, mpfr_exp_t insideSquareExponent, mpfr_t (&scratch)[5], bool& outside) {
    const int reRounded = mpfr_sqr(componentSquares.re, value.re, RND);
    const int imRounded = mpfr_sqr(componentSquares.im, value.im, RND);
    if (mpfr_nan_p(componentSquares.re) || mpfr_nan_p(componentSquares.im)) return false;
    if (mpfr_inf_p(componentSquares.re) || mpfr_inf_p(componentSquares.im)) {
        outside = true;
        return true;
    }
    const bool realClearlyInside = mpfr_zero_p(componentSquares.re) || mpfr_get_exp(componentSquares.re) <= insideSquareExponent;
    const bool imaginaryClearlyInside = mpfr_zero_p(componentSquares.im) || mpfr_get_exp(componentSquares.im) <= insideSquareExponent;
    if (realClearlyInside && imaginaryClearlyInside) {
        outside = false;
        return true;
    }
    mpfr_add(scratch[2], componentSquares.re, componentSquares.im, MPFR_RNDU);
    // Each nearest-rounded square can undershoot by at most half an
    // ulp. Since the nonnegative sum has an ulp at least as large as
    // either operand, one upward step covers both possible errors.
    if (reRounded < 0 || imRounded < 0) mpfr_nextabove(scratch[2]);
    if (mpfr_cmp(scratch[2], bailoutSquared) <= 0) {
        outside = false;
        return true;
    }
    mpfr_hypot(magnitude, value.re, value.im, RND);
    if (mpfr_nan_p(magnitude)) return false;
    outside = mpfr_inf_p(magnitude) || mpfr_cmp_d(magnitude, bailout) > 0;
    return true;
}

bool checkedAddSize(size_t& total, size_t count, size_t bytes) {
    if (bytes != 0 && count > (std::numeric_limits<size_t>::max() - total) / bytes) return false;
    total += count * bytes;
    return true;
}

uint64_t ceilLog2(uint64_t value) {
    uint64_t result = 0;
    uint64_t power = 1;
    while (power < value && result < 63) {
        power <<= 1;
        ++result;
    }
    return result;
}

bool selectAutomaticViewBits(const ExpressionDeepRenderRequest& request, ExpressionReferencePrecisionPolicy& policy, std::string& error) {
    const uint64_t centerBits = request.center.usesMpf() ? std::max(mpfBits(request.center.realMpf), mpfBits(request.center.imaginaryMpf)) : std::max(decimalBits(request.center.realDecimal), decimalBits(request.center.imaginaryDecimal));
    const uint64_t scaleBits = request.scale.usesMpf() ? mpfBits(request.scale.mpf) : decimalBits(request.scale.decimal);
    const uint64_t inputBits = std::max(centerBits, scaleBits);
    if (inputBits > static_cast<uint64_t>(MPFR_PREC_MAX)) {
        error = "exact input precision exceeds MPFR";
        return false;
    }
    policy.requestedBits = std::max(policy.requestedBits, static_cast<mpfr_prec_t>(inputBits));

    const mpfr_prec_t probePrecision = 256;
    ExactGeometry geometry(probePrecision);
    if (!geometry.initialize(request)) {
        error = "failed to parse exact center or positive scale";
        return false;
    }

    uint64_t required = 53;
    const mpfr_exp_t scaleExponent = mpfr_get_exp(geometry.scale);
    if (scaleExponent > 0) { required = std::max<uint64_t>(required, static_cast<uint64_t>(scaleExponent) + ceilLog2(static_cast<uint64_t>(std::max(request.width, request.height))) + 8); }
    auto coverAddition = [&](mpfr_srcptr center, mpfr_srcptr step) {
        if (mpfr_zero_p(center) || mpfr_zero_p(step)) return;
        const mpfr_exp_t difference = mpfr_get_exp(center) - mpfr_get_exp(step);
        if (difference > 0) required = std::max<uint64_t>(required, static_cast<uint64_t>(difference) + 8);
    };
    coverAddition(geometry.center.re, geometry.dxHalf);
    coverAddition(geometry.center.im, geometry.dyHalf);
    if (required > static_cast<uint64_t>(MPFR_PREC_MAX)) {
        error = "view precision exceeds MPFR";
        return false;
    }
    policy.viewBits = std::max(policy.viewBits, static_cast<mpfr_prec_t>(required));
    return true;
}

bool selectDeepPrecision(const ExpressionReferencePrecisionPolicy& policy, mpfr_prec_t& precision, std::string& error) {
    const uint64_t base = std::max<uint64_t>({53, static_cast<uint64_t>(policy.requestedBits), static_cast<uint64_t>(policy.viewBits), static_cast<uint64_t>(policy.minimumBits)});
    const uint64_t guard = static_cast<uint64_t>(policy.guardBits);
    if (base > std::numeric_limits<uint64_t>::max() - guard) {
        error = "MPFR precision calculation overflow";
        return false;
    }
    const uint64_t selected = base + guard;
    if (selected > static_cast<uint64_t>(policy.maximumBits) || selected > static_cast<uint64_t>(ExpressionReferencePrecisionPolicy::ApplicationMaximumBits) || selected > static_cast<uint64_t>(MPFR_PREC_MAX)) {
        error = "required MPFR precision exceeds policy maximum";
        return false;
    }
    precision = static_cast<mpfr_prec_t>(selected);
    return true;
}

void configureFixed(const ExpressionContext& fixed, ExpressionOracleContext& context) {
    context.z.set(fixed.z.real(), fixed.z.imag());
    context.c.set(fixed.c.real(), fixed.c.imag());
    context.z0.set(fixed.z0.real(), fixed.z0.imag());
    for (size_t i = 0; i < fixed.parameters.size(); ++i) { context.parameters[i].set(fixed.parameters[i].real(), fixed.parameters[i].imag()); }
    context.iteration = fixed.iteration;
}

ExpressionDeepFallbackReason reasonForCapability(ExpressionScaledResidualCapability capability) {
    switch (capability) {
    case ExpressionScaledResidualCapability::CertifiedEntireCandidate:
    case ExpressionScaledResidualCapability::CertifiedMeromorphicCandidate:
    case ExpressionScaledResidualCapability::CertifiedBranchCandidate:
    case ExpressionScaledResidualCapability::CertifiedRealCandidate:
    case ExpressionScaledResidualCapability::CertifiedPiecewiseCandidate: return ExpressionDeepFallbackReason::CertificationFailure;
    case ExpressionScaledResidualCapability::UncertifiedSeries: return ExpressionDeepFallbackReason::UncertifiedSeries;
    case ExpressionScaledResidualCapability::BranchSensitive: return ExpressionDeepFallbackReason::BranchSensitive;
    case ExpressionScaledResidualCapability::Unsupported: return ExpressionDeepFallbackReason::UnsupportedOperation;
    case ExpressionScaledResidualCapability::ExactCenteredArithmetic: break;
    }
    return ExpressionDeepFallbackReason::InvalidTape;
}

bool certifiedTaylorCapability(ExpressionScaledResidualCapability capability) {
    return capability == ExpressionScaledResidualCapability::CertifiedEntireCandidate || capability == ExpressionScaledResidualCapability::CertifiedMeromorphicCandidate || capability == ExpressionScaledResidualCapability::CertifiedBranchCandidate || capability == ExpressionScaledResidualCapability::CertifiedRealCandidate || capability == ExpressionScaledResidualCapability::CertifiedPiecewiseCandidate;
}

bool certifiedReferenceCapability(ExpressionScaledResidualCapability capability) {
    return capability == ExpressionScaledResidualCapability::ExactCenteredArithmetic || certifiedTaylorCapability(capability);
}

ExpressionDeepFallbackReason reasonForResidualStatus(ExpressionScaledResidualStatus status) {
    switch (status) {
    case ExpressionScaledResidualStatus::BranchUncertain: return ExpressionDeepFallbackReason::BranchSensitive;
    case ExpressionScaledResidualStatus::Singular: return ExpressionDeepFallbackReason::Singular;
    case ExpressionScaledResidualStatus::Unsupported: return ExpressionDeepFallbackReason::UnsupportedOperation;
    case ExpressionScaledResidualStatus::Nonfinite: return ExpressionDeepFallbackReason::Nonfinite;
    case ExpressionScaledResidualStatus::ExponentRange: return ExpressionDeepFallbackReason::ExponentRange;
    case ExpressionScaledResidualStatus::InvalidTape:
    case ExpressionScaledResidualStatus::InvalidInput: return ExpressionDeepFallbackReason::InvalidTape;
    case ExpressionScaledResidualStatus::Success: break;
    }
    return ExpressionDeepFallbackReason::InvalidTape;
}

ExpressionDeepFallbackReason reasonForArithmeticStatus(ScaledArithmeticStatus status) {
    switch (status) {
    case ScaledArithmeticStatus::Nonfinite: return ExpressionDeepFallbackReason::Nonfinite;
    case ScaledArithmeticStatus::ExponentRange: return ExpressionDeepFallbackReason::ExponentRange;
    case ScaledArithmeticStatus::Singular: return ExpressionDeepFallbackReason::Singular;
    case ScaledArithmeticStatus::Success: break;
    }
    return ExpressionDeepFallbackReason::ReconstructionFailure;
}

enum class BailoutDecision : uint8_t { Inside,
                                       Outside,
                                       Uncertain,
                                       Error };

struct BailoutThreshold {
    ScaledRealValue midpoint;
    ScaledRealValue error;
};

bool makeBailoutThreshold(double bailout, BailoutThreshold& threshold) {
    ScaledComplexValue radius;
    if (makeScaledRealValue(bailout, radius.re) != ScaledArithmeticStatus::Success || certifiedScaledNormSquared(radius, {}, threshold.midpoint, threshold.error) != ScaledArithmeticStatus::Success) return false;
    return true;
}

BailoutDecision decideBailout(const ScaledComplexValue& value, const ScaledRealValue& radius, const BailoutThreshold& threshold, ScaledArithmeticStatus& arithmeticStatus) {
    ScaledComplexBall state{value, radius};
    arithmeticStatus = certifyScaledMpfrExponentRange(state);
    if (arithmeticStatus != ScaledArithmeticStatus::Success) return BailoutDecision::Error;
    ScaledRealValue norm, error;
    arithmeticStatus = certifiedScaledNormSquared(value, radius, norm, error);
    if (arithmeticStatus != ScaledArithmeticStatus::Success) return BailoutDecision::Error;
    arithmeticStatus = certifyScaledMpfrExponentRange(norm);
    if (arithmeticStatus != ScaledArithmeticStatus::Success) return BailoutDecision::Error;
    arithmeticStatus = certifyScaledMpfrExponentRange(error);
    if (arithmeticStatus != ScaledArithmeticStatus::Success) return BailoutDecision::Error;
    ScaledRealValue thresholdUpper, outsideMargin;
    arithmeticStatus = scaledAddUp(threshold.midpoint, threshold.error, thresholdUpper);
    if (arithmeticStatus != ScaledArithmeticStatus::Success) return BailoutDecision::Error;
    arithmeticStatus = certifyScaledMpfrExponentRange(thresholdUpper);
    if (arithmeticStatus != ScaledArithmeticStatus::Success) return BailoutDecision::Error;
    arithmeticStatus = scaledAddUp(thresholdUpper, error, outsideMargin);
    if (arithmeticStatus != ScaledArithmeticStatus::Success) return BailoutDecision::Error;
    arithmeticStatus = certifyScaledMpfrExponentRange(outsideMargin);
    if (arithmeticStatus != ScaledArithmeticStatus::Success) return BailoutDecision::Error;
    if (compareScaledNonnegative(norm, outsideMargin) > 0) return BailoutDecision::Outside;
    ScaledRealValue normUpper, insideMargin;
    arithmeticStatus = scaledAddUp(norm, error, normUpper);
    if (arithmeticStatus != ScaledArithmeticStatus::Success) return BailoutDecision::Error;
    arithmeticStatus = certifyScaledMpfrExponentRange(normUpper);
    if (arithmeticStatus != ScaledArithmeticStatus::Success) return BailoutDecision::Error;
    arithmeticStatus = scaledAddUp(normUpper, threshold.error, insideMargin);
    if (arithmeticStatus != ScaledArithmeticStatus::Success) return BailoutDecision::Error;
    arithmeticStatus = certifyScaledMpfrExponentRange(insideMargin);
    if (arithmeticStatus != ScaledArithmeticStatus::Success) return BailoutDecision::Error;
    if (compareScaledNonnegative(insideMargin, threshold.midpoint) <= 0) return BailoutDecision::Inside;
    return BailoutDecision::Uncertain;
}

struct ScaledOffset {
    ScaledRealValue value;
    ScaledRealValue error;
};

struct TaylorTileJet {
    int xBegin = 0;
    int yBegin = 0;
    int xEnd = 0;
    int yEnd = 0;
    int depth = 0;
    ExpressionTaylorJetResult jet;
};

struct PendingTaylorTile {
    int xBegin = 0;
    int yBegin = 0;
    int xEnd = 0;
    int yEnd = 0;
    int depth = 0;
};

ScaledRealValue absoluteScaled(ScaledRealValue value) {
    value.mantissa = std::abs(value.mantissa);
    return value;
}

uint64_t saturatingAdd(uint64_t left, uint64_t right) {
    return right > std::numeric_limits<uint64_t>::max() - left ? std::numeric_limits<uint64_t>::max() : left + right;
}

void mergeMaximum(ScaledRealValue& target, const ScaledRealValue& value) {
    if (compareScaledNonnegative(value, target) > 0) target = value;
}

void mergeMinimumNonzero(ScaledRealValue& target, const ScaledRealValue& value) {
    if (value.isZero()) return;
    if (target.isZero() || compareScaledNonnegative(value, target) < 0) target = value;
}

bool taylorPersistentBytes(const ExpressionTaylorJetResult& jet, size_t& bytes) {
    bytes = 0;
    return checkedAddSize(bytes, jet.coefficients.capacity(), sizeof(ScaledComplexValue)) && checkedAddSize(bytes, jet.coefficientRadii.capacity(), sizeof(ScaledRealValue)) && checkedAddSize(bytes, jet.intermediateEscapeMargins.capacity(), sizeof(ScaledRealValue));
}

bool makePixelOffsetBall(const std::vector<ScaledOffset>& xOffsets, const std::vector<ScaledOffset>& yOffsets, int x, int y, ScaledComplexBall& offset) {
    if (x < 0 || y < 0 || static_cast<size_t>(x) >= xOffsets.size() || static_cast<size_t>(y) >= yOffsets.size()) return false;
    const ScaledOffset& xOffset = xOffsets[static_cast<size_t>(x)];
    const ScaledOffset& yOffset = yOffsets[static_cast<size_t>(y)];
    offset.value.re = xOffset.value;
    offset.value.im = yOffset.value;
    offset.radius = compareScaledNonnegative(xOffset.error, yOffset.error) >= 0 ? xOffset.error : yOffset.error;
    return certifyScaledMpfrExponentRange(offset) == ScaledArithmeticStatus::Success;
}

#if defined(_MSC_VER)
__declspec(noinline)
#endif
bool evaluateTaylorBatch(const std::vector<ScaledOffset>& xOffsets, const ScaledOffset& yOffset, int x, int xEnd, const ExpressionTaylorJetResult& jet, size_t& count, std::array<ExpressionTaylorJetEvaluation, 4>& landing) {
    count = static_cast<size_t>(std::min(4, xEnd - x));
    std::array<ScaledComplexBall, 4> q{};
    for (size_t lane = 0; lane < count; ++lane) {
        const ScaledOffset& xOffset = xOffsets[static_cast<size_t>(x) + lane];
        ScaledComplexBall offset;
        offset.value.re = xOffset.value;
        offset.value.im = yOffset.value;
        offset.radius = compareScaledNonnegative(xOffset.error, yOffset.error) >= 0 ? xOffset.error : yOffset.error;
        if (!makeExpressionTaylorLocalQ(offset, jet.parameterOffset, jet.parameterScale, q[lane])) return false;
    }
    return ExpressionTaylorJetEvaluator::evaluateBatch(jet, q.data(), count, landing.data());
}

struct TaylorBatchCache {
    int y = -1;
    int xBegin = 0;
    size_t count = 0;
    bool valid = false;
    std::array<ExpressionTaylorJetEvaluation, 4> landing{};
};

#if defined(_MSC_VER)
__declspec(noinline)
#endif
bool getTaylorBatchLanding(TaylorBatchCache& cache, const std::vector<ScaledOffset>& xOffsets, const ScaledOffset& yOffset, int x, int y, int xEnd, const ExpressionTaylorJetResult& jet, ExpressionTaylorJetEvaluation& landing, uint64_t& evaluationNanoseconds) {
    if (cache.y != y || x < cache.xBegin || static_cast<size_t>(x - cache.xBegin) >= cache.count) {
        const Clock::time_point start = Clock::now();
        cache.y = y;
        cache.xBegin = x;
        cache.valid = evaluateTaylorBatch(xOffsets, yOffset, x, xEnd, jet, cache.count, cache.landing);
        evaluationNanoseconds += static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count());
    }
    if (!cache.valid) return false;
    landing = cache.landing[static_cast<size_t>(x - cache.xBegin)];
    return landing.valid;
}

bool makeTaylorTileParameterization(const std::vector<ScaledOffset>& xOffsets, const std::vector<ScaledOffset>& yOffsets, const PendingTaylorTile& tile, const ScaledComplexValue& frameScale, ScaledComplexBall& parameterOffset, ScaledComplexValue& parameterScale) {
    if (tile.xBegin < 0 || tile.yBegin < 0 || tile.xBegin >= tile.xEnd || tile.yBegin >= tile.yEnd || static_cast<size_t>(tile.xEnd) > xOffsets.size() || static_cast<size_t>(tile.yEnd) > yOffsets.size()) return false;
    const int centerX = tile.xBegin + (tile.xEnd - tile.xBegin - 1) / 2;
    const int centerY = tile.yBegin + (tile.yEnd - tile.yBegin - 1) / 2;
    if (!makePixelOffsetBall(xOffsets, yOffsets, centerX, centerY, parameterOffset)) return false;

    ScaledRealValue maximumRealMagnitude;
    ScaledRealValue maximumImaginaryMagnitude;
    ScaledRealValue maximumComponentError;
    for (int y = tile.yBegin; y < tile.yEnd; ++y) {
        for (int x = tile.xBegin; x < tile.xEnd; ++x) {
            ScaledComplexBall pixelOffset;
            ScaledComplexBall localOffset;
            if (!makePixelOffsetBall(xOffsets, yOffsets, x, y, pixelOffset) || certifiedScaledSubtract(pixelOffset, parameterOffset, localOffset) != ScaledArithmeticStatus::Success || certifyScaledMpfrExponentRange(localOffset) != ScaledArithmeticStatus::Success) return false;
            mergeMaximum(maximumRealMagnitude, absoluteScaled(localOffset.value.re));
            mergeMaximum(maximumImaginaryMagnitude, absoluteScaled(localOffset.value.im));
            mergeMaximum(maximumComponentError, localOffset.radius);
        }
    }
    const bool madeScale = makeExpressionTaylorFrameScale(maximumRealMagnitude, maximumImaginaryMagnitude, maximumComponentError, parameterScale);
    if (!madeScale) {
        if (!maximumRealMagnitude.isZero() || !maximumImaginaryMagnitude.isZero() || !maximumComponentError.isZero() || frameScale.re.mantissa != 0.5 || !frameScale.im.isZero()) return false;
        parameterScale = frameScale;
        const int64_t shift = static_cast<int64_t>(tile.depth) + 8;
        if (parameterScale.re.exponent < std::numeric_limits<int64_t>::min() + shift) return false;
        parameterScale.re.exponent -= shift;
        if (certifyScaledMpfrExponentRange(parameterScale) != ScaledArithmeticStatus::Success) return false;
    }
    for (int y = tile.yBegin; y < tile.yEnd; ++y) {
        for (int x = tile.xBegin; x < tile.xEnd; ++x) {
            ScaledComplexBall pixelOffset;
            ScaledComplexBall q;
            if (!makePixelOffsetBall(xOffsets, yOffsets, x, y, pixelOffset) || !makeExpressionTaylorLocalQ(pixelOffset, parameterOffset, parameterScale, q) || !expressionTaylorQInsideUnitDisk(q)) return false;
        }
    }
    return true;
}

uint64_t taylorTilePixels(const PendingTaylorTile& tile) {
    return static_cast<uint64_t>(tile.xEnd - tile.xBegin) * static_cast<uint64_t>(tile.yEnd - tile.yBegin);
}

size_t splitTaylorTile(const PendingTaylorTile& tile, const ExpressionDeepTaylorPolicy& policy, std::array<PendingTaylorTile, 4>& children) {
    const int width = tile.xEnd - tile.xBegin;
    const int height = tile.yEnd - tile.yBegin;
    const bool canSplitX = policy.minimumTileWidth <= width / 2;
    const bool canSplitY = policy.minimumTileHeight <= height / 2;
    if (tile.depth >= policy.maximumDepth || (!canSplitX && !canSplitY)) return 0;
    if (canSplitX && canSplitY) {
        const int middleX = tile.xBegin + width / 2;
        const int middleY = tile.yBegin + height / 2;
        for (PendingTaylorTile& child : children) {
            child = tile;
            child.depth = tile.depth + 1;
        }
        children[0].xEnd = middleX;
        children[0].yEnd = middleY;
        children[1].xBegin = middleX;
        children[1].yEnd = middleY;
        children[2].xEnd = middleX;
        children[2].yBegin = middleY;
        children[3].xBegin = middleX;
        children[3].yBegin = middleY;
        return 4;
    }
    children[0] = tile;
    children[1] = tile;
    children[0].depth = children[1].depth = tile.depth + 1;
    if (canSplitX) {
        const int middle = tile.xBegin + width / 2;
        children[0].xEnd = middle;
        children[1].xBegin = middle;
    } else {
        const int middle = tile.yBegin + height / 2;
        children[0].yEnd = middle;
        children[1].yBegin = middle;
    }
    return 2;
}

bool predictedTaylorBenefit(const ExpressionTaylorJetResult& jet, uint64_t pixelCount, int maximumIterations) {
    const long double fullCost = static_cast<long double>(pixelCount) * static_cast<long double>(maximumIterations) * 16.0L;
    const long double evaluation = static_cast<long double>(pixelCount) * (jet.layout == ExpressionTaylorJetLayout::RealBivariate ? 2.0L * static_cast<long double>(jet.monomialCount) + jet.order + 2.0L : 2.0L * jet.order + 2.0L);
    const long double tail = static_cast<long double>(pixelCount) * static_cast<long double>(std::max(0, maximumIterations - jet.landingIteration)) * 16.0L;
    const long double acceleratedCost = static_cast<long double>(jet.operationCount) + evaluation + tail;
    return acceleratedCost < fullCost;
}

uint64_t mixTaylorTileHash(uint64_t hash, const TaylorTileJet& tile) {
    auto mix = [&](uint64_t value) {
        hash ^= value;
        hash *= 1099511628211ULL;
    };
    mix(static_cast<uint64_t>(static_cast<uint32_t>(tile.xBegin)));
    mix(static_cast<uint64_t>(static_cast<uint32_t>(tile.yBegin)));
    mix(static_cast<uint64_t>(static_cast<uint32_t>(tile.xEnd)));
    mix(static_cast<uint64_t>(static_cast<uint32_t>(tile.yEnd)));
    mix(static_cast<uint64_t>(static_cast<uint32_t>(tile.depth)));
    mix(static_cast<uint64_t>(static_cast<uint32_t>(tile.jet.landingIteration)));
    mix(static_cast<uint64_t>(static_cast<uint32_t>(tile.jet.order)));
    return hash;
}

bool componentDiscrepancy(mpfr_srcptr exact, const ScaledRealValue& base, const ScaledRealValue& offset, mpfr_prec_t precision, ScaledRealValue& error) {
    mpfr_t baseValue, offsetValue, lower, upper;
    mpfr_t differenceLower, differenceUpper, maximum;
    mpfr_inits2(precision, baseValue, offsetValue, lower, upper, differenceLower, differenceUpper, maximum, (mpfr_ptr)0);
    const bool reconstructed = setMpfrFromScaledValue(baseValue, base) && setMpfrFromScaledValue(offsetValue, offset);
    if (reconstructed) {
        mpfr_add(lower, baseValue, offsetValue, MPFR_RNDD);
        mpfr_add(upper, baseValue, offsetValue, MPFR_RNDU);
        mpfr_sub(differenceLower, exact, upper, MPFR_RNDD);
        mpfr_sub(differenceUpper, exact, lower, MPFR_RNDU);
        mpfr_abs(differenceLower, differenceLower, MPFR_RNDU);
        mpfr_abs(differenceUpper, differenceUpper, MPFR_RNDU);
        if (mpfr_cmp(differenceLower, differenceUpper) >= 0)
            mpfr_set(maximum, differenceLower, MPFR_RNDU);
        else
            mpfr_set(maximum, differenceUpper, MPFR_RNDU);
    }
    const bool okay = reconstructed && makeScaledNonnegativeUpward(maximum, error) == ScaledArithmeticStatus::Success;
    mpfr_clears(baseValue, offsetValue, lower, upper, differenceLower, differenceUpper, maximum, (mpfr_ptr)0);
    return okay;
}

bool inflateRadius(ScaledRealValue& radius, int bits) {
    if (bits == 0 || radius.isZero()) return true;
    if (radius.exponent > std::numeric_limits<int64_t>::max() - bits) return false;
    radius.exponent += bits;
    return true;
}

bool inflateReferenceErrors(ExpressionReferenceOrbitResult& reference, int bits) {
    if (!inflateRadius(reference.cError, bits) || !inflateRadius(reference.z0Error, bits) || !inflateRadius(reference.pixelError, bits) || !inflateRadius(reference.initialZError, bits)) return false;
    for (ExpressionReferenceSample& sample : reference.samples)
        if (!inflateRadius(sample.zError, bits) || !inflateRadius(sample.nextError, bits) || !inflateRadius(sample.rootError, bits)) return false;
    for (ExpressionReferenceTapeNode& node : reference.tape)
        if (!inflateRadius(node.outputError, bits) || !inflateRadius(node.auxiliaryError, bits)) return false;
    return true;
}

struct FastPixelResult {
    bool decided = false;
    float output = ExpressionDeepEmptyPixel;
    ExpressionDeepFallbackReason reason = ExpressionDeepFallbackReason::InvalidTape;
    uint64_t iterations = 0;
    uint64_t operationCount = 0;
    uint64_t seriesOperationCount = 0;
    uint64_t foldOperationCount = 0;
    uint64_t uncertainFoldCount = 0;
    uint32_t firstUncertainIteration = 0;
};

size_t uncertaintyHistogramBin(uint32_t iteration) {
    if (iteration == 0) return 0;
    size_t bin = 1;
    while (iteration > 1 && bin + 1 < ExpressionDeepUncertaintyHistogramBins) {
        iteration >>= 1;
        ++bin;
    }
    return bin;
}

uint64_t saturatingMultiply(uint64_t left, uint64_t right) {
    if (left == 0 || right == 0) return 0;
    return left > std::numeric_limits<uint64_t>::max() / right ? std::numeric_limits<uint64_t>::max() : left * right;
}

std::vector<size_t> makePreflightSamples(int width, int height, size_t maximumSamples) {
    std::vector<size_t> samples;
    if (width < 1 || height < 1 || maximumSamples == 0) return samples;
    const size_t pixelCount = static_cast<size_t>(width) * static_cast<size_t>(height);
    maximumSamples = std::min(maximumSamples, pixelCount);
    int columns = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(maximumSamples) * static_cast<double>(width) / static_cast<double>(height))));
    columns = std::max(1, std::min(columns, width));
    int rows = static_cast<int>((maximumSamples + static_cast<size_t>(columns) - 1) / static_cast<size_t>(columns));
    rows = std::max(1, std::min(rows, height));
    samples.reserve(maximumSamples);
    const size_t gridCount = static_cast<size_t>(columns) * static_cast<size_t>(rows);
    for (size_t sample = 0; sample < maximumSamples; ++sample) {
        const size_t slot = maximumSamples == 1 ? gridCount / 2 : sample * (gridCount - 1) / (maximumSamples - 1);
        const int row = static_cast<int>(slot / static_cast<size_t>(columns));
        const int column = static_cast<int>(slot % static_cast<size_t>(columns));
        const int y = rows == 1 ? height / 2 : static_cast<int>(static_cast<int64_t>(row) * (height - 1) / (rows - 1));
        const int x = columns == 1 ? width / 2 : static_cast<int>(static_cast<int64_t>(column) * (width - 1) / (columns - 1));
        const size_t index = static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
        if (samples.empty() || samples.back() != index) samples.push_back(index);
    }
    return samples;
}

template <typename Cancel>
FastPixelResult evaluateCertifiedPerStepPixel(const ExpressionDeepRenderRequest& request, ExpressionScaledResidualCapability capability, bool certifiedPiecewiseCandidate, const ExpressionReferenceOrbitResult& reference, ExpressionScaledResidualEvaluator& evaluator, const ScaledComplexBall& initialReference, bool initialReferenceExponentUnsafe, const BailoutThreshold& bailoutThreshold, const ScaledComplexBall& offset, Cancel&& pollCancellation) {
    FastPixelResult pixel;
    ExpressionScaledResidualInput input;
    ScaledComplexBall stateDelta;
    ScaledArithmeticStatus arithmetic = initialReferenceExponentUnsafe ? ScaledArithmeticStatus::ExponentRange : certifyScaledMpfrExponentRange(offset);
    if (arithmetic != ScaledArithmeticStatus::Success) {
        pixel.reason = reasonForArithmeticStatus(arithmetic);
        return pixel;
    }
    if (request.pixelParameter == FormulaParameter::C) {
        input.c = offset.value;
        input.cError = offset.radius;
        stateDelta.radius = initialReference.radius;
    } else {
        input.z0 = offset.value;
        input.z0Error = offset.radius;
        stateDelta = offset;
    }

    ScaledComplexBall initialBase;
    initialBase.value = initialReference.value;
    ScaledComplexBall initialValue;
    arithmetic = certifyScaledMpfrExponentRange(initialBase);
    if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifyScaledMpfrExponentRange(stateDelta);
    if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifiedScaledAdd(initialBase, stateDelta, initialValue);
    if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifyScaledMpfrExponentRange(initialValue);
    if (arithmetic != ScaledArithmeticStatus::Success) {
        pixel.reason = reasonForArithmeticStatus(arithmetic);
        return pixel;
    }

    ScaledArithmeticStatus gateStatus = ScaledArithmeticStatus::Success;
    BailoutDecision decision = decideBailout(initialValue.value, initialValue.radius, bailoutThreshold, gateStatus);
    if (decision == BailoutDecision::Error) {
        pixel.reason = reasonForArithmeticStatus(gateStatus);
        return pixel;
    }
    if (decision == BailoutDecision::Uncertain) {
        pixel.reason = ExpressionDeepFallbackReason::BailoutUncertain;
        return pixel;
    }
    if (decision == BailoutDecision::Outside) {
        pixel.decided = true;
        pixel.output = 0.0f;
        return pixel;
    }

    for (int iteration = 0; iteration < request.maxIterations; ++iteration) {
        pixel.firstUncertainIteration = static_cast<uint32_t>(iteration);
        if (pollCancellation()) break;
        if (static_cast<size_t>(iteration) >= reference.samples.size()) {
            pixel.reason = ExpressionDeepFallbackReason::ReferenceExhausted;
            break;
        }
        input.z = stateDelta.value;
        input.zError = stateDelta.radius;
        input.iteration = iteration;
        const ExpressionScaledResidualResult evaluated = evaluator.evaluate(static_cast<size_t>(iteration), input);
        ++pixel.iterations;
        pixel.operationCount = saturatingAdd(pixel.operationCount, static_cast<uint64_t>(evaluated.operationCount));
        pixel.seriesOperationCount = saturatingAdd(pixel.seriesOperationCount, static_cast<uint64_t>(evaluated.seriesOperationCount));
        pixel.foldOperationCount = saturatingAdd(pixel.foldOperationCount, static_cast<uint64_t>(evaluated.foldOperationCount));
        pixel.uncertainFoldCount = saturatingAdd(pixel.uncertainFoldCount, static_cast<uint64_t>(evaluated.uncertainFoldCount));
        if (evaluated.status != ExpressionScaledResidualStatus::Success) {
            pixel.reason = reasonForResidualStatus(evaluated.status);
            break;
        }
        if (evaluated.uncertified && !request.allowUncertifiedForBenchmark) {
            pixel.reason = ExpressionDeepFallbackReason::UncertifiedSeries;
            break;
        }
        if ((capability == ExpressionScaledResidualCapability::ExactCenteredArithmetic || certifiedPiecewiseCandidate) && !evaluated.certified) {
            pixel.reason = ExpressionDeepFallbackReason::CertificationFailure;
            break;
        }

        const ExpressionReferenceSample& sample = reference.samples[static_cast<size_t>(iteration)];
        ScaledComplexBall outputBase;
        ScaledComplexBall residualBall;
        ScaledComplexBall actualOutput;
        arithmetic = makeScaledComplexValue(sample.next, sample.rootDefect, outputBase.value);
        if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifyScaledMpfrExponentRange(outputBase);
        residualBall.value = evaluated.residual;
        residualBall.radius = evaluated.radius;
        if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifiedScaledAdd(outputBase, residualBall, actualOutput);
        if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifyScaledMpfrExponentRange(actualOutput);
        if (arithmetic != ScaledArithmeticStatus::Success) {
            pixel.reason = reasonForArithmeticStatus(arithmetic);
            break;
        }

        decision = decideBailout(actualOutput.value, actualOutput.radius, bailoutThreshold, gateStatus);
        if (decision == BailoutDecision::Error) {
            pixel.reason = reasonForArithmeticStatus(gateStatus);
            break;
        }
        if (decision == BailoutDecision::Uncertain) {
            pixel.reason = ExpressionDeepFallbackReason::BailoutUncertain;
            break;
        }
        if (decision == BailoutDecision::Outside) {
            pixel.decided = true;
            pixel.output = static_cast<float>(iteration + 1);
            break;
        }
        if (iteration + 1 == request.maxIterations) {
            pixel.decided = true;
            pixel.output = ExpressionDeepInteriorPixel;
            break;
        }
        if (static_cast<size_t>(iteration + 1) >= reference.samples.size()) {
            pixel.reason = ExpressionDeepFallbackReason::ReferenceExhausted;
            break;
        }
        ScaledComplexBall nextBase;
        arithmetic = makeScaledComplexValue(reference.samples[static_cast<size_t>(iteration + 1)].z, reference.samples[static_cast<size_t>(iteration + 1)].zDefect, nextBase.value);
        if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifyScaledMpfrExponentRange(nextBase);
        if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifiedScaledSubtract(actualOutput, nextBase, stateDelta);
        if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifyScaledMpfrExponentRange(stateDelta);
        if (arithmetic != ScaledArithmeticStatus::Success) {
            pixel.reason = reasonForArithmeticStatus(arithmetic);
            break;
        }
    }
    return pixel;
}

ExpressionDeepRenderStatus mapReferenceStatus(ExpressionReferenceBuildStatus status) {
    switch (status) {
    case ExpressionReferenceBuildStatus::ProgramMismatch: return ExpressionDeepRenderStatus::ProgramMismatch;
    case ExpressionReferenceBuildStatus::PrecisionOutOfRange: return ExpressionDeepRenderStatus::PrecisionOutOfRange;
    case ExpressionReferenceBuildStatus::ResourceLimit: return ExpressionDeepRenderStatus::ResourceLimit;
    case ExpressionReferenceBuildStatus::InvalidRequest:
    case ExpressionReferenceBuildStatus::InputParseError: return ExpressionDeepRenderStatus::InvalidRequest;
    case ExpressionReferenceBuildStatus::Success:
    case ExpressionReferenceBuildStatus::UnsupportedProgram:
    case ExpressionReferenceBuildStatus::CompactionOutOfRange: return ExpressionDeepRenderStatus::ReferenceFailure;
    }
    return ExpressionDeepRenderStatus::ReferenceFailure;
}

} // namespace

const char* expressionDeepRenderStatusName(ExpressionDeepRenderStatus status) {
    switch (status) {
    case ExpressionDeepRenderStatus::Success: return "success";
    case ExpressionDeepRenderStatus::Cancelled: return "cancelled";
    case ExpressionDeepRenderStatus::InvalidRequest: return "invalid-request";
    case ExpressionDeepRenderStatus::ProgramMismatch: return "program-mismatch";
    case ExpressionDeepRenderStatus::PrecisionOutOfRange: return "precision-out-of-range";
    case ExpressionDeepRenderStatus::ResourceLimit: return "resource-limit";
    case ExpressionDeepRenderStatus::ReferenceFailure: return "reference-failure";
    case ExpressionDeepRenderStatus::UndefinedPixel: return "undefined-pixel";
    case ExpressionDeepRenderStatus::InternalError: return "internal-error";
    }
    return "internal-error";
}

const char* expressionDeepFallbackReasonName(ExpressionDeepFallbackReason reason) {
    switch (reason) {
    case ExpressionDeepFallbackReason::UncertifiedSeries: return "uncertified-series";
    case ExpressionDeepFallbackReason::BranchSensitive: return "branch-sensitive";
    case ExpressionDeepFallbackReason::UnsupportedOperation: return "unsupported-operation";
    case ExpressionDeepFallbackReason::Singular: return "singular";
    case ExpressionDeepFallbackReason::Nonfinite: return "nonfinite";
    case ExpressionDeepFallbackReason::ExponentRange: return "exponent-range";
    case ExpressionDeepFallbackReason::InvalidTape: return "invalid-tape";
    case ExpressionDeepFallbackReason::ReferenceExhausted: return "reference-exhausted";
    case ExpressionDeepFallbackReason::CertificationFailure: return "certification-failure";
    case ExpressionDeepFallbackReason::BailoutUncertain: return "bailout-uncertain";
    case ExpressionDeepFallbackReason::ReconstructionFailure: return "reconstruction-failure";
    case ExpressionDeepFallbackReason::Count: break;
    }
    return "invalid";
}

const char* expressionDeepPreflightDecisionName(ExpressionDeepPreflightDecision decision) {
    switch (decision) {
    case ExpressionDeepPreflightDecision::NotRun: return "not-run";
    case ExpressionDeepPreflightDecision::ContinueCertifiedFast: return "certified-fast";
    case ExpressionDeepPreflightDecision::DirectMpfr: return "direct-mpfr";
    }
    return "invalid";
}

bool renderExpressionDeepFrame(const ExpressionDeepRenderRequest& request, ExpressionDeepRenderResult& result) {
    result = {};
    auto fail = [&](ExpressionDeepRenderStatus status, const std::string& error) {
        result.status = status;
        result.error = error;
        result.success = false;
        result.cancelled = status == ExpressionDeepRenderStatus::Cancelled;
        return false;
    };

    try {
        if (!request.runtimeProgram || !request.runtimeProgram->valid()) return fail(ExpressionDeepRenderStatus::InvalidRequest, "runtime program is missing or invalid");
        if (request.canonicalProgram && !request.canonicalProgram->valid()) return fail(ExpressionDeepRenderStatus::InvalidRequest, "canonical program is invalid");
        if (!validCenterRepresentation(request.center) || !validScaleRepresentation(request.scale)) return fail(ExpressionDeepRenderStatus::InvalidRequest, "exact center or scale representation is invalid");
        if (request.pixelParameter != FormulaParameter::C && request.pixelParameter != FormulaParameter::InitialZ) return fail(ExpressionDeepRenderStatus::InvalidRequest, "pixel binding must be c or z0");
        if (request.width < 2 || request.height < 2 || request.width > (LONG_MAX + 1LL) / 2 || request.height > (LONG_MAX + 1LL) / 2) return fail(ExpressionDeepRenderStatus::InvalidRequest, "frame dimensions are invalid");
        if (request.maxIterations < 1 || request.maxIterations > (1 << 24)) return fail(ExpressionDeepRenderStatus::InvalidRequest, "iteration count cannot be represented exactly in float output");
        if (!(request.bailout > 0.0) || !std::isfinite(request.bailout)) return fail(ExpressionDeepRenderStatus::InvalidRequest, "bailout must be finite and positive");
        if (!request.output) return fail(ExpressionDeepRenderStatus::InvalidRequest, "output buffer is null");
        if (request.threading.tileWidth < 1 || request.threading.tileHeight < 1 || request.threading.threads < 0 || request.threading.threads > 1024) return fail(ExpressionDeepRenderStatus::InvalidRequest, "thread or tile policy is invalid");
        if (request.memory.fallbackGuardBits < 0) return fail(ExpressionDeepRenderStatus::InvalidRequest, "fallback guard precision is invalid");
        if (request.preflight.maximumSamples < 1 || request.preflight.maximumSamples > 256 || request.preflight.minimumSamples < 1 || request.preflight.minimumSamples > request.preflight.maximumSamples || request.preflight.earlyRejectMinimumFirstUncertainIteration > static_cast<uint32_t>(1 << 24)) return fail(ExpressionDeepRenderStatus::InvalidRequest, "preflight policy is invalid");
        if (request.taylor.minimumLanding < 1 || request.taylor.minimumOrder < 8 || request.taylor.maximumOrder > 20 || request.taylor.maximumBivariateOrder < 8 || request.taylor.maximumBivariateOrder > ExpressionTaylorMaximumBivariateOrder || request.taylor.minimumOrder > request.taylor.order || request.taylor.order > request.taylor.maximumOrder || request.taylor.maximumCompositionOrder < 1 || request.taylor.maximumCompositionOrder > 24 || request.taylor.maximumCandidateIteration < 0 || request.taylor.maximumDepth < 0 || request.taylor.maximumDepth > 30 || request.taylor.minimumTileWidth < 1 || request.taylor.minimumTileHeight < 1 || request.taylor.maximumJetCount < 1 || request.taylor.maximumJetCount > size_t{65536} || request.taylor.maximumRejectedBeforeFirstAcceptance < 1 || request.taylor.maximumRejectedBeforeFirstAcceptance > request.taylor.maximumJetCount || !(request.taylor.accuracyBudget > 0.0) || !std::isfinite(request.taylor.accuracyBudget))
            return fail(ExpressionDeepRenderStatus::InvalidRequest, "Taylor policy is invalid");
        if (request.verificationErrorInflationBits < 0 || request.verificationErrorInflationBits > 4096) return fail(ExpressionDeepRenderStatus::InvalidRequest, "verification error inflation is invalid");
        if (request.verificationFault < ExpressionDeepVerificationFault::None || request.verificationFault > ExpressionDeepVerificationFault::FallbackIterationAllocation) return fail(ExpressionDeepRenderStatus::InvalidRequest, "verification fault selection is invalid");
        if (request.precision.requestedBits < 0 || request.precision.viewBits < 0 || request.precision.guardBits < 0 || request.precision.minimumBits < 0 || request.precision.maximumBits < MPFR_PREC_MIN || request.precision.maximumBits > MPFR_PREC_MAX) return fail(ExpressionDeepRenderStatus::InvalidRequest, "reference precision policy is invalid");
        if (!finiteComplex(request.fixed.z) || !finiteComplex(request.fixed.c) || !finiteComplex(request.fixed.z0)) return fail(ExpressionDeepRenderStatus::InvalidRequest, "fixed context is nonfinite");
        for (Complex parameter : request.fixed.parameters)
            if (!finiteComplex(parameter)) return fail(ExpressionDeepRenderStatus::InvalidRequest, "fixed parameter is nonfinite");

        ExpressionProgram independentlySpecialized;
        if (request.canonicalProgram) {
            ExpressionError expressionError;
            if (!request.canonicalProgram->specialize(request.fixed, request.pixelParameter, independentlySpecialized, &expressionError)) return fail(ExpressionDeepRenderStatus::InvalidRequest, expressionError.message.empty() ? "canonical specialization failed" : expressionError.message);
            if (request.runtimeProgram->source() != request.canonicalProgram->source() || !request.runtimeProgram->semanticallyEquivalent(independentlySpecialized)) return fail(ExpressionDeepRenderStatus::ProgramMismatch, "runtime expression does not match canonical specialization");
        }

        const uint64_t width = static_cast<uint64_t>(request.width);
        const uint64_t height = static_cast<uint64_t>(request.height);
        if (height != 0 && width > std::numeric_limits<uint64_t>::max() / height) return fail(ExpressionDeepRenderStatus::ResourceLimit, "pixel count overflow");
        const uint64_t pixelCount64 = width * height;
        if (pixelCount64 > static_cast<uint64_t>(std::numeric_limits<size_t>::max()) || request.outputCount < static_cast<size_t>(pixelCount64)) return fail(ExpressionDeepRenderStatus::InvalidRequest, "output buffer is smaller than the frame");
        if (pixelCount64 > std::numeric_limits<uint64_t>::max() / static_cast<uint64_t>(request.maxIterations)) return fail(ExpressionDeepRenderStatus::ResourceLimit, "iteration resource multiplication overflow");
        const size_t pixelCount = static_cast<size_t>(pixelCount64);
        const uint64_t tilesX = (width + static_cast<uint64_t>(request.threading.tileWidth) - 1) / static_cast<uint64_t>(request.threading.tileWidth);
        const uint64_t tilesY = (height + static_cast<uint64_t>(request.threading.tileHeight) - 1) / static_cast<uint64_t>(request.threading.tileHeight);
        if (tilesY != 0 && tilesX > std::numeric_limits<uint64_t>::max() / tilesY) return fail(ExpressionDeepRenderStatus::ResourceLimit, "tile count overflow");
        const uint64_t tileCount = tilesX * tilesY;
        if (tileCount > static_cast<uint64_t>(LLONG_MAX) || pixelCount64 > static_cast<uint64_t>(LLONG_MAX)) return fail(ExpressionDeepRenderStatus::ResourceLimit, "OpenMP loop range is too large");

        int threadCount = request.threading.threads;
#ifdef _OPENMP
        if (threadCount == 0) threadCount = omp_get_max_threads();
#else
        threadCount = 1;
#endif
        threadCount = std::max(1, threadCount);

        result.capability = request.runtimeProgram->scaledResidualCapability();
        const bool certifiedTaylorCandidate = certifiedTaylorCapability(result.capability);
        const bool certifiedPiecewiseCandidate = result.capability == ExpressionScaledResidualCapability::CertifiedPiecewiseCandidate;
        const bool piecewisePerStepEligible = certifiedPiecewiseCandidate && !request.runtimeProgram->scaledResidualRequiresTaylor();
        const bool earlyPiecewisePreflight = piecewisePerStepEligible;
        const bool certifiedTaylorEligible = !certifiedTaylorCandidate || piecewisePerStepEligible || (request.taylor.enableTaylor && request.maxIterations > request.taylor.minimumLanding);
        bool runFast = !request.forceMpfrFallbackForVerification && certifiedTaylorEligible && (certifiedReferenceCapability(result.capability) || (result.capability == ExpressionScaledResidualCapability::UncertifiedSeries && request.allowUncertifiedForBenchmark));
        bool certificationUnavailable = request.forceMpfrFallbackForVerification || (certifiedTaylorCandidate && !certifiedTaylorEligible);

        size_t rendererBaseBytes = 0;
        if ((runFast && (!checkedAddSize(rendererBaseBytes, request.width, sizeof(ScaledOffset)) || !checkedAddSize(rendererBaseBytes, request.height, sizeof(ScaledOffset)) || (request.preflight.enable && piecewisePerStepEligible && !checkedAddSize(rendererBaseBytes, request.preflight.maximumSamples, sizeof(size_t))))) || !checkedAddSize(rendererBaseBytes, pixelCount, sizeof(uint8_t)) || !checkedAddSize(rendererBaseBytes, pixelCount, sizeof(size_t))) return fail(ExpressionDeepRenderStatus::ResourceLimit, "renderer memory calculation overflow");
        size_t fastThreadBytes = 2048;
        if (runFast) {
            const size_t instructionCount = request.runtimeProgram->instructionCount();
            if (!checkedAddSize(fastThreadBytes, instructionCount, sizeof(ScaledComplexValue) + sizeof(ScaledRealValue) + sizeof(uint16_t) + 16) || !checkedAddSize(fastThreadBytes, static_cast<size_t>(request.maxIterations), 2 * sizeof(uint8_t))) return fail(ExpressionDeepRenderStatus::ResourceLimit, "thread workspace calculation overflow");
            if (instructionCount != 0 && static_cast<size_t>(request.maxIterations) > std::numeric_limits<size_t>::max() / instructionCount) return fail(ExpressionDeepRenderStatus::ResourceLimit, "reference tape workspace multiplication overflow");
            const size_t maximumTapeNodes = static_cast<size_t>(request.maxIterations) * instructionCount;
            if (!checkedAddSize(fastThreadBytes, maximumTapeNodes, 2 * sizeof(ScaledComplexValue) + sizeof(ScaledRealValue) + sizeof(uint8_t))) return fail(ExpressionDeepRenderStatus::ResourceLimit, "reference tape workspace calculation overflow");
        }
        size_t fastThreadBytesTotal = 0;
        if (runFast && !checkedAddSize(fastThreadBytesTotal, static_cast<size_t>(threadCount), fastThreadBytes)) return fail(ExpressionDeepRenderStatus::ResourceLimit, "thread workspace multiplication overflow");
        size_t rendererBytes = rendererBaseBytes;
        if (!checkedAddSize(rendererBytes, 1, fastThreadBytesTotal)) return fail(ExpressionDeepRenderStatus::ResourceLimit, "renderer memory calculation overflow");
        result.rendererBytes = rendererBytes;
        if (request.memory.memoryLimitBytes != 0 && rendererBytes > request.memory.memoryLimitBytes) {
            if (!certifiedTaylorCandidate) return fail(ExpressionDeepRenderStatus::ResourceLimit, "renderer exceeds memory limit");
            runFast = false;
            certificationUnavailable = true;
            rendererBaseBytes = 0;
            if (!checkedAddSize(rendererBaseBytes, pixelCount, sizeof(uint8_t)) || !checkedAddSize(rendererBaseBytes, pixelCount, sizeof(size_t))) return fail(ExpressionDeepRenderStatus::ResourceLimit, "fallback renderer memory calculation overflow");
            fastThreadBytesTotal = 0;
            rendererBytes = rendererBaseBytes;
            result.rendererBytes = rendererBytes;
            if (rendererBytes > request.memory.memoryLimitBytes) return fail(ExpressionDeepRenderStatus::ResourceLimit, "fallback renderer exceeds memory limit");
        }

        ExpressionReferencePrecisionPolicy precision = request.precision;
        std::string precisionError;
        if (!selectAutomaticViewBits(request, precision, precisionError)) return fail(ExpressionDeepRenderStatus::InvalidRequest, precisionError);
        if (!selectDeepPrecision(precision, result.selectedPrecision, precisionError)) return fail(ExpressionDeepRenderStatus::PrecisionOutOfRange, precisionError);
        if (result.selectedPrecision <= MPFR_PREC_MAX - request.memory.fallbackGuardBits) result.fallbackPrecision = result.selectedPrecision + request.memory.fallbackGuardBits;

        if (runFast && certifiedReferenceCapability(result.capability)) {
            const mpfr_prec_t certificationGuard = std::max<mpfr_prec_t>(128, precision.guardBits);
            if (certificationGuard > precision.maximumBits || certificationGuard > ExpressionReferencePrecisionPolicy::ApplicationMaximumBits || result.selectedPrecision > MPFR_PREC_MAX - certificationGuard || result.selectedPrecision > precision.maximumBits - certificationGuard || result.selectedPrecision > ExpressionReferencePrecisionPolicy::ApplicationMaximumBits - certificationGuard) {
                runFast = false;
                certificationUnavailable = true;
            } else {
                result.certificationPrecision = result.selectedPrecision + certificationGuard;
            }
        }

        std::fill(request.output, request.output + pixelCount, ExpressionDeepEmptyPixel);

        std::atomic_bool cancelled{false};
        auto pollCancellation = [&]() {
            if (cancelled.load(std::memory_order_acquire)) return true;
            if (request.shouldCancel && request.shouldCancel()) {
                cancelled.store(true);
                return true;
            }
            return false;
        };
        std::mutex progressMutex;
        auto notifyProgress = [&](ExpressionDeepRenderPhase phase, uint64_t completed, uint64_t total) {
            if (!request.progress) return;
            std::lock_guard<std::mutex> lock(progressMutex);
            request.progress(phase, completed, total);
        };

        ExpressionReferenceOrbitResult reference;
        if (runFast) {
            notifyProgress(ExpressionDeepRenderPhase::Reference, 0, pixelCount64);
            ExpressionReferenceBuildRequest referenceRequest;
            referenceRequest.canonicalProgram = request.canonicalProgram;
            referenceRequest.runtimeProgram = request.runtimeProgram;
            referenceRequest.pixelParameter = request.pixelParameter;
            referenceRequest.center = request.center;
            referenceRequest.fixed = request.fixed;
            referenceRequest.bailout = request.bailout;
            referenceRequest.maxIterations = request.maxIterations;
            referenceRequest.precision = precision;
            referenceRequest.certificationPrecision = certifiedReferenceCapability(result.capability) ? result.certificationPrecision : 0;
            referenceRequest.memoryLimitBytes = request.memory.memoryLimitBytes;
            referenceRequest.shouldCancel = pollCancellation;

            const Clock::time_point referenceStart = Clock::now();
            bool referenceBuilt = false;
            {
                struct OracleWorkspaceRelease {
                    ~OracleWorkspaceRelease() { ExpressionOracle::releaseThreadWorkspace(); }
                } releaseOracleWorkspace;
                referenceBuilt = buildExpressionReferenceOrbit(referenceRequest, reference);
            }
            result.referenceSeconds = secondsSince(referenceStart);
            result.referenceBytes = reference.memoryBytes;
            if (cancelled.load(std::memory_order_acquire) || reference.cancelled) return fail(ExpressionDeepRenderStatus::Cancelled, "render cancelled while building the reference");
            if (!referenceBuilt) {
                const bool mayFallback = certifiedReferenceCapability(result.capability) && reference.status != ExpressionReferenceBuildStatus::ProgramMismatch && reference.status != ExpressionReferenceBuildStatus::InvalidRequest && reference.status != ExpressionReferenceBuildStatus::InputParseError;
                if (!mayFallback) return fail(mapReferenceStatus(reference.status), reference.error.empty() ? "reference build failed" : reference.error);
                runFast = false;
                certificationUnavailable = true;
                reference = {};
                result.referenceBytes = 0;
            } else if (!reference.valid) {
                runFast = false;
                certificationUnavailable = true;
                reference = {};
                result.referenceBytes = 0;
            } else if (certifiedReferenceCapability(result.capability) && (!reference.certifiedAgainstHigherPrecision || reference.certificationPrecision != result.certificationPrecision)) {
                runFast = false;
                certificationUnavailable = true;
                reference = {};
                result.referenceBytes = 0;
            }
            if (runFast && !inflateReferenceErrors(reference, request.verificationErrorInflationBits)) return fail(ExpressionDeepRenderStatus::InternalError, "certification radius inflation overflow");
            if (runFast && request.memory.memoryLimitBytes != 0 && (reference.memoryBytes > request.memory.memoryLimitBytes || rendererBytes > request.memory.memoryLimitBytes - reference.memoryBytes)) {
                if (!certifiedTaylorCandidate) return fail(ExpressionDeepRenderStatus::ResourceLimit, "reference plus renderer exceeds memory limit");
                runFast = false;
                certificationUnavailable = true;
                reference = {};
                result.referenceBytes = 0;
                rendererBaseBytes = 0;
                if (!checkedAddSize(rendererBaseBytes, pixelCount, sizeof(uint8_t)) || !checkedAddSize(rendererBaseBytes, pixelCount, sizeof(size_t))) return fail(ExpressionDeepRenderStatus::ResourceLimit, "fallback renderer memory calculation overflow");
                fastThreadBytesTotal = 0;
                rendererBytes = rendererBaseBytes;
                result.rendererBytes = rendererBytes;
            }
        }

        std::vector<uint8_t> fallbackReason(pixelCount, NoFallbackReason);
        std::vector<ScaledOffset> xOffsets;
        std::vector<ScaledOffset> yOffsets;
        ScaledComplexBall initialReference;
        bool initialReferenceExponentUnsafe = false;
        BailoutThreshold bailoutThreshold;
        std::unique_ptr<ExpressionScaledResidualEvaluator> validationEvaluator;
        ExpressionTaylorJetResult taylorJet;
        std::vector<TaylorTileJet> taylorTileJets;
        std::vector<int32_t> taylorPixelJet;
        bool useTaylor = false;
        bool useGlobalTaylor = false;
        size_t retainedTaylorBytes = 0;
        if (runFast) {
            ExactGeometry geometry(reference.precision);
            ExactGeometry certificationGeometry(reference.certificationPrecision != 0 ? reference.certificationPrecision : reference.precision);
            if (!geometry.initialize(request) || !certificationGeometry.initialize(request)) return fail(ExpressionDeepRenderStatus::InternalError, "failed to reconstruct reference geometry");
            if (!makeBailoutThreshold(request.bailout, bailoutThreshold)) return fail(ExpressionDeepRenderStatus::InternalError, "scaled geometry or bailout is not representable");

            xOffsets.resize(static_cast<size_t>(request.width));
            yOffsets.resize(static_cast<size_t>(request.height));
            ScaledComplexValue pixelBase;
            const ScaledComplexShadow& pixelPrimary = request.pixelParameter == FormulaParameter::C ? reference.c : reference.z0;
            const ScaledComplexShadow& pixelDefect = request.pixelParameter == FormulaParameter::C ? reference.cDefect : reference.z0Defect;
            if (makeScaledComplexValue(pixelPrimary, pixelDefect, pixelBase) != ScaledArithmeticStatus::Success) return fail(ExpressionDeepRenderStatus::InternalError, "compact pixel reference is not representable");
            ScaledComplexBall initialPrimary;
            ScaledComplexBall initialDefect;
            if (makeScaledComplexValue(reference.initialZ, initialPrimary.value) != ScaledArithmeticStatus::Success || makeScaledComplexValue(reference.initialZDefect, initialDefect.value) != ScaledArithmeticStatus::Success || certifiedScaledAdd(initialPrimary, initialDefect, initialReference) != ScaledArithmeticStatus::Success || scaledAddUp(initialReference.radius, reference.initialZError, initialReference.radius) != ScaledArithmeticStatus::Success) return fail(ExpressionDeepRenderStatus::InternalError, "compact initial reference is not representable");
            initialReferenceExponentUnsafe = certifyScaledMpfrExponentRange(initialPrimary) != ScaledArithmeticStatus::Success || certifyScaledMpfrExponentRange(initialDefect) != ScaledArithmeticStatus::Success || certifyScaledMpfrExponentRange(initialReference) != ScaledArithmeticStatus::Success;
            MpfrComplex lowOffset(reference.precision);
            MpfrComplex exactCoordinate(certificationGeometry.center.precision());
            for (int x = 0; x < request.width; ++x) {
                const long centered = static_cast<long>(2LL * x - (request.width - 1LL));
                ScaledOffset& offset = xOffsets[static_cast<size_t>(x)];
                mpfr_mul_si(lowOffset.re, geometry.dxHalf, centered, RND);
                if (makeScaledRealValue(lowOffset.re, offset.value) != ScaledArithmeticStatus::Success || !certificationGeometry.coordinate(x, 0, request, exactCoordinate) || !componentDiscrepancy(exactCoordinate.re, pixelBase.re, offset.value, certificationGeometry.center.precision(), offset.error) || !inflateRadius(offset.error, request.verificationErrorInflationBits)) return fail(ExpressionDeepRenderStatus::InternalError, "certified x coordinate construction failed");
            }
            for (int y = 0; y < request.height; ++y) {
                const long centered = static_cast<long>(2LL * y - (request.height - 1LL));
                ScaledOffset& offset = yOffsets[static_cast<size_t>(y)];
                mpfr_mul_si(lowOffset.im, geometry.dyHalf, centered, RND);
                if (makeScaledRealValue(lowOffset.im, offset.value) != ScaledArithmeticStatus::Success || !certificationGeometry.coordinate(0, y, request, exactCoordinate) || !componentDiscrepancy(exactCoordinate.im, pixelBase.im, offset.value, certificationGeometry.center.precision(), offset.error) || !inflateRadius(offset.error, request.verificationErrorInflationBits)) return fail(ExpressionDeepRenderStatus::InternalError, "certified y coordinate construction failed");
            }
            validationEvaluator = std::make_unique<ExpressionScaledResidualEvaluator>();
            if (reference.defectPending || !validationEvaluator->prepare(*request.runtimeProgram, reference)) {
                runFast = false;
                std::fill(fallbackReason.begin(), fallbackReason.end(), static_cast<uint8_t>(ExpressionDeepFallbackReason::InvalidTape));
            } else {
                size_t validationThreadBytesTotal = 0;
                if (!checkedAddSize(validationThreadBytesTotal, static_cast<size_t>(threadCount), validationEvaluator->workspaceBytes())) return fail(ExpressionDeepRenderStatus::ResourceLimit, "validation workspace multiplication overflow");
                size_t validationPeak = rendererBaseBytes;
                if (!checkedAddSize(validationPeak, 1, validationThreadBytesTotal)) return fail(ExpressionDeepRenderStatus::ResourceLimit, "validation workspace calculation overflow");
                rendererBytes = std::max(rendererBytes, validationPeak);
                result.rendererBytes = rendererBytes;
                if (request.memory.memoryLimitBytes != 0 && (reference.memoryBytes > request.memory.memoryLimitBytes || rendererBytes > request.memory.memoryLimitBytes - reference.memoryBytes)) {
                    if (!certifiedTaylorCandidate) return fail(ExpressionDeepRenderStatus::ResourceLimit, "validation workspace exceeds memory limit");
                    runFast = false;
                    certificationUnavailable = true;
                    reference = {};
                    result.referenceBytes = 0;
                    validationEvaluator.reset();
                    std::vector<ScaledOffset>().swap(xOffsets);
                    std::vector<ScaledOffset>().swap(yOffsets);
                    rendererBaseBytes = 0;
                    if (!checkedAddSize(rendererBaseBytes, pixelCount, sizeof(uint8_t)) || !checkedAddSize(rendererBaseBytes, pixelCount, sizeof(size_t))) return fail(ExpressionDeepRenderStatus::ResourceLimit, "fallback renderer memory calculation overflow");
                    fastThreadBytesTotal = 0;
                    rendererBytes = rendererBaseBytes;
                    result.rendererBytes = rendererBytes;
                }
            }
            auto rejectFastFromPreflight = [&]() {
                result.preflightRejectedFast = true;
                result.preflightDecision = ExpressionDeepPreflightDecision::DirectMpfr;
                result.preflightAvoidedFastPixelCount = pixelCount64;
                const uint64_t sampleCount = std::max<uint64_t>(1, result.preflightSampleCount);
                const uint64_t averageIterations = (result.preflightIterationCount + sampleCount - 1) / sampleCount;
                const uint64_t averageOperations = (result.preflightOperationCount + sampleCount - 1) / sampleCount;
                result.preflightAvoidedIterationEstimate = saturatingMultiply(averageIterations, pixelCount64);
                result.preflightAvoidedOperationEstimate = saturatingMultiply(averageOperations, pixelCount64);
                runFast = false;
                certificationUnavailable = true;
                reference = {};
                result.referenceBytes = 0;
                validationEvaluator.reset();
                std::vector<ScaledOffset>().swap(xOffsets);
                std::vector<ScaledOffset>().swap(yOffsets);
                rendererBaseBytes = 0;
                if (!checkedAddSize(rendererBaseBytes, pixelCount, sizeof(uint8_t)) || !checkedAddSize(rendererBaseBytes, pixelCount, sizeof(size_t))) return fail(ExpressionDeepRenderStatus::ResourceLimit, "preflight fallback renderer memory calculation overflow");
                fastThreadBytesTotal = 0;
                rendererBytes = rendererBaseBytes;
                result.rendererBytes = rendererBytes;
                std::fill(fallbackReason.begin(), fallbackReason.end(), static_cast<uint8_t>(ExpressionDeepFallbackReason::CertificationFailure));
                return true;
            };
            auto runCertifiedPreflight = [&]() {
                result.preflightAttempted = true;
                notifyProgress(ExpressionDeepRenderPhase::Preflight, 0, request.preflight.maximumSamples);
                const Clock::time_point preflightStart = Clock::now();
                const std::vector<size_t> samples = makePreflightSamples(request.width, request.height, request.preflight.maximumSamples);
                for (size_t sampleNumber = 0; sampleNumber < samples.size(); ++sampleNumber) {
                    if (pollCancellation()) break;
                    const size_t index = samples[sampleNumber];
                    const int y = static_cast<int>(index / static_cast<size_t>(request.width));
                    const int x = static_cast<int>(index % static_cast<size_t>(request.width));
                    ScaledComplexBall offset;
                    const ScaledOffset& xOffset = xOffsets[static_cast<size_t>(x)];
                    const ScaledOffset& yOffset = yOffsets[static_cast<size_t>(y)];
                    offset.value.re = xOffset.value;
                    offset.value.im = yOffset.value;
                    offset.radius = compareScaledNonnegative(xOffset.error, yOffset.error) >= 0 ? xOffset.error : yOffset.error;
                    const FastPixelResult pixel = evaluateCertifiedPerStepPixel(request, result.capability, certifiedPiecewiseCandidate, reference, *validationEvaluator, initialReference, initialReferenceExponentUnsafe, bailoutThreshold, offset, pollCancellation);
                    ++result.preflightSampleCount;
                    result.preflightIterationCount = saturatingAdd(result.preflightIterationCount, pixel.iterations);
                    result.preflightOperationCount = saturatingAdd(result.preflightOperationCount, pixel.operationCount);
                    result.preflightFoldOperationCount = saturatingAdd(result.preflightFoldOperationCount, pixel.foldOperationCount);
                    result.preflightUncertainFoldCount = saturatingAdd(result.preflightUncertainFoldCount, pixel.uncertainFoldCount);
                    if (pixel.decided) {
                        ++result.preflightFastCount;
                    } else {
                        ++result.preflightFallbackCount;
                        result.preflightMinimumFirstUncertainIteration = std::min(result.preflightMinimumFirstUncertainIteration, pixel.firstUncertainIteration);
                        result.preflightMaximumFirstUncertainIteration = std::max(result.preflightMaximumFirstUncertainIteration, pixel.firstUncertainIteration);
                        const size_t reasonIndex = static_cast<size_t>(pixel.reason);
                        if (reasonIndex < result.preflightFallbackReasonCounts.size()) ++result.preflightFallbackReasonCounts[reasonIndex];
                        ++result.preflightFirstUncertainHistogram[uncertaintyHistogramBin(pixel.firstUncertainIteration)];
                    }
                    notifyProgress(ExpressionDeepRenderPhase::Preflight, result.preflightSampleCount, samples.size());
                }
                result.preflightSeconds = secondsSince(preflightStart);
                if (cancelled.load(std::memory_order_acquire)) return fail(ExpressionDeepRenderStatus::Cancelled, "render cancelled during certified preflight");
                const uint64_t minimumSamples = std::min<uint64_t>(static_cast<uint64_t>(request.preflight.minimumSamples), pixelCount64);
                const bool rejectFast = result.preflightSampleCount >= minimumSamples && result.preflightFallbackCount == result.preflightSampleCount && result.preflightFastCount == 0 && (!request.taylor.enableTaylor || result.preflightMinimumFirstUncertainIteration >= request.preflight.earlyRejectMinimumFirstUncertainIteration);
                if (rejectFast) {
                    if (!rejectFastFromPreflight()) return false;
                } else {
                    result.preflightDecision = ExpressionDeepPreflightDecision::ContinueCertifiedFast;
                }
                return true;
            };
            if (runFast && request.preflight.enable && (!request.taylor.enableTaylor || earlyPiecewisePreflight) && piecewisePerStepEligible && validationEvaluator && validationEvaluator->ready()) {
                if (!runCertifiedPreflight()) return false;
            }
            // Validation is complete; destroy its retained vector capacities
            // before Taylor construction and the parallel render peak.
            validationEvaluator.reset();

            if (runFast && request.taylor.enableTaylor && certifiedReferenceCapability(result.capability) && request.maxIterations > request.taylor.minimumLanding) {
                const Clock::time_point taylorPlanningStart = Clock::now();
                const double taylorBuildSecondsBefore = result.taylorBuildSeconds;
                result.taylorAttempted = true;
                long double weightedLandingTotal = 0.0L;
                size_t taylorBuildPeakBytes = 0;
                auto recordJet = [&](const ExpressionTaylorJetResult& jet, bool accepted, uint64_t pixels, int depth) {
                    ++result.taylorAttemptedJetCount;
                    if (accepted) {
                        ++result.taylorAcceptedJetCount;
                        result.taylorAcceptedPixelCoverage = saturatingAdd(result.taylorAcceptedPixelCoverage, pixels);
                        weightedLandingTotal += static_cast<long double>(pixels) * static_cast<long double>(jet.landingIteration);
                        result.taylorStatus = ExpressionTaylorJetStatus::Success;
                    } else {
                        ++result.taylorRejectedJetCount;
                        if (result.taylorAcceptedJetCount == 0) {
                            result.taylorStatus = jet.status;
                            result.taylorFailureReason = jet.failureReason;
                        }
                    }
                    result.taylorMaximumTileDepth = std::max(result.taylorMaximumTileDepth, depth);
                    result.taylorOrder = std::max(result.taylorOrder, jet.order);
                    result.taylorLayout = jet.layout;
                    result.taylorMonomialCount = std::max(result.taylorMonomialCount, jet.monomialCount);
                    result.taylorBivariateConvolutionOperationCount = saturatingAdd(result.taylorBivariateConvolutionOperationCount, jet.bivariateConvolutionOperationCount);
                    result.taylorCoveredIterations = std::max(result.taylorCoveredIterations, jet.landingIteration);
                    result.taylorMaximumFunctionSeriesOrder = std::max(result.taylorMaximumFunctionSeriesOrder, jet.maximumFunctionSeriesOrder);
                    result.taylorFunctionSeriesCount = saturatingAdd(result.taylorFunctionSeriesCount, jet.functionSeriesCount);
                    result.taylorFunctionSeriesOperationCount = saturatingAdd(result.taylorFunctionSeriesOperationCount, jet.functionSeriesOperationCount);
                    mergeMaximum(result.taylorMaximumFunctionSeriesTail, jet.maximumFunctionSeriesTail);
                    result.taylorMaximumReciprocalOrder = std::max(result.taylorMaximumReciprocalOrder, jet.maximumReciprocalOrder);
                    result.taylorReciprocalCount = saturatingAdd(result.taylorReciprocalCount, jet.reciprocalCount);
                    result.taylorReciprocalOperationCount = saturatingAdd(result.taylorReciprocalOperationCount, jet.reciprocalOperationCount);
                    mergeMinimumNonzero(result.taylorMinimumDenominatorClearance, jet.minimumDenominatorClearance);
                    mergeMaximum(result.taylorMaximumReciprocalTail, jet.maximumReciprocalTail);
                    result.taylorMaximumBranchSeriesOrder = std::max(result.taylorMaximumBranchSeriesOrder, jet.maximumBranchSeriesOrder);
                    result.taylorBranchCompositionCount = saturatingAdd(result.taylorBranchCompositionCount, jet.branchCompositionCount);
                    result.taylorBranchCompositionOperationCount = saturatingAdd(result.taylorBranchCompositionOperationCount, jet.branchCompositionOperationCount);
                    mergeMaximum(result.taylorMaximumBranchSeriesTail, jet.maximumBranchSeriesTail);
                    mergeMinimumNonzero(result.taylorMinimumBranchCutClearance, jet.minimumBranchCutClearance);
                    mergeMinimumNonzero(result.taylorMinimumBranchZeroClearance, jet.minimumBranchZeroClearance);
                    result.taylorArgCompositionCount = saturatingAdd(result.taylorArgCompositionCount, jet.argCompositionCount);
                    if (result.taylorArgRejectionReason.empty()) result.taylorArgRejectionReason = jet.argRejectionReason;
                    result.taylorPolarCompositionCount = saturatingAdd(result.taylorPolarCompositionCount, jet.polarCompositionCount);
                    mergeMinimumNonzero(result.taylorMinimumPolarRadiusClearance, jet.minimumPolarRadiusClearance);
                    if (result.taylorPolarRejectionReason.empty()) result.taylorPolarRejectionReason = jet.polarRejectionReason;
                    result.taylorAbsBranchCount = saturatingAdd(result.taylorAbsBranchCount, jet.absBranchCount);
                    result.taylorAbsPositiveCellCount = saturatingAdd(result.taylorAbsPositiveCellCount, jet.absPositiveCellCount);
                    result.taylorAbsNegativeCellCount = saturatingAdd(result.taylorAbsNegativeCellCount, jet.absNegativeCellCount);
                    mergeMinimumNonzero(result.taylorMinimumFoldClearance, jet.minimumFoldClearance);
                    if (jet.poleRejected || jet.status == ExpressionTaylorJetStatus::PoleRejected) {
                        result.taylorPoleRejected = true;
                        ++result.taylorPoleRejectedJetCount;
                    }
                    if (jet.branchRejected && !jet.foldRejected) {
                        result.taylorBranchRejected = true;
                        ++result.taylorCutRejectedJetCount;
                    }
                    if (jet.polarRejected) result.taylorPolarRejected = true;
                    if (jet.foldRejected) {
                        result.taylorFoldRejected = true;
                        ++result.taylorFoldRejectedJetCount;
                        if (result.taylorFoldRejectionIteration < 0) result.taylorFoldRejectionIteration = jet.foldRejectionIteration;
                        if (result.taylorFoldRejectionReason.empty()) result.taylorFoldRejectionReason = jet.foldRejectionReason;
                    }
                };
                auto makeJetRequest = [&](const ScaledComplexBall& parameterOffset, const ScaledComplexValue& parameterScale, size_t retainedBytes, ExpressionTaylorJetRequest& jetRequest) {
                    jetRequest.program = request.runtimeProgram;
                    jetRequest.reference = &reference;
                    jetRequest.pixelParameter = request.pixelParameter;
                    jetRequest.parameterOffset = parameterOffset;
                    jetRequest.parameterScale = parameterScale;
                    jetRequest.minimumOrder = request.taylor.minimumOrder;
                    jetRequest.preferredOrder = request.taylor.order;
                    jetRequest.maximumOrder = request.taylor.maximumOrder;
                    jetRequest.maximumBivariateOrder = request.taylor.maximumBivariateOrder;
                    jetRequest.maximumCompositionOrder = request.taylor.maximumCompositionOrder;
                    jetRequest.minimumLanding = request.taylor.minimumLanding;
                    jetRequest.maximumCandidateIteration = request.taylor.maximumCandidateIteration;
                    jetRequest.bailout = request.bailout;
                    jetRequest.accuracyBudget = request.taylor.accuracyBudget;
                    jetRequest.shouldCancel = pollCancellation;

                    size_t available = std::numeric_limits<size_t>::max();
                    bool limited = false;
                    if (request.memory.memoryLimitBytes != 0) {
                        limited = true;
                        size_t fixed = reference.memoryBytes;
                        if (!checkedAddSize(fixed, 1, rendererBaseBytes)) return false;
                        if (fixed > request.memory.memoryLimitBytes || retainedBytes > request.memory.memoryLimitBytes - fixed) return false;
                        available = request.memory.memoryLimitBytes - fixed - retainedBytes;
                    }
                    if (request.taylor.maximumJetMemoryBytes != 0) {
                        limited = true;
                        if (retainedBytes > request.taylor.maximumJetMemoryBytes) return false;
                        available = std::min(available, request.taylor.maximumJetMemoryBytes - retainedBytes);
                    }
                    if (limited && available == 0) return false;
                    jetRequest.memoryLimitBytes = limited ? available : 0;
                    return true;
                };
                auto attemptJet = [&](const ScaledComplexBall& parameterOffset, const ScaledComplexValue& parameterScale, uint64_t pixels, int depth, size_t retainedBytes, ExpressionTaylorJetResult& jet, size_t& persistentBytes) {
                    persistentBytes = 0;
                    ExpressionTaylorJetRequest jetRequest;
                    bool built = false;
                    if (!makeJetRequest(parameterOffset, parameterScale, retainedBytes, jetRequest)) {
                        jet.status = ExpressionTaylorJetStatus::ResourceLimit;
                        jet.failureReason = "Taylor planner memory budget is exhausted";
                    } else {
                        built = ExpressionTaylorJetBuilder::build(jetRequest, jet);
                    }
                    result.taylorBuildSeconds += jet.buildSeconds;
                    size_t combinedPeak = retainedBytes;
                    if (!checkedAddSize(combinedPeak, 1, jet.memoryBytes)) combinedPeak = std::numeric_limits<size_t>::max();
                    taylorBuildPeakBytes = std::max(taylorBuildPeakBytes, combinedPeak);
                    result.taylorMemoryBytes = std::max(result.taylorMemoryBytes, combinedPeak);

                    bool accepted = built;
                    if (accepted && certifiedTaylorCapability(result.capability) && (!certifiedPiecewiseCandidate || jet.argCompositionCount > 0 || jet.polarCompositionCount > 0) && !request.allowUncertifiedForBenchmark && jet.landingIteration < request.maxIterations) {
                        accepted = false;
                        jet.status = ExpressionTaylorJetStatus::NoCoverage;
                        jet.failureReason = "certified Taylor jet does not cover the full iteration horizon";
                    }
                    if (accepted && request.taylor.requirePredictedBenefit && !predictedTaylorBenefit(jet, pixels, request.maxIterations)) {
                        accepted = false;
                        jet.status = ExpressionTaylorJetStatus::NoCoverage;
                        jet.failureReason = "Taylor cost predicts no tile benefit";
                    }
                    if (accepted && !taylorPersistentBytes(jet, persistentBytes)) {
                        accepted = false;
                        jet.status = ExpressionTaylorJetStatus::ResourceLimit;
                        jet.failureReason = "Taylor coefficient memory calculation overflow";
                    }
                    size_t retainedAfter = retainedBytes;
                    if (accepted && !checkedAddSize(retainedAfter, 1, persistentBytes)) {
                        accepted = false;
                        jet.status = ExpressionTaylorJetStatus::ResourceLimit;
                        jet.failureReason = "Taylor retained memory calculation overflow";
                    }
                    if (accepted && request.taylor.maximumJetMemoryBytes != 0 && retainedAfter > request.taylor.maximumJetMemoryBytes) {
                        accepted = false;
                        jet.status = ExpressionTaylorJetStatus::ResourceLimit;
                        jet.failureReason = "Taylor retained coefficients exceed jet memory policy";
                    }
                    size_t fastPeak = rendererBytes;
                    if (accepted && !checkedAddSize(fastPeak, 1, retainedAfter)) {
                        accepted = false;
                        jet.status = ExpressionTaylorJetStatus::ResourceLimit;
                        jet.failureReason = "Taylor renderer peak memory calculation overflow";
                    }
                    if (accepted && request.memory.memoryLimitBytes != 0 && (reference.memoryBytes > request.memory.memoryLimitBytes || fastPeak > request.memory.memoryLimitBytes - reference.memoryBytes)) {
                        accepted = false;
                        jet.status = ExpressionTaylorJetStatus::ResourceLimit;
                        jet.failureReason = "Taylor coefficients exceed renderer memory policy";
                    }
                    recordJet(jet, accepted, pixels, depth);
                    return accepted;
                };

                ScaledRealValue maximumRealMagnitude;
                ScaledRealValue maximumImaginaryMagnitude;
                ScaledRealValue maximumComponentError;
                for (const ScaledOffset& offset : xOffsets) {
                    const ScaledRealValue magnitude = absoluteScaled(offset.value);
                    if (compareScaledNonnegative(magnitude, maximumRealMagnitude) > 0) maximumRealMagnitude = magnitude;
                    if (compareScaledNonnegative(offset.error, maximumComponentError) > 0) maximumComponentError = offset.error;
                }
                for (const ScaledOffset& offset : yOffsets) {
                    const ScaledRealValue magnitude = absoluteScaled(offset.value);
                    if (compareScaledNonnegative(magnitude, maximumImaginaryMagnitude) > 0) maximumImaginaryMagnitude = magnitude;
                    if (compareScaledNonnegative(offset.error, maximumComponentError) > 0) maximumComponentError = offset.error;
                }
                ScaledComplexValue frameScale;
                const bool haveFrameScale = makeExpressionTaylorFrameScale(maximumRealMagnitude, maximumImaginaryMagnitude, maximumComponentError, frameScale);
                size_t rootPersistentBytes = 0;
                if (haveFrameScale) {
                    ScaledComplexBall zeroOffset;
                    useTaylor = attemptJet(zeroOffset, frameScale, pixelCount64, 0, 0, taylorJet, rootPersistentBytes);
                } else {
                    ExpressionTaylorJetResult failed;
                    failed.status = ExpressionTaylorJetStatus::ExponentRange;
                    failed.failureReason = "Taylor frame normalization failed";
                    recordJet(failed, false, pixelCount64, 0);
                }
                if (cancelled.load(std::memory_order_acquire) || taylorJet.status == ExpressionTaylorJetStatus::Cancelled) return fail(ExpressionDeepRenderStatus::Cancelled, "render cancelled while building Taylor jet");
                if (useTaylor) {
                    useGlobalTaylor = true;
                    retainedTaylorBytes = rootPersistentBytes;
                    TaylorTileJet rootTile;
                    rootTile.xEnd = request.width;
                    rootTile.yEnd = request.height;
                    rootTile.jet.landingIteration = taylorJet.landingIteration;
                    rootTile.jet.order = taylorJet.order;
                    result.taylorTileMapHash = mixTaylorTileHash(1469598103934665603ULL, rootTile);
                } else if (request.taylor.enableTileTaylor && request.taylor.maximumDepth > 0 && taylorJet.failureReason != "Taylor cost predicts no tile benefit" && result.taylorAttemptedJetCount < request.taylor.maximumJetCount) {
                    std::vector<ScaledComplexValue>().swap(taylorJet.coefficients);
                    std::vector<ScaledRealValue>().swap(taylorJet.coefficientRadii);
                    std::vector<ScaledRealValue>().swap(taylorJet.intermediateEscapeMargins);
                    rootPersistentBytes = 0;
                    std::vector<PendingTaylorTile> pending;
                    bool plannerReady = true;
                    size_t plannerMetadataBytes = 0;
                    if (!checkedAddSize(plannerMetadataBytes, pixelCount, sizeof(int32_t)) || !checkedAddSize(plannerMetadataBytes, request.taylor.maximumJetCount, sizeof(TaylorTileJet)) || !checkedAddSize(plannerMetadataBytes, request.taylor.maximumJetCount, sizeof(PendingTaylorTile))) plannerReady = false;
                    if (plannerReady && request.taylor.maximumJetMemoryBytes != 0 && plannerMetadataBytes > request.taylor.maximumJetMemoryBytes) plannerReady = false;
                    if (plannerReady && request.memory.memoryLimitBytes != 0) {
                        size_t fixed = reference.memoryBytes;
                        if (!checkedAddSize(fixed, 1, rendererBytes)) { plannerReady = false; }
                        if (plannerReady && (fixed > request.memory.memoryLimitBytes || plannerMetadataBytes > request.memory.memoryLimitBytes - fixed)) plannerReady = false;
                    }
                    if (plannerReady) {
                        try {
                            taylorPixelJet.assign(pixelCount, int32_t{-1});
                            taylorTileJets.reserve(request.taylor.maximumJetCount);
                            pending.reserve(request.taylor.maximumJetCount);
                        } catch (const std::bad_alloc&) { plannerReady = false; } catch (const std::length_error&) {
                            plannerReady = false;
                        }
                    }
                    if (plannerReady) {
                        size_t actualMetadataBytes = 0;
                        if (!checkedAddSize(actualMetadataBytes, taylorPixelJet.capacity(), sizeof(int32_t)) || !checkedAddSize(actualMetadataBytes, taylorTileJets.capacity(), sizeof(TaylorTileJet)) || !checkedAddSize(actualMetadataBytes, pending.capacity(), sizeof(PendingTaylorTile))) {
                            plannerReady = false;
                        } else {
                            plannerMetadataBytes = actualMetadataBytes;
                        }
                        if (plannerReady && request.taylor.maximumJetMemoryBytes != 0 && plannerMetadataBytes > request.taylor.maximumJetMemoryBytes) plannerReady = false;
                        if (plannerReady && request.memory.memoryLimitBytes != 0) {
                            size_t fixed = reference.memoryBytes;
                            if (!checkedAddSize(fixed, 1, rendererBytes) || fixed > request.memory.memoryLimitBytes || plannerMetadataBytes > request.memory.memoryLimitBytes - fixed) plannerReady = false;
                        }
                    }
                    if (plannerReady) {
                        retainedTaylorBytes = plannerMetadataBytes;
                        PendingTaylorTile root;
                        root.xEnd = request.width;
                        root.yEnd = request.height;
                        std::array<PendingTaylorTile, 4> children;
                        const size_t childCount = splitTaylorTile(root, request.taylor, children);
                        if (childCount != 0) {
                            ++result.taylorTileSplitCount;
                            const size_t available = request.taylor.maximumJetCount - result.taylorAttemptedJetCount;
                            const size_t retainedChildren = std::min(childCount, available);
                            for (size_t child = retainedChildren; child > 0; --child) pending.push_back(children[child - 1]);
                        }
                        bool stopForBudget = false;
                        uint64_t tileMapHash = 1469598103934665603ULL;
                        while (!pending.empty() && !stopForBudget && result.taylorAttemptedJetCount < request.taylor.maximumJetCount) {
                            if (pollCancellation()) break;
                            const PendingTaylorTile tile = pending.back();
                            pending.pop_back();
                            ScaledComplexBall parameterOffset;
                            ScaledComplexValue parameterScale;
                            ExpressionTaylorJetResult jet;
                            size_t persistentBytes = 0;
                            bool accepted = false;
                            if (makeTaylorTileParameterization(xOffsets, yOffsets, tile, frameScale, parameterOffset, parameterScale)) {
                                accepted = attemptJet(parameterOffset, parameterScale, taylorTilePixels(tile), tile.depth, retainedTaylorBytes, jet, persistentBytes);
                            } else {
                                jet.status = ExpressionTaylorJetStatus::ExponentRange;
                                jet.failureReason = "Taylor tile normalization failed";
                                recordJet(jet, false, taylorTilePixels(tile), tile.depth);
                            }
                            if (cancelled.load(std::memory_order_acquire) || jet.status == ExpressionTaylorJetStatus::Cancelled) break;
                            if (accepted) {
                                TaylorTileJet acceptedTile;
                                acceptedTile.xBegin = tile.xBegin;
                                acceptedTile.yBegin = tile.yBegin;
                                acceptedTile.xEnd = tile.xEnd;
                                acceptedTile.yEnd = tile.yEnd;
                                acceptedTile.depth = tile.depth;
                                acceptedTile.jet = std::move(jet);
                                const int32_t jetIndex = static_cast<int32_t>(taylorTileJets.size());
                                taylorTileJets.push_back(std::move(acceptedTile));
                                retainedTaylorBytes += persistentBytes;
                                const TaylorTileJet& retained = taylorTileJets.back();
                                tileMapHash = mixTaylorTileHash(tileMapHash, retained);
                                for (int y = tile.yBegin; y < tile.yEnd; ++y)
                                    for (int x = tile.xBegin; x < tile.xEnd; ++x) taylorPixelJet[static_cast<size_t>(y) * request.width + x] = jetIndex;
                            } else {
                                if (result.taylorAcceptedJetCount == 0 && result.taylorRejectedJetCount >= request.taylor.maximumRejectedBeforeFirstAcceptance) {
                                    pending.clear();
                                    continue;
                                }
                                if (jet.status == ExpressionTaylorJetStatus::ResourceLimit) {
                                    stopForBudget = true;
                                    continue;
                                }
                                if (jet.failureReason == "Taylor cost predicts no tile benefit") continue;
                                std::array<PendingTaylorTile, 4> children;
                                const size_t childCount = splitTaylorTile(tile, request.taylor, children);
                                if (childCount != 0 && result.taylorAttemptedJetCount + pending.size() < request.taylor.maximumJetCount) {
                                    ++result.taylorTileSplitCount;
                                    const size_t available = request.taylor.maximumJetCount - result.taylorAttemptedJetCount - pending.size();
                                    const size_t retainedChildren = std::min(childCount, available);
                                    for (size_t child = retainedChildren; child > 0; --child) pending.push_back(children[child - 1]);
                                }
                            }
                        }
                        if (cancelled.load(std::memory_order_acquire)) return fail(ExpressionDeepRenderStatus::Cancelled, "render cancelled while planning Taylor tiles");
                        useTaylor = !taylorTileJets.empty();
                        if (useTaylor) result.taylorTileMapHash = tileMapHash;
                    }
                    if (!plannerReady || taylorTileJets.empty()) {
                        std::vector<TaylorTileJet>().swap(taylorTileJets);
                        std::vector<int32_t>().swap(taylorPixelJet);
                        retainedTaylorBytes = 0;
                    }
                }
                if (result.taylorAcceptedPixelCoverage != 0) result.taylorWeightedLanding = static_cast<double>(weightedLandingTotal / static_cast<long double>(result.taylorAcceptedPixelCoverage));
                result.taylorAccepted = useTaylor;
                result.taylorRetainedBytes = retainedTaylorBytes;
                if (useTaylor) {
                    size_t fastWithTaylor = rendererBytes;
                    size_t buildWithTaylor = rendererBaseBytes;
                    if (!checkedAddSize(fastWithTaylor, 1, retainedTaylorBytes) || !checkedAddSize(buildWithTaylor, 1, taylorBuildPeakBytes)) {
                        useTaylor = false;
                        result.taylorFailureReason = "Taylor peak memory calculation overflow";
                    }
                    const size_t withTaylor = std::max(fastWithTaylor, buildWithTaylor);
                    if (useTaylor && request.memory.memoryLimitBytes != 0 && (reference.memoryBytes > request.memory.memoryLimitBytes || withTaylor > request.memory.memoryLimitBytes - reference.memoryBytes)) {
                        useTaylor = false;
                        result.taylorFailureReason = "Taylor coefficients exceed renderer memory policy";
                    } else if (useTaylor) {
                        rendererBytes = withTaylor;
                        result.rendererBytes = rendererBytes;
                    }
                }
                if (!useTaylor) {
                    std::vector<ScaledComplexValue>().swap(taylorJet.coefficients);
                    std::vector<ScaledRealValue>().swap(taylorJet.coefficientRadii);
                    std::vector<ScaledRealValue>().swap(taylorJet.intermediateEscapeMargins);
                    std::vector<TaylorTileJet>().swap(taylorTileJets);
                    std::vector<int32_t>().swap(taylorPixelJet);
                    retainedTaylorBytes = 0;
                    useGlobalTaylor = false;
                }
                result.taylorAccepted = useTaylor;
                result.taylorRetainedBytes = retainedTaylorBytes;
                result.taylorPlanningSeconds = std::max(0.0, secondsSince(taylorPlanningStart) - (result.taylorBuildSeconds - taylorBuildSecondsBefore));
            }
            if (runFast && request.preflight.enable && !result.preflightAttempted && piecewisePerStepEligible && !useTaylor) {
                validationEvaluator = std::make_unique<ExpressionScaledResidualEvaluator>();
                if (validationEvaluator->prepare(*request.runtimeProgram, reference) && validationEvaluator->ready()) {
                    if (!runCertifiedPreflight()) return false;
                }
                validationEvaluator.reset();
            }
            if (runFast && piecewisePerStepEligible && result.preflightAttempted && !result.preflightRejectedFast && !useTaylor && result.preflightSampleCount >= std::min<uint64_t>(request.preflight.minimumSamples, pixelCount64) && result.preflightFallbackCount == result.preflightSampleCount && result.preflightFastCount == 0) {
                if (!rejectFastFromPreflight()) return false;
            }
            if (runFast && certifiedTaylorCapability(result.capability) && !piecewisePerStepEligible && !request.allowUncertifiedForBenchmark && !useTaylor) {
                runFast = false;
                certificationUnavailable = true;
                reference = {};
                result.referenceBytes = 0;
                std::vector<ScaledOffset>().swap(xOffsets);
                std::vector<ScaledOffset>().swap(yOffsets);
                retainedTaylorBytes = 0;
                rendererBaseBytes = 0;
                if (!checkedAddSize(rendererBaseBytes, pixelCount, sizeof(uint8_t)) || !checkedAddSize(rendererBaseBytes, pixelCount, sizeof(size_t))) return fail(ExpressionDeepRenderStatus::ResourceLimit, "fallback renderer memory calculation overflow");
                fastThreadBytesTotal = 0;
                rendererBytes = rendererBaseBytes;
                result.rendererBytes = rendererBytes;
                std::fill(fallbackReason.begin(), fallbackReason.end(), static_cast<uint8_t>(ExpressionDeepFallbackReason::CertificationFailure));
            }
        }
        if (!runFast && fallbackReason[0] == NoFallbackReason) {
            const uint8_t reason = static_cast<uint8_t>(certificationUnavailable ? ExpressionDeepFallbackReason::CertificationFailure : reasonForCapability(result.capability));
            std::fill(fallbackReason.begin(), fallbackReason.end(), reason);
        }

        notifyProgress(ExpressionDeepRenderPhase::Fast, 0, pixelCount64);
        const Clock::time_point fastStart = Clock::now();
        const bool useBatchTaylor = useGlobalTaylor && request.taylor.enableBatchEvaluation && result.capability != ExpressionScaledResidualCapability::ExactCenteredArithmetic && taylorJet.layout == ExpressionTaylorJetLayout::ComplexUnivariate;
        std::atomic<uint64_t> fastPixels{0};
        std::atomic<uint64_t> totalIterations{0};
        std::atomic<uint64_t> taylorAcceptedPixels{0};
        std::atomic<uint64_t> taylorFallbackPixels{0};
        std::atomic<uint64_t> taylorEvaluationNanoseconds{0};
        std::atomic<uint64_t> taylorResidualNanoseconds{0};
        std::atomic<uint64_t> fastOperationCount{0};
        std::atomic<uint64_t> fastSeriesOperationCount{0};
        std::atomic<uint64_t> fastFoldOperationCount{0};
        std::atomic<uint64_t> fastUncertainFoldCount{0};
        std::array<std::atomic<uint64_t>, ExpressionDeepUncertaintyHistogramBins> firstUncertainHistogram{};
        for (auto& count : firstUncertainHistogram) count.store(0, std::memory_order_relaxed);
        std::atomic<uint64_t> fastCompleted{0};
        std::atomic_bool fastResourceError{false};
        std::atomic_bool fastInternalError{false};
        std::atomic_bool fastIterationFaultInjected{false};
        if (runFast) {
#pragma omp parallel num_threads(threadCount)
            {
                ExpressionScaledResidualEvaluator evaluator;
                bool prepared = false;
                bool workerReady = true;
                try {
                    if (request.verificationFault == ExpressionDeepVerificationFault::FastWorkerAllocation) throw std::bad_alloc();
                    prepared = evaluator.prepare(*request.runtimeProgram, reference);
                } catch (const std::bad_alloc&) {
                    fastResourceError.store(true, std::memory_order_release);
                    workerReady = false;
                } catch (const std::length_error&) {
                    fastResourceError.store(true, std::memory_order_release);
                    workerReady = false;
                } catch (...) {
                    fastInternalError.store(true, std::memory_order_release);
                    workerReady = false;
                }
                uint64_t localIterations = 0;
                uint64_t localFastPixels = 0;
                uint64_t localTaylorAcceptedPixels = 0;
                uint64_t localTaylorFallbackPixels = 0;
                uint64_t localTaylorEvaluationNanoseconds = 0;
                uint64_t localTaylorResidualNanoseconds = 0;
                uint64_t localFastOperationCount = 0;
                uint64_t localFastSeriesOperationCount = 0;
                uint64_t localFastFoldOperationCount = 0;
                uint64_t localFastUncertainFoldCount = 0;
                std::array<uint64_t, ExpressionDeepUncertaintyHistogramBins> localFirstUncertainHistogram{};
                TaylorBatchCache batchCache;
#pragma omp for schedule(dynamic, 1)
                for (long long tile = 0; tile < static_cast<long long>(tileCount); ++tile) {
                    try {
                        if (!workerReady || fastResourceError.load(std::memory_order_acquire) || fastInternalError.load(std::memory_order_acquire)) continue;
                        if (request.verificationFault == ExpressionDeepVerificationFault::FastIterationAllocation && !fastIterationFaultInjected.exchange(true, std::memory_order_acq_rel)) throw std::bad_alloc();
                        if (pollCancellation()) continue;
                        const int tileX = static_cast<int>(tile % static_cast<long long>(tilesX));
                        const int tileY = static_cast<int>(tile / static_cast<long long>(tilesX));
                        const int xBegin = tileX * request.threading.tileWidth;
                        const int yBegin = tileY * request.threading.tileHeight;
                        const int xEnd = std::min(request.width, xBegin + request.threading.tileWidth);
                        const int yEnd = std::min(request.height, yBegin + request.threading.tileHeight);
                        uint64_t completedInTile = 0;
                        for (int y = yBegin; y < yEnd; ++y) {
                            for (int x = xBegin; x < xEnd; ++x) {
                                if (pollCancellation()) break;
                                const size_t index = static_cast<size_t>(y) * request.width + x;
                                FastPixelResult pixel;
                                if (!prepared) {
                                    pixel.reason = ExpressionDeepFallbackReason::InvalidTape;
                                } else {
                                    ScaledComplexBall offset;
                                    const ScaledOffset& xOffset = xOffsets[static_cast<size_t>(x)];
                                    const ScaledOffset& yOffset = yOffsets[static_cast<size_t>(y)];
                                    offset.value.re = xOffset.value;
                                    offset.value.im = yOffset.value;
                                    offset.radius = compareScaledNonnegative(xOffset.error, yOffset.error) >= 0 ? xOffset.error : yOffset.error;
                                    ExpressionScaledResidualInput input;
                                    ScaledComplexBall stateDelta;
                                    ScaledArithmeticStatus arithmetic = initialReferenceExponentUnsafe ? ScaledArithmeticStatus::ExponentRange : certifyScaledMpfrExponentRange(offset);
                                    if (arithmetic != ScaledArithmeticStatus::Success) {
                                        pixel.reason = reasonForArithmeticStatus(arithmetic);
                                    } else if (request.pixelParameter == FormulaParameter::C) {
                                        input.c = offset.value;
                                        input.cError = offset.radius;
                                        stateDelta.radius = initialReference.radius;
                                    } else {
                                        input.z0 = offset.value;
                                        input.z0Error = offset.radius;
                                        stateDelta = offset;
                                    }

                                    ScaledComplexBall initialBase;
                                    initialBase.value = initialReference.value;
                                    ScaledComplexBall initialValue;
                                    if (pixel.reason == ExpressionDeepFallbackReason::InvalidTape) arithmetic = certifyScaledMpfrExponentRange(initialBase);
                                    if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifyScaledMpfrExponentRange(stateDelta);
                                    if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifiedScaledAdd(initialBase, stateDelta, initialValue);
                                    if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifyScaledMpfrExponentRange(initialValue);
                                    if (arithmetic != ScaledArithmeticStatus::Success) {
                                        pixel.reason = reasonForArithmeticStatus(arithmetic);
                                    } else {
                                        ScaledArithmeticStatus gateStatus = ScaledArithmeticStatus::Success;
                                        BailoutDecision decision = decideBailout(initialValue.value, initialValue.radius, bailoutThreshold, gateStatus);
                                        if (decision == BailoutDecision::Error) {
                                            pixel.reason = reasonForArithmeticStatus(gateStatus);
                                        } else if (decision == BailoutDecision::Uncertain) {
                                            pixel.reason = ExpressionDeepFallbackReason::BailoutUncertain;
                                        } else if (decision == BailoutDecision::Outside) {
                                            pixel.decided = true;
                                            pixel.output = 0.0f;
                                        } else {
                                            int firstIteration = 0;
                                            const ExpressionTaylorJetResult* pixelTaylorJet = nullptr;
                                            if (useGlobalTaylor) {
                                                pixelTaylorJet = &taylorJet;
                                            } else if (!taylorPixelJet.empty()) {
                                                const int32_t mappedJet = taylorPixelJet[index];
                                                if (mappedJet >= 0 && static_cast<size_t>(mappedJet) < taylorTileJets.size()) pixelTaylorJet = &taylorTileJets[static_cast<size_t>(mappedJet)].jet;
                                            }
                                            if (pixelTaylorJet) {
                                                ScaledComplexBall q;
                                                ScaledComplexBall landingDelta;
                                                bool landed = false;
                                                ExpressionTaylorJetEvaluation landing;
                                                bool evaluatedLanding = false;
                                                bool batchedLanding = false;
                                                if (useBatchTaylor) {
                                                    batchedLanding = getTaylorBatchLanding(batchCache, xOffsets, yOffset, x, y, xEnd, *pixelTaylorJet, landing, localTaylorEvaluationNanoseconds);
                                                    evaluatedLanding = batchedLanding;
                                                }
                                                Clock::time_point evaluationStart;
                                                if (!batchedLanding) { evaluationStart = Clock::now(); }
                                                if (!evaluatedLanding && makeExpressionTaylorLocalQ(offset, pixelTaylorJet->parameterOffset, pixelTaylorJet->parameterScale, q)) { evaluatedLanding = ExpressionTaylorJetEvaluator::evaluate(*pixelTaylorJet, q, landing); }
                                                if (evaluatedLanding) {
                                                    landingDelta = landing.residual;
                                                    const ExpressionReferenceSample& landingSample = reference.samples[pixelTaylorJet->landingSample];
                                                    ScaledComplexBall landingBase;
                                                    ScaledComplexBall landingValue;
                                                    arithmetic = makeScaledComplexValue(pixelTaylorJet->landingUsesSampleOutput ? landingSample.next : landingSample.z, pixelTaylorJet->landingUsesSampleOutput ? landingSample.rootDefect : landingSample.zDefect, landingBase.value);
                                                    if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifiedScaledAdd(landingBase, landingDelta, landingValue);
                                                    if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifyScaledMpfrExponentRange(landingValue);
                                                    if (arithmetic == ScaledArithmeticStatus::Success) {
                                                        decision = decideBailout(landingValue.value, landingValue.radius, bailoutThreshold, gateStatus);
                                                        landed = decision == BailoutDecision::Inside;
                                                    }
                                                }
                                                if (!batchedLanding) localTaylorEvaluationNanoseconds += static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - evaluationStart).count());
                                                if (landed) {
                                                    stateDelta = landingDelta;
                                                    firstIteration = pixelTaylorJet->landingIteration;
                                                    ++localTaylorAcceptedPixels;
                                                    if (firstIteration == request.maxIterations) {
                                                        pixel.decided = true;
                                                        pixel.output = ExpressionDeepInteriorPixel;
                                                    }
                                                } else {
                                                    ++localTaylorFallbackPixels;
                                                    if (certifiedTaylorCapability(result.capability) && !certifiedPiecewiseCandidate && !request.allowUncertifiedForBenchmark) {
                                                        pixel.reason = ExpressionDeepFallbackReason::CertificationFailure;
                                                        firstIteration = request.maxIterations;
                                                    }
                                                }
                                            } else if (useTaylor) {
                                                ++localTaylorFallbackPixels;
                                                if (certifiedTaylorCapability(result.capability) && !piecewisePerStepEligible && !request.allowUncertifiedForBenchmark) {
                                                    pixel.reason = ExpressionDeepFallbackReason::CertificationFailure;
                                                    firstIteration = request.maxIterations;
                                                }
                                            }
                                            const bool timeTaylorResidual = pixelTaylorJet && firstIteration < request.maxIterations;
                                            Clock::time_point residualStart;
                                            if (timeTaylorResidual) residualStart = Clock::now();
                                            for (int iteration = firstIteration; iteration < request.maxIterations; ++iteration) {
                                                pixel.firstUncertainIteration = static_cast<uint32_t>(iteration);
                                                if (pollCancellation()) break;
                                                if (static_cast<size_t>(iteration) >= reference.samples.size()) {
                                                    pixel.reason = ExpressionDeepFallbackReason::ReferenceExhausted;
                                                    break;
                                                }
                                                input.z = stateDelta.value;
                                                input.zError = stateDelta.radius;
                                                input.iteration = iteration;
                                                ExpressionScaledResidualResult evaluated = evaluator.evaluate(static_cast<size_t>(iteration), input);
                                                ++pixel.iterations;
                                                pixel.operationCount = saturatingAdd(pixel.operationCount, static_cast<uint64_t>(evaluated.operationCount));
                                                pixel.seriesOperationCount = saturatingAdd(pixel.seriesOperationCount, static_cast<uint64_t>(evaluated.seriesOperationCount));
                                                pixel.foldOperationCount = saturatingAdd(pixel.foldOperationCount, static_cast<uint64_t>(evaluated.foldOperationCount));
                                                pixel.uncertainFoldCount = saturatingAdd(pixel.uncertainFoldCount, static_cast<uint64_t>(evaluated.uncertainFoldCount));
                                                if (evaluated.status != ExpressionScaledResidualStatus::Success) {
                                                    pixel.reason = reasonForResidualStatus(evaluated.status);
                                                    break;
                                                }
                                                if (evaluated.uncertified && !request.allowUncertifiedForBenchmark) {
                                                    pixel.reason = ExpressionDeepFallbackReason::UncertifiedSeries;
                                                    break;
                                                }
                                                if ((result.capability == ExpressionScaledResidualCapability::ExactCenteredArithmetic || certifiedPiecewiseCandidate) && !evaluated.certified) {
                                                    pixel.reason = ExpressionDeepFallbackReason::CertificationFailure;
                                                    break;
                                                }
                                                const ExpressionReferenceSample& sample = reference.samples[static_cast<size_t>(iteration)];
                                                ScaledComplexBall outputBase;
                                                ScaledComplexBall residualBall;
                                                ScaledComplexBall actualOutput;
                                                arithmetic = makeScaledComplexValue(sample.next, sample.rootDefect, outputBase.value);
                                                if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifyScaledMpfrExponentRange(outputBase);
                                                residualBall.value = evaluated.residual;
                                                residualBall.radius = evaluated.radius;
                                                if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifiedScaledAdd(outputBase, residualBall, actualOutput);
                                                if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifyScaledMpfrExponentRange(actualOutput);
                                                if (arithmetic != ScaledArithmeticStatus::Success) {
                                                    pixel.reason = reasonForArithmeticStatus(arithmetic);
                                                    break;
                                                }
                                                decision = decideBailout(actualOutput.value, actualOutput.radius, bailoutThreshold, gateStatus);
                                                if (decision == BailoutDecision::Error) {
                                                    pixel.reason = reasonForArithmeticStatus(gateStatus);
                                                    break;
                                                }
                                                if (decision == BailoutDecision::Uncertain) {
                                                    pixel.reason = ExpressionDeepFallbackReason::BailoutUncertain;
                                                    break;
                                                }
                                                if (decision == BailoutDecision::Outside) {
                                                    pixel.decided = true;
                                                    pixel.output = static_cast<float>(iteration + 1);
                                                    break;
                                                }
                                                if (iteration + 1 == request.maxIterations) {
                                                    pixel.decided = true;
                                                    pixel.output = ExpressionDeepInteriorPixel;
                                                    break;
                                                }
                                                if (static_cast<size_t>(iteration + 1) >= reference.samples.size()) {
                                                    pixel.reason = ExpressionDeepFallbackReason::ReferenceExhausted;
                                                    break;
                                                }
                                                ScaledComplexBall nextBase;
                                                arithmetic = makeScaledComplexValue(reference.samples[static_cast<size_t>(iteration + 1)].z, reference.samples[static_cast<size_t>(iteration + 1)].zDefect, nextBase.value);
                                                if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifyScaledMpfrExponentRange(nextBase);
                                                if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifiedScaledSubtract(actualOutput, nextBase, stateDelta);
                                                if (arithmetic == ScaledArithmeticStatus::Success) arithmetic = certifyScaledMpfrExponentRange(stateDelta);
                                                if (arithmetic != ScaledArithmeticStatus::Success) {
                                                    pixel.reason = reasonForArithmeticStatus(arithmetic);
                                                    break;
                                                }
                                            }
                                            if (timeTaylorResidual) localTaylorResidualNanoseconds += static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - residualStart).count());
                                        }
                                    }
                                }
                                localIterations += pixel.iterations;
                                localFastOperationCount = saturatingAdd(localFastOperationCount, pixel.operationCount);
                                localFastSeriesOperationCount = saturatingAdd(localFastSeriesOperationCount, pixel.seriesOperationCount);
                                localFastFoldOperationCount = saturatingAdd(localFastFoldOperationCount, pixel.foldOperationCount);
                                localFastUncertainFoldCount = saturatingAdd(localFastUncertainFoldCount, pixel.uncertainFoldCount);
                                if (cancelled.load(std::memory_order_acquire)) break;
                                if (pixel.decided) {
                                    if (!cancelled.load()) {
                                        request.output[index] = pixel.output;
                                        if (cancelled.load()) {
                                            request.output[index] = ExpressionDeepEmptyPixel;
                                        } else {
                                            ++localFastPixels;
                                        }
                                    }
                                } else {
                                    fallbackReason[index] = static_cast<uint8_t>(pixel.reason);
                                    ++localFirstUncertainHistogram[uncertaintyHistogramBin(pixel.firstUncertainIteration)];
                                }
                                ++completedInTile;
                            }
                            if (cancelled.load(std::memory_order_acquire)) break;
                        }
                        const uint64_t completed = fastCompleted.fetch_add(completedInTile, std::memory_order_relaxed) + completedInTile;
                        notifyProgress(ExpressionDeepRenderPhase::Fast, completed, pixelCount64);
                    } catch (const std::bad_alloc&) { fastResourceError.store(true, std::memory_order_release); } catch (const std::length_error&) {
                        fastResourceError.store(true, std::memory_order_release);
                    } catch (...) { fastInternalError.store(true, std::memory_order_release); }
                }
                totalIterations.fetch_add(localIterations, std::memory_order_relaxed);
                fastPixels.fetch_add(localFastPixels, std::memory_order_relaxed);
                taylorAcceptedPixels.fetch_add(localTaylorAcceptedPixels, std::memory_order_relaxed);
                taylorFallbackPixels.fetch_add(localTaylorFallbackPixels, std::memory_order_relaxed);
                taylorEvaluationNanoseconds.fetch_add(localTaylorEvaluationNanoseconds, std::memory_order_relaxed);
                taylorResidualNanoseconds.fetch_add(localTaylorResidualNanoseconds, std::memory_order_relaxed);
                fastOperationCount.fetch_add(localFastOperationCount, std::memory_order_relaxed);
                fastSeriesOperationCount.fetch_add(localFastSeriesOperationCount, std::memory_order_relaxed);
                fastFoldOperationCount.fetch_add(localFastFoldOperationCount, std::memory_order_relaxed);
                fastUncertainFoldCount.fetch_add(localFastUncertainFoldCount, std::memory_order_relaxed);
                for (size_t bin = 0; bin < localFirstUncertainHistogram.size(); ++bin) firstUncertainHistogram[bin].fetch_add(localFirstUncertainHistogram[bin], std::memory_order_relaxed);
            }
        }
        result.fastSeconds = secondsSince(fastStart);
        result.fastPixelCount = fastPixels.load(std::memory_order_relaxed);
        result.totalIterations = totalIterations.load(std::memory_order_relaxed);
        result.fastIterationCount = result.totalIterations;
        result.taylorAcceptedPixelCount = taylorAcceptedPixels.load(std::memory_order_relaxed);
        result.taylorFallbackPixelCount = taylorFallbackPixels.load(std::memory_order_relaxed);
        result.taylorEvaluationSeconds = static_cast<double>(taylorEvaluationNanoseconds.load(std::memory_order_relaxed)) / 1.0e9;
        result.taylorResidualSeconds = static_cast<double>(taylorResidualNanoseconds.load(std::memory_order_relaxed)) / 1.0e9;
        result.fastOperationCount = fastOperationCount.load(std::memory_order_relaxed);
        result.fastSeriesOperationCount = fastSeriesOperationCount.load(std::memory_order_relaxed);
        result.fastFoldOperationCount = fastFoldOperationCount.load(std::memory_order_relaxed);
        result.fastUncertainFoldCount = fastUncertainFoldCount.load(std::memory_order_relaxed);
        for (size_t bin = 0; bin < result.fallbackFirstUncertainHistogram.size(); ++bin) result.fallbackFirstUncertainHistogram[bin] = firstUncertainHistogram[bin].load(std::memory_order_relaxed);
        if (fastResourceError.load(std::memory_order_acquire)) return fail(ExpressionDeepRenderStatus::ResourceLimit, "fast worker allocation failed");
        if (fastInternalError.load(std::memory_order_acquire)) return fail(ExpressionDeepRenderStatus::InternalError, "fast worker failed");

        std::vector<size_t> fallbackQueue;
        fallbackQueue.reserve(pixelCount);
        for (size_t index = 0; index < pixelCount; ++index) {
            const uint8_t encoded = fallbackReason[index];
            if (encoded == NoFallbackReason) continue;
            fallbackQueue.push_back(index);
            const auto reason = static_cast<ExpressionDeepFallbackReason>(encoded);
            const size_t reasonIndex = static_cast<size_t>(reason);
            if (reasonIndex < result.fallbackReasonCounts.size()) ++result.fallbackReasonCounts[reasonIndex];
            if (reason == ExpressionDeepFallbackReason::UncertifiedSeries || reason == ExpressionDeepFallbackReason::BranchSensitive || reason == ExpressionDeepFallbackReason::CertificationFailure || reason == ExpressionDeepFallbackReason::BailoutUncertain) ++result.uncertainPixelCount;
        }
        result.fallbackPixelCount = static_cast<uint64_t>(fallbackQueue.size());
        for (uint64_t count : result.fallbackReasonCounts) result.maxFallbackReasonCount = std::max(result.maxFallbackReasonCount, count);

        for (uint64_t tile = 0; tile < tileCount; ++tile) {
            const int tileX = static_cast<int>(tile % tilesX);
            const int tileY = static_cast<int>(tile / tilesX);
            const int xBegin = tileX * request.threading.tileWidth;
            const int yBegin = tileY * request.threading.tileHeight;
            const int xEnd = std::min(request.width, xBegin + request.threading.tileWidth);
            const int yEnd = std::min(request.height, yBegin + request.threading.tileHeight);
            uint64_t fallbackInTile = 0;
            for (int y = yBegin; y < yEnd; ++y)
                for (int x = xBegin; x < xEnd; ++x) fallbackInTile += fallbackReason[static_cast<size_t>(y) * request.width + x] != NoFallbackReason;
            if (fallbackInTile == 0) continue;
            ++result.fallbackTileCount;
            const uint64_t pixelsInTile = static_cast<uint64_t>(xEnd - xBegin) * static_cast<uint64_t>(yEnd - yBegin);
            result.maxTileFallbackRate = std::max(result.maxTileFallbackRate, static_cast<double>(fallbackInTile) / static_cast<double>(pixelsInTile));
        }

        if (cancelled.load(std::memory_order_acquire)) return fail(ExpressionDeepRenderStatus::Cancelled, "render cancelled during the fast pass");

        const ExpressionPiecewiseQuadraticKind piecewiseQuadraticKind = request.disableSpecializedPiecewiseMpfrForVerification ? ExpressionPiecewiseQuadraticKind::None : request.runtimeProgram->piecewiseQuadraticKind();
        const char* periodicSetting = std::getenv("MANDEL_EXPR_PERIODIC");
        const bool fallbackPeriodicEnabled = (!periodicSetting || std::atoi(periodicSetting) != 0) && request.maxIterations >= 4096 && !request.runtimeProgram->iterationDependent();
        const bool useBigFixedPiecewise = [] {
            const char* value = std::getenv("MANDEL_DEEP_BIGFIXED_PIXELS");
            return !value || atoi(value) != 0;
        }() && piecewiseQuadraticKind == ExpressionPiecewiseQuadraticKind::BurningShip &&
                                          !request.disablePiecewiseBigFixedForVerification && result.fallbackPrecision >= 257 && result.fallbackPrecision <= 960 && request.bailout >= 2.0 && request.bailout <= 0x1p14 && std::floor(request.bailout) == request.bailout;
        const int bigFixedLimbCount = std::max(3, static_cast<int>((result.fallbackPrecision + 63) / 64 + 1));
        if (!fallbackQueue.empty()) {
            if (result.selectedPrecision > MPFR_PREC_MAX - request.memory.fallbackGuardBits) return fail(ExpressionDeepRenderStatus::PrecisionOutOfRange, "fallback precision overflow");
            result.fallbackPrecision = result.selectedPrecision + request.memory.fallbackGuardBits;
            if (result.fallbackPrecision > ExpressionReferencePrecisionPolicy::ApplicationMaximumBits || result.fallbackPrecision > request.precision.maximumBits) return fail(ExpressionDeepRenderStatus::PrecisionOutOfRange, "fallback precision exceeds policy maximum");
            const size_t limbs = (static_cast<size_t>(result.fallbackPrecision) + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
            size_t mpfrValueBytes = sizeof(__mpfr_struct);
            if (!checkedAddSize(mpfrValueBytes, limbs, sizeof(mp_limb_t))) return fail(ExpressionDeepRenderStatus::ResourceLimit, "MPFR workspace calculation overflow");
            // Worker state (34), oracle scratch (12), up to eight simultaneous
            // operation temporaries, and five reusable piecewise temporaries,
            // plus the bytecode stack.
            size_t mpfrValuesPerThread = 62;
            if (!checkedAddSize(mpfrValuesPerThread, request.runtimeProgram->stackDepth(), 2)) return fail(ExpressionDeepRenderStatus::ResourceLimit, "oracle stack calculation overflow");
            size_t fallbackThreadBytes = 0;
            if (!checkedAddSize(fallbackThreadBytes, mpfrValuesPerThread, mpfrValueBytes)) return fail(ExpressionDeepRenderStatus::ResourceLimit, "fallback workspace calculation overflow");
            if (useBigFixedPiecewise && !checkedAddSize(fallbackThreadBytes, 1, 2048 + static_cast<size_t>(bigFixedLimbCount) * 16 * sizeof(uint64_t))) return fail(ExpressionDeepRenderStatus::ResourceLimit, "BigFixed fallback workspace calculation overflow");
            size_t fallbackThreadBytesTotal = 0;
            if (!checkedAddSize(fallbackThreadBytesTotal, static_cast<size_t>(threadCount), fallbackThreadBytes)) return fail(ExpressionDeepRenderStatus::ResourceLimit, "fallback workspace multiplication overflow");
            size_t fallbackRendererBytes = rendererBaseBytes;
            if (!checkedAddSize(fallbackRendererBytes, 1, retainedTaylorBytes) || !checkedAddSize(fallbackRendererBytes, 1, std::max(fastThreadBytesTotal, fallbackThreadBytesTotal))) return fail(ExpressionDeepRenderStatus::ResourceLimit, "renderer peak memory calculation overflow");
            rendererBytes = std::max(rendererBytes, fallbackRendererBytes);
            result.rendererBytes = rendererBytes;
            if (request.memory.memoryLimitBytes != 0 && (reference.memoryBytes > request.memory.memoryLimitBytes || rendererBytes > request.memory.memoryLimitBytes - reference.memoryBytes)) return fail(ExpressionDeepRenderStatus::ResourceLimit, "MPFR fallback exceeds memory limit");
        }

        notifyProgress(ExpressionDeepRenderPhase::Fallback, 0, static_cast<uint64_t>(fallbackQueue.size()));
        const Clock::time_point fallbackStart = Clock::now();
        std::atomic<uint64_t> fallbackCompleted{0};
        std::atomic<uint64_t> undefinedPixels{0};
        std::atomic<uint64_t> specializedPiecewisePixels{0};
        std::atomic<uint64_t> specializedPiecewiseIterations{0};
        std::atomic<uint64_t> specializedPiecewisePeriodicPixels{0};
        std::atomic<uint64_t> genericMpfrPeriodicPixels{0};
        std::atomic<uint64_t> piecewiseBigFixedPixels{0};
        std::atomic<uint64_t> piecewiseBigFixedIterations{0};
        std::atomic_bool fallbackResourceError{false};
        std::atomic_bool fallbackInternalError{false};
        std::atomic_bool fallbackIterationFaultInjected{false};
        if (!fallbackQueue.empty()) {
#pragma omp parallel num_threads(threadCount)
            {
                std::unique_ptr<FallbackWorkerWorkspace> workspace;
                std::unique_ptr<BigFixedPixelWorkspace> bigFixedWorkspace;
                try {
                    if (request.verificationFault == ExpressionDeepVerificationFault::FallbackWorkerAllocation) throw std::bad_alloc();
                    workspace = std::make_unique<FallbackWorkerWorkspace>(result.fallbackPrecision, request);
                    if (useBigFixedPiecewise) bigFixedWorkspace = std::make_unique<BigFixedPixelWorkspace>(bigFixedLimbCount, request.bailout, result.fallbackPrecision);
                    if (!workspace->geometryReady) fallbackInternalError.store(true, std::memory_order_release);
                } catch (const std::bad_alloc&) { fallbackResourceError.store(true, std::memory_order_release); } catch (const std::length_error&) {
                    fallbackResourceError.store(true, std::memory_order_release);
                } catch (...) { fallbackInternalError.store(true, std::memory_order_release); }
                uint64_t localIterations = 0;
                uint64_t localSpecializedPixels = 0;
                uint64_t localSpecializedIterations = 0;
                uint64_t localSpecializedPeriodicPixels = 0;
                uint64_t localGenericPeriodicPixels = 0;
                uint64_t localBigFixedPixels = 0;
                uint64_t localBigFixedIterations = 0;
#pragma omp for schedule(dynamic, 8)
                for (long long queueIndex = 0; queueIndex < static_cast<long long>(fallbackQueue.size()); ++queueIndex) {
                    try {
                        if (!workspace || fallbackResourceError.load(std::memory_order_acquire) || fallbackInternalError.load(std::memory_order_acquire)) continue;
                        if (request.verificationFault == ExpressionDeepVerificationFault::FallbackIterationAllocation && !fallbackIterationFaultInjected.exchange(true, std::memory_order_acq_rel)) throw std::bad_alloc();
                        if (pollCancellation()) continue;
                        ExactGeometry& geometry = workspace->geometry;
                        ExpressionOracleContext& context = workspace->context;
                        MpfrComplex& pixel = workspace->pixel;
                        MpfrComplex& next = workspace->next;
                        mpfr_ptr magnitude = workspace->magnitudeStorage.re;
                        const size_t outputIndex = fallbackQueue[static_cast<size_t>(queueIndex)];
                        const int y = static_cast<int>(outputIndex / static_cast<size_t>(request.width));
                        const int x = static_cast<int>(outputIndex % static_cast<size_t>(request.width));
                        configureFixed(request.fixed, context);
                        bool undefined = !geometry.coordinate(x, y, request, pixel);
                        const bool specializedPiecewise = piecewiseQuadraticKind != ExpressionPiecewiseQuadraticKind::None;
                        if (!undefined) {
                            if (request.pixelParameter == FormulaParameter::C) {
                                context.c.set(pixel);
                            } else {
                                context.z0.set(pixel);
                            }
                            context.z.set(context.z0);
                        }

                        float output = ExpressionDeepInteriorPixel;
                        bool decided = false;
                        uint64_t periodicPower = 1;
                        uint64_t periodicLength = 0;
                        bool periodicReady = false;
                        bool usedBigFixed = false;
                        if (!undefined && bigFixedWorkspace) {
                            uint64_t bigFixedIterations = 0;
                            const BigFixedPixelStatus bigFixedStatus = evaluateBurningShipBigFixedPixel(request, context, *bigFixedWorkspace, output, bigFixedIterations, pollCancellation);
                            if (bigFixedStatus == BigFixedPixelStatus::Success) {
                                localIterations += bigFixedIterations;
                                localBigFixedIterations += bigFixedIterations;
                                ++localBigFixedPixels;
                                usedBigFixed = true;
                                decided = true;
                            } else if (bigFixedStatus == BigFixedPixelStatus::Cancelled) {
                                cancelled.store(true, std::memory_order_release);
                                continue;
                            }
                        }
                        if (!undefined && !decided) {
                            if (mpfr_nan_p(context.z.re) || mpfr_nan_p(context.z.im)) {
                                undefined = true;
                            } else if (mpfr_inf_p(context.z.re) || mpfr_inf_p(context.z.im)) {
                                output = 0.0f;
                                decided = true;
                            } else if (specializedPiecewise) {
                                bool outside = false;
                                if (!piecewiseOutsideBailout(context.z, workspace->piecewiseBailoutSquared, request.bailout, magnitude, workspace->piecewiseSquares, workspace->piecewiseInsideSquareExponent, workspace->piecewiseScratch, outside)) {
                                    undefined = true;
                                } else if (outside) {
                                    output = 0.0f;
                                    decided = true;
                                }
                            } else {
                                bool outside = false;
                                if (!piecewiseOutsideBailout(context.z, workspace->piecewiseBailoutSquared, request.bailout, magnitude, workspace->piecewiseSquares, workspace->piecewiseInsideSquareExponent, workspace->piecewiseScratch, outside)) {
                                    undefined = true;
                                } else if (outside) {
                                    output = 0.0f;
                                    decided = true;
                                }
                            }
                        }
                        if (!undefined && !decided && fallbackPeriodicEnabled) {
                            workspace->periodicState.set(context.z);
                            periodicReady = true;
                        }
                        for (int iteration = 0; !undefined && !decided && iteration < request.maxIterations; ++iteration) {
                            if ((iteration & 15) == 0 && pollCancellation()) break;
                            context.iteration = iteration;
                            bool oracleDefined = true;
                            if (specializedPiecewise) {
                                oracleDefined = evaluatePiecewiseQuadraticMpfr(piecewiseQuadraticKind, context, workspace->piecewiseSquares, next, workspace->piecewiseScratch);
                            } else {
                                oracleDefined = ExpressionOracle::evaluateOrbitStep(*request.runtimeProgram, context, next, &workspace->piecewiseSquares, nullptr);
                            }
                            ++localIterations;
                            if (specializedPiecewise) ++localSpecializedIterations;
                            if (mpfr_nan_p(next.re) || mpfr_nan_p(next.im)) {
                                undefined = true;
                                break;
                            }
                            if (mpfr_inf_p(next.re) || mpfr_inf_p(next.im)) {
                                output = static_cast<float>(iteration + 1);
                                decided = true;
                                break;
                            }
                            if (!oracleDefined || !mpfr_number_p(next.re) || !mpfr_number_p(next.im)) {
                                undefined = true;
                                break;
                            }
                            mpfr_swap(context.z.re, next.re);
                            mpfr_swap(context.z.im, next.im);
                            bool outside = false;
                            if (specializedPiecewise) {
                                if (!piecewiseOutsideBailout(context.z, workspace->piecewiseBailoutSquared, request.bailout, magnitude, workspace->piecewiseSquares, workspace->piecewiseInsideSquareExponent, workspace->piecewiseScratch, outside)) {
                                    undefined = true;
                                    break;
                                }
                            } else {
                                if (!piecewiseOutsideBailout(context.z, workspace->piecewiseBailoutSquared, request.bailout, magnitude, workspace->piecewiseSquares, workspace->piecewiseInsideSquareExponent, workspace->piecewiseScratch, outside)) {
                                    undefined = true;
                                    break;
                                }
                            }
                            if (outside) {
                                output = static_cast<float>(iteration + 1);
                                decided = true;
                            } else {
                                if (periodicReady) {
                                    ++periodicLength;
                                    if (sameMpfrStateComponent(context.z.re, workspace->periodicState.re) && sameMpfrStateComponent(context.z.im, workspace->periodicState.im)) {
                                        output = ExpressionDeepInteriorPixel;
                                        decided = true;
                                        if (specializedPiecewise)
                                            ++localSpecializedPeriodicPixels;
                                        else
                                            ++localGenericPeriodicPixels;
                                    } else if (periodicLength == periodicPower) {
                                        workspace->periodicState.set(context.z);
                                        periodicLength = 0;
                                        periodicPower *= 2;
                                    }
                                }
                                if (!decided && iteration + 1 == request.maxIterations) {
                                    output = ExpressionDeepInteriorPixel;
                                    decided = true;
                                }
                            }
                        }
                        if (cancelled.load(std::memory_order_acquire)) continue;
                        if (undefined) {
                            undefinedPixels.fetch_add(1, std::memory_order_relaxed);
                        } else if (decided) {
                            if (!cancelled.load()) {
                                request.output[outputIndex] = output;
                                if (cancelled.load()) request.output[outputIndex] = ExpressionDeepEmptyPixel;
                            }
                        }
                        if (specializedPiecewise && !usedBigFixed) ++localSpecializedPixels;
                        const uint64_t completed = fallbackCompleted.fetch_add(1, std::memory_order_relaxed) + 1;
                        if ((completed & 31) == 0) notifyProgress(ExpressionDeepRenderPhase::Fallback, completed, static_cast<uint64_t>(fallbackQueue.size()));
                    } catch (const std::bad_alloc&) { fallbackResourceError.store(true, std::memory_order_release); } catch (const std::length_error&) {
                        fallbackResourceError.store(true, std::memory_order_release);
                    } catch (...) { fallbackInternalError.store(true, std::memory_order_release); }
                }
                totalIterations.fetch_add(localIterations, std::memory_order_relaxed);
                specializedPiecewisePixels.fetch_add(localSpecializedPixels, std::memory_order_relaxed);
                specializedPiecewiseIterations.fetch_add(localSpecializedIterations, std::memory_order_relaxed);
                specializedPiecewisePeriodicPixels.fetch_add(localSpecializedPeriodicPixels, std::memory_order_relaxed);
                genericMpfrPeriodicPixels.fetch_add(localGenericPeriodicPixels, std::memory_order_relaxed);
                piecewiseBigFixedPixels.fetch_add(localBigFixedPixels, std::memory_order_relaxed);
                piecewiseBigFixedIterations.fetch_add(localBigFixedIterations, std::memory_order_relaxed);
                ExpressionOracle::releaseThreadWorkspace();
            }
        }
        result.fallbackSeconds = secondsSince(fallbackStart);
        result.totalIterations = totalIterations.load(std::memory_order_relaxed);
        result.undefinedPixelCount = undefinedPixels.load(std::memory_order_relaxed);
        result.specializedPiecewiseMpfrPixelCount = specializedPiecewisePixels.load(std::memory_order_relaxed);
        result.specializedPiecewiseMpfrIterationCount = specializedPiecewiseIterations.load(std::memory_order_relaxed);
        result.specializedPiecewiseMpfrPeriodicPixelCount = specializedPiecewisePeriodicPixels.load(std::memory_order_relaxed);
        result.genericMpfrPeriodicPixelCount = genericMpfrPeriodicPixels.load(std::memory_order_relaxed);
        result.usedSpecializedPiecewiseMpfr = result.specializedPiecewiseMpfrPixelCount != 0;
        result.piecewiseBigFixedPixelCount = piecewiseBigFixedPixels.load(std::memory_order_relaxed);
        result.piecewiseBigFixedIterationCount = piecewiseBigFixedIterations.load(std::memory_order_relaxed);
        result.usedPiecewiseBigFixed = result.piecewiseBigFixedPixelCount != 0;
        if (fallbackResourceError.load(std::memory_order_acquire)) return fail(ExpressionDeepRenderStatus::ResourceLimit, "fallback worker allocation failed");
        if (fallbackInternalError.load(std::memory_order_acquire)) return fail(ExpressionDeepRenderStatus::InternalError, "fallback geometry initialization failed");
        if (cancelled.load(std::memory_order_acquire)) return fail(ExpressionDeepRenderStatus::Cancelled, "render cancelled during MPFR fallback");
        if (result.undefinedPixelCount != 0) return fail(ExpressionDeepRenderStatus::UndefinedPixel, "one or more pixel orbits are undefined; their output remains empty");

        notifyProgress(ExpressionDeepRenderPhase::Complete, pixelCount64, pixelCount64);
        result.status = ExpressionDeepRenderStatus::Success;
        result.success = true;
        result.cancelled = false;
        result.error.clear();
        return true;
    } catch (const std::bad_alloc&) { return fail(ExpressionDeepRenderStatus::ResourceLimit, "renderer allocation failed"); } catch (const std::length_error&) {
        return fail(ExpressionDeepRenderStatus::ResourceLimit, "renderer allocation exceeds container limits");
    } catch (...) { return fail(ExpressionDeepRenderStatus::InternalError, "renderer failed with an unexpected exception"); }
}

} // namespace formula
