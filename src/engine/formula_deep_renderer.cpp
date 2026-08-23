#include "formula_deep_renderer.h"
#include "bigfixed.h"
#include "formula_expression_centered.h"
#include "formula_expression_vector_math.h"

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

class ExpressionDoubleDoubleEvaluator {
  public:
    struct Real {
        double high = 0.0;
        double low = 0.0;
    };

    struct ComplexValue {
        Real real;
        Real imaginary;
    };

    struct OrbitResult {
        bool defined = false;
        int escapeIteration = -1;
        uint64_t iterations = 0;
    };

    static bool supports(const ExpressionProgram& program) {
        if (!program._valid) return false;
        for (const ExpressionProgram::Instruction& instruction : program._code) {
            switch (instruction.op) {
            case ExpressionProgram::Op::Constant:
            case ExpressionProgram::Op::Z:
            case ExpressionProgram::Op::C:
            case ExpressionProgram::Op::Z0:
            case ExpressionProgram::Op::Iteration:
            case ExpressionProgram::Op::Parameter:
            case ExpressionProgram::Op::Negate:
            case ExpressionProgram::Op::Add:
            case ExpressionProgram::Op::Subtract:
            case ExpressionProgram::Op::Multiply:
            case ExpressionProgram::Op::Divide:
            case ExpressionProgram::Op::Square:
            case ExpressionProgram::Op::Sin:
            case ExpressionProgram::Op::Cos:
            case ExpressionProgram::Op::Tan:
            case ExpressionProgram::Op::Sinh:
            case ExpressionProgram::Op::Cosh:
            case ExpressionProgram::Op::Tanh:
            case ExpressionProgram::Op::Exp:
            case ExpressionProgram::Op::Log:
            case ExpressionProgram::Op::Log10:
            case ExpressionProgram::Op::Sqrt:
            case ExpressionProgram::Op::Power:
            case ExpressionProgram::Op::Norm:
            case ExpressionProgram::Op::Abs:
            case ExpressionProgram::Op::Conjugate:
            case ExpressionProgram::Op::Real:
            case ExpressionProgram::Op::Imaginary:
            case ExpressionProgram::Op::MakeComplex: break;
            default: return false;
            }
        }
        return true;
    }

    static bool splitMpfr(mpfr_srcptr value, mpfr_ptr temporary, Real& output) {
        if (!mpfr_number_p(value)) return false;
        output.high = mpfr_get_d(value, MPFR_RNDN);
        if (!std::isfinite(output.high)) return false;
        mpfr_set_d(temporary, output.high, MPFR_RNDN);
        mpfr_sub(temporary, value, temporary, MPFR_RNDN);
        output.low = mpfr_get_d(temporary, MPFR_RNDN);
        return std::isfinite(output.low) && (mpfr_zero_p(value) || output.high != 0.0 || output.low != 0.0);
    }

    template <typename Cancel>
    static OrbitResult evaluateOrbit(const ExpressionProgram& program, const ExpressionContext& fixed, FormulaParameter pixelParameter, ComplexValue pixel, int maxIterations, double bailout, Cancel&& shouldCancel) {
        OrbitResult result;
        if (!supports(program) || (pixelParameter != FormulaParameter::C && pixelParameter != FormulaParameter::InitialZ) || maxIterations < 1 || !std::isfinite(bailout) || bailout <= 0.0 || !finite(pixel)) return result;
        const ComplexValue pixelC = pixelParameter == FormulaParameter::C ? pixel : fromComplex(fixed.c);
        const ComplexValue pixelZ0 = pixelParameter == FormulaParameter::InitialZ ? pixel : fromComplex(fixed.z0);
        ComplexValue state = pixelZ0;
        Real bailoutSquared;
        if (!square(fromDouble(bailout), bailoutSquared)) return result;
        BailoutDecision decision = outsideBailout(state, bailoutSquared);
        if (decision == BailoutDecision::Invalid) return result;
        if (decision == BailoutDecision::Outside) {
            result.defined = true;
            result.escapeIteration = 0;
            return result;
        }
        for (int iteration = 0; iteration < maxIterations; ++iteration) {
            if ((iteration & 63) == 0 && shouldCancel()) return result;
            ComplexValue next;
            if (!evaluateStep(program, fixed, pixelC, pixelZ0, state, iteration, next)) return result;
            ++result.iterations;
            state = next;
            decision = outsideBailout(state, bailoutSquared);
            if (decision == BailoutDecision::Invalid) return result;
            if (decision == BailoutDecision::Outside) {
                result.defined = true;
                result.escapeIteration = iteration + 1;
                return result;
            }
        }
        result.defined = true;
        return result;
    }

    template <typename Cancel>
    static void evaluateOrbit4(const ExpressionProgram& program, const ExpressionContext& fixed, FormulaParameter pixelParameter, const ComplexValue* pixels, uint8_t inputMask, int maxIterations, double bailout, Cancel&& shouldCancel, OrbitResult* results) {
        for (int lane = 0; lane < 4; ++lane) results[lane] = {};
        if (!supports(program) || (pixelParameter != FormulaParameter::C && pixelParameter != FormulaParameter::InitialZ) || !pixels || inputMask == 0 || maxIterations < 1 || !std::isfinite(bailout) || bailout <= 0.0) return;
        Real bailoutSquared;
        if (!square(fromDouble(bailout), bailoutSquared)) return;

        uint8_t activeMask = inputMask;
        for (int lane = 0; lane < 4; ++lane) {
            const uint8_t laneMask = uint8_t{1} << lane;
            if ((activeMask & laneMask) == 0) continue;
            if (!finite(pixels[lane])) {
                activeMask &= static_cast<uint8_t>(~laneMask);
                continue;
            }
            const ComplexValue initialState = pixelParameter == FormulaParameter::InitialZ ? pixels[lane] : fromComplex(fixed.z0);
            const BailoutDecision initialDecision = outsideBailout(initialState, bailoutSquared);
            if (initialDecision == BailoutDecision::Invalid) {
                activeMask &= static_cast<uint8_t>(~laneMask);
            } else if (initialDecision == BailoutDecision::Outside) {
                results[lane].defined = true;
                results[lane].escapeIteration = 0;
                activeMask &= static_cast<uint8_t>(~laneMask);
            }
        }
        if (activeMask == 0) return;

        alignas(32) double cRealHigh[4]{};
        alignas(32) double cRealLow[4]{};
        alignas(32) double cImaginaryHigh[4]{};
        alignas(32) double cImaginaryLow[4]{};
        for (int lane = 0; lane < 4; ++lane) {
            if ((activeMask & (uint8_t{1} << lane)) == 0) continue;
            cRealHigh[lane] = pixels[lane].real.high;
            cRealLow[lane] = pixels[lane].real.low;
            cImaginaryHigh[lane] = pixels[lane].imaginary.high;
            cImaginaryLow[lane] = pixels[lane].imaginary.low;
        }
        VectorComplex pixelCoordinates{{_mm256_load_pd(cRealHigh), _mm256_load_pd(cRealLow)}, {_mm256_load_pd(cImaginaryHigh), _mm256_load_pd(cImaginaryLow)}};
        const VectorComplex pixelC = pixelParameter == FormulaParameter::C ? pixelCoordinates : vectorFromComplex(fixed.c);
        const VectorComplex pixelZ0 = pixelParameter == FormulaParameter::InitialZ ? pixelCoordinates : vectorFromComplex(fixed.z0);
        VectorComplex state = pixelZ0;
        const __m256d bailoutHigh = _mm256_set1_pd(bailoutSquared.high);
        const __m256d bailoutLow = _mm256_set1_pd(bailoutSquared.low);
        for (int iteration = 0; iteration < maxIterations && activeMask != 0; ++iteration) {
            if ((iteration & 63) == 0 && shouldCancel()) return;
            VectorComplex next;
            if (!evaluateStep4(program, fixed, pixelC, pixelZ0, state, iteration, next)) return;
            state = next;
            VectorReal realSquared;
            VectorReal imaginarySquared;
            VectorReal magnitudeSquared;
            if (!vectorSquare(state.real, realSquared) || !vectorSquare(state.imaginary, imaginarySquared) || !vectorAdd(realSquared, imaginarySquared, magnitudeSquared)) return;
            __m256d finiteState = _mm256_and_pd(_mm256_cmp_pd(state.real.high, state.real.high, _CMP_ORD_Q), _mm256_cmp_pd(state.real.low, state.real.low, _CMP_ORD_Q));
            finiteState = _mm256_and_pd(finiteState, _mm256_cmp_pd(state.imaginary.high, state.imaginary.high, _CMP_ORD_Q));
            finiteState = _mm256_and_pd(finiteState, _mm256_cmp_pd(state.imaginary.low, state.imaginary.low, _CMP_ORD_Q));
            finiteState = _mm256_and_pd(finiteState, _mm256_cmp_pd(magnitudeSquared.high, magnitudeSquared.high, _CMP_ORD_Q));
            finiteState = _mm256_and_pd(finiteState, _mm256_cmp_pd(magnitudeSquared.low, magnitudeSquared.low, _CMP_ORD_Q));
            const uint8_t finiteMask = static_cast<uint8_t>(_mm256_movemask_pd(finiteState));
            const __m256d highOutside = _mm256_cmp_pd(magnitudeSquared.high, bailoutHigh, _CMP_GT_OQ);
            const __m256d highEqual = _mm256_cmp_pd(magnitudeSquared.high, bailoutHigh, _CMP_EQ_OQ);
            const __m256d lowOutside = _mm256_cmp_pd(magnitudeSquared.low, bailoutLow, _CMP_GT_OQ);
            const uint8_t outsideMask = activeMask & finiteMask & static_cast<uint8_t>(_mm256_movemask_pd(_mm256_or_pd(highOutside, _mm256_and_pd(highEqual, lowOutside))));
            for (int lane = 0; lane < 4; ++lane) {
                const uint8_t laneMask = uint8_t{1} << lane;
                if ((activeMask & laneMask) == 0) continue;
                ++results[lane].iterations;
                if ((finiteMask & laneMask) == 0) {
                    activeMask &= static_cast<uint8_t>(~laneMask);
                    continue;
                }
                if ((outsideMask & laneMask) != 0) {
                    results[lane].defined = true;
                    results[lane].escapeIteration = iteration + 1;
                    activeMask &= static_cast<uint8_t>(~laneMask);
                } else if (iteration + 1 == maxIterations) {
                    results[lane].defined = true;
                    activeMask &= static_cast<uint8_t>(~laneMask);
                }
            }
        }
    }

  private:
    enum class BailoutDecision : uint8_t {
        Inside,
        Outside,
        Invalid
    };

    struct VectorReal {
        __m256d high;
        __m256d low;
    };

    struct VectorComplex {
        VectorReal real;
        VectorReal imaginary;
    };

    static VectorReal vectorFromDouble(double value) {
        return {_mm256_set1_pd(value), _mm256_setzero_pd()};
    }

    static VectorReal vectorFromReal(Real value) {
        return {_mm256_set1_pd(value.high), _mm256_set1_pd(value.low)};
    }

    static VectorComplex vectorFromComplex(Complex value) {
        return {vectorFromDouble(value.real()), vectorFromDouble(value.imag())};
    }

    static VectorReal vectorNormalize(__m256d high, __m256d low) {
        const __m256d sum = _mm256_add_pd(high, low);
        return {sum, _mm256_sub_pd(low, _mm256_sub_pd(sum, high))};
    }

    static VectorReal vectorTwoSum(__m256d left, __m256d right) {
        const __m256d sum = _mm256_add_pd(left, right);
        const __m256d virtualRight = _mm256_sub_pd(sum, left);
        return {sum, _mm256_add_pd(_mm256_sub_pd(left, _mm256_sub_pd(sum, virtualRight)), _mm256_sub_pd(right, virtualRight))};
    }

    static bool vectorAdd(VectorReal left, VectorReal right, VectorReal& output) {
        VectorReal high = vectorTwoSum(left.high, right.high);
        const VectorReal low = vectorTwoSum(left.low, right.low);
        high.low = _mm256_add_pd(high.low, low.high);
        high = vectorNormalize(high.high, high.low);
        high.low = _mm256_add_pd(high.low, low.low);
        output = vectorNormalize(high.high, high.low);
        return true;
    }

    static bool vectorSubtract(VectorReal left, VectorReal right, VectorReal& output) {
        const __m256d signMask = _mm256_set1_pd(-0.0);
        return vectorAdd(left, {_mm256_xor_pd(right.high, signMask), _mm256_xor_pd(right.low, signMask)}, output);
    }

    static bool vectorMultiply(VectorReal left, VectorReal right, VectorReal& output) {
        const __m256d product = _mm256_mul_pd(left.high, right.high);
        __m256d error = _mm256_fmadd_pd(left.high, right.high, _mm256_xor_pd(product, _mm256_set1_pd(-0.0)));
        error = _mm256_add_pd(error, _mm256_add_pd(_mm256_mul_pd(left.high, right.low), _mm256_mul_pd(left.low, right.high)));
        output = vectorNormalize(product, error);
        output.low = _mm256_add_pd(output.low, _mm256_mul_pd(left.low, right.low));
        output = vectorNormalize(output.high, output.low);
        return true;
    }

    static bool vectorSquare(VectorReal value, VectorReal& output) {
        return vectorMultiply(value, value, output);
    }

    static void vectorRealSinCosDouble(__m256d value, __m256d& sine, __m256d& cosine) {
        int quadrants[4];
        detail::vectorMathRealSinCos(value, sine, cosine, quadrants);
        alignas(32) double sineValues[4];
        alignas(32) double cosineValues[4];
        _mm256_store_pd(sineValues, sine);
        _mm256_store_pd(cosineValues, cosine);
        for (int lane = 0; lane < 4; ++lane) {
            switch (quadrants[lane] & 3) {
            case 1:
                std::swap(sineValues[lane], cosineValues[lane]);
                cosineValues[lane] = -cosineValues[lane];
                break;
            case 2:
                sineValues[lane] = -sineValues[lane];
                cosineValues[lane] = -cosineValues[lane];
                break;
            case 3:
                std::swap(sineValues[lane], cosineValues[lane]);
                sineValues[lane] = -sineValues[lane];
                break;
            }
        }
        sine = _mm256_load_pd(sineValues);
        cosine = _mm256_load_pd(cosineValues);
    }

    static bool vectorSineCosine(VectorReal value, VectorReal& sineOutput, VectorReal& cosineOutput) {
        constexpr Real PiOverTwo{0x1.921fb54442d18p+0, 0x1.1a62633145c07p-54};
        constexpr double TwoOverPi = 0x1.45f306dc9c883p-1;
        constexpr Real SineCoefficients[] = {
            {-0x1.5555555555555p-3, -0x1.5555555555555p-57},
            {0x1.1111111111111p-7, 0x1.1111111111111p-63},
            {-0x1.a01a01a01a01ap-13, -0x1.a01a01a01a01ap-73},
            {0x1.71de3a556c734p-19, -0x1.c154f8ddc6c00p-73},
            {-0x1.ae64567f544e4p-26, 0x1.c062e06d1f209p-80},
            {0x1.6124613a86d09p-33, 0x1.f28e0cc748ebep-87},
            {-0x1.ae7f3e733b81fp-41, -0x1.1d8656b0ee8cbp-97},
            {0x1.952c77030ad4ap-49, 0x1.ac981465ddc6cp-103},
            {-0x1.2f49b46814157p-57, -0x1.2650f61dbdcb4p-112},
            {0x1.71b8ef6dcf572p-66, -0x1.d043ae40c4647p-120},
            {-0x1.761b41316381ap-75, 0x1.3423c7d91404fp-130},
            {0x1.3f3ccdd165fa9p-84, -0x1.58ddadf344487p-139},
            {-0x1.d1ab1c2dccea3p-94, -0x1.054d0c78aea14p-149},
            {0x1.259f98b4358adp-103, 0x1.eaf8c39dd9bc5p-157},
        };
        constexpr Real CosineCoefficients[] = {
            {-0x1.0000000000000p-1, 0x0.0p+0},
            {0x1.5555555555555p-5, 0x1.5555555555555p-59},
            {-0x1.6c16c16c16c17p-10, 0x1.f49f49f49f49fp-65},
            {0x1.a01a01a01a01ap-16, 0x1.a01a01a01a01ap-76},
            {-0x1.27e4fb7789f5cp-22, -0x1.cbbc05b4fa99ap-76},
            {0x1.1eed8eff8d898p-29, -0x1.2aec959e14c06p-83},
            {-0x1.93974a8c07c9dp-37, -0x1.05d6f8a2efd1fp-92},
            {0x1.ae7f3e733b81fp-45, 0x1.1d8656b0ee8cbp-101},
            {-0x1.6827863b97d97p-53, -0x1.eec01221a8b0bp-107},
            {0x1.e542ba4020225p-62, 0x1.ea72b4afe3c2fp-120},
            {-0x1.0ce396db7f853p-70, 0x1.aebcdbd20331cp-124},
            {0x1.f2cf01972f578p-80, -0x1.9ada5fcc1ab14p-135},
            {-0x1.88e85fc6a4e5ap-89, 0x1.71c37ebd16540p-143},
            {0x1.0a18a2635085dp-98, 0x1.b9e2e28e1aa54p-153},
            {-0x1.3932c5047d60ep-108, -0x1.832b7b530a627p-162},
        };

        const __m256d signMask = _mm256_set1_pd(-0.0);
        __m256d valid = _mm256_and_pd(_mm256_cmp_pd(value.high, value.high, _CMP_ORD_Q), _mm256_cmp_pd(value.low, value.low, _CMP_ORD_Q));
        valid = _mm256_and_pd(valid, _mm256_cmp_pd(_mm256_andnot_pd(signMask, value.high), _mm256_set1_pd(0x1p20), _CMP_LE_OQ));
        const VectorReal safeValue{_mm256_and_pd(value.high, valid), _mm256_and_pd(value.low, valid)};
        const __m256d multiple = _mm256_round_pd(_mm256_add_pd(_mm256_mul_pd(safeValue.high, _mm256_set1_pd(TwoOverPi)), _mm256_mul_pd(safeValue.low, _mm256_set1_pd(TwoOverPi))), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
        VectorReal quadrantOffset;
        VectorReal reduced;
        if (!vectorMultiply(vectorFromReal(PiOverTwo), {multiple, _mm256_setzero_pd()}, quadrantOffset) || !vectorSubtract(safeValue, quadrantOffset, reduced)) return false;
        VectorReal squared;
        if (!vectorSquare(reduced, squared)) return false;

        VectorReal sinePolynomial = vectorFromReal(SineCoefficients[std::size(SineCoefficients) - 1]);
        VectorReal cosinePolynomial = vectorFromReal(CosineCoefficients[std::size(CosineCoefficients) - 1]);
        if (!vectorMultiply(squared, cosinePolynomial, cosinePolynomial) || !vectorAdd(vectorFromReal(CosineCoefficients[std::size(CosineCoefficients) - 2]), cosinePolynomial, cosinePolynomial)) return false;
        for (size_t index = std::size(SineCoefficients) - 1; index > 0; --index) {
            VectorReal sineProduct;
            VectorReal cosineProduct;
            if (!vectorMultiply(squared, sinePolynomial, sineProduct) || !vectorMultiply(squared, cosinePolynomial, cosineProduct) || !vectorAdd(vectorFromReal(SineCoefficients[index - 1]), sineProduct, sinePolynomial) || !vectorAdd(vectorFromReal(CosineCoefficients[index - 1]), cosineProduct, cosinePolynomial)) return false;
        }
        VectorReal cubic;
        VectorReal sineCorrection;
        VectorReal cosineCorrection;
        if (!vectorMultiply(reduced, squared, cubic) || !vectorMultiply(cubic, sinePolynomial, sineCorrection) || !vectorAdd(reduced, sineCorrection, sineOutput) || !vectorMultiply(squared, cosinePolynomial, cosineCorrection) || !vectorAdd(vectorFromDouble(1.0), cosineCorrection, cosineOutput)) return false;

        alignas(16) int quadrants[4];
        _mm_store_si128(reinterpret_cast<__m128i*>(quadrants), _mm256_cvttpd_epi32(multiple));
        alignas(32) double sineHigh[4];
        alignas(32) double sineLow[4];
        alignas(32) double cosineHigh[4];
        alignas(32) double cosineLow[4];
        _mm256_store_pd(sineHigh, sineOutput.high);
        _mm256_store_pd(sineLow, sineOutput.low);
        _mm256_store_pd(cosineHigh, cosineOutput.high);
        _mm256_store_pd(cosineLow, cosineOutput.low);
        for (int lane = 0; lane < 4; ++lane) {
            switch (quadrants[lane] & 3) {
            case 1:
                std::swap(sineHigh[lane], cosineHigh[lane]);
                std::swap(sineLow[lane], cosineLow[lane]);
                cosineHigh[lane] = -cosineHigh[lane];
                cosineLow[lane] = -cosineLow[lane];
                break;
            case 2:
                sineHigh[lane] = -sineHigh[lane];
                sineLow[lane] = -sineLow[lane];
                cosineHigh[lane] = -cosineHigh[lane];
                cosineLow[lane] = -cosineLow[lane];
                break;
            case 3:
                std::swap(sineHigh[lane], cosineHigh[lane]);
                std::swap(sineLow[lane], cosineLow[lane]);
                sineHigh[lane] = -sineHigh[lane];
                sineLow[lane] = -sineLow[lane];
                break;
            }
        }
        sineOutput = {_mm256_load_pd(sineHigh), _mm256_load_pd(sineLow)};
        cosineOutput = {_mm256_load_pd(cosineHigh), _mm256_load_pd(cosineLow)};
        const __m256d nan = _mm256_set1_pd(std::numeric_limits<double>::quiet_NaN());
        sineOutput.high = _mm256_blendv_pd(nan, sineOutput.high, valid);
        sineOutput.low = _mm256_blendv_pd(nan, sineOutput.low, valid);
        cosineOutput.high = _mm256_blendv_pd(nan, cosineOutput.high, valid);
        cosineOutput.low = _mm256_blendv_pd(nan, cosineOutput.low, valid);
        return true;
    }

    static bool vectorSine(VectorReal value, VectorReal& output) {
        VectorReal cosineOutput;
        return vectorSineCosine(value, output, cosineOutput);
    }

    static bool vectorCosine(VectorReal value, VectorReal& output) {
        VectorReal sineOutput;
        return vectorSineCosine(value, sineOutput, output);
    }

    static bool vectorHyperbolicSineCosine(VectorReal value, VectorReal& sineOutput, VectorReal& cosineOutput) {
        constexpr Real LnTwo{0x1.62e42fefa39efp-1, 0x1.abc9e3b39803fp-56};
        constexpr double LogTwoE = 0x1.71547652b82fep+0;
        constexpr Real SineCoefficients[] = {
            {0x1.0000000000000p+0, 0x0.0p+0},
            {0x1.5555555555555p-3, 0x1.5555555555555p-57},
            {0x1.1111111111111p-7, 0x1.1111111111111p-63},
            {0x1.a01a01a01a01ap-13, 0x1.a01a01a01a01ap-73},
            {0x1.71de3a556c734p-19, -0x1.c154f8ddc6c00p-73},
            {0x1.ae64567f544e4p-26, -0x1.c062e06d1f209p-80},
            {0x1.6124613a86d09p-33, 0x1.f28e0cc748ebep-87},
            {0x1.ae7f3e733b81fp-41, 0x1.1d8656b0ee8cbp-97},
            {0x1.952c77030ad4ap-49, 0x1.ac981465ddc6cp-103},
            {0x1.2f49b46814157p-57, 0x1.2650f61dbdcb4p-112},
            {0x1.71b8ef6dcf572p-66, -0x1.d043ae40c4647p-120},
        };
        constexpr Real CosineCoefficients[] = {
            {0x1.0000000000000p-1, 0x0.0p+0},
            {0x1.5555555555555p-5, 0x1.5555555555555p-59},
            {0x1.6c16c16c16c17p-10, -0x1.f49f49f49f49fp-65},
            {0x1.a01a01a01a01ap-16, 0x1.a01a01a01a01ap-76},
            {0x1.27e4fb7789f5cp-22, 0x1.cbbc05b4fa99ap-76},
            {0x1.1eed8eff8d898p-29, -0x1.2aec959e14c06p-83},
            {0x1.93974a8c07c9dp-37, 0x1.05d6f8a2efd1fp-92},
            {0x1.ae7f3e733b81fp-45, 0x1.1d8656b0ee8cbp-101},
            {0x1.6827863b97d97p-53, 0x1.eec01221a8b0bp-107},
            {0x1.e542ba4020225p-62, 0x1.ea72b4afe3c2fp-120},
            {0x1.0ce396db7f853p-70, -0x1.aebcdbd20331cp-124},
        };

        const __m256d signMask = _mm256_set1_pd(-0.0);
        __m256d valid = _mm256_and_pd(_mm256_cmp_pd(value.high, value.high, _CMP_ORD_Q), _mm256_cmp_pd(value.low, value.low, _CMP_ORD_Q));
        valid = _mm256_and_pd(valid, _mm256_cmp_pd(_mm256_andnot_pd(signMask, value.high), _mm256_set1_pd(709.0), _CMP_LE_OQ));
        const VectorReal safeValue{_mm256_and_pd(value.high, valid), _mm256_and_pd(value.low, valid)};
        const __m256d multiple = _mm256_round_pd(_mm256_add_pd(_mm256_mul_pd(safeValue.high, _mm256_set1_pd(LogTwoE)), _mm256_mul_pd(safeValue.low, _mm256_set1_pd(LogTwoE))), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
        VectorReal offset;
        VectorReal reduced;
        if (!vectorMultiply(vectorFromReal(LnTwo), {multiple, _mm256_setzero_pd()}, offset) || !vectorSubtract(safeValue, offset, reduced)) return false;
        VectorReal squared;
        if (!vectorSquare(reduced, squared)) return false;

        VectorReal sinePolynomial = vectorFromReal(SineCoefficients[std::size(SineCoefficients) - 1]);
        VectorReal cosinePolynomial = vectorFromReal(CosineCoefficients[std::size(CosineCoefficients) - 1]);
        for (size_t index = std::size(SineCoefficients) - 1; index > 0; --index) {
            VectorReal sineProduct;
            VectorReal cosineProduct;
            if (!vectorMultiply(squared, sinePolynomial, sineProduct) || !vectorMultiply(squared, cosinePolynomial, cosineProduct) || !vectorAdd(vectorFromReal(SineCoefficients[index - 1]), sineProduct, sinePolynomial) || !vectorAdd(vectorFromReal(CosineCoefficients[index - 1]), cosineProduct, cosinePolynomial)) return false;
        }
        VectorReal reducedSine;
        if (!vectorMultiply(reduced, sinePolynomial, reducedSine)) return false;
        VectorReal cosineCorrection;
        VectorReal reducedCosine;
        if (!vectorMultiply(squared, cosinePolynomial, cosineCorrection) || !vectorAdd(vectorFromDouble(1.0), cosineCorrection, reducedCosine)) return false;

        alignas(16) int exponents[4];
        alignas(32) double sumScaleHigh[4];
        alignas(32) double sumScaleLow[4];
        alignas(32) double differenceScaleHigh[4];
        alignas(32) double differenceScaleLow[4];
        _mm_store_si128(reinterpret_cast<__m128i*>(exponents), _mm256_cvttpd_epi32(multiple));
        for (int lane = 0; lane < 4; ++lane) {
            const double positive = std::ldexp(0.5, exponents[lane]);
            const double negative = std::ldexp(0.5, -exponents[lane]);
            const double sum = positive + negative;
            const double sumVirtualNegative = sum - positive;
            sumScaleHigh[lane] = sum;
            sumScaleLow[lane] = (positive - (sum - sumVirtualNegative)) + (negative - sumVirtualNegative);
            const double difference = positive - negative;
            const double differenceVirtualNegative = difference - positive;
            differenceScaleHigh[lane] = difference;
            differenceScaleLow[lane] = (positive - (difference - differenceVirtualNegative)) - (negative + differenceVirtualNegative);
        }
        const VectorReal sum{_mm256_load_pd(sumScaleHigh), _mm256_load_pd(sumScaleLow)};
        const VectorReal difference{_mm256_load_pd(differenceScaleHigh), _mm256_load_pd(differenceScaleLow)};
        VectorReal first;
        VectorReal second;
        if (!vectorMultiply(reducedCosine, difference, first) || !vectorMultiply(reducedSine, sum, second) || !vectorAdd(first, second, sineOutput) || !vectorMultiply(reducedCosine, sum, first) || !vectorMultiply(reducedSine, difference, second) || !vectorAdd(first, second, cosineOutput)) return false;
        const __m256d nan = _mm256_set1_pd(std::numeric_limits<double>::quiet_NaN());
        sineOutput.high = _mm256_blendv_pd(nan, sineOutput.high, valid);
        sineOutput.low = _mm256_blendv_pd(nan, sineOutput.low, valid);
        cosineOutput.high = _mm256_blendv_pd(nan, cosineOutput.high, valid);
        cosineOutput.low = _mm256_blendv_pd(nan, cosineOutput.low, valid);
        return true;
    }

    static bool vectorHyperbolicSine(VectorReal value, VectorReal& output) {
        VectorReal cosineOutput;
        return vectorHyperbolicSineCosine(value, output, cosineOutput);
    }

    static bool vectorHyperbolicCosine(VectorReal value, VectorReal& output) {
        VectorReal sineOutput;
        return vectorHyperbolicSineCosine(value, sineOutput, output);
    }

    static bool vectorExponential(VectorReal value, VectorReal& output) {
        constexpr Real LnTwo{0x1.62e42fefa39efp-1, 0x1.abc9e3b39803fp-56};
        constexpr double LogTwoE = 0x1.71547652b82fep+0;
        constexpr Real Coefficients[] = {
            {0x1.0000000000000p+0, 0x0.0p+0},
            {0x1.0000000000000p-1, 0x0.0p+0},
            {0x1.5555555555555p-3, 0x1.5555555555555p-57},
            {0x1.5555555555555p-5, 0x1.5555555555555p-59},
            {0x1.1111111111111p-7, 0x1.1111111111111p-63},
            {0x1.6c16c16c16c17p-10, -0x1.f49f49f49f49fp-65},
            {0x1.a01a01a01a01ap-13, 0x1.a01a01a01a01ap-73},
            {0x1.a01a01a01a01ap-16, 0x1.a01a01a01a01ap-76},
            {0x1.71de3a556c734p-19, -0x1.c154f8ddc6c00p-73},
            {0x1.27e4fb7789f5cp-22, 0x1.cbbc05b4fa99ap-76},
            {0x1.ae64567f544e4p-26, -0x1.c062e06d1f209p-80},
            {0x1.1eed8eff8d898p-29, -0x1.2aec959e14c06p-83},
            {0x1.6124613a86d09p-33, 0x1.f28e0cc748ebep-87},
            {0x1.93974a8c07c9dp-37, 0x1.05d6f8a2efd1fp-92},
            {0x1.ae7f3e733b81fp-41, 0x1.1d8656b0ee8cbp-97},
            {0x1.ae7f3e733b81fp-45, 0x1.1d8656b0ee8cbp-101},
            {0x1.952c77030ad4ap-49, 0x1.ac981465ddc6cp-103},
            {0x1.6827863b97d97p-53, 0x1.eec01221a8b0bp-107},
            {0x1.2f49b46814157p-57, 0x1.2650f61dbdcb4p-112},
            {0x1.e542ba4020225p-62, 0x1.ea72b4afe3c2fp-120},
            {0x1.71b8ef6dcf572p-66, -0x1.d043ae40c4647p-120},
            {0x1.0ce396db7f853p-70, -0x1.aebcdbd20331cp-124},
        };
        __m256d valid = _mm256_and_pd(_mm256_cmp_pd(value.high, value.high, _CMP_ORD_Q), _mm256_cmp_pd(value.low, value.low, _CMP_ORD_Q));
        valid = _mm256_and_pd(valid, _mm256_cmp_pd(value.high, _mm256_set1_pd(-745.0), _CMP_GE_OQ));
        valid = _mm256_and_pd(valid, _mm256_cmp_pd(value.high, _mm256_set1_pd(709.0), _CMP_LE_OQ));
        const VectorReal safeValue{_mm256_and_pd(value.high, valid), _mm256_and_pd(value.low, valid)};
        const __m256d multiple = _mm256_round_pd(_mm256_add_pd(_mm256_mul_pd(safeValue.high, _mm256_set1_pd(LogTwoE)), _mm256_mul_pd(safeValue.low, _mm256_set1_pd(LogTwoE))), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
        VectorReal offset;
        VectorReal reduced;
        if (!vectorMultiply(vectorFromReal(LnTwo), {multiple, _mm256_setzero_pd()}, offset) || !vectorSubtract(safeValue, offset, reduced)) return false;
        VectorReal polynomial = vectorFromReal(Coefficients[std::size(Coefficients) - 1]);
        for (size_t index = std::size(Coefficients) - 1; index-- > 0;) {
            if (!vectorMultiply(reduced, polynomial, polynomial) || !vectorAdd(vectorFromReal(Coefficients[index]), polynomial, polynomial)) return false;
        }
        if (!vectorMultiply(reduced, polynomial, output) || !vectorAdd(vectorFromDouble(1.0), output, output)) return false;
        alignas(16) int exponents[4];
        alignas(32) double high[4];
        alignas(32) double low[4];
        _mm_store_si128(reinterpret_cast<__m128i*>(exponents), _mm256_cvttpd_epi32(multiple));
        _mm256_store_pd(high, output.high);
        _mm256_store_pd(low, output.low);
        for (int lane = 0; lane < 4; ++lane) {
            high[lane] = std::ldexp(high[lane], exponents[lane]);
            low[lane] = std::ldexp(low[lane], exponents[lane]);
        }
        output = {_mm256_load_pd(high), _mm256_load_pd(low)};
        const __m256d nan = _mm256_set1_pd(std::numeric_limits<double>::quiet_NaN());
        output.high = _mm256_blendv_pd(nan, output.high, valid);
        output.low = _mm256_blendv_pd(nan, output.low, valid);
        return true;
    }

    static bool vectorAdd(VectorComplex left, VectorComplex right, VectorComplex& output) {
        return vectorAdd(left.real, right.real, output.real) && vectorAdd(left.imaginary, right.imaginary, output.imaginary);
    }

    static bool vectorSubtract(VectorComplex left, VectorComplex right, VectorComplex& output) {
        return vectorSubtract(left.real, right.real, output.real) && vectorSubtract(left.imaginary, right.imaginary, output.imaginary);
    }

    static bool vectorMultiply(VectorComplex left, VectorComplex right, VectorComplex& output) {
        VectorReal realLeft;
        VectorReal realRight;
        VectorReal imaginaryLeft;
        VectorReal imaginaryRight;
        return vectorMultiply(left.real, right.real, realLeft) && vectorMultiply(left.imaginary, right.imaginary, realRight) && vectorSubtract(realLeft, realRight, output.real) && vectorMultiply(left.real, right.imaginary, imaginaryLeft) && vectorMultiply(left.imaginary, right.real, imaginaryRight) && vectorAdd(imaginaryLeft, imaginaryRight, output.imaginary);
    }

    static bool vectorSquare(VectorComplex value, VectorComplex& output) {
        VectorReal realSquared;
        VectorReal imaginarySquared;
        VectorReal product;
        return vectorSquare(value.real, realSquared) && vectorSquare(value.imaginary, imaginarySquared) && vectorSubtract(realSquared, imaginarySquared, output.real) && vectorMultiply(value.real, value.imaginary, product) && vectorAdd(product, product, output.imaginary);
    }

    static bool vectorDivide(VectorComplex numerator, VectorComplex denominator, VectorComplex& output) {
        alignas(32) double numeratorRealHigh[4], numeratorRealLow[4], numeratorImaginaryHigh[4], numeratorImaginaryLow[4];
        alignas(32) double denominatorRealHigh[4], denominatorRealLow[4], denominatorImaginaryHigh[4], denominatorImaginaryLow[4];
        alignas(32) double outputRealHigh[4], outputRealLow[4], outputImaginaryHigh[4], outputImaginaryLow[4];
        _mm256_store_pd(numeratorRealHigh, numerator.real.high);
        _mm256_store_pd(numeratorRealLow, numerator.real.low);
        _mm256_store_pd(numeratorImaginaryHigh, numerator.imaginary.high);
        _mm256_store_pd(numeratorImaginaryLow, numerator.imaginary.low);
        _mm256_store_pd(denominatorRealHigh, denominator.real.high);
        _mm256_store_pd(denominatorRealLow, denominator.real.low);
        _mm256_store_pd(denominatorImaginaryHigh, denominator.imaginary.high);
        _mm256_store_pd(denominatorImaginaryLow, denominator.imaginary.low);
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (int lane = 0; lane < 4; ++lane) {
            ComplexValue quotient;
            const ComplexValue laneNumerator{{numeratorRealHigh[lane], numeratorRealLow[lane]}, {numeratorImaginaryHigh[lane], numeratorImaginaryLow[lane]}};
            const ComplexValue laneDenominator{{denominatorRealHigh[lane], denominatorRealLow[lane]}, {denominatorImaginaryHigh[lane], denominatorImaginaryLow[lane]}};
            if (!divide(laneNumerator, laneDenominator, quotient)) quotient = {{nan, nan}, {nan, nan}};
            outputRealHigh[lane] = quotient.real.high;
            outputRealLow[lane] = quotient.real.low;
            outputImaginaryHigh[lane] = quotient.imaginary.high;
            outputImaginaryLow[lane] = quotient.imaginary.low;
        }
        output = {{_mm256_load_pd(outputRealHigh), _mm256_load_pd(outputRealLow)}, {_mm256_load_pd(outputImaginaryHigh), _mm256_load_pd(outputImaginaryLow)}};
        return true;
    }

    static bool vectorSine(VectorComplex value, VectorComplex& output) {
        VectorReal sineReal;
        VectorReal cosineReal;
        VectorReal sineImaginary;
        VectorReal cosineImaginary;
        return vectorSineCosine(value.real, sineReal, cosineReal) && vectorHyperbolicSineCosine(value.imaginary, sineImaginary, cosineImaginary) && vectorMultiply(sineReal, cosineImaginary, output.real) && vectorMultiply(cosineReal, sineImaginary, output.imaginary);
    }

    static bool vectorCosine(VectorComplex value, VectorComplex& output) {
        VectorReal sineReal;
        VectorReal cosineReal;
        VectorReal sineImaginary;
        VectorReal cosineImaginary;
        if (!vectorSineCosine(value.real, sineReal, cosineReal) || !vectorHyperbolicSineCosine(value.imaginary, sineImaginary, cosineImaginary) || !vectorMultiply(cosineReal, cosineImaginary, output.real) || !vectorMultiply(sineReal, sineImaginary, output.imaginary)) return false;
        const __m256d signMask = _mm256_set1_pd(-0.0);
        output.imaginary.high = _mm256_xor_pd(output.imaginary.high, signMask);
        output.imaginary.low = _mm256_xor_pd(output.imaginary.low, signMask);
        return true;
    }

    static bool vectorTangent(VectorComplex value, VectorComplex& output) {
        VectorComplex sineValue;
        VectorComplex cosineValue;
        return vectorSine(value, sineValue) && vectorCosine(value, cosineValue) && vectorDivide(sineValue, cosineValue, output);
    }

    static bool vectorHyperbolicSine(VectorComplex value, VectorComplex& output) {
        VectorReal sineReal;
        VectorReal cosineReal;
        VectorReal sineImaginary;
        VectorReal cosineImaginary;
        return vectorHyperbolicSineCosine(value.real, sineReal, cosineReal) && vectorSineCosine(value.imaginary, sineImaginary, cosineImaginary) && vectorMultiply(sineReal, cosineImaginary, output.real) && vectorMultiply(cosineReal, sineImaginary, output.imaginary);
    }

    static bool vectorHyperbolicCosine(VectorComplex value, VectorComplex& output) {
        VectorReal sineReal;
        VectorReal cosineReal;
        VectorReal sineImaginary;
        VectorReal cosineImaginary;
        return vectorHyperbolicSineCosine(value.real, sineReal, cosineReal) && vectorSineCosine(value.imaginary, sineImaginary, cosineImaginary) && vectorMultiply(cosineReal, cosineImaginary, output.real) && vectorMultiply(sineReal, sineImaginary, output.imaginary);
    }

    static bool vectorHyperbolicTangent(VectorComplex value, VectorComplex& output) {
        VectorComplex sineValue;
        VectorComplex cosineValue;
        return vectorHyperbolicSine(value, sineValue) && vectorHyperbolicCosine(value, cosineValue) && vectorDivide(sineValue, cosineValue, output);
    }

    static bool vectorExponential(VectorComplex value, VectorComplex& output) {
        VectorReal magnitude;
        VectorReal sineImaginary;
        VectorReal cosineImaginary;
        return vectorExponential(value.real, magnitude) && vectorSineCosine(value.imaginary, sineImaginary, cosineImaginary) && vectorMultiply(magnitude, cosineImaginary, output.real) && vectorMultiply(magnitude, sineImaginary, output.imaginary);
    }

    static bool evaluateStep4(const ExpressionProgram& program, const ExpressionContext& fixed, VectorComplex pixelC, VectorComplex pixelZ0, VectorComplex state, int iteration, VectorComplex& output) {
        std::array<VectorComplex, ExpressionProgram::MAX_STACK> stack;
        size_t top = 0;
        auto push = [&](VectorComplex value) {
            if (top >= stack.size()) return false;
            stack[top++] = value;
            return true;
        };
        for (const ExpressionProgram::Instruction& instruction : program._code) {
            VectorComplex value;
            switch (instruction.op) {
            case ExpressionProgram::Op::Constant:
                if (!push(vectorFromComplex(instruction.value))) return false;
                break;
            case ExpressionProgram::Op::Z:
                if (!push(state)) return false;
                break;
            case ExpressionProgram::Op::C:
                if (!push(pixelC)) return false;
                break;
            case ExpressionProgram::Op::Z0:
                if (!push(pixelZ0)) return false;
                break;
            case ExpressionProgram::Op::Iteration:
                if (!push({vectorFromDouble(static_cast<double>(iteration)), vectorFromDouble(0.0)})) return false;
                break;
            case ExpressionProgram::Op::Parameter:
                if (instruction.argument >= fixed.parameters.size() || !push(vectorFromComplex(fixed.parameters[instruction.argument]))) return false;
                break;
            case ExpressionProgram::Op::Negate: {
                if (top < 1) return false;
                const __m256d signMask = _mm256_set1_pd(-0.0);
                stack[top - 1].real.high = _mm256_xor_pd(stack[top - 1].real.high, signMask);
                stack[top - 1].real.low = _mm256_xor_pd(stack[top - 1].real.low, signMask);
                stack[top - 1].imaginary.high = _mm256_xor_pd(stack[top - 1].imaginary.high, signMask);
                stack[top - 1].imaginary.low = _mm256_xor_pd(stack[top - 1].imaginary.low, signMask);
                break;
            }
            case ExpressionProgram::Op::Add:
            case ExpressionProgram::Op::Subtract:
            case ExpressionProgram::Op::Multiply:
            case ExpressionProgram::Op::Divide:
            case ExpressionProgram::Op::Power:
            case ExpressionProgram::Op::MakeComplex:
                if (top < 2) return false;
                if (instruction.op == ExpressionProgram::Op::Add) {
                    if (!vectorAdd(stack[top - 2], stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Subtract) {
                    if (!vectorSubtract(stack[top - 2], stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::MakeComplex) {
                    value = {stack[top - 2].real, stack[top - 1].real};
                } else if (instruction.op == ExpressionProgram::Op::Divide) {
                    if (!vectorDivide(stack[top - 2], stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Power) {
                    alignas(32) double lReH[4], lImH[4], rReH[4], rImH[4], oReH[4], oImH[4];
                    _mm256_store_pd(lReH, stack[top - 2].real.high);
                    _mm256_store_pd(lImH, stack[top - 2].imaginary.high);
                    _mm256_store_pd(rReH, stack[top - 1].real.high);
                    _mm256_store_pd(rImH, stack[top - 1].imaginary.high);
                    for (int lane = 0; lane < 4; ++lane) {
                        Complex res = std::pow(Complex{lReH[lane], lImH[lane]}, Complex{rReH[lane], rImH[lane]});
                        oReH[lane] = res.real();
                        oImH[lane] = res.imag();
                    }
                    value = {{_mm256_load_pd(oReH), _mm256_setzero_pd()}, {_mm256_load_pd(oImH), _mm256_setzero_pd()}};
                } else if (!vectorMultiply(stack[top - 2], stack[top - 1], value)) {
                    return false;
                }
                stack[--top - 1] = value;
                break;
            case ExpressionProgram::Op::Square:
            case ExpressionProgram::Op::Sin:
            case ExpressionProgram::Op::Cos:
            case ExpressionProgram::Op::Tan:
            case ExpressionProgram::Op::Sinh:
            case ExpressionProgram::Op::Cosh:
            case ExpressionProgram::Op::Tanh:
            case ExpressionProgram::Op::Exp:
            case ExpressionProgram::Op::Log:
            case ExpressionProgram::Op::Log10:
            case ExpressionProgram::Op::Sqrt:
            case ExpressionProgram::Op::Norm:
            case ExpressionProgram::Op::Abs:
            case ExpressionProgram::Op::Conjugate:
            case ExpressionProgram::Op::Real:
            case ExpressionProgram::Op::Imaginary:
                if (top < 1) return false;
                if (instruction.op == ExpressionProgram::Op::Square) {
                    if (!vectorSquare(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Sin) {
                    if (!vectorSine(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Cos) {
                    if (!vectorCosine(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Tan) {
                    if (!vectorTangent(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Sinh) {
                    if (!vectorHyperbolicSine(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Cosh) {
                    if (!vectorHyperbolicCosine(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Tanh) {
                    if (!vectorHyperbolicTangent(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Exp) {
                    if (!vectorExponential(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Log || instruction.op == ExpressionProgram::Op::Log10 || instruction.op == ExpressionProgram::Op::Sqrt) {
                    alignas(32) double reH[4], imH[4], oReH[4], oImH[4];
                    _mm256_store_pd(reH, stack[top - 1].real.high);
                    _mm256_store_pd(imH, stack[top - 1].imaginary.high);
                    for (int lane = 0; lane < 4; ++lane) {
                        Complex c{reH[lane], imH[lane]};
                        Complex res = instruction.op == ExpressionProgram::Op::Sqrt ? std::sqrt(c) : (instruction.op == ExpressionProgram::Op::Log10 ? std::log(c) / std::log(10.0) : std::log(c));
                        oReH[lane] = res.real();
                        oImH[lane] = res.imag();
                    }
                    value = {{_mm256_load_pd(oReH), _mm256_setzero_pd()}, {_mm256_load_pd(oImH), _mm256_setzero_pd()}};
                } else if (instruction.op == ExpressionProgram::Op::Norm || instruction.op == ExpressionProgram::Op::Abs) {
                    VectorReal realSquared;
                    VectorReal imaginarySquared;
                    if (!vectorSquare(stack[top - 1].real, realSquared) || !vectorSquare(stack[top - 1].imaginary, imaginarySquared) || !vectorAdd(realSquared, imaginarySquared, value.real)) return false;
                    if (instruction.op == ExpressionProgram::Op::Abs) {
                        alignas(32) double nH[4];
                        _mm256_store_pd(nH, value.real.high);
                        alignas(32) double oH[4];
                        for (int lane = 0; lane < 4; ++lane) oH[lane] = std::sqrt(nH[lane]);
                        value.real = {_mm256_load_pd(oH), _mm256_setzero_pd()};
                    }
                    value.imaginary = vectorFromDouble(0.0);
                } else if (instruction.op == ExpressionProgram::Op::Conjugate) {
                    const __m256d signMask = _mm256_set1_pd(-0.0);
                    value = stack[top - 1];
                    value.imaginary.high = _mm256_xor_pd(value.imaginary.high, signMask);
                    value.imaginary.low = _mm256_xor_pd(value.imaginary.low, signMask);
                } else if (instruction.op == ExpressionProgram::Op::Real) {
                    value = {stack[top - 1].real, vectorFromDouble(0.0)};
                } else {
                    value = {stack[top - 1].imaginary, vectorFromDouble(0.0)};
                }
                stack[top - 1] = value;
                break;
            default: return false;
            }
        }
        if (top != 1) return false;
        output = stack[0];
        return true;
    }

    static Real fromDouble(double value) {
        return {value, 0.0};
    }

    static ComplexValue fromComplex(Complex value) {
        return {fromDouble(value.real()), fromDouble(value.imag())};
    }

    static bool finite(Real value) {
        return std::isfinite(value.high) && std::isfinite(value.low);
    }

    static bool finite(ComplexValue value) {
        return finite(value.real) && finite(value.imaginary);
    }

    static Real normalize(double high, double low) {
        const double sum = high + low;
        return {sum, low - (sum - high)};
    }

    static Real twoSum(double left, double right) {
        const double sum = left + right;
        const double virtualRight = sum - left;
        return {sum, (left - (sum - virtualRight)) + (right - virtualRight)};
    }

    static bool add(Real left, Real right, Real& output) {
        if (!finite(left) || !finite(right)) return false;
        Real high = twoSum(left.high, right.high);
        const Real low = twoSum(left.low, right.low);
        high.low += low.high;
        high = normalize(high.high, high.low);
        high.low += low.low;
        output = normalize(high.high, high.low);
        return finite(output);
    }

    static bool subtract(Real left, Real right, Real& output) {
        return add(left, {-right.high, -right.low}, output);
    }

    static bool multiply(Real left, Real right, Real& output) {
        if (!finite(left) || !finite(right)) return false;
        const double product = left.high * right.high;
        double error = std::fma(left.high, right.high, -product);
        error += left.high * right.low + left.low * right.high;
        output = normalize(product, error);
        output.low += left.low * right.low;
        output = normalize(output.high, output.low);
        return finite(output);
    }

    static bool square(Real value, Real& output) {
        return multiply(value, value, output);
    }

    static bool divide(Real numerator, Real denominator, Real& output) {
        if (!finite(numerator) || !finite(denominator)) return false;
        numerator = normalize(numerator.high, numerator.low);
        denominator = normalize(denominator.high, denominator.low);
        if (denominator.high == 0.0) return false;
        const double firstQuotient = numerator.high / denominator.high;
        Real product;
        Real remainder;
        if (!std::isfinite(firstQuotient) || !multiply(denominator, fromDouble(firstQuotient), product) || !subtract(numerator, product, remainder)) return false;
        const double secondQuotient = remainder.high / denominator.high;
        if (!std::isfinite(secondQuotient) || !multiply(denominator, fromDouble(secondQuotient), product) || !subtract(remainder, product, remainder)) return false;
        const double thirdQuotient = (remainder.high + remainder.low) / denominator.high;
        Real first;
        return std::isfinite(thirdQuotient) && add(fromDouble(firstQuotient), fromDouble(secondQuotient), first) && add(first, fromDouble(thirdQuotient), output);
    }

    static bool divide(Real numerator, double denominator, Real& output) {
        return divide(numerator, fromDouble(denominator), output);
    }

    static bool scale(Real value, int exponent, Real& output) {
        if (!finite(value)) return false;
        output = {std::scalbn(value.high, exponent), std::scalbn(value.low, exponent)};
        return finite(output) && ((value.high == 0.0 && value.low == 0.0) || output.high != 0.0 || output.low != 0.0);
    }

    static bool maximumExponent(ComplexValue value, int& exponent) {
        if (!finite(value)) return false;
        const double magnitude = std::max({std::fabs(value.real.high), std::fabs(value.real.low), std::fabs(value.imaginary.high), std::fabs(value.imaginary.low)});
        if (magnitude == 0.0) return false;
        exponent = std::ilogb(magnitude);
        return exponent != FP_ILOGB0 && exponent != FP_ILOGBNAN;
    }

    static bool sineCosine(Real value, Real& sineOutput, Real& cosineOutput) {
        if (!finite(value) || std::fabs(value.high) > 0x1p20) return false;
        constexpr Real PiOverTwo{0x1.921fb54442d18p+0, 0x1.1a62633145c07p-54};
        constexpr double TwoOverPi = 0x1.45f306dc9c883p-1;
        const double quadrantValue = std::nearbyint(value.high * TwoOverPi + value.low * TwoOverPi);
        if (!std::isfinite(quadrantValue) || std::fabs(quadrantValue) > 0x1p52) return false;
        const long long quadrant = static_cast<long long>(quadrantValue);
        Real quadrantOffset;
        Real reduced;
        if (!multiply(PiOverTwo, fromDouble(quadrantValue), quadrantOffset) || !subtract(value, quadrantOffset, reduced)) return false;
        Real reducedSquared;
        if (!square(reduced, reducedSquared)) return false;

        Real sine = reduced;
        Real sineTerm = reduced;
        Real cosine = fromDouble(1.0);
        Real cosineTerm = fromDouble(1.0);
        Real negativeSquared{-reducedSquared.high, -reducedSquared.low};
        for (int term = 1; term <= 18; ++term) {
            Real next;
            if (!multiply(sineTerm, negativeSquared, next) || !divide(next, static_cast<double>((2 * term) * (2 * term + 1)), sineTerm) || !add(sine, sineTerm, sine) || !multiply(cosineTerm, negativeSquared, next) || !divide(next, static_cast<double>((2 * term - 1) * (2 * term)), cosineTerm) || !add(cosine, cosineTerm, cosine)) return false;
        }
        switch ((quadrant % 4 + 4) % 4) {
        case 0:
            sineOutput = sine;
            cosineOutput = cosine;
            break;
        case 1:
            sineOutput = cosine;
            cosineOutput = {-sine.high, -sine.low};
            break;
        case 2:
            sineOutput = {-sine.high, -sine.low};
            cosineOutput = {-cosine.high, -cosine.low};
            break;
        case 3:
            sineOutput = {-cosine.high, -cosine.low};
            cosineOutput = sine;
            break;
        }
        return finite(sineOutput) && finite(cosineOutput);
    }

    static bool sine(Real value, Real& output) {
        Real cosineOutput;
        return sineCosine(value, output, cosineOutput);
    }

    static bool cosine(Real value, Real& output) {
        Real sineOutput;
        return sineCosine(value, sineOutput, output);
    }

    static bool exponential(Real value, Real& output) {
        if (!finite(value) || value.high > 709.0 || value.high < -745.0) return false;
        constexpr Real LnTwo{0x1.62e42fefa39efp-1, 0x1.abc9e3b39803fp-56};
        constexpr double LogTwoE = 0x1.71547652b82fep+0;
        const double exponentValue = std::nearbyint(value.high * LogTwoE + value.low * LogTwoE);
        if (!std::isfinite(exponentValue) || std::fabs(exponentValue) > 0x1p20) return false;
        const int exponent = static_cast<int>(exponentValue);
        Real exponentOffset;
        Real reduced;
        if (!multiply(LnTwo, fromDouble(exponentValue), exponentOffset) || !subtract(value, exponentOffset, reduced)) return false;
        Real sum = fromDouble(1.0);
        Real term = fromDouble(1.0);
        for (int order = 1; order <= 40; ++order) {
            if (!multiply(term, reduced, term) || !divide(term, static_cast<double>(order), term) || !add(sum, term, sum)) return false;
        }
        output = {std::ldexp(sum.high, exponent), std::ldexp(sum.low, exponent)};
        return finite(output) && (sum.high == 0.0 || output.high != 0.0 || output.low != 0.0);
    }

    static bool hyperbolicSine(Real value, Real& output) {
        Real positive;
        Real negative;
        Real difference;
        const Real negated{-value.high, -value.low};
        return exponential(value, positive) && exponential(negated, negative) && subtract(positive, negative, difference) && multiply(difference, fromDouble(0.5), output);
    }

    static bool hyperbolicCosine(Real value, Real& output) {
        Real positive;
        Real negative;
        Real sum;
        const Real negated{-value.high, -value.low};
        return exponential(value, positive) && exponential(negated, negative) && add(positive, negative, sum) && multiply(sum, fromDouble(0.5), output);
    }

    static bool add(ComplexValue left, ComplexValue right, ComplexValue& output) {
        return add(left.real, right.real, output.real) && add(left.imaginary, right.imaginary, output.imaginary);
    }

    static bool subtract(ComplexValue left, ComplexValue right, ComplexValue& output) {
        return subtract(left.real, right.real, output.real) && subtract(left.imaginary, right.imaginary, output.imaginary);
    }

    static bool multiply(ComplexValue left, ComplexValue right, ComplexValue& output) {
        Real realLeft;
        Real realRight;
        Real imaginaryLeft;
        Real imaginaryRight;
        return multiply(left.real, right.real, realLeft) && multiply(left.imaginary, right.imaginary, realRight) && subtract(realLeft, realRight, output.real) && multiply(left.real, right.imaginary, imaginaryLeft) && multiply(left.imaginary, right.real, imaginaryRight) && add(imaginaryLeft, imaginaryRight, output.imaginary);
    }

    static bool square(ComplexValue value, ComplexValue& output) {
        Real realSquared;
        Real imaginarySquared;
        Real product;
        return square(value.real, realSquared) && square(value.imaginary, imaginarySquared) && subtract(realSquared, imaginarySquared, output.real) && multiply(value.real, value.imaginary, product) && add(product, product, output.imaginary);
    }

    static bool divide(ComplexValue numerator, ComplexValue denominator, ComplexValue& output) {
        int numeratorExponent = 0;
        int denominatorExponent = 0;
        if (!maximumExponent(denominator, denominatorExponent)) return false;
        const bool zeroNumerator = numerator.real.high == 0.0 && numerator.real.low == 0.0 && numerator.imaginary.high == 0.0 && numerator.imaginary.low == 0.0;
        if (!zeroNumerator && !maximumExponent(numerator, numeratorExponent)) return false;
        ComplexValue scaledNumerator;
        ComplexValue scaledDenominator;
        if (zeroNumerator) {
            scaledNumerator = {};
            numeratorExponent = denominatorExponent;
        } else if (!scale(numerator.real, -numeratorExponent, scaledNumerator.real) || !scale(numerator.imaginary, -numeratorExponent, scaledNumerator.imaginary)) {
            return false;
        }
        if (!scale(denominator.real, -denominatorExponent, scaledDenominator.real) || !scale(denominator.imaginary, -denominatorExponent, scaledDenominator.imaginary)) return false;

        Real denominatorRealSquared;
        Real denominatorImaginarySquared;
        Real denominatorNorm;
        if (!square(scaledDenominator.real, denominatorRealSquared) || !square(scaledDenominator.imaginary, denominatorImaginarySquared) || !add(denominatorRealSquared, denominatorImaginarySquared, denominatorNorm)) return false;
        denominatorNorm = normalize(denominatorNorm.high, denominatorNorm.low);
        if (!(denominatorNorm.high > std::fabs(denominatorNorm.low))) return false;

        Real first;
        Real second;
        Real realNumerator;
        Real imaginaryNumerator;
        Real normalizedReal;
        Real normalizedImaginary;
        if (!multiply(scaledNumerator.real, scaledDenominator.real, first) || !multiply(scaledNumerator.imaginary, scaledDenominator.imaginary, second) || !add(first, second, realNumerator) || !multiply(scaledNumerator.imaginary, scaledDenominator.real, first) || !multiply(scaledNumerator.real, scaledDenominator.imaginary, second) || !subtract(first, second, imaginaryNumerator) || !divide(realNumerator, denominatorNorm, normalizedReal) || !divide(imaginaryNumerator, denominatorNorm, normalizedImaginary)) return false;
        const int resultExponent = numeratorExponent - denominatorExponent;
        return scale(normalizedReal, resultExponent, output.real) && scale(normalizedImaginary, resultExponent, output.imaginary) && finite(output);
    }

    static bool sine(ComplexValue value, ComplexValue& output) {
        Real sineReal;
        Real cosineReal;
        Real sineImaginary;
        Real cosineImaginary;
        return sine(value.real, sineReal) && cosine(value.real, cosineReal) && hyperbolicSine(value.imaginary, sineImaginary) && hyperbolicCosine(value.imaginary, cosineImaginary) && multiply(sineReal, cosineImaginary, output.real) && multiply(cosineReal, sineImaginary, output.imaginary);
    }

    static bool cosine(ComplexValue value, ComplexValue& output) {
        Real sineReal;
        Real cosineReal;
        Real sineImaginary;
        Real cosineImaginary;
        if (!sine(value.real, sineReal) || !cosine(value.real, cosineReal) || !hyperbolicSine(value.imaginary, sineImaginary) || !hyperbolicCosine(value.imaginary, cosineImaginary) || !multiply(cosineReal, cosineImaginary, output.real) || !multiply(sineReal, sineImaginary, output.imaginary)) return false;
        output.imaginary.high = -output.imaginary.high;
        output.imaginary.low = -output.imaginary.low;
        return true;
    }

    static bool tangent(ComplexValue value, ComplexValue& output) {
        ComplexValue sineValue;
        ComplexValue cosineValue;
        return sine(value, sineValue) && cosine(value, cosineValue) && divide(sineValue, cosineValue, output);
    }

    static bool hyperbolicSine(ComplexValue value, ComplexValue& output) {
        Real sineReal;
        Real cosineReal;
        Real sineImaginary;
        Real cosineImaginary;
        return hyperbolicSine(value.real, sineReal) && hyperbolicCosine(value.real, cosineReal) && sine(value.imaginary, sineImaginary) && cosine(value.imaginary, cosineImaginary) && multiply(sineReal, cosineImaginary, output.real) && multiply(cosineReal, sineImaginary, output.imaginary);
    }

    static bool hyperbolicCosine(ComplexValue value, ComplexValue& output) {
        Real sineReal;
        Real cosineReal;
        Real sineImaginary;
        Real cosineImaginary;
        return hyperbolicSine(value.real, sineReal) && hyperbolicCosine(value.real, cosineReal) && sine(value.imaginary, sineImaginary) && cosine(value.imaginary, cosineImaginary) && multiply(cosineReal, cosineImaginary, output.real) && multiply(sineReal, sineImaginary, output.imaginary);
    }

    static bool hyperbolicTangent(ComplexValue value, ComplexValue& output) {
        ComplexValue sineValue;
        ComplexValue cosineValue;
        return hyperbolicSine(value, sineValue) && hyperbolicCosine(value, cosineValue) && divide(sineValue, cosineValue, output);
    }

    static bool exponential(ComplexValue value, ComplexValue& output) {
        Real magnitude;
        Real sineImaginary;
        Real cosineImaginary;
        return exponential(value.real, magnitude) && sine(value.imaginary, sineImaginary) && cosine(value.imaginary, cosineImaginary) && multiply(magnitude, cosineImaginary, output.real) && multiply(magnitude, sineImaginary, output.imaginary);
    }

    static bool logarithm(ComplexValue value, ComplexValue& output) {
        if (!finite(value)) return false;
        Complex c{value.real.high, value.imaginary.high};
        if (c.real() == 0.0 && c.imag() == 0.0) return false;
        output = fromComplex(std::log(c));
        return finite(output);
    }

    static bool squareRoot(ComplexValue value, ComplexValue& output) {
        if (!finite(value)) return false;
        Complex c{value.real.high, value.imaginary.high};
        output = fromComplex(std::sqrt(c));
        return finite(output);
    }

    static bool power(ComplexValue left, ComplexValue right, ComplexValue& output) {
        if (!finite(left) || !finite(right)) return false;
        Complex l{left.real.high, left.imaginary.high};
        Complex r{right.real.high, right.imaginary.high};
        output = fromComplex(std::pow(l, r));
        return finite(output);
    }

    static BailoutDecision outsideBailout(ComplexValue value, Real bailoutSquared) {
        Real realSquared;
        Real imaginarySquared;
        Real magnitudeSquared;
        if (!finite(value) || !square(value.real, realSquared) || !square(value.imaginary, imaginarySquared) || !add(realSquared, imaginarySquared, magnitudeSquared)) return BailoutDecision::Invalid;
        if (magnitudeSquared.high > bailoutSquared.high || (magnitudeSquared.high == bailoutSquared.high && magnitudeSquared.low > bailoutSquared.low)) return BailoutDecision::Outside;
        return BailoutDecision::Inside;
    }

    static bool evaluateStep(const ExpressionProgram& program, const ExpressionContext& fixed, ComplexValue pixelC, ComplexValue pixelZ0, ComplexValue state, int iteration, ComplexValue& output) {
        std::array<ComplexValue, ExpressionProgram::MAX_STACK> stack;
        size_t top = 0;
        auto push = [&](ComplexValue value) {
            if (top >= stack.size() || !finite(value)) return false;
            stack[top++] = value;
            return true;
        };
        for (const ExpressionProgram::Instruction& instruction : program._code) {
            ComplexValue value;
            switch (instruction.op) {
            case ExpressionProgram::Op::Constant:
                if (!push(fromComplex(instruction.value))) return false;
                break;
            case ExpressionProgram::Op::Z:
                if (!push(state)) return false;
                break;
            case ExpressionProgram::Op::C:
                if (!push(pixelC)) return false;
                break;
            case ExpressionProgram::Op::Z0:
                if (!push(pixelZ0)) return false;
                break;
            case ExpressionProgram::Op::Iteration:
                if (!push({fromDouble(static_cast<double>(iteration)), fromDouble(0.0)})) return false;
                break;
            case ExpressionProgram::Op::Parameter:
                if (instruction.argument >= fixed.parameters.size() || !push(fromComplex(fixed.parameters[instruction.argument]))) return false;
                break;
            case ExpressionProgram::Op::Negate:
                if (top < 1) return false;
                stack[top - 1].real = {-stack[top - 1].real.high, -stack[top - 1].real.low};
                stack[top - 1].imaginary = {-stack[top - 1].imaginary.high, -stack[top - 1].imaginary.low};
                break;
            case ExpressionProgram::Op::Add:
            case ExpressionProgram::Op::Subtract:
            case ExpressionProgram::Op::Multiply:
            case ExpressionProgram::Op::Divide:
            case ExpressionProgram::Op::Power:
            case ExpressionProgram::Op::MakeComplex:
                if (top < 2) return false;
                if (instruction.op == ExpressionProgram::Op::Add) {
                    if (!add(stack[top - 2], stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Subtract) {
                    if (!subtract(stack[top - 2], stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::MakeComplex) {
                    value = {stack[top - 2].real, stack[top - 1].real};
                } else if (instruction.op == ExpressionProgram::Op::Divide) {
                    if (!divide(stack[top - 2], stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Power) {
                    if (!power(stack[top - 2], stack[top - 1], value)) return false;
                } else if (!multiply(stack[top - 2], stack[top - 1], value)) {
                    return false;
                }
                stack[--top - 1] = value;
                break;
            case ExpressionProgram::Op::Square:
            case ExpressionProgram::Op::Sin:
            case ExpressionProgram::Op::Cos:
            case ExpressionProgram::Op::Tan:
            case ExpressionProgram::Op::Sinh:
            case ExpressionProgram::Op::Cosh:
            case ExpressionProgram::Op::Tanh:
            case ExpressionProgram::Op::Exp:
            case ExpressionProgram::Op::Log:
            case ExpressionProgram::Op::Log10:
            case ExpressionProgram::Op::Sqrt:
            case ExpressionProgram::Op::Norm:
            case ExpressionProgram::Op::Abs:
            case ExpressionProgram::Op::Conjugate:
            case ExpressionProgram::Op::Real:
            case ExpressionProgram::Op::Imaginary:
                if (top < 1) return false;
                if (instruction.op == ExpressionProgram::Op::Square) {
                    if (!square(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Sin) {
                    if (!sine(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Cos) {
                    if (!cosine(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Tan) {
                    if (!tangent(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Sinh) {
                    if (!hyperbolicSine(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Cosh) {
                    if (!hyperbolicCosine(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Tanh) {
                    if (!hyperbolicTangent(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Exp) {
                    if (!exponential(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Log) {
                    if (!logarithm(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Log10) {
                    if (!logarithm(stack[top - 1], value)) return false;
                    Real log10Val = fromDouble(std::log(10.0));
                    if (!divide(value.real, log10Val, value.real) || !divide(value.imaginary, log10Val, value.imaginary)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Sqrt) {
                    if (!squareRoot(stack[top - 1], value)) return false;
                } else if (instruction.op == ExpressionProgram::Op::Norm || instruction.op == ExpressionProgram::Op::Abs) {
                    Real realSquared;
                    Real imaginarySquared;
                    if (!square(stack[top - 1].real, realSquared) || !square(stack[top - 1].imaginary, imaginarySquared) || !add(realSquared, imaginarySquared, value.real)) return false;
                    if (instruction.op == ExpressionProgram::Op::Abs) {
                        value.real = fromDouble(std::sqrt(value.real.high));
                    }
                    value.imaginary = fromDouble(0.0);
                } else if (instruction.op == ExpressionProgram::Op::Conjugate) {
                    value = stack[top - 1];
                    value.imaginary = {-value.imaginary.high, -value.imaginary.low};
                } else if (instruction.op == ExpressionProgram::Op::Real) {
                    value = {stack[top - 1].real, fromDouble(0.0)};
                } else {
                    value = {stack[top - 1].imaginary, fromDouble(0.0)};
                }
                stack[top - 1] = value;
                break;
            default: return false;
            }
        }
        if (top != 1 || !finite(stack[0])) return false;
        output = stack[0];
        return true;
    }
};

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
        const int fullHeight = request.fullHeight > 0 ? request.fullHeight : request.height;

        // dx/2 = 2/(scale*(width-1)).
        mpfr_mul_ui(temporary, scale, static_cast<unsigned long>(request.width - 1), RND);
        mpfr_ui_div(dxHalf, 2, temporary, RND);

        // dy/2 = 2*height/(scale*width*(height-1)).
        mpfr_mul_ui(temporary, scale, static_cast<unsigned long>(request.width), RND);
        mpfr_mul_ui(temporary, temporary, static_cast<unsigned long>(fullHeight - 1), RND);
        mpfr_ui_div(dyHalf, static_cast<unsigned long>(fullHeight), temporary, RND);
        mpfr_mul_ui(dyHalf, dyHalf, 2, RND);
        return mpfr_number_p(dxHalf) && mpfr_number_p(dyHalf) && !mpfr_zero_p(dxHalf) && !mpfr_zero_p(dyHalf);
    }

    bool coordinate(int x, int y, const ExpressionDeepRenderRequest& request, MpfrComplex& output) {
        const int fullHeight = request.fullHeight > 0 ? request.fullHeight : request.height;
        const int globalY = y + request.rowBase;
        const long centeredX = static_cast<long>(2LL * x - (request.width - 1LL));
        const long centeredY = static_cast<long>(2LL * globalY - (fullHeight - 1LL));
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

struct TransferComplexInterval {
    explicit TransferComplexInterval(mpfr_prec_t precision) : lower(precision), upper(precision) {}

    MpfrComplex lower;
    MpfrComplex upper;
};

struct TransferIntervalWorkspace {
    explicit TransferIntervalWorkspace(mpfr_prec_t precision) : first(precision), second(precision), centerLower(precision), centerUpper(precision), coordinate(precision) {
        mpfr_inits2(precision, radius, candidate, error, temporary, bailoutSquared, realMagnitude, imaginaryMagnitude, magnitudeSquared, (mpfr_ptr)0);
        for (mpfr_ptr value : products) mpfr_init2(value, precision);
    }

    ~TransferIntervalWorkspace() {
        for (mpfr_ptr value : products) mpfr_clear(value);
        mpfr_clears(radius, candidate, error, temporary, bailoutSquared, realMagnitude, imaginaryMagnitude, magnitudeSquared, (mpfr_ptr)0);
    }

    TransferComplexInterval first;
    TransferComplexInterval second;
    MpfrComplex centerLower;
    MpfrComplex centerUpper;
    MpfrComplex coordinate;
    mpfr_t products[8];
    mpfr_t radius;
    mpfr_t candidate;
    mpfr_t error;
    mpfr_t temporary;
    mpfr_t bailoutSquared;
    mpfr_t realMagnitude;
    mpfr_t imaginaryMagnitude;
    mpfr_t magnitudeSquared;
};

enum class CertifiedTransferBuildStatus : uint8_t {
    Accepted,
    Rejected,
    Cancelled,
    ResourceLimit,
    InvalidReference
};

struct CertifiedTransferBuildResult {
    CertifiedTransferBuildStatus status = CertifiedTransferBuildStatus::Rejected;
    int coveredIterations = 0;
    size_t memoryBytes = 0;
    ScaledRealValue finalRadius;
};

void copyTransferInterval(const TransferComplexInterval& input, TransferComplexInterval& output) {
    mpfr_set(output.lower.re, input.lower.re, MPFR_RNDD);
    mpfr_set(output.upper.re, input.upper.re, MPFR_RNDU);
    mpfr_set(output.lower.im, input.lower.im, MPFR_RNDD);
    mpfr_set(output.upper.im, input.upper.im, MPFR_RNDU);
}

bool finiteTransferInterval(const TransferComplexInterval& value) {
    return mpfr_number_p(value.lower.re) && mpfr_number_p(value.upper.re) && mpfr_number_p(value.lower.im) && mpfr_number_p(value.upper.im) && mpfr_cmp(value.lower.re, value.upper.re) <= 0 && mpfr_cmp(value.lower.im, value.upper.im) <= 0;
}

void addTransferInterval(const TransferComplexInterval& left, const TransferComplexInterval& right, TransferComplexInterval& output) {
    mpfr_add(output.lower.re, left.lower.re, right.lower.re, MPFR_RNDD);
    mpfr_add(output.upper.re, left.upper.re, right.upper.re, MPFR_RNDU);
    mpfr_add(output.lower.im, left.lower.im, right.lower.im, MPFR_RNDD);
    mpfr_add(output.upper.im, left.upper.im, right.upper.im, MPFR_RNDU);
}

void subtractTransferInterval(const TransferComplexInterval& left, const TransferComplexInterval& right, TransferComplexInterval& output) {
    mpfr_sub(output.lower.re, left.lower.re, right.upper.re, MPFR_RNDD);
    mpfr_sub(output.upper.re, left.upper.re, right.lower.re, MPFR_RNDU);
    mpfr_sub(output.lower.im, left.lower.im, right.upper.im, MPFR_RNDD);
    mpfr_sub(output.upper.im, left.upper.im, right.lower.im, MPFR_RNDU);
}

void negateTransferInterval(const TransferComplexInterval& input, TransferComplexInterval& output) {
    mpfr_neg(output.lower.re, input.upper.re, MPFR_RNDD);
    mpfr_neg(output.upper.re, input.lower.re, MPFR_RNDU);
    mpfr_neg(output.lower.im, input.upper.im, MPFR_RNDD);
    mpfr_neg(output.upper.im, input.lower.im, MPFR_RNDU);
}

void multiplyTransferReal(mpfr_srcptr leftLower, mpfr_srcptr leftUpper, mpfr_srcptr rightLower, mpfr_srcptr rightUpper, mpfr_ptr outputLower, mpfr_ptr outputUpper, TransferIntervalWorkspace& workspace) {
    mpfr_mul(workspace.products[0], leftLower, rightLower, MPFR_RNDD);
    mpfr_mul(workspace.products[1], leftLower, rightUpper, MPFR_RNDD);
    mpfr_mul(workspace.products[2], leftUpper, rightLower, MPFR_RNDD);
    mpfr_mul(workspace.products[3], leftUpper, rightUpper, MPFR_RNDD);
    mpfr_set(outputLower, workspace.products[0], MPFR_RNDD);
    for (size_t index = 1; index < 4; ++index)
        if (mpfr_cmp(workspace.products[index], outputLower) < 0) mpfr_set(outputLower, workspace.products[index], MPFR_RNDD);

    mpfr_mul(workspace.products[4], leftLower, rightLower, MPFR_RNDU);
    mpfr_mul(workspace.products[5], leftLower, rightUpper, MPFR_RNDU);
    mpfr_mul(workspace.products[6], leftUpper, rightLower, MPFR_RNDU);
    mpfr_mul(workspace.products[7], leftUpper, rightUpper, MPFR_RNDU);
    mpfr_set(outputUpper, workspace.products[4], MPFR_RNDU);
    for (size_t index = 5; index < 8; ++index)
        if (mpfr_cmp(workspace.products[index], outputUpper) > 0) mpfr_set(outputUpper, workspace.products[index], MPFR_RNDU);
}

void multiplyTransferInterval(const TransferComplexInterval& left, const TransferComplexInterval& right, TransferComplexInterval& output, TransferIntervalWorkspace& workspace) {
    multiplyTransferReal(left.lower.re, left.upper.re, right.lower.re, right.upper.re, workspace.first.lower.re, workspace.first.upper.re, workspace);
    multiplyTransferReal(left.lower.im, left.upper.im, right.lower.im, right.upper.im, workspace.first.lower.im, workspace.first.upper.im, workspace);
    multiplyTransferReal(left.lower.re, left.upper.re, right.lower.im, right.upper.im, workspace.second.lower.re, workspace.second.upper.re, workspace);
    multiplyTransferReal(left.lower.im, left.upper.im, right.lower.re, right.upper.re, workspace.second.lower.im, workspace.second.upper.im, workspace);
    mpfr_sub(output.lower.re, workspace.first.lower.re, workspace.first.upper.im, MPFR_RNDD);
    mpfr_sub(output.upper.re, workspace.first.upper.re, workspace.first.lower.im, MPFR_RNDU);
    mpfr_add(output.lower.im, workspace.second.lower.re, workspace.second.lower.im, MPFR_RNDD);
    mpfr_add(output.upper.im, workspace.second.upper.re, workspace.second.upper.im, MPFR_RNDU);
}

bool setTransferReferenceInterval(const ScaledComplexShadow& primary, const ScaledComplexShadow& defect, const ScaledComplexExpansionTail& tail, const ScaledRealValue& referenceError, TransferComplexInterval& output, TransferIntervalWorkspace& workspace) {
    if (!reconstructMpfrFromExpansion(workspace.centerLower, primary, defect, &tail, MPFR_RNDD) || !reconstructMpfrFromExpansion(workspace.centerUpper, primary, defect, &tail, MPFR_RNDU) || !setMpfrFromScaledValue(workspace.error, referenceError, MPFR_RNDU)) return false;
    mpfr_sub(output.lower.re, workspace.centerLower.re, workspace.error, MPFR_RNDD);
    mpfr_add(output.upper.re, workspace.centerUpper.re, workspace.error, MPFR_RNDU);
    mpfr_sub(output.lower.im, workspace.centerLower.im, workspace.error, MPFR_RNDD);
    mpfr_add(output.upper.im, workspace.centerUpper.im, workspace.error, MPFR_RNDU);
    return finiteTransferInterval(output);
}

bool recenterTransferInterval(const TransferComplexInterval& computed, const ScaledComplexShadow& primary, const ScaledComplexShadow& defect, const ScaledComplexExpansionTail& tail, const ScaledRealValue& referenceError, TransferComplexInterval& output, TransferIntervalWorkspace& workspace) {
    if (!finiteTransferInterval(computed) || !reconstructMpfrFromExpansion(workspace.centerLower, primary, defect, &tail, MPFR_RNDD) || !reconstructMpfrFromExpansion(workspace.centerUpper, primary, defect, &tail, MPFR_RNDU) || !setMpfrFromScaledValue(workspace.radius, referenceError, MPFR_RNDU)) return false;
    auto includeDistance = [&](mpfr_srcptr endpoint, mpfr_srcptr center) {
        if (mpfr_cmp(endpoint, center) >= 0)
            mpfr_sub(workspace.candidate, endpoint, center, MPFR_RNDU);
        else
            mpfr_sub(workspace.candidate, center, endpoint, MPFR_RNDU);
        if (mpfr_cmp(workspace.candidate, workspace.radius) > 0) mpfr_set(workspace.radius, workspace.candidate, MPFR_RNDU);
    };
    includeDistance(computed.lower.re, workspace.centerUpper.re);
    includeDistance(computed.upper.re, workspace.centerLower.re);
    includeDistance(computed.lower.im, workspace.centerUpper.im);
    includeDistance(computed.upper.im, workspace.centerLower.im);
    mpfr_sub(output.lower.re, workspace.centerLower.re, workspace.radius, MPFR_RNDD);
    mpfr_add(output.upper.re, workspace.centerUpper.re, workspace.radius, MPFR_RNDU);
    mpfr_sub(output.lower.im, workspace.centerLower.im, workspace.radius, MPFR_RNDD);
    mpfr_add(output.upper.im, workspace.centerUpper.im, workspace.radius, MPFR_RNDU);
    return finiteTransferInterval(output);
}

bool transferIntervalInsideBailout(const TransferComplexInterval& value, TransferIntervalWorkspace& workspace) {
    mpfr_abs(workspace.realMagnitude, value.lower.re, MPFR_RNDU);
    mpfr_abs(workspace.temporary, value.upper.re, MPFR_RNDU);
    if (mpfr_cmp(workspace.temporary, workspace.realMagnitude) > 0) mpfr_set(workspace.realMagnitude, workspace.temporary, MPFR_RNDU);
    mpfr_abs(workspace.imaginaryMagnitude, value.lower.im, MPFR_RNDU);
    mpfr_abs(workspace.temporary, value.upper.im, MPFR_RNDU);
    if (mpfr_cmp(workspace.temporary, workspace.imaginaryMagnitude) > 0) mpfr_set(workspace.imaginaryMagnitude, workspace.temporary, MPFR_RNDU);
    mpfr_sqr(workspace.realMagnitude, workspace.realMagnitude, MPFR_RNDU);
    mpfr_sqr(workspace.imaginaryMagnitude, workspace.imaginaryMagnitude, MPFR_RNDU);
    mpfr_add(workspace.magnitudeSquared, workspace.realMagnitude, workspace.imaginaryMagnitude, MPFR_RNDU);
    return mpfr_number_p(workspace.magnitudeSquared) && mpfr_cmp(workspace.magnitudeSquared, workspace.bailoutSquared) <= 0;
}

bool estimateCertifiedTransferBytes(size_t nodeCount, mpfr_prec_t precision, size_t& bytes) {
    const size_t limbs = (static_cast<size_t>(precision) + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
    size_t mpfrBytes = sizeof(__mpfr_struct);
    if (!checkedAddSize(mpfrBytes, limbs + 1, sizeof(mp_limb_t))) return false;
    size_t intervalBytes = 0;
    if (!checkedAddSize(intervalBytes, 4, mpfrBytes)) return false;
    bytes = 0;
    return checkedAddSize(bytes, nodeCount, intervalBytes) && checkedAddSize(bytes, 46, mpfrBytes) && checkedAddSize(bytes, 1, sizeof(std::vector<TransferComplexInterval>));
}

template <typename Cancel>
CertifiedTransferBuildResult buildCertifiedTerminalTransfer(const ExpressionDeepRenderRequest& request, const ExpressionReferenceOrbitResult& reference, ExactGeometry& certificationGeometry, uint64_t pixelCount, size_t availableMemoryBytes, Cancel&& shouldCancel) {
    CertifiedTransferBuildResult result;
    if (request.runtimeProgram->scaledResidualCapability() != ExpressionScaledResidualCapability::ExactCenteredArithmetic || request.pixelParameter != FormulaParameter::C || reference.compaction != ExpressionReferenceCompaction::FourTermCertifiedTransfer || !reference.fourTerm || !reference.certifiedAgainstHigherPrecision || reference.certificationPrecision == 0 || reference.samples.size() < static_cast<size_t>(request.maxIterations) || reference.fourTerm->samples.size() != reference.samples.size() || reference.fourTerm->tape.size() != reference.tape.size() || reference.escaped || reference.undefined) {
        result.status = CertifiedTransferBuildStatus::InvalidReference;
        return result;
    }
    const uint64_t instructionCount = static_cast<uint64_t>(request.runtimeProgram->instructionCount());
    if (instructionCount == 0 || pixelCount < request.transfer.minimumPixelCount || request.maxIterations < request.transfer.minimumIterations) return result;
    const long double baselineWork = static_cast<long double>(pixelCount) * request.maxIterations * instructionCount;
    const long double segmentWork = static_cast<long double>(request.maxIterations) * instructionCount * 48.0L + pixelCount;
    if (request.transfer.requirePredictedBenefit && !(segmentWork < baselineWork)) return result;

    if (!estimateCertifiedTransferBytes(request.runtimeProgram->instructionCount(), reference.certificationPrecision, result.memoryBytes)) {
        result.status = CertifiedTransferBuildStatus::ResourceLimit;
        return result;
    }
    if (availableMemoryBytes != 0 && result.memoryBytes > availableMemoryBytes) {
        result.status = CertifiedTransferBuildStatus::ResourceLimit;
        return result;
    }

    TransferIntervalWorkspace workspace(reference.certificationPrecision);
    mpfr_set_d(workspace.bailoutSquared, request.bailout, MPFR_RNDD);
    mpfr_sqr(workspace.bailoutSquared, workspace.bailoutSquared, MPFR_RNDD);
    if (!mpfr_number_p(workspace.bailoutSquared) || mpfr_sgn(workspace.bailoutSquared) <= 0) {
        result.status = CertifiedTransferBuildStatus::InvalidReference;
        return result;
    }

    TransferComplexInterval pixel(reference.certificationPrecision);
    bool haveReal = false;
    for (int x = 0; x < request.width; ++x) {
        if ((x & 63) == 0 && shouldCancel()) {
            result.status = CertifiedTransferBuildStatus::Cancelled;
            return result;
        }
        if (!certificationGeometry.coordinate(x, 0, request, workspace.coordinate)) {
            result.status = CertifiedTransferBuildStatus::InvalidReference;
            return result;
        }
        if (!haveReal) {
            mpfr_set(pixel.lower.re, workspace.coordinate.re, MPFR_RNDD);
            mpfr_set(pixel.upper.re, workspace.coordinate.re, MPFR_RNDU);
            haveReal = true;
        } else {
            if (mpfr_cmp(workspace.coordinate.re, pixel.lower.re) < 0) mpfr_set(pixel.lower.re, workspace.coordinate.re, MPFR_RNDD);
            if (mpfr_cmp(workspace.coordinate.re, pixel.upper.re) > 0) mpfr_set(pixel.upper.re, workspace.coordinate.re, MPFR_RNDU);
        }
    }
    bool haveImaginary = false;
    for (int y = 0; y < request.height; ++y) {
        if ((y & 63) == 0 && shouldCancel()) {
            result.status = CertifiedTransferBuildStatus::Cancelled;
            return result;
        }
        if (!certificationGeometry.coordinate(0, y, request, workspace.coordinate)) {
            result.status = CertifiedTransferBuildStatus::InvalidReference;
            return result;
        }
        if (!haveImaginary) {
            mpfr_set(pixel.lower.im, workspace.coordinate.im, MPFR_RNDD);
            mpfr_set(pixel.upper.im, workspace.coordinate.im, MPFR_RNDU);
            haveImaginary = true;
        } else {
            if (mpfr_cmp(workspace.coordinate.im, pixel.lower.im) < 0) mpfr_set(pixel.lower.im, workspace.coordinate.im, MPFR_RNDD);
            if (mpfr_cmp(workspace.coordinate.im, pixel.upper.im) > 0) mpfr_set(pixel.upper.im, workspace.coordinate.im, MPFR_RNDU);
        }
    }
    if (!finiteTransferInterval(pixel)) {
        result.status = CertifiedTransferBuildStatus::InvalidReference;
        return result;
    }

    TransferComplexInterval fixedC(reference.certificationPrecision);
    TransferComplexInterval fixedZ0(reference.certificationPrecision);
    TransferComplexInterval state(reference.certificationPrecision);
    if (!setTransferReferenceInterval(reference.c, reference.cDefect, reference.fourTerm->c, reference.cError, fixedC, workspace) || !setTransferReferenceInterval(reference.z0, reference.z0Defect, reference.fourTerm->z0, reference.z0Error, fixedZ0, workspace) || !setTransferReferenceInterval(reference.initialZ, reference.initialZDefect, reference.fourTerm->initialZ, reference.initialZError, state, workspace)) {
        result.status = CertifiedTransferBuildStatus::InvalidReference;
        return result;
    }
    if (!transferIntervalInsideBailout(state, workspace)) return result;

    std::vector<TransferComplexInterval> nodes;
    nodes.reserve(request.runtimeProgram->instructionCount());
    for (size_t index = 0; index < request.runtimeProgram->instructionCount(); ++index) nodes.emplace_back(reference.certificationPrecision);

    for (int iteration = 0; iteration < request.maxIterations; ++iteration) {
        if ((iteration & 7) == 0 && shouldCancel()) {
            result.status = CertifiedTransferBuildStatus::Cancelled;
            return result;
        }
        const ExpressionReferenceSample& sample = reference.samples[static_cast<size_t>(iteration)];
        if (sample.iteration != iteration || sample.tapeOffset > reference.tape.size() || sample.tapeCount != nodes.size() || sample.tapeCount > reference.tape.size() - static_cast<size_t>(sample.tapeOffset) || sample.rootNode >= sample.tapeCount) {
            result.status = CertifiedTransferBuildStatus::InvalidReference;
            return result;
        }
        for (size_t local = 0; local < sample.tapeCount; ++local) {
            const size_t tapeIndex = static_cast<size_t>(sample.tapeOffset) + local;
            const ExpressionReferenceTapeNode& node = reference.tape[tapeIndex];
            const ExpressionReferenceTapeNodeFourTerm& fourTermNode = reference.fourTerm->tape[tapeIndex];
            TransferComplexInterval& output = nodes[local];
            auto validUnary = [&] { return node.leftNode < local; };
            auto validBinary = [&] { return node.leftNode < local && node.rightNode < local; };
            bool computed = true;
            switch (node.operation) {
            case ExpressionOracleOperation::Z:
                computed = recenterTransferInterval(state, node.output, node.outputDefect, fourTermNode.output, node.outputError, output, workspace);
                break;
            case ExpressionOracleOperation::C:
                computed = recenterTransferInterval(request.pixelParameter == FormulaParameter::C ? pixel : fixedC, node.output, node.outputDefect, fourTermNode.output, node.outputError, output, workspace);
                break;
            case ExpressionOracleOperation::Z0:
                computed = recenterTransferInterval(request.pixelParameter == FormulaParameter::InitialZ ? pixel : fixedZ0, node.output, node.outputDefect, fourTermNode.output, node.outputError, output, workspace);
                break;
            case ExpressionOracleOperation::Constant:
            case ExpressionOracleOperation::Iteration:
            case ExpressionOracleOperation::Parameter:
                computed = setTransferReferenceInterval(node.output, node.outputDefect, fourTermNode.output, node.outputError, output, workspace);
                break;
            case ExpressionOracleOperation::Negate:
                if (!validUnary()) {
                    computed = false;
                    break;
                }
                negateTransferInterval(nodes[node.leftNode], workspace.first);
                computed = recenterTransferInterval(workspace.first, node.output, node.outputDefect, fourTermNode.output, node.outputError, output, workspace);
                break;
            case ExpressionOracleOperation::Add:
            case ExpressionOracleOperation::Subtract:
                if (!validBinary()) {
                    computed = false;
                    break;
                }
                if (node.operation == ExpressionOracleOperation::Add)
                    addTransferInterval(nodes[node.leftNode], nodes[node.rightNode], workspace.first);
                else
                    subtractTransferInterval(nodes[node.leftNode], nodes[node.rightNode], workspace.first);
                computed = recenterTransferInterval(workspace.first, node.output, node.outputDefect, fourTermNode.output, node.outputError, output, workspace);
                break;
            case ExpressionOracleOperation::Multiply:
            case ExpressionOracleOperation::Square:
                if (!validUnary() || (node.operation == ExpressionOracleOperation::Multiply && !validBinary())) {
                    computed = false;
                    break;
                }
                multiplyTransferInterval(nodes[node.leftNode], nodes[node.operation == ExpressionOracleOperation::Square ? node.leftNode : node.rightNode], workspace.first, workspace);
                computed = recenterTransferInterval(workspace.first, node.output, node.outputDefect, fourTermNode.output, node.outputError, output, workspace);
                break;
            default: computed = false; break;
            }
            if (!computed || !finiteTransferInterval(output)) {
                result.status = CertifiedTransferBuildStatus::InvalidReference;
                return result;
            }
        }
        copyTransferInterval(nodes[sample.rootNode], state);
        if (!transferIntervalInsideBailout(state, workspace)) return result;
        result.coveredIterations = iteration + 1;
    }
    if (result.coveredIterations != request.maxIterations || makeScaledNonnegativeUpward(workspace.radius, result.finalRadius) != ScaledArithmeticStatus::Success) return result;
    result.status = CertifiedTransferBuildStatus::Accepted;
    return result;
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
    const int fullHeight = request.fullHeight > 0 ? request.fullHeight : request.height;
    if (scaleExponent > 0) { required = std::max<uint64_t>(required, static_cast<uint64_t>(scaleExponent) + ceilLog2(static_cast<uint64_t>(std::max(request.width, fullHeight))) + 8); }
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

enum class Centered4ExperimentalStatus : uint8_t {
    NotApplicable,
    Success,
    Cancelled,
    ResourceLimit,
    UndefinedPixel,
    Failed
};

template <typename Cancel, typename Progress>
Centered4ExperimentalStatus renderAdaptiveCenteredExperimental(const ExpressionDeepRenderRequest& request, ExpressionDeepRenderResult& result, mpfr_prec_t selectedPrecision, mpfr_prec_t fallbackPrecision, int threadCount, size_t outerLiveBytes, Cancel&& shouldCancel, Progress&& progress, std::string& error) {
    struct OracleWorkspaceRelease {
        ~OracleWorkspaceRelease() { ExpressionOracle::releaseThreadWorkspace(); }
    } callerOracleWorkspaceRelease;
    const ExpressionScaledResidualCapability capability = request.runtimeProgram->scaledResidualCapability();
    const bool adaptiveCapability = capability == ExpressionScaledResidualCapability::ExactCenteredArithmetic || capability == ExpressionScaledResidualCapability::CertifiedEntireCandidate || capability == ExpressionScaledResidualCapability::CertifiedMeromorphicCandidate || capability == ExpressionScaledResidualCapability::CertifiedBranchCandidate || capability == ExpressionScaledResidualCapability::CertifiedRealCandidate || capability == ExpressionScaledResidualCapability::CertifiedPiecewiseCandidate;
    const size_t pixelCount = static_cast<size_t>(request.width) * request.height;
    if (!request.centered.enableAdaptiveCandidate || request.forceMpfrFallbackForVerification || (request.pixelParameter != FormulaParameter::C && request.pixelParameter != FormulaParameter::InitialZ) || !adaptiveCapability || !ExpressionDoubleDoubleEvaluator::supports(*request.runtimeProgram) || selectedPrecision > 512 || selectedPrecision > MPFR_PREC_MAX - request.memory.fallbackGuardBits || fallbackPrecision > request.precision.maximumBits || fallbackPrecision > ExpressionReferencePrecisionPolicy::ApplicationMaximumBits || pixelCount < request.centered.minimumPixelCount || request.maxIterations < request.centered.minimumIterations || request.centered.maximumReferences < 2) return Centered4ExperimentalStatus::NotApplicable;
    if (request.memory.memoryLimitBytes != 0 && outerLiveBytes >= request.memory.memoryLimitBytes) return Centered4ExperimentalStatus::NotApplicable;
    const size_t centeredMemoryBudget = request.memory.memoryLimitBytes == 0 ? 0 : request.memory.memoryLimitBytes - outerLiveBytes;
    constexpr mpfr_prec_t ReferencePrecision = 512;
    constexpr double DefaultStateDisagreementThreshold = 0x1p-11;
    constexpr double DefaultPrimaryErrorThreshold = 1e10;
    const double StateDisagreementThreshold = std::isfinite(request.verificationCenteredStateDisagreementThreshold) && request.verificationCenteredStateDisagreementThreshold > 0.0 ? request.verificationCenteredStateDisagreementThreshold : DefaultStateDisagreementThreshold;
    const double PrimaryErrorThreshold = std::isfinite(request.verificationCenteredPrimaryErrorThreshold) && request.verificationCenteredPrimaryErrorThreshold > 0.0 ? request.verificationCenteredPrimaryErrorThreshold : DefaultPrimaryErrorThreshold;
    constexpr double PrimaryMarginThreshold = 0x1p-16;
    constexpr double PrimaryInitialError = 1e-15;
    constexpr double PrimaryOperationErrorScale = 64.0 * std::numeric_limits<double>::epsilon();
    constexpr size_t TraceStateCount = 5;
    constexpr size_t DefaultFallbackValidationSamples = 128;
    constexpr size_t DefaultMandatoryFallbackValidationSamples = 32;
    constexpr mpfr_prec_t MandatoryFallbackCandidatePrecision = 224;
    constexpr std::array<mpfr_prec_t, 4> FallbackPrecisionLadder = {160, 192, 224, 256};
    const mpfr_prec_t requiredCandidatePrecision = std::max<mpfr_prec_t>(MPFR_PREC_MIN, selectedPrecision - request.precision.guardBits);
    mpfr_prec_t candidateFallbackPrecision = fallbackPrecision;
    if (request.verificationCenteredFallbackCandidatePrecision > 0) {
        candidateFallbackPrecision = request.verificationCenteredFallbackCandidatePrecision;
    } else {
        for (mpfr_prec_t precision : FallbackPrecisionLadder) {
            if (precision < requiredCandidatePrecision) continue;
            candidateFallbackPrecision = precision;
            break;
        }
    }
    const bool useLowPrecisionFallback = candidateFallbackPrecision >= MPFR_PREC_MIN && candidateFallbackPrecision < fallbackPrecision;
    if (!useLowPrecisionFallback) candidateFallbackPrecision = fallbackPrecision;
    const size_t maximumFallbackValidationSamples = useLowPrecisionFallback ? std::min(pixelCount, request.verificationCenteredFallbackValidationSamples != 0 ? request.verificationCenteredFallbackValidationSamples : DefaultFallbackValidationSamples) : 0;
    result.centeredFallbackCandidatePrecision = candidateFallbackPrecision;
    result.centeredFallbackFullPrecision = fallbackPrecision;
    result.centeredFallbackValidationIsEmpirical = useLowPrecisionFallback;
    const bool mandatoryPrecisionOverridden = request.verificationCenteredMandatoryFallbackCandidatePrecision > 0;
    mpfr_prec_t mandatoryCandidatePrecision = fallbackPrecision;
    if (mandatoryPrecisionOverridden) {
        mandatoryCandidatePrecision = request.verificationCenteredMandatoryFallbackCandidatePrecision;
    } else if (requiredCandidatePrecision <= MandatoryFallbackCandidatePrecision) {
        mandatoryCandidatePrecision = MandatoryFallbackCandidatePrecision;
    }
    if (mandatoryCandidatePrecision < MPFR_PREC_MIN || mandatoryCandidatePrecision >= fallbackPrecision) mandatoryCandidatePrecision = fallbackPrecision;
    const bool useMandatoryCandidate = mandatoryCandidatePrecision < fallbackPrecision;
    const size_t maximumMandatoryValidationSamples = useMandatoryCandidate ? std::min(pixelCount, mandatoryPrecisionOverridden ? request.verificationCenteredMandatoryFallbackValidationSamples : DefaultMandatoryFallbackValidationSamples) : 0;
    result.centeredFallbackMandatoryCandidatePrecision = mandatoryCandidatePrecision;
    ExactGeometry geometry(ReferencePrecision);
    if (!geometry.initialize(request)) {
        error = "adaptive centered geometry initialization failed";
        return Centered4ExperimentalStatus::Failed;
    }
    result.centeredAttempted = true;

    struct Candidate {
        int x = 0;
        int y = 0;
    };
    std::vector<Candidate> candidates;
    constexpr int GridWidth = 9;
    constexpr int GridHeight = 9;
    for (int gridY = 0; gridY < GridHeight; ++gridY) {
        const int y = static_cast<int>(static_cast<long long>(gridY) * (request.height - 1) / (GridHeight - 1));
        for (int gridX = 0; gridX < GridWidth; ++gridX) {
            const int x = static_cast<int>(static_cast<long long>(gridX) * (request.width - 1) / (GridWidth - 1));
            candidates.push_back({x, y});
        }
    }
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& left, const Candidate& right) { return left.x < right.x || (left.x == right.x && left.y < right.y); });
    candidates.erase(std::unique(candidates.begin(), candidates.end(), [](const Candidate& left, const Candidate& right) { return left.x == right.x && left.y == right.y; }), candidates.end());

    struct AdaptiveReference {
        ExpressionReferenceOrbitResult orbit;
        int x = 0;
        int y = 0;
        Complex c{};
        Complex z0{};
        std::vector<Complex> z;
        std::vector<Complex> nodes;
        std::vector<Complex> nodeAuxiliaries;
    };
    std::vector<AdaptiveReference> references;
    const size_t referenceLimit = std::min(request.centered.maximumReferences, candidates.size());
    references.reserve(referenceLimit);
    MpfrComplex candidatePixel(ReferencePrecision);
    mpf_t candidateReal;
    mpf_t candidateImaginary;
    mpf_init2(candidateReal, ReferencePrecision);
    mpf_init2(candidateImaginary, ReferencePrecision);
    const Clock::time_point referenceStart = Clock::now();
    size_t retainedOrbitBytes = 0;
    while (!candidates.empty() && references.size() < referenceLimit) {
        if (shouldCancel()) {
            mpf_clears(candidateReal, candidateImaginary, (mpf_ptr)0);
            return Centered4ExperimentalStatus::Cancelled;
        }
        size_t selected = 0;
        long long selectedScore = std::numeric_limits<long long>::min();
        for (size_t candidateIndex = 0; candidateIndex < candidates.size(); ++candidateIndex) {
            const Candidate& candidate = candidates[candidateIndex];
            long long score = 0;
            if (references.empty()) {
                const long long dx = 2LL * candidate.x - (request.width - 1LL);
                const long long dy = 2LL * candidate.y - (request.height - 1LL);
                score = -(dx * dx + dy * dy);
            } else {
                score = std::numeric_limits<long long>::max();
                for (const AdaptiveReference& reference : references) {
                    const long long dx = candidate.x - reference.x;
                    const long long dy = candidate.y - reference.y;
                    score = std::min(score, dx * dx + dy * dy);
                }
            }
            if (score > selectedScore) {
                selectedScore = score;
                selected = candidateIndex;
            }
        }
        const Candidate candidate = candidates[selected];
        candidates.erase(candidates.begin() + static_cast<std::ptrdiff_t>(selected));
        ++result.centeredReferenceAttemptCount;
        if (!geometry.coordinate(candidate.x, candidate.y, request, candidatePixel)) continue;
        mpfr_get_f(candidateReal, candidatePixel.re, MPFR_RNDN);
        mpfr_get_f(candidateImaginary, candidatePixel.im, MPFR_RNDN);
        ExpressionReferenceBuildRequest referenceRequest;
        referenceRequest.canonicalProgram = request.canonicalProgram;
        referenceRequest.runtimeProgram = request.runtimeProgram;
        referenceRequest.pixelParameter = request.pixelParameter;
        referenceRequest.center.realMpf = candidateReal;
        referenceRequest.center.imaginaryMpf = candidateImaginary;
        referenceRequest.fixed = request.fixed;
        referenceRequest.bailout = request.bailout;
        referenceRequest.maxIterations = request.maxIterations;
        referenceRequest.precision.requestedBits = ReferencePrecision;
        referenceRequest.precision.minimumBits = 53;
        referenceRequest.precision.guardBits = 0;
        referenceRequest.precision.maximumBits = request.precision.maximumBits;
        if (centeredMemoryBudget != 0 && retainedOrbitBytes >= centeredMemoryBudget) {
            mpf_clears(candidateReal, candidateImaginary, (mpf_ptr)0);
            error = "adaptive centered references exceed the remaining memory budget";
            return Centered4ExperimentalStatus::ResourceLimit;
        }
        referenceRequest.memoryLimitBytes = centeredMemoryBudget == 0 ? 0 : centeredMemoryBudget - retainedOrbitBytes;
        referenceRequest.shouldCancel = shouldCancel;
        AdaptiveReference reference;
        if (!buildExpressionReferenceOrbit(referenceRequest, reference.orbit)) {
            if (reference.orbit.cancelled) {
                mpf_clears(candidateReal, candidateImaginary, (mpf_ptr)0);
                return Centered4ExperimentalStatus::Cancelled;
            }
            continue;
        }
        if (!reference.orbit.valid || reference.orbit.escaped || reference.orbit.undefined || reference.orbit.samples.size() != static_cast<size_t>(request.maxIterations)) continue;
        if (reference.orbit.memoryBytes > std::numeric_limits<size_t>::max() - retainedOrbitBytes || (centeredMemoryBudget != 0 && reference.orbit.memoryBytes > centeredMemoryBudget - retainedOrbitBytes)) {
            error = "adaptive centered retained references exceed the remaining memory budget";
            return Centered4ExperimentalStatus::ResourceLimit;
        }
        retainedOrbitBytes += reference.orbit.memoryBytes;
        reference.x = candidate.x;
        reference.y = candidate.y;
        references.push_back(std::move(reference));
    }
    mpf_clears(candidateReal, candidateImaginary, (mpf_ptr)0);
    result.referenceSeconds += secondsSince(referenceStart);
    ExpressionOracle::releaseThreadWorkspace();
    if (references.size() < 2) {
        result.centeredPreflightRejected = true;
        return Centered4ExperimentalStatus::NotApplicable;
    }
    result.centeredReferenceCount = references.size();

    auto reconstruct = [](const ScaledComplexShadow& primary, const ScaledComplexShadow& defect, Complex& output) {
        ScaledComplexValue scaled;
        return makeScaledComplexValue(primary, defect, scaled) == ScaledArithmeticStatus::Success && scaledValueToDouble(scaled, output);
    };
    size_t referenceBytes = 0;
    for (AdaptiveReference& reference : references) {
        if (!reconstruct(reference.orbit.c, reference.orbit.cDefect, reference.c) || !reconstruct(reference.orbit.z0, reference.orbit.z0Defect, reference.z0)) return Centered4ExperimentalStatus::NotApplicable;
        reference.z.resize(reference.orbit.samples.size());
        reference.nodes.resize(reference.orbit.tape.size());
        reference.nodeAuxiliaries.resize(reference.orbit.tape.size());
        for (size_t index = 0; index < reference.nodes.size(); ++index) {
            if (!reconstruct(reference.orbit.tape[index].output, reference.orbit.tape[index].outputDefect, reference.nodes[index])) return Centered4ExperimentalStatus::NotApplicable;
            if (!reconstruct(reference.orbit.tape[index].auxiliary, reference.orbit.tape[index].auxiliaryDefect, reference.nodeAuxiliaries[index])) return Centered4ExperimentalStatus::NotApplicable;
        }
        for (size_t iteration = 0; iteration < reference.z.size(); ++iteration)
            if (!reconstruct(reference.orbit.samples[iteration].z, reference.orbit.samples[iteration].zDefect, reference.z[iteration])) return Centered4ExperimentalStatus::NotApplicable;
        if (!checkedAddSize(referenceBytes, 1, reference.orbit.memoryBytes) || !checkedAddSize(referenceBytes, reference.z.size() + reference.nodes.size() + reference.nodeAuxiliaries.size(), sizeof(Complex))) {
            error = "adaptive centered reference size overflow";
            return Centered4ExperimentalStatus::ResourceLimit;
        }
    }
    size_t coordinateAndAssignmentBytes = 0;
    if (!checkedAddSize(coordinateAndAssignmentBytes, static_cast<size_t>(request.width + request.height), sizeof(double) + sizeof(ExpressionDoubleDoubleEvaluator::Real)) || !checkedAddSize(coordinateAndAssignmentBytes, pixelCount, 3 * sizeof(size_t)) || !checkedAddSize(coordinateAndAssignmentBytes, references.size(), sizeof(std::vector<size_t>)) || !checkedAddSize(coordinateAndAssignmentBytes, references.capacity(), sizeof(AdaptiveReference))) {
        error = "adaptive centered coordinate and assignment size overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    if (request.memory.memoryLimitBytes != 0 && (outerLiveBytes > request.memory.memoryLimitBytes || referenceBytes > request.memory.memoryLimitBytes - outerLiveBytes || coordinateAndAssignmentBytes > request.memory.memoryLimitBytes - outerLiveBytes - referenceBytes)) return Centered4ExperimentalStatus::NotApplicable;

    std::vector<double> xOffsets;
    std::vector<double> yOffsets;
    std::vector<ExpressionDoubleDoubleEvaluator::Real> xCoordinates;
    std::vector<ExpressionDoubleDoubleEvaluator::Real> yCoordinates;
    try {
        xOffsets.resize(static_cast<size_t>(request.width));
        yOffsets.resize(static_cast<size_t>(request.height));
        xCoordinates.resize(static_cast<size_t>(request.width));
        yCoordinates.resize(static_cast<size_t>(request.height));
    } catch (const std::bad_alloc&) {
        error = "adaptive centered coordinate allocation failed";
        return Centered4ExperimentalStatus::ResourceLimit;
    } catch (const std::length_error&) {
        error = "adaptive centered coordinate length overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    mpfr_t offset;
    mpfr_t coordinate;
    mpfr_inits2(ReferencePrecision, offset, coordinate, (mpfr_ptr)0);
    bool coordinatesReady = true;
    for (int x = 0; x < request.width; ++x) {
        const long centered = static_cast<long>(2LL * x - (request.width - 1LL));
        mpfr_mul_si(offset, geometry.dxHalf, centered, MPFR_RNDN);
        xOffsets[static_cast<size_t>(x)] = mpfr_get_d(offset, MPFR_RNDN);
        mpfr_add(coordinate, geometry.center.re, offset, MPFR_RNDN);
        coordinatesReady = coordinatesReady && std::isfinite(xOffsets[static_cast<size_t>(x)]) && (mpfr_zero_p(offset) || xOffsets[static_cast<size_t>(x)] != 0.0) && ExpressionDoubleDoubleEvaluator::splitMpfr(coordinate, geometry.temporary, xCoordinates[static_cast<size_t>(x)]);
        if (!coordinatesReady) break;
    }
    const int fullHeight = request.fullHeight > 0 ? request.fullHeight : request.height;
    for (int y = 0; coordinatesReady && y < request.height; ++y) {
        const int globalY = y + request.rowBase;
        const long centered = static_cast<long>(2LL * globalY - (fullHeight - 1LL));
        mpfr_mul_si(offset, geometry.dyHalf, centered, MPFR_RNDN);
        yOffsets[static_cast<size_t>(y)] = mpfr_get_d(offset, MPFR_RNDN);
        mpfr_add(coordinate, geometry.center.im, offset, MPFR_RNDN);
        coordinatesReady = coordinatesReady && std::isfinite(yOffsets[static_cast<size_t>(y)]) && (mpfr_zero_p(offset) || yOffsets[static_cast<size_t>(y)] != 0.0) && ExpressionDoubleDoubleEvaluator::splitMpfr(coordinate, geometry.temporary, yCoordinates[static_cast<size_t>(y)]);
    }
    auto sameDoubleDouble = [](const ExpressionDoubleDoubleEvaluator::Real& left, const ExpressionDoubleDoubleEvaluator::Real& right) { return left.high == right.high && left.low == right.low; };
    for (size_t x = 1; coordinatesReady && x < xCoordinates.size(); ++x) coordinatesReady = !sameDoubleDouble(xCoordinates[x - 1], xCoordinates[x]);
    for (size_t y = 1; coordinatesReady && y < yCoordinates.size(); ++y) coordinatesReady = !sameDoubleDouble(yCoordinates[y - 1], yCoordinates[y]);
    mpfr_clears(offset, coordinate, (mpfr_ptr)0);
    if (!coordinatesReady) {
        result.centeredPreflightRejected = true;
        return Centered4ExperimentalStatus::NotApplicable;
    }

    std::vector<size_t> nearest(pixelCount);
    std::vector<size_t> secondNearest(pixelCount);
    std::vector<std::vector<size_t>> primaryAssignments(references.size());
    for (size_t index = 0; index < pixelCount; ++index) {
        const int x = static_cast<int>(index % static_cast<size_t>(request.width));
        const int y = static_cast<int>(index / static_cast<size_t>(request.width));
        size_t first = 0;
        size_t second = 0;
        long long firstDistance = std::numeric_limits<long long>::max();
        long long secondDistance = std::numeric_limits<long long>::max();
        for (size_t referenceIndex = 0; referenceIndex < references.size(); ++referenceIndex) {
            const long long dx = x - references[referenceIndex].x;
            const long long dy = y - references[referenceIndex].y;
            const long long distance = dx * dx + dy * dy;
            if (distance < firstDistance || (distance == firstDistance && referenceIndex < first)) {
                secondDistance = firstDistance;
                second = first;
                firstDistance = distance;
                first = referenceIndex;
            } else if (distance < secondDistance || (distance == secondDistance && referenceIndex < second)) {
                secondDistance = distance;
                second = referenceIndex;
            }
        }
        nearest[index] = first;
        secondNearest[index] = second;
        primaryAssignments[first].push_back(index);
    }

    struct CandidateTrace {
        std::array<Complex, TraceStateCount> states{};
        uint8_t validMask = 0;
    };
    enum PrimaryRiskFlag : uint8_t {
        PrimaryRiskFailure = 1 << 0,
        PrimaryRiskNonfinite = 1 << 1,
        PrimaryRiskPerimeter = 1 << 2,
        PrimaryRiskError = 1 << 3,
        PrimaryRiskMargin = 1 << 4,
        PrimaryRiskInterior = 1 << 5,
        PrimaryRiskDenominator = 1 << 6,
    };
    struct PrimaryRiskTelemetry {
        double currentError = PrimaryInitialError;
        double maximumError = PrimaryInitialError;
        double minimumBailoutMargin = std::numeric_limits<double>::infinity();
        double escapeMargin = std::numeric_limits<double>::infinity();
        uint8_t flags = 0;
    };
    struct PixelEvaluation {
        float output = ExpressionDeepInteriorPixel;
        int escapeIteration = -1;
        CandidateTrace trace;
        PrimaryRiskTelemetry risk;
        bool failed = false;
        bool denominatorInstability = false;
        uint64_t iterations = 0;
    };
    auto initialCandidateStateComplete = [&](const AdaptiveReference& reference, const Complex& stateDelta, auto& candidate) {
        const Complex actual = reference.z[0] + stateDelta;
        if (!finiteComplex(actual)) {
            candidate.failed = true;
            candidate.risk.flags |= PrimaryRiskFailure | PrimaryRiskNonfinite;
            return true;
        }
        const double magnitude = std::abs(actual);
        if (magnitude <= request.bailout) return false;
        candidate.output = 0.0f;
        candidate.escapeIteration = 0;
        candidate.trace.states[4] = actual;
        candidate.trace.validMask |= uint8_t{1} << 4;
        const double margin = (magnitude - request.bailout) / (1.0 + request.bailout);
        candidate.risk.minimumBailoutMargin = margin;
        candidate.risk.escapeMargin = margin;
        return true;
    };
    const std::array<int, 4> probeIterations = {std::max(1, request.maxIterations / 4), std::max(1, request.maxIterations / 2), std::max(1, 3 * request.maxIterations / 4), request.maxIterations};
    auto evaluatePixel = [&](const AdaptiveReference& reference, size_t outputIndex) {
        PixelEvaluation pixel;
        const int x = static_cast<int>(outputIndex % static_cast<size_t>(request.width));
        const int y = static_cast<int>(outputIndex / static_cast<size_t>(request.width));
        ExpressionContext referenceContext = request.fixed;
        referenceContext.c = reference.c;
        referenceContext.z0 = reference.z0;
        ExpressionDeltaContext delta;
        const Complex pixelOffset{xOffsets[static_cast<size_t>(x)] - xOffsets[static_cast<size_t>(reference.x)], yOffsets[static_cast<size_t>(y)] - yOffsets[static_cast<size_t>(reference.y)]};
        if (request.pixelParameter == FormulaParameter::C)
            delta.c = pixelOffset;
        else
            delta.z0 = pixelOffset;
        Complex stateDelta = request.pixelParameter == FormulaParameter::InitialZ ? pixelOffset : Complex{};
        if (initialCandidateStateComplete(reference, stateDelta, pixel)) return pixel;
        for (int iteration = 0; iteration < request.maxIterations; ++iteration) {
            if ((iteration & 63) == 0 && shouldCancel()) {
                pixel.failed = true;
                break;
            }
            referenceContext.iteration = iteration;
            referenceContext.z = reference.z[static_cast<size_t>(iteration)];
            delta.z = stateDelta;
            const ExpressionReferenceSample& sample = reference.orbit.samples[static_cast<size_t>(iteration)];
            const Complex* nodeBases = reference.nodes.data() + static_cast<size_t>(sample.tapeOffset);
            const Complex* nodeAuxiliaries = reference.nodeAuxiliaries.data() + static_cast<size_t>(sample.tapeOffset);
            ExpressionDeltaContext deltas[4]{};
            ExpressionCenteredResult evaluated[4];
            double jacobianNorms[4];
            std::fill(std::begin(deltas), std::end(deltas), delta);
            if (!ExpressionCenteredEvaluator::evaluate4WithNodeBasesFastEntireDerivative(*request.runtimeProgram, referenceContext, nodeBases, nodeAuxiliaries, deltas, evaluated, jacobianNorms)) {
                pixel.failed = true;
                pixel.risk.flags |= PrimaryRiskFailure;
                break;
            }
            ++pixel.iterations;
            if (evaluated[0].denominatorInstability) {
                pixel.denominatorInstability = true;
                pixel.risk.flags |= PrimaryRiskDenominator;
            }
            if (!evaluated[0].success()) {
                pixel.failed = true;
                pixel.risk.flags |= PrimaryRiskFailure;
                break;
            }
            const Complex nextReference = nodeBases[sample.rootNode];
            const Complex actual = nextReference + evaluated[0].delta;
            stateDelta = iteration + 1 < request.maxIterations ? evaluated[0].delta + (nextReference - reference.z[static_cast<size_t>(iteration + 1)]) : evaluated[0].delta;
            const double magnitude = std::abs(actual);
            const double jacobianNorm = jacobianNorms[0];
            if (std::isfinite(jacobianNorm)) {
                pixel.risk.currentError = jacobianNorm * pixel.risk.currentError + PrimaryOperationErrorScale * std::max(1.0, magnitude);
                pixel.risk.maximumError = std::max(pixel.risk.maximumError, pixel.risk.currentError);
            } else {
                pixel.risk.currentError = pixel.risk.maximumError = std::numeric_limits<double>::infinity();
                pixel.risk.flags |= PrimaryRiskNonfinite;
            }
            pixel.trace.states[4] = actual;
            pixel.trace.validMask |= uint8_t{1} << 4;
            for (size_t probe = 0; probe < probeIterations.size(); ++probe) {
                if (iteration + 1 == probeIterations[probe]) {
                    pixel.trace.states[probe] = actual;
                    pixel.trace.validMask |= uint8_t{1} << probe;
                }
            }
            if (std::isfinite(magnitude)) {
                pixel.risk.minimumBailoutMargin = std::min(pixel.risk.minimumBailoutMargin, std::fabs(magnitude - request.bailout) / (1.0 + request.bailout));
            } else {
                pixel.risk.flags |= PrimaryRiskNonfinite;
            }
            if (!finiteComplex(actual) || std::hypot(actual.real(), actual.imag()) > request.bailout) {
                pixel.output = static_cast<float>(iteration + 1);
                pixel.escapeIteration = iteration + 1;
                if (!finiteComplex(actual)) {
                    pixel.failed = true;
                    pixel.risk.flags |= PrimaryRiskFailure | PrimaryRiskNonfinite;
                } else {
                    pixel.risk.escapeMargin = (magnitude - request.bailout) / (1.0 + request.bailout);
                }
                break;
            }
        }
        return pixel;
    };
    auto primaryRiskFlags = [&](const PrimaryRiskTelemetry& telemetry, int escapeIteration, bool failed, size_t index) {
        uint8_t flags = telemetry.flags;
        if (failed) flags |= PrimaryRiskFailure;
        const size_t x = index % static_cast<size_t>(request.width);
        const size_t y = index / static_cast<size_t>(request.width);
        if (x == 0 || x + 1 == static_cast<size_t>(request.width) || y == 0 || y + 1 == static_cast<size_t>(request.height)) flags |= PrimaryRiskPerimeter;
        if (!std::isfinite(telemetry.currentError) || !std::isfinite(telemetry.maximumError) || !std::isfinite(telemetry.minimumBailoutMargin)) flags |= PrimaryRiskNonfinite;
        if (telemetry.maximumError > PrimaryErrorThreshold) flags |= PrimaryRiskError;
        if (escapeIteration >= 0) {
            if (!std::isfinite(telemetry.escapeMargin) || telemetry.escapeMargin < PrimaryMarginThreshold || telemetry.minimumBailoutMargin < PrimaryMarginThreshold) flags |= PrimaryRiskMargin;
        } else {
            const bool independentlySafe = (flags & (PrimaryRiskFailure | PrimaryRiskNonfinite | PrimaryRiskError | PrimaryRiskDenominator)) == 0 && telemetry.currentError <= PrimaryErrorThreshold && telemetry.minimumBailoutMargin >= PrimaryMarginThreshold;
            if (!independentlySafe) flags |= PrimaryRiskInterior;
        }
        return flags;
    };
    auto traceGap = [](const CandidateTrace& primary, const CandidateTrace& secondary) {
        double gap = 0.0;
        const uint8_t common = primary.validMask & secondary.validMask;
        for (size_t state = 0; state < primary.states.size(); ++state) {
            if ((common & (uint8_t{1} << state)) == 0) continue;
            const double scale = 1.0 + std::max(std::abs(primary.states[state]), std::abs(secondary.states[state]));
            gap = std::max(gap, std::abs(primary.states[state] - secondary.states[state]) / scale);
        }
        return gap;
    };
    auto finiteTrace = [](const CandidateTrace& trace) {
        if (trace.validMask == 0) return false;
        for (size_t state = 0; state < trace.states.size(); ++state)
            if ((trace.validMask & (uint8_t{1} << state)) != 0 && !finiteComplex(trace.states[state])) return false;
        return true;
    };
    auto exactPixel = [&](FallbackWorkerWorkspace& workspace, size_t outputIndex, uint64_t& iterations, float& output) {
        const int y = static_cast<int>(outputIndex / static_cast<size_t>(request.width));
        const int x = static_cast<int>(outputIndex % static_cast<size_t>(request.width));
        configureFixed(request.fixed, workspace.context);
        if (!workspace.geometry.coordinate(x, y, request, workspace.pixel)) return false;
        if (request.pixelParameter == FormulaParameter::C)
            workspace.context.c.set(workspace.pixel);
        else
            workspace.context.z0.set(workspace.pixel);
        workspace.context.z.set(workspace.context.z0);
        bool outside = false;
        if (!piecewiseOutsideBailout(workspace.context.z, workspace.piecewiseBailoutSquared, request.bailout, workspace.magnitudeStorage.re, workspace.piecewiseSquares, workspace.piecewiseInsideSquareExponent, workspace.piecewiseScratch, outside)) return false;
        if (outside) {
            output = 0.0f;
            return true;
        }
        output = ExpressionDeepInteriorPixel;
        for (int iteration = 0; iteration < request.maxIterations; ++iteration) {
            if ((iteration & 15) == 0 && shouldCancel()) return false;
            workspace.context.iteration = iteration;
            const bool defined = ExpressionOracle::evaluateOrbitStep(*request.runtimeProgram, workspace.context, workspace.next, &workspace.piecewiseSquares, nullptr);
            if (mpfr_nan_p(workspace.next.re) || mpfr_nan_p(workspace.next.im)) return false;
            ++iterations;
            if (mpfr_inf_p(workspace.next.re) || mpfr_inf_p(workspace.next.im)) {
                output = static_cast<float>(iteration + 1);
                return true;
            }
            if (!defined || !mpfr_number_p(workspace.next.re) || !mpfr_number_p(workspace.next.im)) return false;
            mpfr_swap(workspace.context.z.re, workspace.next.re);
            mpfr_swap(workspace.context.z.im, workspace.next.im);
            if (!piecewiseOutsideBailout(workspace.context.z, workspace.piecewiseBailoutSquared, request.bailout, workspace.magnitudeStorage.re, workspace.piecewiseSquares, workspace.piecewiseInsideSquareExponent, workspace.piecewiseScratch, outside)) return false;
            if (outside) {
                output = static_cast<float>(iteration + 1);
                break;
            }
        }
        return true;
    };
    struct FallbackValidationResult {
        size_t index = 0;
        float lowOutput = ExpressionDeepEmptyPixel;
        float fullOutput = ExpressionDeepEmptyPixel;
        uint64_t lowIterations = 0;
        uint64_t fullIterations = 0;
        bool lowDefined = false;
        bool fullDefined = false;
    };

    const size_t limbs = (static_cast<size_t>(fallbackPrecision) + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
    size_t mpfrValueBytes = sizeof(__mpfr_struct);
    if (!checkedAddSize(mpfrValueBytes, limbs, sizeof(mp_limb_t))) {
        error = "adaptive centered MPFR workspace size overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    size_t mpfrValuesPerThread = 62;
    if (!checkedAddSize(mpfrValuesPerThread, request.runtimeProgram->stackDepth(), 2)) {
        error = "adaptive centered oracle stack size overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    size_t fallbackThreadBytes = 0;
    if (!checkedAddSize(fallbackThreadBytes, mpfrValuesPerThread, mpfrValueBytes)) {
        error = "adaptive centered thread workspace overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    size_t estimatedThreadBytes = fallbackThreadBytes;
    if (!checkedAddSize(estimatedThreadBytes, 1, size_t{64} << 10)) {
        error = "adaptive centered thread allocation size overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    struct PreflightSampleResult {
        size_t index = 0;
        float output = ExpressionDeepEmptyPixel;
        int centeredEscapeIteration = -1;
        bool accepted = false;
        bool doubleDoubleRequired = false;
    };
    const size_t maximumPreflightSamples = std::min(request.centered.preflightSamples, pixelCount);
    size_t preflightTemporaryBytes = 0;
    if (!checkedAddSize(preflightTemporaryBytes, maximumPreflightSamples, sizeof(size_t) + sizeof(PreflightSampleResult))) {
        error = "adaptive centered preflight workspace size overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    size_t preflightLiveBytes = preflightTemporaryBytes;
    if (!checkedAddSize(preflightLiveBytes, static_cast<size_t>(request.width + request.height), sizeof(double) + sizeof(ExpressionDoubleDoubleEvaluator::Real)) || !checkedAddSize(preflightLiveBytes, pixelCount, 3 * sizeof(size_t)) || !checkedAddSize(preflightLiveBytes, references.capacity(), sizeof(AdaptiveReference)) || !checkedAddSize(preflightLiveBytes, references.size(), sizeof(std::vector<size_t>))) {
        error = "adaptive centered live preflight workspace size overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    int preflightThreadCount = threadCount;
    if (request.memory.memoryLimitBytes != 0) {
        if (outerLiveBytes > request.memory.memoryLimitBytes || referenceBytes > request.memory.memoryLimitBytes - outerLiveBytes) return Centered4ExperimentalStatus::NotApplicable;
        const size_t available = request.memory.memoryLimitBytes - outerLiveBytes - referenceBytes;
        if (preflightLiveBytes > available) return Centered4ExperimentalStatus::NotApplicable;
        const size_t affordableThreads = estimatedThreadBytes == 0 ? static_cast<size_t>(threadCount) : (available - preflightLiveBytes) / estimatedThreadBytes;
        if (affordableThreads == 0) return Centered4ExperimentalStatus::NotApplicable;
        preflightThreadCount = static_cast<int>(std::min<size_t>(static_cast<size_t>(threadCount), affordableThreads));
    }

    const Clock::time_point preflightStart = Clock::now();
    size_t preflightFlags = 0;
    std::vector<size_t> preflightSamples;
    std::vector<PreflightSampleResult> preflightEvaluated;
    try {
        preflightSamples = makePreflightSamples(request.width, request.height, request.centered.preflightSamples);
        preflightEvaluated.resize(preflightSamples.size());
    } catch (const std::bad_alloc&) {
        error = "adaptive centered preflight allocation failed";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    std::atomic<uint64_t> preflightSampleCount{0};
    std::atomic<uint64_t> preflightPrimaryRiskCount{0};
    std::atomic<uint64_t> preflightSecondaryCount{0};
    std::atomic<uint64_t> preflightInitialFlagCount{0};
    std::atomic<uint64_t> preflightAdditionalCount{0};
#pragma omp parallel for num_threads(preflightThreadCount) schedule(dynamic, 1)
    for (long long rawSample = 0; rawSample < static_cast<long long>(preflightSamples.size()); ++rawSample) {
        if (shouldCancel()) continue;
        const size_t index = preflightSamples[static_cast<size_t>(rawSample)];
        const PixelEvaluation primary = evaluatePixel(references[nearest[index]], index);
        const uint8_t riskFlags = primaryRiskFlags(primary.risk, primary.escapeIteration, primary.failed, index);
        const bool secondaryRequired = riskFlags != 0;
        bool mandatoryFallback = primary.failed || !finiteTrace(primary.trace);
        bool initiallyFlagged = false;
        double maximumStateDisagreement = 0.0;
        uint64_t additionalEvaluationCount = 0;
        const size_t x = index % static_cast<size_t>(request.width);
        const size_t y = index / static_cast<size_t>(request.width);
        const bool perimeter = x == 0 || x + 1 == static_cast<size_t>(request.width) || y == 0 || y + 1 == static_cast<size_t>(request.height);
        if (secondaryRequired) {
            const PixelEvaluation secondary = evaluatePixel(references[secondNearest[index]], index);
            maximumStateDisagreement = traceGap(primary.trace, secondary.trace);
            mandatoryFallback = mandatoryFallback || secondary.failed || !finiteTrace(secondary.trace) || primary.escapeIteration != secondary.escapeIteration;
            initiallyFlagged = perimeter || mandatoryFallback || !std::isfinite(maximumStateDisagreement) || maximumStateDisagreement >= StateDisagreementThreshold;
        }
        if (initiallyFlagged && !mandatoryFallback) {
            for (size_t referenceIndex = 0; referenceIndex < references.size(); ++referenceIndex) {
                if (referenceIndex == nearest[index] || referenceIndex == secondNearest[index]) continue;
                const PixelEvaluation additional = evaluatePixel(references[referenceIndex], index);
                ++additionalEvaluationCount;
                mandatoryFallback = mandatoryFallback || additional.failed || !finiteTrace(additional.trace) || additional.escapeIteration != primary.escapeIteration;
                const double stateDisagreement = traceGap(primary.trace, additional.trace);
                if (!std::isfinite(stateDisagreement)) {
                    maximumStateDisagreement = std::numeric_limits<double>::infinity();
                } else {
                    maximumStateDisagreement = std::max(maximumStateDisagreement, stateDisagreement);
                }
            }
        }
        const bool highStateDisagreement = secondaryRequired && (!std::isfinite(maximumStateDisagreement) || maximumStateDisagreement >= StateDisagreementThreshold);
        const bool doubleDoubleRequired = !mandatoryFallback && highStateDisagreement && primary.escapeIteration >= 0;
        preflightEvaluated[static_cast<size_t>(rawSample)] = {index, primary.output, primary.escapeIteration, !mandatoryFallback && !highStateDisagreement, doubleDoubleRequired};
        preflightSampleCount.fetch_add(1, std::memory_order_relaxed);
        preflightPrimaryRiskCount.fetch_add(secondaryRequired, std::memory_order_relaxed);
        preflightSecondaryCount.fetch_add(secondaryRequired, std::memory_order_relaxed);
        preflightInitialFlagCount.fetch_add(initiallyFlagged, std::memory_order_relaxed);
        preflightAdditionalCount.fetch_add(additionalEvaluationCount, std::memory_order_relaxed);
    }
    result.centeredPreflightSampleCount = preflightSampleCount.load(std::memory_order_relaxed);
    result.centeredPreflightPrimaryRiskFlagCount = preflightPrimaryRiskCount.load(std::memory_order_relaxed);
    result.centeredPreflightSecondaryEvaluationCount = preflightSecondaryCount.load(std::memory_order_relaxed);
    result.centeredPreflightInitialFlagCount = preflightInitialFlagCount.load(std::memory_order_relaxed);
    result.centeredPreflightAdditionalReferenceEvaluationCount = preflightAdditionalCount.load(std::memory_order_relaxed);
    if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;

    std::atomic<uint64_t> preflightDoubleDoubleVerified{0};
    std::atomic<uint64_t> preflightDoubleDoubleAgreed{0};
    std::atomic<uint64_t> preflightDoubleDoubleIterations{0};
    const Clock::time_point preflightDoubleDoubleStart = Clock::now();
#pragma omp parallel for num_threads(preflightThreadCount) schedule(dynamic, 1)
    for (long long rawIndex = 0; rawIndex < static_cast<long long>(preflightEvaluated.size()); ++rawIndex) {
        PreflightSampleResult& sample = preflightEvaluated[static_cast<size_t>(rawIndex)];
        if (!sample.doubleDoubleRequired || shouldCancel()) continue;
        const size_t x = sample.index % static_cast<size_t>(request.width);
        const size_t y = sample.index / static_cast<size_t>(request.width);
        ExpressionDoubleDoubleEvaluator::ComplexValue pixels[4]{};
        ExpressionDoubleDoubleEvaluator::OrbitResult verified[4];
        pixels[0] = {xCoordinates[x], yCoordinates[y]};
        ExpressionDoubleDoubleEvaluator::evaluateOrbit4(*request.runtimeProgram, request.fixed, request.pixelParameter, pixels, 1, request.maxIterations, request.bailout, shouldCancel, verified);
        const bool agreed = verified[0].defined && verified[0].escapeIteration == sample.centeredEscapeIteration;
        sample.accepted = agreed;
        preflightDoubleDoubleVerified.fetch_add(1, std::memory_order_relaxed);
        preflightDoubleDoubleAgreed.fetch_add(agreed, std::memory_order_relaxed);
        preflightDoubleDoubleIterations.fetch_add(verified[0].iterations, std::memory_order_relaxed);
    }
    const double preflightDoubleDoubleSeconds = secondsSince(preflightDoubleDoubleStart);
    const uint64_t preflightVerified = preflightDoubleDoubleVerified.load(std::memory_order_relaxed);
    const uint64_t preflightAgreed = preflightDoubleDoubleAgreed.load(std::memory_order_relaxed);
    result.centeredDoubleDoubleVerifiedPixelCount = saturatingAdd(result.centeredDoubleDoubleVerifiedPixelCount, preflightVerified);
    result.centeredDoubleDoubleAgreementCount = saturatingAdd(result.centeredDoubleDoubleAgreementCount, preflightAgreed);
    result.centeredDoubleDoubleRejectedCount = saturatingAdd(result.centeredDoubleDoubleRejectedCount, preflightVerified - preflightAgreed);
    result.centeredDoubleDoubleIterationCount = saturatingAdd(result.centeredDoubleDoubleIterationCount, preflightDoubleDoubleIterations.load(std::memory_order_relaxed));
    result.centeredDoubleDoubleSeconds += preflightDoubleDoubleSeconds;
    for (const PreflightSampleResult& sample : preflightEvaluated) preflightFlags += !sample.accepted;
    if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;

    std::atomic_bool preflightMismatch{false};
    std::atomic_bool preflightResourceError{false};
    std::atomic_bool preflightInternalError{false};
#pragma omp parallel num_threads(preflightThreadCount)
    {
        OracleWorkspaceRelease workerOracleWorkspaceRelease;
        std::unique_ptr<FallbackWorkerWorkspace> workspace;
        try {
            workspace = std::make_unique<FallbackWorkerWorkspace>(fallbackPrecision, request);
            if (!workspace->geometryReady) preflightInternalError.store(true, std::memory_order_relaxed);
        } catch (const std::bad_alloc&) {
            preflightResourceError.store(true, std::memory_order_relaxed);
        } catch (const std::length_error&) {
            preflightResourceError.store(true, std::memory_order_relaxed);
        } catch (...) {
            preflightInternalError.store(true, std::memory_order_relaxed);
        }
#pragma omp for schedule(dynamic, 1)
        for (long long rawIndex = 0; rawIndex < static_cast<long long>(preflightEvaluated.size()); ++rawIndex) {
            if (shouldCancel() || !workspace || preflightResourceError.load(std::memory_order_relaxed) || preflightInternalError.load(std::memory_order_relaxed)) continue;
            const PreflightSampleResult& sample = preflightEvaluated[static_cast<size_t>(rawIndex)];
            float exact = ExpressionDeepEmptyPixel;
            uint64_t iterations = 0;
            if (!exactPixel(*workspace, sample.index, iterations, exact) || (sample.accepted && sample.output != exact)) preflightMismatch.store(true, std::memory_order_relaxed);
        }
    }
    result.centeredPreflightFlagCount = preflightFlags;
    result.centeredPreflightSeconds = secondsSince(preflightStart);
    if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;
    if (preflightResourceError.load(std::memory_order_relaxed)) {
        error = "adaptive centered preflight worker allocation failed";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    if (preflightInternalError.load(std::memory_order_relaxed)) {
        error = "adaptive centered preflight worker workspace failed";
        return Centered4ExperimentalStatus::Failed;
    }
    const double predictedFallbackRate = result.centeredPreflightSampleCount == 0 ? 1.0 : static_cast<double>(preflightFlags) / result.centeredPreflightSampleCount;
    if (preflightMismatch.load(std::memory_order_relaxed) || (request.centered.requirePredictedBenefit && predictedFallbackRate > request.centered.maximumPredictedFallbackRate)) {
        result.centeredPreflightRejected = true;
        return Centered4ExperimentalStatus::NotApplicable;
    }
    std::vector<size_t>().swap(preflightSamples);
    std::vector<PreflightSampleResult>().swap(preflightEvaluated);

    size_t rendererBytes = 0;
    if (!checkedAddSize(rendererBytes, pixelCount, sizeof(size_t) * 2 + sizeof(CandidateTrace) + sizeof(PrimaryRiskTelemetry) + sizeof(double) + 5 * sizeof(uint8_t) + sizeof(int) * 2) || !checkedAddSize(rendererBytes, pixelCount, sizeof(size_t) * 5) || !checkedAddSize(rendererBytes, pixelCount, 2 * sizeof(uint8_t)) || !checkedAddSize(rendererBytes, pixelCount, sizeof(size_t) + sizeof(uint8_t)) || !checkedAddSize(rendererBytes, maximumFallbackValidationSamples, sizeof(FallbackValidationResult)) || !checkedAddSize(rendererBytes, maximumMandatoryValidationSamples, sizeof(FallbackValidationResult)) || !checkedAddSize(rendererBytes, 2 * references.size(), sizeof(std::vector<size_t>)) || !checkedAddSize(rendererBytes, static_cast<size_t>(request.width + request.height), sizeof(double) + sizeof(ExpressionDoubleDoubleEvaluator::Real))) {
        error = "adaptive centered workspace size overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    if (!checkedAddSize(rendererBytes, static_cast<size_t>(threadCount), fallbackThreadBytes) || !checkedAddSize(rendererBytes, static_cast<size_t>(threadCount), size_t{64} << 10)) {
        error = "adaptive centered thread workspace overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    if (request.memory.memoryLimitBytes != 0 && (outerLiveBytes > request.memory.memoryLimitBytes || referenceBytes > request.memory.memoryLimitBytes - outerLiveBytes || rendererBytes > request.memory.memoryLimitBytes - outerLiveBytes - referenceBytes)) return Centered4ExperimentalStatus::NotApplicable;
    result.referenceBytes += referenceBytes;
    result.rendererBytes += rendererBytes;

    std::vector<CandidateTrace> primaryTraces(pixelCount);
    std::vector<PrimaryRiskTelemetry> primaryRisk(pixelCount);
    std::vector<int> primaryEscapeIteration(pixelCount, -1);
    std::vector<int> secondaryEscapeIteration(pixelCount, -1);
    std::vector<double> stateDisagreement(pixelCount, std::numeric_limits<double>::quiet_NaN());
    std::vector<uint8_t> candidateFailure(pixelCount, 0);
    std::vector<uint8_t> candidateNonfinite(pixelCount, 0);
    std::vector<uint8_t> referenceOutputDisagreement(pixelCount, 0);
    std::vector<uint8_t> denominatorInstability(pixelCount, 0);
    struct CandidateLaneRecord {
        size_t outputIndex = 0;
        int iteration = 0;
        ExpressionContext referenceContext;
        ExpressionDeltaContext delta;
        Complex stateDelta{};
        CandidateTrace trace;
        PrimaryRiskTelemetry risk;
        float output = ExpressionDeepInteriorPixel;
        int escapeIteration = -1;
        bool failed = false;
        bool nonfinite = false;
        bool denominatorInstability = false;
        bool active = false;
    };
    enum class CandidateLaneLoadStatus : uint8_t {
        Exhausted,
        Ready,
        Complete
    };
    auto loadCandidateLane = [&](CandidateLaneRecord& lane, const AdaptiveReference& reference, const std::vector<size_t>& indices, std::atomic<size_t>& nextIndex) {
        const size_t position = nextIndex.fetch_add(1, std::memory_order_relaxed);
        if (position >= indices.size()) {
            lane.active = false;
            return CandidateLaneLoadStatus::Exhausted;
        }
        lane.outputIndex = indices[position];
        lane.iteration = 0;
        lane.delta = {};
        lane.stateDelta = {};
        lane.trace = {};
        lane.risk.currentError = PrimaryInitialError;
        lane.risk.maximumError = PrimaryInitialError;
        lane.risk.minimumBailoutMargin = std::numeric_limits<double>::infinity();
        lane.risk.escapeMargin = std::numeric_limits<double>::infinity();
        lane.risk.flags = 0;
        lane.output = ExpressionDeepInteriorPixel;
        lane.escapeIteration = -1;
        lane.failed = false;
        lane.nonfinite = false;
        lane.denominatorInstability = false;
        const int x = static_cast<int>(lane.outputIndex % static_cast<size_t>(request.width));
        const int y = static_cast<int>(lane.outputIndex / static_cast<size_t>(request.width));
        lane.referenceContext = request.fixed;
        lane.referenceContext.c = reference.c;
        lane.referenceContext.z0 = reference.z0;
        const Complex pixelOffset{xOffsets[static_cast<size_t>(x)] - xOffsets[static_cast<size_t>(reference.x)], yOffsets[static_cast<size_t>(y)] - yOffsets[static_cast<size_t>(reference.y)]};
        if (request.pixelParameter == FormulaParameter::C)
            lane.delta.c = pixelOffset;
        else {
            lane.delta.z0 = pixelOffset;
            lane.stateDelta = pixelOffset;
        }
        if (initialCandidateStateComplete(reference, lane.stateDelta, lane)) {
            lane.nonfinite = (lane.risk.flags & PrimaryRiskNonfinite) != 0;
            lane.active = false;
            return CandidateLaneLoadStatus::Complete;
        }
        lane.active = true;
        return CandidateLaneLoadStatus::Ready;
    };
    auto refillCandidateLane = [&](CandidateLaneRecord& lane, const AdaptiveReference& reference, const std::vector<size_t>& indices, std::atomic<size_t>& nextIndex, auto&& recordComplete) {
        for (;;) {
            const CandidateLaneLoadStatus status = loadCandidateLane(lane, reference, indices, nextIndex);
            if (status == CandidateLaneLoadStatus::Exhausted) return false;
            if (status == CandidateLaneLoadStatus::Ready) return true;
            recordComplete(lane);
        }
    };
    auto recordPrimaryLane = [&](const CandidateLaneRecord& lane) {
        request.output[lane.outputIndex] = lane.output;
        primaryEscapeIteration[lane.outputIndex] = lane.escapeIteration;
        primaryTraces[lane.outputIndex] = lane.trace;
        primaryRisk[lane.outputIndex] = lane.risk;
        if (lane.failed) candidateFailure[lane.outputIndex] = 1;
        if (lane.nonfinite || (lane.risk.flags & PrimaryRiskNonfinite) != 0) candidateNonfinite[lane.outputIndex] = 1;
        if (lane.denominatorInstability || (lane.risk.flags & PrimaryRiskDenominator) != 0) denominatorInstability[lane.outputIndex] = 1;
    };
    auto recordSecondaryLane = [&](const CandidateLaneRecord& lane) {
        if (lane.failed) candidateFailure[lane.outputIndex] = 1;
        if (lane.nonfinite) candidateNonfinite[lane.outputIndex] = 1;
        if (lane.denominatorInstability) denominatorInstability[lane.outputIndex] = 1;
        secondaryEscapeIteration[lane.outputIndex] = lane.escapeIteration;
        if (lane.escapeIteration != primaryEscapeIteration[lane.outputIndex]) referenceOutputDisagreement[lane.outputIndex] = 1;
        stateDisagreement[lane.outputIndex] = traceGap(primaryTraces[lane.outputIndex], lane.trace);
    };
    auto recordAdditionalLane = [&](const CandidateLaneRecord& lane) {
        const bool finiteCandidateTrace = finiteTrace(lane.trace);
        candidateFailure[lane.outputIndex] = candidateFailure[lane.outputIndex] != 0 || lane.failed || !finiteCandidateTrace;
        candidateNonfinite[lane.outputIndex] = candidateNonfinite[lane.outputIndex] != 0 || lane.nonfinite;
        denominatorInstability[lane.outputIndex] = denominatorInstability[lane.outputIndex] != 0 || lane.denominatorInstability;
        referenceOutputDisagreement[lane.outputIndex] = referenceOutputDisagreement[lane.outputIndex] != 0 || lane.escapeIteration != primaryEscapeIteration[lane.outputIndex];
        const double gap = traceGap(primaryTraces[lane.outputIndex], lane.trace);
        if (!std::isfinite(gap)) {
            stateDisagreement[lane.outputIndex] = std::numeric_limits<double>::infinity();
        } else {
            stateDisagreement[lane.outputIndex] = std::max(stateDisagreement[lane.outputIndex], gap);
        }
    };
    std::atomic<uint64_t> vmIterations{0};
    std::atomic<uint64_t> vectorSteps{0};
    std::atomic<uint64_t> vectorActiveLanes{0};
    progress(ExpressionDeepRenderPhase::Fast, 0, pixelCount);
    const Clock::time_point primaryStart = Clock::now();
    for (size_t referenceIndex = 0; referenceIndex < references.size(); ++referenceIndex) {
        const AdaptiveReference& reference = references[referenceIndex];
        const std::vector<size_t>& indices = primaryAssignments[referenceIndex];
        std::atomic<size_t> nextIndex{0};
#pragma omp parallel num_threads(threadCount)
        {
            std::array<CandidateLaneRecord, 4> lanes;
            int activeCount = 0;
            for (CandidateLaneRecord& lane : lanes) activeCount += refillCandidateLane(lane, reference, indices, nextIndex, recordPrimaryLane);
            uint64_t localIterations = 0;
            uint64_t localVectorSteps = 0;
            uint64_t localActiveLanes = 0;
            while (activeCount > 0) {
                uint8_t activeMask = 0;
                bool pollCancellation = false;
                std::array<const ExpressionContext*, 4> referenceContexts{};
                std::array<const Complex*, 4> nodeBasePointers{};
                std::array<const Complex*, 4> nodeAuxiliaryPointers{};
                std::array<const ExpressionDeltaContext*, 4> deltas{};
                for (int laneIndex = 0; laneIndex < 4; ++laneIndex) {
                    CandidateLaneRecord& lane = lanes[static_cast<size_t>(laneIndex)];
                    if (!lane.active) continue;
                    activeMask |= uint8_t{1} << laneIndex;
                    pollCancellation = pollCancellation || (lane.iteration & 63) == 0;
                    lane.referenceContext.iteration = lane.iteration;
                    lane.referenceContext.z = reference.z[static_cast<size_t>(lane.iteration)];
                    lane.delta.z = lane.stateDelta;
                    const ExpressionReferenceSample& sample = reference.orbit.samples[static_cast<size_t>(lane.iteration)];
                    referenceContexts[laneIndex] = &lane.referenceContext;
                    nodeBasePointers[laneIndex] = reference.nodes.data() + static_cast<size_t>(sample.tapeOffset);
                    nodeAuxiliaryPointers[laneIndex] = reference.nodeAuxiliaries.data() + static_cast<size_t>(sample.tapeOffset);
                    deltas[laneIndex] = &lane.delta;
                }
                if (pollCancellation && shouldCancel()) break;
                ExpressionCenteredResult evaluated[4];
                double jacobianNorms[4];
                const bool evaluatedBatch = ExpressionCenteredEvaluator::evaluate4WithLaneNodeBasesFastEntireDerivative(*request.runtimeProgram, referenceContexts.data(), nodeBasePointers.data(), nodeAuxiliaryPointers.data(), deltas.data(), activeMask, evaluated, jacobianNorms);
                ++localVectorSteps;
                localActiveLanes += static_cast<uint64_t>((activeMask & 1) != 0) + static_cast<uint64_t>((activeMask & 2) != 0) + static_cast<uint64_t>((activeMask & 4) != 0) + static_cast<uint64_t>((activeMask & 8) != 0);
                for (int laneIndex = 0; laneIndex < 4; ++laneIndex) {
                    CandidateLaneRecord& lane = lanes[static_cast<size_t>(laneIndex)];
                    if (!lane.active) continue;
                    if (evaluatedBatch) ++localIterations;
                    bool complete = false;
                    if (evaluatedBatch && evaluated[laneIndex].denominatorInstability) {
                        lane.denominatorInstability = true;
                        lane.risk.flags |= PrimaryRiskDenominator;
                    }
                    if (!evaluatedBatch || !evaluated[laneIndex].success()) {
                        lane.failed = true;
                        lane.risk.flags |= PrimaryRiskFailure;
                        complete = true;
                    } else {
                        const ExpressionReferenceSample& sample = reference.orbit.samples[static_cast<size_t>(lane.iteration)];
                        const Complex* laneNodeBases = nodeBasePointers[laneIndex];
                        const Complex nextReference = laneNodeBases[sample.rootNode];
                        const Complex actual = nextReference + evaluated[laneIndex].delta;
                        lane.stateDelta = lane.iteration + 1 < request.maxIterations ? evaluated[laneIndex].delta + (nextReference - reference.z[static_cast<size_t>(lane.iteration + 1)]) : evaluated[laneIndex].delta;
                        const double magnitude = std::abs(actual);
                        const double jacobianNorm = jacobianNorms[laneIndex];
                        if (std::isfinite(jacobianNorm)) {
                            lane.risk.currentError = jacobianNorm * lane.risk.currentError + PrimaryOperationErrorScale * std::max(1.0, magnitude);
                            lane.risk.maximumError = std::max(lane.risk.maximumError, lane.risk.currentError);
                        } else {
                            lane.risk.currentError = lane.risk.maximumError = std::numeric_limits<double>::infinity();
                            lane.risk.flags |= PrimaryRiskNonfinite;
                            lane.nonfinite = true;
                        }
                        lane.trace.states[4] = actual;
                        lane.trace.validMask |= uint8_t{1} << 4;
                        for (size_t probe = 0; probe < probeIterations.size(); ++probe) {
                            if (lane.iteration + 1 == probeIterations[probe]) {
                                lane.trace.states[probe] = actual;
                                lane.trace.validMask |= uint8_t{1} << probe;
                            }
                        }
                        if (std::isfinite(magnitude)) {
                            lane.risk.minimumBailoutMargin = std::min(lane.risk.minimumBailoutMargin, std::fabs(magnitude - request.bailout) / (1.0 + request.bailout));
                        } else {
                            lane.risk.flags |= PrimaryRiskNonfinite;
                        }
                        if (!finiteComplex(actual) || magnitude > request.bailout) {
                            lane.output = static_cast<float>(lane.iteration + 1);
                            lane.escapeIteration = lane.iteration + 1;
                            if (!finiteComplex(actual)) {
                                lane.failed = true;
                                lane.nonfinite = true;
                                lane.risk.flags |= PrimaryRiskFailure | PrimaryRiskNonfinite;
                            } else {
                                lane.risk.escapeMargin = (magnitude - request.bailout) / (1.0 + request.bailout);
                            }
                            complete = true;
                        } else if (lane.iteration + 1 >= request.maxIterations) {
                            complete = true;
                        }
                    }
                    if (!complete) {
                        ++lane.iteration;
                        continue;
                    }
                    recordPrimaryLane(lane);
                    lane.active = false;
                    --activeCount;
                    if (refillCandidateLane(lane, reference, indices, nextIndex, recordPrimaryLane)) ++activeCount;
                }
            }
            vmIterations.fetch_add(localIterations, std::memory_order_relaxed);
            vectorSteps.fetch_add(localVectorSteps, std::memory_order_relaxed);
            vectorActiveLanes.fetch_add(localActiveLanes, std::memory_order_relaxed);
        }
    }
    result.centeredPrimarySeconds = secondsSince(primaryStart);
    result.fastSeconds += result.centeredPrimarySeconds;
    if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;

    const Clock::time_point riskSelectorStart = Clock::now();
    std::vector<std::vector<size_t>> secondaryAssignments;
    try {
        secondaryAssignments.resize(references.size());
        for (size_t index = 0; index < pixelCount; ++index) {
            uint8_t flags = primaryRiskFlags(primaryRisk[index], primaryEscapeIteration[index], candidateFailure[index] != 0, index);
            primaryRisk[index].flags = flags;
            if (flags == 0) continue;
            secondaryAssignments[secondNearest[index]].push_back(index);
            ++result.centeredPrimaryRiskFlagCount;
        }
    } catch (const std::bad_alloc&) {
        error = "adaptive centered risk selector allocation failed";
        return Centered4ExperimentalStatus::ResourceLimit;
    } catch (const std::length_error&) {
        error = "adaptive centered risk selector length overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    result.centeredSecondaryEvaluationCount = result.centeredPrimaryRiskFlagCount;
    result.centeredSelectorSeconds = secondsSince(riskSelectorStart);

    const Clock::time_point secondaryStart = Clock::now();
    for (size_t referenceIndex = 0; referenceIndex < references.size(); ++referenceIndex) {
        const AdaptiveReference& reference = references[referenceIndex];
        const std::vector<size_t>& indices = secondaryAssignments[referenceIndex];
        std::atomic<size_t> nextIndex{0};
#pragma omp parallel num_threads(threadCount)
        {
            std::array<CandidateLaneRecord, 4> lanes;
            int activeCount = 0;
            for (CandidateLaneRecord& lane : lanes) activeCount += refillCandidateLane(lane, reference, indices, nextIndex, recordSecondaryLane);
            uint64_t localIterations = 0;
            uint64_t localVectorSteps = 0;
            uint64_t localActiveLanes = 0;
            while (activeCount > 0) {
                uint8_t activeMask = 0;
                bool pollCancellation = false;
                std::array<const ExpressionContext*, 4> referenceContexts{};
                std::array<const Complex*, 4> nodeBasePointers{};
                std::array<const Complex*, 4> nodeAuxiliaryPointers{};
                std::array<const ExpressionDeltaContext*, 4> deltas{};
                for (int laneIndex = 0; laneIndex < 4; ++laneIndex) {
                    CandidateLaneRecord& lane = lanes[static_cast<size_t>(laneIndex)];
                    if (!lane.active) continue;
                    activeMask |= uint8_t{1} << laneIndex;
                    pollCancellation = pollCancellation || (lane.iteration & 63) == 0;
                    lane.referenceContext.iteration = lane.iteration;
                    lane.referenceContext.z = reference.z[static_cast<size_t>(lane.iteration)];
                    lane.delta.z = lane.stateDelta;
                    const ExpressionReferenceSample& sample = reference.orbit.samples[static_cast<size_t>(lane.iteration)];
                    referenceContexts[laneIndex] = &lane.referenceContext;
                    nodeBasePointers[laneIndex] = reference.nodes.data() + static_cast<size_t>(sample.tapeOffset);
                    nodeAuxiliaryPointers[laneIndex] = reference.nodeAuxiliaries.data() + static_cast<size_t>(sample.tapeOffset);
                    deltas[laneIndex] = &lane.delta;
                }
                if (pollCancellation && shouldCancel()) break;
                ExpressionCenteredResult evaluated[4];
                const bool evaluatedBatch = ExpressionCenteredEvaluator::evaluate4WithLaneNodeBasesFastEntire(*request.runtimeProgram, referenceContexts.data(), nodeBasePointers.data(), nodeAuxiliaryPointers.data(), deltas.data(), activeMask, evaluated);
                ++localVectorSteps;
                localActiveLanes += static_cast<uint64_t>((activeMask & 1) != 0) + static_cast<uint64_t>((activeMask & 2) != 0) + static_cast<uint64_t>((activeMask & 4) != 0) + static_cast<uint64_t>((activeMask & 8) != 0);
                for (int laneIndex = 0; laneIndex < 4; ++laneIndex) {
                    CandidateLaneRecord& lane = lanes[static_cast<size_t>(laneIndex)];
                    if (!lane.active) continue;
                    if (evaluatedBatch) ++localIterations;
                    bool complete = false;
                    if (evaluatedBatch && evaluated[laneIndex].denominatorInstability) {
                        lane.denominatorInstability = true;
                        lane.risk.flags |= PrimaryRiskDenominator;
                    }
                    if (!evaluatedBatch || !evaluated[laneIndex].success()) {
                        lane.failed = true;
                        complete = true;
                    } else {
                        const ExpressionReferenceSample& sample = reference.orbit.samples[static_cast<size_t>(lane.iteration)];
                        const Complex* laneNodeBases = nodeBasePointers[laneIndex];
                        const Complex nextReference = laneNodeBases[sample.rootNode];
                        const Complex actual = nextReference + evaluated[laneIndex].delta;
                        lane.stateDelta = lane.iteration + 1 < request.maxIterations ? evaluated[laneIndex].delta + (nextReference - reference.z[static_cast<size_t>(lane.iteration + 1)]) : evaluated[laneIndex].delta;
                        lane.trace.states[4] = actual;
                        lane.trace.validMask |= uint8_t{1} << 4;
                        for (size_t probe = 0; probe < probeIterations.size(); ++probe) {
                            if (lane.iteration + 1 == probeIterations[probe]) {
                                lane.trace.states[probe] = actual;
                                lane.trace.validMask |= uint8_t{1} << probe;
                            }
                        }
                        if (!finiteComplex(actual) || std::hypot(actual.real(), actual.imag()) > request.bailout) {
                            lane.escapeIteration = lane.iteration + 1;
                            lane.failed = !finiteComplex(actual);
                            lane.nonfinite = lane.failed;
                            complete = true;
                        } else if (lane.iteration + 1 >= request.maxIterations) {
                            complete = true;
                        }
                    }
                    if (!complete) {
                        ++lane.iteration;
                        continue;
                    }
                    recordSecondaryLane(lane);
                    lane.active = false;
                    --activeCount;
                    if (refillCandidateLane(lane, reference, indices, nextIndex, recordSecondaryLane)) ++activeCount;
                }
            }
            vmIterations.fetch_add(localIterations, std::memory_order_relaxed);
            vectorSteps.fetch_add(localVectorSteps, std::memory_order_relaxed);
            vectorActiveLanes.fetch_add(localActiveLanes, std::memory_order_relaxed);
        }
    }
    result.centeredSecondarySeconds = secondsSince(secondaryStart);
    result.centeredCandidateSeconds = result.centeredPrimarySeconds + result.centeredSecondarySeconds;
    result.fastSeconds += result.centeredSecondarySeconds;
    if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;

    const Clock::time_point disagreementSelectorStart = Clock::now();
    std::vector<size_t> fallback;
    try {
        fallback.reserve(static_cast<size_t>(result.centeredSecondaryEvaluationCount));
        for (size_t index = 0; index < pixelCount; ++index) {
            if (primaryRisk[index].flags == 0) continue;
            const size_t x = index % static_cast<size_t>(request.width);
            const size_t y = index / static_cast<size_t>(request.width);
            // A perimeter pixel must first be checked against every retained
            // reference before the approximate candidate can be accepted.
            const bool perimeter = x == 0 || x + 1 == static_cast<size_t>(request.width) || y == 0 || y + 1 == static_cast<size_t>(request.height);
            const bool outputDisagreement = referenceOutputDisagreement[index] != 0;
            const bool interiorConservatism = primaryEscapeIteration[index] < 0;
            const bool flagged = perimeter || interiorConservatism || candidateFailure[index] != 0 || outputDisagreement || !std::isfinite(stateDisagreement[index]) || stateDisagreement[index] >= StateDisagreementThreshold;
            if (!flagged) continue;
            fallback.push_back(index);
        }
    } catch (const std::bad_alloc&) {
        error = "adaptive centered selector allocation failed";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    result.centeredSelectorSeconds += secondsSince(disagreementSelectorStart);
    result.centeredSelectorFlagCount = fallback.size();
    const size_t initialFlagCount = fallback.size();

    const Clock::time_point additionalReferenceStart = Clock::now();
    std::vector<size_t> additionalIndices;
    try {
        additionalIndices.reserve(initialFlagCount);
        for (size_t referenceIndex = 0; referenceIndex < references.size(); ++referenceIndex) {
            additionalIndices.clear();
            for (size_t outputIndex : fallback)
                if (candidateFailure[outputIndex] == 0 && referenceIndex != nearest[outputIndex] && referenceIndex != secondNearest[outputIndex]) additionalIndices.push_back(outputIndex);
            result.centeredAdditionalReferenceEvaluationCount = saturatingAdd(result.centeredAdditionalReferenceEvaluationCount, additionalIndices.size());
            result.centeredHierarchyEvaluationCount = saturatingAdd(result.centeredHierarchyEvaluationCount, additionalIndices.size());
            if (additionalIndices.empty()) continue;
            const AdaptiveReference& reference = references[referenceIndex];
            std::atomic<size_t> nextIndex{0};
#pragma omp parallel num_threads(threadCount)
            {
                std::array<CandidateLaneRecord, 4> lanes;
                int activeCount = 0;
                for (CandidateLaneRecord& lane : lanes) activeCount += refillCandidateLane(lane, reference, additionalIndices, nextIndex, recordAdditionalLane);
                uint64_t localIterations = 0;
                uint64_t localVectorSteps = 0;
                uint64_t localActiveLanes = 0;
                while (activeCount > 0) {
                    uint8_t activeMask = 0;
                    bool pollCancellation = false;
                    std::array<const ExpressionContext*, 4> referenceContexts{};
                    std::array<const Complex*, 4> nodeBasePointers{};
                    std::array<const Complex*, 4> nodeAuxiliaryPointers{};
                    std::array<const ExpressionDeltaContext*, 4> deltas{};
                    for (int laneIndex = 0; laneIndex < 4; ++laneIndex) {
                        CandidateLaneRecord& lane = lanes[static_cast<size_t>(laneIndex)];
                        if (!lane.active) continue;
                        activeMask |= uint8_t{1} << laneIndex;
                        pollCancellation = pollCancellation || (lane.iteration & 63) == 0;
                        lane.referenceContext.iteration = lane.iteration;
                        lane.referenceContext.z = reference.z[static_cast<size_t>(lane.iteration)];
                        lane.delta.z = lane.stateDelta;
                        const ExpressionReferenceSample& sample = reference.orbit.samples[static_cast<size_t>(lane.iteration)];
                        referenceContexts[laneIndex] = &lane.referenceContext;
                        nodeBasePointers[laneIndex] = reference.nodes.data() + static_cast<size_t>(sample.tapeOffset);
                        nodeAuxiliaryPointers[laneIndex] = reference.nodeAuxiliaries.data() + static_cast<size_t>(sample.tapeOffset);
                        deltas[laneIndex] = &lane.delta;
                    }
                    if (pollCancellation && shouldCancel()) break;
                    ExpressionCenteredResult evaluated[4];
                    const bool evaluatedBatch = ExpressionCenteredEvaluator::evaluate4WithLaneNodeBasesFastEntire(*request.runtimeProgram, referenceContexts.data(), nodeBasePointers.data(), nodeAuxiliaryPointers.data(), deltas.data(), activeMask, evaluated);
                    ++localVectorSteps;
                    localActiveLanes += static_cast<uint64_t>((activeMask & 1) != 0) + static_cast<uint64_t>((activeMask & 2) != 0) + static_cast<uint64_t>((activeMask & 4) != 0) + static_cast<uint64_t>((activeMask & 8) != 0);
                    for (int laneIndex = 0; laneIndex < 4; ++laneIndex) {
                        CandidateLaneRecord& lane = lanes[static_cast<size_t>(laneIndex)];
                        if (!lane.active) continue;
                        if (evaluatedBatch) ++localIterations;
                        bool complete = false;
                        if (evaluatedBatch && evaluated[laneIndex].denominatorInstability) {
                            lane.denominatorInstability = true;
                            lane.risk.flags |= PrimaryRiskDenominator;
                        }
                        if (!evaluatedBatch || !evaluated[laneIndex].success()) {
                            lane.failed = true;
                            complete = true;
                        } else {
                            const ExpressionReferenceSample& sample = reference.orbit.samples[static_cast<size_t>(lane.iteration)];
                            const Complex* laneNodeBases = nodeBasePointers[laneIndex];
                            const Complex nextReference = laneNodeBases[sample.rootNode];
                            const Complex actual = nextReference + evaluated[laneIndex].delta;
                            lane.stateDelta = lane.iteration + 1 < request.maxIterations ? evaluated[laneIndex].delta + (nextReference - reference.z[static_cast<size_t>(lane.iteration + 1)]) : evaluated[laneIndex].delta;
                            lane.trace.states[4] = actual;
                            lane.trace.validMask |= uint8_t{1} << 4;
                            for (size_t probe = 0; probe < probeIterations.size(); ++probe) {
                                if (lane.iteration + 1 == probeIterations[probe]) {
                                    lane.trace.states[probe] = actual;
                                    lane.trace.validMask |= uint8_t{1} << probe;
                                }
                            }
                            if (!finiteComplex(actual) || std::hypot(actual.real(), actual.imag()) > request.bailout) {
                                lane.escapeIteration = lane.iteration + 1;
                                lane.failed = !finiteComplex(actual);
                                lane.nonfinite = lane.failed;
                                complete = true;
                            } else if (lane.iteration + 1 >= request.maxIterations) {
                                complete = true;
                            }
                        }
                        if (!complete) {
                            ++lane.iteration;
                            continue;
                        }
                        recordAdditionalLane(lane);
                        lane.active = false;
                        --activeCount;
                        if (refillCandidateLane(lane, reference, additionalIndices, nextIndex, recordAdditionalLane)) ++activeCount;
                    }
                }
                vmIterations.fetch_add(localIterations, std::memory_order_relaxed);
                vectorSteps.fetch_add(localVectorSteps, std::memory_order_relaxed);
                vectorActiveLanes.fetch_add(localActiveLanes, std::memory_order_relaxed);
            }
        }
    } catch (const std::bad_alloc&) {
        error = "adaptive centered additional-reference allocation failed";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    result.centeredAdditionalReferenceSeconds = secondsSince(additionalReferenceStart);
    result.fastSeconds += result.centeredAdditionalReferenceSeconds;
    if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;

    const Clock::time_point finalSelectorStart = Clock::now();
    std::vector<size_t> doubleDoubleCandidates;
    std::vector<uint8_t> doubleDoubleAgreement;
    try {
        doubleDoubleCandidates.reserve(fallback.size());
        for (size_t outputIndex : fallback) {
            const bool finiteHighStateDisagreement = std::isfinite(stateDisagreement[outputIndex]) && stateDisagreement[outputIndex] >= StateDisagreementThreshold;
            if (candidateFailure[outputIndex] == 0 && candidateNonfinite[outputIndex] == 0 && denominatorInstability[outputIndex] == 0 && referenceOutputDisagreement[outputIndex] == 0 && finiteHighStateDisagreement && primaryEscapeIteration[outputIndex] >= 0) doubleDoubleCandidates.push_back(outputIndex);
        }
        doubleDoubleAgreement.assign(doubleDoubleCandidates.size(), 0);
    } catch (const std::bad_alloc&) {
        error = "adaptive centered double-double selector allocation failed";
        return Centered4ExperimentalStatus::ResourceLimit;
    } catch (const std::length_error&) {
        error = "adaptive centered double-double selector length overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }

    std::atomic<uint64_t> doubleDoubleIterations{0};
    std::atomic<uint64_t> doubleDoubleAgreed{0};
    const Clock::time_point doubleDoubleStart = Clock::now();
#pragma omp parallel for num_threads(threadCount) schedule(dynamic, 2)
    for (long long rawGroup = 0; rawGroup < static_cast<long long>((doubleDoubleCandidates.size() + 3) / 4); ++rawGroup) {
        if (shouldCancel()) continue;
        const size_t groupBegin = static_cast<size_t>(rawGroup) * 4;
        ExpressionDoubleDoubleEvaluator::ComplexValue pixels[4]{};
        ExpressionDoubleDoubleEvaluator::OrbitResult verified[4];
        uint8_t activeMask = 0;
        for (int lane = 0; lane < 4 && groupBegin + static_cast<size_t>(lane) < doubleDoubleCandidates.size(); ++lane) {
            const size_t outputIndex = doubleDoubleCandidates[groupBegin + static_cast<size_t>(lane)];
            const size_t x = outputIndex % static_cast<size_t>(request.width);
            const size_t y = outputIndex / static_cast<size_t>(request.width);
            pixels[lane] = {xCoordinates[x], yCoordinates[y]};
            activeMask |= uint8_t{1} << lane;
        }
        ExpressionDoubleDoubleEvaluator::evaluateOrbit4(*request.runtimeProgram, request.fixed, request.pixelParameter, pixels, activeMask, request.maxIterations, request.bailout, shouldCancel, verified);
        uint64_t localIterations = 0;
        uint64_t localAgreed = 0;
        for (int lane = 0; lane < 4 && groupBegin + static_cast<size_t>(lane) < doubleDoubleCandidates.size(); ++lane) {
            const size_t rawIndex = groupBegin + static_cast<size_t>(lane);
            const size_t outputIndex = doubleDoubleCandidates[rawIndex];
            const bool agreed = verified[lane].defined && verified[lane].escapeIteration == primaryEscapeIteration[outputIndex];
            doubleDoubleAgreement[rawIndex] = agreed;
            localIterations += verified[lane].iterations;
            localAgreed += agreed;
        }
        doubleDoubleIterations.fetch_add(localIterations, std::memory_order_relaxed);
        doubleDoubleAgreed.fetch_add(localAgreed, std::memory_order_relaxed);
    }
    const double doubleDoubleSeconds = secondsSince(doubleDoubleStart);
    if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;
    const uint64_t verifiedCount = doubleDoubleCandidates.size();
    const uint64_t agreedCount = doubleDoubleAgreed.load(std::memory_order_relaxed);
    result.centeredDoubleDoubleVerifiedPixelCount = saturatingAdd(result.centeredDoubleDoubleVerifiedPixelCount, verifiedCount);
    result.centeredDoubleDoubleAgreementCount = saturatingAdd(result.centeredDoubleDoubleAgreementCount, agreedCount);
    result.centeredDoubleDoubleRejectedCount = saturatingAdd(result.centeredDoubleDoubleRejectedCount, verifiedCount - agreedCount);
    result.centeredDoubleDoubleIterationCount = saturatingAdd(result.centeredDoubleDoubleIterationCount, doubleDoubleIterations.load(std::memory_order_relaxed));
    result.centeredDoubleDoubleSeconds += doubleDoubleSeconds;
    result.fastSeconds += doubleDoubleSeconds;

    const auto centeredFallbackReasonBit = [](ExpressionDeepCenteredFallbackReason reason) { return expressionDeepCenteredFallbackReasonMask(reason); };
    size_t finalFallbackCount = 0;
    size_t doubleDoubleIndex = 0;
    std::vector<size_t> mandatoryFullFallback;
    std::vector<size_t> lowEligibleFallback;
    std::vector<uint8_t> fallbackMask;
    try {
        mandatoryFullFallback.reserve(fallback.size());
        lowEligibleFallback.reserve(fallback.size());
        fallbackMask.assign(pixelCount, 0);
        result.centeredFallbackPrecisionReasonMask.assign(pixelCount, 0);
        for (size_t outputIndex : fallback) {
            uint8_t reasonMask = 0;
            if (candidateFailure[outputIndex] != 0) reasonMask |= centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::CandidateFailure);
            if (candidateNonfinite[outputIndex] != 0) reasonMask |= centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::Nonfinite);
            if (denominatorInstability[outputIndex] != 0) reasonMask |= centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::DenominatorInstability);
            if (referenceOutputDisagreement[outputIndex] != 0) reasonMask |= centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::ReferenceOutputDisagreement);

            const bool finiteStateDisagreement = std::isfinite(stateDisagreement[outputIndex]);
            const bool highStateDisagreement = !finiteStateDisagreement || stateDisagreement[outputIndex] >= StateDisagreementThreshold;
            if (highStateDisagreement) reasonMask |= centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::StateDisagreement);
            if (!finiteStateDisagreement) reasonMask |= centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::StructuralInconsistency);

            bool finalFallback = (reasonMask & ExpressionDeepCenteredFallbackMandatoryMask) != 0;
            const bool doubleDoubleCandidate = candidateFailure[outputIndex] == 0 && candidateNonfinite[outputIndex] == 0 && denominatorInstability[outputIndex] == 0 && referenceOutputDisagreement[outputIndex] == 0 && finiteStateDisagreement && highStateDisagreement && primaryEscapeIteration[outputIndex] >= 0;
            if (!finalFallback && highStateDisagreement) {
                if (primaryEscapeIteration[outputIndex] < 0) {
                    reasonMask |= centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::InteriorConservatism);
                    finalFallback = true;
                } else if (!doubleDoubleCandidate || doubleDoubleIndex >= doubleDoubleAgreement.size()) {
                    reasonMask |= centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::StructuralInconsistency);
                    finalFallback = true;
                } else {
                    if (doubleDoubleAgreement[doubleDoubleIndex] == 0) {
                        reasonMask |= centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::DoubleDoubleRejected);
                        finalFallback = true;
                    }
                    ++doubleDoubleIndex;
                }
            } else if (doubleDoubleCandidate) {
                ++doubleDoubleIndex;
            }
            if (!finalFallback) continue;

            if ((reasonMask & ExpressionDeepCenteredFallbackMandatoryMask) == 0) {
                const uint8_t lowEligibleReasons = centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::StateDisagreement) | centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::InteriorConservatism);
                if (reasonMask != lowEligibleReasons) reasonMask |= centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::StructuralInconsistency);
            }
            fallback[finalFallbackCount++] = outputIndex;
            result.centeredFallbackPrecisionReasonMask[outputIndex] = reasonMask;
            for (size_t reason = 0; reason < static_cast<size_t>(ExpressionDeepCenteredFallbackReason::Count); ++reason)
                if ((reasonMask & (uint8_t{1} << reason)) != 0) ++result.centeredFallbackPrecisionReasonCounts[reason];
            if ((reasonMask & ExpressionDeepCenteredFallbackMandatoryMask) != 0) {
                fallbackMask[outputIndex] = 1;
                mandatoryFullFallback.push_back(outputIndex);
            } else {
                fallbackMask[outputIndex] = 2;
                lowEligibleFallback.push_back(outputIndex);
            }
        }
        fallback.resize(finalFallbackCount);
        finalFallbackCount = fallback.size();
    } catch (const std::bad_alloc&) {
        std::fill(request.output, request.output + pixelCount, ExpressionDeepEmptyPixel);
        error = "adaptive centered fallback classification allocation failed";
        return Centered4ExperimentalStatus::ResourceLimit;
    } catch (const std::length_error&) {
        std::fill(request.output, request.output + pixelCount, ExpressionDeepEmptyPixel);
        error = "adaptive centered fallback classification length overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    fallback.resize(finalFallbackCount);
    result.centeredFinalSelectorSeconds = secondsSince(finalSelectorStart);
    result.centeredFinalFallbackFlagCount = fallback.size();
    result.centeredFallbackMandatoryFullPixelCount = mandatoryFullFallback.size();
    result.centeredFallbackLowEligiblePixelCount = lowEligibleFallback.size();
    result.fastPixelCount = static_cast<uint64_t>(pixelCount - fallback.size());
    result.fallbackPixelCount = static_cast<uint64_t>(fallback.size());
    result.fallbackReasonCounts.fill(0);
    result.fallbackReasonCounts[static_cast<size_t>(ExpressionDeepFallbackReason::CertificationFailure)] = result.fallbackPixelCount;
    result.maxFallbackReasonCount = result.fallbackPixelCount;
    result.uncertainPixelCount = result.fallbackPixelCount;
    result.undefinedPixelCount = 0;
    result.fallbackTileCount = 0;
    result.maxTileFallbackRate = 0.0;
    const int tilesX = (request.width + request.threading.tileWidth - 1) / request.threading.tileWidth;
    const int tilesY = (request.height + request.threading.tileHeight - 1) / request.threading.tileHeight;
    for (int tileY = 0; tileY < tilesY; ++tileY) {
        for (int tileX = 0; tileX < tilesX; ++tileX) {
            const int xBegin = tileX * request.threading.tileWidth;
            const int yBegin = tileY * request.threading.tileHeight;
            const int xEnd = std::min(request.width, xBegin + request.threading.tileWidth);
            const int yEnd = std::min(request.height, yBegin + request.threading.tileHeight);
            uint64_t fallbackInTile = 0;
            for (int y = yBegin; y < yEnd; ++y)
                for (int x = xBegin; x < xEnd; ++x) fallbackInTile += fallbackMask[static_cast<size_t>(y) * request.width + x] != 0;
            if (fallbackInTile == 0) continue;
            ++result.fallbackTileCount;
            const uint64_t pixelsInTile = static_cast<uint64_t>(xEnd - xBegin) * static_cast<uint64_t>(yEnd - yBegin);
            result.maxTileFallbackRate = std::max(result.maxTileFallbackRate, static_cast<double>(fallbackInTile) / pixelsInTile);
        }
    }

    progress(ExpressionDeepRenderPhase::Fallback, 0, fallback.size());
    const Clock::time_point fallbackStart = Clock::now();
    std::atomic_bool resourceError{false};
    std::atomic_bool internalError{false};
    std::vector<FallbackValidationResult> validation;
    std::vector<FallbackValidationResult> mandatoryValidation;
    if (maximumMandatoryValidationSamples != 0 && !mandatoryFullFallback.empty()) {
        try {
            mandatoryValidation.reserve(std::min(maximumMandatoryValidationSamples, mandatoryFullFallback.size()));
            auto addValidation = [&](size_t index) {
                if (mandatoryValidation.size() >= maximumMandatoryValidationSamples || fallbackMask[index] != 1) return;
                fallbackMask[index] = 4;
                mandatoryValidation.push_back({index});
            };
            auto riskFlagCount = [](uint8_t flags) {
                int count = 0;
                for (; flags != 0; flags &= static_cast<uint8_t>(flags - 1)) ++count;
                return count;
            };
            auto moreRisky = [&](size_t left, size_t right) {
                const bool leftFailed = candidateFailure[left] != 0;
                const bool rightFailed = candidateFailure[right] != 0;
                if (leftFailed != rightFailed) return leftFailed;
                const bool leftNonfinite = !std::isfinite(stateDisagreement[left]);
                const bool rightNonfinite = !std::isfinite(stateDisagreement[right]);
                if (leftNonfinite != rightNonfinite) return leftNonfinite;
                const int leftFlags = riskFlagCount(primaryRisk[left].flags);
                const int rightFlags = riskFlagCount(primaryRisk[right].flags);
                if (leftFlags != rightFlags) return leftFlags > rightFlags;
                if (stateDisagreement[left] != stateDisagreement[right]) return stateDisagreement[left] > stateDisagreement[right];
                return left < right;
            };
            auto addEvenly = [&](const std::vector<size_t>& indices, size_t count) {
                count = std::min(count, indices.size());
                for (size_t sample = 0; sample < count && mandatoryValidation.size() < maximumMandatoryValidationSamples; ++sample) {
                    const size_t position = count == 1 ? 0 : sample * (indices.size() - 1) / (count - 1);
                    addValidation(indices[position]);
                }
            };
            const size_t categorySamples = std::max<size_t>(1, maximumMandatoryValidationSamples / 5);
            const uint8_t structuralMask = centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::CandidateFailure) | centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::Nonfinite) | centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::StructuralInconsistency) | centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::DenominatorInstability);
            for (uint8_t reasonMask : {structuralMask, centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::ReferenceOutputDisagreement), centeredFallbackReasonBit(ExpressionDeepCenteredFallbackReason::DoubleDoubleRejected)}) {
                additionalIndices.clear();
                for (size_t index : mandatoryFullFallback)
                    if ((result.centeredFallbackPrecisionReasonMask[index] & reasonMask) != 0) additionalIndices.push_back(index);
                std::sort(additionalIndices.begin(), additionalIndices.end(), moreRisky);
                addEvenly(additionalIndices, categorySamples);
            }

            additionalIndices.clear();
            for (size_t index : mandatoryFullFallback) {
                const size_t x = index % static_cast<size_t>(request.width);
                const size_t y = index / static_cast<size_t>(request.width);
                if (x == 0 || x + 1 == static_cast<size_t>(request.width) || y == 0 || y + 1 == static_cast<size_t>(request.height)) additionalIndices.push_back(index);
            }
            std::sort(additionalIndices.begin(), additionalIndices.end(), moreRisky);
            addEvenly(additionalIndices, categorySamples);

            additionalIndices = mandatoryFullFallback;
            std::sort(additionalIndices.begin(), additionalIndices.end(), moreRisky);
            addEvenly(additionalIndices, maximumMandatoryValidationSamples - mandatoryValidation.size());
            addEvenly(mandatoryFullFallback, maximumMandatoryValidationSamples - mandatoryValidation.size());
            for (size_t index : mandatoryFullFallback) addValidation(index);
        } catch (const std::bad_alloc&) {
            error = "adaptive centered mandatory fallback validation allocation failed";
            return Centered4ExperimentalStatus::ResourceLimit;
        } catch (const std::length_error&) {
            error = "adaptive centered mandatory fallback validation length overflow";
            return Centered4ExperimentalStatus::ResourceLimit;
        }
    }
    result.centeredFallbackMandatoryValidationSampleCount = mandatoryValidation.size();
    if (maximumFallbackValidationSamples != 0 && !lowEligibleFallback.empty()) {
        try {
            validation.reserve(std::min(maximumFallbackValidationSamples, lowEligibleFallback.size()));
            auto addValidation = [&](size_t index) {
                if (validation.size() >= maximumFallbackValidationSamples || fallbackMask[index] != 2) return;
                fallbackMask[index] = 3;
                validation.push_back({index});
            };
            auto riskFlagCount = [](uint8_t flags) {
                int count = 0;
                for (; flags != 0; flags &= static_cast<uint8_t>(flags - 1)) ++count;
                return count;
            };
            auto moreRisky = [&](size_t left, size_t right) {
                const bool leftFailed = candidateFailure[left] != 0;
                const bool rightFailed = candidateFailure[right] != 0;
                if (leftFailed != rightFailed) return leftFailed;
                const bool leftNonfinite = !std::isfinite(stateDisagreement[left]);
                const bool rightNonfinite = !std::isfinite(stateDisagreement[right]);
                if (leftNonfinite != rightNonfinite) return leftNonfinite;
                const int leftFlags = riskFlagCount(primaryRisk[left].flags);
                const int rightFlags = riskFlagCount(primaryRisk[right].flags);
                if (leftFlags != rightFlags) return leftFlags > rightFlags;
                if (stateDisagreement[left] != stateDisagreement[right]) return stateDisagreement[left] > stateDisagreement[right];
                return left < right;
            };
            auto addEvenly = [&](const std::vector<size_t>& indices, size_t count) {
                count = std::min(count, indices.size());
                for (size_t sample = 0; sample < count && validation.size() < maximumFallbackValidationSamples; ++sample) {
                    const size_t position = count == 1 ? 0 : sample * (indices.size() - 1) / (count - 1);
                    addValidation(indices[position]);
                }
            };
            const size_t categorySamples = std::max<size_t>(1, maximumFallbackValidationSamples / 4);

            additionalIndices.clear();
            for (size_t index : lowEligibleFallback) {
                const size_t x = index % static_cast<size_t>(request.width);
                const size_t y = index / static_cast<size_t>(request.width);
                if (x == 0 || x + 1 == static_cast<size_t>(request.width) || y == 0 || y + 1 == static_cast<size_t>(request.height)) additionalIndices.push_back(index);
            }
            std::sort(additionalIndices.begin(), additionalIndices.end(), moreRisky);
            addEvenly(additionalIndices, categorySamples);

            additionalIndices.clear();
            for (size_t index : lowEligibleFallback)
                if (primaryEscapeIteration[index] < 0) additionalIndices.push_back(index);
            std::sort(additionalIndices.begin(), additionalIndices.end(), moreRisky);
            addEvenly(additionalIndices, categorySamples);

            additionalIndices = lowEligibleFallback;
            std::sort(additionalIndices.begin(), additionalIndices.end(), moreRisky);
            addEvenly(additionalIndices, categorySamples);
            for (uint8_t flag = PrimaryRiskFailure; flag <= PrimaryRiskDenominator && validation.size() < maximumFallbackValidationSamples; flag <<= 1) {
                const auto found = std::find_if(additionalIndices.begin(), additionalIndices.end(), [&](size_t index) { return (primaryRisk[index].flags & flag) != 0; });
                if (found != additionalIndices.end()) addValidation(*found);
            }

            addEvenly(lowEligibleFallback, maximumFallbackValidationSamples - validation.size());
            for (size_t index : lowEligibleFallback) addValidation(index);
        } catch (const std::bad_alloc&) {
            error = "adaptive centered fallback validation allocation failed";
            return Centered4ExperimentalStatus::ResourceLimit;
        } catch (const std::length_error&) {
            error = "adaptive centered fallback validation length overflow";
            return Centered4ExperimentalStatus::ResourceLimit;
        }
    }
    result.centeredFallbackValidationSampleCount = validation.size();

    auto evaluateValidation = [&](std::vector<FallbackValidationResult>& samples, mpfr_prec_t precision, bool candidatePrecision) {
        if (samples.empty()) return;
#pragma omp parallel num_threads(threadCount)
        {
            OracleWorkspaceRelease workerOracleWorkspaceRelease;
            std::unique_ptr<FallbackWorkerWorkspace> workspace;
            try {
                workspace = std::make_unique<FallbackWorkerWorkspace>(precision, request);
                if (!workspace->geometryReady && !candidatePrecision) internalError.store(true, std::memory_order_relaxed);
            } catch (const std::bad_alloc&) {
                resourceError.store(true, std::memory_order_relaxed);
            } catch (const std::length_error&) {
                resourceError.store(true, std::memory_order_relaxed);
            } catch (...) {
                internalError.store(true, std::memory_order_relaxed);
            }
#pragma omp for schedule(dynamic, 1)
            for (long long rawIndex = 0; rawIndex < static_cast<long long>(samples.size()); ++rawIndex) {
                try {
                    if (shouldCancel() || !workspace || resourceError.load(std::memory_order_relaxed) || internalError.load(std::memory_order_relaxed)) continue;
                    FallbackValidationResult& sample = samples[static_cast<size_t>(rawIndex)];
                    float output = ExpressionDeepEmptyPixel;
                    uint64_t iterations = 0;
                    const bool defined = workspace->geometryReady && exactPixel(*workspace, sample.index, iterations, output);
                    if (candidatePrecision) {
                        sample.lowDefined = defined;
                        sample.lowOutput = output;
                        sample.lowIterations = iterations;
                    } else {
                        sample.fullDefined = defined;
                        sample.fullOutput = output;
                        sample.fullIterations = iterations;
                    }
                } catch (const std::bad_alloc&) {
                    resourceError.store(true, std::memory_order_relaxed);
                } catch (const std::length_error&) {
                    resourceError.store(true, std::memory_order_relaxed);
                } catch (...) {
                    internalError.store(true, std::memory_order_relaxed);
                }
            }
        }
    };

    if (!mandatoryValidation.empty()) {
        const Clock::time_point candidateValidationStart = Clock::now();
        evaluateValidation(mandatoryValidation, mandatoryCandidatePrecision, true);
        result.centeredFallbackMandatoryCandidateSeconds += secondsSince(candidateValidationStart);
        if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;
        if (resourceError.load(std::memory_order_relaxed)) {
            error = "adaptive centered mandatory candidate validation allocation failed";
            return Centered4ExperimentalStatus::ResourceLimit;
        }
        if (internalError.load(std::memory_order_relaxed)) {
            error = "adaptive centered mandatory candidate validation failed";
            return Centered4ExperimentalStatus::Failed;
        }
        const Clock::time_point fullValidationStart = Clock::now();
        evaluateValidation(mandatoryValidation, fallbackPrecision, false);
        result.centeredFallbackMandatoryFullPrecisionSeconds += secondsSince(fullValidationStart);
        if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;
        if (resourceError.load(std::memory_order_relaxed)) {
            error = "adaptive centered mandatory full validation allocation failed";
            return Centered4ExperimentalStatus::ResourceLimit;
        }
        if (internalError.load(std::memory_order_relaxed)) {
            error = "adaptive centered mandatory full validation failed";
            return Centered4ExperimentalStatus::Failed;
        }
        result.centeredFallbackMandatoryCandidatePixelCount = mandatoryValidation.size();
        result.centeredFallbackMandatoryFullPrecisionPixelCount = mandatoryValidation.size();
        for (const FallbackValidationResult& sample : mandatoryValidation) {
            result.centeredFallbackMandatoryCandidateIterationCount = saturatingAdd(result.centeredFallbackMandatoryCandidateIterationCount, sample.lowIterations);
            result.centeredFallbackMandatoryFullPrecisionIterationCount = saturatingAdd(result.centeredFallbackMandatoryFullPrecisionIterationCount, sample.fullIterations);
            if (!sample.lowDefined || !sample.fullDefined || sample.lowOutput != sample.fullOutput) ++result.centeredFallbackMandatoryValidationMismatchCount;
        }
        result.centeredFallbackMandatoryUpgraded = result.centeredFallbackMandatoryValidationMismatchCount != 0;
    }

    if (!validation.empty()) {
        const Clock::time_point lowValidationStart = Clock::now();
        evaluateValidation(validation, candidateFallbackPrecision, true);
        result.centeredFallbackLowPrecisionSeconds += secondsSince(lowValidationStart);
        if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;
        if (resourceError.load(std::memory_order_relaxed)) {
            error = "adaptive centered low-precision validation allocation failed";
            return Centered4ExperimentalStatus::ResourceLimit;
        }
        if (internalError.load(std::memory_order_relaxed)) {
            error = "adaptive centered low-precision validation failed";
            return Centered4ExperimentalStatus::Failed;
        }
        const Clock::time_point fullValidationStart = Clock::now();
        evaluateValidation(validation, fallbackPrecision, false);
        result.centeredFallbackFullPrecisionSeconds += secondsSince(fullValidationStart);
        if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;
        if (resourceError.load(std::memory_order_relaxed)) {
            error = "adaptive centered full-precision validation allocation failed";
            return Centered4ExperimentalStatus::ResourceLimit;
        }
        if (internalError.load(std::memory_order_relaxed)) {
            error = "adaptive centered full-precision validation failed";
            return Centered4ExperimentalStatus::Failed;
        }
        result.centeredFallbackLowPrecisionPixelCount = validation.size();
        result.centeredFallbackFullPrecisionPixelCount = validation.size();
        for (const FallbackValidationResult& sample : validation) {
            result.centeredFallbackLowPrecisionIterationCount = saturatingAdd(result.centeredFallbackLowPrecisionIterationCount, sample.lowIterations);
            result.centeredFallbackFullPrecisionIterationCount = saturatingAdd(result.centeredFallbackFullPrecisionIterationCount, sample.fullIterations);
            if (!sample.lowDefined || !sample.fullDefined || sample.lowOutput != sample.fullOutput) ++result.centeredFallbackValidationMismatchCount;
        }
        result.centeredFallbackUpgraded = result.centeredFallbackValidationMismatchCount != 0;
    }

    struct FallbackBatchResult {
        uint64_t iterations = 0;
        uint64_t undefinedPixels = 0;
        double seconds = 0.0;
    };
    auto evaluateFallbackPixels = [&](const std::vector<size_t>& indices, mpfr_prec_t precision, uint8_t skippedMask) {
        std::atomic<uint64_t> batchIterations{0};
        std::atomic<uint64_t> batchUndefinedPixels{0};
        const Clock::time_point start = Clock::now();
        if (indices.empty()) return FallbackBatchResult{0, 0, secondsSince(start)};
#pragma omp parallel num_threads(threadCount)
        {
            OracleWorkspaceRelease workerOracleWorkspaceRelease;
            std::unique_ptr<FallbackWorkerWorkspace> workspace;
            try {
                workspace = std::make_unique<FallbackWorkerWorkspace>(precision, request);
                if (!workspace->geometryReady) internalError.store(true, std::memory_order_relaxed);
            } catch (const std::bad_alloc&) {
                resourceError.store(true, std::memory_order_relaxed);
            } catch (const std::length_error&) {
                resourceError.store(true, std::memory_order_relaxed);
            } catch (...) {
                internalError.store(true, std::memory_order_relaxed);
            }
            uint64_t localIterations = 0;
#pragma omp for schedule(dynamic, 1)
            for (long long rawIndex = 0; rawIndex < static_cast<long long>(indices.size()); ++rawIndex) {
                try {
                    if (shouldCancel() || !workspace || resourceError.load(std::memory_order_relaxed) || internalError.load(std::memory_order_relaxed)) continue;
                    const size_t outputIndex = indices[static_cast<size_t>(rawIndex)];
                    if (skippedMask != 0 && fallbackMask[outputIndex] == skippedMask) continue;
                    float output = ExpressionDeepEmptyPixel;
                    uint64_t iterations = 0;
                    const bool defined = exactPixel(*workspace, outputIndex, iterations, output);
                    localIterations += iterations;
                    if (!defined) {
                        batchUndefinedPixels.fetch_add(1, std::memory_order_relaxed);
                        continue;
                    }
                    request.output[outputIndex] = output;
                } catch (const std::bad_alloc&) {
                    resourceError.store(true, std::memory_order_relaxed);
                } catch (const std::length_error&) {
                    resourceError.store(true, std::memory_order_relaxed);
                } catch (...) {
                    internalError.store(true, std::memory_order_relaxed);
                }
            }
            batchIterations.fetch_add(localIterations, std::memory_order_relaxed);
        }
        return FallbackBatchResult{batchIterations.load(std::memory_order_relaxed), batchUndefinedPixels.load(std::memory_order_relaxed), secondsSince(start)};
    };

    const size_t remainingMandatoryPixelCount = mandatoryFullFallback.size() - mandatoryValidation.size();
    mpfr_prec_t mandatoryRemainingPrecision = result.centeredFallbackMandatoryUpgraded || !useMandatoryCandidate ? fallbackPrecision : mandatoryCandidatePrecision;
    FallbackBatchResult mandatoryResult = evaluateFallbackPixels(mandatoryFullFallback, mandatoryRemainingPrecision, 4);
    if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;
    if (resourceError.load(std::memory_order_relaxed)) {
        error = "adaptive centered mandatory-full fallback allocation failed";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    if (internalError.load(std::memory_order_relaxed)) {
        error = "adaptive centered mandatory-full fallback workspace failed";
        return Centered4ExperimentalStatus::Failed;
    }
    if (mandatoryRemainingPrecision != fallbackPrecision && mandatoryResult.undefinedPixels != 0) {
        result.centeredFallbackMandatoryUpgraded = true;
        result.centeredFallbackMandatoryCandidatePixelCount = saturatingAdd(result.centeredFallbackMandatoryCandidatePixelCount, remainingMandatoryPixelCount);
        result.centeredFallbackMandatoryCandidateIterationCount = saturatingAdd(result.centeredFallbackMandatoryCandidateIterationCount, mandatoryResult.iterations);
        result.centeredFallbackMandatoryCandidateSeconds += mandatoryResult.seconds;
        mandatoryRemainingPrecision = fallbackPrecision;
        mandatoryResult = evaluateFallbackPixels(mandatoryFullFallback, mandatoryRemainingPrecision, 4);
        if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;
        if (resourceError.load(std::memory_order_relaxed)) {
            error = "adaptive centered upgraded mandatory fallback allocation failed";
            return Centered4ExperimentalStatus::ResourceLimit;
        }
        if (internalError.load(std::memory_order_relaxed)) {
            error = "adaptive centered upgraded mandatory fallback workspace failed";
            return Centered4ExperimentalStatus::Failed;
        }
    }
    if (mandatoryRemainingPrecision == fallbackPrecision) {
        result.centeredFallbackMandatoryFullPrecisionPixelCount = saturatingAdd(result.centeredFallbackMandatoryFullPrecisionPixelCount, remainingMandatoryPixelCount);
        result.centeredFallbackMandatoryFullPrecisionIterationCount = saturatingAdd(result.centeredFallbackMandatoryFullPrecisionIterationCount, mandatoryResult.iterations);
        result.centeredFallbackMandatoryFullPrecisionSeconds += mandatoryResult.seconds;
    } else {
        result.centeredFallbackMandatoryCandidatePixelCount = saturatingAdd(result.centeredFallbackMandatoryCandidatePixelCount, remainingMandatoryPixelCount);
        result.centeredFallbackMandatoryCandidateIterationCount = saturatingAdd(result.centeredFallbackMandatoryCandidateIterationCount, mandatoryResult.iterations);
        result.centeredFallbackMandatoryCandidateSeconds += mandatoryResult.seconds;
    }
    for (const FallbackValidationResult& sample : mandatoryValidation) {
        const bool useFull = result.centeredFallbackMandatoryUpgraded;
        const bool defined = useFull ? sample.fullDefined : sample.lowDefined;
        if (!defined) continue;
        request.output[sample.index] = useFull ? sample.fullOutput : sample.lowOutput;
    }
    if (result.centeredFallbackMandatoryUpgraded) result.centeredFallbackMandatoryUpgradedPixelCount = mandatoryFullFallback.size();
    result.centeredFallbackMandatoryFullIterationCount = mandatoryRemainingPrecision == fallbackPrecision ? result.centeredFallbackMandatoryFullPrecisionIterationCount : result.centeredFallbackMandatoryCandidateIterationCount;
    result.centeredFallbackMandatoryFullSeconds = mandatoryRemainingPrecision == fallbackPrecision ? result.centeredFallbackMandatoryFullPrecisionSeconds : result.centeredFallbackMandatoryCandidateSeconds;
    result.centeredFallbackFullPrecisionPixelCount = saturatingAdd(result.centeredFallbackFullPrecisionPixelCount, result.centeredFallbackMandatoryFullPrecisionPixelCount);
    result.centeredFallbackFullPrecisionIterationCount = saturatingAdd(result.centeredFallbackFullPrecisionIterationCount, result.centeredFallbackMandatoryFullPrecisionIterationCount);
    result.centeredFallbackFullPrecisionSeconds += result.centeredFallbackMandatoryFullPrecisionSeconds;

    const size_t remainingLowEligiblePixelCount = lowEligibleFallback.size() - validation.size();
    mpfr_prec_t remainingPrecision = result.centeredFallbackUpgraded || !useLowPrecisionFallback ? fallbackPrecision : candidateFallbackPrecision;
    FallbackBatchResult remainingResult = evaluateFallbackPixels(lowEligibleFallback, remainingPrecision, 3);
    if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;
    if (resourceError.load(std::memory_order_relaxed)) {
        error = "adaptive centered low-eligible fallback allocation failed";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    if (internalError.load(std::memory_order_relaxed)) {
        error = "adaptive centered low-eligible fallback workspace failed";
        return Centered4ExperimentalStatus::Failed;
    }
    if (remainingPrecision != fallbackPrecision && remainingResult.undefinedPixels != 0) {
        result.centeredFallbackUpgraded = true;
        result.centeredFallbackLowPrecisionPixelCount = saturatingAdd(result.centeredFallbackLowPrecisionPixelCount, remainingLowEligiblePixelCount);
        result.centeredFallbackLowPrecisionIterationCount = saturatingAdd(result.centeredFallbackLowPrecisionIterationCount, remainingResult.iterations);
        result.centeredFallbackLowPrecisionSeconds += remainingResult.seconds;
        remainingPrecision = fallbackPrecision;
        remainingResult = evaluateFallbackPixels(lowEligibleFallback, remainingPrecision, 3);
        if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;
        if (resourceError.load(std::memory_order_relaxed)) {
            error = "adaptive centered upgraded fallback allocation failed";
            return Centered4ExperimentalStatus::ResourceLimit;
        }
        if (internalError.load(std::memory_order_relaxed)) {
            error = "adaptive centered upgraded fallback workspace failed";
            return Centered4ExperimentalStatus::Failed;
        }
    }
    if (remainingPrecision == fallbackPrecision) {
        result.centeredFallbackFullPrecisionPixelCount = saturatingAdd(result.centeredFallbackFullPrecisionPixelCount, remainingLowEligiblePixelCount);
        result.centeredFallbackFullPrecisionIterationCount = saturatingAdd(result.centeredFallbackFullPrecisionIterationCount, remainingResult.iterations);
        result.centeredFallbackFullPrecisionSeconds += remainingResult.seconds;
    } else {
        result.centeredFallbackLowPrecisionPixelCount = saturatingAdd(result.centeredFallbackLowPrecisionPixelCount, remainingLowEligiblePixelCount);
        result.centeredFallbackLowPrecisionIterationCount = saturatingAdd(result.centeredFallbackLowPrecisionIterationCount, remainingResult.iterations);
        result.centeredFallbackLowPrecisionSeconds += remainingResult.seconds;
    }

    uint64_t undefinedPixels = mandatoryResult.undefinedPixels + remainingResult.undefinedPixels;
    for (const FallbackValidationResult& sample : mandatoryValidation) {
        const bool useFull = result.centeredFallbackMandatoryUpgraded;
        if (!(useFull ? sample.fullDefined : sample.lowDefined)) ++undefinedPixels;
    }
    for (const FallbackValidationResult& sample : validation) {
        const bool useFull = result.centeredFallbackUpgraded;
        const bool defined = useFull ? sample.fullDefined : sample.lowDefined;
        if (!defined) {
            ++undefinedPixels;
            continue;
        }
        request.output[sample.index] = useFull ? sample.fullOutput : sample.lowOutput;
    }
    if (result.centeredFallbackUpgraded) result.centeredFallbackUpgradedPixelCount = lowEligibleFallback.size();
    result.fallbackSeconds += secondsSince(fallbackStart);
    if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;
    result.undefinedPixelCount = undefinedPixels;
    if (result.undefinedPixelCount != 0) {
        std::fill(request.output, request.output + pixelCount, ExpressionDeepEmptyPixel);
        error = "adaptive centered fallback produced an undefined pixel";
        return Centered4ExperimentalStatus::UndefinedPixel;
    }
    const uint64_t centeredFastIterations = saturatingAdd(vmIterations.load(std::memory_order_relaxed), doubleDoubleIterations.load(std::memory_order_relaxed));
    result.centeredVectorStepCount = saturatingAdd(result.centeredVectorStepCount, vectorSteps.load(std::memory_order_relaxed));
    result.centeredVectorActiveLaneCount = saturatingAdd(result.centeredVectorActiveLaneCount, vectorActiveLanes.load(std::memory_order_relaxed));
    result.fastIterationCount = saturatingAdd(result.fastIterationCount, centeredFastIterations);
    const uint64_t fallbackIterations = saturatingAdd(result.centeredFallbackMandatoryCandidateIterationCount, saturatingAdd(result.centeredFallbackLowPrecisionIterationCount, result.centeredFallbackFullPrecisionIterationCount));
    result.totalIterations = saturatingAdd(result.totalIterations, saturatingAdd(centeredFastIterations, fallbackIterations));
    result.selectedPrecision = selectedPrecision;
    result.fallbackPrecision = fallbackPrecision;
    result.centeredAccepted = true;
    result.status = ExpressionDeepRenderStatus::Success;
    result.success = true;
    progress(ExpressionDeepRenderPhase::Complete, pixelCount, pixelCount);
    return Centered4ExperimentalStatus::Success;
}

template <typename Cancel, typename Progress>
Centered4ExperimentalStatus renderCentered4Experimental(const ExpressionDeepRenderRequest& request, ExpressionDeepRenderResult& result, mpfr_prec_t selectedPrecision, mpfr_prec_t fallbackPrecision, int threadCount, size_t outerLiveBytes, Cancel&& shouldCancel, Progress&& progress, std::string& error) {
    struct OracleWorkspaceRelease {
        ~OracleWorkspaceRelease() { ExpressionOracle::releaseThreadWorkspace(); }
    } callerOracleWorkspaceRelease;
    const char* setting = std::getenv("MANDEL_DEEP_CENTERED4");
    const ExpressionScaledResidualCapability capability = request.runtimeProgram->scaledResidualCapability();
    const bool analyticCapability = capability == ExpressionScaledResidualCapability::ExactCenteredArithmetic || capability == ExpressionScaledResidualCapability::CertifiedEntireCandidate;
    if (!setting || std::atoi(setting) == 0 || request.forceMpfrFallbackForVerification || request.pixelParameter != FormulaParameter::C || !request.runtimeProgram->derivativeCompatible() || !analyticCapability || selectedPrecision > 512 || selectedPrecision > MPFR_PREC_MAX - request.memory.fallbackGuardBits || fallbackPrecision > request.precision.maximumBits || fallbackPrecision > ExpressionReferencePrecisionPolicy::ApplicationMaximumBits) return Centered4ExperimentalStatus::NotApplicable;
    if (request.memory.memoryLimitBytes != 0 && outerLiveBytes >= request.memory.memoryLimitBytes) return Centered4ExperimentalStatus::NotApplicable;
    const size_t centeredMemoryBudget = request.memory.memoryLimitBytes == 0 ? 0 : request.memory.memoryLimitBytes - outerLiveBytes;
    constexpr mpfr_prec_t ReferencePrecision = 512;
    const size_t pixelCount = static_cast<size_t>(request.width) * request.height;
    ExactGeometry geometry(ReferencePrecision);
    if (!geometry.initialize(request)) {
        error = "centered4 geometry initialization failed";
        return Centered4ExperimentalStatus::Failed;
    }

    struct Candidate {
        long long distance = 0;
        int x = 0;
        int y = 0;
    };
    std::vector<Candidate> candidates;
    constexpr int GridWidth = 9;
    constexpr int GridHeight = 9;
    for (int gridY = 0; gridY < GridHeight; ++gridY) {
        const int y = static_cast<int>(static_cast<long long>(gridY) * (request.height - 1) / (GridHeight - 1));
        for (int gridX = 0; gridX < GridWidth; ++gridX) {
            const int x = static_cast<int>(static_cast<long long>(gridX) * (request.width - 1) / (GridWidth - 1));
            const long long dx = 2LL * x - (request.width - 1LL);
            const long long dy = 2LL * y - (request.height - 1LL);
            candidates.push_back({dx * dx + dy * dy, x, y});
        }
    }
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& left, const Candidate& right) { return left.distance < right.distance || (left.distance == right.distance && (left.y < right.y || (left.y == right.y && left.x < right.x))); });
    candidates.erase(std::unique(candidates.begin(), candidates.end(), [](const Candidate& left, const Candidate& right) { return left.x == right.x && left.y == right.y; }), candidates.end());

    const Clock::time_point referenceStart = Clock::now();
    ExpressionReferenceOrbitResult reference;
    int referenceX = -1;
    int referenceY = -1;
    MpfrComplex candidatePixel(ReferencePrecision);
    mpf_t candidateReal;
    mpf_t candidateImaginary;
    mpf_init2(candidateReal, ReferencePrecision);
    mpf_init2(candidateImaginary, ReferencePrecision);
    for (const Candidate& candidate : candidates) {
        if (shouldCancel()) {
            mpf_clears(candidateReal, candidateImaginary, (mpf_ptr)0);
            return Centered4ExperimentalStatus::Cancelled;
        }
        if (!geometry.coordinate(candidate.x, candidate.y, request, candidatePixel)) continue;
        mpfr_get_f(candidateReal, candidatePixel.re, MPFR_RNDN);
        mpfr_get_f(candidateImaginary, candidatePixel.im, MPFR_RNDN);
        ExpressionReferenceBuildRequest referenceRequest;
        referenceRequest.canonicalProgram = request.canonicalProgram;
        referenceRequest.runtimeProgram = request.runtimeProgram;
        referenceRequest.pixelParameter = request.pixelParameter;
        referenceRequest.center.realMpf = candidateReal;
        referenceRequest.center.imaginaryMpf = candidateImaginary;
        referenceRequest.fixed = request.fixed;
        referenceRequest.bailout = request.bailout;
        referenceRequest.maxIterations = request.maxIterations;
        referenceRequest.precision.requestedBits = ReferencePrecision;
        referenceRequest.precision.minimumBits = 53;
        referenceRequest.precision.guardBits = 0;
        referenceRequest.precision.maximumBits = request.precision.maximumBits;
        referenceRequest.memoryLimitBytes = centeredMemoryBudget;
        referenceRequest.shouldCancel = shouldCancel;
        ExpressionReferenceOrbitResult candidateReference;
        if (!buildExpressionReferenceOrbit(referenceRequest, candidateReference)) {
            if (candidateReference.cancelled) {
                mpf_clears(candidateReal, candidateImaginary, (mpf_ptr)0);
                return Centered4ExperimentalStatus::Cancelled;
            }
            continue;
        }
        if (candidateReference.valid && !candidateReference.escaped && !candidateReference.undefined && candidateReference.samples.size() == static_cast<size_t>(request.maxIterations)) {
            reference = std::move(candidateReference);
            referenceX = candidate.x;
            referenceY = candidate.y;
            break;
        }
    }
    mpf_clears(candidateReal, candidateImaginary, (mpf_ptr)0);
    const double centeredReferenceSeconds = secondsSince(referenceStart);
    ExpressionOracle::releaseThreadWorkspace();
    if (referenceX < 0 || referenceY < 0) return Centered4ExperimentalStatus::NotApplicable;

    size_t rendererBytes = 0;
    if (!checkedAddSize(rendererBytes, pixelCount, sizeof(double) + sizeof(size_t) + sizeof(uint8_t)) || !checkedAddSize(rendererBytes, static_cast<size_t>(request.width + request.height), sizeof(double)) || !checkedAddSize(rendererBytes, reference.samples.size(), sizeof(Complex)) || !checkedAddSize(rendererBytes, reference.tape.size(), sizeof(Complex))) {
        error = "centered4 workspace size overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    const size_t limbs = (static_cast<size_t>(fallbackPrecision) + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
    size_t mpfrValueBytes = sizeof(__mpfr_struct);
    if (!checkedAddSize(mpfrValueBytes, limbs, sizeof(mp_limb_t))) {
        error = "centered4 MPFR workspace size overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    size_t mpfrValuesPerThread = 62;
    if (!checkedAddSize(mpfrValuesPerThread, request.runtimeProgram->stackDepth(), 2)) {
        error = "centered4 oracle stack size overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    size_t fallbackThreadBytes = 0;
    if (!checkedAddSize(fallbackThreadBytes, mpfrValuesPerThread, mpfrValueBytes)) {
        error = "centered4 fallback workspace size overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    size_t fallbackThreadBytesTotal = 0;
    if (!checkedAddSize(fallbackThreadBytesTotal, static_cast<size_t>(threadCount), fallbackThreadBytes) || !checkedAddSize(rendererBytes, 1, fallbackThreadBytesTotal)) {
        error = "centered4 fallback workspace multiplication overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    if (!checkedAddSize(rendererBytes, static_cast<size_t>(threadCount), size_t{64} << 10)) {
        error = "centered4 VM thread workspace multiplication overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    if (request.memory.memoryLimitBytes != 0 && (outerLiveBytes > request.memory.memoryLimitBytes || reference.memoryBytes > request.memory.memoryLimitBytes - outerLiveBytes || rendererBytes > request.memory.memoryLimitBytes - outerLiveBytes - reference.memoryBytes)) {
        return Centered4ExperimentalStatus::NotApplicable;
    }
    auto reconstruct = [](const ScaledComplexShadow& primary, const ScaledComplexShadow& defect, Complex& output) {
        ScaledComplexValue scaled;
        return makeScaledComplexValue(primary, defect, scaled) == ScaledArithmeticStatus::Success && scaledValueToDouble(scaled, output);
    };
    Complex referenceC;
    Complex referenceZ0;
    if (!reconstruct(reference.c, reference.cDefect, referenceC) || !reconstruct(reference.z0, reference.z0Defect, referenceZ0)) return Centered4ExperimentalStatus::NotApplicable;
    std::vector<Complex> referenceZ(reference.samples.size());
    std::vector<Complex> referenceNodes(reference.tape.size());
    for (size_t index = 0; index < reference.tape.size(); ++index)
        if (!reconstruct(reference.tape[index].output, reference.tape[index].outputDefect, referenceNodes[index])) return Centered4ExperimentalStatus::NotApplicable;
    for (size_t iteration = 0; iteration < reference.samples.size(); ++iteration)
        if (!reconstruct(reference.samples[iteration].z, reference.samples[iteration].zDefect, referenceZ[iteration])) return Centered4ExperimentalStatus::NotApplicable;

    MpfrComplex referencePixel(ReferencePrecision);
    if (!geometry.coordinate(referenceX, referenceY, request, referencePixel)) return Centered4ExperimentalStatus::NotApplicable;
    std::vector<double> xOffsets(static_cast<size_t>(request.width));
    std::vector<double> yOffsets(static_cast<size_t>(request.height));
    MpfrComplex coordinate(ReferencePrecision);
    for (int x = 0; x < request.width; ++x) {
        if (!geometry.coordinate(x, referenceY, request, coordinate)) return Centered4ExperimentalStatus::NotApplicable;
        mpfr_sub(coordinate.re, coordinate.re, referencePixel.re, MPFR_RNDN);
        const double offset = mpfr_get_d(coordinate.re, MPFR_RNDN);
        if (!std::isfinite(offset) || (!mpfr_zero_p(coordinate.re) && offset == 0.0)) return Centered4ExperimentalStatus::NotApplicable;
        xOffsets[static_cast<size_t>(x)] = offset;
    }
    for (int y = 0; y < request.height; ++y) {
        if (!geometry.coordinate(referenceX, y, request, coordinate)) return Centered4ExperimentalStatus::NotApplicable;
        mpfr_sub(coordinate.im, coordinate.im, referencePixel.im, MPFR_RNDN);
        const double offset = mpfr_get_d(coordinate.im, MPFR_RNDN);
        if (!std::isfinite(offset) || (!mpfr_zero_p(coordinate.im) && offset == 0.0)) return Centered4ExperimentalStatus::NotApplicable;
        yOffsets[static_cast<size_t>(y)] = offset;
    }
    result.referenceSeconds += centeredReferenceSeconds;
    result.referenceBytes += reference.memoryBytes;
    result.rendererBytes += rendererBytes;

    double errorThreshold = 1e5;
    if (const char* value = std::getenv("MANDEL_DEEP_CENTERED4_ERROR")) {
        const double parsed = std::strtod(value, nullptr);
        if (std::isfinite(parsed) && parsed >= 0.0) errorThreshold = parsed;
    }
    std::vector<double> errors(pixelCount, 0.0);
    std::atomic<uint64_t> vmIterations{0};
    progress(ExpressionDeepRenderPhase::Fast, 0, pixelCount);
    const Clock::time_point fastStart = Clock::now();
#pragma omp parallel for num_threads(threadCount) schedule(dynamic, 1)
    for (int y = 0; y < request.height; ++y) {
        if (shouldCancel()) continue;
        ExpressionContext referenceContext = request.fixed;
        referenceContext.c = referenceC;
        referenceContext.z0 = referenceZ0;
        uint64_t localIterations = 0;
        for (int xBase = 0; xBase < request.width; xBase += 4) {
            const int lanes = std::min(4, request.width - xBase);
            std::array<ExpressionDeltaContext, 4> deltas{};
            std::array<Complex, 4> stateDelta{};
            std::array<double, 4> stateError{};
            std::array<bool, 4> active{};
            std::array<size_t, 4> outputIndices{};
            int activeCount = lanes;
            for (int lane = 0; lane < lanes; ++lane) {
                const int x = xBase + lane;
                outputIndices[lane] = static_cast<size_t>(y) * request.width + x;
                deltas[lane].c = {xOffsets[static_cast<size_t>(x)], yOffsets[static_cast<size_t>(y)]};
                stateError[lane] = 1e-15;
                active[lane] = true;
                request.output[outputIndices[lane]] = ExpressionDeepInteriorPixel;
            }
            for (int iteration = 0; iteration < request.maxIterations && activeCount > 0; ++iteration) {
                if ((iteration & 63) == 0 && shouldCancel()) break;
                referenceContext.iteration = iteration;
                referenceContext.z = referenceZ[static_cast<size_t>(iteration)];
                for (int lane = 0; lane < 4; ++lane) deltas[lane].z = stateDelta[lane];
                const ExpressionReferenceSample& sample = reference.samples[static_cast<size_t>(iteration)];
                const Complex* nodeBases = referenceNodes.data() + static_cast<size_t>(sample.tapeOffset);
                ExpressionCenteredResult evaluated[4];
                if (!ExpressionCenteredEvaluator::evaluate4WithNodeBases(*request.runtimeProgram, referenceContext, nodeBases, deltas.data(), evaluated)) {
                    for (int lane = 0; lane < lanes; ++lane)
                        if (active[lane]) {
                            errors[outputIndices[lane]] = std::numeric_limits<double>::infinity();
                            active[lane] = false;
                            --activeCount;
                        }
                    break;
                }
                const Complex nextReference = nodeBases[sample.rootNode];
                for (int lane = 0; lane < lanes; ++lane) {
                    if (!active[lane]) continue;
                    ++localIterations;
                    if (!evaluated[lane].success()) {
                        errors[outputIndices[lane]] = std::numeric_limits<double>::infinity();
                        active[lane] = false;
                        --activeCount;
                        continue;
                    }
                    const Complex currentActual = referenceZ[static_cast<size_t>(iteration)] + stateDelta[lane];
                    const Complex outputDelta = evaluated[lane].delta;
                    const Complex actual = nextReference + outputDelta;
                    stateDelta[lane] = iteration + 1 < request.maxIterations ? outputDelta + (nextReference - referenceZ[static_cast<size_t>(iteration + 1)]) : outputDelta;
                    ExpressionContext actualContext = request.fixed;
                    actualContext.z = currentActual;
                    actualContext.c = referenceC + deltas[lane].c;
                    actualContext.iteration = iteration;
                    ExpressionDerivativeSeed seed;
                    seed.z = {1.0, 0.0};
                    Complex value;
                    Complex derivative;
                    if (request.runtimeProgram->evaluateWithDerivative(actualContext, seed, value, derivative) && finiteComplex(derivative)) {
                        stateError[lane] = std::abs(derivative) * stateError[lane] + 64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(actual));
                    } else {
                        stateError[lane] = std::numeric_limits<double>::infinity();
                    }
                    errors[outputIndices[lane]] = stateError[lane];
                    if (!std::isfinite(stateError[lane]) || stateError[lane] > errorThreshold) {
                        active[lane] = false;
                        --activeCount;
                        continue;
                    }
                    const double magnitude = std::hypot(actual.real(), actual.imag());
                    const double magnitudeError = 2.0 * stateError[lane];
                    if (!finiteComplex(actual)) {
                        errors[outputIndices[lane]] = std::numeric_limits<double>::infinity();
                        active[lane] = false;
                        --activeCount;
                    } else if (magnitude > request.bailout && magnitude - magnitudeError > request.bailout) {
                        request.output[outputIndices[lane]] = static_cast<float>(iteration + 1);
                        active[lane] = false;
                        --activeCount;
                    } else if (magnitude + magnitudeError >= request.bailout) {
                        errors[outputIndices[lane]] = std::numeric_limits<double>::infinity();
                        active[lane] = false;
                        --activeCount;
                    }
                }
            }
        }
        vmIterations.fetch_add(localIterations, std::memory_order_relaxed);
    }
    result.fastSeconds += secondsSince(fastStart);
    if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;

    std::vector<size_t> fallback;
    std::vector<uint8_t> fallbackMask;
    try {
        fallback.reserve(pixelCount);
        for (size_t index = 0; index < pixelCount; ++index)
            if (!std::isfinite(errors[index]) || errors[index] > errorThreshold) fallback.push_back(index);
        fallbackMask.assign(pixelCount, 0);
    } catch (const std::bad_alloc&) {
        error = "centered4 fallback queue allocation failed";
        return Centered4ExperimentalStatus::ResourceLimit;
    } catch (const std::length_error&) {
        error = "centered4 fallback queue length overflow";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    result.fastPixelCount = static_cast<uint64_t>(pixelCount - fallback.size());
    result.fallbackPixelCount = static_cast<uint64_t>(fallback.size());
    result.fallbackReasonCounts.fill(0);
    result.uncertainPixelCount = 0;
    result.maxFallbackReasonCount = 0;
    result.fallbackTileCount = 0;
    result.maxTileFallbackRate = 0.0;
    result.undefinedPixelCount = 0;
    result.fallbackReasonCounts[static_cast<size_t>(ExpressionDeepFallbackReason::CertificationFailure)] = result.fallbackPixelCount;
    result.maxFallbackReasonCount = result.fallbackPixelCount;
    result.uncertainPixelCount = result.fallbackPixelCount;
    for (size_t index : fallback) fallbackMask[index] = 1;
    const int tilesX = (request.width + request.threading.tileWidth - 1) / request.threading.tileWidth;
    const int tilesY = (request.height + request.threading.tileHeight - 1) / request.threading.tileHeight;
    for (int tileY = 0; tileY < tilesY; ++tileY) {
        for (int tileX = 0; tileX < tilesX; ++tileX) {
            const int xBegin = tileX * request.threading.tileWidth;
            const int yBegin = tileY * request.threading.tileHeight;
            const int xEnd = std::min(request.width, xBegin + request.threading.tileWidth);
            const int yEnd = std::min(request.height, yBegin + request.threading.tileHeight);
            uint64_t fallbackInTile = 0;
            for (int y = yBegin; y < yEnd; ++y)
                for (int x = xBegin; x < xEnd; ++x) fallbackInTile += fallbackMask[static_cast<size_t>(y) * request.width + x];
            if (fallbackInTile == 0) continue;
            ++result.fallbackTileCount;
            const uint64_t pixelsInTile = static_cast<uint64_t>(xEnd - xBegin) * static_cast<uint64_t>(yEnd - yBegin);
            result.maxTileFallbackRate = std::max(result.maxTileFallbackRate, static_cast<double>(fallbackInTile) / pixelsInTile);
        }
    }

    progress(ExpressionDeepRenderPhase::Fallback, 0, fallback.size());
    const Clock::time_point fallbackStart = Clock::now();
    std::atomic<uint64_t> fallbackIterations{0};
    std::atomic<uint64_t> undefinedPixels{0};
    std::atomic_bool resourceError{false};
    std::atomic_bool internalError{false};
#pragma omp parallel num_threads(threadCount)
    {
        OracleWorkspaceRelease workerOracleWorkspaceRelease;
        std::unique_ptr<FallbackWorkerWorkspace> workspace;
        try {
            workspace = std::make_unique<FallbackWorkerWorkspace>(fallbackPrecision, request);
            if (!workspace->geometryReady) internalError.store(true, std::memory_order_relaxed);
        } catch (const std::bad_alloc&) {
            resourceError.store(true, std::memory_order_relaxed);
        } catch (const std::length_error&) {
            resourceError.store(true, std::memory_order_relaxed);
        } catch (...) {
            internalError.store(true, std::memory_order_relaxed);
        }
        uint64_t localIterations = 0;
#pragma omp for schedule(dynamic, 8)
        for (long long rawIndex = 0; rawIndex < static_cast<long long>(fallback.size()); ++rawIndex) {
            try {
                if (shouldCancel() || !workspace || resourceError.load(std::memory_order_relaxed) || internalError.load(std::memory_order_relaxed)) continue;
                const size_t outputIndex = fallback[static_cast<size_t>(rawIndex)];
                const int y = static_cast<int>(outputIndex / static_cast<size_t>(request.width));
                const int x = static_cast<int>(outputIndex % static_cast<size_t>(request.width));
                configureFixed(request.fixed, workspace->context);
                if (!workspace->geometry.coordinate(x, y, request, workspace->pixel)) {
                    request.output[outputIndex] = ExpressionDeepEmptyPixel;
                    undefinedPixels.fetch_add(1, std::memory_order_relaxed);
                    continue;
                }
                workspace->context.c.set(workspace->pixel);
                workspace->context.z.set(workspace->context.z0);
                bool outside = false;
                if (!piecewiseOutsideBailout(workspace->context.z, workspace->piecewiseBailoutSquared, request.bailout, workspace->magnitudeStorage.re, workspace->piecewiseSquares, workspace->piecewiseInsideSquareExponent, workspace->piecewiseScratch, outside)) {
                    request.output[outputIndex] = ExpressionDeepEmptyPixel;
                    undefinedPixels.fetch_add(1, std::memory_order_relaxed);
                    continue;
                }
                if (outside) {
                    request.output[outputIndex] = 0.0f;
                    continue;
                }
                request.output[outputIndex] = ExpressionDeepInteriorPixel;
                for (int iteration = 0; iteration < request.maxIterations; ++iteration) {
                    if ((iteration & 15) == 0 && shouldCancel()) break;
                    workspace->context.iteration = iteration;
                    const bool defined = ExpressionOracle::evaluateOrbitStep(*request.runtimeProgram, workspace->context, workspace->next, &workspace->piecewiseSquares, nullptr);
                    if (mpfr_nan_p(workspace->next.re) || mpfr_nan_p(workspace->next.im)) {
                        request.output[outputIndex] = ExpressionDeepEmptyPixel;
                        undefinedPixels.fetch_add(1, std::memory_order_relaxed);
                        break;
                    }
                    ++localIterations;
                    if (mpfr_inf_p(workspace->next.re) || mpfr_inf_p(workspace->next.im)) {
                        request.output[outputIndex] = static_cast<float>(iteration + 1);
                        break;
                    }
                    if (!defined || !mpfr_number_p(workspace->next.re) || !mpfr_number_p(workspace->next.im)) {
                        request.output[outputIndex] = ExpressionDeepEmptyPixel;
                        undefinedPixels.fetch_add(1, std::memory_order_relaxed);
                        break;
                    }
                    mpfr_swap(workspace->context.z.re, workspace->next.re);
                    mpfr_swap(workspace->context.z.im, workspace->next.im);
                    if (!piecewiseOutsideBailout(workspace->context.z, workspace->piecewiseBailoutSquared, request.bailout, workspace->magnitudeStorage.re, workspace->piecewiseSquares, workspace->piecewiseInsideSquareExponent, workspace->piecewiseScratch, outside)) {
                        request.output[outputIndex] = ExpressionDeepEmptyPixel;
                        undefinedPixels.fetch_add(1, std::memory_order_relaxed);
                        break;
                    }
                    if (outside) {
                        request.output[outputIndex] = static_cast<float>(iteration + 1);
                        break;
                    }
                }
            } catch (const std::bad_alloc&) {
                resourceError.store(true, std::memory_order_relaxed);
            } catch (const std::length_error&) {
                resourceError.store(true, std::memory_order_relaxed);
            } catch (...) {
                internalError.store(true, std::memory_order_relaxed);
            }
        }
        fallbackIterations.fetch_add(localIterations, std::memory_order_relaxed);
    }
    result.fallbackSeconds += secondsSince(fallbackStart);
    if (shouldCancel()) return Centered4ExperimentalStatus::Cancelled;
    if (resourceError.load(std::memory_order_relaxed)) {
        error = "centered4 fallback allocation failed";
        return Centered4ExperimentalStatus::ResourceLimit;
    }
    if (internalError.load(std::memory_order_relaxed)) {
        error = "centered4 fallback workspace initialization failed";
        return Centered4ExperimentalStatus::Failed;
    }
    result.undefinedPixelCount = undefinedPixels.load(std::memory_order_relaxed);
    if (result.undefinedPixelCount != 0) {
        error = "centered4 fallback produced an undefined pixel";
        return Centered4ExperimentalStatus::UndefinedPixel;
    }
    const uint64_t centeredFastIterations = vmIterations.load(std::memory_order_relaxed);
    result.fastIterationCount = saturatingAdd(result.fastIterationCount, centeredFastIterations);
    result.totalIterations = saturatingAdd(result.totalIterations, saturatingAdd(centeredFastIterations, fallbackIterations.load(std::memory_order_relaxed)));
    result.selectedPrecision = selectedPrecision;
    result.fallbackPrecision = fallbackPrecision;
    result.status = ExpressionDeepRenderStatus::Success;
    result.success = true;
    progress(ExpressionDeepRenderPhase::Complete, pixelCount, pixelCount);
    return Centered4ExperimentalStatus::Success;
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

const char* expressionDeepCenteredFallbackReasonName(ExpressionDeepCenteredFallbackReason reason) {
    switch (reason) {
    case ExpressionDeepCenteredFallbackReason::CandidateFailure: return "candidate-failure";
    case ExpressionDeepCenteredFallbackReason::Nonfinite: return "nonfinite";
    case ExpressionDeepCenteredFallbackReason::ReferenceOutputDisagreement: return "reference-output-disagreement";
    case ExpressionDeepCenteredFallbackReason::StateDisagreement: return "state-disagreement";
    case ExpressionDeepCenteredFallbackReason::InteriorConservatism: return "interior-conservatism";
    case ExpressionDeepCenteredFallbackReason::DoubleDoubleRejected: return "double-double-rejected";
    case ExpressionDeepCenteredFallbackReason::StructuralInconsistency: return "structural-inconsistency";
    case ExpressionDeepCenteredFallbackReason::DenominatorInstability: return "denominator-instability";
    case ExpressionDeepCenteredFallbackReason::Count: break;
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
        const int fullHeight = request.fullHeight > 0 ? request.fullHeight : request.height;
        if (request.width < 2 || request.height < 1 || fullHeight < 2 || request.rowBase < 0 || request.rowBase > fullHeight - request.height || request.width > (LONG_MAX + 1LL) / 2 || fullHeight > (LONG_MAX + 1LL) / 2) return fail(ExpressionDeepRenderStatus::InvalidRequest, "frame dimensions are invalid");
        if (request.maxIterations < 1 || request.maxIterations > (1 << 24)) return fail(ExpressionDeepRenderStatus::InvalidRequest, "iteration count cannot be represented exactly in float output");
        if (!(request.bailout > 0.0) || !std::isfinite(request.bailout)) return fail(ExpressionDeepRenderStatus::InvalidRequest, "bailout must be finite and positive");
        if (!request.output) return fail(ExpressionDeepRenderStatus::InvalidRequest, "output buffer is null");
        if (request.threading.tileWidth < 1 || request.threading.tileHeight < 1 || request.threading.threads < 0 || request.threading.threads > 1024) return fail(ExpressionDeepRenderStatus::InvalidRequest, "thread or tile policy is invalid");
        if (request.memory.fallbackGuardBits < 0) return fail(ExpressionDeepRenderStatus::InvalidRequest, "fallback guard precision is invalid");
        if (request.preflight.maximumSamples < 1 || request.preflight.maximumSamples > 256 || request.preflight.minimumSamples < 1 || request.preflight.minimumSamples > request.preflight.maximumSamples || request.preflight.earlyRejectMinimumFirstUncertainIteration > static_cast<uint32_t>(1 << 24)) return fail(ExpressionDeepRenderStatus::InvalidRequest, "preflight policy is invalid");
        if (request.transfer.minimumPixelCount < 1 || request.transfer.minimumIterations < 1 || request.transfer.minimumIterations > (1 << 24)) return fail(ExpressionDeepRenderStatus::InvalidRequest, "transfer policy is invalid");
        if (request.taylor.minimumLanding < 1 || request.taylor.minimumOrder < 8 || request.taylor.maximumOrder > 20 || request.taylor.maximumBivariateOrder < 8 || request.taylor.maximumBivariateOrder > ExpressionTaylorMaximumBivariateOrder || request.taylor.minimumOrder > request.taylor.order || request.taylor.order > request.taylor.maximumOrder || request.taylor.maximumCompositionOrder < 1 || request.taylor.maximumCompositionOrder > 24 || request.taylor.maximumCandidateIteration < 0 || request.taylor.maximumDepth < 0 || request.taylor.maximumDepth > 30 || request.taylor.minimumTileWidth < 1 || request.taylor.minimumTileHeight < 1 || request.taylor.maximumJetCount < 1 || request.taylor.maximumJetCount > size_t{65536} || request.taylor.maximumRejectedBeforeFirstAcceptance < 1 || request.taylor.maximumRejectedBeforeFirstAcceptance > request.taylor.maximumJetCount || !(request.taylor.accuracyBudget > 0.0) || !std::isfinite(request.taylor.accuracyBudget))
            return fail(ExpressionDeepRenderStatus::InvalidRequest, "Taylor policy is invalid");
        if (request.verificationErrorInflationBits < 0 || request.verificationErrorInflationBits > 4096) return fail(ExpressionDeepRenderStatus::InvalidRequest, "verification error inflation is invalid");
        if (request.verificationCenteredFallbackCandidatePrecision < 0 || request.verificationCenteredFallbackCandidatePrecision > request.precision.maximumBits || request.verificationCenteredFallbackCandidatePrecision > ExpressionReferencePrecisionPolicy::ApplicationMaximumBits) return fail(ExpressionDeepRenderStatus::InvalidRequest, "verification centered fallback precision is invalid");
        if (request.verificationCenteredMandatoryFallbackCandidatePrecision < 0 || request.verificationCenteredMandatoryFallbackCandidatePrecision > request.precision.maximumBits || request.verificationCenteredMandatoryFallbackCandidatePrecision > ExpressionReferencePrecisionPolicy::ApplicationMaximumBits) return fail(ExpressionDeepRenderStatus::InvalidRequest, "verification centered mandatory fallback precision is invalid");
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
        const bool certifiedTransferCandidate = request.transfer.enableCertifiedSegments && result.capability == ExpressionScaledResidualCapability::ExactCenteredArithmetic && request.pixelParameter == FormulaParameter::C && !request.forceMpfrFallbackForVerification;
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
            size_t referenceCacheBytes = 0;
            if (!ExpressionScaledResidualEvaluator::estimateReferenceCacheBytes(maximumTapeNodes, referenceCacheBytes) || !checkedAddSize(fastThreadBytes, maximumTapeNodes, sizeof(uint8_t)) || !checkedAddSize(fastThreadBytes, 1, referenceCacheBytes)) return fail(ExpressionDeepRenderStatus::ResourceLimit, "reference tape workspace calculation overflow");
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
            referenceRequest.compaction = certifiedTransferCandidate ? ExpressionReferenceCompaction::FourTermCertifiedTransfer : ExpressionReferenceCompaction::TwoTerm;
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

            if (certifiedTransferCandidate) {
                result.transferAttempted = true;
                size_t availableTransferBytes = 0;
                if (request.memory.memoryLimitBytes != 0) {
                    if (reference.memoryBytes >= request.memory.memoryLimitBytes || rendererBaseBytes >= request.memory.memoryLimitBytes - reference.memoryBytes) {
                        availableTransferBytes = 1;
                    } else {
                        availableTransferBytes = request.memory.memoryLimitBytes - reference.memoryBytes - rendererBaseBytes;
                    }
                }
                const Clock::time_point transferBuildStart = Clock::now();
                const CertifiedTransferBuildResult transfer = buildCertifiedTerminalTransfer(request, reference, certificationGeometry, pixelCount64, availableTransferBytes, pollCancellation);
                result.transferBuildSeconds = secondsSince(transferBuildStart);
                result.transferMemoryBytes = transfer.memoryBytes;
                result.transferCoveredIterations = transfer.coveredIterations;
                result.transferFinalRadius = transfer.finalRadius;
                size_t transferPeak = rendererBaseBytes;
                if (!checkedAddSize(transferPeak, 1, transfer.memoryBytes)) return fail(ExpressionDeepRenderStatus::ResourceLimit, "transfer memory calculation overflow");
                result.rendererBytes = std::max(result.rendererBytes, transferPeak);
                if (transfer.status == CertifiedTransferBuildStatus::Cancelled || cancelled.load(std::memory_order_acquire)) return fail(ExpressionDeepRenderStatus::Cancelled, "render cancelled while building certified transfer segment");
                if (transfer.status == CertifiedTransferBuildStatus::Accepted) {
                    notifyProgress(ExpressionDeepRenderPhase::Fast, 0, pixelCount64);
                    const Clock::time_point transferApplyStart = Clock::now();
                    std::fill(request.output, request.output + pixelCount, ExpressionDeepInteriorPixel);
                    result.transferApplySeconds = secondsSince(transferApplyStart);
                    result.fastSeconds = result.transferApplySeconds;
                    result.transferAccepted = true;
                    result.transferAcceptedSegmentCount = 1;
                    result.transferSkippedIterationCount = saturatingMultiply(pixelCount64, static_cast<uint64_t>(request.maxIterations));
                    result.fastPixelCount = pixelCount64;
                    result.totalIterations = result.transferSkippedIterationCount;
                    notifyProgress(ExpressionDeepRenderPhase::Fast, pixelCount64, pixelCount64);
                    notifyProgress(ExpressionDeepRenderPhase::Complete, pixelCount64, pixelCount64);
                    result.status = ExpressionDeepRenderStatus::Success;
                    result.success = true;
                    result.error.clear();
                    return true;
                }
                runFast = false;
                certificationUnavailable = true;
                reference = {};
                result.referenceBytes = 0;
                rendererBaseBytes = 0;
                if (!checkedAddSize(rendererBaseBytes, pixelCount, sizeof(uint8_t)) || !checkedAddSize(rendererBaseBytes, pixelCount, sizeof(size_t))) return fail(ExpressionDeepRenderStatus::ResourceLimit, "transfer fallback renderer memory calculation overflow");
                fastThreadBytesTotal = 0;
                rendererBytes = rendererBaseBytes;
                result.rendererBytes = std::max(result.rendererBytes, transferPeak);
                goto certifiedTransferFallback;
            }

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
                const long centered = static_cast<long>(2LL * (y + request.rowBase) - (fullHeight - 1LL));
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
            if (runFast && certifiedTaylorCapability(result.capability) && !request.allowUncertifiedForBenchmark && !useTaylor) {
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
    certifiedTransferFallback:
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
        if (result.fastPixelCount == 0 && result.fallbackPixelCount == pixelCount64) {
            size_t centered4OuterLiveBytes = result.rendererBytes;
            if (!checkedAddSize(centered4OuterLiveBytes, 1, result.referenceBytes)) return fail(ExpressionDeepRenderStatus::ResourceLimit, "centered4 outer workspace size overflow");
            std::string centered4Error;
            const Centered4ExperimentalStatus adaptiveCenteredStatus = renderAdaptiveCenteredExperimental(request, result, result.selectedPrecision, result.fallbackPrecision, threadCount, centered4OuterLiveBytes, pollCancellation, notifyProgress, centered4Error);
            if (adaptiveCenteredStatus == Centered4ExperimentalStatus::Success) return true;
            if (adaptiveCenteredStatus == Centered4ExperimentalStatus::Cancelled) {
                std::fill(request.output, request.output + pixelCount, ExpressionDeepEmptyPixel);
                return fail(ExpressionDeepRenderStatus::Cancelled, "render cancelled during adaptive centered experiment");
            }
            if (adaptiveCenteredStatus == Centered4ExperimentalStatus::ResourceLimit) {
                std::fill(request.output, request.output + pixelCount, ExpressionDeepEmptyPixel);
                return fail(ExpressionDeepRenderStatus::ResourceLimit, centered4Error.empty() ? "adaptive centered resource limit" : centered4Error);
            }
            if (adaptiveCenteredStatus == Centered4ExperimentalStatus::UndefinedPixel) return fail(ExpressionDeepRenderStatus::UndefinedPixel, centered4Error.empty() ? "one or more adaptive centered pixels are undefined" : centered4Error);
            if (adaptiveCenteredStatus == Centered4ExperimentalStatus::Failed) {
                std::fill(request.output, request.output + pixelCount, ExpressionDeepEmptyPixel);
                return fail(ExpressionDeepRenderStatus::InternalError, centered4Error.empty() ? "adaptive centered experiment failed" : centered4Error);
            }
            centered4Error.clear();
            const Centered4ExperimentalStatus centered4Status = renderCentered4Experimental(request, result, result.selectedPrecision, result.fallbackPrecision, threadCount, centered4OuterLiveBytes, pollCancellation, notifyProgress, centered4Error);
            if (centered4Status == Centered4ExperimentalStatus::Success) return true;
            if (centered4Status == Centered4ExperimentalStatus::Cancelled) {
                std::fill(request.output, request.output + pixelCount, ExpressionDeepEmptyPixel);
                return fail(ExpressionDeepRenderStatus::Cancelled, "render cancelled during centered4 experiment");
            }
            if (centered4Status == Centered4ExperimentalStatus::ResourceLimit) {
                std::fill(request.output, request.output + pixelCount, ExpressionDeepEmptyPixel);
                return fail(ExpressionDeepRenderStatus::ResourceLimit, centered4Error.empty() ? "centered4 resource limit" : centered4Error);
            }
            if (centered4Status == Centered4ExperimentalStatus::UndefinedPixel) return fail(ExpressionDeepRenderStatus::UndefinedPixel, centered4Error.empty() ? "one or more centered4 pixel orbits are undefined" : centered4Error);
            if (centered4Status == Centered4ExperimentalStatus::Failed) {
                std::fill(request.output, request.output + pixelCount, ExpressionDeepEmptyPixel);
                return fail(ExpressionDeepRenderStatus::InternalError, centered4Error.empty() ? "centered4 experiment failed" : centered4Error);
            }
        }

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
