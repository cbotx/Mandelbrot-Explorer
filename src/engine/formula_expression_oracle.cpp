#include "formula_expression_oracle.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace formula {

namespace {

constexpr mpfr_rnd_t RND = MPFR_RNDN;

struct Scratch {
    explicit Scratch(mpfr_prec_t precision) {
        for (mpfr_ptr value : values) mpfr_init2(value, precision);
    }
    ~Scratch() {
        for (mpfr_ptr value : values) mpfr_clear(value);
    }
    mpfr_t values[12];
};

void setNan(MpfrComplex& value) {
    mpfr_set_nan(value.re);
    mpfr_set_nan(value.im);
}

void add(MpfrComplex& out, const MpfrComplex& a, const MpfrComplex& b) {
    mpfr_add(out.re, a.re, b.re, RND);
    mpfr_add(out.im, a.im, b.im, RND);
}

void sub(MpfrComplex& out, const MpfrComplex& a, const MpfrComplex& b) {
    mpfr_sub(out.re, a.re, b.re, RND);
    mpfr_sub(out.im, a.im, b.im, RND);
}

bool scaleExponent(mpfr_ptr out, mpfr_srcptr in, mpfr_exp_t exponent) {
    if (mpfr_zero_p(in)) {
        mpfr_set_zero(out, mpfr_signbit(in) ? -1 : 1);
        return true;
    }
    mpfr_exp_t target = mpfr_get_exp(in) - exponent;
    if (target < mpfr_get_emin()) {
        mpfr_set_zero(out, mpfr_signbit(in) ? -1 : 1);
        return true;
    }
    if (target > mpfr_get_emax()) return false;
    mpfr_set(out, in, RND);
    return mpfr_set_exp(out, target) == 0;
}

bool shiftExponent(mpfr_ptr value, mpfr_exp_t exponent) {
    if (mpfr_zero_p(value)) return true;
    mpfr_exp_t target = mpfr_get_exp(value) + exponent;
    if (target < mpfr_get_emin()) {
        mpfr_set_zero(value, mpfr_signbit(value) ? -1 : 1);
        return true;
    }
    if (target > mpfr_get_emax()) return false;
    return mpfr_set_exp(value, target) == 0;
}

void mul(MpfrComplex& out, const MpfrComplex& a, const MpfrComplex& b, Scratch& s) {
    mpfr_mul(s.values[0], a.re, b.re, RND);
    mpfr_mul(s.values[1], a.im, b.im, RND);
    mpfr_sub(s.values[2], s.values[0], s.values[1], RND);
    mpfr_mul(s.values[0], a.re, b.im, RND);
    mpfr_mul(s.values[1], a.im, b.re, RND);
    mpfr_add(s.values[3], s.values[0], s.values[1], RND);
    mpfr_set(out.re, s.values[2], RND);
    mpfr_set(out.im, s.values[3], RND);
}

bool divide(MpfrComplex& out, const MpfrComplex& a, const MpfrComplex& b, Scratch& s) {
    if (mpfr_zero_p(b.re) && mpfr_zero_p(b.im)) {
        setNan(out);
        return false;
    }
    mpfr_exp_t ea = mpfr_zero_p(a.re) ? mpfr_get_exp(a.im) : mpfr_get_exp(a.re);
    if (!mpfr_zero_p(a.im)) ea = std::max(ea, mpfr_get_exp(a.im));
    mpfr_exp_t eb = mpfr_zero_p(b.re) ? mpfr_get_exp(b.im) : mpfr_get_exp(b.re);
    if (!mpfr_zero_p(b.im)) eb = std::max(eb, mpfr_get_exp(b.im));
    ++ea; ++eb; // each normalized vector has max component below 0.5
    if (!scaleExponent(s.values[0], a.re, ea) ||
        !scaleExponent(s.values[1], a.im, ea) ||
        !scaleExponent(s.values[2], b.re, eb) ||
        !scaleExponent(s.values[3], b.im, eb)) {
        setNan(out);
        return false;
    }
    mpfr_sqr(s.values[4], s.values[2], RND);
    mpfr_sqr(s.values[5], s.values[3], RND);
    mpfr_add(s.values[4], s.values[4], s.values[5], RND); // denominator
    mpfr_mul(s.values[5], s.values[0], s.values[2], RND);
    mpfr_mul(s.values[6], s.values[1], s.values[3], RND);
    mpfr_add(s.values[5], s.values[5], s.values[6], RND);
    mpfr_div(s.values[5], s.values[5], s.values[4], RND);
    mpfr_mul(s.values[6], s.values[1], s.values[2], RND);
    mpfr_mul(s.values[7], s.values[0], s.values[3], RND);
    mpfr_sub(s.values[6], s.values[6], s.values[7], RND);
    mpfr_div(s.values[6], s.values[6], s.values[4], RND);
    mpfr_set(out.re, s.values[5], RND);
    mpfr_set(out.im, s.values[6], RND);
    mpfr_exp_t quotientExponent = ea - eb;
    if (!shiftExponent(out.re, quotientExponent) ||
        !shiftExponent(out.im, quotientExponent)) {
        setNan(out);
        return false;
    }
    if (!mpfr_number_p(out.re) || !mpfr_number_p(out.im)) {
        setNan(out);
        return false;
    }
    return true;
}

void absolute(MpfrComplex& out, const MpfrComplex& a) {
    mpfr_hypot(out.re, a.re, a.im, RND);
    mpfr_set_zero(out.im, 0);
}

void norm(MpfrComplex& out, const MpfrComplex& a, Scratch& s) {
    mpfr_sqr(s.values[0], a.re, RND);
    mpfr_sqr(s.values[1], a.im, RND);
    mpfr_add(out.re, s.values[0], s.values[1], RND);
    mpfr_set_zero(out.im, 0);
}

void argument(MpfrComplex& out, const MpfrComplex& a) {
    mpfr_atan2(out.re, a.im, a.re, RND);
    mpfr_set_zero(out.im, 0);
}

void exponential(MpfrComplex& out, const MpfrComplex& a, Scratch& s) {
    mpfr_exp(s.values[0], a.re, RND);
    mpfr_cos(s.values[1], a.im, RND);
    mpfr_sin(s.values[2], a.im, RND);
    mpfr_mul(out.re, s.values[0], s.values[1], RND);
    mpfr_mul(out.im, s.values[0], s.values[2], RND);
}

bool logarithm(MpfrComplex& out, const MpfrComplex& a, Scratch& s) {
    mpfr_hypot(s.values[0], a.re, a.im, RND);
    if (mpfr_zero_p(s.values[0])) {
        setNan(out);
        return false;
    }
    mpfr_log(out.re, s.values[0], RND);
    mpfr_atan2(out.im, a.im, a.re, RND);
    return true;
}

void sine(MpfrComplex& out, const MpfrComplex& a, Scratch& s) {
    mpfr_sin(s.values[0], a.re, RND);
    mpfr_cosh(s.values[1], a.im, RND);
    mpfr_mul(out.re, s.values[0], s.values[1], RND);
    mpfr_cos(s.values[0], a.re, RND);
    mpfr_sinh(s.values[1], a.im, RND);
    mpfr_mul(out.im, s.values[0], s.values[1], RND);
}

void cosine(MpfrComplex& out, const MpfrComplex& a, Scratch& s) {
    mpfr_cos(s.values[0], a.re, RND);
    mpfr_cosh(s.values[1], a.im, RND);
    mpfr_mul(out.re, s.values[0], s.values[1], RND);
    mpfr_sin(s.values[0], a.re, RND);
    mpfr_sinh(s.values[1], a.im, RND);
    mpfr_mul(out.im, s.values[0], s.values[1], RND);
    mpfr_neg(out.im, out.im, RND);
}

void hyperbolicSine(MpfrComplex& out, const MpfrComplex& a, Scratch& s) {
    mpfr_sinh(s.values[0], a.re, RND);
    mpfr_cos(s.values[1], a.im, RND);
    mpfr_mul(out.re, s.values[0], s.values[1], RND);
    mpfr_cosh(s.values[0], a.re, RND);
    mpfr_sin(s.values[1], a.im, RND);
    mpfr_mul(out.im, s.values[0], s.values[1], RND);
}

void hyperbolicCosine(MpfrComplex& out, const MpfrComplex& a, Scratch& s) {
    mpfr_cosh(s.values[0], a.re, RND);
    mpfr_cos(s.values[1], a.im, RND);
    mpfr_mul(out.re, s.values[0], s.values[1], RND);
    mpfr_sinh(s.values[0], a.re, RND);
    mpfr_sin(s.values[1], a.im, RND);
    mpfr_mul(out.im, s.values[0], s.values[1], RND);
}

bool squareRoot(MpfrComplex& out, const MpfrComplex& a, Scratch& s) {
    if (mpfr_nan_p(a.re) || mpfr_nan_p(a.im)) {
        setNan(out);
        return false;
    }
    if (mpfr_zero_p(a.re) && mpfr_zero_p(a.im)) {
        out.set(0.0, 0.0);
        return true;
    }
    mpfr_exp_t e = mpfr_zero_p(a.re) ? mpfr_get_exp(a.im) : mpfr_get_exp(a.re);
    if (!mpfr_zero_p(a.im)) e = std::max(e, mpfr_get_exp(a.im));
    if (e & 1) --e; // even exponent; scaled inputs stay below 2
    if (!scaleExponent(s.values[10], a.re, e) ||
        !scaleExponent(s.values[11], a.im, e)) {
        setNan(out);
        return false;
    }
    mpfr_hypot(s.values[0], s.values[10], s.values[11], RND);
    if (mpfr_sgn(s.values[10]) >= 0) {
        // Real part is the large component; derive the small imaginary component
        // from y/(2*real) instead of subtracting two nearly equal magnitudes.
        mpfr_add(s.values[1], s.values[0], s.values[10], RND);
        mpfr_div_2ui(s.values[1], s.values[1], 1, RND);
        mpfr_sqrt(s.values[2], s.values[1], RND);
        mpfr_set(out.re, s.values[2], RND);
        mpfr_set_zero(out.im, mpfr_signbit(a.im) ? -1 : 1);
    } else {
        // Imaginary part is the large component; derive real=|y|/(2*|imag|).
        mpfr_sub(s.values[1], s.values[0], s.values[10], RND);
        mpfr_div_2ui(s.values[1], s.values[1], 1, RND);
        mpfr_sqrt(s.values[2], s.values[1], RND);
        mpfr_copysign(out.im, s.values[2], s.values[11], RND);
        mpfr_set_zero(out.re, 0);
    }
    mpfr_exp_t half = e / 2;
    if (!shiftExponent(out.re, half) || !shiftExponent(out.im, half)) {
        setNan(out);
        return false;
    }
    if (mpfr_sgn(a.re) >= 0 && !mpfr_zero_p(out.re)) {
        // Recover the small component at its FINAL exponent.
        mpfr_mul_2ui(s.values[3], out.re, 1, RND);
        mpfr_div(out.im, a.im, s.values[3], RND);
    } else if (mpfr_sgn(a.re) < 0 && !mpfr_zero_p(out.im)) {
        mpfr_abs(s.values[3], a.im, RND);
        mpfr_abs(s.values[4], out.im, RND);
        mpfr_mul_2ui(s.values[4], s.values[4], 1, RND);
        mpfr_div(out.re, s.values[3], s.values[4], RND);
    }
    return mpfr_number_p(out.re) && mpfr_number_p(out.im);
}

bool tangent(MpfrComplex& out, const MpfrComplex& a, Scratch& s) {
    if (mpfr_cmpabs_ui(a.im, 16) < 0) {
        MpfrComplex sn(out.precision()), cs(out.precision());
        sine(sn, a, s);
        cosine(cs, a, s);
        return divide(out, sn, cs, s);
    }
    // Multiply the standard formula by 2*exp(-2|y|):
    // tan(x+iy) = [2q sin(2x) + i sign(y)(1-q^2)] /
    //             [1+q^2+2q cos(2x)], q=exp(-2|y|).
    mpfr_abs(s.values[0], a.im, RND);
    mpfr_mul_si(s.values[0], s.values[0], -2, RND);
    mpfr_exp(s.values[1], s.values[0], RND);             // q
    mpfr_sqr(s.values[2], s.values[1], RND);             // q^2
    mpfr_sin_cos(s.values[3], s.values[4], a.re, RND);   // sin x, cos x
    mpfr_mul(s.values[5], s.values[3], s.values[4], RND);
    mpfr_mul_2ui(s.values[5], s.values[5], 1, RND);       // sin(2x)
    mpfr_sqr(s.values[6], s.values[4], RND);
    mpfr_sqr(s.values[7], s.values[3], RND);
    mpfr_sub(s.values[6], s.values[6], s.values[7], RND); // cos(2x)
    mpfr_mul(s.values[7], s.values[6], s.values[1], RND);
    mpfr_mul_2ui(s.values[7], s.values[7], 1, RND);
    mpfr_add_ui(s.values[8], s.values[2], 1, RND);
    mpfr_add(s.values[8], s.values[8], s.values[7], RND); // denominator
    if (mpfr_zero_p(s.values[8])) { setNan(out); return false; }
    mpfr_mul(s.values[7], s.values[5], s.values[1], RND);
    mpfr_mul_2ui(s.values[7], s.values[7], 1, RND);
    mpfr_div(out.re, s.values[7], s.values[8], RND);
    mpfr_ui_sub(s.values[7], 1, s.values[2], RND);
    mpfr_copysign(s.values[7], s.values[7], a.im, RND);
    mpfr_div(out.im, s.values[7], s.values[8], RND);
    return mpfr_number_p(out.re) && mpfr_number_p(out.im);
}

bool hyperbolicTangent(MpfrComplex& out, const MpfrComplex& a, Scratch& s) {
    if (mpfr_cmpabs_ui(a.re, 16) < 0) {
        MpfrComplex sh(out.precision()), ch(out.precision());
        hyperbolicSine(sh, a, s);
        hyperbolicCosine(ch, a, s);
        return divide(out, sh, ch, s);
    }
    // tanh(x+iy) is the same scaled identity with x/y exchanged.
    mpfr_abs(s.values[0], a.re, RND);
    mpfr_mul_si(s.values[0], s.values[0], -2, RND);
    mpfr_exp(s.values[1], s.values[0], RND);              // q
    mpfr_sqr(s.values[2], s.values[1], RND);
    mpfr_sin_cos(s.values[3], s.values[4], a.im, RND);    // sin y, cos y
    mpfr_mul(s.values[5], s.values[3], s.values[4], RND);
    mpfr_mul_2ui(s.values[5], s.values[5], 1, RND);        // sin(2y)
    mpfr_sqr(s.values[6], s.values[4], RND);
    mpfr_sqr(s.values[7], s.values[3], RND);
    mpfr_sub(s.values[6], s.values[6], s.values[7], RND);  // cos(2y)
    mpfr_mul(s.values[7], s.values[6], s.values[1], RND);
    mpfr_mul_2ui(s.values[7], s.values[7], 1, RND);
    mpfr_add_ui(s.values[8], s.values[2], 1, RND);
    mpfr_add(s.values[8], s.values[8], s.values[7], RND);
    if (mpfr_zero_p(s.values[8])) { setNan(out); return false; }
    mpfr_ui_sub(s.values[7], 1, s.values[2], RND);
    mpfr_copysign(s.values[7], s.values[7], a.re, RND);
    mpfr_div(out.re, s.values[7], s.values[8], RND);
    mpfr_mul(s.values[7], s.values[5], s.values[1], RND);
    mpfr_mul_2ui(s.values[7], s.values[7], 1, RND);
    mpfr_div(out.im, s.values[7], s.values[8], RND);
    return mpfr_number_p(out.re) && mpfr_number_p(out.im);
}

bool power(MpfrComplex& out, const MpfrComplex& a, const MpfrComplex& b, Scratch& s) {
    if (mpfr_zero_p(a.re) && mpfr_zero_p(a.im)) {
        if (mpfr_zero_p(b.im) && mpfr_sgn(b.re) > 0) {
            out.set(0.0, 0.0);
            return true;
        }
        if (mpfr_zero_p(b.re) && mpfr_zero_p(b.im)) {
            out.set(1.0, 0.0);
            return true;
        }
        setNan(out);
        return false;
    }
    MpfrComplex logA(out.precision()), product(out.precision());
    if (!logarithm(logA, a, s)) return false;
    mul(product, b, logA, s);
    exponential(out, product, s);
    return true;
}

} // namespace

MpfrComplex::MpfrComplex(mpfr_prec_t precision) {
    mpfr_init2(re, precision);
    mpfr_init2(im, precision);
}

MpfrComplex::MpfrComplex(const MpfrComplex& other)
    : MpfrComplex(other.precision()) {
    set(other);
}

MpfrComplex& MpfrComplex::operator=(const MpfrComplex& other) {
    if (this != &other) {
        mpfr_set_prec(re, other.precision());
        mpfr_set_prec(im, other.precision());
        set(other);
    }
    return *this;
}

MpfrComplex::~MpfrComplex() {
    mpfr_clear(re);
    mpfr_clear(im);
}

void MpfrComplex::set(double real, double imaginary) {
    mpfr_set_d(re, real, RND);
    mpfr_set_d(im, imaginary, RND);
}

bool MpfrComplex::set(const std::string& real, const std::string& imaginary) {
    MpfrComplex parsed(precision());
    if (mpfr_set_str(parsed.re, real.c_str(), 10, RND) != 0 ||
        mpfr_set_str(parsed.im, imaginary.c_str(), 10, RND) != 0)
        return false;
    set(parsed);
    return true;
}

void MpfrComplex::set(const MpfrComplex& other) {
    mpfr_set(re, other.re, RND);
    mpfr_set(im, other.im, RND);
}

std::complex<double> MpfrComplex::toDouble() const {
    return { mpfr_get_d(re, RND), mpfr_get_d(im, RND) };
}

mpfr_prec_t MpfrComplex::precision() const {
    return std::min(mpfr_get_prec(re), mpfr_get_prec(im));
}

ExpressionOracleContext::ExpressionOracleContext(mpfr_prec_t precision)
    : z(precision), c(precision), z0(precision),
      parameters{ MpfrComplex(precision), MpfrComplex(precision),
                  MpfrComplex(precision), MpfrComplex(precision),
                  MpfrComplex(precision), MpfrComplex(precision),
                  MpfrComplex(precision), MpfrComplex(precision) } {}

bool ExpressionOracle::evaluate(const ExpressionProgram& program,
                                const ExpressionOracleContext& context,
                                MpfrComplex& output,
                                std::string* error) {
    if (error) error->clear();
    if (!program.valid()) {
        if (error) *error = "program is not compiled";
        setNan(output);
        return false;
    }
    const mpfr_prec_t precision = output.precision();
    std::vector<MpfrComplex> stack;
    stack.reserve(program.stackDepth());
    for (size_t i = 0; i < program.stackDepth(); ++i)
        stack.emplace_back(precision);
    Scratch scratch(precision);
    size_t top = 0;
    bool exactDomain = true;

    auto unary = [&](auto function) {
        function(stack[top - 1], stack[top - 1], scratch);
    };
    auto binary = [&](auto function) {
        MpfrComplex right(stack[top - 1]);
        --top;
        function(stack[top - 1], stack[top - 1], right, scratch);
    };

    for (const ExpressionProgram::Instruction& instruction : program._code) {
        switch (instruction.op) {
        case ExpressionProgram::Op::Constant:
            stack[top++].set(instruction.value.real(), instruction.value.imag()); break;
        case ExpressionProgram::Op::Z: stack[top++].set(context.z); break;
        case ExpressionProgram::Op::C: stack[top++].set(context.c); break;
        case ExpressionProgram::Op::Z0: stack[top++].set(context.z0); break;
        case ExpressionProgram::Op::Iteration:
            mpfr_set_si(stack[top].re, context.iteration, RND);
            mpfr_set_zero(stack[top].im, 0); ++top; break;
        case ExpressionProgram::Op::Parameter:
            stack[top++].set(context.parameters[instruction.argument]); break;
        case ExpressionProgram::Op::Negate:
            mpfr_neg(stack[top - 1].re, stack[top - 1].re, RND);
            mpfr_neg(stack[top - 1].im, stack[top - 1].im, RND); break;
        case ExpressionProgram::Op::Add:
            binary([](MpfrComplex& o, const MpfrComplex& a, const MpfrComplex& b, Scratch&) { add(o, a, b); }); break;
        case ExpressionProgram::Op::Subtract:
            binary([](MpfrComplex& o, const MpfrComplex& a, const MpfrComplex& b, Scratch&) { sub(o, a, b); }); break;
        case ExpressionProgram::Op::Multiply:
            binary([](MpfrComplex& o, const MpfrComplex& a, const MpfrComplex& b, Scratch& s) { mul(o, a, b, s); }); break;
        case ExpressionProgram::Op::Divide:
            binary([&](MpfrComplex& o, const MpfrComplex& a, const MpfrComplex& b, Scratch& s) {
                exactDomain = divide(o, a, b, s) && exactDomain;
            }); break;
        case ExpressionProgram::Op::Power:
            binary([&](MpfrComplex& o, const MpfrComplex& a, const MpfrComplex& b, Scratch& s) {
                exactDomain = power(o, a, b, s) && exactDomain;
            }); break;
        case ExpressionProgram::Op::Square:
            unary([](MpfrComplex& o, const MpfrComplex& a, Scratch& s) {
                MpfrComplex copy(a); mul(o, copy, copy, s);
            }); break;
        case ExpressionProgram::Op::Sin:
            unary([](MpfrComplex& o, const MpfrComplex& a, Scratch& s) { MpfrComplex copy(a); sine(o, copy, s); }); break;
        case ExpressionProgram::Op::Cos:
            unary([](MpfrComplex& o, const MpfrComplex& a, Scratch& s) { MpfrComplex copy(a); cosine(o, copy, s); }); break;
        case ExpressionProgram::Op::Tan:
            unary([&](MpfrComplex& o, const MpfrComplex& a, Scratch& s) {
                MpfrComplex copy(a);
                exactDomain = tangent(o, copy, s) && exactDomain;
            }); break;
        case ExpressionProgram::Op::Sinh:
            unary([](MpfrComplex& o, const MpfrComplex& a, Scratch& s) { MpfrComplex copy(a); hyperbolicSine(o, copy, s); }); break;
        case ExpressionProgram::Op::Cosh:
            unary([](MpfrComplex& o, const MpfrComplex& a, Scratch& s) { MpfrComplex copy(a); hyperbolicCosine(o, copy, s); }); break;
        case ExpressionProgram::Op::Tanh:
            unary([&](MpfrComplex& o, const MpfrComplex& a, Scratch& s) {
                MpfrComplex copy(a);
                exactDomain = hyperbolicTangent(o, copy, s) && exactDomain;
            }); break;
        case ExpressionProgram::Op::Exp:
            unary([](MpfrComplex& o, const MpfrComplex& a, Scratch& s) { MpfrComplex copy(a); exponential(o, copy, s); }); break;
        case ExpressionProgram::Op::Log:
            unary([&](MpfrComplex& o, const MpfrComplex& a, Scratch& s) {
                MpfrComplex copy(a); exactDomain = logarithm(o, copy, s) && exactDomain;
            }); break;
        case ExpressionProgram::Op::Log10:
            unary([&](MpfrComplex& o, const MpfrComplex& a, Scratch& s) {
                MpfrComplex copy(a);
                bool valid = logarithm(o, copy, s);
                exactDomain = valid && exactDomain;
                if (valid) {
                    mpfr_const_log2(s.values[0], RND);
                    mpfr_log_ui(s.values[1], 5, RND);
                    mpfr_add(s.values[0], s.values[0], s.values[1], RND);
                    mpfr_div(o.re, o.re, s.values[0], RND);
                    mpfr_div(o.im, o.im, s.values[0], RND);
                }
            }); break;
        case ExpressionProgram::Op::Sqrt:
            unary([&](MpfrComplex& o, const MpfrComplex& a, Scratch& s) {
                MpfrComplex copy(a); exactDomain = squareRoot(o, copy, s) && exactDomain;
            }); break;
        case ExpressionProgram::Op::Abs:
            unary([](MpfrComplex& o, const MpfrComplex& a, Scratch&) { MpfrComplex copy(a); absolute(o, copy); }); break;
        case ExpressionProgram::Op::Norm:
            unary([](MpfrComplex& o, const MpfrComplex& a, Scratch& s) { MpfrComplex copy(a); norm(o, copy, s); }); break;
        case ExpressionProgram::Op::Arg:
            unary([](MpfrComplex& o, const MpfrComplex& a, Scratch&) { MpfrComplex copy(a); argument(o, copy); }); break;
        case ExpressionProgram::Op::Conjugate:
            mpfr_neg(stack[top - 1].im, stack[top - 1].im, RND); break;
        case ExpressionProgram::Op::Real:
            mpfr_set_zero(stack[top - 1].im, 0); break;
        case ExpressionProgram::Op::Imaginary:
            mpfr_set(stack[top - 1].re, stack[top - 1].im, RND);
            mpfr_set_zero(stack[top - 1].im, 0); break;
        case ExpressionProgram::Op::MakeComplex: {
            MpfrComplex right(stack[top - 1]); --top;
            mpfr_set(stack[top - 1].im, right.re, RND);
            break;
        }
        case ExpressionProgram::Op::Polar: {
            MpfrComplex angle(stack[top - 1]); --top;
            MpfrComplex& radius = stack[top - 1];
            if (!mpfr_zero_p(radius.im) || !mpfr_zero_p(angle.im) ||
                mpfr_nan_p(radius.re) || mpfr_inf_p(radius.re) ||
                mpfr_sgn(radius.re) < 0 || mpfr_nan_p(angle.re) ||
                mpfr_inf_p(angle.re)) {
                setNan(radius); exactDomain = false;
            } else {
                mpfr_cos(scratch.values[0], angle.re, RND);
                mpfr_sin(scratch.values[1], angle.re, RND);
                mpfr_mul(radius.im, radius.re, scratch.values[1], RND);
                mpfr_mul(radius.re, radius.re, scratch.values[0], RND);
            }
            break;
        }
        }
    }
    output.set(stack[0]);
    if (!exactDomain && error) *error = "expression evaluated outside a function domain";
    return exactDomain;
}

} // namespace formula
