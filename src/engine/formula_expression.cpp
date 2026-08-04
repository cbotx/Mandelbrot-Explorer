#include "formula_expression.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string_view>
#include <immintrin.h>

namespace formula {

class ExpressionParser {
public:
    ExpressionParser(ExpressionProgram& program, const std::string& source,
                     ExpressionError* error)
        : _program(program), _source(source), _error(error) {}

    bool parse() {
        if (_source.empty()) return fail(0, "expression is empty");
        if (_source.size() > ExpressionProgram::MAX_SOURCE)
            return fail(ExpressionProgram::MAX_SOURCE, "expression is too long");
        if (!parseExpression()) return false;
        skipSpace();
        if (_pos != _source.size()) return fail(_pos, "unexpected token");
        if (_stack != 1) return fail(_pos, "invalid expression stack");
        _program._stackDepth = _maxStack;
        return true;
    }

private:
    ExpressionProgram& _program;
    const std::string& _source;
    ExpressionError* _error;
    size_t _pos = 0;
    size_t _depth = 0;
    size_t _stack = 0;
    size_t _maxStack = 0;

    struct DepthGuard {
        ExpressionParser& parser;
        bool ok;
        explicit DepthGuard(ExpressionParser& p) : parser(p) {
            ok = ++parser._depth <= ExpressionProgram::MAX_PARSE_DEPTH;
        }
        ~DepthGuard() { --parser._depth; }
    };

    bool fail(size_t position, const char* message) {
        if (_error) {
            _error->position = position;
            _error->message = message;
        }
        return false;
    }

    void skipSpace() {
        while (_pos < _source.size() &&
               (_source[_pos] == ' ' || _source[_pos] == '\t' ||
                _source[_pos] == '\r' || _source[_pos] == '\n'))
            ++_pos;
    }

    bool consume(char token) {
        skipSpace();
        if (_pos >= _source.size() || _source[_pos] != token) return false;
        ++_pos;
        return true;
    }

    bool emit(ExpressionProgram::Op op, int pops, int pushes,
              uint8_t argument = 0, Complex value = {}) {
        if (_program._code.size() >= ExpressionProgram::MAX_INSTRUCTIONS)
            return fail(_pos, "expression has too many operations");
        if (_stack < (size_t)pops) return fail(_pos, "invalid operand count");
        _stack = _stack - pops + pushes;
        _maxStack = std::max(_maxStack, _stack);
        if (_maxStack > ExpressionProgram::MAX_STACK)
            return fail(_pos, "expression stack is too deep");
        _program._code.push_back({ op, argument, value });
        return true;
    }

    bool parseExpression() {
        DepthGuard guard(*this);
        if (!guard.ok) return fail(_pos, "expression nesting is too deep");
        return parseAddSubtract();
    }

    bool parseAddSubtract() {
        if (!parseMultiplyDivide()) return false;
        for (;;) {
            skipSpace();
            if (_pos >= _source.size()) return true;
            char op = _source[_pos];
            if (op != '+' && op != '-') return true;
            ++_pos;
            if (!parseMultiplyDivide()) return false;
            if (!emit(op == '+' ? ExpressionProgram::Op::Add
                                : ExpressionProgram::Op::Subtract, 2, 1))
                return false;
        }
    }

    bool parseMultiplyDivide() {
        if (!parseUnary()) return false;
        for (;;) {
            skipSpace();
            if (_pos >= _source.size()) return true;
            char op = _source[_pos];
            if (op != '*' && op != '/') return true;
            ++_pos;
            if (!parseUnary()) return false;
            if (!emit(op == '*' ? ExpressionProgram::Op::Multiply
                                : ExpressionProgram::Op::Divide, 2, 1))
                return false;
        }
    }

    bool parseUnary() {
        skipSpace();
        bool negate = false;
        while (_pos < _source.size() &&
               (_source[_pos] == '+' || _source[_pos] == '-')) {
            if (_source[_pos] == '-') negate = !negate;
            ++_pos;
            skipSpace();
        }
        if (!parsePower()) return false;
        return !negate || emit(ExpressionProgram::Op::Negate, 1, 1);
    }

    bool parsePower() {
        DepthGuard guard(*this);
        if (!guard.ok) return fail(_pos, "power nesting is too deep");
        if (!parsePrimary()) return false;
        skipSpace();
        if (!consume('^')) return true;
        if (!parseUnary()) return false; // right associative: z^2^3
        return emit(ExpressionProgram::Op::Power, 2, 1);
    }

    bool parsePrimary() {
        skipSpace();
        if (_pos >= _source.size()) return fail(_pos, "expected expression");
        if (consume('(')) {
            if (!parseExpression()) return false;
            if (!consume(')')) return fail(_pos, "expected ')'");
            return true;
        }
        unsigned char first = (unsigned char)_source[_pos];
        if ((first >= '0' && first <= '9') || first == '.')
            return parseNumber();
        if ((first >= 'A' && first <= 'Z') ||
            (first >= 'a' && first <= 'z') || first == '_')
            return parseIdentifier();
        return fail(_pos, "expected number, variable, function, or '('");
    }

    bool parseNumber() {
        const char* begin = _source.c_str() + _pos;
        char* end = nullptr;
        errno = 0;
        double value = std::strtod(begin, &end);
        if (end == begin || errno == ERANGE)
            return fail(_pos, "invalid numeric literal");
        _pos += (size_t)(end - begin);
        skipSpace();
        bool imaginary = _pos < _source.size() && _source[_pos] == 'i';
        if (imaginary) ++_pos;
        return emit(ExpressionProgram::Op::Constant, 0, 1, 0,
                    imaginary ? Complex{ 0.0, value } : Complex{ value, 0.0 });
    }

    std::string parseName() {
        skipSpace();
        size_t begin = _pos;
        while (_pos < _source.size()) {
            unsigned char ch = (unsigned char)_source[_pos];
            if (!((ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') ||
                  (ch >= '0' && ch <= '9') || ch == '_'))
                break;
            ++_pos;
        }
        std::string name = _source.substr(begin, _pos - begin);
        for (char& ch : name)
            if (ch >= 'A' && ch <= 'Z') ch = (char)(ch - 'A' + 'a');
        return name;
    }

    bool parseIdentifier() {
        size_t namePosition = _pos;
        std::string name = parseName();
        skipSpace();
        if (consume('(')) return parseFunction(name, namePosition);

        if (name == "z") return emit(ExpressionProgram::Op::Z, 0, 1);
        if (name == "c") return emit(ExpressionProgram::Op::C, 0, 1);
        if (name == "z0") return emit(ExpressionProgram::Op::Z0, 0, 1);
        if (name == "n") return emit(ExpressionProgram::Op::Iteration, 0, 1);
        if (name == "i") return emit(ExpressionProgram::Op::Constant, 0, 1, 0,
                                     Complex{ 0.0, 1.0 });
        if (name == "pi") return emit(ExpressionProgram::Op::Constant, 0, 1, 0,
                                      Complex{ 3.14159265358979323846, 0.0 });
        if (name == "e") return emit(ExpressionProgram::Op::Constant, 0, 1, 0,
                                     Complex{ 2.71828182845904523536, 0.0 });
        if (name.size() == 2 && name[0] == 'p' && name[1] >= '0' && name[1] <= '7')
            return emit(ExpressionProgram::Op::Parameter, 0, 1,
                        (uint8_t)(name[1] - '0'));
        return fail(namePosition, "unknown variable or constant");
    }

    bool parseFunction(const std::string& name, size_t namePosition) {
        struct Unary {
            std::string_view name;
            ExpressionProgram::Op op;
        };
        static constexpr Unary unary[] = {
            { "sqr", ExpressionProgram::Op::Square },
            { "sin", ExpressionProgram::Op::Sin },
            { "cos", ExpressionProgram::Op::Cos },
            { "tan", ExpressionProgram::Op::Tan },
            { "sinh", ExpressionProgram::Op::Sinh },
            { "cosh", ExpressionProgram::Op::Cosh },
            { "tanh", ExpressionProgram::Op::Tanh },
            { "exp", ExpressionProgram::Op::Exp },
            { "log", ExpressionProgram::Op::Log },
            { "log10", ExpressionProgram::Op::Log10 },
            { "sqrt", ExpressionProgram::Op::Sqrt },
            { "abs", ExpressionProgram::Op::Abs },
            { "norm", ExpressionProgram::Op::Norm },
            { "arg", ExpressionProgram::Op::Arg },
            { "conj", ExpressionProgram::Op::Conjugate },
            { "real", ExpressionProgram::Op::Real },
            { "re", ExpressionProgram::Op::Real },
            { "imag", ExpressionProgram::Op::Imaginary },
            { "im", ExpressionProgram::Op::Imaginary }
        };
        for (const Unary& function : unary) {
            if (name == function.name) {
                if (!parseExpression()) return false;
                if (consume(',')) return fail(_pos - 1, "function expects one argument");
                if (!consume(')')) return fail(_pos, "expected ')'");
                return emit(function.op, 1, 1);
            }
        }
        ExpressionProgram::Op op;
        if (name == "pow") op = ExpressionProgram::Op::Power;
        else if (name == "complex") op = ExpressionProgram::Op::MakeComplex;
        else if (name == "polar") op = ExpressionProgram::Op::Polar;
        else return fail(namePosition, "unknown function");

        if (!parseExpression()) return false;
        if (!consume(',')) return fail(_pos, "function expects two arguments");
        if (!parseExpression()) return false;
        if (consume(',')) return fail(_pos - 1, "function expects two arguments");
        if (!consume(')')) return fail(_pos, "expected ')'");
        return emit(op, 2, 1);
    }
};

bool ExpressionProgram::compile(const std::string& source, ExpressionError* error) {
    _valid = false;
    _source.clear();
    _code.clear();
    _stackDepth = 0;
    _fastPath = FastPath::None;
    _fastIntegerPower = 0;
    _avx2Compatible = false;
    _derivativeCompatible = false;
    if (error) *error = {};
    ExpressionParser parser(*this, source, error);
    if (!parser.parse()) {
        _code.clear();
        return false;
    }
    _source = source;
    const auto is = [&](size_t index, Op op) {
        return index < _code.size() && _code[index].op == op;
    };
    int degree = 0;
    if (_code.size() >= 5 && is(0, Op::Z) && is(1, Op::Z) &&
        is(2, Op::Multiply)) {
        degree = 2;
        size_t cursor = 3;
        while (degree < 16 && cursor + 1 < _code.size() &&
               is(cursor, Op::Z) && is(cursor + 1, Op::Multiply)) {
            ++degree;
            cursor += 2;
        }
        if (!(cursor + 2 == _code.size() && is(cursor, Op::C) &&
              is(cursor + 1, Op::Add)))
            degree = 0;
    } else if (_code.size() == 4 && is(0, Op::Z) && is(1, Op::Square) &&
               is(2, Op::C) && is(3, Op::Add)) {
        degree = 2;
    }
    if (degree >= 2) {
        _fastPath = FastPath::IntegerPowerPlusC;
        _fastIntegerPower = (uint8_t)degree;
    }
    _avx2Compatible = std::all_of(_code.begin(), _code.end(),
        [](const Instruction& instruction) {
            switch (instruction.op) {
            case Op::Constant:
            case Op::Z:
            case Op::C:
            case Op::Z0:
            case Op::Iteration:
            case Op::Parameter:
            case Op::Negate:
            case Op::Add:
            case Op::Subtract:
            case Op::Multiply:
            case Op::Square:
            case Op::Conjugate:
            case Op::Real:
            case Op::Imaginary:
            case Op::MakeComplex:
                return true;
            default:
                return false;
            }
        });
    _derivativeCompatible = std::all_of(_code.begin(), _code.end(),
        [](const Instruction& instruction) {
            switch (instruction.op) {
            case Op::Constant:
            case Op::Z:
            case Op::C:
            case Op::Z0:
            case Op::Iteration:
            case Op::Parameter:
            case Op::Negate:
            case Op::Add:
            case Op::Subtract:
            case Op::Multiply:
            case Op::Divide:
            case Op::Power:
            case Op::Square:
            case Op::Sin:
            case Op::Cos:
            case Op::Tan:
            case Op::Sinh:
            case Op::Cosh:
            case Op::Tanh:
            case Op::Exp:
            case Op::Log:
            case Op::Log10:
            case Op::Sqrt:
                return true;
            default:
                return false;
            }
        });
    _valid = true;
    return true;
}

Complex ExpressionProgram::evaluate(const ExpressionContext& context) const {
    std::array<Complex, MAX_STACK> stack;
    return evaluate(context, stack.data(), stack.size());
}

Complex ExpressionProgram::evaluate(const ExpressionContext& context,
                                    Complex* stack, size_t capacity) const {
    if (!_valid) return {
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()
    };
    if (!stack || capacity < _stackDepth) return {
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()
    };
    size_t top = 0;
    auto unary = [&](auto function) { stack[top - 1] = function(stack[top - 1]); };
    auto binary = [&](auto function) {
        Complex right = stack[--top];
        stack[top - 1] = function(stack[top - 1], right);
    };

    for (const Instruction& instruction : _code) {
        switch (instruction.op) {
        case Op::Constant: stack[top++] = instruction.value; break;
        case Op::Z: stack[top++] = context.z; break;
        case Op::C: stack[top++] = context.c; break;
        case Op::Z0: stack[top++] = context.z0; break;
        case Op::Iteration: stack[top++] = Complex{ (double)context.iteration, 0.0 }; break;
        case Op::Parameter: stack[top++] = context.parameters[instruction.argument]; break;
        case Op::Negate: unary([](Complex a) { return -a; }); break;
        case Op::Add: binary([](Complex a, Complex b) { return a + b; }); break;
        case Op::Subtract: binary([](Complex a, Complex b) { return a - b; }); break;
        case Op::Multiply: binary([](Complex a, Complex b) { return a * b; }); break;
        case Op::Divide: binary([](Complex a, Complex b) { return a / b; }); break;
        case Op::Power: binary([](Complex a, Complex b) { return std::pow(a, b); }); break;
        case Op::Square: unary([](Complex a) { return a * a; }); break;
        case Op::Sin: unary([](Complex a) { return std::sin(a); }); break;
        case Op::Cos: unary([](Complex a) { return std::cos(a); }); break;
        case Op::Tan: unary([](Complex a) { return std::tan(a); }); break;
        case Op::Sinh: unary([](Complex a) { return std::sinh(a); }); break;
        case Op::Cosh: unary([](Complex a) { return std::cosh(a); }); break;
        case Op::Tanh: unary([](Complex a) { return std::tanh(a); }); break;
        case Op::Exp: unary([](Complex a) { return std::exp(a); }); break;
        case Op::Log: unary([](Complex a) { return std::log(a); }); break;
        case Op::Log10: unary([](Complex a) { return std::log10(a); }); break;
        case Op::Sqrt: unary([](Complex a) { return std::sqrt(a); }); break;
        case Op::Abs: unary([](Complex a) { return Complex{ std::abs(a), 0.0 }; }); break;
        case Op::Norm: unary([](Complex a) { return Complex{ std::norm(a), 0.0 }; }); break;
        case Op::Arg: unary([](Complex a) { return Complex{ std::arg(a), 0.0 }; }); break;
        case Op::Conjugate: unary([](Complex a) { return std::conj(a); }); break;
        case Op::Real: unary([](Complex a) { return Complex{ a.real(), 0.0 }; }); break;
        case Op::Imaginary: unary([](Complex a) { return Complex{ a.imag(), 0.0 }; }); break;
        case Op::MakeComplex:
            binary([](Complex a, Complex b) { return Complex{ a.real(), b.real() }; });
            break;
        case Op::Polar:
            binary([](Complex a, Complex b) {
                const double nan = std::numeric_limits<double>::quiet_NaN();
                if (a.imag() != 0.0 || b.imag() != 0.0 ||
                    !std::isfinite(a.real()) || a.real() < 0.0 ||
                    !std::isfinite(b.real()))
                    return Complex{ nan, nan };
                return std::polar(a.real(), b.real());
            });
            break;
        }
    }
    return stack[0];
}

bool ExpressionProgram::evaluate4(const ExpressionContext* contexts,
                                  Complex* outputs) const {
    if (!_valid || !_avx2Compatible || !contexts || !outputs) return false;
    struct VecComplex { __m256d re, im; };
    std::array<VecComplex, MAX_STACK> stack;
    size_t top = 0;
    auto real = [&](auto getter) {
        return _mm256_set_pd(getter(contexts[3]).real(), getter(contexts[2]).real(),
                             getter(contexts[1]).real(), getter(contexts[0]).real());
    };
    auto imag = [&](auto getter) {
        return _mm256_set_pd(getter(contexts[3]).imag(), getter(contexts[2]).imag(),
                             getter(contexts[1]).imag(), getter(contexts[0]).imag());
    };
    for (const Instruction& instruction : _code) {
        switch (instruction.op) {
        case Op::Constant:
            stack[top++] = {
                _mm256_set1_pd(instruction.value.real()),
                _mm256_set1_pd(instruction.value.imag())
            };
            break;
        case Op::Z:
            stack[top++] = { real([](const ExpressionContext& c) { return c.z; }),
                             imag([](const ExpressionContext& c) { return c.z; }) };
            break;
        case Op::C:
            stack[top++] = { real([](const ExpressionContext& c) { return c.c; }),
                             imag([](const ExpressionContext& c) { return c.c; }) };
            break;
        case Op::Z0:
            stack[top++] = { real([](const ExpressionContext& c) { return c.z0; }),
                             imag([](const ExpressionContext& c) { return c.z0; }) };
            break;
        case Op::Iteration:
            stack[top++] = {
                _mm256_set_pd((double)contexts[3].iteration, (double)contexts[2].iteration,
                              (double)contexts[1].iteration, (double)contexts[0].iteration),
                _mm256_setzero_pd()
            };
            break;
        case Op::Parameter: {
            int p = instruction.argument;
            stack[top++] = {
                _mm256_set_pd(contexts[3].parameters[p].real(),
                              contexts[2].parameters[p].real(),
                              contexts[1].parameters[p].real(),
                              contexts[0].parameters[p].real()),
                _mm256_set_pd(contexts[3].parameters[p].imag(),
                              contexts[2].parameters[p].imag(),
                              contexts[1].parameters[p].imag(),
                              contexts[0].parameters[p].imag())
            };
            break;
        }
        case Op::Negate:
            stack[top - 1].re = _mm256_xor_pd(stack[top - 1].re, _mm256_set1_pd(-0.0));
            stack[top - 1].im = _mm256_xor_pd(stack[top - 1].im, _mm256_set1_pd(-0.0));
            break;
        case Op::Add:
        case Op::Subtract:
        case Op::Multiply: {
            VecComplex right = stack[--top];
            VecComplex& left = stack[top - 1];
            if (instruction.op == Op::Add) {
                left.re = _mm256_add_pd(left.re, right.re);
                left.im = _mm256_add_pd(left.im, right.im);
            } else if (instruction.op == Op::Subtract) {
                left.re = _mm256_sub_pd(left.re, right.re);
                left.im = _mm256_sub_pd(left.im, right.im);
            } else {
                __m256d re = _mm256_sub_pd(_mm256_mul_pd(left.re, right.re),
                                           _mm256_mul_pd(left.im, right.im));
                __m256d im = _mm256_add_pd(_mm256_mul_pd(left.re, right.im),
                                           _mm256_mul_pd(left.im, right.re));
                left.re = re; left.im = im;
            }
            break;
        }
        case Op::Square: {
            VecComplex& value = stack[top - 1];
            __m256d re = _mm256_sub_pd(_mm256_mul_pd(value.re, value.re),
                                       _mm256_mul_pd(value.im, value.im));
            __m256d im = _mm256_add_pd(_mm256_mul_pd(value.re, value.im),
                                       _mm256_mul_pd(value.im, value.re));
            value.re = re; value.im = im;
            break;
        }
        case Op::Conjugate:
            stack[top - 1].im = _mm256_xor_pd(stack[top - 1].im, _mm256_set1_pd(-0.0));
            break;
        case Op::Real:
            stack[top - 1].im = _mm256_setzero_pd();
            break;
        case Op::Imaginary:
            stack[top - 1].re = stack[top - 1].im;
            stack[top - 1].im = _mm256_setzero_pd();
            break;
        case Op::MakeComplex: {
            VecComplex right = stack[--top];
            stack[top - 1].im = right.re;
            break;
        }
        default:
            return false;
        }
    }
    alignas(32) double re[4], im[4];
    _mm256_store_pd(re, stack[0].re);
    _mm256_store_pd(im, stack[0].im);
    for (int lane = 0; lane < 4; ++lane) outputs[lane] = { re[lane], im[lane] };
    return true;
}

bool ExpressionProgram::evaluateWithDerivative(
        const ExpressionContext& context, const ExpressionDerivativeSeed& seed,
        Complex& value, Complex& derivative) const {
    value = derivative = {
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()
    };
    if (!_valid || !_derivativeCompatible) return false;
    struct Dual { Complex value, derivative; };
    std::array<Dual, MAX_STACK> stack;
    size_t top = 0;
    auto push = [&](Complex v, Complex d = {}) { stack[top++] = { v, d }; };
    auto unary = [&](auto function) {
        Dual input = stack[top - 1];
        stack[top - 1] = function(input);
    };
    auto binary = [&](auto function) {
        Dual right = stack[--top];
        Dual left = stack[top - 1];
        stack[top - 1] = function(left, right);
    };
    bool domain = true;
    const double ln10 = std::log(10.0);

    for (const Instruction& instruction : _code) {
        switch (instruction.op) {
        case Op::Constant: push(instruction.value); break;
        case Op::Z: push(context.z, seed.z); break;
        case Op::C: push(context.c, seed.c); break;
        case Op::Z0: push(context.z0, seed.z0); break;
        case Op::Iteration: push({ (double)context.iteration, 0.0 }); break;
        case Op::Parameter:
            push(context.parameters[instruction.argument],
                 seed.parameters[instruction.argument]);
            break;
        case Op::Negate:
            unary([](Dual a) { return Dual{ -a.value, -a.derivative }; });
            break;
        case Op::Add:
            binary([](Dual a, Dual b) {
                return Dual{ a.value + b.value, a.derivative + b.derivative };
            });
            break;
        case Op::Subtract:
            binary([](Dual a, Dual b) {
                return Dual{ a.value - b.value, a.derivative - b.derivative };
            });
            break;
        case Op::Multiply:
            binary([](Dual a, Dual b) {
                return Dual{ a.value * b.value,
                             a.derivative * b.value + a.value * b.derivative };
            });
            break;
        case Op::Divide:
            binary([&](Dual a, Dual b) {
                if (b.value == Complex{}) {
                    domain = false;
                    return Dual{};
                }
                Complex result = a.value / b.value;
                return Dual{ result,
                             (a.derivative - result * b.derivative) / b.value };
            });
            break;
        case Op::Power:
            binary([&](Dual a, Dual b) {
                bool constantInteger =
                    b.derivative == Complex{} && b.value.imag() == 0.0 &&
                    std::isfinite(b.value.real()) &&
                    std::trunc(b.value.real()) == b.value.real();
                if (constantInteger && a.value == Complex{} &&
                    b.value.real() >= 0.0) {
                    if (b.value.real() == 0.0)
                        return Dual{ Complex{ 1.0, 0.0 }, {} };
                    int exponent = b.value.real() <= 1.0 ? 1 : 2;
                    return Dual{ {},
                        exponent == 1 ? a.derivative : Complex{} };
                }
                if (a.value == Complex{}) {
                    domain = false;
                    return Dual{};
                }
                Complex result = std::pow(a.value, b.value);
                Complex d;
                if (constantInteger) {
                    d = b.value * std::pow(a.value, b.value - Complex{ 1.0, 0.0 }) *
                        a.derivative;
                } else {
                    if (a.value.imag() == 0.0 && a.value.real() < 0.0) {
                        domain = false;
                        return Dual{};
                    }
                    d = result * (b.derivative * std::log(a.value) +
                                  b.value * a.derivative / a.value);
                }
                return Dual{ result, d };
            });
            break;
        case Op::Square:
            unary([](Dual a) {
                return Dual{ a.value * a.value,
                             a.derivative * a.value + a.value * a.derivative };
            });
            break;
        case Op::Sin:
            unary([](Dual a) {
                return Dual{ std::sin(a.value), std::cos(a.value) * a.derivative };
            });
            break;
        case Op::Cos:
            unary([](Dual a) {
                return Dual{ std::cos(a.value), -std::sin(a.value) * a.derivative };
            });
            break;
        case Op::Tan:
            unary([](Dual a) {
                Complex result = std::tan(a.value);
                Complex factor;
                if (std::abs(a.value.imag()) < 300.0) {
                    Complex c = std::cos(a.value);
                    factor = 1.0 / (c * c);
                } else {
                    double q = std::exp(-2.0 * std::abs(a.value.imag()));
                    double sign = std::signbit(a.value.imag()) ? -1.0 : 1.0;
                    Complex scaled{
                        std::cos(a.value.real()) * (1.0 + q),
                        -sign * std::sin(a.value.real()) * (1.0 - q)
                    };
                    factor = q == 0.0 ? Complex{} : 4.0 * q / (scaled * scaled);
                }
                return Dual{ result, factor * a.derivative };
            });
            break;
        case Op::Sinh:
            unary([](Dual a) {
                return Dual{ std::sinh(a.value), std::cosh(a.value) * a.derivative };
            });
            break;
        case Op::Cosh:
            unary([](Dual a) {
                return Dual{ std::cosh(a.value), std::sinh(a.value) * a.derivative };
            });
            break;
        case Op::Tanh:
            unary([](Dual a) {
                Complex result = std::tanh(a.value);
                Complex factor;
                if (std::abs(a.value.real()) < 300.0) {
                    Complex c = std::cosh(a.value);
                    factor = 1.0 / (c * c);
                } else {
                    double q = std::exp(-2.0 * std::abs(a.value.real()));
                    double sign = std::signbit(a.value.real()) ? -1.0 : 1.0;
                    Complex scaled{
                        std::cos(a.value.imag()) * (1.0 + q),
                        sign * std::sin(a.value.imag()) * (1.0 - q)
                    };
                    factor = q == 0.0 ? Complex{} : 4.0 * q / (scaled * scaled);
                }
                return Dual{ result, factor * a.derivative };
            });
            break;
        case Op::Exp:
            unary([](Dual a) {
                Complex result = std::exp(a.value);
                return Dual{ result, result * a.derivative };
            });
            break;
        case Op::Log:
            unary([&](Dual a) {
                if (a.value == Complex{} ||
                    (a.value.imag() == 0.0 && a.value.real() < 0.0)) {
                    domain = false; return Dual{};
                }
                return Dual{ std::log(a.value), a.derivative / a.value };
            });
            break;
        case Op::Log10:
            unary([&](Dual a) {
                if (a.value == Complex{} ||
                    (a.value.imag() == 0.0 && a.value.real() < 0.0)) {
                    domain = false; return Dual{};
                }
                return Dual{ std::log10(a.value),
                             a.derivative / (a.value * ln10) };
            });
            break;
        case Op::Sqrt:
            unary([&](Dual a) {
                if (a.value.imag() == 0.0 && a.value.real() < 0.0) {
                    domain = false; return Dual{};
                }
                Complex result = std::sqrt(a.value);
                if (result == Complex{}) { domain = false; return Dual{}; }
                return Dual{ result, a.derivative / (2.0 * result) };
            });
            break;
        default:
            return false;
        }
        if (!domain) return false;
    }
    value = stack[0].value;
    derivative = stack[0].derivative;
    return std::isfinite(value.real()) && std::isfinite(value.imag()) &&
           std::isfinite(derivative.real()) && std::isfinite(derivative.imag());
}

} // namespace formula
