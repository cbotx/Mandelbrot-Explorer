#ifndef MANDEL_FORMULA_EXPRESSION_H
#define MANDEL_FORMULA_EXPRESSION_H

#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "formula_spec.h"

namespace formula {

using Complex = std::complex<double>;
class ExpressionOracle;
class ExpressionJit4;

struct ExpressionError {
    size_t position = 0;
    std::string message;
};

struct ExpressionContext {
    Complex z{};
    Complex c{};
    Complex z0{};
    std::array<Complex, 8> parameters{};
    int iteration = 0;
};

struct ExpressionDerivativeSeed {
    Complex z{};
    Complex c{};
    Complex z0{};
    std::array<Complex, 8> parameters{};
};

enum class ExpressionColoring : uint8_t {
    Raw,
    Smooth,
    Distance
};

class ExpressionProgram {
public:
    static constexpr size_t MAX_SOURCE = 4096;
    static constexpr size_t MAX_INSTRUCTIONS = 512;
    static constexpr size_t MAX_STACK = 128;
    static constexpr size_t MAX_PARSE_DEPTH = 64;

    enum class FastPath : uint8_t {
        None,
        IntegerPowerPlusC
    };

    bool compile(const std::string& source, ExpressionError* error = nullptr);
    bool specialize(const ExpressionContext& fixed,
                    FormulaParameter pixelParameter,
                    ExpressionProgram& output,
                    ExpressionError* error = nullptr) const;
    Complex evaluate(const ExpressionContext& context) const;
    Complex evaluate(const ExpressionContext& context, Complex* stack,
                     size_t capacity) const;
    bool evaluate4(const ExpressionContext* contexts, Complex* outputs) const;
    bool evaluate4Hybrid(const ExpressionContext* contexts,
                         Complex* outputs) const;
    bool evaluateWithDerivative(const ExpressionContext& context,
                                const ExpressionDerivativeSeed& seed,
                                Complex& value, Complex& derivative) const;

    bool valid() const { return _valid; }
    const std::string& source() const { return _source; }
    size_t instructionCount() const { return _code.size(); }
    size_t stackDepth() const { return _stackDepth; }
    FastPath fastPath() const { return _fastPath; }
    int fastIntegerPower() const { return _fastIntegerPower; }
    bool avx2Compatible() const { return _avx2Compatible; }
    bool batchCompatible() const { return _batchCompatible; }
    bool derivativeCompatible() const { return _derivativeCompatible; }

private:
    enum class Op : uint8_t {
        Constant,
        Z, C, Z0, Iteration, Parameter,
        Negate, Add, Subtract, Multiply, Divide, Power,
        Square,
        Sin, Cos, Tan, Sinh, Cosh, Tanh,
        Exp, Log, Log10, Sqrt,
        Abs, Norm, Arg, Conjugate, Real, Imaginary,
        MakeComplex, Polar
    };

    struct Instruction {
        Op op{};
        uint8_t argument = 0;
        Complex value{};
    };

    static constexpr uint8_t CONSTANT_FIXED_C = 1;

    std::vector<Instruction> _code;
    std::string _source;
    size_t _stackDepth = 0;
    FastPath _fastPath = FastPath::None;
    uint8_t _fastIntegerPower = 0;
    bool _avx2Compatible = false;
    bool _batchCompatible = false;
    bool _derivativeCompatible = false;
    bool _valid = false;

    bool analyze(ExpressionError* error);
    static int operandCount(Op op);
    static Complex evaluateUnary(Op op, Complex value);
    static Complex evaluateBinary(Op op, Complex left, Complex right);

    friend class ExpressionParser;
    friend class ExpressionOracle;
    friend class ExpressionJit4;
};

} // namespace formula

#endif
