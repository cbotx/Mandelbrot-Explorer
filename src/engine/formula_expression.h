#ifndef MANDEL_FORMULA_EXPRESSION_H
#define MANDEL_FORMULA_EXPRESSION_H

#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace formula {

using Complex = std::complex<double>;
class ExpressionOracle;

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

class ExpressionProgram {
public:
    static constexpr size_t MAX_SOURCE = 4096;
    static constexpr size_t MAX_INSTRUCTIONS = 512;
    static constexpr size_t MAX_STACK = 128;
    static constexpr size_t MAX_PARSE_DEPTH = 64;

    enum class FastPath : uint8_t {
        None,
        QuadraticPlusC
    };

    bool compile(const std::string& source, ExpressionError* error = nullptr);
    Complex evaluate(const ExpressionContext& context) const;
    Complex evaluate(const ExpressionContext& context, Complex* stack,
                     size_t capacity) const;

    bool valid() const { return _valid; }
    const std::string& source() const { return _source; }
    size_t instructionCount() const { return _code.size(); }
    size_t stackDepth() const { return _stackDepth; }
    FastPath fastPath() const { return _fastPath; }

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

    std::vector<Instruction> _code;
    std::string _source;
    size_t _stackDepth = 0;
    FastPath _fastPath = FastPath::None;
    bool _valid = false;

    friend class ExpressionParser;
    friend class ExpressionOracle;
};

} // namespace formula

#endif
