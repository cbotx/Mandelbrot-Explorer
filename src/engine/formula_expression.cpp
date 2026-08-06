#include "formula_expression.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iterator>
#include <limits>
#include <string_view>
#include <unordered_map>
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

int ExpressionProgram::operandCount(Op op) {
    switch (op) {
    case Op::Constant:
    case Op::Z:
    case Op::C:
    case Op::Z0:
    case Op::Iteration:
    case Op::Parameter:
    case Op::OrbitInvariant:
        return 0;
    case Op::Negate:
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
    case Op::Abs:
    case Op::Norm:
    case Op::Arg:
    case Op::Conjugate:
    case Op::Real:
    case Op::Imaginary:
        return 1;
    case Op::Add:
    case Op::Subtract:
    case Op::Multiply:
    case Op::Divide:
    case Op::Power:
    case Op::MakeComplex:
    case Op::Polar:
        return 2;
    }
    return -1;
}

Complex ExpressionProgram::evaluateUnary(Op op, Complex value) {
    switch (op) {
    case Op::Negate: return -value;
    case Op::Square: return value * value;
    case Op::Sin: return std::sin(value);
    case Op::Cos: return std::cos(value);
    case Op::Tan: return std::tan(value);
    case Op::Sinh: return std::sinh(value);
    case Op::Cosh: return std::cosh(value);
    case Op::Tanh: return std::tanh(value);
    case Op::Exp: return std::exp(value);
    case Op::Log: return std::log(value);
    case Op::Log10: return std::log10(value);
    case Op::Sqrt: return std::sqrt(value);
    case Op::Abs: return { std::abs(value), 0.0 };
    case Op::Norm: return { std::norm(value), 0.0 };
    case Op::Arg: return { std::arg(value), 0.0 };
    case Op::Conjugate: return std::conj(value);
    case Op::Real: return { value.real(), 0.0 };
    case Op::Imaginary: return { value.imag(), 0.0 };
    default:
        break;
    }
    const double nan = std::numeric_limits<double>::quiet_NaN();
    return { nan, nan };
}

Complex ExpressionProgram::evaluateBinary(Op op, Complex left, Complex right) {
    switch (op) {
    case Op::Add: return left + right;
    case Op::Subtract: return left - right;
    case Op::Multiply: return left * right;
    case Op::Divide: return left / right;
    case Op::Power: return std::pow(left, right);
    case Op::MakeComplex: return { left.real(), right.real() };
    case Op::Polar: {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        if (left.imag() != 0.0 || right.imag() != 0.0 ||
            !std::isfinite(left.real()) || left.real() < 0.0 ||
            !std::isfinite(right.real()))
            return { nan, nan };
        return std::polar(left.real(), right.real());
    }
    default:
        break;
    }
    const double nan = std::numeric_limits<double>::quiet_NaN();
    return { nan, nan };
}

bool ExpressionProgram::analyze(ExpressionError* error) {
    _stackDepth = 0;
    _fastPath = FastPath::None;
    _fastIntegerPower = 0;
    _avx2Compatible = false;
    _batchCompatible = false;
    _derivativeCompatible = false;
    _valid = false;
    auto fail = [&](const char* message) {
        if (error) {
            error->position = _source.size();
            error->message = message;
        }
        return false;
    };
    if (_code.empty()) return fail("expression bytecode is empty");
    if (_code.size() > MAX_INSTRUCTIONS)
        return fail("expression has too many operations");

    size_t stack = 0;
    for (const Instruction& instruction : _code) {
        int operands = operandCount(instruction.op);
        if (operands < 0) return fail("expression contains an invalid operation");
        if (instruction.op == Op::Parameter && instruction.argument >= 8)
            return fail("expression contains an invalid parameter");
        if (stack < (size_t)operands)
            return fail("invalid expression stack");
        stack = stack - (size_t)operands + 1;
        _stackDepth = std::max(_stackDepth, stack);
        if (_stackDepth > MAX_STACK)
            return fail("expression stack is too deep");
    }
    if (stack != 1) return fail("invalid expression stack");

    const auto is = [&](size_t index, Op op) {
        return index < _code.size() && _code[index].op == op;
    };
    const auto isRecurrenceC = [&](size_t index) {
        return is(index, Op::C) ||
               (is(index, Op::Constant) &&
                _code[index].argument == CONSTANT_FIXED_C);
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
        if (!(cursor + 2 == _code.size() && isRecurrenceC(cursor) &&
              is(cursor + 1, Op::Add)))
            degree = 0;
    } else if (_code.size() == 4 && is(0, Op::Z) && is(1, Op::Square) &&
               isRecurrenceC(2) && is(3, Op::Add)) {
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
            case Op::OrbitInvariant:
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
    _batchCompatible = std::all_of(_code.begin(), _code.end(),
        [](const Instruction& instruction) {
            switch (instruction.op) {
            case Op::Constant:
            case Op::Z:
            case Op::C:
            case Op::Z0:
            case Op::Iteration:
            case Op::Parameter:
            case Op::OrbitInvariant:
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
            case Op::Abs:
            case Op::Norm:
            case Op::Arg:
            case Op::Conjugate:
            case Op::Real:
            case Op::Imaginary:
            case Op::MakeComplex:
            case Op::Polar:
                return true;
            }
            return false;
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

bool ExpressionProgram::compile(const std::string& source, ExpressionError* error) {
    _valid = false;
    _source.clear();
    _code.clear();
    _stackDepth = 0;
    _fastPath = FastPath::None;
    _fastIntegerPower = 0;
    _avx2Compatible = false;
    _batchCompatible = false;
    _derivativeCompatible = false;
    if (error) *error = {};
    ExpressionParser parser(*this, source, error);
    if (!parser.parse()) {
        _code.clear();
        return false;
    }
    _source = source;
    if (!analyze(error)) {
        _source.clear();
        _code.clear();
        return false;
    }
    return true;
}

bool ExpressionProgram::specialize(
        const ExpressionContext& fixed, FormulaParameter pixelParameter,
        ExpressionProgram& output, ExpressionError* error) const {
    if (error) *error = {};
    auto fail = [&](const char* message) {
        if (error) {
            error->position = 0;
            error->message = message;
        }
        return false;
    };
    if (!_valid) return fail("expression program is invalid");
    if (pixelParameter != FormulaParameter::C &&
        pixelParameter != FormulaParameter::InitialZ)
        return fail("unsupported pixel parameter");

    struct Node {
        bool constant = false;
        Complex value{};
        uint8_t constantArgument = 0;
        std::vector<Instruction> code;
    };
    std::vector<Node> stack;
    stack.reserve(_stackDepth);

    auto materialize = [](Node& node) {
        if (node.constant) {
            node.code.push_back({
                Op::Constant, node.constantArgument, node.value
            });
            node.constant = false;
        }
    };

    ExpressionProgram candidate;
    candidate._source = _source;
    for (const Instruction& instruction : _code) {
        int operands = operandCount(instruction.op);
        if (operands == 0) {
            Node node;
            switch (instruction.op) {
            case Op::Constant:
                node.constant = true;
                node.value = instruction.value;
                node.constantArgument = instruction.argument;
                break;
            case Op::C:
                if (pixelParameter == FormulaParameter::InitialZ) {
                    node.constant = true;
                    node.value = fixed.c;
                    node.constantArgument = CONSTANT_FIXED_C;
                } else {
                    node.code.push_back(instruction);
                }
                break;
            case Op::Z0:
                if (pixelParameter == FormulaParameter::C) {
                    node.constant = true;
                    node.value = fixed.z0;
                } else {
                    node.code.push_back(instruction);
                }
                break;
            case Op::Parameter:
                if (instruction.argument >= fixed.parameters.size())
                    return fail("expression contains an invalid parameter");
                node.constant = true;
                node.value = fixed.parameters[instruction.argument];
                break;
            case Op::Z:
            case Op::Iteration:
                node.code.push_back(instruction);
                break;
            default:
                return fail("expression contains an invalid operation");
            }
            stack.push_back(std::move(node));
            continue;
        }
        if (stack.size() < (size_t)operands)
            return fail("invalid expression stack");
        if (operands == 1) {
            Node node = std::move(stack.back());
            stack.pop_back();
            if (node.constant) {
                node.value = evaluateUnary(instruction.op, node.value);
                node.constantArgument = 0;
            } else {
                node.code.push_back(instruction);
            }
            stack.push_back(std::move(node));
            continue;
        }
        if (operands != 2)
            return fail("expression contains an invalid operation");
        Node right = std::move(stack.back());
        stack.pop_back();
        Node left = std::move(stack.back());
        stack.pop_back();
        Node node;
        if (left.constant && right.constant) {
            node.constant = true;
            node.value = evaluateBinary(
                instruction.op, left.value, right.value);
        } else {
            materialize(left);
            materialize(right);
            node.code.reserve(
                left.code.size() + right.code.size() + 1);
            node.code.insert(
                node.code.end(),
                std::make_move_iterator(left.code.begin()),
                std::make_move_iterator(left.code.end()));
            node.code.insert(
                node.code.end(),
                std::make_move_iterator(right.code.begin()),
                std::make_move_iterator(right.code.end()));
            node.code.push_back(instruction);
        }
        stack.push_back(std::move(node));
    }
    if (stack.size() != 1) return fail("invalid expression stack");
    materialize(stack.back());
    candidate._code = std::move(stack.back().code);
    if (!candidate.analyze(error)) return false;
    output = std::move(candidate);
    return true;
}

bool ExpressionOrbitPlan::build(
        const ExpressionProgram& program, ExpressionError* error) {
    if (error) *error = {};
    auto fail = [&](const char* message) {
        if (error) {
            error->position = program.source().size();
            error->message = message;
        }
        return false;
    };
    if (!program.valid()) return fail("expression program is invalid");

    struct Node {
        ExpressionProgram::Instruction instruction;
        int operands = 0;
        int left = -1;
        int right = -1;
        size_t instructionCount = 1;
        uint8_t dependencies = DependencyNone;
        std::string key;
    };
    std::vector<Node> nodes;
    std::vector<int> stack;
    nodes.reserve(program._code.size());
    stack.reserve(program._stackDepth);

    auto appendBytes = [](std::string& key, const void* data, size_t size) {
        key.append(static_cast<const char*>(data), size);
    };
    for (const ExpressionProgram::Instruction& instruction : program._code) {
        Node node;
        node.instruction = instruction;
        node.operands = ExpressionProgram::operandCount(instruction.op);
        if (node.operands < 0 || stack.size() < (size_t)node.operands)
            return fail("invalid expression stack");
        if (node.operands == 1) {
            node.left = stack.back();
            stack.pop_back();
        } else if (node.operands == 2) {
            node.right = stack.back();
            stack.pop_back();
            node.left = stack.back();
            stack.pop_back();
        }

        switch (instruction.op) {
        case ExpressionProgram::Op::Z:
            node.dependencies = DependencyZ;
            break;
        case ExpressionProgram::Op::Iteration:
            node.dependencies = DependencyIteration;
            break;
        case ExpressionProgram::Op::C:
            node.dependencies = DependencyC;
            break;
        case ExpressionProgram::Op::Z0:
            node.dependencies = DependencyZ0;
            break;
        case ExpressionProgram::Op::Parameter:
            node.dependencies = DependencyParameter;
            break;
        case ExpressionProgram::Op::OrbitInvariant:
            return fail("nested orbit plan is not supported");
        default:
            break;
        }
        if (node.left >= 0) {
            node.dependencies |= nodes[(size_t)node.left].dependencies;
            node.instructionCount +=
                nodes[(size_t)node.left].instructionCount;
        }
        if (node.right >= 0) {
            node.dependencies |= nodes[(size_t)node.right].dependencies;
            node.instructionCount +=
                nodes[(size_t)node.right].instructionCount;
        }

        const uint8_t op = static_cast<uint8_t>(instruction.op);
        node.key.push_back(static_cast<char>(op));
        node.key.push_back(static_cast<char>(instruction.argument));
        uint64_t realBits = 0, imaginaryBits = 0;
        const double real = instruction.value.real();
        const double imaginary = instruction.value.imag();
        std::memcpy(&realBits, &real, sizeof(realBits));
        std::memcpy(&imaginaryBits, &imaginary, sizeof(imaginaryBits));
        appendBytes(node.key, &realBits, sizeof(realBits));
        appendBytes(node.key, &imaginaryBits, sizeof(imaginaryBits));
        node.key.push_back(static_cast<char>(node.operands));
        auto appendChildKey = [&](int child) {
            const std::string& childKey = nodes[(size_t)child].key;
            const uint32_t size = static_cast<uint32_t>(childKey.size());
            appendBytes(node.key, &size, sizeof(size));
            appendBytes(node.key, childKey.data(), childKey.size());
        };
        if (node.left >= 0) appendChildKey(node.left);
        if (node.right >= 0) appendChildKey(node.right);

        nodes.push_back(std::move(node));
        stack.push_back(static_cast<int>(nodes.size() - 1));
    }
    if (stack.size() != 1) return fail("invalid expression stack");
    const int root = stack.back();

    ExpressionOrbitPlan candidate;
    candidate._source = program._source;
    candidate._programKey = nodes[(size_t)root].key;
    candidate._originalCode = program._code;
    candidate._dependencyMask = nodes[(size_t)root].dependencies;
    candidate._programAvx2Compatible = program._avx2Compatible;
    candidate._body._source = program._source;

    const uint8_t iterationDependencies =
        DependencyZ | DependencyIteration;
    auto isInvariant = [&](const Node& node) {
        return (node.dependencies & iterationDependencies) == 0 &&
               node.instructionCount > 1;
    };
    std::vector<int> maximalInvariantRoots;
    std::function<void(int)> collectMaximal;
    collectMaximal = [&](int index) {
        const Node& node = nodes[(size_t)index];
        if (isInvariant(node)) {
            maximalInvariantRoots.push_back(index);
            return;
        }
        if (node.left >= 0) collectMaximal(node.left);
        if (node.right >= 0) collectMaximal(node.right);
    };
    collectMaximal(root);

    std::unordered_map<std::string, size_t> invariantOccurrences;
    std::function<void(int)> countInvariantNodes;
    countInvariantNodes = [&](int index) {
        const Node& node = nodes[(size_t)index];
        if (isInvariant(node)) ++invariantOccurrences[node.key];
        if (node.left >= 0) countInvariantNodes(node.left);
        if (node.right >= 0) countInvariantNodes(node.right);
    };
    for (int index : maximalInvariantRoots)
        countInvariantNodes(index);

    std::unordered_map<std::string, uint8_t> invariantIndices;
    bool okay = true;
    std::function<uint8_t(int)> internInvariant;
    std::function<void(
        int, std::vector<ExpressionProgram::Instruction>&)> emitInvariant;
    emitInvariant = [&](int index,
                        std::vector<ExpressionProgram::Instruction>& code) {
        if (!okay) return;
        const Node& node = nodes[(size_t)index];
        auto emitChild = [&](int child) {
            if (child < 0 || !okay) return;
            const Node& childNode = nodes[(size_t)child];
            auto occurrence = invariantOccurrences.find(childNode.key);
            if (isInvariant(childNode) &&
                occurrence != invariantOccurrences.end() &&
                occurrence->second > 1) {
                uint8_t childIndex = internInvariant(child);
                code.push_back({
                    ExpressionProgram::Op::OrbitInvariant,
                    childIndex, {}
                });
            } else {
                emitInvariant(child, code);
            }
        };
        emitChild(node.left);
        emitChild(node.right);
        if (okay) code.push_back(node.instruction);
    };
    internInvariant = [&](int index) -> uint8_t {
        const Node& node = nodes[(size_t)index];
        auto found = invariantIndices.find(node.key);
        if (found != invariantIndices.end()) return found->second;
        if (candidate._invariantPrograms.size() >= MAX_INVARIANTS) {
            okay = false;
            return 0;
        }
        ExpressionProgram invariantProgram;
        invariantProgram._source = program._source;
        emitInvariant(index, invariantProgram._code);
        if (!okay || !invariantProgram.analyze(error)) {
            okay = false;
            return 0;
        }
        uint8_t invariantIndex = static_cast<uint8_t>(
            candidate._invariantPrograms.size());
        invariantIndices.emplace(node.key, invariantIndex);
        candidate._invariantDependencies.push_back(node.dependencies);
        candidate._invariantOperationCount +=
            std::count_if(
                invariantProgram._code.begin(),
                invariantProgram._code.end(),
                [](const ExpressionProgram::Instruction& item) {
                    return ExpressionProgram::operandCount(item.op) > 0;
                });
        candidate._invariantPrograms.push_back(
            std::move(invariantProgram));
        return invariantIndex;
    };
    std::function<void(int)> emitBody;
    emitBody = [&](int index) {
        if (!okay) return;
        const Node& node = nodes[(size_t)index];
        if (isInvariant(node)) {
            uint8_t invariantIndex = internInvariant(index);
            candidate._body._code.push_back({
                ExpressionProgram::Op::OrbitInvariant,
                invariantIndex, {}
            });
            return;
        }
        if (node.left >= 0) emitBody(node.left);
        if (node.right >= 0) emitBody(node.right);
        candidate._body._code.push_back(node.instruction);
    };
    emitBody(root);
    if (!okay) {
        if (!error || error->message.empty())
            return fail("expression has too many orbit invariants");
        return false;
    }
    if (!candidate._body.analyze(error)) return false;
    candidate._bodyOperationCount = std::count_if(
        candidate._body._code.begin(), candidate._body._code.end(),
        [](const ExpressionProgram::Instruction& instruction) {
            return ExpressionProgram::operandCount(instruction.op) > 0;
        });
    candidate._valid = true;
    *this = std::move(candidate);
    return true;
}

void ExpressionOrbitPlan::reset() {
    *this = ExpressionOrbitPlan{};
}

bool ExpressionOrbitPlan::matches(
        const ExpressionProgram& program) const {
    if (!_valid || !program.valid() ||
        _source != program._source ||
        _originalCode.size() != program._code.size())
        return false;
    for (size_t i = 0; i < _originalCode.size(); ++i) {
        const ExpressionProgram::Instruction& left = _originalCode[i];
        const ExpressionProgram::Instruction& right = program._code[i];
        if (left.op != right.op || left.argument != right.argument)
            return false;
        uint64_t leftReal = 0, leftImaginary = 0;
        uint64_t rightReal = 0, rightImaginary = 0;
        double value = left.value.real();
        std::memcpy(&leftReal, &value, sizeof(leftReal));
        value = left.value.imag();
        std::memcpy(&leftImaginary, &value, sizeof(leftImaginary));
        value = right.value.real();
        std::memcpy(&rightReal, &value, sizeof(rightReal));
        value = right.value.imag();
        std::memcpy(&rightImaginary, &value, sizeof(rightImaginary));
        if (leftReal != rightReal || leftImaginary != rightImaginary)
            return false;
    }
    return true;
}

bool ExpressionOrbitPlan::prepare(
        const ExpressionContext& context, Prepared& prepared) const {
    if (!_valid) return false;
    std::array<Complex, ExpressionProgram::MAX_STACK> stack;
    for (size_t i = 0; i < _invariantPrograms.size(); ++i) {
        prepared.values[i] = _invariantPrograms[i].evaluatePrepared(
            context, prepared.values.data(), i,
            stack.data(), stack.size());
    }
    return true;
}

bool ExpressionOrbitPlan::prepare4(
        const ExpressionContext* contexts, uint8_t activeMask,
        Prepared* prepared) const {
    if (!_valid || !contexts || !prepared) return false;
    activeMask &= 0x0f;
    if (activeMask == 0) return true;
    if (activeMask != 0x0f) {
        for (int lane = 0; lane < 4; ++lane) {
            if ((activeMask & (1u << lane)) &&
                !prepare(contexts[lane], prepared[lane]))
                return false;
        }
        return true;
    }
    for (size_t index = 0;
         index < _invariantPrograms.size(); ++index) {
        Complex values[4];
        const ExpressionProgram& invariant =
            _invariantPrograms[index];
        const Complex* invariants[4] = {
            prepared[0].values.data(), prepared[1].values.data(),
            prepared[2].values.data(), prepared[3].values.data()
        };
        bool evaluated = _programAvx2Compatible
            ? invariant.evaluate4Prepared(
                contexts, invariants, index, values)
            : invariant.evaluate4HybridPrepared(
                contexts, invariants, index, values);
        if (!evaluated) return false;
        for (int lane = 0; lane < 4; ++lane)
            prepared[lane].values[index] = values[lane];
    }
    return true;
}

Complex ExpressionOrbitPlan::evaluate(
        const ExpressionContext& context,
        const Prepared& prepared) const {
    std::array<Complex, ExpressionProgram::MAX_STACK> stack;
    return evaluate(
        context, prepared, stack.data(), stack.size());
}

Complex ExpressionOrbitPlan::evaluate(
        const ExpressionContext& context, const Prepared& prepared,
        Complex* stack, size_t capacity) const {
    if (!_valid) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return { nan, nan };
    }
    return _body.evaluatePrepared(
        context, prepared.values.data(), _invariantPrograms.size(),
        stack, capacity);
}

bool ExpressionOrbitPlan::evaluate4(
        const ExpressionContext* contexts, const Prepared* prepared,
        Complex* outputs) const {
    if (!_valid || !prepared) return false;
    const Complex* invariants[4] = {
        prepared[0].values.data(), prepared[1].values.data(),
        prepared[2].values.data(), prepared[3].values.data()
    };
    return _body.evaluate4Prepared(
        contexts, invariants, _invariantPrograms.size(), outputs);
}

bool ExpressionOrbitPlan::evaluate4Hybrid(
        const ExpressionContext* contexts, const Prepared* prepared,
        Complex* outputs) const {
    if (!_valid || !prepared) return false;
    const Complex* invariants[4] = {
        prepared[0].values.data(), prepared[1].values.data(),
        prepared[2].values.data(), prepared[3].values.data()
    };
    return _body.evaluate4HybridPrepared(
        contexts, invariants, _invariantPrograms.size(), outputs);
}

Complex ExpressionProgram::evaluate(const ExpressionContext& context) const {
    std::array<Complex, MAX_STACK> stack;
    return evaluate(context, stack.data(), stack.size());
}

Complex ExpressionProgram::evaluate(const ExpressionContext& context,
                                    Complex* stack, size_t capacity) const {
    return evaluatePrepared(context, nullptr, 0, stack, capacity);
}

Complex ExpressionProgram::evaluatePrepared(
        const ExpressionContext& context, const Complex* invariants,
        size_t invariantCount, Complex* stack, size_t capacity) const {
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
        case Op::OrbitInvariant:
            if (!invariants || instruction.argument >= invariantCount)
                return {
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN()
                };
            stack[top++] = invariants[instruction.argument];
            break;
        case Op::Negate:
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
        case Op::Abs:
        case Op::Norm:
        case Op::Arg:
        case Op::Conjugate:
        case Op::Real:
        case Op::Imaginary:
            unary([&](Complex value) {
                return evaluateUnary(instruction.op, value);
            });
            break;
        case Op::Add:
        case Op::Subtract:
        case Op::Multiply:
        case Op::Divide:
        case Op::Power:
        case Op::MakeComplex:
        case Op::Polar:
            binary([&](Complex left, Complex right) {
                return evaluateBinary(instruction.op, left, right);
            });
            break;
        }
    }
    return stack[0];
}

bool ExpressionProgram::evaluate4(const ExpressionContext* contexts,
                                  Complex* outputs) const {
    return evaluate4Prepared(contexts, nullptr, 0, outputs);
}

bool ExpressionProgram::evaluate4Prepared(
        const ExpressionContext* contexts,
        const Complex* const* invariants, size_t invariantCount,
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
        case Op::OrbitInvariant: {
            if (!invariants || instruction.argument >= invariantCount)
                return false;
            const size_t index = instruction.argument;
            stack[top++] = {
                _mm256_set_pd(invariants[3][index].real(),
                              invariants[2][index].real(),
                              invariants[1][index].real(),
                              invariants[0][index].real()),
                _mm256_set_pd(invariants[3][index].imag(),
                              invariants[2][index].imag(),
                              invariants[1][index].imag(),
                              invariants[0][index].imag())
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

bool ExpressionProgram::evaluate4Hybrid(const ExpressionContext* contexts,
                                        Complex* outputs) const {
    return evaluate4HybridPrepared(contexts, nullptr, 0, outputs);
}

bool ExpressionProgram::evaluate4HybridPrepared(
        const ExpressionContext* contexts,
        const Complex* const* invariants, size_t invariantCount,
        Complex* outputs) const {
    if (!_valid || !_batchCompatible || !contexts || !outputs) return false;
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
    auto unary = [&](auto function) {
        VecComplex& value = stack[top - 1];
        alignas(32) double re[4], im[4];
        _mm256_store_pd(re, value.re);
        _mm256_store_pd(im, value.im);
        for (int lane = 0; lane < 4; ++lane) {
            Complex result = function(Complex{ re[lane], im[lane] });
            re[lane] = result.real();
            im[lane] = result.imag();
        }
        value.re = _mm256_load_pd(re);
        value.im = _mm256_load_pd(im);
    };
    auto binary = [&](auto function) {
        VecComplex right = stack[--top];
        VecComplex& left = stack[top - 1];
        alignas(32) double leftRe[4], leftIm[4], rightRe[4], rightIm[4];
        _mm256_store_pd(leftRe, left.re);
        _mm256_store_pd(leftIm, left.im);
        _mm256_store_pd(rightRe, right.re);
        _mm256_store_pd(rightIm, right.im);
        for (int lane = 0; lane < 4; ++lane) {
            Complex result = function(Complex{ leftRe[lane], leftIm[lane] },
                                      Complex{ rightRe[lane], rightIm[lane] });
            leftRe[lane] = result.real();
            leftIm[lane] = result.imag();
        }
        left.re = _mm256_load_pd(leftRe);
        left.im = _mm256_load_pd(leftIm);
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
        case Op::OrbitInvariant: {
            if (!invariants || instruction.argument >= invariantCount)
                return false;
            const size_t index = instruction.argument;
            stack[top++] = {
                _mm256_set_pd(invariants[3][index].real(),
                              invariants[2][index].real(),
                              invariants[1][index].real(),
                              invariants[0][index].real()),
                _mm256_set_pd(invariants[3][index].imag(),
                              invariants[2][index].imag(),
                              invariants[1][index].imag(),
                              invariants[0][index].imag())
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
        case Op::Divide:
            binary([](Complex a, Complex b) { return a / b; });
            break;
        case Op::Power:
            binary([](Complex a, Complex b) { return std::pow(a, b); });
            break;
        case Op::Square: {
            VecComplex& value = stack[top - 1];
            __m256d re = _mm256_sub_pd(_mm256_mul_pd(value.re, value.re),
                                       _mm256_mul_pd(value.im, value.im));
            __m256d im = _mm256_add_pd(_mm256_mul_pd(value.re, value.im),
                                       _mm256_mul_pd(value.im, value.re));
            value.re = re; value.im = im;
            break;
        }
        case Op::Sin:
            unary([](Complex a) { return std::sin(a); });
            break;
        case Op::Cos:
            unary([](Complex a) { return std::cos(a); });
            break;
        case Op::Tan:
            unary([](Complex a) { return std::tan(a); });
            break;
        case Op::Sinh:
            unary([](Complex a) { return std::sinh(a); });
            break;
        case Op::Cosh:
            unary([](Complex a) { return std::cosh(a); });
            break;
        case Op::Tanh:
            unary([](Complex a) { return std::tanh(a); });
            break;
        case Op::Exp:
            unary([](Complex a) { return std::exp(a); });
            break;
        case Op::Log:
            unary([](Complex a) { return std::log(a); });
            break;
        case Op::Log10:
            unary([](Complex a) { return std::log10(a); });
            break;
        case Op::Sqrt:
            unary([](Complex a) { return std::sqrt(a); });
            break;
        case Op::Abs:
            unary([](Complex a) { return Complex{ std::abs(a), 0.0 }; });
            break;
        case Op::Norm:
            unary([](Complex a) { return Complex{ std::norm(a), 0.0 }; });
            break;
        case Op::Arg:
            unary([](Complex a) { return Complex{ std::arg(a), 0.0 }; });
            break;
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
