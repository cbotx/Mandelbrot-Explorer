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
class ExpressionOrbitPlan;
class ExpressionCenteredEvaluator;
class ExpressionScaledResidualEvaluator;
class ExpressionTaylorJetBuilder;
class ExpressionRealTaylorJetBuilder;

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

enum class ExpressionColoring : uint8_t { Raw,
                                          Smooth,
                                          Distance,
                                          Feather,
                                          OrbitTrap };

enum class ExpressionScaledResidualCapability : uint8_t { ExactCenteredArithmetic,
                                                          CertifiedEntireCandidate,
                                                          CertifiedMeromorphicCandidate,
                                                          CertifiedBranchCandidate,
                                                          CertifiedRealCandidate,
                                                          CertifiedPiecewiseCandidate,
                                                          UncertifiedSeries,
                                                          BranchSensitive,
                                                          Unsupported };

enum class ExpressionPiecewiseQuadraticKind : uint8_t { None,
                                                        BurningShip };

class ExpressionProgram {
  public:
    static constexpr size_t MAX_SOURCE = 4096;
    static constexpr size_t MAX_INSTRUCTIONS = 512;
    static constexpr size_t MAX_STACK = 128;
    static constexpr size_t MAX_PARSE_DEPTH = 64;

    enum class FastPath : uint8_t { None,
                                    IntegerPowerPlusC };

    ExpressionProgram();
    ExpressionProgram(const ExpressionProgram& other);
    ExpressionProgram(ExpressionProgram&& other) noexcept;
    ExpressionProgram& operator=(const ExpressionProgram& other);
    ExpressionProgram& operator=(ExpressionProgram&& other) noexcept;

    bool compile(const std::string& source, ExpressionError* error = nullptr);
    bool specialize(const ExpressionContext& fixed, FormulaParameter pixelParameter, ExpressionProgram& output, ExpressionError* error = nullptr) const;
    Complex evaluate(const ExpressionContext& context) const;
    Complex evaluate(const ExpressionContext& context, Complex* stack, size_t capacity) const;
    bool evaluate4(const ExpressionContext* contexts, Complex* outputs) const;
    bool evaluate4Hybrid(const ExpressionContext* contexts, Complex* outputs, int vectorTranscendentalMode = 0) const;
    bool evaluateWithDerivative(const ExpressionContext& context, const ExpressionDerivativeSeed& seed, Complex& value, Complex& derivative) const;

    bool valid() const { return _valid; }
    const std::string& source() const { return _source; }
    size_t instructionCount() const { return _code.size(); }
    size_t stackDepth() const { return _stackDepth; }
    uint64_t semanticHash() const;
    bool semanticallyEquivalent(const ExpressionProgram& other) const;
    bool containsOrbitInvariant() const;
    bool iterationDependent() const;
    FastPath fastPath() const { return _fastPath; }
    int fastIntegerPower() const { return _fastIntegerPower; }
    bool isCanonicalQuadraticPlusC() const;
    bool avx2Compatible() const { return _avx2Compatible; }
    bool batchCompatible() const { return _batchCompatible; }
    bool derivativeCompatible() const { return _derivativeCompatible; }
    ExpressionScaledResidualCapability scaledResidualCapability() const;
    bool scaledResidualRequiresTaylor() const;
    ExpressionPiecewiseQuadraticKind piecewiseQuadraticKind() const;

  private:
    enum class Op : uint8_t { Constant,
                              Z,
                              C,
                              Z0,
                              Iteration,
                              Parameter,
                              OrbitInvariant,
                              Negate,
                              Add,
                              Subtract,
                              Multiply,
                              Divide,
                              Power,
                              Square,
                              Sin,
                              Cos,
                              Tan,
                              Sinh,
                              Cosh,
                              Tanh,
                              Exp,
                              Log,
                              Log10,
                              Sqrt,
                              Abs,
                              Norm,
                              Arg,
                              Conjugate,
                              Real,
                              Imaginary,
                              MakeComplex,
                              Polar };

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
    uint64_t _identity = 0;
    uint64_t _revision = 1;

    bool analyze(ExpressionError* error);
    static int operandCount(Op op);
    static Complex evaluateUnary(Op op, Complex value);
    static Complex evaluateBinary(Op op, Complex left, Complex right);
    Complex evaluatePrepared(const ExpressionContext& context, const Complex* invariants, size_t invariantCount, Complex* stack, size_t capacity) const;
    bool evaluate4Prepared(const ExpressionContext* contexts, const Complex* const* invariants, size_t invariantCount, Complex* outputs) const;
    bool evaluate4HybridPrepared(const ExpressionContext* contexts, const Complex* const* invariants, size_t invariantCount, Complex* outputs, int vectorTranscendentalMode = 0) const;

    friend class ExpressionParser;
    friend class ExpressionOracle;
    friend class ExpressionJit4;
    friend class ExpressionOrbitPlan;
    friend class ExpressionCenteredEvaluator;
    friend class ExpressionScaledResidualEvaluator;
    friend class ExpressionTaylorJetBuilder;
    friend class ExpressionRealTaylorJetBuilder;
};

class ExpressionOrbitPlan {
  public:
    static constexpr size_t MAX_INVARIANTS = 255;

    enum Dependency : uint8_t { DependencyNone = 0,
                                DependencyZ = 1 << 0,
                                DependencyIteration = 1 << 1,
                                DependencyC = 1 << 2,
                                DependencyZ0 = 1 << 3,
                                DependencyParameter = 1 << 4 };

    struct Prepared {
        std::array<Complex, MAX_INVARIANTS> values{};
    };

    bool build(const ExpressionProgram& program, ExpressionError* error = nullptr);
    void reset();

    bool valid() const { return _valid; }
    bool profitable() const { return _valid && !_invariantPrograms.empty(); }
    bool matches(const ExpressionProgram& program) const;
    const std::string& source() const { return _source; }
    uint8_t dependencyMask() const { return _dependencyMask; }
    size_t invariantCount() const { return _invariantPrograms.size(); }
    uint8_t invariantDependencyMask(size_t index) const { return index < _invariantDependencies.size() ? _invariantDependencies[index] : DependencyNone; }
    size_t invariantInstructionCount(size_t index) const { return index < _invariantPrograms.size() ? _invariantPrograms[index].instructionCount() : 0; }
    size_t invariantOperationCount() const { return _invariantOperationCount; }
    size_t bodyInstructionCount() const { return _body.instructionCount(); }
    const ExpressionProgram& bodyProgram() const { return _body; }
    size_t bodyOperationCount() const { return _bodyOperationCount; }
    size_t bodyStackDepth() const { return _body.stackDepth(); }
    bool avx2Compatible() const { return _valid && _body.avx2Compatible(); }
    bool batchCompatible() const { return _valid && _body.batchCompatible(); }

    bool prepare(const ExpressionContext& context, Prepared& prepared) const;
    bool prepare4(const ExpressionContext* contexts, uint8_t activeMask, Prepared* prepared) const;
    Complex evaluate(const ExpressionContext& context, const Prepared& prepared) const;
    Complex evaluate(const ExpressionContext& context, const Prepared& prepared, Complex* stack, size_t capacity) const;
    bool evaluate4(const ExpressionContext* contexts, const Prepared* prepared, Complex* outputs) const;
    bool evaluate4Hybrid(const ExpressionContext* contexts, const Prepared* prepared, Complex* outputs, int vectorTranscendentalMode = 0) const;

  private:
    ExpressionProgram _body;
    std::vector<ExpressionProgram> _invariantPrograms;
    std::vector<uint8_t> _invariantDependencies;
    std::vector<ExpressionProgram::Instruction> _originalCode;
    std::string _source;
    std::string _programKey;
    size_t _invariantOperationCount = 0;
    size_t _bodyOperationCount = 0;
    uint8_t _dependencyMask = DependencyNone;
    bool _programAvx2Compatible = false;
    bool _valid = false;

    friend class ExpressionJit4;
};

} // namespace formula

#endif
