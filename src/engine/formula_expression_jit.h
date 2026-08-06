#ifndef MANDEL_FORMULA_EXPRESSION_JIT_H
#define MANDEL_FORMULA_EXPRESSION_JIT_H

#include <array>
#include <complex>
#include <memory>
#include <string>

#include "formula_expression.h"

namespace formula {

struct alignas(32) ExpressionJitInput4 {
    enum VectorIndex : size_t {
        Z_RE = 0, Z_IM,
        C_RE, C_IM,
        Z0_RE, Z0_IM,
        N_RE, N_IM,
        PARAM_RE = 8
    };
    static constexpr size_t VECTOR_COUNT = 24;
    double vectors[VECTOR_COUNT][4];

    void setContexts(const ExpressionContext* contexts);
    void setContextLane(int lane, const ExpressionContext& context);
};

struct alignas(32) ExpressionJitInvariantInput4 {
    static constexpr size_t VECTOR_COUNT =
        2 * ExpressionOrbitPlan::MAX_INVARIANTS;
    double vectors[VECTOR_COUNT][4];

    void setPreparedLane(int lane, const ExpressionOrbitPlan& plan,
                         const ExpressionOrbitPlan::Prepared& prepared);
};

struct alignas(32) ExpressionJitOutput4 {
    double re[4]{};
    double im[4]{};
};

class ExpressionJit4 {
public:
    ExpressionJit4();
    ~ExpressionJit4();
    ExpressionJit4(ExpressionJit4&&) noexcept;
    ExpressionJit4& operator=(ExpressionJit4&&) noexcept;
    ExpressionJit4(const ExpressionJit4&) = delete;
    ExpressionJit4& operator=(const ExpressionJit4&) = delete;

    bool compile(const ExpressionProgram& program, std::string* error = nullptr);
    bool compile(const ExpressionOrbitPlan& plan,
                 std::string* error = nullptr);
    void reset();
    bool valid() const;
    bool supports(const ExpressionOrbitPlan& plan) const;
    bool usesDualMapping() const;
    const void* codeAddress() const;
    void evaluate(const ExpressionJitInput4& input, ExpressionJitOutput4& output) const;
    void evaluate(const ExpressionJitInput4& input,
                  const ExpressionJitInvariantInput4* invariants,
                  ExpressionJitOutput4& output) const;
    bool evaluate(const ExpressionContext* contexts, Complex* outputs) const;
    bool evaluateOrbit(const ExpressionContext* contexts, int lanes,
                       int mxit, double bailout, float* results,
                       const volatile bool* halt = nullptr,
                       const ExpressionOrbitPlan* plan = nullptr) const;

private:
    struct Impl;
    std::unique_ptr<Impl> _impl;
    bool compileProgram(const ExpressionProgram& program,
                        size_t invariantCount,
                        const std::string* planKey,
                        std::string* error);
};

} // namespace formula

#endif
