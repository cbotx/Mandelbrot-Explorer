#include "formula_expression_jit.h"

#if !defined(MANDEL_VERIFY_NO_JIT)
#error formula_expression_jit_stub.cpp is only for the non-JIT verifier
#endif

namespace formula {

struct ExpressionJit4::Impl {};

void ExpressionJitInput4::setContexts(
        const ExpressionContext* contexts) {
    if (!contexts) return;
    for (int lane = 0; lane < 4; ++lane)
        setContextLane(lane, contexts[lane]);
}

void ExpressionJitInput4::setContextLane(
        int lane, const ExpressionContext& context) {
    if (lane < 0 || lane >= 4) return;
    vectors[Z_RE][lane] = context.z.real();
    vectors[Z_IM][lane] = context.z.imag();
    vectors[C_RE][lane] = context.c.real();
    vectors[C_IM][lane] = context.c.imag();
    vectors[Z0_RE][lane] = context.z0.real();
    vectors[Z0_IM][lane] = context.z0.imag();
    vectors[N_RE][lane] =
        static_cast<double>(context.iteration);
    vectors[N_IM][lane] = 0.0;
    for (int parameter = 0; parameter < 8; ++parameter) {
        vectors[PARAM_RE + 2 * parameter][lane] =
            context.parameters[parameter].real();
        vectors[PARAM_RE + 2 * parameter + 1][lane] =
            context.parameters[parameter].imag();
    }
}

void ExpressionJitInvariantInput4::setPreparedLane(
        int lane, const ExpressionOrbitPlan& plan,
        const ExpressionOrbitPlan::Prepared& prepared) {
    if (lane < 0 || lane >= 4) return;
    for (size_t index = 0;
         index < plan.invariantCount(); ++index) {
        vectors[2 * index][lane] =
            prepared.values[index].real();
        vectors[2 * index + 1][lane] =
            prepared.values[index].imag();
    }
}

ExpressionJit4::ExpressionJit4()
    : _impl(std::make_unique<Impl>()) {}
ExpressionJit4::~ExpressionJit4() = default;
ExpressionJit4::ExpressionJit4(
    ExpressionJit4&&) noexcept = default;
ExpressionJit4& ExpressionJit4::operator=(
    ExpressionJit4&&) noexcept = default;

bool ExpressionJit4::compile(
        const ExpressionProgram&, std::string* error) {
    if (error) *error = "JIT is disabled in this build";
    return false;
}

bool ExpressionJit4::compile(
        const ExpressionOrbitPlan&, std::string* error) {
    if (error) *error = "JIT is disabled in this build";
    return false;
}

bool ExpressionJit4::compileProgram(
        const ExpressionProgram&, size_t,
        const std::string*, std::string* error) {
    if (error) *error = "JIT is disabled in this build";
    return false;
}

void ExpressionJit4::reset() {}
bool ExpressionJit4::valid() const { return false; }
bool ExpressionJit4::supports(
    const ExpressionOrbitPlan&) const { return false; }
bool ExpressionJit4::usesDualMapping() const { return false; }
const void* ExpressionJit4::codeAddress() const { return nullptr; }

void ExpressionJit4::evaluate(
        const ExpressionJitInput4&,
        ExpressionJitOutput4&) const {}
void ExpressionJit4::evaluate(
        const ExpressionJitInput4&,
        const ExpressionJitInvariantInput4*,
        ExpressionJitOutput4&) const {}
bool ExpressionJit4::evaluate(
        const ExpressionContext*, Complex*) const {
    return false;
}
bool ExpressionJit4::evaluateOrbit(
        const ExpressionContext*, int, int, double,
        float*, const volatile bool*,
        const ExpressionOrbitPlan*) const {
    return false;
}

} // namespace formula
