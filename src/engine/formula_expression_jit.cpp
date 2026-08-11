#include "formula_expression_jit.h"

#include <asmjit/x86.h>

#include <algorithm>
#include <cmath>
#include <cstring>

namespace formula {

using JitFunction = void (*)(const double*, const double*, double*);

struct ExpressionJit4::Impl {
    asmjit::JitAllocator::CreateParams params;
    asmjit::JitRuntime runtime;
    JitFunction function = nullptr;
    std::string planKey;
    size_t invariantCount = 0;

    Impl()
        : params(makeParams()),
          runtime(&params) {}

    ~Impl() {
        if (function) runtime.release(function);
    }

    static asmjit::JitAllocator::CreateParams makeParams() {
        asmjit::JitAllocator::CreateParams result;
        result.options = asmjit::JitAllocatorOptions::kUseDualMapping |
                         asmjit::JitAllocatorOptions::kFillUnusedMemory |
                         asmjit::JitAllocatorOptions::kImmediateRelease;
        return result;
    }
};

void ExpressionJitInput4::setContexts(const ExpressionContext* contexts) {
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
    vectors[N_RE][lane] = (double)context.iteration;
    vectors[N_IM][lane] = 0.0;
    for (int p = 0; p < 8; ++p) {
        vectors[PARAM_RE + 2 * p][lane] =
            context.parameters[p].real();
        vectors[PARAM_RE + 2 * p + 1][lane] =
            context.parameters[p].imag();
    }
}

void ExpressionJitInvariantInput4::setPreparedLane(
        int lane, const ExpressionOrbitPlan& plan,
        const ExpressionOrbitPlan::Prepared& prepared) {
    if (lane < 0 || lane >= 4) return;
    for (size_t index = 0; index < plan.invariantCount(); ++index) {
        vectors[2 * index][lane] =
            prepared.values[index].real();
        vectors[2 * index + 1][lane] =
            prepared.values[index].imag();
    }
}

ExpressionJit4::ExpressionJit4() : _impl(std::make_unique<Impl>()) {}
ExpressionJit4::~ExpressionJit4() = default;
ExpressionJit4::ExpressionJit4(ExpressionJit4&&) noexcept = default;
ExpressionJit4& ExpressionJit4::operator=(ExpressionJit4&&) noexcept = default;

bool ExpressionJit4::valid() const {
    return _impl && _impl->function;
}

bool ExpressionJit4::supports(const ExpressionOrbitPlan& plan) const {
    return valid() && plan.valid() &&
           _impl->invariantCount == plan.invariantCount() &&
           _impl->planKey == plan._programKey;
}

bool ExpressionJit4::usesDualMapping() const {
    return _impl && _impl->runtime.allocator()->hasOption(
        asmjit::JitAllocatorOptions::kUseDualMapping);
}

const void* ExpressionJit4::codeAddress() const {
    return _impl ? reinterpret_cast<const void*>(_impl->function) : nullptr;
}

bool ExpressionJit4::compile(const ExpressionProgram& program, std::string* error) {
    return compileProgram(program, 0, nullptr, error);
}

bool ExpressionJit4::compile(
        const ExpressionOrbitPlan& plan, std::string* error) {
    if (!plan.valid() || !plan.profitable()) {
        reset();
        if (error) *error = "orbit plan is not profitable";
        return false;
    }
    return compileProgram(
        plan._body, plan.invariantCount(), &plan._programKey, error);
}

bool ExpressionJit4::compileProgram(
        const ExpressionProgram& program, size_t invariantCount,
        const std::string* planKey, std::string* error) {
    if (error) error->clear();
    if (!_impl) _impl = std::make_unique<Impl>();
    if (_impl->function) {
        _impl->runtime.release(_impl->function);
        _impl->function = nullptr;
    }
    _impl->planKey.clear();
    _impl->invariantCount = 0;
    if (!program.valid() || !program.avx2Compatible()) {
        if (error) *error = "program is not AVX2-JIT compatible";
        return false;
    }
    if (program.stackDepth() == 0 || program.stackDepth() > 32) {
        if (error) *error = "JIT stack depth exceeds 32";
        return false;
    }

    asmjit::CodeHolder code;
    asmjit::Error err = code.init(_impl->runtime.environment(),
                                  _impl->runtime.cpuFeatures());
    if (err) {
        if (error) *error = asmjit::DebugUtils::errorAsString(err);
        return false;
    }
    asmjit::x86::Assembler a(&code);
    using namespace asmjit::x86;
    const uint32_t frameSize = (uint32_t)((program.stackDepth() * 64 + 31) & ~31u);
    a.sub(rsp, frameSize);
    size_t top = 0;
    auto stackRe = [&](size_t index) { return ptr(rsp, (int32_t)(index * 64)); };
    auto stackIm = [&](size_t index) { return ptr(rsp, (int32_t)(index * 64 + 32)); };
    auto store = [&](size_t index, const Ymm& re, const Ymm& im) {
        a.vmovupd(stackRe(index), re);
        a.vmovupd(stackIm(index), im);
    };
    auto load = [&](size_t index, const Ymm& re, const Ymm& im) {
        a.vmovupd(re, stackRe(index));
        a.vmovupd(im, stackIm(index));
    };
    auto loadInput = [&](size_t vectorIndex, const Ymm& destination) {
        a.vmovupd(destination, ptr(rcx, (int32_t)(vectorIndex * 32)));
    };
    auto loadInvariant = [&](size_t vectorIndex,
                             const Ymm& destination) {
        a.vmovupd(destination, ptr(rdx, (int32_t)(vectorIndex * 32)));
    };
    auto broadcast = [&](double value, const Ymm& destination) {
        uint64_t bits;
        std::memcpy(&bits, &value, sizeof(bits));
        a.mov(rax, bits);
        a.vmovq(destination.xmm(), rax);
        a.vbroadcastsd(destination, destination.xmm());
    };
    Ymm yr0 = ymm0, yi0 = ymm1, yr1 = ymm2, yi1 = ymm3;
    Ymm tmp0 = ymm4, tmp1 = ymm5;

    for (const ExpressionProgram::Instruction& instruction : program._code) {
        switch (instruction.op) {
        case ExpressionProgram::Op::Constant:
            broadcast(instruction.value.real(), yr0);
            broadcast(instruction.value.imag(), yi0);
            store(top++, yr0, yi0);
            break;
        case ExpressionProgram::Op::Z:
            loadInput(ExpressionJitInput4::Z_RE, yr0);
            loadInput(ExpressionJitInput4::Z_IM, yi0);
            store(top++, yr0, yi0);
            break;
        case ExpressionProgram::Op::C:
            loadInput(ExpressionJitInput4::C_RE, yr0);
            loadInput(ExpressionJitInput4::C_IM, yi0);
            store(top++, yr0, yi0);
            break;
        case ExpressionProgram::Op::Z0:
            loadInput(ExpressionJitInput4::Z0_RE, yr0);
            loadInput(ExpressionJitInput4::Z0_IM, yi0);
            store(top++, yr0, yi0);
            break;
        case ExpressionProgram::Op::Iteration:
            loadInput(ExpressionJitInput4::N_RE, yr0);
            loadInput(ExpressionJitInput4::N_IM, yi0);
            store(top++, yr0, yi0);
            break;
        case ExpressionProgram::Op::Parameter: {
            size_t index = ExpressionJitInput4::PARAM_RE + 2 * instruction.argument;
            loadInput(index, yr0);
            loadInput(index + 1, yi0);
            store(top++, yr0, yi0);
            break;
        }
        case ExpressionProgram::Op::OrbitInvariant: {
            if (instruction.argument >= invariantCount) {
                if (error) *error = "invalid orbit invariant input";
                return false;
            }
            size_t index = 2 * instruction.argument;
            loadInvariant(index, yr0);
            loadInvariant(index + 1, yi0);
            store(top++, yr0, yi0);
            break;
        }
        case ExpressionProgram::Op::Negate:
            load(top - 1, yr0, yi0);
            broadcast(-0.0, tmp0);
            a.vxorpd(yr0, yr0, tmp0);
            a.vxorpd(yi0, yi0, tmp0);
            store(top - 1, yr0, yi0);
            break;
        case ExpressionProgram::Op::Add:
        case ExpressionProgram::Op::Subtract:
        case ExpressionProgram::Op::Multiply:
            load(top - 2, yr0, yi0);
            load(top - 1, yr1, yi1);
            --top;
            if (instruction.op == ExpressionProgram::Op::Add) {
                a.vaddpd(yr0, yr0, yr1);
                a.vaddpd(yi0, yi0, yi1);
            } else if (instruction.op == ExpressionProgram::Op::Subtract) {
                a.vsubpd(yr0, yr0, yr1);
                a.vsubpd(yi0, yi0, yi1);
            } else {
                a.vmulpd(tmp0, yr0, yr1);
                a.vmulpd(tmp1, yi0, yi1);
                a.vsubpd(tmp0, tmp0, tmp1);
                a.vmulpd(tmp1, yr0, yi1);
                a.vmulpd(yr0, yi0, yr1);
                a.vaddpd(yi0, tmp1, yr0);
                a.vmovapd(yr0, tmp0);
            }
            store(top - 1, yr0, yi0);
            break;
        case ExpressionProgram::Op::Square:
            load(top - 1, yr0, yi0);
            a.vmulpd(tmp0, yr0, yr0);
            a.vmulpd(tmp1, yi0, yi0);
            a.vsubpd(tmp0, tmp0, tmp1);
            a.vmulpd(tmp1, yr0, yi0);
            a.vmulpd(yr0, yi0, yr0);
            a.vaddpd(yi0, tmp1, yr0);
            a.vmovapd(yr0, tmp0);
            store(top - 1, yr0, yi0);
            break;
        case ExpressionProgram::Op::Conjugate:
            load(top - 1, yr0, yi0);
            broadcast(-0.0, tmp0);
            a.vxorpd(yi0, yi0, tmp0);
            store(top - 1, yr0, yi0);
            break;
        case ExpressionProgram::Op::Real:
            load(top - 1, yr0, yi0);
            a.vxorpd(yi0, yi0, yi0);
            store(top - 1, yr0, yi0);
            break;
        case ExpressionProgram::Op::Imaginary:
            load(top - 1, yr0, yi0);
            a.vmovapd(yr0, yi0);
            a.vxorpd(yi0, yi0, yi0);
            store(top - 1, yr0, yi0);
            break;
        case ExpressionProgram::Op::MakeComplex:
            load(top - 2, yr0, yi0);
            load(top - 1, yr1, yi1);
            --top;
            a.vmovapd(yi0, yr1);
            store(top - 1, yr0, yi0);
            break;
        default:
            if (error) *error = "unsupported JIT opcode";
            return false;
        }
    }

    load(0, yr0, yi0);
    a.vmovupd(ptr(r8, 0), yr0);
    a.vmovupd(ptr(r8, 32), yi0);
    a.vzeroupper();
    a.add(rsp, frameSize);
    a.ret();
    err = _impl->runtime.add(&_impl->function, &code);
    if (err) {
        if (error) *error = asmjit::DebugUtils::errorAsString(err);
        _impl->function = nullptr;
        return false;
    }
    _impl->invariantCount = invariantCount;
    if (planKey) _impl->planKey = *planKey;
    return true;
}

void ExpressionJit4::reset() {
    if (_impl && _impl->function) {
        _impl->runtime.release(_impl->function);
        _impl->function = nullptr;
    }
    if (_impl) {
        _impl->planKey.clear();
        _impl->invariantCount = 0;
    }
}

void ExpressionJit4::evaluate(const ExpressionJitInput4& input,
                              ExpressionJitOutput4& output) const {
    evaluate(input, nullptr, output);
}

void ExpressionJit4::evaluate(
        const ExpressionJitInput4& input,
        const ExpressionJitInvariantInput4* invariants,
        ExpressionJitOutput4& output) const {
    if (!_impl || !_impl->function ||
        (!_impl->planKey.empty() && !invariants))
        return;
    _impl->function(
        &input.vectors[0][0],
        invariants ? &invariants->vectors[0][0] : nullptr,
        &output.re[0]);
}

bool ExpressionJit4::evaluate(const ExpressionContext* contexts,
                              Complex* outputs) const {
    if (!valid() || !_impl->planKey.empty() || !contexts || !outputs)
        return false;
    ExpressionJitInput4 input;
    ExpressionJitOutput4 output;
    input.setContexts(contexts);
    evaluate(input, output);
    for (int lane = 0; lane < 4; ++lane)
        outputs[lane] = { output.re[lane], output.im[lane] };
    return true;
}

bool ExpressionJit4::evaluateOrbit(
        const ExpressionContext* contexts, int lanes,
        int mxit, double bailout, float* results,
        const std::atomic_bool* halt,
        const ExpressionOrbitPlan* plan) const {
    if (!valid() || !contexts || !results ||
        lanes < 1 || lanes > 4 || mxit < 1 ||
        !(bailout > 0.0) || !std::isfinite(bailout))
        return false;
    if ((!_impl->planKey.empty() &&
         (!plan || !supports(*plan))) ||
        (_impl->planKey.empty() && plan))
        return false;

    ExpressionJitInput4 input;
    ExpressionJitOutput4 output;
    std::unique_ptr<ExpressionJitInvariantInput4> invariantInput;
    std::unique_ptr<ExpressionOrbitPlan::Prepared[]> prepared;
    input.setContexts(contexts);
    bool active[4] = { false, false, false, false };
    uint8_t activeMask = 0;
    int activeCount = 0;
    for (int lane = 0; lane < 4; ++lane) {
        results[lane] = -2.0f;
        if (lane >= lanes) continue;
        double re = input.vectors[ExpressionJitInput4::Z_RE][lane];
        double im = input.vectors[ExpressionJitInput4::Z_IM][lane];
        if (!std::isfinite(re) || !std::isfinite(im) ||
            std::hypot(re, im) > bailout) {
            results[lane] = 0.0f;
        } else {
            active[lane] = true;
            activeMask |= static_cast<uint8_t>(1u << lane);
            ++activeCount;
        }
    }
    if (plan) {
        invariantInput =
            std::make_unique<ExpressionJitInvariantInput4>();
        prepared =
            std::make_unique<ExpressionOrbitPlan::Prepared[]>(4);
        if (!plan->prepare4(
                contexts, activeMask, prepared.get()))
            return false;
        for (int lane = 0; lane < 4; ++lane)
            invariantInput->setPreparedLane(
                lane, *plan, prepared[(size_t)lane]);
    }

    for (int iteration = 0;
         iteration < mxit && activeCount > 0; ++iteration) {
        if ((iteration & 255) == 0 && halt && *halt)
            return false;
        for (int lane = 0; lane < 4; ++lane)
            input.vectors[ExpressionJitInput4::N_RE][lane] =
                (double)iteration;
        evaluate(input, invariantInput.get(), output);
        for (int lane = 0; lane < lanes; ++lane) {
            if (!active[lane]) continue;
            double re = output.re[lane];
            double im = output.im[lane];
            input.vectors[ExpressionJitInput4::Z_RE][lane] = re;
            input.vectors[ExpressionJitInput4::Z_IM][lane] = im;
            if (!std::isfinite(re) || !std::isfinite(im) ||
                std::hypot(re, im) > bailout) {
                results[lane] = (float)(iteration + 1);
                active[lane] = false;
                --activeCount;
            }
        }
    }
    return true;
}

} // namespace formula
