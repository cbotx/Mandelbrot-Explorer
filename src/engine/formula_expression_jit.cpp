#include "formula_expression_jit.h"

#include <asmjit/x86.h>

#include <algorithm>
#include <cstring>

namespace formula {

using JitFunction = void (*)(const double*, double*);

struct ExpressionJit4::Impl {
    asmjit::JitAllocator::CreateParams params;
    asmjit::JitRuntime runtime;
    JitFunction function = nullptr;

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
    for (int lane = 0; lane < 4; ++lane) {
        vectors[Z_RE][lane] = contexts[lane].z.real();
        vectors[Z_IM][lane] = contexts[lane].z.imag();
        vectors[C_RE][lane] = contexts[lane].c.real();
        vectors[C_IM][lane] = contexts[lane].c.imag();
        vectors[Z0_RE][lane] = contexts[lane].z0.real();
        vectors[Z0_IM][lane] = contexts[lane].z0.imag();
        vectors[N_RE][lane] = (double)contexts[lane].iteration;
        vectors[N_IM][lane] = 0.0;
        for (int p = 0; p < 8; ++p) {
            vectors[PARAM_RE + 2 * p][lane] = contexts[lane].parameters[p].real();
            vectors[PARAM_RE + 2 * p + 1][lane] = contexts[lane].parameters[p].imag();
        }
    }
}

ExpressionJit4::ExpressionJit4() : _impl(std::make_unique<Impl>()) {}
ExpressionJit4::~ExpressionJit4() = default;
ExpressionJit4::ExpressionJit4(ExpressionJit4&&) noexcept = default;
ExpressionJit4& ExpressionJit4::operator=(ExpressionJit4&&) noexcept = default;

bool ExpressionJit4::valid() const {
    return _impl && _impl->function;
}

bool ExpressionJit4::usesDualMapping() const {
    return _impl && _impl->runtime.allocator()->hasOption(
        asmjit::JitAllocatorOptions::kUseDualMapping);
}

const void* ExpressionJit4::codeAddress() const {
    return _impl ? reinterpret_cast<const void*>(_impl->function) : nullptr;
}

bool ExpressionJit4::compile(const ExpressionProgram& program, std::string* error) {
    if (error) error->clear();
    if (!_impl) _impl = std::make_unique<Impl>();
    if (_impl->function) {
        _impl->runtime.release(_impl->function);
        _impl->function = nullptr;
    }
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
    a.vmovupd(ptr(rdx, 0), yr0);
    a.vmovupd(ptr(rdx, 32), yi0);
    a.vzeroupper();
    a.add(rsp, frameSize);
    a.ret();
    err = _impl->runtime.add(&_impl->function, &code);
    if (err) {
        if (error) *error = asmjit::DebugUtils::errorAsString(err);
        _impl->function = nullptr;
        return false;
    }
    return true;
}

void ExpressionJit4::reset() {
    if (_impl && _impl->function) {
        _impl->runtime.release(_impl->function);
        _impl->function = nullptr;
    }
}

void ExpressionJit4::evaluate(const ExpressionJitInput4& input,
                              ExpressionJitOutput4& output) const {
    if (_impl && _impl->function)
        _impl->function(&input.vectors[0][0], &output.re[0]);
}

bool ExpressionJit4::evaluate(const ExpressionContext* contexts,
                              Complex* outputs) const {
    if (!valid() || !contexts || !outputs) return false;
    ExpressionJitInput4 input;
    ExpressionJitOutput4 output;
    input.setContexts(contexts);
    evaluate(input, output);
    for (int lane = 0; lane < 4; ++lane)
        outputs[lane] = { output.re[lane], output.im[lane] };
    return true;
}

} // namespace formula
