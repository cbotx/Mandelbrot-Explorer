#ifndef MANDEL_FORMULA_SPEC_H
#define MANDEL_FORMULA_SPEC_H

#include <cstdint>

enum class FormulaId : uint8_t {
    PowerPlusC
};

enum class FormulaParameter : uint8_t {
    InitialZ,
    C,
    Power
};

enum FormulaCapability : uint32_t {
    FORMULA_HOLOMORPHIC       = 1u << 0,
    FORMULA_POLYNOMIAL        = 1u << 1,
    FORMULA_PERTURBATION      = 1u << 2,
    FORMULA_SERIES_APPROX     = 1u << 3,
    FORMULA_BLA               = 1u << 4,
    FORMULA_PERIODIC_REF      = 1u << 5,
    FORMULA_DISTANCE_ESTIMATE = 1u << 6,
    FORMULA_FLOATEXP          = 1u << 7,
    FORMULA_BIGFIXED          = 1u << 8
};

struct FormulaSpec {
    FormulaId id = FormulaId::PowerPlusC;
    int power = 2;
    uint32_t capabilities =
        FORMULA_HOLOMORPHIC | FORMULA_POLYNOMIAL | FORMULA_PERTURBATION |
        FORMULA_SERIES_APPROX | FORMULA_BLA | FORMULA_PERIODIC_REF |
        FORMULA_DISTANCE_ESTIMATE | FORMULA_FLOATEXP | FORMULA_BIGFIXED;

    constexpr bool supports(FormulaCapability capability) const {
        return (capabilities & (uint32_t)capability) != 0;
    }
};

struct RenderSlice {
    FormulaParameter pixel = FormulaParameter::C;
};

struct FormulaContext {
    FormulaSpec formula{};
    RenderSlice slice{};
    uint32_t enabledCapabilities = formula.capabilities;

    constexpr bool supports(FormulaCapability capability) const {
        return formula.supports(capability) &&
               (enabledCapabilities & (uint32_t)capability) != 0;
    }
};

template <class Complex>
struct OrbitInput {
    Complex z0;
    Complex c;
};

constexpr FormulaContext quadraticMandelbrot() {
    return FormulaContext{};
}

constexpr FormulaContext quadraticJulia() {
    FormulaContext context{};
    context.slice.pixel = FormulaParameter::InitialZ;
    context.enabledCapabilities = FORMULA_DISTANCE_ESTIMATE;
    return context;
}

static_assert(quadraticMandelbrot().supports(FORMULA_BLA),
              "Quadratic Mandelbrot must retain its optimized kernel");
static_assert(!quadraticJulia().supports(FORMULA_BLA),
              "Julia accelerations are enabled only after their recurrences are implemented");

#endif
