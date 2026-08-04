#ifndef MANDEL_FORMULA_DIALOG_H
#define MANDEL_FORMULA_DIALOG_H

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

#include <array>
#include <complex>
#include <string>

#include "formula_spec.h"

struct FormulaDialogConfig {
    std::string source = "z*z+c";
    FormulaParameter pixelParameter = FormulaParameter::C;
    std::complex<double> fixedZ0{};
    std::complex<double> fixedC{};
    std::array<std::complex<double>, 8> parameters{};
    double bailout = 4.0;
};

#endif
