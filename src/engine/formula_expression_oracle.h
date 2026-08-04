#ifndef MANDEL_FORMULA_EXPRESSION_ORACLE_H
#define MANDEL_FORMULA_EXPRESSION_ORACLE_H

#include <mpfr.h>

#include <array>
#include <complex>
#include <string>

#include "formula_expression.h"

namespace formula {

class MpfrComplex {
public:
    explicit MpfrComplex(mpfr_prec_t precision = 256);
    MpfrComplex(const MpfrComplex& other);
    MpfrComplex& operator=(const MpfrComplex& other);
    ~MpfrComplex();

    void set(double real, double imaginary = 0.0);
    bool set(const std::string& real, const std::string& imaginary = "0");
    void set(const MpfrComplex& other);
    std::complex<double> toDouble() const;
    mpfr_prec_t precision() const;

    mpfr_t re;
    mpfr_t im;
};

struct ExpressionOracleContext {
    explicit ExpressionOracleContext(mpfr_prec_t precision = 256);

    MpfrComplex z;
    MpfrComplex c;
    MpfrComplex z0;
    std::array<MpfrComplex, 8> parameters;
    int iteration = 0;
};

class ExpressionOracle {
public:
    static bool evaluate(const ExpressionProgram& program,
                         const ExpressionOracleContext& context,
                         MpfrComplex& output,
                         std::string* error = nullptr);
};

} // namespace formula

#endif
