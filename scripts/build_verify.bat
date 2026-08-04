@echo off
REM Build the headless verification/benchmark harness (MSVC + static GMP).
setlocal
cd /d "%~dp0.."
if not exist build mkdir build
call "%~dp0msvcenv.bat" || goto :err

cl /nologo /O2 /EHsc /openmp /std:c++17 /arch:AVX2 /MT /DMANDEL_ENABLE_ASMJIT ^
   /I "%VCPKG%\include" /I src\engine ^
   /Fo:build\ ^
   src\tools\verify.cpp src\engine\compute_backend.cpp src\engine\compute_backend_d3d11.cpp src\engine\mandel_perturbation.cpp src\engine\mandel_expression.cpp src\engine\float_math.cpp src\engine\formula_expression.cpp src\engine\formula_expression_jit.cpp src\engine\formula_expression_oracle.cpp ^
   /Fe:build\verify.exe ^
   /link /LIBPATH:"%VCPKG%\lib" asmjit.lib mpfr.lib gmp.lib d3d11.lib dxgi.lib d3dcompiler.lib || goto :err

echo Build OK: build\verify.exe
exit /b 0
:err
echo BUILD FAILED
exit /b 1
