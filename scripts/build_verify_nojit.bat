@echo off
REM Build the Defender-safe non-JIT verification harness.
setlocal
cd /d "%~dp0.."
if not exist build mkdir build
call "%~dp0msvcenv.bat" || goto :err

cl /nologo /O2 /EHsc /openmp /std:c++17 /arch:AVX2 /MT /DMANDEL_VERIFY_NO_JIT ^
   /I "%VCPKG%\include" /I src\engine /I src\gui ^
   /Fo:build\ ^
   src\tools\verify.cpp src\gui\mandel_navigator.cpp src\gui\navigator.cpp src\gui\orbit_overlay.cpp src\gui\orbit_thumbnail.cpp src\gui\gui_color.cpp ^
   src\engine\compute_backend.cpp src\engine\compute_backend_d3d11.cpp src\engine\mandel_perturbation.cpp src\engine\mandel_expression.cpp src\engine\float_math.cpp src\engine\formula_expression.cpp src\engine\formula_expression_centered.cpp src\engine\formula_expression_orbit.cpp src\engine\formula_reference_orbit.cpp src\engine\formula_scaled_residual.cpp src\engine\formula_expression_jit_stub.cpp src\engine\formula_expression_oracle.cpp src\engine\interpolate.cpp src\engine\color.cpp ^
   /Fe:build\verify_nojit.exe ^
   /link /LIBPATH:"%VCPKG%\lib" mpfr.lib gmp.lib d3d11.lib dxgi.lib d3dcompiler.lib || goto :err

echo Build OK: build\verify_nojit.exe
exit /b 0
:err
echo BUILD FAILED
exit /b 1
