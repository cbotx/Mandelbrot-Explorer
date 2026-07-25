@echo off
REM Build the headless verification/benchmark harness (MSVC + static GMP).
setlocal
cd /d "%~dp0"
set VCPKG=C:\Users\yuxiangchi\repo\vcpkg\installed\x64-windows-static
if not exist build mkdir build
call "%~dp0msvcenv.bat" || goto :err

cl /nologo /O2 /EHsc /openmp /std:c++17 /arch:AVX2 /MT ^
   /I "%VCPKG%\include" ^
   /Fo:build\ ^
   verify.cpp mandel_perturbation.cpp float_math.cpp ^
   /Fe:build\verify.exe ^
   /link /LIBPATH:"%VCPKG%\lib" gmp.lib || goto :err

echo Build OK: build\verify.exe
exit /b 0
:err
echo BUILD FAILED
exit /b 1
