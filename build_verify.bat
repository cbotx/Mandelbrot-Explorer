@echo off
REM Build the headless verification/benchmark harness with MSVC + GMP (vcpkg).
setlocal
cd /d "%~dp0"
set VCVARS="C:\Program Files\Microsoft Visual Studio\18\Enterprise\VC\Auxiliary\Build\vcvars64.bat"
set VCPKG=C:\Users\yuxiangchi\repo\vcpkg\installed\x64-windows
call %VCVARS% >nul || goto :err

cl /nologo /O2 /EHsc /openmp /std:c++17 /arch:AVX2 ^
   /I "%VCPKG%\include" ^
   verify.cpp mandel_perturbation.cpp float_math.cpp ^
   /Fe:verify.exe ^
   /link /LIBPATH:"%VCPKG%\lib" gmp.lib || goto :err

REM Make the GMP DLL available next to the exe.
copy /y "%VCPKG%\bin\gmp-10.dll" . >nul
echo Build OK: verify.exe
exit /b 0
:err
echo BUILD FAILED
exit /b 1
