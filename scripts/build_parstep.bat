@echo off
setlocal
cd /d "%~dp0.."
if not exist build mkdir build
call "%~dp0msvcenv.bat" || goto :err
cl /nologo /O2 /EHsc /std:c++17 /arch:AVX2 /openmp /MT ^
   /I "%VCPKG%\include" /I src\engine /Fo:build\ src\bench\bench_parstep.cpp ^
   /Fe:build\bench_parstep.exe /link /LIBPATH:"%VCPKG%\lib" gmp.lib || goto :err
echo Build OK: build\bench_parstep.exe
exit /b 0
:err
echo BUILD FAILED
exit /b 1
