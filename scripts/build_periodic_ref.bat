@echo off
setlocal
cd /d "%~dp0.."
if not exist build mkdir build
call "%~dp0msvcenv.bat" || goto :err

cl /nologo /O2 /EHsc /std:c++17 /arch:AVX2 /MT ^
   /I "%VCPKG%\include" /I src\engine /Fo:build\ ^
   src\bench\periodic_ref_bench.cpp ^
   /Fe:build\periodic_ref_bench.exe ^
   /link /LIBPATH:"%VCPKG%\lib" gmp.lib || goto :err

echo Build OK: build\periodic_ref_bench.exe
exit /b 0
:err
echo BUILD FAILED
exit /b 1
