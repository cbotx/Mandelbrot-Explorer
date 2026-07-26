@echo off
setlocal
cd /d "%~dp0.."
if not exist build mkdir build
call "%~dp0msvcenv.bat" || goto :err

cl /nologo /O2 /EHsc /openmp /std:c++17 /arch:AVX2 /MT ^
   /I "%VCPKG%\include" /I src\engine ^
   /Fo:build\ ^
   src\bench\palette_preview.cpp src\engine\mandel_perturbation.cpp src\engine\float_math.cpp src\engine\interpolate.cpp src\engine\color.cpp ^
   /Fe:build\palette_preview.exe ^
   /link /LIBPATH:"%VCPKG%\lib" gmp.lib || goto :err

echo Build OK: build\palette_preview.exe
exit /b 0
:err
echo BUILD FAILED
exit /b 1
