@echo off
setlocal
cd /d "%~dp0.."
if not exist build mkdir build
call "%~dp0msvcenv.bat" || goto :err
cl /nologo /O2 /EHsc /openmp /std:c++17 /arch:AVX2 /MT ^
   /I src\engine /Fo:build\ src\bench\coloring_demo.cpp src\engine\interpolate.cpp src\engine\color.cpp ^
   /Fe:build\coloring_demo.exe || goto :err
echo Build OK: build\coloring_demo.exe
exit /b 0
:err
echo BUILD FAILED
exit /b 1
