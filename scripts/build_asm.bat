@echo off
setlocal
cd /d "%~dp0.."
if not exist build mkdir build
call "%~dp0msvcenv.bat" || goto :err
ml64 /nologo /c /Fo build\mm_mul.obj src\engine\mm_mul.asm || goto :err
cl /nologo /O2 /EHsc /std:c++17 /arch:AVX2 /MT ^
   /I "%VCPKG%\include" /I src\engine /Fo:build\ src\bench\bench_asm.cpp build\mm_mul.obj ^
   /Fe:build\bench_asm.exe /link /LIBPATH:"%VCPKG%\lib" gmp.lib || goto :err
echo Build OK: build\bench_asm.exe
exit /b 0
:err
echo BUILD FAILED
exit /b 1
