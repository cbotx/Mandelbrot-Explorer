@echo off
setlocal
cd /d "%~dp0.."
if not exist build mkdir build
call "%~dp0msvcenv.bat" || goto :err
set "LLVM=C:\Program Files\Microsoft Visual Studio\18\Enterprise\VC\Tools\Llvm\x64\bin"
"%LLVM%\clang-cl.exe" /O2 /EHsc /std:c++17 /clang:-march=native /MT ^
   /I "%VCPKG%\include" /I src\engine /Fo:build\ src\bench\bench_mulhigh.cpp ^
   /Fe:build\bench_mulhigh_clang.exe /link /LIBPATH:"%VCPKG%\lib" gmp.lib || goto :err
echo Build OK: build\bench_mulhigh_clang.exe
exit /b 0
:err
echo BUILD FAILED
exit /b 1
