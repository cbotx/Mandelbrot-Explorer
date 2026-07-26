@echo off
setlocal
cd /d "%~dp0.."
if not exist build mkdir build
call "%~dp0msvcenv.bat" || goto :err

cl /nologo /O2 /EHsc /std:c++17 /MT ^
   /I "%VCPKG%\include" /I src\engine ^
   /Fo:build\ ^
   src\bench\bench_bigfixed.cpp ^
   /Fe:build\bench_bigfixed.exe ^
   /link /LIBPATH:"%VCPKG%\lib" gmp.lib || goto :err

echo Build OK: build\bench_bigfixed.exe
exit /b 0
:err
echo BUILD FAILED
exit /b 1
