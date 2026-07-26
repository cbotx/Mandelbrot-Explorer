@echo off
setlocal
cd /d "%~dp0"
set VCPKG=C:\Users\yuxiangchi\repo\vcpkg\installed\x64-windows-static
if not exist build mkdir build
call "%~dp0msvcenv.bat" || goto :err
cl /nologo /O2 /EHsc /std:c++17 /arch:AVX2 /MT ^
   /I "%VCPKG%\include" /Fo:build\ bench_refbuild.cpp ^
   /Fe:build\bench_refbuild.exe /link /LIBPATH:"%VCPKG%\lib" gmp.lib || goto :err
echo Build OK: build\bench_refbuild.exe
exit /b 0
:err
echo BUILD FAILED
exit /b 1
