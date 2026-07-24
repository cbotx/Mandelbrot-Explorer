@echo off
setlocal
cd /d "%~dp0"
set VCVARS="C:\Program Files\Microsoft Visual Studio\18\Enterprise\VC\Auxiliary\Build\vcvars64.bat"
set VCPKG=C:\Users\yuxiangchi\repo\vcpkg\installed\x64-windows
call %VCVARS% >nul || goto :err

cl /nologo /O2 /EHsc /std:c++17 ^
   /I "%VCPKG%\include" ^
   pbench.cpp ^
   /Fe:pbench.exe ^
   /link /LIBPATH:"%VCPKG%\lib" gmp.lib || goto :err

copy /y "%VCPKG%\bin\gmp-10.dll" . >nul
echo Build OK: pbench.exe
exit /b 0
:err
echo BUILD FAILED
exit /b 1
