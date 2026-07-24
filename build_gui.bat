@echo off
setlocal
cd /d "%~dp0"
set VCVARS="C:\Program Files\Microsoft Visual Studio\18\Enterprise\VC\Auxiliary\Build\vcvars64.bat"
set VCPKG=C:\Users\yuxiangchi\repo\vcpkg\installed\x64-windows-static
if not exist build mkdir build
call %VCVARS% >nul || goto :err

REM Static GMP (x64-windows-static, /MT) -> single-file exe, no gmp-10.dll.
cl /nologo /O2 /EHsc /openmp /std:c++17 /arch:AVX2 /MT ^
   /DUNICODE /D_UNICODE ^
   /I "%VCPKG%\include" ^
   /Fo:build\ ^
   win32_main.cpp mandel_navigator.cpp navigator.cpp ^
   mandel_perturbation.cpp float_math.cpp interpolate.cpp color.cpp ^
   /Fe:build\mandel_gui.exe ^
   /link /SUBSYSTEM:WINDOWS /LIBPATH:"%VCPKG%\lib" ^
   gmp.lib user32.lib gdi32.lib comctl32.lib comdlg32.lib || goto :err

REM Ship the OpenMP runtime next to the exe so build\ is copy-and-run portable.
for /f "delims=" %%D in ('dir /b /s "C:\Program Files\Microsoft Visual Studio\18\Enterprise\VC\Redist\*\x64\*\vcomp140.dll" 2^>nul') do copy /y "%%D" build\ >nul
echo Build OK: build\mandel_gui.exe  (deps: build\vcomp140.dll + Windows system DLLs)
exit /b 0
:err
echo BUILD FAILED
exit /b 1
