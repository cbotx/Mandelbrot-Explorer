@echo off
setlocal
cd /d "%~dp0.."
call "%~dp0msvcenv.bat" || goto :err

for %%I in ("%VCPKG%\..\x64-windows") do set "VCPKG_RELEASE=%%~fI"
if not exist "%VCPKG_RELEASE%\include\mpfr.h" (
  echo Missing dynamic MPFR package. Install mpfr:x64-windows.
  goto :err
)
if not exist "%VCPKG_RELEASE%\lib\asmjit.lib" (
  echo Missing dynamic-triplet AsmJit package. Install asmjit:x64-windows.
  goto :err
)
if not exist "%VCPKG_RELEASE%\bin\asmjit.dll" (
  echo Missing AsmJit runtime DLL. Install asmjit:x64-windows.
  goto :err
)
if not exist "%VCPKG_RELEASE%\bin\gmp-10.dll" (
  echo Missing dynamic GMP package. Install gmp:x64-windows.
  goto :err
)

if not exist build\release mkdir build\release
if exist build\release\obj rmdir /s /q build\release\obj
del /q build\release\mandel_gui.exe build\release\mandel_gui.res build\release\*.dll 2>nul
mkdir build\release\obj

rc /nologo /fo build\release\mandel_gui.res src\gui\mandel_gui.rc || goto :err

cl /nologo /O2 /GL /EHsc /openmp /std:c++17 /arch:AVX2 /MD /DMANDEL_ENABLE_ASMJIT ^
   /DUNICODE /D_UNICODE ^
   /I "%VCPKG_RELEASE%\include" /I src\engine /I src\gui ^
   /Fo:build\release\obj\ ^
   src\gui\win32_main.cpp src\gui\gui_color.cpp src\gui\gui_theme.cpp src\gui\gui_export.cpp src\gui\ui_framework.cpp src\gui\ui_controls.cpp src\gui\formula_editor_accessibility.cpp src\gui\formula_editor_panel.cpp src\gui\mandel_navigator.cpp src\gui\navigator.cpp src\gui\orbit_overlay.cpp src\gui\orbit_thumbnail.cpp ^
   src\engine\compute_backend.cpp src\engine\compute_backend_d3d11.cpp src\engine\mandel_perturbation.cpp src\engine\mandel_expression.cpp src\engine\formula_expression.cpp src\engine\formula_expression_centered.cpp src\engine\formula_expression_orbit.cpp src\engine\formula_reference_orbit.cpp src\engine\formula_scaled_residual.cpp src\engine\formula_taylor_jet.cpp src\engine\formula_deep_renderer.cpp src\engine\formula_expression_oracle.cpp src\engine\formula_expression_jit.cpp src\engine\float_math.cpp src\engine\interpolate.cpp src\engine\color.cpp ^
   build\release\mandel_gui.res ^
   /Fe:build\release\mandel_gui.exe ^
   /link /LTCG /OPT:REF /OPT:ICF /DYNAMICBASE /NXCOMPAT /HIGHENTROPYVA /MANIFEST:NO /SUBSYSTEM:WINDOWS ^
   /LIBPATH:"%VCPKG_RELEASE%\lib" ^
   asmjit.lib mpfr.lib gmp.lib d3d11.lib dxgi.lib d3dcompiler.lib user32.lib gdi32.lib comctl32.lib comdlg32.lib windowscodecs.lib ole32.lib oleaut32.lib oleacc.lib || goto :err

copy /y "%VCPKG_RELEASE%\bin\gmp-10.dll" build\release\ >nul || goto :err
copy /y "%VCPKG_RELEASE%\bin\asmjit.dll" build\release\ >nul || goto :err
for /f "delims=" %%D in ('dir /b "%VCPKG_RELEASE%\bin\mpfr-*.dll" 2^>nul') do copy /y "%VCPKG_RELEASE%\bin\%%D" build\release\ >nul
if not exist build\release\mpfr-*.dll (
  echo MPFR runtime DLL was not copied.
  goto :err
)

set "REDISTVER="
for /f "delims=" %%D in ('dir /b /ad /o-n "%VCREDIST%" 2^>nul') do if not defined REDISTVER if exist "%VCREDIST%\%%D\x64\Microsoft.VC143.CRT" set "REDISTVER=%%D"
if not defined REDISTVER (
  echo Visual C++ redistributable files were not found.
  goto :err
)
copy /y "%VCREDIST%\%REDISTVER%\x64\Microsoft.VC143.CRT\*.dll" build\release\ >nul || goto :err
copy /y "%VCREDIST%\%REDISTVER%\x64\Microsoft.VC143.OpenMP\*.dll" build\release\ >nul || goto :err

echo Build OK: build\release\mandel_gui.exe
exit /b 0
:err
echo RELEASE BUILD FAILED
exit /b 1
