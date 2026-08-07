@echo off
setlocal
cd /d "%~dp0.."
if not exist build mkdir build
call "%~dp0msvcenv.bat" || goto :err

REM Static GMP (x64-windows-static, /MT) -> single-file exe, no gmp-10.dll.
cl /nologo /O2 /EHsc /openmp /std:c++17 /arch:AVX2 /MT /DMANDEL_ENABLE_ASMJIT ^
   /DUNICODE /D_UNICODE ^
   /I "%VCPKG%\include" /I src\engine /I src\gui ^
   /Fo:build\ ^
   src\gui\win32_main.cpp src\gui\gui_color.cpp src\gui\gui_theme.cpp src\gui\gui_export.cpp src\gui\ui_framework.cpp src\gui\ui_controls.cpp src\gui\formula_editor_accessibility.cpp src\gui\formula_editor_panel.cpp src\gui\mandel_navigator.cpp src\gui\navigator.cpp src\gui\orbit_overlay.cpp ^
   src\engine\compute_backend.cpp src\engine\compute_backend_d3d11.cpp src\engine\mandel_perturbation.cpp src\engine\mandel_expression.cpp src\engine\formula_expression.cpp src\engine\formula_expression_centered.cpp src\engine\formula_expression_orbit.cpp src\engine\formula_expression_jit.cpp src\engine\float_math.cpp src\engine\interpolate.cpp src\engine\color.cpp ^
   /Fe:build\mandel_gui.exe ^
   /link /SUBSYSTEM:WINDOWS /LIBPATH:"%VCPKG%\lib" ^
   asmjit.lib gmp.lib d3d11.lib dxgi.lib d3dcompiler.lib user32.lib gdi32.lib comctl32.lib comdlg32.lib windowscodecs.lib ole32.lib oleaut32.lib oleacc.lib || goto :err

REM Ship the OpenMP runtime next to the exe so build\ is copy-and-run portable.
for /f "delims=" %%D in ('dir /b /s "%VCREDIST%\*\x64\*\vcomp140.dll" 2^>nul') do copy /y "%%D" build\ >nul
echo Build OK: build\mandel_gui.exe  (deps: build\vcomp140.dll + Windows system DLLs)
exit /b 0
:err
echo BUILD FAILED
exit /b 1
