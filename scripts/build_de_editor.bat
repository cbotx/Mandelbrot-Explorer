@echo off
setlocal
cd /d "%~dp0.."
if not exist build mkdir build
call "%~dp0msvcenv.bat" || goto :err

REM Static GMP (x64-windows-static, /MT) -> single-file exe, no gmp-10.dll.
cl /nologo /O2 /EHsc /openmp /std:c++17 /arch:AVX2 /MT ^
   /DUNICODE /D_UNICODE ^
   /I "%VCPKG%\include" /I src\engine ^
   /Fo:build\ ^
   src\tools\de_editor.cpp src\engine\mandel_perturbation.cpp src\engine\float_math.cpp src\engine\interpolate.cpp ^
   /Fe:build\de_editor.exe ^
   /link /SUBSYSTEM:WINDOWS /LIBPATH:"%VCPKG%\lib" ^
   gmp.lib user32.lib gdi32.lib comctl32.lib || goto :err

for /f "delims=" %%D in ('dir /b /s "%VCREDIST%\*\x64\*\vcomp140.dll" 2^>nul') do copy /y "%%D" build\ >nul
echo Build OK: build\de_editor.exe
exit /b 0
:err
echo BUILD FAILED
exit /b 1
