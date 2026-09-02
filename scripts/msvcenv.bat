@echo off
REM Set up an MSVC x64 build environment and locate vcpkg + the OpenMP redist

if defined WindowsSdkDir goto :findvcpkg
if defined VCINSTALLDIR goto :findvcpkg

if exist "%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" (
  for /f "usebackq delims=" %%i in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -property installationPath 2^>nul`) do (
    if exist "%%i\VC\Auxiliary\Build\vcvarsall.bat" (
      call "%%i\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul
    ) else if exist "%%i\VC\Auxiliary\Build\vcvars64.bat" (
      call "%%i\VC\Auxiliary\Build\vcvars64.bat" >nul
    )
  )
)
if defined WindowsSdkDir goto :findvcpkg
if defined VCINSTALLDIR goto :findvcpkg

REM ---- fallback: manual search ---------------------------------------------
set "VSROOT="
set "MSVCVER="
for %%P in ("%ProgramFiles%\Microsoft Visual Studio" "%ProgramFiles(x86)%\Microsoft Visual Studio") do (
  for /f "delims=" %%y in ('dir /b /ad /o-n "%%~P" 2^>nul') do (
    for %%e in (Enterprise Professional Community BuildTools Preview) do (
      if not defined MSVCVER if exist "%%~P\%%y\%%e\VC\Tools\MSVC" (
        for /f "delims=" %%d in ('dir /b /ad /o-n "%%~P\%%y\%%e\VC\Tools\MSVC" 2^>nul') do (
          if not defined MSVCVER if exist "%%~P\%%y\%%e\VC\Tools\MSVC\%%d\bin\HostX64\x64\cl.exe" (
            set "VSROOT=%%~P\%%y\%%e"
            set "MSVCVER=%%d"
          )
        )
      )
    )
  )
)
if not defined MSVCVER ( echo msvcenv: no MSVC toolset with cl.exe found & exit /b 1 )
set "MSVC=%VSROOT%\VC\Tools\MSVC\%MSVCVER%"

set "SDKROOT=%ProgramFiles(x86)%\Windows Kits\10"
set "SDKVER="
for /f "delims=" %%d in ('dir /b /ad /o-n "%SDKROOT%\Include" 2^>nul') do (
  if not defined SDKVER if exist "%SDKROOT%\Lib\%%d\um\x64\kernel32.lib" set "SDKVER=%%d"
)
if not defined SDKVER ( echo msvcenv: no usable Windows SDK found & exit /b 1 )

set "INCLUDE=%MSVC%\include;%SDKROOT%\Include\%SDKVER%\ucrt;%SDKROOT%\Include\%SDKVER%\shared;%SDKROOT%\Include\%SDKVER%\um;%SDKROOT%\Include\%SDKVER%\winrt"
set "LIB=%MSVC%\lib\x64;%SDKROOT%\Lib\%SDKVER%\ucrt\x64;%SDKROOT%\Lib\%SDKVER%\um\x64"
set "PATH=%MSVC%\bin\HostX64\x64;%SDKROOT%\bin\%SDKVER%\x64;%PATH%"

:findvcpkg
REM ---- vcpkg static triplet (VCPKG_ROOT, then PATH, then common locations) --
set "VCPKG="
if defined VCPKG_ROOT call :trypkg "%VCPKG_ROOT%"
if not defined VCPKG for /f "delims=" %%i in ('where vcpkg 2^>nul') do if not defined VCPKG call :trypkg "%%~dpi"
if not defined VCPKG for %%R in (
  "%USERPROFILE%\repo\vcpkg" "%USERPROFILE%\vcpkg" "%USERPROFILE%\source\repos\vcpkg"
  "C:\vcpkg" "C:\src\vcpkg" "C:\dev\vcpkg"
) do if not defined VCPKG call :trypkg %%R
if not defined VCPKG echo msvcenv: WARNING vcpkg (x64-windows-static) not found - set VCPKG_ROOT if a target needs GMP

echo msvcenv: MSVC environment ready
if defined VCPKG echo msvcenv: vcpkg %VCPKG%
exit /b 0

:trypkg
set "_r=%~1"
if "%_r:~-1%"=="\" set "_r=%_r:~0,-1%"
if exist "%_r%\installed\x64-windows-static\lib\gmp.lib" set "VCPKG=%_r%\installed\x64-windows-static"
exit /b 0
