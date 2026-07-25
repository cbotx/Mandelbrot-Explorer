@echo off
REM Set up an MSVC x64 build environment WITHOUT vcvars64.bat.
REM Needed because on this machine vcvars is broken: vswhere does not index the
REM VS 18 ("2026") preview, and the default toolset can be mid-install (no cl.exe).
REM This auto-selects the newest installed MSVC toolset that actually has cl.exe,
REM and the newest Windows 10/11 SDK that has the core libs.
set "VCROOT=C:\Program Files\Microsoft Visual Studio\18\Enterprise\VC\Tools\MSVC"
set "MSVCVER="
for /f "delims=" %%d in ('dir /b /ad /o-n "%VCROOT%" 2^>nul') do (
  if not defined MSVCVER if exist "%VCROOT%\%%d\bin\HostX64\x64\cl.exe" set "MSVCVER=%%d"
)
set "SDKROOT=C:\Program Files (x86)\Windows Kits\10"
set "SDKVER="
for /f "delims=" %%d in ('dir /b /ad /o-n "%SDKROOT%\Include" 2^>nul') do (
  if not defined SDKVER if exist "%SDKROOT%\Lib\%%d\um\x64\kernel32.lib" set "SDKVER=%%d"
)
if not defined MSVCVER ( echo msvcenv: no MSVC toolset with cl.exe found & exit /b 1 )
if not defined SDKVER  ( echo msvcenv: no usable Windows SDK found & exit /b 1 )
set "MSVC=%VCROOT%\%MSVCVER%"
set "INCLUDE=%MSVC%\include;%SDKROOT%\Include\%SDKVER%\ucrt;%SDKROOT%\Include\%SDKVER%\shared;%SDKROOT%\Include\%SDKVER%\um;%SDKROOT%\Include\%SDKVER%\winrt"
set "LIB=%MSVC%\lib\x64;%SDKROOT%\Lib\%SDKVER%\ucrt\x64;%SDKROOT%\Lib\%SDKVER%\um\x64"
set "PATH=%MSVC%\bin\HostX64\x64;%PATH%"
echo msvcenv: MSVC %MSVCVER%, SDK %SDKVER%
exit /b 0
