@echo off
setlocal
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0format_cpp.ps1" -Check
exit /b %errorlevel%
