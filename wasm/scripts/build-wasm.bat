@echo off
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0build-wasm.ps1" %*
exit /b %errorlevel%
