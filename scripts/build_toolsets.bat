@echo off
setlocal

set "CONFIG=Release"
set "BUILD_DIR=build\ninja-msvc"
set "ZSPACE_CORE_DIR=%~1"
set "CMAKE_EXE=cmake"
set "BUILD_TOOL_ARGS=/m /clp:ErrorsOnly"

if "%ZSPACE_CORE_DIR%"=="" (
    set "ZSPACE_CORE_DIR=%~dp0..\..\zspace_core"
)

where cmake >nul 2>nul
if errorlevel 1 (
    if exist "C:\Program Files\Microsoft Visual Studio\2022\Professional\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" (
        set "CMAKE_EXE=C:\Program Files\Microsoft Visual Studio\2022\Professional\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
    ) else if exist "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" (
        set "CMAKE_EXE=C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
    ) else if exist "C:\Program Files\Microsoft Visual Studio\2022\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" (
        set "CMAKE_EXE=C:\Program Files\Microsoft Visual Studio\2022\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
    ) else (
        echo.
        echo [zspace_toolsets] CMake was not found on PATH or in Visual Studio 2022.
        goto :fail
    )
)

echo [zspace_toolsets] Building with zspace_core source
echo [zspace_toolsets] zspace_core: "%ZSPACE_CORE_DIR%"
echo [zspace_toolsets] cmake: "%CMAKE_EXE%"
echo.

pushd "%~dp0.."
if errorlevel 1 goto :fail

if not exist "%BUILD_DIR%" mkdir "%BUILD_DIR%"

"%CMAKE_EXE%" -S . -B "%BUILD_DIR%" -DCMAKE_BUILD_TYPE=%CONFIG% -DZSPACE_TOOLSETS_USE_ZSPACE_SOURCE=ON -DZSPACE_CORE_DIR="%ZSPACE_CORE_DIR%" -DZSPACE_TOOLSETS_BUILD_TESTS=ON
if errorlevel 1 goto :fail_pop

"%CMAKE_EXE%" --build "%BUILD_DIR%" --config %CONFIG% -- %BUILD_TOOL_ARGS%
if errorlevel 1 goto :fail_pop

echo.
echo [zspace_toolsets] build finished successfully.
popd
if not defined ZSPACE_TOOLSETS_NO_PAUSE pause
exit /b 0

:fail_pop
popd

:fail
echo.
echo [zspace_toolsets] build failed.
if not defined ZSPACE_TOOLSETS_NO_PAUSE pause
exit /b 1
