@echo off
REM run_qml.bat -- Launch the CSHORE.jl QML GUI on Windows.
REM Requires Julia 1.10+ on PATH. First launch installs QML.jl + Qt6_jll
REM (~300 MB download), so be patient — subsequent launches are fast.

cd /d "%~dp0\.."

where julia >nul 2>nul
if not %ERRORLEVEL%==0 (
    echo.
    echo Julia not found on PATH. Install from https://julialang.org/downloads/
    echo and re-run this file.
    echo.
    pause
    goto :eof
)

julia --threads=auto --project=qml_gui qml_gui/cshore_qml.jl
