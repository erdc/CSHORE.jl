#!/usr/bin/env bash
# run_qml.command -- Double-click launcher for the CSHORE.jl QML GUI on macOS.
# Requires Julia 1.10+ on PATH. First launch installs QML.jl + Qt6_jll
# (~300 MB download), so be patient -- subsequent launches are fast.
set -e
cd "$(dirname "$0")/.."

if ! command -v julia >/dev/null 2>&1; then
    osascript -e 'display dialog "Julia not found.\n\nInstall Julia 1.10 LTS from https://julialang.org/downloads/ (or via juliaup) and re-run this file." buttons {"OK"} default button 1 with icon stop'
    exit 1
fi

exec julia --threads=auto --project=qml_gui qml_gui/cshore_qml.jl
