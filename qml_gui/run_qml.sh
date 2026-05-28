#!/usr/bin/env bash
# run_qml.sh -- Launcher for the CSHORE.jl QML GUI on Linux.
# Requires Julia 1.10+ on PATH. First launch installs QML.jl + Qt6_jll
# (~300 MB download); subsequent launches are fast.
#
# On Linux you also need the system XCB libraries Qt links against:
#   sudo apt install libxcb-cursor0 libxcb-icccm4 libxcb-keysyms1   (Debian/Ubuntu)
#   sudo dnf install xcb-util-cursor xcb-util-keysyms                (Fedora)
set -e
cd "$(dirname "$0")/.."

if ! command -v julia >/dev/null 2>&1; then
    echo "Julia not found on PATH." >&2
    echo "Install Julia 1.10 LTS from https://julialang.org/downloads/ or via juliaup." >&2
    exit 1
fi

exec julia --threads=auto --project=qml_gui qml_gui/cshore_qml.jl
