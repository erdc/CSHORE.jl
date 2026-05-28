#!/usr/bin/env bash
# Launcher for the CSHORE.jl web GUI (HTTP.jl + Plotly stack). Mirrors
# qml_gui/run_qml.sh.
#
# 1. Kills any pre-existing julia processes running web_gui/app.jl — they
#    otherwise hold port 8000 and the new process fails to bind.
# 2. Runs Julia with --threads=auto so simulations don't block the HTTP
#    event loop.
set -euo pipefail
cd "$(dirname "$0")/.."

PORT="${PORT:-8000}"

# Identify any stale instances and terminate them. pkill is best-effort —
# we still continue if it finds nothing.
if pgrep -f "julia.*web_gui/app.jl" >/dev/null 2>&1; then
    echo "[run_web] found stale web_gui Julia process(es); terminating…"
    pkill -9 -f "julia.*web_gui/app.jl" || true
    sleep 1
fi

# If port is still occupied (something else holding it), warn loudly.
if lsof -nP -iTCP:"$PORT" -sTCP:LISTEN >/dev/null 2>&1; then
    echo "[run_web] WARNING: port $PORT still in use after pkill; the new"
    echo "[run_web]          process will fail to bind. Kill the holder:"
    lsof -nP -iTCP:"$PORT" -sTCP:LISTEN
fi

exec julia --threads=auto --project=web_gui web_gui/app.jl
