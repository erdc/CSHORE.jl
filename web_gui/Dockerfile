# =============================================================================
# CSHORE.jl Web GUI — Fly.io deploy image.
#
# Lean stack: HTTP.jl + vanilla HTML/JS + Plotly.js (client-side). No
# CairoMakie / Makie / Genie / Stipple — those used to live here and got
# dropped because together they peaked at ~3-4 GB during precompile and
# blew up Fly's 2 GB remote builder.
#
# Build strategy: install packages but DON'T precompile at build time
# (JULIA_PKG_PRECOMPILE_AUTO=0). Precompile + JIT specialization happen
# at server boot in app.jl, where a synthetic _warmup() call pre-pays the
# specialization cost so the user's first /run is fast.
#
# Boot time: ~5-10 min after deploy (package precompile + warmup). The
# port doesn't open until that finishes. Subsequent /run requests are
# ~10-60 s (just the simulation itself).
# =============================================================================

FROM julia:1.12-bookworm

# Runtime libs only: TLS for any HTTP fetches + libgomp for the HDF5 /
# NetCDF JLLs that NCDatasets uses to read the simulation output. No
# Cairo / Pango / fontconfig because all plotting happens in the browser.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl \
        libstdc++6 libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . /app

RUN test -f /app/Project.toml \
 && test -f /app/web_gui/Project.toml \
 && test -f /app/web_gui/app.jl

ENV JULIA_DEPOT_PATH=/app/.julia
ENV JULIA_PROJECT=/app/web_gui
ENV HOST=0.0.0.0

# Install only — precompile is deferred to first server boot to keep
# build memory under Fly's 2 GB cap. Pkg.develop wires the parent
# CSHORE package at /app into the web_gui project.
RUN julia --project=/app/web_gui -e ' \
        ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0"; \
        using Pkg; \
        Pkg.develop(path="/app"); \
        Pkg.instantiate(); \
        @info "instantiated WITHOUT precompile — first boot will compile"; \
    '

EXPOSE 8000

CMD ["julia", "--threads=auto", "--project=/app/web_gui", "/app/web_gui/app.jl"]
