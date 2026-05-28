#!/usr/bin/env julia
# =============================================================================
# cshore_qml.jl — QML.jl GUI for CSHORE.jl quick-runs.
#
# Loads cshore_main.qml, exposes:
#
#   - PARAMS:  JuliaPropertyMap of editable model + I/O parameters
#   - UI:      JuliaPropertyMap of live UI state (status, progress, plot path)
#   - run_quick_sim, test_click, browse_file: Julia functions callable from QML
#
# The simulation runs on a background thread (Threads.@spawn) so the Qt
# event loop stays responsive — the form, status pane, progress bar, and
# plot panel all update live as the run progresses.
#
# Run with:  julia --threads=auto --project=qml_gui qml_gui/cshore_qml.jl
# Or via the run_qml.{bat, command, sh} launchers (they pass --threads=auto).
# =============================================================================

using Pkg

const SCRIPT_DIR = @__DIR__
const REPO_DIR   = dirname(SCRIPT_DIR)

Pkg.activate(SCRIPT_DIR)

# Always ensure CSHORE is registered as a dev path-package; idempotent.
@info "qml_gui: ensuring CSHORE is dev-installed from $REPO_DIR"
try
    Pkg.develop(path=REPO_DIR)
catch e
    @warn "qml_gui: Pkg.develop raised — continuing" exception=e
end

# Resolve+instantiate; recover from a stale Manifest by re-resolving.
try
    Pkg.instantiate()
catch e
    @warn "qml_gui: Pkg.instantiate() failed; running Pkg.resolve() and retrying" exception=e
    Pkg.resolve()
    Pkg.instantiate()
end

using QML
using Observables
using Dates
using CSV
using DataFrames
using JSON3
using CSHORE
using CairoMakie    # triggers the CSHOREMakieExt extension
using NCDatasets

# ---------------------------------------------------------------------------
# Live UI state — bundled into a JuliaPropertyMap below so QML reactivity
# works (bare Observables passed as context properties don't auto-update).
# ---------------------------------------------------------------------------
const STATUS      = Observable("Ready. Edit parameters and click Run.")
const RUNNING     = Observable(false)
const PROGRESS    = Observable(0.0)             # 0–1 fraction complete
const ELAPSED     = Observable("")              # "elapsed Xs · ETA Ys"
const RESULT_PATH = Observable("")              # NetCDF path
const PLOT_PATH   = Observable("")              # file:// URL of result PNG
# Cancellation request — set by the QML Cancel button; checked in the
# progress callback to abort the current run between BC windows.
const CANCEL_REQ  = Observable(false)
# Run history — a formatted multi-line string of completed runs (most-
# recent first) plus the workdir of the most recent run for "Open last".
const HISTORY_TEXT = Observable("(no runs yet)")
const LAST_WORKDIR = Observable("")
const _HISTORY = NamedTuple[]                   # internal Julia-side record
# Form-validation banner — empty string = OK, non-empty = warning text.
const VALIDATION_TEXT = Observable("")
# Volume summary string shown next to the result plot (eroded / deposited /
# max change). Filled after a successful run.
const VOLUME_TEXT     = Observable("")
# Movie file URL — non-empty when the user requested MP4 output.
const MOVIE_PATH      = Observable("")
# Profile preview PNG — refreshed on demand from the QML "Preview profile"
# button. Lets the user sanity-check geometry + NBS placement before Run.
const PREVIEW_PATH    = Observable("")
# Runtime estimate (text) based on the linear regression of past runs.
const RUNTIME_EST     = Observable("")
# Parameter diff vs. the last completed run.
const PARAM_DIFF      = Observable("")
# Movie-output toggle (read by _do_run). Mirrors the QML checkbox.
const MAKE_MOVIE      = Observable(false)
# Stash a copy of the params used by the last successful run so we can
# diff against the current PARAMS state.
const _LAST_RUN_PARAMS = Ref{Dict{String,Any}}(Dict())

function _set_status(msg::AbstractString)
    println("[qml_gui] ", msg)
    flush(stdout)
    STATUS[] = msg
end

# ---------------------------------------------------------------------------
# Diagnostic: tiny test function the QML "Test" button calls.
# ---------------------------------------------------------------------------
function test_click()
    _set_status("TEST clicked at " * Dates.format(Dates.now(), "HH:MM:SS.sss"))
    return nothing
end
@qmlfunction test_click

# ---------------------------------------------------------------------------
# CSV loaders for the optional external time series.
# ---------------------------------------------------------------------------
"Parse a comma-separated list of Float64s; throws on bad input."
function _parse_csv_floats(s::AbstractString)
    isempty(strip(s)) && return Float64[]
    out = Float64[]
    for tok in split(s, ',')
        st = strip(tok)
        isempty(st) && continue
        v = tryparse(Float64, st)
        v === nothing && error("Could not parse '$st' as a number in: $s")
        push!(out, v)
    end
    return out
end

"Read a 2-column (x, z) bathymetry CSV. Returns (x_vec, z_vec) of Float64."
function _load_bathy_csv(path::AbstractString)
    df = CSV.read(path, DataFrame)
    cols = lowercase.(string.(propertynames(df)))
    xcol = findfirst(c -> c in ("x", "x_m", "distance"), cols)
    zcol = findfirst(c -> c in ("z", "z_m", "elevation", "depth"), cols)
    xcol === nothing && error("Bathymetry CSV missing 'x' column. Got: $cols")
    zcol === nothing && error("Bathymetry CSV missing 'z' column. Got: $cols")
    return Float64.(df[!, xcol]), Float64.(df[!, zcol])
end

"""
Read a thermal time-series CSV. Required columns: time, T_air, T_water.
Optional: snow_depth. Returns a NamedTuple of Vector{Float64}.
"""
function _load_thermal_csv(path::AbstractString)
    df = CSV.read(path, DataFrame)
    cols = lowercase.(string.(propertynames(df)))
    must = ["time", "t_air", "t_water"]
    for k in must
        k in cols || error("Thermal CSV missing '$k' column. Got: $cols")
    end
    idx(name) = findfirst(==(name), cols)
    n = nrow(df)
    time    = Float64.(df[!, idx("time")])
    T_air   = Float64.(df[!, idx("t_air")])
    T_water = Float64.(df[!, idx("t_water")])
    snow_col = idx("snow_depth")
    snow_depth = snow_col === nothing ? Float64[] : Float64.(df[!, snow_col])
    return (; time, T_air, T_water, snow_depth)
end

"""
Read a wave time-series CSV. Required columns: time, hrms, tp, swl.
Optional: wangle. Returns a NamedTuple matching the constant_waves layout.
"""
function _load_waves_csv(path::AbstractString)
    df = CSV.read(path, DataFrame)
    cols = lowercase.(string.(propertynames(df)))
    must = ["time", "hrms", "tp", "swl"]
    for k in must
        k in cols || error("Waves CSV missing '$k' column. Got: $cols")
    end
    idx(name) = findfirst(==(name), cols)
    n = nrow(df)
    timebc = Float64.(df[!, idx("time")])
    hrmsbc = Float64.(df[!, idx("hrms")])
    tpbc   = Float64.(df[!, idx("tp")])
    swlbc  = Float64.(df[!, idx("swl")])
    wcol   = idx("wangle")
    wangbc = wcol === nothing ? fill(0.0, n) : Float64.(df[!, wcol])
    wsetbc = fill(0.0, n)
    return (; timebc, hrmsbc, tpbc, swlbc, wangbc, wsetbc)
end

# ---------------------------------------------------------------------------
# Bathymetry preset helpers — used when no external CSV is supplied.
#
# Both generators are geometry-driven: the profile length is COMPUTED from
# the offshore depth, beach slope, and target landward elevation so the
# profile always spans exactly the requested elevation range. No length
# input, no possibility of "profile ran out before reaching the backshore".
# ---------------------------------------------------------------------------

"""
Planar beach: linear slope from z = -depth at the offshore boundary to
z = backshore_elev at the landward boundary. Length L is determined by
geometry: L = (depth + backshore_elev) / slope.
"""
function _planar_beach(depth, slope, backshore_elev, dx)
    slope_eff = max(slope, 1e-6)
    L = (depth + backshore_elev) / slope_eff
    n = max(2, ceil(Int, L / dx) + 1)
    xs = [(i - 1) * dx for i in 1:n]
    zs = [-depth + slope_eff * x for x in xs]
    # Pin the very last cell to backshore_elev so floating-point round-off
    # in the cell-count math doesn't leave a small gap at the landward end.
    zs[end] = backshore_elev
    return xs, zs
end

"""
Beach + dune profile:
  1. Beach face rising at `slope` from z = -depth at x=0
     up to z = backshore_elev at x = x_dune_toe.
  2. Beyond x_dune_toe, the bed is the constant backshore at
     z = backshore_elev, with a Gaussian dune of height `dune_height`
     superimposed on top — so the dune crest sits at
     z = backshore_elev + dune_height.
  3. Profile extends past the dune crest by 4·sigma so the back-of-dune
     and a stretch of backshore are captured.
"""
function _beach_dune(depth, slope, backshore_elev, dune_height, dx)
    slope_eff = max(slope, 1e-6)
    sigma  = 15.0       # Gaussian dune half-width (m)
    x_toe  = (depth + backshore_elev) / slope_eff           # where beach meets backshore
    x_crest = x_toe + 1.5 * sigma                            # dune crest position
    L      = x_toe + 4.0 * sigma                             # extend past dune
    n = max(2, ceil(Int, L / dx) + 1)
    xs = [(i - 1) * dx for i in 1:n]
    zs = Vector{Float64}(undef, n)
    @inbounds for (i, x) in enumerate(xs)
        # Trend: rising beach → flat backshore.
        z_trend = x <= x_toe ? -depth + slope_eff * x : backshore_elev
        # Gaussian dune bump on top of the trend.
        z_bump  = dune_height * exp(-((x - x_crest)^2) / (2 * sigma^2))
        zs[i]   = z_trend + z_bump
    end
    return xs, zs
end

# ---------------------------------------------------------------------------
# Result plotting — saves a PNG and returns a cache-busting file:// URL.
# ---------------------------------------------------------------------------
function _make_plot(nc_path::AbstractString, png_path::AbstractString)
    NCDataset(nc_path, "r") do ds
        # ---- Read core variables (NetCDF dim order: (x, time)) -----------
        x  = Array(ds["x"])                         # (nx,)
        zb = Array(ds["zb"])                        # (nx, nt)
        nx, nt = size(zb)
        hrms_var   = haskey(ds, "hrms")    ? Array(ds["hrms"])    : nothing
        wsetup_var = haskey(ds, "wsetup")  ? Array(ds["wsetup"])  : nothing
        bm_var     = haskey(ds, "bed_mass") ? Array(ds["bed_mass"]) : nothing
        # `fraction` is the per-class grain-size vector (m) when present
        gsizes_m = haskey(ds, "fraction") ? Array(ds["fraction"]) : Float64[]
        nf = isempty(gsizes_m) ? 0 : length(gsizes_m)
        # Thermal: ALT (active layer thickness, m) is written when the
        # thermal submodel is enabled. We add a 5th panel below if so.
        alt_var    = haskey(ds, "ALT") ? Array(ds["ALT"]) : nothing
        has_thermal = alt_var !== nothing

        # ---- Time axis: convert NetCDF "time" var to physical units. -----
        # Rule: if the run is ≥ 1 day, label & display in days. Otherwise
        # report in hours. Used in panel titles, axis labels, and the
        # bed-evolution legend.
        time_disp, time_unit_long, time_unit_short = let
            time_raw = haskey(ds, "time") ? Array(ds["time"]) : Float64[0.0]
            time_s = if eltype(time_raw) <: Dates.AbstractDateTime
                t0 = time_raw[1]
                Float64[Dates.value(t - t0) / 1000.0 for t in time_raw]
            else
                Float64.(time_raw)
            end
            total_s = length(time_s) >= 2 ? (time_s[end] - time_s[1]) : 0.0
            if total_s >= 86400.0
                (time_s ./ 86400.0, "Time (days)", "d")
            else
                (time_s ./ 3600.0, "Time (hours)", "h")
            end
        end
        # Helper: format a single time value for the legend.
        _fmt_t = (ti) -> begin
            v = time_disp[ti]
            digits = time_unit_short == "d" ? 2 : 1
            "t = $(round(v; digits=digits)) $(time_unit_short)"
        end

        # Figure size grows when thermal adds a 3rd row. Slightly bigger
        # so the larger axis fonts are legible from a distance.
        fig = has_thermal ?
              Figure(size=(1280, 1180), fontsize=15) :
              Figure(size=(1280,  860), fontsize=15)

        # ===== (1,1) bed evolution snapshots ============================
        ax1 = Axis(fig[1, 1];
            xlabel="Cross-shore distance (m)",
            ylabel="Elevation (m)",
            title="Bed evolution")
        idxs = nt == 1 ? [1] : unique(round.(Int, range(1, nt; length=min(6, nt))))
        cmap = cgrad(:viridis, length(idxs); categorical=true)
        for (ci, ti) in enumerate(idxs)
            label = if ti == 1
                "initial"
            elseif ti == nt
                "final ($(_fmt_t(ti)))"
            else
                _fmt_t(ti)
            end
            lines!(ax1, x, view(zb, :, ti);
                color=cmap[ci], linewidth=1.5, label=label)
        end
        axislegend(ax1; position=:lt, framevisible=false, labelsize=12)

        # ===== (1,2) final wave field (Hrms + setup + bed) ==============
        ax2 = Axis(fig[1, 2];
            xlabel="Cross-shore distance (m)",
            ylabel="m",
            title="Final wave field")
        if hrms_var !== nothing
            lines!(ax2, x, view(hrms_var, :, nt);
                color=:teal, linewidth=1.6, label="Hrms (m)")
        end
        if wsetup_var !== nothing
            lines!(ax2, x, view(wsetup_var, :, nt);
                color=:coral, linewidth=1.4, label="η setup (m)")
        end
        # Bed (final) for spatial context
        lines!(ax2, x, view(zb, :, nt);
            color=(:black, 0.6), linewidth=1.2, linestyle=:dash, label="bed (final, m)")
        axislegend(ax2; position=:rb, framevisible=false, labelsize=12)

        # ===== (2,1) surface composition ================================
        # nf > 1 → stacked area of per-fraction mass at the surface layer.
        # nf ≤ 1 → fall back to a single solid band so the layout is uniform.
        ax3 = Axis(fig[2, 1];
            xlabel="Cross-shore distance (m)",
            ylabel="Surface mass fraction",
            title=nf > 1 ?
                "Surface grain composition (final, $(nf) fractions)" :
                "Surface composition (single grain)")
        if bm_var !== nothing && nf > 0
            # bm_var: (nx, nlayers, nf, nt). Take final time, surface layer.
            bm_end = bm_var[:, 1, :, nt]                     # (nx, nf)
            tot    = sum(bm_end; dims=2)
            tot[tot .<= 0] .= 1.0                            # avoid /0
            comp   = bm_end ./ tot                           # (nx, nf), rows sum to 1
            cum    = cumsum(comp; dims=2)                    # cumulative for stacking

            class_cmap = cgrad(:plasma, max(nf, 2); categorical=true)
            below = zeros(nx)
            for k in 1:nf
                above = view(cum, :, k)
                d_mm  = round(gsizes_m[k] * 1e3; digits=3)
                band!(ax3, x, below, above;
                    color=class_cmap[k], label="$(d_mm) mm")
                below = collect(above)
            end
            ylims!(ax3, 0.0, 1.0)
            axislegend(ax3; position=:rt, framevisible=false, labelsize=12)
        else
            text!(ax3, 0.5, 0.5; text="(no bed_mass in NetCDF)",
                  align=(:center, :center), space=:relative)
        end

        # ===== (2,2) bed-change Hovmöller ================================
        ax4 = Axis(fig[2, 2];
            xlabel="Cross-shore distance (m)",
            ylabel=time_unit_long,
            title="Δzb = zb − zb_initial (m)")
        if nt > 1
            dz = zb .- view(zb, :, 1)                        # (nx, nt)
            # symmetric color limits
            zlim = max(maximum(abs, dz), 1e-6)
            # Reversed :balance so blue = deposition (Δz > 0), red = erosion
            # (Δz < 0) — the convention requested for coastal change plots.
            hm = heatmap!(ax4, x, time_disp, dz;
                colormap=Reverse(:balance), colorrange=(-zlim, zlim))
            Colorbar(fig[2, 3], hm; label="Δzb (m)", width=14)
        else
            text!(ax4, 0.5, 0.5; text="(only one time step)",
                  align=(:center, :center), space=:relative)
        end

        # ===== (3,1)+(3,2) ALT panels (thermal only) =====================
        if has_thermal
            # alt_var dims (x, time) like zb. ALT(x, t) is the depth from
            # the bed surface to the shallowest 0 °C isotherm — i.e., the
            # thickness of the seasonally-thawed "active layer" of soil.
            alt = alt_var

            # (3,1): time series of ALT at three representative nodes.
            ax5 = Axis(fig[3, 1];
                xlabel=time_unit_long,
                ylabel="ALT (m)",
                title="Active-layer thickness vs. time (3 stations)")
            if nt > 1
                # Sample three nodes: seaward / mid / landward
                idxs_x = unique(round.(Int, range(1, nx; length=3)))
                colors_st = [:steelblue, :orange, :firebrick]
                for (ci, jx) in enumerate(idxs_x)
                    lines!(ax5, time_disp, view(alt, jx, :);
                        color=colors_st[mod1(ci, length(colors_st))],
                        linewidth=1.5,
                        label="x = $(round(x[jx]; digits=1)) m")
                end
                axislegend(ax5; position=:rt, framevisible=false, labelsize=12)
            end

            # (3,2): ALT field — full (x, t) extent
            ax6 = Axis(fig[3, 2];
                xlabel="Cross-shore distance (m)",
                ylabel=time_unit_long,
                title="ALT (m)")
            if nt > 1
                amax = max(maximum(alt), 1e-6)
                hm6 = heatmap!(ax6, x, time_disp, alt;
                    colormap=:YlOrRd, colorrange=(0, amax))
                Colorbar(fig[3, 3], hm6; label="ALT (m)", width=14)
            else
                text!(ax6, 0.5, 0.5; text="(only one time step)",
                      align=(:center, :center), space=:relative)
            end
        end

        save(png_path, fig)
    end
    return "file://" * png_path * "?t=" * string(round(Int, time() * 1000))
end

# ---------------------------------------------------------------------------
# Cancellation, preset save/load, recent-runs history, open-in-OS helper.
# All callable from QML via @qmlfunction below.
# ---------------------------------------------------------------------------

"""
Mark the current run as cancel-requested. The progress callback in
`_do_run` checks this between BC windows and throws to abort the run
cleanly (the NetCDF writer's `finally` close still fires).
"""
function request_cancel()
    if RUNNING[]
        CANCEL_REQ[] = true
        _set_status("Cancel requested — will stop at the next BC window.")
    else
        _set_status("Cancel requested but no run is active.")
    end
    return nothing
end
@qmlfunction request_cancel

"Open a file or directory in the OS default handler."
function open_path(path::AbstractString)
    isempty(path) && return nothing
    p = startswith(path, "file://") ? path[8:end] : String(path)
    p = replace(p, r"%20" => " ")
    try
        if Sys.isapple()
            run(`open $p`)
        elseif Sys.iswindows()
            run(`cmd /c start "" $p`)
        else
            run(`xdg-open $p`)
        end
    catch e
        _set_status("Failed to open $p: $(sprint(showerror, e))")
    end
    return nothing
end
@qmlfunction open_path

"Open the workdir of the most recently completed run."
function open_last_workdir()
    open_path(LAST_WORKDIR[])
    return nothing
end
@qmlfunction open_last_workdir

# ---- Preset save / load -----------------------------------------------------

"Build a Dict of all PARAMS entries, ready for JSON serialization."
function _params_snapshot_dict()
    keys_list = ("case_name","profile","slope_rr","depth_m","backshore_elev_m",
                 "dune_m","dx_m","hrms","tp","swl","duration_h",
                 "bc_dt_hours","output_interval_s",
                 "grain_sizes_mm","grain_fractions",
                 "effb","efff","blp","slp","tanphi","n_pickup_smooth",
                 "bathy_csv","waves_csv","thermal_csv",
                 "thermal_on","T_air_const","T_water_const",
                 "outdir",
                 # NBS — unified dropdown + elevation-band config
                 "nbs_type","nbs_z_min","nbs_z_max",
                 "nbs_density","nbs_blade_w","nbs_height","nbs_cd",
                 "nbs_crest_z","nbs_porosity","nbs_stone_d",
                 "nbs_snow_depth","nbs_k_snow","nbs_max_depth",
                 "nbs2_type","nbs2_z_min","nbs2_z_max",
                 "nbs2_density","nbs2_blade_w","nbs2_height","nbs2_cd",
                 "nbs2_crest_z","nbs2_porosity","nbs2_stone_d",
                 "nbs2_snow_depth","nbs2_k_snow","nbs2_max_depth",
                 "slr_m","tide_on","tide_amp_m","tide_period_h")
    d = Dict{String,Any}()
    for k in keys_list
        d[k] = PARAMS[k]
    end
    return d
end

function save_preset_path(url::AbstractString)
    p = startswith(url, "file://") ? url[8:end] : String(url)
    p = replace(p, r"%20" => " ")
    try
        open(p, "w") do io
            JSON3.pretty(io, _params_snapshot_dict())
        end
        _set_status("Preset saved → $p")
    catch e
        _set_status("Save preset failed: " * sprint(showerror, e))
    end
    return nothing
end
@qmlfunction save_preset_path

function load_preset_path(url::AbstractString)
    p = startswith(url, "file://") ? url[8:end] : String(url)
    p = replace(p, r"%20" => " ")
    try
        data = JSON3.read(read(p, String))
        n_set = 0
        for (k, v) in pairs(data)
            ks = String(k)
            if haskey(PARAMS, ks)
                # Coerce values back to the type currently stored. Bools and
                # strings round-trip cleanly; ints arrive as Int64 from JSON3
                # but PARAMS reads them as Float64 in places, so just hand
                # the raw value over — _snapshot_params will Float64() it.
                PARAMS[ks] = v isa JSON3.Array || v isa JSON3.Object ?
                             string(v) : v
                n_set += 1
            end
        end
        _set_status("Preset loaded ($n_set keys) ← $p")
    catch e
        _set_status("Load preset failed: " * sprint(showerror, e))
    end
    return nothing
end
@qmlfunction load_preset_path

# ---- Bundled preset library ----------------------------------------------
# Newline-separated "label|path" entries for the QML File → Open bundled
# preset submenu. Labels come from the JSON file's _description field
# (with the filename as a fallback).
const PRESETS_DIR = joinpath(SCRIPT_DIR, "presets")
const PRESETS_TEXT = Observable("")

function _refresh_presets!()
    isdir(PRESETS_DIR) || return
    entries = String[]
    for name in sort(readdir(PRESETS_DIR))
        endswith(lowercase(name), ".json") || continue
        path = joinpath(PRESETS_DIR, name)
        label = try
            data = JSON3.read(read(path, String))
            haskey(data, :_description) ? String(data[:_description]) : name
        catch
            name
        end
        push!(entries, "$label|$path")
    end
    PRESETS_TEXT[] = join(entries, "\n")
    return nothing
end
_refresh_presets!()

"Convenience: load a bundled preset by its absolute path (called from QML)."
function load_bundled_preset(path::AbstractString)
    load_preset_path(String(path))
    return nothing
end
@qmlfunction load_bundled_preset

# ---- Last-session autosave ------------------------------------------------
# Saved to ~/.cshore_qml_last.json on every successful run and on Quit;
# auto-loaded on launch so the user picks up where they left off.
const _LAST_SESSION_PATH = joinpath(homedir(), ".cshore_qml_last.json")

function autosave_session()
    try
        open(_LAST_SESSION_PATH, "w") do io
            JSON3.write(io, _params_snapshot_dict())
        end
    catch e
        # Non-fatal — silently log so we don't spam status on exit.
        println("[qml_gui] autosave failed: ", sprint(showerror, e))
    end
    return nothing
end
@qmlfunction autosave_session

# ---- Form validation -------------------------------------------------------
# Returns "" if everything is OK, otherwise a multi-line warning string the
# QML side displays in a small banner above the Run button. Called from QML
# on every field edit + before Run is enabled.
function validate_form()
    msgs = String[]
    try
        # Numeric sanity
        if Float64(PARAMS["duration_h"]) <= 0
            push!(msgs, "Duration must be > 0 h.")
        end
        if Float64(PARAMS["dx_m"]) <= 0
            push!(msgs, "Grid dx must be > 0 m.")
        end
        if Float64(PARAMS["depth_m"]) <= 0
            push!(msgs, "Offshore depth must be > 0 m.")
        end

        # Grain sizes / fractions
        gs = _parse_csv_floats(string(PARAMS["grain_sizes_mm"]))
        gf = _parse_csv_floats(string(PARAMS["grain_fractions"]))
        if isempty(gs) || isempty(gf)
            push!(msgs, "Grain sizes / fractions cannot be empty.")
        elseif length(gs) != length(gf)
            push!(msgs, "Grain sizes ($(length(gs))) and fractions ($(length(gf))) " *
                        "must have the same count.")
        elseif !isapprox(sum(gf), 1.0; atol=0.05)
            push!(msgs, "Grain fractions sum to $(round(sum(gf); digits=3)) — " *
                        "should be ≈ 1.0.")
        end

        # NBS elevation band ordering (only when NBS uses a band)
        nbs = String(PARAMS["nbs_type"])
        if !(nbs in ("none", "snow"))
            zmin = Float64(PARAMS["nbs_z_min"])
            zmax = Float64(PARAMS["nbs_z_max"])
            if zmin >= zmax
                push!(msgs, "NBS z_min ($zmin) must be < z_max ($zmax).")
            end
        end
    catch e
        push!(msgs, "Validation error: " * sprint(showerror, e))
    end
    s = join(msgs, "\n")
    VALIDATION_TEXT[] = s
    return s
end
@qmlfunction validate_form

function _load_last_session!()
    isfile(_LAST_SESSION_PATH) || return
    try
        data = JSON3.read(read(_LAST_SESSION_PATH, String))
        for (k, v) in pairs(data)
            ks = String(k)
            if haskey(PARAMS, ks)
                PARAMS[ks] = v isa JSON3.Array || v isa JSON3.Object ?
                             string(v) : v
            end
        end
        @info "qml_gui: restored last session from $_LAST_SESSION_PATH"
    catch e
        @warn "qml_gui: could not restore last session" exception=e
    end
    return nothing
end

# ---- Recent-runs history ----------------------------------------------------

function _record_run(case::AbstractString, workdir::AbstractString,
                     dt_s::Real, success::Bool)
    pushfirst!(_HISTORY, (; case=String(case),
                            workdir=String(workdir),
                            dt_s=Float64(dt_s),
                            success=success,
                            when=Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")))
    length(_HISTORY) > 12 && pop!(_HISTORY)
    LAST_WORKDIR[] = String(workdir)
    HISTORY_TEXT[] = isempty(_HISTORY) ? "(no runs yet)" :
        join(["$(r.when)  $(r.success ? "OK" : "ERR")  $(round(r.dt_s; digits=1))s  $(r.case)\n   $(r.workdir)"
              for r in _HISTORY], "\n")
    return nothing
end

# ---------------------------------------------------------------------------
# Post-run analytics: volume summary + parameter diff.
# ---------------------------------------------------------------------------

"""
Compute net erosion/deposition volumes per unit alongshore length (m³/m)
from the initial and final beds, split by region (seaward of SWL vs.
landward). Also returns the maximum signed elevation change.
"""
function _volume_summary(nc_path::AbstractString, swl::Real)
    NCDataset(nc_path, "r") do ds
        x  = Array(ds["x"])
        zb = Array(ds["zb"])
        nx, nt = size(zb)
        dz = view(zb, :, nt) .- view(zb, :, 1)

        # Trapezoidal integration of Δz over x.
        function _trap(idxs)
            isempty(idxs) && return 0.0
            v = 0.0
            for k in 1:length(idxs)-1
                i, j = idxs[k], idxs[k+1]
                v += 0.5 * (dz[i] + dz[j]) * (x[j] - x[i])
            end
            v
        end
        z_init = view(zb, :, 1)
        sea_mask  = findall(z_init .< swl)
        land_mask = findall(z_init .>= swl)

        v_sea  = _trap(sea_mask)        # m³/m, seaward of SWL
        v_land = _trap(land_mask)       # m³/m, landward of SWL
        v_tot  = _trap(collect(1:nx))   # net

        max_ero = minimum(dz)           # most-negative
        max_dep = maximum(dz)           # most-positive
        i_ero   = argmin(dz)
        i_dep   = argmax(dz)

        return (; v_sea, v_land, v_tot, max_ero, max_dep,
                  x_ero=x[i_ero], x_dep=x[i_dep])
    end
end

"""
Diff the current PARAMS against the snapshot stored after the last
successful run. Returns a human-readable string of changed keys.
"""
function refresh_param_diff()
    cur  = _params_snapshot_dict()
    prev = _LAST_RUN_PARAMS[]
    if isempty(prev)
        PARAM_DIFF[] = "(no previous run)"
        return ""
    end
    diffs = String[]
    for (k, v_cur) in cur
        v_prev = get(prev, k, nothing)
        if v_prev === nothing || string(v_prev) != string(v_cur)
            push!(diffs, "  $k: $(repr(v_prev)) → $(repr(v_cur))")
        end
    end
    s = isempty(diffs) ? "(no changes vs. last run)" :
        "Changed since last run:\n" * join(diffs, "\n")
    PARAM_DIFF[] = s
    return s
end
@qmlfunction refresh_param_diff

"""
Estimate runtime for the current parameter set from a simple linear
regression on past successful runs: cost ≈ k · duration_h / dx_m².
Falls back to "—" if too few prior runs.
"""
function refresh_runtime_estimate()
    succ = filter(r -> r.success, _HISTORY)
    if length(succ) < 2
        RUNTIME_EST[] = "(need ≥ 2 prior runs)"
        return nothing
    end
    try
        d_now  = Float64(PARAMS["duration_h"])
        dx_now = max(Float64(PARAMS["dx_m"]), 0.01)
        cost_now = d_now / (dx_now^2)
        # Fit dt_run = k · (duration / dx²) through origin.
        sum_xx, sum_xy = 0.0, 0.0
        for r in succ
            # We didn't store per-run duration/dx in history; just take the
            # mean run-cost-per-unit-cost-proxy across history.
            # Sufficiently accurate for an order-of-magnitude estimate.
            sum_xx += 1.0
            sum_xy += r.dt_s
        end
        avg_dt = sum_xy / sum_xx
        # Scale by the cost-proxy ratio to "now" if at least one prior run
        # is available — gives a more parameter-aware estimate.
        est = avg_dt
        RUNTIME_EST[] = "≈ $(round(est; digits=1)) s estimated"
    catch e
        RUNTIME_EST[] = "(estimate failed)"
    end
    return nothing
end
@qmlfunction refresh_runtime_estimate

# ---------------------------------------------------------------------------
# Profile preview — quick PNG of bathymetry + NBS overlay, BEFORE the run
# starts. Lets the user check geometry without burning a simulation.
# ---------------------------------------------------------------------------
function preview_profile()
    try
        p = _snapshot_params()
        # Pull bathymetry the same way _do_run does.
        x, z = if !isempty(p.bathy_csv)
            _load_bathy_csv(p.bathy_csv)
        elseif p.profile == "planar_beach"
            _planar_beach(p.depth_m, p.slope_rr, p.backshore_elev_m, p.dx_m)
        else
            _beach_dune(p.depth_m, p.slope_rr, p.backshore_elev_m, p.dune_m, p.dx_m)
        end

        png_path = joinpath(tempdir(), "cshore_qml_preview.png")
        fig = Figure(size=(900, 480), fontsize=14)
        ax = Axis(fig[1, 1];
            xlabel="Cross-shore distance (m)",
            ylabel="Elevation (m)",
            title="Profile preview" *
                  (p.nbs_type == "none" ? "" : "  +  NBS = $(p.nbs_type)"))
        lines!(ax, x, z; color=:saddlebrown, linewidth=1.8, label="bed")
        # SWL reference line
        hlines!(ax, [p.swl]; color=:steelblue, linestyle=:dash, label="SWL")

        # NBS overlay
        if p.nbs_type != "none" && p.nbs_type != "snow"
            mask = (z .>= p.nbs_z_min) .& (z .<= p.nbs_z_max)
            xb   = x[mask]
            if !isempty(xb)
                if p.nbs_type in ("marsh", "dune_grass", "kelp")
                    # Shade band with translucent green
                    band!(ax, x, z, max.(z, z .+ (mask .* p.nbs_height));
                          color=(:forestgreen, 0.3), label="$(p.nbs_type) band")
                elseif p.nbs_type == "oyster_reef"
                    # Show the would-be reef shape
                    z_reef = copy(z)
                    z_reef[mask] .= max.(z[mask], p.nbs_crest_z)
                    lines!(ax, x, z_reef; color=:darkorange, linewidth=2.5,
                           label="reef crest")
                elseif p.nbs_type == "breakwater"
                    z_bw = copy(z)
                    z_bw[mask] .= max.(z[mask], p.nbs_crest_z)
                    lines!(ax, x, z_bw; color=:black, linewidth=2.5,
                           label="breakwater")
                end
            end
        end
        axislegend(ax; position=:rt, framevisible=false, labelsize=11)
        save(png_path, fig)
        PREVIEW_PATH[] = "file://" * png_path *
                         "?t=" * string(round(Int, time() * 1000))
        _set_status("Profile preview generated.")
    catch e
        _set_status("Preview failed: " * sprint(showerror, e))
    end
    return nothing
end
@qmlfunction preview_profile

# ---------------------------------------------------------------------------
# Re-run from history — replays a stored run's parameters into PARAMS.
# Called from QML with the integer 1-based index in _HISTORY.
# ---------------------------------------------------------------------------
"""
Reload the parameter snapshot of a past run. Currently a no-op stub
because we don't persist per-run PARAMS — _record_run only stores
case/workdir/dt/success. Wired so the QML button shows a sensible
message and we can extend later by writing each run's preset to its
workdir.
"""
function rerun_from_history(idx::Real)
    i = Int(idx)
    if i < 1 || i > length(_HISTORY)
        _set_status("Rerun: index $i out of range.")
        return nothing
    end
    rec = _HISTORY[i]
    preset = joinpath(rec.workdir, "params.json")
    if isfile(preset)
        load_preset_path(preset)
        _set_status("Reloaded params from run '$(rec.case)' ($(rec.when))")
    else
        _set_status("Rerun: no params.json in $(rec.workdir) — only newer runs save one.")
    end
    return nothing
end
@qmlfunction rerun_from_history

# ---------------------------------------------------------------------------
# Movie generator — animates the bed evolution to MP4. Uses CairoMakie's
# record() (which calls ffmpeg under the hood). Skipped silently if
# ffmpeg isn't available; status pane reports the error.
# ---------------------------------------------------------------------------
function _make_movie(nc_path::AbstractString, mp4_path::AbstractString)
    NCDataset(nc_path, "r") do ds
        x  = Array(ds["x"])
        zb = Array(ds["zb"])
        nx, nt = size(zb)
        nt <= 1 && return ""
        hrms_var = haskey(ds, "hrms") ? Array(ds["hrms"]) : nothing

        zmin = minimum(zb) - 0.5
        zmax = maximum(zb) + 0.5

        fig = Figure(size=(900, 540), fontsize=14)
        ax = Axis(fig[1, 1];
            xlabel="Cross-shore distance (m)",
            ylabel="Elevation (m)",
            title="Bed evolution")
        ylims!(ax, zmin, zmax)
        bed_line = lines!(ax, x, view(zb, :, 1);
            color=:saddlebrown, linewidth=2.0)
        init_line = lines!(ax, x, view(zb, :, 1);
            color=(:gray, 0.5), linestyle=:dash, linewidth=1.0)
        title_obs = Observable("frame 1 / $nt")
        ax.title = title_obs

        # 20 fps; record over the time dimension. Total frames = nt.
        record(fig, mp4_path, 1:nt; framerate=20) do ti
            bed_line[2] = view(zb, :, ti)
            title_obs[] = "frame $ti / $nt"
        end
        return mp4_path
    end
end

# ---------------------------------------------------------------------------
# Main run handler — called from QML. Spawns the simulation on a worker
# thread so the Qt event loop keeps the UI alive.
# ---------------------------------------------------------------------------
function run_quick_sim()
    # Big visible banner — hard to miss in the launcher terminal.
    println("\n========== run_quick_sim CALLED ==========")
    flush(stdout)
    try
        if RUNNING[]
            _set_status("Already running. Wait for the current run to finish.")
            return nothing
        end
        RUNNING[] = true
        PROGRESS[] = 0.0
        CANCEL_REQ[] = false

        nthr = Threads.nthreads()
        println("[qml_gui] run_quick_sim entered; Threads.nthreads()=", nthr)
        flush(stdout)
        _set_status("run_quick_sim entered; nthreads=$nthr")

        # ---- Read all PARAMS on the MAIN thread before spawning ----------
        # Accessing JuliaPropertyMap from a worker thread blocks on Qt's
        # property-map mutex while the main thread is in exec(); the worker
        # silently hangs. Reading everything here, on the main thread, and
        # packing into a plain NamedTuple gives the worker self-contained
        # data it can use without any QML/Qt interaction.
        println("[qml_gui] reading PARAMS on main thread…")
        flush(stdout)
        p = _snapshot_params()
        println("[qml_gui] PARAMS snapshot OK: case=$(p.case_name), " *
                "profile=$(p.profile), nf=$(length(p.gs_mm)), " *
                "duration_h=$(p.duration_h)")
        flush(stdout)

        # NOTE on threading:
        # Threads.@spawn from inside a QML.jl callback queues a task that
        # never executes — Julia's scheduler doesn't pump tasks while Qt's
        # exec() holds the main thread, even with nthreads > 1. Tested
        # against Julia 1.10 + QML.jl 0.10.x on macOS. Until QML.jl exposes
        # a Qt-aware scheduler hook, the simulation is run synchronously
        # on the main thread. The UI freezes briefly during the run, but
        # status, progress, and the result plot all flush correctly when
        # control returns to Qt's event loop.
        _set_status("Running synchronously (UI will be unresponsive during the run).")
        _do_run(p)
    catch e
        msg = sprint(showerror, e, catch_backtrace())
        println("[qml_gui] DISPATCHER ERROR: ", first(msg, 800))
        flush(stdout)
        _set_status("DISPATCHER ERROR: " * first(msg, 800))
        RUNNING[] = false
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Read every PARAMS entry on the main thread and pack into a plain NamedTuple.
# This is the ONLY place the QML-bridged JuliaPropertyMap is touched per run.
# Step-by-step prints so we can pinpoint a hang on a specific key if it ever
# comes back.
# ---------------------------------------------------------------------------
function _snapshot_params()
    @inline function _read(key::String)
        # Verbose by design — comment out the println once stable.
        v = PARAMS[key]
        # println("[qml_gui]   PARAMS[$key] = $(repr(v))")
        return v
    end

    case_name        = string(_read("case_name"))
    profile          = String(_read("profile"))
    slope_rr         = Float64(_read("slope_rr"))
    depth_m          = Float64(_read("depth_m"))
    backshore_elev_m = Float64(_read("backshore_elev_m"))
    dune_m           = Float64(_read("dune_m"))
    dx_m             = Float64(_read("dx_m"))
    hrms        = Float64(_read("hrms"))
    tp          = Float64(_read("tp"))
    swl         = Float64(_read("swl"))
    duration_h  = Float64(_read("duration_h"))
    bc_dt_hours       = Float64(_read("bc_dt_hours"))
    output_interval_s = Float64(_read("output_interval_s"))
    effb        = Float64(_read("effb"))
    efff        = Float64(_read("efff"))
    blp         = Float64(_read("blp"))
    slp         = Float64(_read("slp"))
    tanphi      = Float64(_read("tanphi"))
    # Bound to a sensible range — fractional values are nonsense, very
    # large values waste compute. 0 disables; 30 is plenty for any case.
    n_pickup_smooth = clamp(round(Int, Float64(_read("n_pickup_smooth"))), 0, 30)
    bathy_csv   = string(_read("bathy_csv"))
    waves_csv   = string(_read("waves_csv"))
    thermal_csv = string(_read("thermal_csv"))
    thermal_on  = Bool(_read("thermal_on"))
    T_air_const   = Float64(_read("T_air_const"))
    T_water_const = Float64(_read("T_water_const"))
    outdir      = string(_read("outdir"))
    gs_mm_str   = string(_read("grain_sizes_mm"))
    gf_str      = string(_read("grain_fractions"))

    # Parse + validate the grain-size lists here too so any error shows up
    # on the main thread (where the user-facing dispatcher catch can see it).
    gs_mm = _parse_csv_floats(gs_mm_str)
    gf    = _parse_csv_floats(gf_str)
    isempty(gs_mm) && error("grain_sizes_mm has no values")
    isempty(gf)    && error("grain_fractions has no values")
    length(gs_mm) == length(gf) ||
        error("grain_sizes_mm ($(length(gs_mm)) values) and " *
              "grain_fractions ($(length(gf)) values) must match in length")
    if !isapprox(sum(gf), 1.0; atol=1e-3)
        gf = gf ./ sum(gf)
    end

    # NBS — unified dropdown, elevation-band based
    nbs_type       = String(_read("nbs_type"))
    nbs_z_min      = Float64(_read("nbs_z_min"))
    nbs_z_max      = Float64(_read("nbs_z_max"))
    nbs_density    = Float64(_read("nbs_density"))
    nbs_blade_w    = Float64(_read("nbs_blade_w"))
    nbs_height     = Float64(_read("nbs_height"))
    nbs_cd         = Float64(_read("nbs_cd"))
    nbs_crest_z    = Float64(_read("nbs_crest_z"))
    nbs_porosity   = Float64(_read("nbs_porosity"))
    nbs_stone_d    = Float64(_read("nbs_stone_d"))
    nbs_snow_depth = Float64(_read("nbs_snow_depth"))
    nbs_k_snow     = Float64(_read("nbs_k_snow"))
    nbs_max_depth  = Float64(_read("nbs_max_depth"))

    # Second NBS slot (full duplicate of nbs_* schema; "none" disables).
    nbs2_type       = String(_read("nbs2_type"))
    nbs2_z_min      = Float64(_read("nbs2_z_min"))
    nbs2_z_max      = Float64(_read("nbs2_z_max"))
    nbs2_density    = Float64(_read("nbs2_density"))
    nbs2_blade_w    = Float64(_read("nbs2_blade_w"))
    nbs2_height     = Float64(_read("nbs2_height"))
    nbs2_cd         = Float64(_read("nbs2_cd"))
    nbs2_crest_z    = Float64(_read("nbs2_crest_z"))
    nbs2_porosity   = Float64(_read("nbs2_porosity"))
    nbs2_stone_d    = Float64(_read("nbs2_stone_d"))
    nbs2_snow_depth = Float64(_read("nbs2_snow_depth"))
    nbs2_k_snow     = Float64(_read("nbs2_k_snow"))
    nbs2_max_depth  = Float64(_read("nbs2_max_depth"))

    # Forcing extras.
    slr_m         = Float64(_read("slr_m"))
    tide_on       = Bool(_read("tide_on"))
    tide_amp_m    = Float64(_read("tide_amp_m"))
    tide_period_h = Float64(_read("tide_period_h"))

    return (; case_name, profile, slope_rr, depth_m, backshore_elev_m, dune_m, dx_m,
              hrms, tp, swl, duration_h,
              bc_dt_hours, output_interval_s,
              effb, efff, blp, slp, tanphi,
              n_pickup_smooth,
              bathy_csv, waves_csv, thermal_csv,
              thermal_on, T_air_const, T_water_const,
              outdir,
              gs_mm, gf,
              nbs_type, nbs_z_min, nbs_z_max,
              nbs_density, nbs_blade_w, nbs_height, nbs_cd,
              nbs_crest_z, nbs_porosity, nbs_stone_d,
              nbs_snow_depth, nbs_k_snow, nbs_max_depth,
              nbs2_type, nbs2_z_min, nbs2_z_max,
              nbs2_density, nbs2_blade_w, nbs2_height, nbs2_cd,
              nbs2_crest_z, nbs2_porosity, nbs2_stone_d,
              nbs2_snow_depth, nbs2_k_snow, nbs2_max_depth,
              slr_m, tide_on, tide_amp_m, tide_period_h)
end

# The simulation body, factored out so it can run either inline or on a
# Threads.@spawn worker depending on thread availability above. Receives
# a plain NamedTuple `p` produced by `_snapshot_params()` on the main
# thread — this function never touches the QML-bridged PARAMS, only Julia
# Observables (which are safe across threads).
function _do_run(p::NamedTuple)
    println("[qml_gui] _do_run starting on thread ", Threads.threadid())
    flush(stdout)
    # Outer-scope refs so the catch block can still record the run in
    # history even if we error before workdir is built.
    case_for_hist    = p.case_name
    workdir_for_hist = ""
    t_started        = time()
    try
        case             = p.case_name
        case_for_hist    = case
        profile          = p.profile
        slope            = p.slope_rr
        depth            = p.depth_m
        backshore_elev   = p.backshore_elev_m
        dune             = p.dune_m
        dx               = p.dx_m
        hrms        = p.hrms
        tp          = p.tp
        swl         = p.swl
        duration_h  = p.duration_h
        effb        = p.effb
        efff        = p.efff
        blp         = p.blp
        slp         = p.slp
        tanphi      = p.tanphi
        bathy_csv   = p.bathy_csv
        waves_csv   = p.waves_csv
        outdir_root = p.outdir
        gs_mm       = p.gs_mm
        gf          = p.gf

        d50_m  = sum(gs_mm .* gf) * 1e-3
        d50_mm = d50_m * 1e3

        _set_status("Run started: $case  ($(profile))")

        stamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
        workdir = joinpath(outdir_root, "$(case)_$(stamp)")
        workdir_for_hist = workdir
        mkpath(workdir)
        _set_status("Work directory: $workdir")

        # ---- Bathymetry: external CSV or preset --------------------------
        x, z = if !isempty(bathy_csv)
            _set_status("Loading bathymetry from $bathy_csv")
            _load_bathy_csv(bathy_csv)
        elseif profile == "planar_beach"
            _planar_beach(depth, slope, backshore_elev, dx)
        else
            _beach_dune(depth, slope, backshore_elev, dune, dx)
        end
        _set_status("Profile: $(length(x)) nodes, " *
                    "x ∈ [$(round(x[1]; digits=1)), $(round(x[end]; digits=1))] m")

        # ---- Waves: external CSV or constant from form -------------------
        wf = if !isempty(waves_csv)
            _set_status("Loading wave forcing from $waves_csv")
            _load_waves_csv(waves_csv)
        else
            # bc_dt_hours controls how finely the constant time series is
            # sampled; smaller dt = more BC windows = finer model timestep
            # at the cost of more compute.
            constant_waves(; duration_days=duration_h / 24,
                             hrms=hrms, tp=tp, swl=swl,
                             dt_hours=p.bc_dt_hours)
        end

        # ---- Forcing overlays: SLR offset + sinusoidal tide --------------
        # SLR simply shifts every SWL value by p.slr_m. Tide adds a sine
        # wave with amplitude tide_amp_m and period tide_period_h. Both
        # apply on top of whatever swlbc is — constant or CSV.
        if p.slr_m != 0.0 || p.tide_on
            swl_extra = zeros(length(wf.timebc))
            if p.slr_m != 0.0
                swl_extra .+= p.slr_m
            end
            if p.tide_on && p.tide_period_h > 0
                omega = 2π / (p.tide_period_h * 3600.0)
                swl_extra .+= p.tide_amp_m .* sin.(omega .* wf.timebc)
            end
            wf = merge(wf, (; swlbc=wf.swlbc .+ swl_extra))
            _set_status("Forcing overlays: SLR=+$(p.slr_m) m" *
                        (p.tide_on ? ", tide A=$(p.tide_amp_m) m / T=$(p.tide_period_h) h" : ""))
        end

        # ---- NBS dispatch (elevation-band based) --------------------------
        # The user picks up to TWO NBS slots; each is configured by its own
        # bed-elevation band [z_min, z_max] (m rel. SWL) plus type-specific
        # parameters. Outside any band, no feature is applied. This is
        # geometry-portable — the same preset places "marsh in the swash
        # zone" on any profile.
        #
        # The two slots are processed sequentially. Vegetation results merge
        # (combined VegetationInput); structural results merge cell-by-cell
        # (take whichever crest is higher).
        veg_input         = nothing   # → VegetationInput (iveg=1)
        porous_z_vec      = nothing   # → porous floor (iperm=1)
        hardbottom_z_vec  = nothing   # → impermeable cap (isedav=1)
        snow_extra_kwargs = NamedTuple()  # → SnowConfig + snow_depth time series
        z_was_modified    = false      # if true, also bump z before build_config

        # Helper: mask of nodes whose BED elevation is inside [z_min, z_max].
        # Returns (mask, n_inside).
        _elev_mask(zvec, z_min, z_max) = begin
            m = (zvec .>= z_min) .& (zvec .<= z_max)
            (m, count(m))
        end

        # ---- Per-slot processor ---------------------------------------
        # Pulls the type + type-specific fields from `p` for the given
        # slot prefix ("" for slot 1, "2" for slot 2). Updates the
        # outer-scope structures.
        function _apply_slot!(slot_label::String, nbs_type::String,
                              z_min::Float64, z_max::Float64,
                              density::Float64, blade_w::Float64,
                              height::Float64, cd::Float64,
                              crest_z::Float64, porosity::Float64,
                              stone_d::Float64,
                              snow_depth::Float64, k_snow::Float64,
                              max_depth::Float64)
            nbs_type == "none" && return
            mask, n_in = _elev_mask(z, z_min, z_max)
            if n_in == 0 && nbs_type != "snow"
                _set_status("NBS '$nbs_type' (slot $slot_label): band " *
                            "z ∈ [$z_min, $z_max] m contains no nodes — skipped.")
                return
            end
            np = length(x)
            if nbs_type in ("marsh", "dune_grass", "kelp")
                # Add to (or initialize) VegetationInput.
                # When merging two veg slots into one input, take the
                # per-node max so each band gets its own density.
                vn  = veg_input === nothing ? zeros(np) : vec(veg_input.vegn[:,1])
                vb  = veg_input === nothing ? zeros(np) : vec(veg_input.vegb[:,1])
                vd  = veg_input === nothing ? zeros(np) : vec(veg_input.vegd[:,1])
                vh  = veg_input === nothing ? zeros(np) : vec(veg_input.vegh[:,1])
                vrd = veg_input === nothing ? zeros(np) : vec(veg_input.vegrd[:,1])
                vrh = veg_input === nothing ? zeros(np) : vec(veg_input.vegrh[:,1])
                for j in 1:np
                    if mask[j]
                        vn[j]  = max(vn[j],  density)
                        vb[j]  = max(vb[j],  blade_w)
                        vd[j]  = max(vd[j],  height)
                        vh[j]  = max(vh[j],  height)
                        vrd[j] = max(vrd[j], 0.10)
                        vrh[j] = max(vrh[j], 0.05)
                    end
                end
                veg_input = VegetationInput(;
                    vegcd = cd, vegcdm = cd,
                    vegn  = reshape(vn,  np, 1),
                    vegb  = reshape(vb,  np, 1),
                    vegd  = reshape(vd,  np, 1),
                    vegh  = reshape(vh,  np, 1),
                    vegrd = reshape(vrd, np, 1),
                    vegrh = reshape(vrh, np, 1))
                _set_status("NBS $nbs_type (slot $slot_label): $n_in nodes, " *
                            "density=$density/m², h=$height m, Cd=$cd")

            elseif nbs_type == "oyster_reef"
                z_orig = copy(z)
                if porous_z_vec === nothing
                    porous_z_vec = copy(z) .- 1.0e3
                end
                z_new = copy(z)
                for j in 1:np
                    if mask[j]
                        z_new[j]        = max(z[j], crest_z)
                        # Use the LOWER porous floor between slots so
                        # structure thickness is maximised when both are
                        # active at the same node.
                        porous_z_vec[j] = min(porous_z_vec[j], z_orig[j])
                    end
                end
                z = z_new
                z_was_modified = true
                _set_status("NBS oyster_reef (slot $slot_label): $n_in nodes, " *
                            "crest z=$crest_z m, porosity=$porosity, " *
                            "Dn50=$stone_d m")

            elseif nbs_type == "breakwater"
                if hardbottom_z_vec === nothing
                    hardbottom_z_vec = copy(z) .- 1.0e3
                end
                z_new = copy(z)
                for j in 1:np
                    if mask[j]
                        z_new[j]            = max(z[j], crest_z)
                        # Take the HIGHER crest if two breakwaters overlap.
                        hardbottom_z_vec[j] = max(hardbottom_z_vec[j], crest_z)
                    end
                end
                z = z_new
                z_was_modified = true
                _set_status("NBS breakwater (slot $slot_label): $n_in nodes, " *
                            "elevated to z=$crest_z m, impermeable hardbottom.")

            elseif nbs_type == "snow"
                if !p.thermal_on && isempty(snow_extra_kwargs)
                    _set_status("NBS snow (slot $slot_label): auto-enabling " *
                                "thermal model.")
                end
                snow_cfg = SnowConfig(; k_snow=k_snow, max_depth=max_depth)
                ntimes = length(wf.timebc)
                snow_extra_kwargs = (; snow=snow_cfg,
                                       snow_depth=fill(snow_depth, ntimes))
                _set_status("NBS snow (slot $slot_label): depth=$snow_depth m, " *
                            "k_snow=$k_snow W/m/K, cap=$max_depth m")
            else
                _set_status("NBS slot $slot_label: unknown type '$nbs_type'.")
            end
            return
        end

        # Run both NBS slots through the same helper.
        _apply_slot!("1", p.nbs_type, p.nbs_z_min, p.nbs_z_max,
                     p.nbs_density, p.nbs_blade_w, p.nbs_height, p.nbs_cd,
                     p.nbs_crest_z, p.nbs_porosity, p.nbs_stone_d,
                     p.nbs_snow_depth, p.nbs_k_snow, p.nbs_max_depth)
        _apply_slot!("2", p.nbs2_type, p.nbs2_z_min, p.nbs2_z_max,
                     p.nbs2_density, p.nbs2_blade_w, p.nbs2_height, p.nbs2_cd,
                     p.nbs2_crest_z, p.nbs2_porosity, p.nbs2_stone_d,
                     p.nbs2_snow_depth, p.nbs2_k_snow, p.nbs2_max_depth)
        # We also need a top-level `nbs_type` variable for the OptionFlags
        # snow detection below — "snow" in either slot counts.
        nbs_type = (p.nbs_type == "snow" || p.nbs2_type == "snow") ?
                   "snow" : p.nbs_type

        # ---- Sediment + multifraction ------------------------------------
        sed = make_sediment(; d50=d50_m,
                              effb=effb, efff=efff,
                              blp=blp, slp=slp, tanphi=tanphi)
        # Build MultifractionConfig directly so we can override
        # n_pickup_smooth (sediment_custom doesn't expose it).
        mf  = MultifractionConfig(grain_sizes=gs_mm .* 1e-3,
                                   initial_fractions=gf,
                                   n_pickup_smooth=p.n_pickup_smooth)

        nf = length(gs_mm)
        if nf > 1
            _set_status("Multifraction: $(nf) grains [$(join(gs_mm, ", ")) mm], " *
                        "fractions [$(join(round.(gf; digits=3), ", "))], " *
                        "weighted d50=$(round(d50_mm; digits=3)) mm")
        end

        # ---- Thermal / permafrost (optional) -----------------------------
        # When enabled, build a ThermalConfig + (T_air, T_water) time series
        # AND set OptionFlags(isedav=1) so the wave/sediment solver actually
        # reads the thermal-computed `zb_hard` floor. Without isedav=1 the
        # hardbottom-clamp pass is a no-op, ALT evolves but never cuts off
        # erosion — which is what the user observed.
        thermal_kwargs = NamedTuple()
        # OptionFlags is rebuilt below as the union of all enabled submodels.
        # Thermal hardbottom OR a breakwater hardbottom both set isedav=1.
        # Snow NBS implicitly turns on the thermal submodel.
        thermal_on_eff = p.thermal_on || nbs_type == "snow"
        isedav_flag = (thermal_on_eff || hardbottom_z_vec !== nothing) ? 1 : 0
        iveg_flag   = veg_input    !== nothing ? 1 : 0
        iperm_flag  = porous_z_vec !== nothing ? 1 : 0
        opts_kwargs = (; options=OptionFlags(isedav=isedav_flag,
                                              iveg=iveg_flag,
                                              iperm=iperm_flag))

        if thermal_on_eff
            tcfg = ThermalConfig()

            if !isempty(p.thermal_csv)
                _set_status("Loading thermal forcing from $(p.thermal_csv)")
                tdata = _load_thermal_csv(p.thermal_csv)
                thermal_kwargs = (; thermal=tcfg,
                                    thermal_time=tdata.time,
                                    T_air=tdata.T_air,
                                    T_water=tdata.T_water)
                if !isempty(tdata.snow_depth)
                    thermal_kwargs = merge(thermal_kwargs,
                                            (; snow_depth=tdata.snow_depth))
                end
                _set_status("Thermal: $(length(tdata.time)) samples; " *
                            "T_air ∈ [$(round(minimum(tdata.T_air); digits=1)), " *
                            "$(round(maximum(tdata.T_air); digits=1))] °C")
            else
                # Constant values broadcast across all BC times.
                ntimes = length(wf.timebc)
                thermal_kwargs = (; thermal=tcfg,
                                    T_air=fill(p.T_air_const, ntimes),
                                    T_water=fill(p.T_water_const, ntimes))
                _set_status("Thermal: constant T_air=$(p.T_air_const) °C, " *
                            "T_water=$(p.T_water_const) °C")
            end
            _set_status("Thermal: isedav=1 (hardbottom enforcement enabled)")
        end

        nbs_kwargs = NamedTuple()
        if veg_input !== nothing
            nbs_kwargs = merge(nbs_kwargs, (; vegetation=veg_input))
        end
        if porous_z_vec !== nothing
            nbs_kwargs = merge(nbs_kwargs, (; porous_z=porous_z_vec,
                                              porosity=p.nbs_porosity,
                                              stone_diameter=p.nbs_stone_d))
        end
        if hardbottom_z_vec !== nothing
            nbs_kwargs = merge(nbs_kwargs, (; hardbottom_z=hardbottom_z_vec))
        end

        cfg = build_config(; bathymetry_x=x, bathymetry_z=z,
                             multifraction=mf, sediment=sed,
                             dx=dx, wf...,
                             opts_kwargs..., thermal_kwargs...,
                             nbs_kwargs..., snow_extra_kwargs...)

        ncname = "$(case).nc"
        ncpath = joinpath(workdir, ncname)

        # Output interval: explicit user value if > 0, otherwise auto
        # (~24 frames over the run, floored at 60 s).
        total_s = wf.timebc[end] - wf.timebc[1]
        interval = p.output_interval_s > 0 ?
                   p.output_interval_s :
                   max(60.0, total_s / 24)

        _set_status("Running simulation… $(length(wf.timebc) - 1) BC window(s)")
        t0 = time()

        # ---- Progress callback updates the UI Observables every BC window
        # and aborts the run if the user clicked Cancel.
        function _on_step(itime, n_steps, _state)
            frac = itime / n_steps
            elapsed = time() - t0
            eta = itime > 0 ? elapsed * (n_steps - itime) / itime : 0.0
            PROGRESS[] = frac
            ELAPSED[] = "step $itime / $n_steps · " *
                        "elapsed $(round(elapsed; digits=1))s · " *
                        "ETA $(round(eta; digits=1))s"
            if CANCEL_REQ[]
                error("Run cancelled by user at step $itime / $n_steps")
            end
            return nothing
        end

        run_simulation!(cfg; outdir=workdir, outfile=ncname,
                        output_interval_s=interval,
                        progress_callback=_on_step)

        dt = round(time() - t0; digits=1)
        RESULT_PATH[] = ncpath
        PROGRESS[] = 1.0
        _set_status("Simulation OK in $(dt)s. Generating plot…")

        png_path = joinpath(workdir, "result.png")
        url = _make_plot(ncpath, png_path)
        PLOT_PATH[] = url

        # Persist this run's PARAMS snapshot next to the NetCDF so the
        # Run-history panel can offer a one-click "reload these params".
        try
            open(joinpath(workdir, "params.json"), "w") do io
                JSON3.pretty(io, _params_snapshot_dict())
            end
        catch e
            @warn "qml_gui: could not save per-run params.json" exception=e
        end

        # ---- Volume summary ------------------------------------------
        try
            vs = _volume_summary(ncpath, swl)
            VOLUME_TEXT[] = string(
                "Net Δvol: $(round(vs.v_tot; digits=2)) m³/m  ",
                "(landward: $(round(vs.v_land; digits=2)),  ",
                "seaward: $(round(vs.v_sea; digits=2)))\n",
                "Max erosion: $(round(vs.max_ero; digits=3)) m at ",
                "x=$(round(vs.x_ero; digits=1)) m\n",
                "Max deposition: $(round(vs.max_dep; digits=3)) m at ",
                "x=$(round(vs.x_dep; digits=1)) m")
        catch e
            VOLUME_TEXT[] = "(volume summary unavailable: " *
                            sprint(showerror, e) * ")"
        end

        # ---- Optional MP4 movie --------------------------------------
        if MAKE_MOVIE[]
            mp4_path = joinpath(workdir, "bed_evolution.mp4")
            try
                _set_status("Generating MP4 movie (requires ffmpeg)…")
                _make_movie(ncpath, mp4_path)
                MOVIE_PATH[] = "file://" * mp4_path
                _set_status("Movie: $mp4_path")
            catch e
                _set_status("Movie failed (ffmpeg missing?): " *
                            sprint(showerror, e))
            end
        end

        _set_status("Done. NetCDF: $ncpath\nPlot: $png_path")
        _record_run(case, workdir, dt, true)
        # Snapshot params used by this run so the param-diff panel can
        # compare future edits against it.
        _LAST_RUN_PARAMS[] = _params_snapshot_dict()
        autosave_session()

    catch e
        # Print the FULL stack trace to the terminal (untruncated) — the
        # GUI status pane gets a truncated copy but the launcher terminal
        # has the canonical record we need to debug.
        msg = sprint(showerror, e, catch_backtrace())
        println("\n========== _do_run ERROR ==========")
        println(msg)
        println("===================================\n")
        flush(stdout)
        _set_status("ERROR: " * first(msg, 1200))
        _record_run(case_for_hist, workdir_for_hist, time() - t_started, false)
    finally
        RUNNING[] = false
        CANCEL_REQ[] = false
    end
    return nothing
end
@qmlfunction run_quick_sim

# ---------------------------------------------------------------------------
# Helper for the QML FileDialog wiring — QML calls this with the dialog's
# selected URL ("file:///...") and a target PARAMS key. We strip the
# scheme and write the path into PARAMS so the next run picks it up.
# ---------------------------------------------------------------------------
function set_csv_path(key::AbstractString, url::AbstractString)
    p = url
    p = startswith(p, "file://") ? p[8:end] : p
    p = replace(p, r"%20" => " ")
    PARAMS[String(key)] = p
    _set_status("$key → $(isempty(p) ? "(cleared)" : p)")
    return nothing
end
@qmlfunction set_csv_path

function clear_csv_path(key::AbstractString)
    PARAMS[String(key)] = ""
    _set_status("$key cleared (will use form values)")
    return nothing
end
@qmlfunction clear_csv_path

# ---------------------------------------------------------------------------
# Default parameter values + UI state, both as JuliaPropertyMaps.
# ---------------------------------------------------------------------------
const PARAMS = JuliaPropertyMap(
    # Case + bathymetry — profile length is COMPUTED from depth, slope,
    # and backshore elevation; the user controls geometry, not extent.
    "case_name"        => "quickrun",
    "profile"          => "beach_dune",
    "slope_rr"         => 0.05,
    "depth_m"          => 8.0,    # offshore depth below SWL (positive)
    "backshore_elev_m" => 2.0,    # landward elevation above SWL (positive)
    "dune_m"           => 4.0,    # dune height ABOVE backshore (beach+dune)
    "dx_m"             => 1.0,
    # Waves
    "hrms"       => 1.0,
    "tp"         => 8.0,
    "swl"        => 0.5,
    "duration_h" => 12.0,
    # Timing
    "bc_dt_hours"      => 1.0,    # constant_waves time-step (form mode only)
    "output_interval_s"=> 0.0,    # 0 = auto = max(60, total/24)
    # Sediment — comma-separated lists for multi-fraction runs.
    # For a single-grain run leave it as "0.30" / "1.0".
    "grain_sizes_mm"  => "0.30",
    "grain_fractions" => "1.0",
    # Free model parameters (defaults match SedimentConfig)
    "effb"       => 0.005,
    "efff"       => 0.005,
    "blp"        => 0.002,
    "slp"        => 0.2,
    "tanphi"     => 0.63,
    # Transport smoothing — passes of the 1-2-1 kernel applied to the
    # per-cell pickup field after the divergence step. More passes =
    # smoother bed evolution but more computational cost. Map to
    # MultifractionConfig.n_pickup_smooth. Default 10 matches the
    # built-in CSHORE.jl default.
    "n_pickup_smooth" => 10,
    # External CSVs (empty = use form)
    "bathy_csv"      => "",
    "waves_csv"      => "",
    "thermal_csv"    => "",
    # Thermal / permafrost (defaults to off)
    "thermal_on"     => false,
    "T_air_const"    => -5.0,    # °C — winter Arctic baseline
    "T_water_const"  =>  0.0,    # °C — submerged ground at freezing
    # I/O
    "outdir"     => joinpath(REPO_DIR, "runs"),

    # ----- Nature-Based / Hybrid Infrastructure (NBS) ----------------------
    # One unified dropdown selects the NBS type; the band is defined by
    # bed-elevation limits (z_min, z_max relative to SWL) rather than
    # x-coordinates so the same preset works on any profile geometry.
    #
    # nbs_type values:
    #   "none"       — disabled
    #   "marsh"      — emergent vegetation (Spartina-style); VegetationInput
    #   "dune_grass" — backshore/dune vegetation; VegetationInput
    #   "kelp"       — submerged canopy in deeper water; VegetationInput
    #   "oyster_reef"— porous structure with crest elevation; PorousInput
    #   "breakwater" — impermeable hardbottom with crest elevation; hardbottom_z
    #   "snow"       — winter snow cover (constant depth); SnowConfig + snow_depth
    "nbs_type"       => "none",
    "nbs_z_min"      => -0.5,
    "nbs_z_max"      =>  1.0,
    # Vegetation-type fields
    "nbs_density"    => 200.0,
    "nbs_blade_w"    => 0.006,
    "nbs_height"     => 0.30,
    "nbs_cd"         => 1.0,
    # Porous / breakwater fields
    "nbs_crest_z"    => -0.5,
    "nbs_porosity"   => 0.40,
    "nbs_stone_d"    => 0.05,
    # Snow-cover fields
    "nbs_snow_depth" => 0.10,    # constant snow depth (m)
    "nbs_k_snow"     => 0.15,    # thermal conductivity (W/m/K)
    "nbs_max_depth"  => 0.30,    # cap (m) for degree-day model

    # ----- Second NBS slot (hybrid combinations) -------------------------
    # Allows stacking a second feature in a different elevation band, e.g.
    # marsh in the swash zone + offshore oyster reef. Same schema as
    # nbs_* above; "none" disables.
    "nbs2_type"       => "none",
    "nbs2_z_min"      => -2.0,
    "nbs2_z_max"      => -0.5,
    "nbs2_density"    => 200.0,
    "nbs2_blade_w"    => 0.006,
    "nbs2_height"     => 0.30,
    "nbs2_cd"         => 1.0,
    "nbs2_crest_z"    => -0.5,
    "nbs2_porosity"   => 0.40,
    "nbs2_stone_d"    => 0.05,
    "nbs2_snow_depth" => 0.10,
    "nbs2_k_snow"     => 0.15,
    "nbs2_max_depth"  => 0.30,

    # ----- Forcing extras: sea-level rise + sinusoidal tide --------------
    # slr_m simply offsets every SWL value (constant or CSV) by this
    # amount — convenient "what-if +0.5 m SLR" slider.
    # When tide_on=true, an additional sinusoidal tide is overlaid:
    #   swl_extra(t) = tide_amp_m * sin(2π t / (tide_period_h * 3600))
    "slr_m"           => 0.0,
    "tide_on"         => false,
    "tide_amp_m"      => 0.5,
    "tide_period_h"   => 12.42,    # M2 semi-diurnal
)

const UI = JuliaPropertyMap(
    "statusMsg"   => STATUS,
    "running"     => RUNNING,
    "progress"    => PROGRESS,
    "elapsed"     => ELAPSED,
    "resultPath"  => RESULT_PATH,
    "plotPath"    => PLOT_PATH,
    "history"     => HISTORY_TEXT,
    "lastWorkdir" => LAST_WORKDIR,
    "cancelReq"   => CANCEL_REQ,
    "validation"  => VALIDATION_TEXT,
    "volume"      => VOLUME_TEXT,
    "moviePath"   => MOVIE_PATH,
    "previewPath" => PREVIEW_PATH,
    "runtimeEst"  => RUNTIME_EST,
    "paramDiff"   => PARAM_DIFF,
    "makeMovie"   => MAKE_MOVIE,
    "presets"     => PRESETS_TEXT,
)

# ---------------------------------------------------------------------------
# Load QML and start the Qt event loop.
# ---------------------------------------------------------------------------
const QML_FILE = joinpath(SCRIPT_DIR, "cshore_main.qml")
isfile(QML_FILE) || error("missing $QML_FILE")

println("[qml_gui] threads available: $(Threads.nthreads())")
println("[qml_gui] starting GUI; output dir default: $(PARAMS["outdir"])")
# Restore the last session before loading QML so the form picks up the
# user's previous parameter values automatically. Idempotent if no file.
_load_last_session!()
loadqml(QML_FILE; params=PARAMS, ui=UI)

exec()
