#!/usr/bin/env julia
# =============================================================================
# web_gui/app.jl — LEAN variant
#
# A bare HTTP.jl server with handwritten HTML/JS instead of Stipple/Genie/
# StippleUI. The client posts form parameters as JSON to /run; the server
# runs the CSHORE simulation, reads the result NetCDF, and replies with a
# JSON blob containing the arrays the browser needs to draw Plotly charts.
#
#   GET  /          → static index.html (Plotly via CDN, plain HTML form)
#   POST /run       → run simulation, return JSON {plot_x, plot_z_snaps,
#                     plot_dz_rows, mix_*, volume_text, status, ...}
#   GET  /healthz   → 200 OK, used by Fly's healthcheck
#
# Run locally:
#   julia --project=web_gui --threads=auto web_gui/app.jl
#   then open http://localhost:8000
# =============================================================================

const SCRIPT_DIR = @__DIR__
const REPO_DIR   = dirname(SCRIPT_DIR)

# In the Fly container the Dockerfile already ran `Pkg.develop` + `Pkg.instantiate`
# during the build, so the project is pinned, all artifacts are downloaded, and
# the manifest is current. Running them again at every server start just pays
# the cost of loading Pkg, walking the resolver, and re-stat'ing the manifest —
# ~10-20 s on each boot. Skip it when a manifest is already present (the
# container case). Locally, a fresh `julia --project=web_gui app.jl` from a
# clean clone has no manifest, so we fall back to the original behavior to
# preserve the documented dev workflow.
if isfile(joinpath(SCRIPT_DIR, "Manifest.toml"))
    # Pkg.activate is implicit via JULIA_PROJECT (set in the Dockerfile + when
    # invoked with `--project=web_gui`), so no further setup needed here.
else
    using Pkg
    Pkg.activate(SCRIPT_DIR)
    try
        Pkg.develop(path=REPO_DIR)
    catch e
        @warn "web_gui: Pkg.develop raised — continuing" exception=e
    end
    try
        Pkg.instantiate()
    catch
        Pkg.resolve(); Pkg.instantiate()
    end
end

using HTTP
using JSON3
using NCDatasets
using CSHORE
using Dates

# ---------------------------------------------------------------------------
# Bathymetry preset helpers (same shapes as qml_gui).
# ---------------------------------------------------------------------------
function _planar_beach(depth, slope, backshore_elev, dx)
    slope_eff = max(slope, 1e-6)
    L = (depth + backshore_elev) / slope_eff
    n = max(2, ceil(Int, L / dx) + 1)
    xs = [(i - 1) * dx for i in 1:n]
    zs = [-depth + slope_eff * x for x in xs]
    zs[end] = backshore_elev
    return xs, zs
end

function _beach_dune(depth, slope, backshore_elev, dune_height, dx)
    slope_eff = max(slope, 1e-6)
    sigma  = 15.0
    x_toe  = (depth + backshore_elev) / slope_eff
    x_crest = x_toe + 1.5 * sigma
    L      = x_toe + 4.0 * sigma
    n = max(2, ceil(Int, L / dx) + 1)
    xs = [(i - 1) * dx for i in 1:n]
    zs = Vector{Float64}(undef, n)
    @inbounds for (i, x) in enumerate(xs)
        z_trend = x <= x_toe ? -depth + slope_eff * x : backshore_elev
        z_bump  = dune_height * exp(-((x - x_crest)^2) / (2 * sigma^2))
        zs[i]   = z_trend + z_bump
    end
    return xs, zs
end

# ---------------------------------------------------------------------------
# Static asset: the HTML/JS UI. Served from /public/index.html on disk.
# ---------------------------------------------------------------------------
const PUBLIC_DIR = joinpath(SCRIPT_DIR, "public")
isdir(PUBLIC_DIR) || mkpath(PUBLIC_DIR)
const INDEX_HTML_PATH = joinpath(PUBLIC_DIR, "index.html")

function load_index_html()
    if isfile(INDEX_HTML_PATH)
        return read(INDEX_HTML_PATH, String)
    end
    return "<h1>web_gui/public/index.html missing — see app.jl docstring</h1>"
end

# ---------------------------------------------------------------------------
# Hard caps to prevent compute-bombs. Reject any run whose cost-proxy
# (duration_h / dx_m²) exceeds a threshold — that's roughly the wall-time
# floor on the smallest VM.
# ---------------------------------------------------------------------------
const MAX_COST_PROXY = 50_000.0   # tuned for shared-cpu-2x; ~10 min worst case

function _check_cost(params)
    # Effective duration: from wave CSV if present, otherwise from form.
    wave_csv = get(params, :wave_csv, nothing)
    dh = if wave_csv !== nothing && haskey(wave_csv, :t)
        t = wave_csv[:t]
        n = length(t)
        n > 10_000 && error("wave_csv has $n samples (max 10000)")
        n < 2 ? 0.0 : (Float64(last(t)) - Float64(first(t))) / 3600.0
    else
        Float64(get(params, :duration_h, 0))
    end

    # Effective dx: from profile CSV spacing if present, otherwise from form.
    prof_csv = get(params, :profile_csv, nothing)
    dx = if prof_csv !== nothing && haskey(prof_csv, :x)
        xc = prof_csv[:x]
        n = length(xc)
        n > 5_000 && error("profile_csv has $n points (max 5000)")
        n < 2 ? 1.0 : (Float64(last(xc)) - Float64(first(xc))) / max(1, n - 1)
    else
        Float64(get(params, :dx_m, 1.0))
    end
    dx = max(dx, 0.01)

    cost = dh / (dx^2)
    if cost > MAX_COST_PROXY
        error("Run rejected: cost proxy duration_h / dx² = $(round(cost; digits=1)) " *
              "exceeds MAX_COST_PROXY ($MAX_COST_PROXY). " *
              "Pick a shorter duration_h or a coarser dx_m.")
    end
end

# ---------------------------------------------------------------------------
# Build the CSHORE config from the form params + run the simulation +
# read the NetCDF + pack the arrays the browser needs into a plain Dict
# that JSON3 will serialize.
# ---------------------------------------------------------------------------
function run_simulation_json(params)
    _check_cost(params)

    # Snapshot params with defaults for any missing field.
    g(k, d) = get(params, k, d)
    case_       = String(g(:case_name, "webrun"))
    profile_    = String(g(:profile,   "beach_dune"))
    depth_      = Float64(g(:depth_m,           8.0))
    slope_      = Float64(g(:slope_rr,          0.05))
    bs_         = Float64(g(:backshore_elev_m,  2.0))
    dune_       = Float64(g(:dune_m,            4.0))
    dx_         = Float64(g(:dx_m,              1.0))
    hrms_       = Float64(g(:hrms,              1.0))
    tp_         = Float64(g(:tp,                8.0))
    swl_        = Float64(g(:swl,               0.5))
    dur_h_      = Float64(g(:duration_h,       12.0))
    tide_on_    = Bool(g(:tide_on,             false))
    tide_amp_   = Float64(g(:tide_amp_m,        0.5))
    tide_per_   = Float64(g(:tide_period_h,    12.42))
    storm_on_     = Bool(g(:storm_on,              false))
    storm_hrms_   = Float64(g(:storm_hrms_peak,    3.0))
    storm_recur_d = Float64(g(:storm_recurrence_d, 30.0))
    gs_str_     = String(g(:grain_sizes_mm,    "0.30"))
    gf_str_     = String(g(:grain_fractions,   "1.0"))
    thermal_on_ = Bool(g(:thermal_on,          false))
    T_air_      = Float64(g(:T_air_C,          -5.0))
    T_water_    = Float64(g(:T_water_C,         0.0))
    T_init_     = Float64(g(:T_init_C,         -2.0))
    nbs_kind    = String(g(:nbs_type,          "none"))
    nbs_zmin_   = Float64(g(:nbs_z_min,        -0.5))
    nbs_zmax_   = Float64(g(:nbs_z_max,         1.0))
    nbs_crest_  = Float64(g(:nbs_crest_z,      -0.5))
    nbs_dens_   = Float64(g(:nbs_density,     200.0))
    nbs_h_      = Float64(g(:nbs_height,        0.3))
    nbs_cd_     = Float64(g(:nbs_cd,            1.0))
    nbs_por_    = Float64(g(:nbs_porosity,      0.4))
    nbs_stone_  = Float64(g(:nbs_stone_d,       0.05))
    nbs_thick_  = Float64(g(:nbs_thickness,     1.0))   # gravel-revetment layer thickness (m)
    nbs_vol_    = Float64(g(:nbs_volume_m3,    20.0))   # beach-nourishment volume per alongshore meter
    effb_       = Float64(g(:effb,              0.005))
    efff_       = Float64(g(:efff,              0.005))
    blp_        = Float64(g(:blp,               0.002))
    slp_        = Float64(g(:slp,               0.2))
    tanphi_     = Float64(g(:tanphi,            0.63))

    stamp   = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    workdir = joinpath(REPO_DIR, "runs", "$(case_)_$(stamp)_lean")
    mkpath(workdir)

    # ---- Profile: CSV upload overrides the preset --------------------
    prof_csv = g(:profile_csv, nothing)
    x, z = if prof_csv !== nothing && haskey(prof_csv, :x) && haskey(prof_csv, :z)
        xc = Float64.(collect(prof_csv[:x]))
        zc = Float64.(collect(prof_csv[:z]))
        length(xc) == length(zc) ||
            error("profile_csv: x and z column lengths differ ($(length(xc)) vs $(length(zc)))")
        length(xc) >= 2 || error("profile_csv: need at least 2 points")
        all(isfinite, xc) && all(isfinite, zc) ||
            error("profile_csv: x and z must be finite numbers")
        # Override the form's dx with the CSV's actual mean spacing so
        # build_config's metadata matches the supplied grid. The user
        # is responsible for providing uniformly-spaced data.
        dx_ = (xc[end] - xc[1]) / (length(xc) - 1)
        xc, zc
    elseif profile_ == "planar_beach"
        _planar_beach(depth_, slope_, bs_, dx_)
    else
        _beach_dune(depth_, slope_, bs_, dune_, dx_)
    end

    # ---- Structure dispatch -----------------------------------------
    veg_input        = nothing
    porous_z_vec     = nothing
    hardbottom_z_vec = nothing
    if nbs_kind != "none"
        mask = (z .>= nbs_zmin_) .& (z .<= nbs_zmax_)
        np   = length(x); n_in = count(mask)
        if n_in > 0
            if nbs_kind in ("low_marsh", "high_marsh", "marsh_shrub",
                            "phragmites", "mangrove",
                            "dune_grass", "dune_shrub",
                            "beach_grass", "beach_forbs",
                            "eelgrass", "kelp")
                # All vegetation kinds share the same physics — the difference
                # is the defaults pre-filled in the form (density, height,
                # Cd, z-band). Per-species defaults come from
                # `NBS_VEGETATION_PARAMS` in src/presets/nbs.jl. Stem
                # diameter is held constant here since it has minor
                # influence relative to density × height.
                vn  = zeros(np); vb = fill(0.006, np); vd = fill(nbs_h_, np)
                vh  = fill(nbs_h_, np); vrd = fill(0.10, np); vrh = fill(0.05, np)
                vn[mask]  .= nbs_dens_
                veg_input = VegetationInput(;
                    vegcd=nbs_cd_, vegcdm=nbs_cd_,
                    vegn=reshape(vn, np, 1), vegb=reshape(vb, np, 1),
                    vegd=reshape(vd, np, 1), vegh=reshape(vh, np, 1),
                    vegrd=reshape(vrd, np, 1), vegrh=reshape(vrh, np, 1))
            elseif nbs_kind in ("gravel_revetment", "rock_revetment")
                # Porous layer of `nbs_thickness` below the bed surface
                # across the z-band. Same physics for gravel (small Dn50,
                # thinner layer) and rock revetment (large Dn50, thicker);
                # only the default parameters differ at the UI layer.
                # Bed elevation unchanged; CSHORE treats the porous_z
                # surface as the underside of the permeable layer.
                # Outside the band porous_z is set far below the bed so
                # the layer effectively has zero thickness there.
                porous_z_vec = copy(z) .- 1.0e3
                for j in 1:np
                    if mask[j]
                        porous_z_vec[j] = z[j] - nbs_thick_
                    end
                end
            elseif nbs_kind == "breakwater"
                hardbottom_z_vec = copy(z) .- 1.0e3
                z_new = copy(z)
                for j in 1:np
                    if mask[j]
                        z_new[j] = max(z[j], nbs_crest_)
                        hardbottom_z_vec[j] = nbs_crest_
                    end
                end
                z = z_new
            elseif nbs_kind == "beach_nourishment"
                # Half-cosine taper that's full thickness at the seaward
                # edge of the band and goes smoothly to zero at the
                # landward edge, so the placed sand merges into the
                # upper beach without a step. ∫ 0.5*(1+cos(π t)) dt
                # from 0→1 = 0.5, so dz_max = 2 * vol / width gives the
                # same cross-shore-integrated added area as the
                # requested volume.
                idx = findall(mask)
                x_lo = x[first(idx)]; x_hi = x[last(idx)]
                width = x_hi - x_lo
                if width > 0
                    dz_max = 2 * nbs_vol_ / width
                    z = copy(z)
                    for j in idx
                        t = (x[j] - x_lo) / width
                        z[j] += dz_max * 0.5 * (1 + cos(pi * t))
                    end
                end
            end
        end
    end

    # ---- Waves: CSV upload overrides constant_waves + overlays -------
    wave_csv = g(:wave_csv, nothing)
    custom_waves = wave_csv !== nothing && haskey(wave_csv, :t) && haskey(wave_csv, :hrms)
    wf = if custom_waves
        tcsv = Float64.(collect(wave_csv[:t]))
        hcsv = Float64.(collect(wave_csv[:hrms]))
        pcsv = haskey(wave_csv, :tp)  ? Float64.(collect(wave_csv[:tp]))  : fill(tp_,  length(tcsv))
        scsv = haskey(wave_csv, :swl) ? Float64.(collect(wave_csv[:swl])) : fill(swl_, length(tcsv))
        length(tcsv) >= 2 || error("wave_csv: need at least 2 time samples")
        length(hcsv) == length(tcsv) ||
            error("wave_csv: time and hrms column lengths differ")
        all(isfinite, tcsv) && all(isfinite, hcsv) && all(isfinite, pcsv) && all(isfinite, scsv) ||
            error("wave_csv: all values must be finite numbers")
        # Tide and storm overlays are intentionally ignored when the
        # user provides their own time series — the CSV already encodes
        # whatever forcing they want.
        (; timebc=tcsv, hrmsbc=hcsv, tpbc=pcsv, swlbc=scsv, wangbc=zeros(length(tcsv)))
    else
        wf0 = constant_waves(; duration_days=dur_h_/24,
                              hrms=hrms_, tp=tp_, swl=swl_, dt_hours=0.5)
        timebc = wf0.timebc
        swlbc  = collect(wf0.swlbc); hrmsbc = collect(wf0.hrmsbc)
        if tide_on_ && tide_per_ > 0
            omega = 2π / (tide_per_ * 3600.0)
            swlbc .+= tide_amp_ .* sin.(omega .* timebc)
        end
        # Storm envelope: sinusoidal recurrence with period = storm_recur_d
        # (days). sin² peak shape; baseline is recovered exactly between
        # events. For a 12 h sim with the 30-day default the envelope
        # stays near baseline — shorten the recurrence for visible storm
        # impact in a short run.
        if storm_on_ && storm_recur_d > 0
            period_s = storm_recur_d * 86400.0
            for (i, t) in enumerate(timebc)
                rel = t - timebc[1]
                env = sin(pi * rel / period_s)^2
                hrmsbc[i] = hrms_ + (storm_hrms_ - hrms_) * env
            end
        end
        merge(wf0, (; swlbc=swlbc, hrmsbc=hrmsbc))
    end

    # ---- Multifraction sediment --------------------------------------
    gs_mm = [parse(Float64, strip(t)) for t in split(gs_str_, ",") if !isempty(strip(t))]
    gf    = [parse(Float64, strip(t)) for t in split(gf_str_, ",") if !isempty(strip(t))]
    # MultifractionConfig.validate requires sum(initial_fractions) ≈ 1.0 at
    # atol=1e-9, so always renormalize. The previous `if !isapprox(...; atol=1e-3)`
    # skipped normalization for inputs like "0.333,0.333,0.333" (sums to 0.999),
    # which then crashed validate() at the tight tolerance and silently failed
    # every multifraction run.
    s = sum(gf)
    s > 0 || error("grain_fractions must contain at least one positive value")
    gf = gf ./ s
    d50_m = sum(gs_mm .* gf) * 1e-3
    sed = make_sediment(; d50=d50_m, effb=effb_, efff=efff_, blp=blp_, slp=slp_, tanphi=tanphi_)
    mf  = MultifractionConfig(grain_sizes=gs_mm .* 1e-3, initial_fractions=gf)
    nf  = length(gs_mm)
    has_mixing = nf > 1

    # ---- Thermal ------------------------------------------------------
    thermal_kwargs = NamedTuple()
    if thermal_on_
        tcfg = ThermalConfig(T_init=T_init_)
        ntimes = length(wf.timebc)
        thermal_kwargs = (; thermal=tcfg,
                            T_air   = fill(T_air_,   ntimes),
                            T_water = fill(T_water_, ntimes))
    end

    # ---- OptionFlags + nbs kwargs ------------------------------------
    isedav = (hardbottom_z_vec !== nothing || thermal_on_) ? 1 : 0
    iveg   = veg_input         !== nothing ? 1 : 0
    iperm  = porous_z_vec      !== nothing ? 1 : 0
    opts_kwargs = (; options=OptionFlags(isedav=isedav, iveg=iveg, iperm=iperm))
    nbs_kwargs = NamedTuple()
    veg_input        !== nothing && (nbs_kwargs = merge(nbs_kwargs, (; vegetation=veg_input)))
    porous_z_vec     !== nothing && (nbs_kwargs = merge(nbs_kwargs, (; porous_z=porous_z_vec,
                                                                       porosity=nbs_por_, stone_diameter=nbs_stone_)))
    hardbottom_z_vec !== nothing && (nbs_kwargs = merge(nbs_kwargs, (; hardbottom_z=hardbottom_z_vec)))

    cfg = build_config(; bathymetry_x=x, bathymetry_z=z,
                         multifraction=mf, sediment=sed, dx=dx_, wf...,
                         opts_kwargs..., nbs_kwargs..., thermal_kwargs...)

    ncname = "$(case_).nc"
    ncpath = joinpath(workdir, ncname)
    total_s  = wf.timebc[end] - wf.timebc[1]
    interval = max(60.0, total_s / 24)

    t0 = Base.time()
    run_simulation!(cfg; outdir=workdir, outfile=ncname,
                    output_interval_s=interval,
                    progress_callback=nothing)
    dt = round(Base.time() - t0; digits=1)

    # ---- Pack arrays for the browser ---------------------------------
    out = Dict{Symbol,Any}(:status      => "ok",
                            :workdir     => workdir,
                            :ncpath      => ncpath,
                            :elapsed_s   => dt,
                            :has_mixing  => has_mixing)

    NCDataset(ncpath, "r") do ds
        xs   = Float64.(Array(ds["x"]))
        zb   = Array(ds["zb"])
        nx, nt = size(zb)
        times_s = haskey(ds, "time") ? Array(ds["time"]) : Float64[0.0]
        time_h  = if eltype(times_s) <: Dates.AbstractDateTime
            t0_ = times_s[1]
            [Dates.value(t - t0_)/3.6e6 for t in times_s]
        else
            Float64.(times_s) ./ 3600.0
        end

        idxs = nt == 1 ? [1] : unique(round.(Int, range(1, nt; length=min(6, nt))))
        out[:plot_x]           = xs
        out[:plot_z_snapshots] = [Float64.(@view zb[:, ti]) for ti in idxs]
        out[:plot_times_h]     = Float64.(time_h[idxs])
        out[:plot_swl]         = swl_
        out[:plot_dz_rows]     = nt > 1 ? [Float64.(zb[:, t] .- zb[:, 1]) for t in 1:nt] : Vector{Float64}[]
        out[:plot_times_full]  = Float64.(time_h)

        dz = Float64.(@view zb[:, nt]) .- Float64.(@view zb[:, 1])
        v_tot = 0.0
        for k in 1:length(xs)-1
            v_tot += 0.5 * (dz[k] + dz[k+1]) * (xs[k+1] - xs[k])
        end
        out[:volume_text] = "Net Δvol: $(round(v_tot; digits=2)) m³/m  ·  " *
                            "max ero: $(round(minimum(dz); digits=3)) m  ·  " *
                            "max dep: $(round(maximum(dz); digits=3)) m"

        # ---- Hydrodynamic envelope: max/min Hrms and water surface ---
        # Surface elevation η(x,t) = zb(x,t) + h(x,t). The max envelope
        # is the highest the water reached at each x across the run; the
        # min envelope is the lowest. Add ±Hrms_max around the max-WL
        # line to suggest the wave-crest envelope.
        if haskey(ds, "hrms") && haskey(ds, "h")
            hrms_var = Float64.(Array(ds["hrms"]))   # (nx, nt)
            h_var    = Float64.(Array(ds["h"]))      # (nx, nt) water depth
            wl_var   = Float64.(zb) .+ max.(h_var, 0.0)
            out[:hydro_x]       = xs
            out[:hydro_zb]      = Float64.(@view zb[:, 1])
            out[:hydro_wl_max]  = vec(maximum(wl_var;   dims=2))
            out[:hydro_wl_min]  = vec(minimum(wl_var;   dims=2))
            out[:hydro_hrms_max] = vec(maximum(hrms_var; dims=2))
            out[:hydro_hrms_min] = vec(minimum(hrms_var; dims=2))
        end

        if has_mixing && haskey(ds, "bed_mass") && haskey(ds, "fraction")
            bm   = Array(ds["bed_mass"])     # (nx, nlayers, nf, nt)
            gs_m = Float64.(Array(ds["fraction"]))
            nf_  = length(gs_m)
            bm_end = bm[:, 1, :, nt]
            tot    = sum(bm_end; dims=2); tot[tot .<= 0] .= 1.0
            comp   = bm_end ./ tot
            out[:mix_x]                   = xs
            out[:mix_grain_sizes_mm]      = gs_m .* 1e3
            out[:mix_surface_composition] = [Float64.(comp[:, k]) for k in 1:nf_]

            if nt > 1
                bm_surf = bm[:, 1, :, :]
                tot_t   = sum(bm_surf; dims=2); tot_t[tot_t .<= 0] .= 1.0
                d50_surf = dropdims(sum(bm_surf .* reshape(gs_m .* 1e3, 1, :, 1); dims=2) ./ tot_t; dims=2)
                out[:mix_d50_rows] = [Float64.(d50_surf[:, t]) for t in 1:nt]
                out[:mix_times_h]  = Float64.(time_h)
            end
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Ensemble: vary one knob at a time, run baseline + N members, return
# the final bed profile + an erosion-volume metric for each.
# The client builds the member-override list (see public/index.html), so
# the server is a thin runner — it just merges each member's overrides
# into the base params and dispatches run_simulation_json.
# ---------------------------------------------------------------------------

# Erosion volume above a threshold elevation. Counts only cells where
# the bed dropped (z_final < z_initial) and only over the cross-shore
# region where the initial bed was above z_threshold. Returns m³/m.
function _erosion_above_threshold(x::AbstractVector, z_init::AbstractVector,
                                   z_final::AbstractVector, z_thresh::Float64)
    v = 0.0
    @inbounds for i in 1:length(x)-1
        z0a = (z_init[i]   + z_init[i+1]  ) / 2
        z0a >= z_thresh || continue
        ero_l = max(z_init[i]   - z_final[i],   0.0)
        ero_r = max(z_init[i+1] - z_final[i+1], 0.0)
        v += 0.5 * (ero_l + ero_r) * (x[i+1] - x[i])
    end
    return v
end

# Walk a JSON-bound value and replace any non-finite Float (NaN, +Inf, -Inf)
# with `nothing` so JSON3 doesn't reject the response. JSON3 follows the
# strict JSON spec, which has no NaN/Inf literal — without this, a single
# divergent sim or an empty-erosion baseline would crash the whole response.
function _sanitize_json(v)
    if v isa AbstractFloat
        return isfinite(v) ? v : nothing
    elseif v isa AbstractArray
        return [_sanitize_json(x) for x in v]
    elseif v isa AbstractDict
        return Dict(k => _sanitize_json(val) for (k, val) in pairs(v))
    else
        return v
    end
end

const ENSEMBLE_MAX_MEMBERS = 10

function run_ensemble_json(spec)
    base = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in pairs(spec[:base_params]))
    members_in = collect(spec[:members])
    isempty(members_in) && error("ensemble: need at least 1 member")
    length(members_in) > ENSEMBLE_MAX_MEMBERS &&
        error("ensemble: $(length(members_in)) members exceeds cap of $(ENSEMBLE_MAX_MEMBERS)")

    # Effectiveness threshold defaults to SWL (water-line erosion).
    z_thresh = Float64(get(spec, :effectiveness_threshold_z,
                            Float64(get(base, :swl, 0.0))))

    base_case = String(get(base, :case_name, "webrun"))
    sweep_param = String(get(spec, :sweep_param, ""))

    # Job 0 = baseline (base params unchanged). Jobs 1..N = sweep members
    # with the client-supplied override dict merged on top.
    jobs = Tuple{String, Dict{Symbol,Any}}[]
    push!(jobs, ("baseline",
                 merge(base, Dict{Symbol,Any}(:case_name => base_case * "_baseline"))))
    for (i, m) in enumerate(members_in)
        lbl = String(get(m, :label, "member_$(i)"))
        ovr_raw = get(m, :overrides, nothing)
        ovr = Dict{Symbol,Any}()
        if ovr_raw !== nothing
            for (k, v) in pairs(ovr_raw)
                ovr[Symbol(k)] = v
            end
        end
        p = merge(base, ovr)
        p[:case_name] = base_case * "_m$(i)"
        push!(jobs, (lbl, p))
    end

    # Run sequentially. CSHORE is mostly single-threaded inside a sim, so
    # parallel @spawn could win 2-3x on the 4-CPU machine — but it also
    # multiplies the peak memory footprint, and on the 8 GB tier that's
    # tight. Keep it sequential for v1; revisit if wall time bites.
    out_members = Vector{Dict{Symbol,Any}}(undef, length(jobs))
    for (i, (lbl, p)) in enumerate(jobs)
        @info "ensemble member $(i)/$(length(jobs)): $(lbl)"
        try
            r = run_simulation_json(p)
            snaps = r[:plot_z_snapshots]
            z_init  = snaps[1]
            z_final = snaps[end]
            ero = _erosion_above_threshold(r[:plot_x], z_init, z_final, z_thresh)
            out_members[i] = Dict{Symbol,Any}(
                :label           => lbl,
                :status          => "ok",
                :x               => r[:plot_x],
                :z_initial       => z_init,
                :z_final         => z_final,
                :erosion_above_z => ero,
                :elapsed_s       => r[:elapsed_s],
            )
        catch e
            @error "ensemble member failed" label=lbl exception=(e, catch_backtrace())
            out_members[i] = Dict{Symbol,Any}(
                :label   => lbl,
                :status  => "error",
                :message => first(sprint(showerror, e), 400),
            )
        end
    end

    # Reduction % vs baseline. Use `nothing` instead of NaN where the
    # comparison is undefined (baseline failed, or eroded nothing above
    # the threshold) — JSON3 rejects NaN.
    base_ok = out_members[1][:status] == "ok"
    base_ero = base_ok ? Float64(out_members[1][:erosion_above_z]) : 0.0
    for m in out_members
        if m[:status] != "ok" || !base_ok || base_ero <= 0
            m[:reduction_pct] = nothing
        else
            m[:reduction_pct] = (1 - Float64(m[:erosion_above_z]) / base_ero) * 100
        end
    end

    return Dict{Symbol,Any}(
        :status      => "ok",
        :sweep_param => sweep_param,
        :threshold_z => z_thresh,
        :members     => out_members,
    )
end

# ---------------------------------------------------------------------------
# HTTP request router.
# ---------------------------------------------------------------------------

# CORS headers — open to all origins so a GitHub-Pages-hosted UI (or any
# other static-front-end deployment) can call this backend cross-origin.
# Adjust `Access-Control-Allow-Origin` to a specific domain if you want
# to tighten down later.
const _CORS_HEADERS = [
    "Access-Control-Allow-Origin"  => "*",
    "Access-Control-Allow-Methods" => "GET, POST, OPTIONS",
    "Access-Control-Allow-Headers" => "Content-Type",
    "Access-Control-Max-Age"       => "86400",
]

function _with_cors(headers::Vector)
    vcat(headers, _CORS_HEADERS)
end

function _json(status::Int, body)
    HTTP.Response(status,
        _with_cors(["Content-Type" => "application/json"]),
        JSON3.write(body))
end

function handler(req::HTTP.Request)
    target = HTTP.URIs.URI(req.target).path
    # CORS preflight: respond 204 to any OPTIONS request with the CORS
    # headers. Browsers send these before cross-origin POSTs.
    if req.method == "OPTIONS"
        return HTTP.Response(204, _CORS_HEADERS)
    end
    if req.method == "GET" && target == "/"
        return HTTP.Response(200,
            _with_cors(["Content-Type" => "text/html; charset=utf-8"]),
            load_index_html())
    elseif req.method == "GET" && target == "/healthz"
        return HTTP.Response(200, _CORS_HEADERS, "ok")
    elseif req.method == "POST" && target == "/run"
        try
            params = JSON3.read(String(req.body))
            result = run_simulation_json(params)
            return _json(200, result)
        catch e
            msg = sprint(showerror, e)
            @error "run failed" exception=(e, catch_backtrace())
            return _json(500, Dict(:status=>"error",
                                    :message=>first(msg, 800)))
        end
    elseif req.method == "POST" && target == "/run_ensemble"
        try
            spec = JSON3.read(String(req.body))
            result = _sanitize_json(run_ensemble_json(spec))
            return _json(200, result)
        catch e
            msg = sprint(showerror, e)
            @error "ensemble run failed" exception=(e, catch_backtrace())
            return _json(500, Dict(:status=>"error",
                                    :message=>first(msg, 800)))
        end
    end
    return HTTP.Response(404, _CORS_HEADERS, "Not found")
end

const PORT = parse(Int, get(ENV, "PORT", "8000"))
const HOST = get(ENV, "HOST", "0.0.0.0")

# ---------------------------------------------------------------------------
# Server-startup warmup. Without precompile baked into the image (we're on
# Dockerfile.lazy-precompile), the *first* /run after deploy pays for JIT
# specialization of the entire CSHORE call tree on top of the simulation
# itself — easily several minutes. That exceeds Fly's edge-proxy idle
# timeout, so the browser sees the request hang forever and no plots
# render. We pre-pay the JIT cost here by running a tiny synthetic case
# before HTTP.serve binds. After this, every real /run is just the
# simulation work (~10–60 s) and returns promptly.
#
# Set WEB_GUI_SKIP_WARMUP=1 to skip (useful for fast local dev iteration).
# ---------------------------------------------------------------------------
function _warmup()
    if get(ENV, "WEB_GUI_SKIP_WARMUP", "0") == "1"
        @info "web_gui: WEB_GUI_SKIP_WARMUP=1 — skipping startup warmup"
        return
    end
    @info "web_gui: warming up simulation paths (one-time JIT, ~2–5 min)..."
    t0 = Base.time()
    base = Dict{Symbol,Any}(
        :profile          => "planar_beach",
        :depth_m          => 2.0,
        :slope_rr         => 0.05,
        :backshore_elev_m => 1.0,
        :dune_m           => 1.0,
        :dx_m             => 1.0,
        :hrms             => 0.5,
        :tp               => 6.0,
        :swl              => 0.0,
        :duration_h       => 0.5,
    )
    # Two cases so both the single-fraction and the multifraction
    # (per-fraction transport / mixing / hiding-exposure) code branches
    # get JIT-specialized before the user's first click.
    cases = [
        ("warmup_1frac", "0.30",        "1.0"),
        ("warmup_mfx",  "0.20,0.50",   "0.5,0.5"),
    ]
    for (name, gs, gf) in cases
        try
            run_simulation_json(merge(base, Dict{Symbol,Any}(
                :case_name       => name,
                :grain_sizes_mm  => gs,
                :grain_fractions => gf,
            )))
            @info "web_gui: warmup case $(name) ok"
        catch e
            @warn("web_gui: warmup case $(name) failed — first real /run on " *
                  "that branch will pay full JIT cost",
                  exception=(e, catch_backtrace()))
        end
    end
    @info "web_gui: warmup complete in $(round(Base.time() - t0; digits=1))s"
end
_warmup()

@info "web_gui (lean): starting HTTP server on $HOST:$PORT (threads=$(Threads.nthreads()))"
HTTP.serve(handler, HOST, PORT)
