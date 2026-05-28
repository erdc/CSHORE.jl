#!/usr/bin/env julia
#=============================================================================
revetment_demo.jl — Gravel revetment at the toe of a fine-grained bluff.

10 grain-size classes (0.5–6 mm), size-adaptive transport (CSHORE for
sand, MPM for gravel), seasonal wave + surge forcing that reaches the
bluff several times per year, 1-year simulation.

Produces:
  - NetCDF output with 12-hourly snapshots
  - Comprehensive 8-panel MP4 animation (nearshore zoom, d50-colored profile)
  - Static diagnostic PNGs
=============================================================================#

import Pkg
const CSHORE_PATH = get(ENV, "CSHORE_PATH", joinpath(@__DIR__, ".."))
Pkg.activate(@__DIR__)
try
    using CSHORE
    import CairoMakie, NCDatasets, Dates
catch
    Pkg.develop(path=CSHORE_PATH)
    Pkg.add(["CairoMakie", "NCDatasets", "Dates", "Printf"])
    Pkg.instantiate()
end
using CSHORE
import CairoMakie, NCDatasets, Dates
using Printf

OUTDIR = joinpath(@__DIR__, "movie_output", "revetment_demo")
mkpath(OUTDIR)

# ============================================================================
# 1. PROFILE: gentle shelf → sand beach → wide gravel revetment → silt bluff
# ============================================================================
dx = 2.0

x_pts = [
    0.0,     # offshore
    300.0,   # mid-shelf break
    450.0,   # beach toe (MSL) — short, steep nearshore
    470.0,   # top of beach berm (+0.8 m)
    470.0,   # base of revetment
    510.0,   # top of revetment (+2.8 m) — 40 m wide
    512.0,   # bluff toe behind revetment
    530.0,   # bluff crest (+6 m)
    580.0,   # upland
]
z_pts = [
    -8.0,    # offshore depth
    -3.0,    # steep shoreface
    0.0,     # MSL shoreline
    0.8,     # berm crest
    0.8,     # revetment base
    2.8,     # revetment top (armored, 40 m wide)
    2.8,     # bluff toe
    6.0,     # bluff crest
    6.0,     # flat upland
]

function linterp(xp, zp, xi)
    xi ≤ xp[1] && return zp[1]
    xi ≥ xp[end] && return zp[end]
    for k in 1:(length(xp)-1)
        if xp[k] ≤ xi ≤ xp[k+1]
            t = (xi - xp[k]) / (xp[k+1] - xp[k])
            return zp[k] + t * (zp[k+1] - zp[k])
        end
    end
    return zp[end]
end

x = collect(0.0:dx:x_pts[end])
z = [linterp(x_pts, z_pts, xi) for xi in x]
np = length(x)
println("Profile: $np nodes, x ∈ [$(x[1]), $(x[end])] m")

# ============================================================================
# 2. HARDBOTTOM: revetment is non-erodible (x = 800–840 m, 40 m wide)
# ============================================================================
# Hardbottom: thin veneer over bedrock everywhere, EXCEPT the revetment zone
# where the gravel armor layer is up to 1.5 m thick, tapering to zero at the
# waterline (x=450, MSL). This gives the revetment a wedge of mobile gravel
# sitting on top of bedrock that waves can rework but not erode through.
veneer_sand = 0.15       # thin sand veneer on bedrock outside revetment
revet_max_thick = 1.5    # max gravel thickness at landward end of revetment
revet_x_start = 450.0    # waterline (taper starts, thickness = 0)
revet_x_end   = 510.0    # landward end (full thickness)

zh = copy(z)
for i in 1:np
    if x[i] ≥ revet_x_start && x[i] ≤ revet_x_end
        # Linear taper: 0 at waterline → revet_max_thick at landward end
        frac = (x[i] - revet_x_start) / (revet_x_end - revet_x_start)
        thickness = revet_max_thick * frac
        zh[i] = z[i] - thickness
    else
        zh[i] = z[i] - veneer_sand
    end
end
println("Hardbottom: revetment gravel wedge 0–$(revet_max_thick)m (x=$(revet_x_start)–$(revet_x_end)m), $(veneer_sand)m sand elsewhere")

# ============================================================================
# 3. SEDIMENT: 10 grain classes from 0.5 mm to 6 mm
# ============================================================================
grain_sizes = [
    0.50e-3,   # 1: medium sand
    0.71e-3,   # 2: coarse sand
    1.00e-3,   # 3: very coarse sand
    1.41e-3,   # 4: granule
    2.00e-3,   # 5: fine gravel (MPM threshold)
    2.83e-3,   # 6: fine gravel
    3.36e-3,   # 7: medium gravel
    4.00e-3,   # 8: medium gravel
    5.00e-3,   # 9: coarse gravel
    6.00e-3,   # 10: coarse gravel
]
nf = length(grain_sizes)
println("Sediment: $nf fractions — $(grain_sizes .* 1e3) mm")

# Spatial composition: 4 zones
fracs = zeros(np, nf)
for i in 1:np
    if x[i] < 350.0
        # Offshore shelf: medium-coarse sand dominated
        fracs[i, :] = [0.30, 0.25, 0.20, 0.10, 0.05, 0.04, 0.03, 0.02, 0.005, 0.005]
    elseif x[i] < 470.0
        # Beach/shoreface: coarser sand + some gravel
        fracs[i, :] = [0.15, 0.20, 0.20, 0.15, 0.10, 0.08, 0.05, 0.04, 0.02, 0.01]
    elseif x[i] ≤ 510.0
        # Revetment: coarse gravel dominated
        fracs[i, :] = [0.01, 0.01, 0.02, 0.03, 0.05, 0.10, 0.15, 0.23, 0.25, 0.15]
    else
        # Bluff: finer material (sandy silt/fine sand)
        fracs[i, :] = [0.35, 0.25, 0.15, 0.10, 0.05, 0.04, 0.03, 0.02, 0.005, 0.005]
    end
    fracs[i, :] ./= sum(fracs[i, :])
end

mf = MultifractionConfig(
    grain_sizes             = grain_sizes,
    initial_fractions       = fracs[1, :],
    initial_fractions_spatial = fracs,
    nlayers                 = 5,
    layer_thickness         = 0.25,
    porosity                = 0.4,
    transport_formula       = :size_adaptive,
    n_face_flux_smooth      = 3,
    n_pickup_smooth         = 15,
    n_composition_smooth    = 3,
)

# ============================================================================
# 4. WAVE FORCING: seasonal with storms + surge
# ============================================================================
duration_days = 730
dt_hours = 1.0
ntimes = round(Int, duration_days * 24 / dt_hours) + 1
timebc = collect(range(0.0, duration_days * 86400.0; length=ntimes))
year_s = 365.25 * 86400.0

hrmsbc = [0.5 + 2.5 * (0.5 + 0.5 * cos(2π * t / year_s))^2 for t in timebc]
tpbc   = [8.0 + 4.0 * (0.5 + 0.5 * cos(2π * t / year_s)) for t in timebc]

swlbc = [0.5 * sin(2π * t / (12.42 * 3600)) +
         1.2 * max(0, cos(2π * t / year_s))^3 +
         0.8 * max(0, sin(2π * t / (14 * 86400)))^6 *
               max(0, cos(2π * t / year_s))^2
         for t in timebc]
for i in 1:ntimes
    day = timebc[i] / 86400
    for storm_day in [30.0, 75.0, 105.0, 320.0]
        dt_storm = abs(day - storm_day)
        if dt_storm < 2.5
            spike = 2.2 * exp(-dt_storm^2 / 0.4)
            swlbc[i] += spike
            hrmsbc[i] = max(hrmsbc[i], 3.0 * exp(-dt_storm^2 / 0.6))
        end
    end
end

println("Waves: Hrms ∈ [$(round(minimum(hrmsbc),digits=2)), $(round(maximum(hrmsbc),digits=2))] m")
println("SWL:   ∈ [$(round(minimum(swlbc),digits=2)), $(round(maximum(swlbc),digits=2))] m")

# ============================================================================
# 5. BUILD & RUN
# ============================================================================
opts = OptionFlags(iprofl=1, isedav=1)
output_interval_s = 12.0 * 3600.0
nc_path = joinpath(OUTDIR, "revetment_1yr.nc")

cfg = build_config(
    dx=dx, bathymetry_x=x, bathymetry_z=z, friction=0.003,
    timebc=timebc, tpbc=tpbc, hrmsbc=hrmsbc, swlbc=swlbc,
    options=opts, sediment=make_sediment(d50=1.0e-3),
    multifraction=mf, hardbottom_z=zh,
)

println("\nRunning 1-year simulation...")
@time state = run_simulation!(cfg; outfile=nc_path, output_interval_s=output_interval_s)
jmax = state.jmax[1]
println("  jmax=$jmax, jr=$(state.jr), NetCDF: $nc_path")

# ============================================================================
# 6. STATIC PLOTS
# ============================================================================
println("\nGenerating static plots...")
save_figure(plot_profile(x, z; title="Initial: Beach + Gravel Revetment + Bluff", swl=0.0),
    joinpath(OUTDIR, "01_initial_profile.png"))
save_figure(plot_profile_evolution(nc_path; title="Profile Evolution (1 year)"),
    joinpath(OUTDIR, "02_profile_evolution.png"))
save_figure(plot_hovmoller(nc_path; var=:zb, title="Bed Elevation"),
    joinpath(OUTDIR, "03_pcolor_zb.png"))
save_figure(plot_hovmoller(nc_path; var=:d50_surface, title="Surface d50"),
    joinpath(OUTDIR, "04_pcolor_d50.png"))
save_figure(plot_mass_balance(nc_path),
    joinpath(OUTDIR, "05_mass_balance.png"))

# Composition comparison
frac_colors = CairoMakie.cgrad(:turbo, nf; categorical=true)
fig_comp = CairoMakie.Figure(size=(1100, 700))
for (row, (ttl, get_frac)) in enumerate([
    ("Initial composition", j -> fracs[j, :]),
    ("Final surface composition", j -> begin
        m = [max(0.0, state.bed_mass[j, 1, k]) for k in 1:nf]
        s = sum(m); s > 0 ? m ./ s : mf.initial_fractions
    end),
])
    ax = CairoMakie.Axis(fig_comp[row, 1], ylabel="Fraction", title=ttl,
        xlabel= row==2 ? "x (m)" : "")
    nn_plot = min(jmax, np)
    bottom = zeros(nn_plot)
    for k in 1:nf
        fk = [clamp(get_frac(j)[k], 0, 1) for j in 1:nn_plot]
        CairoMakie.band!(ax, x[1:nn_plot], bottom, bottom .+ fk, color=(frac_colors[k], 0.8))
        bottom .+= fk
    end
    CairoMakie.ylims!(ax, 0, 1.05)
end
CairoMakie.Label(fig_comp[0, 1],
    "Grain sizes: $(join([@sprintf("%.2g",d*1e3) for d in grain_sizes], ", ")) mm",
    fontsize=11, font=:bold)
save_figure(fig_comp, joinpath(OUTDIR, "06_composition.png"))

# ============================================================================
# 7. COMPREHENSIVE 8-PANEL ANIMATION — zoomed nearshore, d50-colored profile
# ============================================================================
println("\nGenerating animation (nearshore zoom, d50-colored profile)...")

NCDatasets.NCDataset(nc_path, "r") do ds
    x_nc   = Array(ds["x"])
    zb_nc  = Array(ds["zb"])
    hrms_nc = Array(ds["hrms"])
    h_nc   = Array(ds["h"])
    d50_nc = Array(ds["d50_surface"])
    bm_nc  = Array(ds["bed_mass"])
    time_raw = Array(ds["time"])

    time_s = if eltype(time_raw) <: Dates.AbstractDateTime
        t0 = time_raw[1]; Float64[Dates.value(t - t0) / 1000.0 for t in time_raw]
    else Float64.(time_raw) end
    # Pick a display unit for the time axis based on run length.
    total = time_s[end]
    divisor = total >= 2 * 365.25 * 86400 ? 365.25 * 86400 :
              total >= 86400              ? 86400.0         :
              total >= 3600               ? 3600.0          : 1.0
    t_unit  = divisor == 365.25 * 86400 ? "yr" :
              divisor == 86400.0        ? "d"  :
              divisor == 3600.0         ? "h"  : "s"
    time_disp = time_s ./ divisor

    nx_nc = size(zb_nc, 1)
    nt_nc = size(zb_nc, 2)

    # Nearshore zoom limits (focus on shoreface → bluff)
    x_zoom_lo = 300.0
    x_zoom_hi = x_nc[end]
    z_zoom_lo = -3.0
    z_zoom_hi = 7.5

    # Forcing interpolated to output times
    hrms_f = [hrmsbc[clamp(round(Int, time_s[ti]/(timebc[end]/(ntimes-1)))+1, 1, ntimes)] for ti in 1:nt_nc]
    swl_f  = [swlbc[clamp(round(Int, time_s[ti]/(timebc[end]/(ntimes-1)))+1, 1, ntimes)] for ti in 1:nt_nc]

    # d50 colormap range
    d50_min_mm = minimum(grain_sizes) * 1e3
    d50_max_mm = maximum(grain_sizes) * 1e3
    gs_mm = grain_sizes .* 1e3

    # ---- Figure layout: 4 rows × 2 cols ----
    fig = CairoMakie.Figure(size=(1500, 1100))
    title_obs = CairoMakie.Observable("")
    CairoMakie.Label(fig[0, 1:2], title_obs, fontsize=14, font=:bold)

    # R1C1: Profile colored by d50
    ax_prof = CairoMakie.Axis(fig[1, 1], ylabel="z (m)", title="Profile (colored by d50)")
    CairoMakie.lines!(ax_prof, x_nc, zb_nc[:, 1], color=(:black, 0.25), linestyle=:dash, linewidth=0.7)
    CairoMakie.band!(ax_prof, [470.0, 510.0], [z_zoom_lo, z_zoom_lo], [2.8, 2.8],
        color=(:gray, 0.2))  # revetment
    # d50-colored profile line — use linesegments for per-segment color
    # Initialize with first frame's data so Observable lengths match
    _init_pts = CairoMakie.Point2f[]
    _init_cols = Float32[]
    for j in 1:(nx_nc-1)
        push!(_init_pts, CairoMakie.Point2f(x_nc[j], zb_nc[j, 1]))
        push!(_init_pts, CairoMakie.Point2f(x_nc[j+1], zb_nc[j+1, 1]))
        push!(_init_cols, Float32(d50_nc[j, 1] * 1e3))
        push!(_init_cols, Float32(d50_nc[j+1, 1] * 1e3))
    end
    prof_segs = CairoMakie.Observable(_init_pts)
    prof_cols = CairoMakie.Observable(_init_cols)
    seg_plot = CairoMakie.linesegments!(ax_prof, prof_segs; color=prof_cols,
        colormap=CairoMakie.Reverse(:RdYlBu), colorrange=(d50_min_mm, d50_max_mm),
        linewidth=3)
    CairoMakie.Colorbar(fig[1, 3], seg_plot; label="d50 (mm)", width=12)
    swl_line = CairoMakie.hlines!(ax_prof, [0.0], color=(:dodgerblue, 0.5), linestyle=:dot)
    CairoMakie.xlims!(ax_prof, x_zoom_lo, x_zoom_hi)
    CairoMakie.ylims!(ax_prof, z_zoom_lo, z_zoom_hi)

    # R1C2: Forcing time series
    ax_force = CairoMakie.Axis(fig[1, 2], ylabel="m", title="Wave + Surge Forcing")
    CairoMakie.lines!(ax_force, time_disp, hrms_f, color=:teal, linewidth=0.7, label="Hrms")
    CairoMakie.lines!(ax_force, time_disp, swl_f, color=:coral, linewidth=0.7, label="SWL")
    cursor_f = CairoMakie.vlines!(ax_force, [0.0], color=:black, linewidth=1)
    CairoMakie.axislegend(ax_force, position=:rt, framevisible=false, labelsize=8)

    # R2C1: Hrms (zoomed)
    ax_hrms = CairoMakie.Axis(fig[2, 1], ylabel="Hrms (m)", title="Wave Height")
    line_hrms = CairoMakie.lines!(ax_hrms, x_nc, hrms_nc[:, 1], color=:teal, linewidth=1.5)
    CairoMakie.xlims!(ax_hrms, x_zoom_lo, x_zoom_hi)
    CairoMakie.ylims!(ax_hrms, 0, maximum(hrms_nc)*1.1+0.01)

    # R2C2: d50 profile
    ax_d50 = CairoMakie.Axis(fig[2, 2], ylabel="d50 (mm)", title="Surface d50")
    for gs in gs_mm
        CairoMakie.hlines!(ax_d50, [gs], color=(:gray, 0.12), linewidth=0.4)
    end
    line_d50 = CairoMakie.lines!(ax_d50, x_nc, d50_nc[:, 1] .* 1e3, color=:darkorange, linewidth=1.5)
    CairoMakie.xlims!(ax_d50, x_zoom_lo, x_zoom_hi)
    CairoMakie.ylims!(ax_d50, 0, d50_max_mm * 1.15)

    # R3: Stacked composition (full width, zoomed) with legend
    ax_comp = CairoMakie.Axis(fig[3, 1:2], ylabel="Fraction",
        title="Surface Composition ($(nf) classes)")
    CairoMakie.xlims!(ax_comp, x_zoom_lo, x_zoom_hi)
    CairoMakie.ylims!(ax_comp, 0, 1.05)
    comp_bands = []
    # Static legend entries (dummy polys for the legend)
    for k in 1:nf
        CairoMakie.poly!(ax_comp, CairoMakie.Point2f[(0,0),(0,0),(0,0)],
            color=(frac_colors[k], 0.8),
            label=@sprintf("%.2g mm", grain_sizes[k]*1e3))
    end
    CairoMakie.axislegend(ax_comp; position=:rt, framevisible=true,
        labelsize=7, nbanks=2, patchsize=(8,8))

    # R4C1: Bluff detail zoom
    ax_bluff = CairoMakie.Axis(fig[4, 1], xlabel="x (m)", ylabel="z (m)",
        title="Revetment + Bluff Detail")
    CairoMakie.lines!(ax_bluff, x_nc, zb_nc[:, 1], color=(:black, 0.25), linestyle=:dash, linewidth=0.7)
    CairoMakie.band!(ax_bluff, [470.0, 510.0], [0.0, 0.0], [2.8, 2.8], color=(:gray, 0.25))
    line_zb_zoom = CairoMakie.lines!(ax_bluff, x_nc, zb_nc[:, 1], color=:steelblue, linewidth=2)
    CairoMakie.xlims!(ax_bluff, 400.0, 560.0)
    CairoMakie.ylims!(ax_bluff, -2, 7)

    # R4C2: Per-class transport
    ax_qt = CairoMakie.Axis(fig[4, 2], xlabel="x (m)", ylabel="Q (m²/s)",
        title="Per-fraction Transport")
    CairoMakie.xlims!(ax_qt, x_zoom_lo, x_zoom_hi)
    qt_lines = []

    # ---- Animate ----
    mp4_path = joinpath(OUTDIR, "revetment_1yr.mp4")
    frame_skip = max(1, nt_nc ÷ 500)

    CairoMakie.record(fig, mp4_path, 1:frame_skip:nt_nc; framerate=20) do ti
        t_val = time_disp[ti]
        title_obs[] = @sprintf("Gravel Revetment — t = %.1f %s / %.1f %s  |  %d grain classes (%.1f–%.1f mm)  |  size-adaptive transport",
            t_val, t_unit, time_disp[end], t_unit, nf, d50_min_mm, d50_max_mm)

        # Profile colored by d50: build line segments
        zb_t = zb_nc[:, ti]
        d50_t = d50_nc[:, ti] .* 1e3
        pts = CairoMakie.Point2f[]
        cols = Float32[]
        for j in 1:(nx_nc-1)
            push!(pts, CairoMakie.Point2f(x_nc[j], zb_t[j]))
            push!(pts, CairoMakie.Point2f(x_nc[j+1], zb_t[j+1]))
            push!(cols, Float32(d50_t[j]))
            push!(cols, Float32(d50_t[j+1]))
        end
        prof_segs[] = pts
        prof_cols[] = cols

        # SWL
        delete!(ax_prof, swl_line)
        swl_line = CairoMakie.hlines!(ax_prof, [swl_f[ti]], color=(:dodgerblue, 0.5), linestyle=:dot)

        # Forcing cursor
        delete!(ax_force, cursor_f)
        cursor_f = CairoMakie.vlines!(ax_force, [t_val], color=:black, linewidth=1)

        # Hrms
        line_hrms[2] = hrms_nc[:, ti]

        # d50
        line_d50[2] = d50_t

        # Composition
        for b in comp_bands; delete!(ax_comp, b); end
        empty!(comp_bands)
        bottom_arr = zeros(nx_nc)
        for k in 1:nf
            surf_k = [max(0.0, bm_nc[j, 1, k, ti]) for j in 1:nx_nc]
            total_s = [sum(max(0.0, bm_nc[j, 1, kk, ti]) for kk in 1:nf) for j in 1:nx_nc]
            fk = clamp.(surf_k ./ max.(total_s, 1e-12), 0.0, 1.0)
            b = CairoMakie.band!(ax_comp, x_nc, bottom_arr, bottom_arr .+ fk,
                color=(frac_colors[k], 0.8))
            push!(comp_bands, b)
            bottom_arr .+= fk
        end

        # Bluff zoom
        line_zb_zoom[2] = zb_t

        # Per-fraction transport
        for ln in qt_lines; delete!(ax_qt, ln); end
        empty!(qt_lines)
        qbx_t = ds["qbx"][:, :, ti]
        qsx_t = ds["qsx"][:, :, ti]
        for k in 1:nf
            qt_k = qbx_t[:, k] .+ qsx_t[:, k]
            ln = CairoMakie.lines!(ax_qt, x_nc, qt_k, color=(frac_colors[k], 0.8), linewidth=0.7)
            push!(qt_lines, ln)
        end
    end
    println("  Movie: $mp4_path ($(nt_nc ÷ frame_skip) frames)")
end

println("\nAll outputs in: $OUTDIR")
println("Done!")
