#=============================================================================
breakwater_demo.jl — Effect of an offshore breakwater on shoreline + bluff
                     erosion, single grain size.

Compares two configurations sharing identical wave / SWL forcing and
identical initial bathymetry except for the presence of a submerged
offshore breakwater modeled as a hardbottom block:

  Case A (no breakwater): clean shoreface → beach → bluff
  Case B (with breakwater): same profile + hardbottom block at x=180–220 m,
                            crest at z = -1.0 m (submerged)

Single-grain-size sand bed (d50 = 0.3 mm), 7-day storm event with one peak.

Output: examples/movie_output/breakwater_demo/breakwater_compare.png
=============================================================================#

import Pkg
const CSHORE_PATH = get(ENV, "CSHORE_PATH", joinpath(@__DIR__, ".."))
Pkg.activate(@__DIR__)
try
    using CSHORE
    import CairoMakie, NCDatasets
catch
    Pkg.develop(path=CSHORE_PATH)
    Pkg.add(["CairoMakie", "NCDatasets", "Printf"])
    Pkg.instantiate()
end
using CSHORE
import CairoMakie, NCDatasets
using Printf

const OUTDIR = joinpath(@__DIR__, "movie_output", "breakwater_demo")
mkpath(OUTDIR)

# ============================================================================
# 1. Shared geometry: shoreface → beach → bluff
# ============================================================================
const dx = 1.0
x = collect(0.0:dx:500.0)

# Profile (three-piece shoreface + beach + bluff):
#   x = 0      — z = -10.5 m  (deeper offshore boundary than the previous
#                              -8 m, gives the gentler slope room to work)
#   x = 0..200 — slope 1:80  (gentle offshore — eases wave shoaling onto
#                              the breakwater seaward face / apron)
#   x = 200..350 — slope 1:21 (mid-shoreface transition — steeper here so
#                              the inshore section stays close to 1:50)
#   x = 350..400 — slope 1:50 (inshore shoreface — same as the original
#                              full-domain slope, preserves beach-front
#                              wave climate near MSL)
#   x = 400..420 — narrow beach, +0.5 m at x=420
#   x = 420..440 — bluff face, +5 m at x=440
#   x = 440+    — plateau
function build_z(x)
    z = similar(x)
    for (i, xi) in pairs(x)
        z[i] = if xi < 200.0
            -10.5 + xi * (2.5 / 200.0)               # 1:80 slope, -10.5 → -8 m
        elseif xi < 350.0
            -8.0 + (xi - 200.0) * (7.0 / 150.0)      # 1:21 slope, -8 → -1 m
        elseif xi < 400.0
            -1.0 + (xi - 350.0) * (1.0 / 50.0)       # 1:50 slope, -1 → 0 m
        elseif xi < 420.0
            0.0 + (xi - 400.0) * (0.5 / 20.0)        # short beach, +0.5 m
        elseif xi < 440.0
            0.5 + (xi - 420.0) * (4.5 / 20.0)        # bluff face, +5 m
        else
            5.0                                      # plateau
        end
    end
    return z
end
z = build_z(x)

# ============================================================================
# 2. Hardbottom configurations
# ============================================================================
# Case A: only the bluff is hardbottom-protected (so the test focuses on
#         shoreline + lower-bluff erosion, not whole-cliff collapse).
# Case B: same as A + a submerged offshore breakwater modeled as a
#         hardbottom block at x ∈ [180, 220] m with crest at z = −1.0 m.

# Bluff hardbottom: face is erodible (zh well below z) but plateau is hard
# to keep simulation stable — a real cliff wouldn't fully collapse in 7 days.
zh_base = fill(-1e30, length(x))
for (i, xi) in pairs(x)
    if xi >= 440.0
        zh_base[i] = 4.5      # plateau hardbottom 0.5 m below crest
    end
end

# Breakwater addition for Case B — submerged offshore breakwater modeled as
# a hardbottom block. Geometry:
#   - 40 m wide crest at z = -1 m
#   - 10 m wide sloping flanks ramping down to the natural shoreface
#   - 20 m wide TOE-PROTECTION APRON on each side (rip-rap / geotextile)
#     extending the hardbottom along the natural seabed to suppress the
#     unrealistic toe scour we'd otherwise see at the abrupt edge of an
#     unprotected breakwater foundation. This mimics real engineering
#     practice (e.g. CIRIA, ASCE rubble-mound design guidance).
zh_break = copy(zh_base)
# Genuine offshore breakwater: in the surf-zone outer edge at ~3 m depth,
# crest emergent at +0.5 m so it overtops only at the storm-surge peak.
# Toe-protection apron (defined below in `break_z`) extends 20 m on each
# flank to suppress the unrealistic toe scour that an unprotected
# hardbottom edge produces in a single-grain CSHORE setup.
const BREAK_CREST = 0.5       # crest 0.5 m above MSL (emergent rubble-mound)
const CREST_X_LO = 220.0
const CREST_X_HI = 270.0
const FLANK_W    = 10.0
const APRON_W    = 20.0
const FLANK_X_LO = CREST_X_LO - FLANK_W
const FLANK_X_HI = CREST_X_HI + FLANK_W
const APRON_X_LO = FLANK_X_LO - APRON_W
const APRON_X_HI = FLANK_X_HI + APRON_W

# Trapezoidal block + toe-apron skirt:
#   - inside [APRON_X_LO, FLANK_X_LO): hardbottom along natural seabed (apron)
#   - inside [FLANK_X_LO, CREST_X_LO): linear ramp up to crest
#   - inside [CREST_X_LO, CREST_X_HI]: flat crest at -1 m
#   - inside (CREST_X_HI, FLANK_X_HI]: linear ramp down to natural seabed
#   - inside (FLANK_X_HI, APRON_X_HI]: hardbottom along natural seabed (apron)
function break_z(xi, z_natural)
    if xi < APRON_X_LO || xi > APRON_X_HI
        return -1e30                     # no hardbottom
    elseif xi < FLANK_X_LO
        return z_natural                 # seaward toe apron
    elseif xi < CREST_X_LO
        f = (xi - FLANK_X_LO) / FLANK_W
        return f * BREAK_CREST + (1 - f) * z_natural
    elseif xi <= CREST_X_HI
        return BREAK_CREST
    elseif xi <= FLANK_X_HI
        f = (FLANK_X_HI - xi) / FLANK_W
        return f * BREAK_CREST + (1 - f) * z_natural
    else
        return z_natural                 # landward toe apron
    end
end

for (i, xi) in pairs(x)
    bk = break_z(xi, z[i])
    if bk > -1e29 && bk > zh_break[i]
        zh_break[i] = bk
    end
end

# Initial bathymetry: sit on top of the breakwater (bed = max of natural
# profile and breakwater hardbottom, so the bed never starts below zh).
z_break = copy(z)
for (i, xi) in pairs(x)
    if FLANK_X_LO <= xi <= FLANK_X_HI && zh_break[i] > -1e29
        z_break[i] = max(z_break[i], zh_break[i])
    end
end

# ============================================================================
# 3. Storm forcing — single 7-day event with one peak
# ============================================================================
const NDAYS = 7
const NTIMES = NDAYS * 24 + 1
tbc    = collect(range(0.0, NDAYS * 86400.0; length=NTIMES))
# Stronger storm: peak Hrms = 4.0 m, peak surge = +1.0 m. With an emergent
# breakwater crest at +0.5 m, this places the still-water level above the
# crest only briefly at the storm peak — enough to drive overtopping but
# small enough that the breakwater still dominates the wave field.
hrmsbc = [0.5 + 2.5 * exp(-((t - 0.5*NDAYS*86400.0)/86400.0)^2) for t in tbc]
tpbc   = fill(10.0, NTIMES)
swlbc  = [0.0 + 0.6 * exp(-((t - 0.5*NDAYS*86400.0)/86400.0)^2) for t in tbc]

# ============================================================================
# 4. Build & run both cases
# ============================================================================
function run_case(z0, zh, label, ncfile)
    cfg = build_config(
        dx=dx, bathymetry_x=x, bathymetry_z=copy(z0), friction=0.003,
        timebc=tbc, tpbc=tpbc, hrmsbc=hrmsbc, swlbc=swlbc,
        options=OptionFlags(iprofl=1, isedav=1),
        sediment=make_sediment(d50=0.3e-3),
        hardbottom_z=zh,
    )
    println("Running case: $label")
    @time st = run_simulation!(cfg; outfile=ncfile, output_interval_s=3600.0)
    return st
end

const NC_A = joinpath(OUTDIR, "no_break.nc")
const NC_B = joinpath(OUTDIR, "with_break.nc")
st_A = run_case(z,       zh_base,  "no breakwater",   NC_A)
st_B = run_case(z_break, zh_break, "with breakwater", NC_B)

# ============================================================================
# 5. Diagnostics — shoreline position, bluff toe retreat, volume change
# ============================================================================
function shoreline_pos(zb, x)
    # Position where bed crosses MSL (z = 0). Linear interpolate between
    # adjacent nodes that bracket the crossing.
    for j in 1:(length(x) - 1)
        if zb[j] < 0.0 && zb[j+1] >= 0.0
            return x[j] + dx * (-zb[j]) / (zb[j+1] - zb[j])
        end
    end
    return NaN
end

function bluff_toe_pos(zb, x; threshold = 1.0)
    # Position where bed first exceeds threshold (default +1 m above MSL).
    # Walks from offshore landward.
    for j in 1:length(x)
        if zb[j] >= threshold
            return x[j]
        end
    end
    return NaN
end

zb0_A = z
zb0_B = z_break
zb1_A = st_A.zb[1:length(x), 1]
zb1_B = st_B.zb[1:length(x), 1]

shore_A0 = shoreline_pos(zb0_A, x);  shore_A1 = shoreline_pos(zb1_A, x)
shore_B0 = shoreline_pos(zb0_B, x);  shore_B1 = shoreline_pos(zb1_B, x)
toe_A0   = bluff_toe_pos(zb0_A, x);  toe_A1   = bluff_toe_pos(zb1_A, x)
toe_B0   = bluff_toe_pos(zb0_B, x);  toe_B1   = bluff_toe_pos(zb1_B, x)

# Volume above MSL (proxy for sub-aerial sand budget, integrating max(zb,0))
function vol_above_msl(zb)
    v = 0.0
    for zj in zb; if zj > 0; v += zj * dx; end; end
    return v
end
vol_A0 = vol_above_msl(zb0_A); vol_A1 = vol_above_msl(zb1_A)
vol_B0 = vol_above_msl(zb0_B); vol_B1 = vol_above_msl(zb1_B)

println()
println("="^70)
println("RESULTS (7-day storm, d50 = 0.3 mm)")
println("="^70)
@printf("%-25s %12s %12s\n", "metric", "no-break", "with-break")
@printf("%-25s %12.2f %12.2f m\n",   "Shoreline retreat",
    shore_A0 - shore_A1, shore_B0 - shore_B1)
@printf("%-25s %12.2f %12.2f m\n",   "Bluff toe retreat",
    toe_A1 - toe_A0,   toe_B1 - toe_B0)
@printf("%-25s %12.2f %12.2f m³/m\n", "ΔVol above MSL",
    vol_A1 - vol_A0,   vol_B1 - vol_B0)
println()
@printf("Breakwater protective benefit:\n")
@printf("  Shoreline retreat reduced by %+.2f m  (%.1f%%)\n",
    (shore_A0 - shore_A1) - (shore_B0 - shore_B1),
    100 * ((shore_A0 - shore_A1) - (shore_B0 - shore_B1)) / max((shore_A0 - shore_A1), 1e-9))
toe_A_retreat = toe_A1 - toe_A0
toe_B_retreat = toe_B1 - toe_B0
toe_pct = abs(toe_A_retreat) > 0.05 ? 100*(toe_A_retreat - toe_B_retreat)/toe_A_retreat : 0.0
@printf("  Bluff toe retreat reduced by %+.2f m  (%.1f%%)\n",
    toe_A_retreat - toe_B_retreat, toe_pct)
@printf("  Sub-aerial vol delta:        %+.2f m³/m\n",
    (vol_B1 - vol_B0) - (vol_A1 - vol_A0))

# ============================================================================
# 6. Comparison plot
# ============================================================================
const CM = CairoMakie
fig = CM.Figure(size = (1300, 800))

# --- Top: full profile overlay ---
ax1 = CM.Axis(fig[1, 1:2],
    xlabel="Cross-shore distance (m)",
    ylabel="Elevation (m)",
    title="Bed profile — initial vs after 7-day storm")
CM.lines!(ax1, x, zb0_A, color=:black, linewidth=2.0, linestyle=:dash, label="Initial (no break)")
CM.lines!(ax1, x, zb0_B, color=:gray,  linewidth=1.5, linestyle=:dot,  label="Initial (with break)")
CM.lines!(ax1, x, zb1_A, color=:tomato, linewidth=2.0, label="Final — no breakwater")
CM.lines!(ax1, x, zb1_B, color=:dodgerblue, linewidth=2.0, label="Final — with breakwater")
CM.hlines!(ax1, [0.0], color=(:steelblue, 0.4), linestyle=:dash, label="MSL")
# Highlight breakwater footprint
CM.text!(ax1, (CREST_X_LO + CREST_X_HI)/2, BREAK_CREST + 0.6,
    text="breakwater", color=:dodgerblue, fontsize=10, align=(:center, :bottom))
CM.ylims!(ax1, -10, 7)
CM.axislegend(ax1, position=:lt)

# --- Bottom-left: zoomed nearshore comparison ---
ax2 = CM.Axis(fig[2, 1],
    xlabel="Cross-shore distance (m)", ylabel="Elevation (m)",
    title="Beach + bluff zoom (final)")
CM.lines!(ax2, x, zb0_A, color=:black, linewidth=1.0, linestyle=:dash, label="Initial")
CM.lines!(ax2, x, zb1_A, color=:tomato, linewidth=2.0, label="No breakwater")
CM.lines!(ax2, x, zb1_B, color=:dodgerblue, linewidth=2.0, label="With breakwater")
CM.hlines!(ax2, [0.0], color=(:steelblue, 0.4), linestyle=:dash)
CM.xlims!(ax2, 380, 460)
CM.ylims!(ax2, -1, 6)
CM.axislegend(ax2, position=:lt)

# --- Bottom-right: bar chart of metrics ---
ax3 = CM.Axis(fig[2, 2],
    title="Erosion metrics", ylabel="Retreat (m)", xticks=([1, 2], ["Shoreline", "Bluff toe"]))
metrics_A = [shore_A0 - shore_A1, toe_A1 - toe_A0]
metrics_B = [shore_B0 - shore_B1, toe_B1 - toe_B0]
CM.barplot!(ax3, [1, 2] .- 0.18, metrics_A, width=0.34, color=:tomato, label="No breakwater")
CM.barplot!(ax3, [1, 2] .+ 0.18, metrics_B, width=0.34, color=:dodgerblue, label="With breakwater")
CM.axislegend(ax3, position=:rt)

CM.Label(fig[0, :],
    "Breakwater effect on shoreline + bluff erosion (single grain, 7-day storm)",
    fontsize=14, font=:bold)

png_path = joinpath(OUTDIR, "breakwater_compare.png")
CM.save(png_path, fig)
println()
println("Plot saved to: ", png_path)
