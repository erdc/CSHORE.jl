#==============================================================================
lip_11d_test_1c.jl — LIP 11D Delta Flume, Test 1C (accretive case)

The LIP 11D experiments (Roelvink & Reniers 1995) in Delft Hydraulics'
Delta Flume are THE standard benchmark for cross-shore morphodynamic
models. Three tests covered stable (1A), erosive (1B), and accretive (1C)
beach states. Test 1C is the classic onshore bar migration case.

This script uses the **actual XBeach skillbed input files**, fetched from
https://svn.oss.deltares.nl/repos/xbeach/skillbed/input/DeltaflumeLIP11D/1C/
  - h_1C.dep       : water depths (201 nodes, 1 m spacing) — positive down
  - jonswap.inp    : Hm0=0.58 m, fp=0.125 Hz (Tp=8s), mainang=270°
  - params_original.txt : grid, flow, sed options used by XBeach

Running from the published inputs rather than a synthesised profile lets
us benchmark CSHORE.jl against the same forcing XBeach itself uses.

CITATION:
  Roelvink, J.A. & Reniers, A.J.H.M. (1995). *LIP 11D Delta Flume
  experiments: A dataset for profile model validation*. Delft Hydraulics
  report H2130.

TEST 1C CONDITIONS (from params_original.txt + jonswap.inp):
  - Flume: 201 nodes × 1 m = 200 m long working section
  - Initial profile: h_1C.dep (digitised LIP-1B final state)
  - Waves: Hm0 = 0.58 m, Tp = 1/fp = 8.0 s, JONSWAP γ=3.3
  - SWL: zs0 = -0.008 m ≈ 0 (at flume datum)
  - Duration: tstop-tstart = 29800-1000 = 28800 s = 8 hr waves
  - Sediment: d50=0.22 mm, d90=0.33 mm, porosity=0.4

EXPECTED (from published XBeach validation):
  - Bar migrates ~5-10 m onshore over the 8 hr run
  - Classic bed-change dipole: erosion seaward, accretion shoreward
==============================================================================#

import Pkg; Pkg.activate(joinpath(@__DIR__, "..", "..", ".."), io=devnull)
using CSHORE, Printf, DelimitedFiles
using CSHORE: build_config, OptionFlags, make_sediment, MultifractionConfig,
              initialize_state, apply_initial_bathymetry!,
              compute_derived_constants!, compute_bed_slope!, step_bc_window!
const CM = CSHORE.CairoMakie

const OUTDIR = joinpath(@__DIR__, "..", "output")
mkpath(OUTDIR)

const DATADIR = joinpath(@__DIR__, "data", "DeltaflumeLIP11D_1C")

# ------------------------------------------------------------------
# 1. Load bathymetry from XBeach .dep file
# ------------------------------------------------------------------
# The .dep file stores positive-down depths (XBeach convention: "posdwn = -1"
# in this case actually means z-coordinates, negative below datum → we see
# negative values in the file → those are actual bed elevations).
# From the file inspection: values are negative (~-4.1 to +0.8), consistent
# with posdwn=-1 meaning the file contains bed ELEVATIONS, not depths.
# So zb = values directly.
bathy_raw = vec(readdlm(joinpath(DATADIR, "h_1C.dep")))

# XBeach grid: 201 nodes, dx=1 m, x from 0 to 200 m
# But XBeach's convention has x increasing from offshore → onshore,
# and CSHORE uses the same convention, so no flip needed.
DX = 1.0
x = collect(0.0:DX:(length(bathy_raw) - 1) * DX)
z0 = collect(bathy_raw)

# Sanity check: offshore end should be deep (~-4 m), onshore shallow (>0)
println("=== Loaded LIP 11D Test 1C bathymetry ===")
@printf("Grid:       %d nodes × %.1f m = %.1f m flume\n",
        length(x), DX, last(x))
@printf("Elevation:  %.2f (offshore) → %.2f (onshore) m\n", z0[1], z0[end])

# ------------------------------------------------------------------
# 2. Wave conditions from jonswap.inp and params
# ------------------------------------------------------------------
# jonswap.inp: Hm0=0.58 m, fp=0.125 Hz, γ=3.3
# CSHORE uses Hrms (= Hm0/√2 for Rayleigh-distributed heights)
Hm0 = 0.58
Hrms = Hm0 / sqrt(2.0)
fp = 0.125
Tp = 1.0 / fp      # = 8.0 s

# params_original.txt: zs0 = -0.008 m
SWL = -0.008

# Run duration (waves active): tstart=1000, tstop=29800 → 28800 s = 8 hr
DURATION_S = 29800.0 - 1000.0
DURATION_H = DURATION_S / 3600.0

# Use hourly BC windows for morphodynamic stepping
ntimes = round(Int, DURATION_H) + 1
timebc = collect(range(0.0, DURATION_S; length=ntimes))

@printf("\n=== Wave conditions (from jonswap.inp) ===\n")
@printf("Hm0:        %.2f m  →  Hrms=%.3f m\n", Hm0, Hrms)
@printf("Tp:         %.2f s\n", Tp)
@printf("SWL:        %+.3f m\n", SWL)
@printf("Duration:   %.1f hr\n", DURATION_H)

# ------------------------------------------------------------------
# 3. Sediment (params_original.txt: d50=0.22 mm, porosity=0.4)
# ------------------------------------------------------------------
d50 = 0.22e-3
poro = 0.4

mf = MultifractionConfig(
    grain_sizes = [d50],
    nlayers = 3,
    layer_thickness = 0.3,
    porosity = poro,
    initial_fractions = [1.0],
)

# ------------------------------------------------------------------
# 4. Run three configurations to bracket the physics
# ------------------------------------------------------------------
# 1. Linear waves, no roller            — baseline
# 2. + Ruessink skewness (IASYM=1)      — standard XBeach-like setup
# 3. + Roller (IROLL=1) + Skewness       — full nonlinear physics

function run_case(label; iasym=0, iroll=0, facSK=1.0)
    cfg = build_config(
        dx = DX,
        bathymetry_x = x,
        bathymetry_z = copy(z0),
        friction = 0.003,
        timebc = timebc,
        tpbc = fill(Tp, ntimes),
        hrmsbc = fill(Hrms, ntimes),
        swlbc = fill(SWL, ntimes),
        options = OptionFlags(iprofl=1, iasym=iasym, iroll=iroll),
        sediment = make_sediment(d50=d50, sporo=poro),
        multifraction = mf,
        facSK = facSK,
        max_dzb_per_step = 0.05,
    )
    state = initialize_state(cfg)
    apply_initial_bathymetry!(state, cfg)
    compute_derived_constants!(state, cfg)
    for l in 1:cfg.options.iline; compute_bed_slope!(state, cfg, l); end
    for itime in 1:(ntimes-1); step_bc_window!(state, cfg, itime); end
    return state, state.jmax[1]
end

println("\n=== Running LIP 11D Test 1C ($(DURATION_H) hr) ===")
print("  Linear baseline... "); t0 = time()
s_lin, jm = run_case("Linear"; iasym=0, iroll=0)
@printf("%.1fs\n", time() - t0)

print("  + Ruessink skewness... "); t0 = time()
s_sk, _   = run_case("Skewness"; iasym=1, iroll=0, facSK=1.0)
@printf("%.1fs\n", time() - t0)

print("  + Roller + Skewness... "); t0 = time()
s_all, _  = run_case("Roller+Skew"; iasym=1, iroll=1, facSK=1.0)
@printf("%.1fs\n", time() - t0)

# ------------------------------------------------------------------
# 5. Bar crest tracking + diagnostics
# ------------------------------------------------------------------
xp = s_lin.xb[1:jm]

# Locate bar between x=120 m and x=160 m
function bar_crest(xb, zb; xrange=(110.0, 165.0))
    js = findall(xi -> xrange[1] <= xi <= xrange[2], xb)
    isempty(js) && return NaN, NaN, -1
    imax = js[argmax(zb[js])]
    return xb[imax], zb[imax], imax
end

# Bar-region centroid of accretion — more robust bar-migration metric.
# Computed as mass-weighted centroid of positive bed change within window.
function bar_centroid(xb, dzb; xrange=(115.0, 165.0))
    js = findall(xi -> xrange[1] <= xi <= xrange[2], xb)
    isempty(js) && return NaN
    # Only use accreted mass (positive dzb)
    mass = max.(dzb[js], 0.0)
    total = sum(mass)
    total < 1e-6 && return NaN
    return sum(xb[js] .* mass) / total
end

# Bar is centered around x=135-140 in the LIP-1B final profile; use a
# tight window that excludes the beach-face climb (starts around x=160).
bar_window = (115.0, 155.0)
xbar0, zbar0, _ = bar_crest(xp, z0[1:jm]; xrange=bar_window)
xbar_lin, _, _  = bar_crest(xp, s_lin.zb[1:jm, 1]; xrange=bar_window)
xbar_sk, _, _   = bar_crest(xp, s_sk.zb[1:jm, 1]; xrange=bar_window)
xbar_all, _, _  = bar_crest(xp, s_all.zb[1:jm, 1]; xrange=bar_window)

@printf("\n=== Bar Crest Tracking ===\n")
@printf("  Initial bar crest:    x = %.1f m,  z = %+.3f m\n", xbar0, zbar0)
@printf("  Linear (iasym=0):     crest x = %.1f m  (ΔX = %+.1f m)\n",
        xbar_lin, xbar_lin - xbar0)
@printf("  + Skewness:           crest x = %.1f m  (ΔX = %+.1f m)\n",
        xbar_sk, xbar_sk - xbar0)
@printf("  + Roller + Skew:      crest x = %.1f m  (ΔX = %+.1f m)\n",
        xbar_all, xbar_all - xbar0)

# Centroid-based onshore migration metric (more robust than crest position)
dzb_lin_arr = s_lin.zb[1:jm,1] .- z0[1:jm]
dzb_sk_arr  = s_sk.zb[1:jm,1]  .- z0[1:jm]
dzb_all_arr = s_all.zb[1:jm,1] .- z0[1:jm]
xcen_lin = bar_centroid(xp, dzb_lin_arr)
xcen_sk  = bar_centroid(xp, dzb_sk_arr)
xcen_all = bar_centroid(xp, dzb_all_arr)

@printf("\n=== Accretion Centroid (bar-migration metric) ===\n")
@printf("  Initial bar crest location:  %.1f m\n", xbar0)
@printf("  Linear accretion centroid:   %.1f m  (ΔX = %+.1f m onshore)\n",
        xcen_lin, xcen_lin - xbar0)
@printf("  Skewness accretion centroid: %.1f m  (ΔX = %+.1f m onshore)\n",
        xcen_sk, xcen_sk - xbar0)
@printf("  Full physics centroid:       %.1f m  (ΔX = %+.1f m onshore)\n",
        xcen_all, xcen_all - xbar0)
@printf("  (Published: ~5-10 m onshore migration over 13 hr run)\n")

# Volume conservation check
dx_ = DX
vol_lin = sum(s_lin.zb[1:jm, 1] .- z0[1:jm]) * dx_
vol_sk  = sum(s_sk.zb[1:jm, 1]  .- z0[1:jm]) * dx_
vol_all = sum(s_all.zb[1:jm, 1] .- z0[1:jm]) * dx_
@printf("\n=== Volume Conservation (should be ~0) ===\n")
@printf("  Linear:      %+.5f m²\n", vol_lin)
@printf("  Skewness:    %+.5f m²\n", vol_sk)
@printf("  Roller+Skew: %+.5f m²\n", vol_all)

# Max bed change
@printf("\n=== Max |Δzb| ===\n")
@printf("  Linear:      %.3f m\n", maximum(abs, s_lin.zb[1:jm,1].-z0[1:jm]))
@printf("  Skewness:    %.3f m\n", maximum(abs, s_sk.zb[1:jm,1] .-z0[1:jm]))
@printf("  Roller+Skew: %.3f m\n", maximum(abs, s_all.zb[1:jm,1].-z0[1:jm]))

# ------------------------------------------------------------------
# 6. Plot
# ------------------------------------------------------------------
fig = CM.Figure(size=(1200, 1100), backgroundcolor=:white)

ax1 = CM.Axis(fig[1,1];
    ylabel = "Elevation (m)",
    title  = "LIP 11D Test 1C — $(DURATION_H)hr, Hm0=$(Hm0)m, Tp=$(Tp)s (XBeach skillbed data)")
CM.lines!(ax1, xp, z0[1:jm];
          color=:gray, linewidth=1.5, linestyle=:dash, label="Initial (XBeach h_1C.dep)")
CM.lines!(ax1, xp, s_lin.zb[1:jm,1]; color=:black,     linewidth=2, label="Linear")
CM.lines!(ax1, xp, s_sk.zb[1:jm,1];  color=:blue,      linewidth=2, label="+Skewness")
CM.lines!(ax1, xp, s_all.zb[1:jm,1]; color=:red,       linewidth=2, label="+Roller+Skew")
CM.vlines!(ax1, [xbar0];    color=:gray,  linestyle=:dot, linewidth=1)
CM.vlines!(ax1, [xbar_all]; color=:red,   linestyle=:dot, linewidth=1)
CM.hlines!(ax1, [SWL]; color=:cyan, linewidth=1, linestyle=:dot, label="SWL")
CM.xlims!(ax1, 100, 200)
CM.axislegend(ax1; position=:lt, labelsize=10)
CM.hidexdecorations!(ax1; grid=false)

# Bed change
dzb_lin = s_lin.zb[1:jm,1] .- z0[1:jm]
dzb_sk  = s_sk.zb[1:jm,1]  .- z0[1:jm]
dzb_all = s_all.zb[1:jm,1] .- z0[1:jm]
ax2 = CM.Axis(fig[2,1]; ylabel="Bed change Δzb (m)")
CM.lines!(ax2, xp, dzb_lin; color=:black, linewidth=2, label="Linear")
CM.lines!(ax2, xp, dzb_sk;  color=:blue,  linewidth=2, label="+Skewness")
CM.lines!(ax2, xp, dzb_all; color=:red,   linewidth=2, label="+Roller+Skew")
CM.hlines!(ax2, [0.0]; color=:gray, linewidth=0.5, linestyle=:dash)
CM.xlims!(ax2, 100, 200)
CM.axislegend(ax2; position=:lt, labelsize=10)
CM.hidexdecorations!(ax2; grid=false)

# Wave height
ax3 = CM.Axis(fig[3,1]; ylabel="Hrms (m)")
CM.lines!(ax3, xp, s_lin.hrms[1:jm]; color=:black, linewidth=2, label="Linear")
CM.lines!(ax3, xp, s_sk.hrms[1:jm];  color=:blue,  linewidth=2, label="+Skewness")
CM.lines!(ax3, xp, s_all.hrms[1:jm]; color=:red,   linewidth=2, label="+Roller+Skew")
CM.xlims!(ax3, 100, 200)
CM.axislegend(ax3; position=:lb, labelsize=10)
CM.hidexdecorations!(ax3; grid=false)

# Undertow
ax4 = CM.Axis(fig[4,1]; ylabel="umean (m/s)")
CM.lines!(ax4, xp, s_lin.umean[1:jm]; color=:black, linewidth=2, label="Linear")
CM.lines!(ax4, xp, s_sk.umean[1:jm];  color=:blue,  linewidth=2, label="+Skewness")
CM.lines!(ax4, xp, s_all.umean[1:jm]; color=:red,   linewidth=2, label="+Roller+Skew")
CM.hlines!(ax4, [0.0]; color=:gray, linewidth=0.5, linestyle=:dash)
CM.xlims!(ax4, 100, 200)
CM.axislegend(ax4; position=:lb, labelsize=10)
CM.hidexdecorations!(ax4; grid=false)

# Full profile inset
ax5 = CM.Axis(fig[5,1]; xlabel="Cross-shore distance (m)",
              ylabel="Elevation (m)",
              title="Full flume profile (200 m)")
CM.lines!(ax5, xp, z0[1:jm];        color=:gray, linewidth=1.5, linestyle=:dash)
CM.lines!(ax5, xp, s_all.zb[1:jm,1]; color=:red, linewidth=2)
CM.hlines!(ax5, [SWL]; color=:cyan, linewidth=1, linestyle=:dot)

CM.rowgap!(fig.layout, 5)
outpath = joinpath(OUTDIR, "benchmark_xbeach_lip_11d_1c.png")
CM.save(outpath, fig; px_per_unit=2)
println("\nSaved figure: $outpath")
