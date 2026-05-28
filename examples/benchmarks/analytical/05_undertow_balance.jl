#==============================================================================
05_undertow_balance.jl — Mean flow mass balance (undertow).

In a closed 1D flume (no net alongshore flow, no overtopping), the
time-averaged cross-shore volume flux must be zero when integrated from
the seaward boundary to the shoreline. The Stokes drift (wave-induced
onshore mass flux above trough) must be balanced by a return flow
(undertow) below trough level.

PHYSICS:
  Mass conservation (depth-averaged):
      Q_total = ∫_{-h}^{η} u(z) dz = Q_stokes + Q_return ≈ 0

  CSHORE captures this through `umean[j]` which includes the depth-averaged
  return flow. The Stokes-drift mass flux per unit width is:
      q_stokes ≈ (1/8) g · Hrms² / c

  The resulting undertow velocity should be:
      U_mean ≈ -q_stokes / h = -(1/8) g Hrms² / (c · h)

SETUP:
  - Steady waves on a planar slope
  - No overtopping (IOVER=0) so total depth-averaged Q_total = 0
  - Fixed bed (IPROFL=0) — just checking the hydrodynamic closure

EXPECTED:
  - |umean| increases toward the shoreline where h decreases
  - CSHORE's undertow within factor of ~2 of the analytical q_stokes/h
    (exact match not expected because CSHORE uses a different depth
     partitioning based on setup gradient, not Stokes drift directly)
==============================================================================#

import Pkg; Pkg.activate(joinpath(@__DIR__, "..", "..", ".."), io=devnull)
using CSHORE, Printf
using CSHORE: build_config, OptionFlags, make_sediment, MultifractionConfig,
              initialize_state, apply_initial_bathymetry!,
              compute_derived_constants!, compute_bed_slope!, step_bc_window!
const CM = CSHORE.CairoMakie

const GRAV = 9.81
const OUTDIR = joinpath(@__DIR__, "..", "output")
mkpath(OUTDIR)

# --- Setup ---
DX = 1.0
x = collect(0.0:DX:400.0)
slope = 1.0 / 50.0
z0 = [-5.0 + slope * xi for xi in x]

Hrms0 = 0.8
Tp = 8.0
swl = 0.0

ntimes = 5
timebc = collect(range(0.0, 3600.0; length = ntimes))

cfg = build_config(
    dx = DX, bathymetry_x = x, bathymetry_z = copy(z0),
    friction = 0.003,
    timebc = timebc,
    tpbc = fill(Tp, ntimes),
    hrmsbc = fill(Hrms0, ntimes),
    swlbc = fill(swl, ntimes),
    options = OptionFlags(iprofl = 0),
    sediment = make_sediment(d50 = 0.3e-3),
    multifraction = MultifractionConfig(
        grain_sizes = [0.3e-3], nlayers = 1,
        layer_thickness = 0.1, porosity = 0.4,
        initial_fractions = [1.0]),
)

state = initialize_state(cfg)
apply_initial_bathymetry!(state, cfg)
compute_derived_constants!(state, cfg)
for l in 1:cfg.options.iline; compute_bed_slope!(state, cfg, l); end
step_bc_window!(state, cfg, 1)
jmax = state.jmax[1]

# --- Analytical Stokes drift & expected undertow ---
# Stokes drift volume flux per unit width (2nd-order Airy):
#   q_stokes ≈ (1/8) g Hrms² / c    (for shallow-water approx where
#   Hm0 ~ √2 Hrms). In fact depth-integrated Stokes mass flux is
#   E/(ρc) = (1/8) g Hrms² / c.
q_stokes = zeros(jmax)
umean_analytical = zeros(jmax)
for j in 1:jmax
    c = state.cp[j]
    h = state.h[j]
    if c > 0 && h > 0
        q_stokes[j] = 0.125 * GRAV * state.hrms[j]^2 / c
        umean_analytical[j] = -q_stokes[j] / h
    end
end

umean_model = state.umean[1:jmax]
good = (state.hrms[1:jmax] .> 1e-3) .& (state.h[1:jmax] .> 0.2)

# Compute ratio where analytical is non-zero
ratios = Float64[]
for j in 1:jmax
    if good[j] && abs(umean_analytical[j]) > 0.01
        push!(ratios, umean_model[j] / umean_analytical[j])
    end
end
mean_ratio = isempty(ratios) ? NaN : sum(ratios) / length(ratios)
median_ratio = isempty(ratios) ? NaN : sort(ratios)[length(ratios)÷2+1]

# Sign check: undertow should be negative (seaward)
max_umean = maximum(umean_model[good])
min_umean = minimum(umean_model[good])

@printf("\n=== Benchmark 05: Undertow / Mass Balance ===\n")
@printf("Hrms0, Tp:                 %.2f m, %.1f s\n", Hrms0, Tp)
@printf("Max umean (should be ≤0):  %+.3f m/s\n", max_umean)
@printf("Min umean (peak undertow): %+.3f m/s  @ x=%.1f m\n",
        min_umean, state.xb[argmin(umean_model)])
@printf("Stokes-drift q_max:        %.4f m²/s\n", maximum(q_stokes[good]))
@printf("Analytical umean peak:     %+.3f m/s\n", minimum(umean_analytical[good]))
@printf("Ratio model/analytical:    mean=%.2f  median=%.2f\n",
        mean_ratio, median_ratio)
# Sanity: undertow sign must be negative in the breaking zone
sign_ok = min_umean < 0.0
ratio_ok = !isnan(median_ratio) && 0.3 < median_ratio < 3.0
@printf("Sign check (undertow<0):   %s\n", sign_ok ? "PASS ✓" : "FAIL ✗")
@printf("Magnitude (factor of 3):   %s\n", ratio_ok ? "PASS ✓" : "FAIL ✗")

# --- Plot ---
fig = CM.Figure(size = (1100, 900), backgroundcolor = :white)

ax1 = CM.Axis(fig[1, 1]; ylabel = "Hrms (m)",
              title = "Undertow Mass Balance: CSHORE vs Stokes-drift analytical")
CM.lines!(ax1, state.xb[1:jmax], state.hrms[1:jmax];
          color = :blue, linewidth = 2)
CM.hidexdecorations!(ax1; grid = false)

ax2 = CM.Axis(fig[2, 1]; ylabel = "Stokes flux q (m²/s)")
CM.lines!(ax2, state.xb[1:jmax], q_stokes;
          color = :darkgreen, linewidth = 2, label = "Analytical q_stokes")
CM.axislegend(ax2; position = :lt)
CM.hidexdecorations!(ax2; grid = false)

ax3 = CM.Axis(fig[3, 1]; ylabel = "Undertow umean (m/s)")
CM.lines!(ax3, state.xb[1:jmax], umean_model;
          color = :red, linewidth = 2, label = "CSHORE.jl")
CM.lines!(ax3, state.xb[1:jmax], umean_analytical;
          color = :black, linestyle = :dash, linewidth = 2,
          label = "Analytical -q_stokes/h")
CM.hlines!(ax3, [0.0]; color = :gray, linestyle = :dash)
CM.axislegend(ax3; position = :lb)
CM.hidexdecorations!(ax3; grid = false)

ax4 = CM.Axis(fig[4, 1]; xlabel = "Cross-shore distance (m)",
              ylabel = "Elevation (m)")
CM.lines!(ax4, state.xb[1:jmax], state.zb[1:jmax, 1];
          color = :black, linewidth = 2, label = "Bed")
CM.hlines!(ax4, [swl]; color = :cyan, linestyle = :dot, label = "SWL")
CM.axislegend(ax4; position = :lt)

CM.rowgap!(fig.layout, 5)
outpath = joinpath(OUTDIR, "benchmark_05_undertow_balance.png")
CM.save(outpath, fig; px_per_unit = 2)
println("Saved figure: $outpath")
