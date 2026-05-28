#==============================================================================
02_snell_refraction.jl — Wave angle refraction via Snell's Law.

Validates that CSHORE.jl correctly tracks wave direction through shoaling
using Snell's Law:
      sin(θ) / c = constant
where θ is the wave angle relative to the shore-normal and c is the phase
velocity. Waves bend toward shore-normal as they shoal (since c decreases
in shallow water).

SETUP:
  - Mild slope (1:150), small non-breaking waves
  - Oblique incidence: offshore angle = 30°
  - Zero friction, fixed bed
  - IANGLE = 1 (oblique refraction enabled)

EXPECTED: sin(θ)/c matches the offshore value at every node within the
shoaling zone. The wave angle θ decreases from ~30° offshore to a few
degrees near the breaker line.
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
DX = 2.0
x = collect(0.0:DX:1000.0)
slope = 1.0 / 150.0
z0 = [-10.0 + slope * xi for xi in x]

Hrms0 = 0.2
Tp = 8.0
swl = 1.0
theta0_deg = 30.0         # offshore angle in degrees from shore-normal
ntimes = 5
timebc = collect(range(0.0, 8640.0; length = ntimes))

cfg = build_config(
    dx = DX, bathymetry_x = x, bathymetry_z = copy(z0),
    friction = 0.0,
    timebc = timebc,
    tpbc = fill(Tp, ntimes),
    hrmsbc = fill(Hrms0, ntimes),
    swlbc = fill(swl, ntimes),
    wangbc = fill(theta0_deg, ntimes),
    options = OptionFlags(iprofl = 0, iangle = 1),
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

# --- Snell invariant ---
# sin(θ)/c should be constant; equivalently k·sin(θ) is the invariant
# (since c = ω/k and ω is constant, sin(θ)/c ~ k·sin(θ)/ω).
# CSHORE stores wkpsin = wkp*sin(theta_offshore) (Snell invariant).
# We compute sin(θ)/c at each node and compare to the offshore value.
snell_now = zeros(jmax)
for j in 1:jmax
    cj = state.cp[j]
    sj = state.stheta[j]
    snell_now[j] = cj > 0 ? sj / cj : 0.0
end
snell_off = snell_now[1]

# Relative error in the Snell invariant
rel_err = zeros(jmax)
for j in 1:jmax
    if state.hrms[j] > 1e-6 && abs(snell_off) > 1e-10
        rel_err[j] = 100 * abs(snell_now[j] - snell_off) / abs(snell_off)
    end
end
valid = state.hrms[1:jmax] .> 1e-6
max_rel = maximum(rel_err[valid])
mean_rel = sum(rel_err[valid]) / count(valid)

# Wave angles in degrees
theta_deg = [rad2deg(asin(clamp(state.stheta[j], -1.0, 1.0))) for j in 1:jmax]

@printf("\n=== Benchmark 02: Snell Refraction ===\n")
@printf("Offshore angle:       %.1f°\n", theta0_deg)
@printf("Offshore c:           %.3f m/s\n", state.cp[1])
@printf("Offshore sin(θ)/c:    %.5f\n", snell_off)
@printf("Max Snell error:      %.3f %%\n", max_rel)
@printf("Mean Snell error:     %.3f %%\n", mean_rel)
@printf("Threshold (<1%%):      %s\n", max_rel < 1.0 ? "PASS ✓" : "FAIL ✗")
@printf("Wave angle range:     %.2f° (offshore) → %.2f° (shoreward)\n",
        theta_deg[1], theta_deg[maximum(findall(valid))])

# --- Plot ---
fig = CM.Figure(size = (1100, 900), backgroundcolor = :white)

ax1 = CM.Axis(fig[1, 1]; ylabel = "Wave angle θ (°)",
              title = "Snell's Law: sin(θ)/c = const")
CM.lines!(ax1, state.xb[1:jmax], theta_deg;
          color = :red, linewidth = 2, label = "CSHORE.jl θ")
CM.hlines!(ax1, [theta0_deg]; color = :gray, linestyle = :dash,
           label = "Offshore θ₀ = $(theta0_deg)°")
CM.axislegend(ax1; position = :rb)
CM.hidexdecorations!(ax1; grid = false)

ax2 = CM.Axis(fig[2, 1]; ylabel = "sin(θ)/c (s/m)")
CM.lines!(ax2, state.xb[1:jmax], snell_now;
          color = :red, linewidth = 2, label = "CSHORE.jl")
CM.hlines!(ax2, [snell_off]; color = :black, linestyle = :dash, linewidth = 2,
           label = "Offshore invariant")
CM.axislegend(ax2; position = :rb)
CM.hidexdecorations!(ax2; grid = false)

ax3 = CM.Axis(fig[3, 1]; ylabel = "Snell error (%)")
CM.lines!(ax3, state.xb[1:jmax], rel_err; color = :red, linewidth = 2)
CM.hlines!(ax3, [1.0]; color = :gray, linestyle = :dash,
           label = "1% threshold")
CM.axislegend(ax3; position = :lt)
CM.hidexdecorations!(ax3; grid = false)

ax4 = CM.Axis(fig[4, 1]; xlabel = "Cross-shore distance (m)",
              ylabel = "Elevation (m)")
CM.lines!(ax4, state.xb[1:jmax], state.zb[1:jmax, 1];
          color = :black, linewidth = 2, label = "Bed")
CM.hlines!(ax4, [swl]; color = :cyan, linestyle = :dot, label = "SWL")
CM.axislegend(ax4; position = :lt)

CM.rowgap!(fig.layout, 5)
outpath = joinpath(OUTDIR, "benchmark_02_snell_refraction.png")
CM.save(outpath, fig; px_per_unit = 2)
println("Saved figure: $outpath")
