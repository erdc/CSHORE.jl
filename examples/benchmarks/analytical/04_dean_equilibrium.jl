#==============================================================================
04_dean_equilibrium.jl — Dean equilibrium profile stability.

The Dean equilibrium profile (Dean 1977, 1991) is the cross-shore bathymetry
that minimises uniform wave energy dissipation per unit area:

      h(x) = A · x^(2/3)

where A is a grain-size-dependent shape parameter (0.10-0.15 for medium sand).
A beach initialised exactly on this profile should be morphodynamically
close to steady under moderate wave forcing — the onshore (bedload/
asymmetry) and offshore (undertow/suspended) transport components are
approximately balanced by construction.

PHYSICS:
  In the wet (non-swash) zone, the equilibrium condition is:
      d(Ax^(2/3)) / dx · (onshore flux) = (offshore flux)
  Dean's derivation: equal energy dissipation per unit volume → profile
  shape that is stable under monochromatic constant waves.

SETUP:
  - Dean profile with A = 0.125 (typical medium sand ~0.3 mm)
  - Constant moderate waves (Hrms=0.8, Tp=8 s)
  - 2-day run to accumulate morphodynamic change
  - Full morphodynamic run (iprofl=1) to see drift from equilibrium

EXPECTED:
  - Profile drift is SMALL (RMSE < 0.1 m over the beach face)
  - No runaway bar/trough formation
  - Any drift is systematic (wave over/under-shoaling) not numerical noise
==============================================================================#

import Pkg; Pkg.activate(joinpath(@__DIR__, "..", "..", ".."), io=devnull)
using CSHORE, Printf
using CSHORE: build_config, OptionFlags, make_sediment, MultifractionConfig,
              initialize_state, apply_initial_bathymetry!,
              compute_derived_constants!, compute_bed_slope!, step_bc_window!
const CM = CSHORE.CairoMakie

const OUTDIR = joinpath(@__DIR__, "..", "output")
mkpath(OUTDIR)

# --- Setup: Dean profile h = A·x^(2/3) ---
A = 0.125       # shape parameter, medium sand
DX = 1.0
# Place x=0 at the shoreline, extend offshore; CSHORE wants increasing x
# We'll invert: build the profile in "distance offshore from shore", then
# flip so x increases shoreward.
x_offshore = collect(5.0:DX:500.0)   # start 5 m offshore (avoid singular h=0)
h_dean = A .* x_offshore.^(2/3)
# CSHORE coordinates: x=0 offshore, x increases toward shore
x = collect(0.0:DX:(length(x_offshore)-1)*DX)
z0 = -reverse(h_dean)    # bed elevation = -water depth; flip to CSHORE orientation

Hrms0 = 0.8
Tp = 8.0
swl = 0.0

DAYS = 2.0
ntimes = DAYS * 24 + 1 |> round |> Int
timebc = collect(range(0.0, DAYS*86400.0; length = ntimes))

cfg = build_config(
    dx = DX, bathymetry_x = x, bathymetry_z = copy(z0),
    friction = 0.003,
    timebc = timebc,
    tpbc = fill(Tp, ntimes),
    hrmsbc = fill(Hrms0, ntimes),
    swlbc = fill(swl, ntimes),
    options = OptionFlags(iprofl = 1),
    sediment = make_sediment(d50 = 0.3e-3),
    multifraction = MultifractionConfig(
        grain_sizes = [0.3e-3], nlayers = 3,
        layer_thickness = 0.3, porosity = 0.4,
        initial_fractions = [1.0]),
    max_dzb_per_step = 0.05,
)

state = initialize_state(cfg)
apply_initial_bathymetry!(state, cfg)
compute_derived_constants!(state, cfg)
for l in 1:cfg.options.iline; compute_bed_slope!(state, cfg, l); end

for itime in 1:(ntimes-1)
    step_bc_window!(state, cfg, itime)
end
jmax = state.jmax[1]

# Compare final profile to initial Dean profile
zb_init = z0[1:jmax]
zb_final = state.zb[1:jmax, 1]
dzb = zb_final .- zb_init

# Metric: RMSE on the beach face (non-trivial depth)
deep_enough = state.h[1:jmax] .> 0.2
rmse = sqrt(sum(dzb[deep_enough].^2) / count(deep_enough))
max_drift = maximum(abs, dzb[deep_enough])
mean_drift = sum(abs, dzb[deep_enough]) / count(deep_enough)

@printf("\n=== Benchmark 04: Dean Equilibrium Profile ===\n")
@printf("Dean parameter A:     %.3f (medium sand)\n", A)
@printf("Duration:             %.1f days\n", DAYS)
@printf("Offshore Hrms/Tp:     %.2f m / %.1f s\n", Hrms0, Tp)
@printf("RMSE bed drift:       %.4f m\n", rmse)
@printf("Max drift:            %.4f m  @  x = %.1f m\n",
        max_drift, state.xb[argmax(abs.(dzb))])
@printf("Mean |drift|:         %.4f m\n", mean_drift)
# A stable-ish Dean profile should have drift << local depth (say < 10% of H)
threshold = 0.2
@printf("Stability threshold (RMSE < %.2f m): %s\n", threshold,
        rmse < threshold ? "PASS ✓" : "FAIL ✗")

# --- Plot ---
fig = CM.Figure(size = (1100, 900), backgroundcolor = :white)

ax1 = CM.Axis(fig[1, 1]; ylabel = "Elevation (m)",
              title = "Dean Equilibrium Profile: $(DAYS)-day morphodynamic drift")
CM.lines!(ax1, state.xb[1:jmax], zb_init;
          color = :black, linewidth = 2, linestyle = :dash,
          label = "Initial Dean h = $(A)·x^(2/3)")
CM.lines!(ax1, state.xb[1:jmax], zb_final;
          color = :red, linewidth = 2, label = "After $(DAYS) days")
CM.hlines!(ax1, [swl]; color = :cyan, linestyle = :dot, label = "SWL")
CM.axislegend(ax1; position = :lt)
CM.hidexdecorations!(ax1; grid = false)

ax2 = CM.Axis(fig[2, 1]; ylabel = "Bed drift Δzb (m)")
CM.lines!(ax2, state.xb[1:jmax], dzb; color = :red, linewidth = 2)
CM.hlines!(ax2, [0.0]; color = :gray, linestyle = :dash)
CM.hidexdecorations!(ax2; grid = false)

ax3 = CM.Axis(fig[3, 1]; ylabel = "Hrms (m)")
CM.lines!(ax3, state.xb[1:jmax], state.hrms[1:jmax];
          color = :blue, linewidth = 2)
CM.hidexdecorations!(ax3; grid = false)

ax4 = CM.Axis(fig[4, 1]; xlabel = "Cross-shore distance (m)",
              ylabel = "q_total (m²/s)")
CM.lines!(ax4, state.xb[1:jmax], state.q_total[1:jmax];
          color = :darkgreen, linewidth = 2, label = "Net transport")
CM.hlines!(ax4, [0.0]; color = :gray, linestyle = :dash)
CM.axislegend(ax4; position = :lt)

CM.rowgap!(fig.layout, 5)
outpath = joinpath(OUTDIR, "benchmark_04_dean_equilibrium.png")
CM.save(outpath, fig; px_per_unit = 2)
println("Saved figure: $outpath")
