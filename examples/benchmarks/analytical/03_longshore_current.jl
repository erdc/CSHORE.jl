#==============================================================================
03_longshore_current.jl — Longuet-Higgins (1970) longshore current.

Validates the longshore momentum closure by comparing CSHORE.jl's `vmean`
to Longuet-Higgins' analytical solution for wave-driven longshore current
on a planar beach with oblique wave incidence.

PHYSICS:
  Longuet-Higgins 1970 showed that, neglecting lateral mixing, the steady
  longshore current at the breaker line is:

      V_max ≈ (5π/16) · tan(θ_b) · √(g·h_b)      (no mixing)

  where θ_b is the wave angle and h_b the depth at breaking. Lateral mixing
  smooths the profile but the peak V occurs near the break point.

  CSHORE uses a simpler, integrated closure via longshore_vstgby() that
  captures the same qualitative behaviour — peak V at the break, decay
  offshore and onshore.

SETUP:
  - Planar beach (slope 1:50)
  - Oblique waves: H=1.0 m, T=8 s, θ₀=20°
  - Breaking-induced Sxy gradient drives V in the surf zone
  - No wind, no tide, no roller (clean test)

EXPECTED:
  - V peaks near the breaker line
  - Peak magnitude O(0.5-1.0) m/s (consistent with Longuet-Higgins)
  - V decays to zero offshore and at the shoreline
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

# --- Setup: planar beach ---
DX = 1.0
x = collect(0.0:DX:400.0)
slope = 1.0 / 50.0
z0 = [-5.0 + slope * xi for xi in x]

Hrms0 = 1.0 / sqrt(2.0)   # so that Hm0 = Hrms·√2 = 1.0 m
Tp = 8.0
swl = 0.0
theta0_deg = 20.0

ntimes = 5
timebc = collect(range(0.0, 3600.0; length = ntimes))

cfg = build_config(
    dx = DX, bathymetry_x = x, bathymetry_z = copy(z0),
    friction = 0.003,
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

# --- Longuet-Higgins analytical peak ---
# Locate the breaker node: first node where qbreak exceeds a small threshold
jb = findfirst(q -> q > 0.2, state.qbreak[1:jmax])
if jb === nothing
    error("No breaking detected; increase Hrms0 or check setup")
end
h_b = state.h[jb]
theta_b = asin(state.stheta[jb])
# Longuet-Higgins no-mixing peak
V_LH = (5π / 16) * tan(theta_b) * sqrt(GRAV * h_b)

# Model peak longshore current
vmean = state.vmean[1:jmax]
V_max_model = maximum(abs, vmean)
j_peak = argmax(abs.(vmean))

@printf("\n=== Benchmark 03: Longuet-Higgins Longshore Current ===\n")
@printf("Offshore H:          %.2f m (Hrms=%.3f m)\n", Hrms0*sqrt(2), Hrms0)
@printf("Offshore angle:      %.1f°\n", theta0_deg)
@printf("Breaker depth h_b:   %.3f m\n", h_b)
@printf("Breaker angle θ_b:   %.2f°\n", rad2deg(theta_b))
@printf("V_max (LH 1970 no mixing):  %.3f m/s\n", V_LH)
@printf("V_max (CSHORE.jl):          %.3f m/s  @ x = %.1f m\n",
        V_max_model, state.xb[j_peak])
@printf("Ratio model/LH:       %.2f\n", V_max_model / V_LH)

# Sanity check: CSHORE should be in the same order of magnitude as LH.
# With lateral mixing (always present implicitly), the real peak is usually
# 50-100% of V_LH.
ratio_ok = 0.3 < V_max_model / V_LH < 1.5
@printf("Order-of-magnitude check: %s\n", ratio_ok ? "PASS ✓" : "FAIL ✗")

# --- Plot ---
fig = CM.Figure(size = (1100, 900), backgroundcolor = :white)

ax1 = CM.Axis(fig[1, 1]; ylabel = "Hrms (m)",
              title = "Longuet-Higgins 1970: Longshore Current on a Plane Beach")
CM.lines!(ax1, state.xb[1:jmax], state.hrms[1:jmax];
          color = :blue, linewidth = 2, label = "Hrms")
CM.vlines!(ax1, [state.xb[jb]]; color = :orange, linestyle = :dash,
           label = "Break point")
CM.axislegend(ax1; position = :lt)
CM.hidexdecorations!(ax1; grid = false)

ax2 = CM.Axis(fig[2, 1]; ylabel = "Qbreak (fraction)")
CM.lines!(ax2, state.xb[1:jmax], state.qbreak[1:jmax];
          color = :orange, linewidth = 2)
CM.vlines!(ax2, [state.xb[jb]]; color = :orange, linestyle = :dash)
CM.hidexdecorations!(ax2; grid = false)

ax3 = CM.Axis(fig[3, 1]; ylabel = "Longshore V (m/s)")
CM.lines!(ax3, state.xb[1:jmax], vmean;
          color = :red, linewidth = 2, label = "CSHORE.jl vmean")
CM.hlines!(ax3, [V_LH]; color = :black, linestyle = :dash,
           label = "LH 1970 no-mixing peak = $(@sprintf("%.2f", V_LH)) m/s")
CM.vlines!(ax3, [state.xb[jb]]; color = :orange, linestyle = :dash,
           label = "Break point")
CM.axislegend(ax3; position = :lt)
CM.hidexdecorations!(ax3; grid = false)

ax4 = CM.Axis(fig[4, 1]; xlabel = "Cross-shore distance (m)",
              ylabel = "Elevation (m)")
CM.lines!(ax4, state.xb[1:jmax], state.zb[1:jmax, 1];
          color = :black, linewidth = 2, label = "Bed")
CM.hlines!(ax4, [swl]; color = :cyan, linestyle = :dot, label = "SWL")
CM.vlines!(ax4, [state.xb[jb]]; color = :orange, linestyle = :dash,
           label = "Break point")
CM.axislegend(ax4; position = :lt)

CM.rowgap!(fig.layout, 5)
outpath = joinpath(OUTDIR, "benchmark_03_longshore_current.png")
CM.save(outpath, fig; px_per_unit = 2)
println("Saved figure: $outpath")
