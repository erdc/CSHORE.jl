#==============================================================================
01_linear_shoaling.jl — Linear (non-breaking) wave shoaling.

Validates LWAVE and the wave action equation by comparing CSHORE.jl's Hrms
profile to the analytical Green's Law solution for an Airy wave shoaling
over a mild slope in the absence of breaking and friction.

PHYSICS:
  For normal incidence, linear shoaling preserves the energy flux:
      E · cg · n = const     ⇒   Hrms(x) = Hrms0 · √(cg0·n0 / (cg·n))
  where cg = group velocity, n = cg/c.

SETUP:
  - Mild slope 1:200 (gentle enough to avoid steepening effects)
  - Small offshore wave (Hrms=0.2 m) so no breaking occurs anywhere
  - Zero friction (fb2=0) to isolate the shoaling physics
  - Fixed bed (iprofl=0)
  - Normal incidence (iangle=0)

EXPECTED: CSHORE Hrms matches Green's Law to within discretization error
(typically <1% in the shoaling zone).
==============================================================================#

import Pkg; Pkg.activate(joinpath(@__DIR__, "..", "..", ".."), io=devnull)
using CSHORE, Printf
using CSHORE: build_config, OptionFlags, make_sediment, MultifractionConfig,
              initialize_state, apply_initial_bathymetry!,
              compute_derived_constants!, compute_bed_slope!, step_bc_window!,
              lwave_dispersion
const CM = CSHORE.CairoMakie

const GRAV = 9.81
const OUTDIR = joinpath(@__DIR__, "..", "output")
mkpath(OUTDIR)

# --- Setup ---
DX = 2.0
x = collect(0.0:DX:1000.0)
slope = 1.0 / 200.0
z0 = [-10.0 + slope * xi for xi in x]    # from -10 m offshore to +0 m onshore

Hrms0 = 0.2      # small offshore wave — stays below breaker limit everywhere
Tp = 8.0
swl = 1.0        # 1 m SWL so the profile has a clear intersection
DAYS = 0.1       # short run — fixed bed, just check wave transform
ntimes = 5
timebc = collect(range(0.0, DAYS * 86400.0; length = ntimes))

cfg = build_config(
    dx = DX, bathymetry_x = x, bathymetry_z = copy(z0),
    friction = 0.0,                              # zero bottom friction
    timebc = timebc,
    tpbc = fill(Tp, ntimes),
    hrmsbc = fill(Hrms0, ntimes),
    swlbc = fill(swl, ntimes),
    options = OptionFlags(iprofl = 0),           # fixed bed
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

# Run a single BC step (iprofl=0 runs the wave transform once)
step_bc_window!(state, cfg, 1)
jmax = state.jmax[1]

# --- Analytical Green's Law ---
# For each node where waves have not broken, compute the exact Hrms.
# At the offshore boundary: cg0 and n0 from Airy dispersion
function wave_numbers(T::Float64, h::Float64)
    # Newton iteration for kh from ω² = gk·tanh(kh)
    ω = 2π / T
    x0 = ω^2 * h / GRAV
    kh = x0
    for _ in 1:50
        f = kh * tanh(kh) - x0
        df = tanh(kh) + kh * (1 - tanh(kh)^2)
        knew = kh - f / df
        abs(knew - kh) < 1e-10 && (kh = knew; break)
        kh = knew
    end
    k = kh / h
    c = ω / k
    n = 0.5 * (1 + 2*kh / sinh(2*kh))
    cg = n * c
    return k, c, cg, n
end

h0 = state.h[1]
_, c0, cg0, n0 = wave_numbers(Tp, h0)

Hrms_exact = zeros(jmax)
valid = trues(jmax)
for j in 1:jmax
    hj = state.h[j]
    if hj <= 0
        valid[j] = false
        continue
    end
    _, cj, cgj, nj = wave_numbers(Tp, hj)
    # Shoaling coefficient Ks = sqrt(cg0 / cgj)
    Ks = sqrt(cg0 / cgj)
    Hrms_exact[j] = Hrms0 * Ks
end

# Identify non-breaking zone: where CSHORE's qbreak is essentially zero
non_breaking = state.qbreak[1:jmax] .< 0.01
# Restrict to valid region where waves exist
good = valid .& non_breaking .& (state.hrms[1:jmax] .> 1e-6)

# Error metrics (on the non-breaking zone only)
Hrms_model = state.hrms[1:jmax]
err = abs.(Hrms_model[good] .- Hrms_exact[good])
rel_err = 100 * err ./ Hrms_exact[good]
rmse = sqrt(sum(err.^2) / count(good))
max_rel = maximum(rel_err)
mean_rel = sum(rel_err) / count(good)

@printf("\n=== Benchmark 01: Linear Shoaling ===\n")
@printf("Non-breaking nodes: %d / %d\n", count(good), jmax)
@printf("Offshore:  h0=%.2f m  Hrms0=%.3f m  cg0=%.3f m/s\n", h0, Hrms0, cg0)
@printf("RMSE (Hrms):  %.5f m\n", rmse)
@printf("Max relative error:  %.2f %%\n", max_rel)
@printf("Mean relative error: %.2f %%\n", mean_rel)
@printf("Threshold (<2%%): %s\n", max_rel < 2.0 ? "PASS ✓" : "FAIL ✗")

# --- Plot ---
fig = CM.Figure(size = (1100, 800), backgroundcolor = :white)

ax1 = CM.Axis(fig[1, 1]; ylabel = "Hrms (m)",
              title = "Linear Shoaling: CSHORE.jl vs Green's Law (analytical)")
CM.lines!(ax1, state.xb[1:jmax], Hrms_exact;
          color = :black, linewidth = 2.5, label = "Analytical (Green's Law)")
CM.lines!(ax1, state.xb[1:jmax], Hrms_model;
          color = :red, linewidth = 2, linestyle = :dash, label = "CSHORE.jl")
# Shade any breaking region
breaking_x = state.xb[1:jmax][state.qbreak[1:jmax] .>= 0.01]
if !isempty(breaking_x)
    CM.vspan!(ax1, breaking_x[1], breaking_x[end];
              color = (:orange, 0.15), label = "Breaking zone (excluded)")
end
CM.axislegend(ax1; position = :lt)
CM.hidexdecorations!(ax1; grid = false)

ax2 = CM.Axis(fig[2, 1]; ylabel = "Relative error (%)")
relerr_full = zeros(jmax)
@inbounds for j in 1:jmax
    if good[j]
        relerr_full[j] = 100 * abs(Hrms_model[j] - Hrms_exact[j]) / Hrms_exact[j]
    end
end
CM.lines!(ax2, state.xb[1:jmax], relerr_full; color = :red, linewidth = 2)
CM.hlines!(ax2, [2.0]; color = :gray, linestyle = :dash, linewidth = 1,
           label = "2% threshold")
CM.axislegend(ax2; position = :lt)
CM.hidexdecorations!(ax2; grid = false)

ax3 = CM.Axis(fig[3, 1]; xlabel = "Cross-shore distance (m)",
              ylabel = "Elevation (m)")
CM.lines!(ax3, state.xb[1:jmax], state.zb[1:jmax, 1];
          color = :black, linewidth = 2, label = "Bed")
CM.hlines!(ax3, [swl]; color = :cyan, linestyle = :dot, label = "SWL")
CM.axislegend(ax3; position = :lt)

CM.rowgap!(fig.layout, 5)
outpath = joinpath(OUTDIR, "benchmark_01_linear_shoaling.png")
CM.save(outpath, fig; px_per_unit = 2)
println("Saved figure: $outpath")
