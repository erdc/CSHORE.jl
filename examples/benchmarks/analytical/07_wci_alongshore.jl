#==============================================================================
07_wci_alongshore.jl — Iterative wave-current interaction with alongshore V

Validates the new IWCINT_ALONG=1 path: when the alongshore-current contribution
is added to the Doppler-shifted dispersion via outer Picard iteration, oblique
waves are modified by the imposed alongshore current.

PHYSICS:
  General Doppler-shifted dispersion: (σ - k⃗·U⃗)² = g k tanh(k h).
  With wave vector k⃗ = k(cos θ, sin θ) and current U⃗ = (u, v):
      σ_relative = σ - k(u cos θ + v sin θ)
  CSHORE-FORTRAN keeps only the cross-shore projection (u cos θ).
  IWCINT_ALONG=1 adds the alongshore projection (v sin θ) via an outer
  iteration: each pass uses the previous-pass vmean to compute an effective
  qdisp = qwx + vmean·sin(θ)·h that is fed back into LWAVE.

EXPECTED:
  - For oblique waves + strong following alongshore current: wave heights
    DROP slightly relative to baseline (action conservation; co-flow
    elongates wavelength, lowers height).
  - For opposing current: wave heights GROW.
  - For shore-normal waves (θ → 0): no effect (sin θ → 0).
  - Convergence in 2-4 outer iterations.

SETUP:
  Planar 1:50 beach, oblique waves at 30°, sweep imposed vbc ∈ {-1, 0, +1} m/s.
  Compare three configurations:
      A. IWCINT=0 baseline (no current effect on dispersion at all)
      B. IWCINT=1 alone (cross-shore-only Doppler — current behavior)
      C. IWCINT=1 + IWCINT_ALONG=1 (full Doppler with iteration)
==============================================================================#

import Pkg; Pkg.activate(joinpath(@__DIR__, "..", "..", ".."), io=devnull)
using CSHORE, Printf
using CSHORE: build_config, OptionFlags, make_sediment, MultifractionConfig,
              initialize_state, apply_initial_bathymetry!,
              compute_derived_constants!, compute_bed_slope!, step_bc_window!,
              CurrentInput
const CM = CSHORE.CairoMakie

const OUTDIR = joinpath(@__DIR__, "..", "output")
mkpath(OUTDIR)

# ------------------------------------------------------------------
# 1. Common setup — planar beach, strongly oblique waves
# ------------------------------------------------------------------
DX = 1.0
x  = collect(0.0:DX:400.0)
slope = 1.0 / 50.0
z0 = [-5.0 + slope * xi for xi in x]

Hrms0 = 1.0 / sqrt(2.0)
Tp = 8.0
swl = 0.0
theta0_deg = 30.0   # strongly oblique to maximise sin θ projection

ntimes = 5
timebc = collect(range(0.0, 3600.0; length = ntimes))

function run_case(; iwcint::Int, iwcint_along::Int, vbc_value::Float64)
    cfg = build_config(
        dx = DX, bathymetry_x = x, bathymetry_z = copy(z0),
        friction = 0.003,
        timebc = timebc,
        tpbc = fill(Tp, ntimes),
        hrmsbc = fill(Hrms0, ntimes),
        swlbc = fill(swl, ntimes),
        wangbc = fill(theta0_deg, ntimes),
        options = OptionFlags(
            iprofl = 0, iangle = 1, icurrent = 1,
            iwcint = iwcint, iwcint_along = iwcint_along,
        ),
        sediment = make_sediment(d50 = 0.3e-3),
        multifraction = MultifractionConfig(
            grain_sizes = [0.3e-3], nlayers = 1,
            layer_thickness = 0.1, porosity = 0.4,
            initial_fractions = [1.0]),
        current = CurrentInput(time = [timebc[1], timebc[end]],
                               vbc  = [vbc_value, vbc_value]),
    )
    state = initialize_state(cfg)
    apply_initial_bathymetry!(state, cfg)
    compute_derived_constants!(state, cfg)
    for l in 1:cfg.options.iline; compute_bed_slope!(state, cfg, l); end
    step_bc_window!(state, cfg, 1)
    jm = state.jmax[1]
    return (
        x      = state.xb[1:jm],
        h      = state.h[1:jm],
        hrms   = state.hrms[1:jm],
        vmean  = state.vmean[1:jm],
        cp     = state.cp[1:jm],
        wt     = state.wt[1:jm],
        detady = state.detady_now,
    )
end

# ------------------------------------------------------------------
# 2. Sweep
# ------------------------------------------------------------------
vbc_sweep = (-1.0, 0.0, 1.0)
configs = (
    (label = "A: IWCINT=0 (no Doppler)",                 iwcint = 0, iwcint_along = 0),
    (label = "B: IWCINT=1 (cross-shore Doppler only)",   iwcint = 1, iwcint_along = 0),
    (label = "C: IWCINT=1 + IWCINT_ALONG=1 (full)",      iwcint = 1, iwcint_along = 1),
)

results = Dict{Tuple{String,Float64}, Any}()
for cfg_spec in configs
    for vb in vbc_sweep
        r = run_case(iwcint = cfg_spec.iwcint,
                     iwcint_along = cfg_spec.iwcint_along,
                     vbc_value = vb)
        results[(cfg_spec.label, vb)] = r
    end
end

# ------------------------------------------------------------------
# 3. Diagnostics — compare Hrms in offshore region (away from breaking)
# ------------------------------------------------------------------
println("=== WCI alongshore-current effect on wave heights ===")
println("Setup: planar 1:50 beach, θ₀=30°, Hrms₀=$(round(Hrms0,digits=3)) m, Tp=$(Tp)s")
println()
@printf("%-50s   %s\n", "Configuration",
        join([@sprintf("vbc=%+.1f m/s", v) for v in vbc_sweep], "    "))
@printf("%-50s   %s\n", "─"^50,
        join(["─"^11 for _ in vbc_sweep], "    "))

# Show Hrms at a fixed shoaling node (offshore of breaking) for each case
i_probe = 100   # x = 100 m, h ≈ 3 m, well outside surf zone
for cfg_spec in configs
    row = String[]
    for vb in vbc_sweep
        r = results[(cfg_spec.label, vb)]
        push!(row, @sprintf("%.4f m", r.hrms[i_probe]))
    end
    @printf("%-50s   %s\n", "Hrms @ x=$(i_probe-1) m for $(cfg_spec.label)",
            join(row, "    "))
end

println()
println("Wave period at probe (full Doppler shifts the apparent period):")
for cfg_spec in configs
    row = String[]
    for vb in vbc_sweep
        r = results[(cfg_spec.label, vb)]
        push!(row, @sprintf("%.3f s", r.wt[i_probe]))
    end
    @printf("%-50s   %s\n", "wt @ x=$(i_probe-1) m for $(cfg_spec.label)",
            join(row, "    "))
end

# Compute height ratio (full vs cross-shore-only) at probe to highlight effect
println()
println("Hrms ratio (full IWCINT_ALONG=1 / IWCINT=1 baseline) at probe:")
for vb in vbc_sweep
    r_full = results[(configs[3].label, vb)]
    r_base = results[(configs[2].label, vb)]
    ratio  = r_full.hrms[i_probe] / r_base.hrms[i_probe]
    @printf("  vbc = %+.1f m/s :  ratio = %.4f   (%+.1f %% Δ Hrms)\n",
            vb, ratio, 100*(ratio-1))
end

# ------------------------------------------------------------------
# 4. Plot — Hrms profiles for each config × vbc combo
# ------------------------------------------------------------------
fig = CM.Figure(size=(1300, 1000), backgroundcolor=:white)

ax1 = CM.Axis(fig[1,1]; ylabel = "Hrms (m)",
              title = "Wave heights with imposed alongshore current — θ₀=30° oblique waves")
colors_cfg = (CM.RGB(0.5,0.5,0.5),  # baseline gray
              CM.RGB(0.0,0.4,0.8),  # cross-shore only (blue)
              CM.RGB(0.85,0.2,0.2)) # full (red)
linestyles = ((:solid, :dash, :dot))  # one per vbc

for (ci, cfg_spec) in enumerate(configs)
    for (vi, vb) in enumerate(vbc_sweep)
        r = results[(cfg_spec.label, vb)]
        lbl = vi == 1 ? cfg_spec.label : nothing
        CM.lines!(ax1, r.x, r.hrms;
                  color = colors_cfg[ci], linewidth = 1.5,
                  linestyle = linestyles[vi],
                  label = lbl)
    end
end
CM.axislegend(ax1; position = :rt, labelsize = 9)
CM.hidexdecorations!(ax1; grid = false)

# Panel 2: Δ Hrms relative to IWCINT=1 baseline (highlights the WCI_ALONG effect)
ax2 = CM.Axis(fig[2,1]; ylabel = "ΔHrms vs IWCINT=1 baseline (m)",
              title = "Effect of IWCINT_ALONG=1 (full Doppler) on wave heights")
for (vi, vb) in enumerate(vbc_sweep)
    r_full = results[(configs[3].label, vb)]
    r_base = results[(configs[2].label, vb)]
    delta  = r_full.hrms .- r_base.hrms
    CM.lines!(ax2, r_full.x, delta;
              color = vb > 0 ? :red : (vb < 0 ? :blue : :black),
              linewidth = 2.0,
              label = @sprintf("vbc = %+.1f m/s", vb))
end
CM.hlines!(ax2, [0.0]; color = :gray, linewidth = 0.5, linestyle = :dash)
CM.axislegend(ax2; position = :rt, labelsize = 9)
CM.hidexdecorations!(ax2; grid = false)

# Panel 3: V profile for the +1.0 m/s case across all three configurations
ax3 = CM.Axis(fig[3,1]; xlabel = "Cross-shore distance (m)",
              ylabel = "V (m/s)",
              title = "Longshore current V (vbc=+1.0 m/s, all configurations)")
for (ci, cfg_spec) in enumerate(configs)
    r = results[(cfg_spec.label, 1.0)]
    CM.lines!(ax3, r.x, r.vmean;
              color = colors_cfg[ci], linewidth = 1.8,
              label = cfg_spec.label)
end
CM.axislegend(ax3; position = :lt, labelsize = 9)

CM.rowgap!(fig.layout, 5)
outpath = joinpath(OUTDIR, "benchmark_07_wci_alongshore.png")
CM.save(outpath, fig; px_per_unit = 2)
println("\nSaved figure: $outpath")
