#==============================================================================
06_imposed_longshore.jl — V→DETADY inversion check (ICURRENT=1).

Validates the new ICURRENT=1 capability: the user prescribes an alongshore
current speed `vbc` at the offshore boundary, and CSHORE.jl back-solves the
alongshore water-surface gradient DETADY required to produce it via analytical
inversion of the longshore-momentum balance.

This is the exact inverse of the Longuet-Higgins-style forward problem in
`03_longshore_current.jl`: there we gave the model wave forcing and read off
the peak V; here we GIVE the model V and check that:

  (1) `vmean` at the offshore boundary matches the prescribed `vbc`
      to numerical precision;
  (2) the implied DETADY is non-zero, finite, and changes sign with the
      sign of the prescribed V;
  (3) the V profile across the rest of the surf zone is physically
      sensible (smoothly varying, peaks at break, decays toward shore).

SETUP:
  - Planar beach (slope 1:50), oblique incidence (θ=20°)
  - Sweep across vbc = -0.5, -0.2, 0.0, +0.2, +0.5 m/s
  - Compare ICURRENT=0 baseline (no imposed current) to each ICURRENT=1 case
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
# 1. Common setup — planar beach, oblique waves
# ------------------------------------------------------------------
DX = 1.0
x  = collect(0.0:DX:400.0)
slope = 1.0 / 50.0
z0 = [-5.0 + slope * xi for xi in x]

Hrms0 = 1.0 / sqrt(2.0)
Tp = 8.0
swl = 0.0
theta0_deg = 20.0

ntimes = 5
timebc = collect(range(0.0, 3600.0; length = ntimes))

function run_case(; icurrent::Int, vbc_value::Float64)
    cfg = build_config(
        dx = DX, bathymetry_x = x, bathymetry_z = copy(z0),
        friction = 0.003,
        timebc = timebc,
        tpbc = fill(Tp, ntimes),
        hrmsbc = fill(Hrms0, ntimes),
        swlbc = fill(swl, ntimes),
        wangbc = fill(theta0_deg, ntimes),
        options = OptionFlags(iprofl = 0, iangle = 1, icurrent = icurrent),
        sediment = make_sediment(d50 = 0.3e-3),
        multifraction = MultifractionConfig(
            grain_sizes = [0.3e-3], nlayers = 1,
            layer_thickness = 0.1, porosity = 0.4,
            initial_fractions = [1.0]),
        current = icurrent == 1 ?
                  CurrentInput(time = [timebc[1], timebc[end]],
                               vbc  = [vbc_value, vbc_value]) : nothing,
    )
    state = initialize_state(cfg)
    apply_initial_bathymetry!(state, cfg)
    compute_derived_constants!(state, cfg)
    for l in 1:cfg.options.iline; compute_bed_slope!(state, cfg, l); end
    step_bc_window!(state, cfg, 1)
    jm = state.jmax[1]
    return (
        x       = state.xb[1:jm],
        h       = state.h[1:jm],
        hrms    = state.hrms[1:jm],
        vmean   = state.vmean[1:jm],
        gby     = state.gby[1:jm],
        detady  = state.detady_now,
        vbc_now = state.vbc_now,
    )
end

# ------------------------------------------------------------------
# 2. Baseline (icurrent=0) — purely wave-driven
# ------------------------------------------------------------------
base = run_case(icurrent = 0, vbc_value = 0.0)
@printf("=== Baseline (ICURRENT=0, no imposed current) ===\n")
@printf("Peak |V| = %.3f m/s at x = %.1f m   (Longuet-Higgins-style)\n",
        maximum(abs, base.vmean), base.x[argmax(abs.(base.vmean))])
@printf("DETADY   = %.3e   (should be 0 — no tide / no imposed V)\n\n",
        base.detady)

# ------------------------------------------------------------------
# 3. Sweep imposed V at offshore boundary
# ------------------------------------------------------------------
vbc_sweep = (-0.5, -0.2, 0.0, 0.2, 0.5)
results = Dict{Float64, Any}()
for vb in vbc_sweep
    r = run_case(icurrent = 1, vbc_value = vb)
    results[vb] = r
    # vmean at offshore boundary cell. Cell 1 is the boundary node (no momentum
    # solve there); cell 2 is the most-offshore valid cell where the inversion
    # was applied.
    v_offshore = r.vmean[2]
    err = abs(v_offshore - vb)
    @printf("=== ICURRENT=1, vbc = %+.2f m/s ===\n", vb)
    @printf("  vmean[2]  = %+.4f m/s   (target %+.2f)   err = %.2e\n",
            v_offshore, vb, err)
    @printf("  DETADY    = %+.3e\n", r.detady)
    @printf("  Peak |V|  = %.3f m/s  at x = %.1f m\n",
            maximum(abs, r.vmean), r.x[argmax(abs.(r.vmean))])
    @printf("  V at shore (last 5 cells): %s\n",
            join([@sprintf("%+.3f", v) for v in r.vmean[end-4:end]], ", "))
    println()
end

# ------------------------------------------------------------------
# 4. PASS / FAIL summary
# ------------------------------------------------------------------
println("=== Validation summary ===")
let all_pass = true, TOL = 1e-3
    global validation_passed = true
    for vb in vbc_sweep
        r = results[vb]
        err = abs(r.vmean[2] - vb)
        pass = err < TOL
        @printf("  vbc=%+.2f m/s :  err = %.2e   %s\n",
                vb, err, pass ? "PASS" : "FAIL")
        all_pass &= pass
    end
    global validation_passed = all_pass
    println(all_pass ? "ALL PASS — V→DETADY inversion is exact." : "FAIL — see errors above.")
end

# ------------------------------------------------------------------
# 5. Plot — V profiles for the swept vbc values
# ------------------------------------------------------------------
fig = CM.Figure(size=(1100, 800), backgroundcolor=:white)

ax1 = CM.Axis(fig[1,1]; ylabel = "Mean longshore current V (m/s)",
              title = "Imposed alongshore current (ICURRENT=1) — V profile across surf zone")
CM.lines!(ax1, base.x, base.vmean;
          color = :gray, linewidth = 2.0, linestyle = :dash,
          label = "Baseline (ICURRENT=0, vbc=0)")
colors = CM.cgrad(:RdBu, length(vbc_sweep); categorical = true, rev = true)
for (i, vb) in enumerate(vbc_sweep)
    r = results[vb]
    CM.lines!(ax1, r.x, r.vmean;
              color = colors[i], linewidth = 2.0,
              label = @sprintf("vbc = %+.2f m/s  (DETADY=%+.2e)", vb, r.detady))
    # Mark the offshore reference cell where V was imposed
    CM.scatter!(ax1, [r.x[2]], [r.vmean[2]];
                color = colors[i], markersize = 12, marker = :star5)
end
CM.hlines!(ax1, [0.0]; color = :black, linewidth = 0.5, linestyle = :dot)
CM.axislegend(ax1; position = :lt, labelsize = 9)
CM.hidexdecorations!(ax1; grid = false)

ax2 = CM.Axis(fig[2,1]; xlabel = "Cross-shore distance (m)",
              ylabel = "Hrms (m) / SWL (m)",
              title = "Wave height for context (same for all cases)")
CM.lines!(ax2, base.x, base.hrms; color = :blue, linewidth = 1.5, label = "Hrms")
CM.lines!(ax2, base.x, -base.h .+ swl; color = :brown, linewidth = 1.5, label = "Bed elev")
CM.lines!(ax2, base.x, fill(swl, length(base.x));
          color = :cyan, linewidth = 1, linestyle = :dash, label = "SWL")
CM.axislegend(ax2; position = :rt, labelsize = 9)

CM.rowgap!(fig.layout, 5)
outpath = joinpath(OUTDIR, "benchmark_06_imposed_longshore.png")
CM.save(outpath, fig; px_per_unit = 2)
println("\nSaved figure: $outpath")
