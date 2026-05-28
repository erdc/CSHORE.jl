"""
Plot cross-shore profile comparisons for all ERDC USACE benchmark cases.

Usage:
  julia --project=../../.. examples/benchmarks/erdc_usace/plot_benchmarks.jl

Produces: examples/benchmarks/erdc_usace/erdc_comparison.png
Each row shows one case with four panels:
  (a) Hrms   (b) wave setup   (c) bed-level change   (d) bed profile
"""

include(joinpath(@__DIR__, "compare_utils.jl"))
using CairoMakie

outpath = joinpath(@__DIR__, "erdc_comparison.png")

ncases = length(CASES)
fig = Figure(size=(1400, 290 * ncases), fontsize=12)
Label(fig[0, 1:4];
    text="Julia CSHORE.jl  vs  FORTRAN CSHORE_USACE — ERDC testbed cases",
    fontsize=15, font=:bold)

ax_kw = (xgridvisible=true, ygridvisible=true,
         xgridcolor=:lightgray, ygridcolor=:lightgray)

for (row, case) in enumerate(CASES)
    # ── Load reference ────────────────────────────────────────────────────────
    ref_os = read_ref_block(joinpath(case.dir, "ref_OSETUP.txt"); ncols=4)
    ref_ob = read_ref_block(joinpath(case.dir, "ref_OBPROF.txt"); ncols=2)

    ft_x_wave = ref_os[:, 1]
    ft_hrms   = sqrt(8.0) .* ref_os[:, 4]
    ft_setup  = ref_os[:, 2] .- case.swl
    ft_x_prof = ref_ob[:, 1]
    ft_zb     = ref_ob[:, 2]

    # Initial bed (from infile)
    cfg = read_infile(joinpath(case.dir, "infile"); strict=false)
    nb0   = cfg.bathymetry.nbinp[1]
    x0    = cfg.bathymetry.xbinp[1:nb0, 1]
    zb0   = cfg.bathymetry.zbinp[1:nb0, 1]
    ft_dzb = ft_zb .- zb0[1:length(ft_zb)]

    # ── Run Julia ─────────────────────────────────────────────────────────────
    state = run_simulation!(cfg)
    jmax = state.jmax[1]; jr = state.jr
    x_jl    = state.xb[1:jmax]
    hrms_jl = state.hrms[1:jmax]
    setup_jl = state.wsetup[1:jmax]
    zb_jl   = state.zb[1:jmax, 1]
    dzb_jl  = zb_jl .- zb0[1:jmax]

    n_w = min(length(ft_x_wave), jr)
    n_p = min(length(ft_x_prof), jmax)

    rms_h = sqrt(sum((ft_hrms[1:n_w] .- hrms_jl[1:n_w]).^2) / n_w)
    rms_z = sqrt(sum((ft_dzb[1:n_p]  .- dzb_jl[1:n_p]).^2)  / n_p)

    # ── Panel A: Hrms ─────────────────────────────────────────────────────────
    ax1 = Axis(fig[row, 1];
        title="$(case.name) — Hrms [RMS=$(round(rms_h,digits=3)) m]",
        xlabel="x (m)", ylabel="Hrms (m)", ax_kw...)
    lines!(ax1, ft_x_wave[1:n_w], ft_hrms[1:n_w];
           color=:steelblue, linewidth=2.5, label="FORTRAN")
    lines!(ax1, x_jl[1:n_w], hrms_jl[1:n_w];
           color=:orangered, linewidth=2, linestyle=:dash, label="Julia")
    axislegend(ax1; position=:lb, labelsize=10, framevisible=false)

    # ── Panel B: Wave setup ───────────────────────────────────────────────────
    ax2 = Axis(fig[row, 2];
        title="$(case.name) — Wave setup",
        xlabel="x (m)", ylabel="setup (m)", ax_kw...)
    hlines!(ax2, [0.0]; color=:gray70, linewidth=0.8)
    lines!(ax2, ft_x_wave[1:n_w], ft_setup[1:n_w];
           color=:steelblue, linewidth=2.5, label="FORTRAN")
    lines!(ax2, x_jl[1:n_w], setup_jl[1:n_w];
           color=:orangered, linewidth=2, linestyle=:dash, label="Julia")
    axislegend(ax2; position=:lt, labelsize=10, framevisible=false)

    # ── Panel C: Bed-level change ─────────────────────────────────────────────
    ax3 = Axis(fig[row, 3];
        title="$(case.name) — Δzb [RMS=$(round(rms_z,digits=3)) m]",
        xlabel="x (m)", ylabel="Δzb (m)", ax_kw...)
    hlines!(ax3, [0.0]; color=:gray70, linewidth=0.8)
    lines!(ax3, ft_x_prof[1:n_p], ft_dzb[1:n_p];
           color=:steelblue, linewidth=2.5, label="FORTRAN")
    lines!(ax3, x_jl[1:n_p], dzb_jl[1:n_p];
           color=:orangered, linewidth=2, linestyle=:dash, label="Julia")
    axislegend(ax3; position=:lt, labelsize=10, framevisible=false)

    # ── Panel D: Bed profile ──────────────────────────────────────────────────
    ax4 = Axis(fig[row, 4];
        title="$(case.name) — Bed profile",
        xlabel="x (m)", ylabel="z (m)", ax_kw...)
    lines!(ax4, x0[1:nb0], zb0[1:nb0];
           color=:gray55, linewidth=1.2, linestyle=:dot, label="Initial")
    lines!(ax4, ft_x_prof[1:n_p], ft_zb[1:n_p];
           color=:steelblue, linewidth=2.5, label="FORTRAN final")
    lines!(ax4, x_jl[1:n_p], zb_jl[1:n_p];
           color=:orangered, linewidth=2, linestyle=:dash, label="Julia final")
    hlines!(ax4, [case.swl]; color=:dodgerblue, linewidth=0.8, linestyle=:dot)
    axislegend(ax4; position=:lt, labelsize=10, framevisible=false)
end

rowgap!(fig.layout, 8)
colgap!(fig.layout, 8)
save(outpath, fig; px_per_unit=2)
println("Saved → $outpath")
