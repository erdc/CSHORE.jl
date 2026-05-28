#==============================================================================
plot_benchmarks_plots.jl — Plots.jl version of plot_benchmarks.jl.

Same content as plot_benchmarks.jl (which depends on CairoMakie that isn't
in the default Project.toml) but using Plots.jl which is already a dep.

Usage:
  julia --project=. examples/benchmarks/erdc_usace/plot_benchmarks_plots.jl [outfile]

Output: a 4-column × 8-row PNG grid, one row per ERDC case.
Columns: (a) Hrms   (b) wave setup   (c) bed-level change   (d) bed profile
==============================================================================#

include(joinpath(@__DIR__, "compare_utils.jl"))
using Plots

outpath = length(ARGS) >= 1 ? ARGS[1] :
          joinpath(@__DIR__, "erdc_comparison.png")

ncases = length(CASES)
panels = Vector{Plots.Plot}(undef, 4 * ncases)

for (row, case) in enumerate(CASES)
    ref_os = read_ref_block(joinpath(case.dir, "ref_OSETUP.txt"); ncols=4)
    ref_ob = read_ref_block(joinpath(case.dir, "ref_OBPROF.txt"); ncols=2)

    ft_x_wave = ref_os[:, 1]
    ft_hrms   = sqrt(8.0) .* ref_os[:, 4]
    ft_setup  = ref_os[:, 2] .- case.swl
    ft_x_prof = ref_ob[:, 1]
    ft_zb     = ref_ob[:, 2]

    cfg = read_infile(joinpath(case.dir, "infile"); strict=false)
    nb0  = cfg.bathymetry.nbinp[1]
    x0   = cfg.bathymetry.xbinp[1:nb0, 1]
    zb0  = cfg.bathymetry.zbinp[1:nb0, 1]
    ft_dzb = ft_zb .- zb0[1:length(ft_zb)]

    state = run_simulation!(cfg)
    jmax = state.jmax[1]; jr = state.jr
    x_jl     = state.xb[1:jmax]
    hrms_jl  = state.hrms[1:jmax]
    setup_jl = state.wsetup[1:jmax]
    zb_jl    = state.zb[1:jmax, 1]
    dzb_jl   = zb_jl .- zb0[1:jmax]

    n_w = min(length(ft_x_wave), jr)
    n_p = min(length(ft_x_prof), jmax)

    rms_h = sqrt(sum((ft_hrms[1:n_w] .- hrms_jl[1:n_w]).^2) / n_w)
    rms_z = sqrt(sum((ft_dzb[1:n_p]  .- dzb_jl[1:n_p]).^2)  / n_p)
    rms_s = sqrt(sum((ft_setup[1:n_w] .- setup_jl[1:n_w]).^2) / n_w)

    # Panel A — Hrms
    p1 = plot(ft_x_wave[1:n_w], ft_hrms[1:n_w]; lw=2.5, color=:steelblue,
        label="FORTRAN", xlabel="x (m)", ylabel="Hrms (m)",
        title="$(case.name) — Hrms\n[RMS=$(round(rms_h,digits=3)) m]",
        titlefontsize=9, legend=:topleft, legendfontsize=7)
    plot!(p1, x_jl[1:n_w], hrms_jl[1:n_w]; lw=2, color=:orangered,
        linestyle=:dash, label="Julia")

    # Panel B — Wave setup
    p2 = plot(ft_x_wave[1:n_w], ft_setup[1:n_w]; lw=2.5, color=:steelblue,
        label="FORTRAN", xlabel="x (m)", ylabel="setup (m)",
        title="$(case.name) — Setup\n[RMS=$(round(rms_s,digits=3)) m]",
        titlefontsize=9, legend=:topleft, legendfontsize=7)
    plot!(p2, x_jl[1:n_w], setup_jl[1:n_w]; lw=2, color=:orangered,
        linestyle=:dash, label="Julia")
    hline!(p2, [0.0]; color=:gray70, lw=0.6, label="")

    # Panel C — Δzb
    p3 = plot(ft_x_prof[1:n_p], ft_dzb[1:n_p]; lw=2.5, color=:steelblue,
        label="FORTRAN", xlabel="x (m)", ylabel="Δzb (m)",
        title="$(case.name) — Δzb\n[RMS=$(round(rms_z,digits=3)) m]",
        titlefontsize=9, legend=:topleft, legendfontsize=7)
    plot!(p3, x_jl[1:n_p], dzb_jl[1:n_p]; lw=2, color=:orangered,
        linestyle=:dash, label="Julia")
    hline!(p3, [0.0]; color=:gray70, lw=0.6, label="")

    # Panel D — Bed profile
    p4 = plot(x0[1:nb0], zb0[1:nb0]; lw=1.2, color=:gray55, linestyle=:dot,
        label="Initial", xlabel="x (m)", ylabel="z (m)",
        title="$(case.name) — Profile",
        titlefontsize=9, legend=:topleft, legendfontsize=7)
    plot!(p4, ft_x_prof[1:n_p], ft_zb[1:n_p]; lw=2.5, color=:steelblue,
        label="FORTRAN final")
    plot!(p4, x_jl[1:n_p], zb_jl[1:n_p]; lw=2, color=:orangered,
        linestyle=:dash, label="Julia final")
    hline!(p4, [case.swl]; color=:dodgerblue, lw=0.6, linestyle=:dot, label="")

    panels[4 * (row - 1) + 1] = p1
    panels[4 * (row - 1) + 2] = p2
    panels[4 * (row - 1) + 3] = p3
    panels[4 * (row - 1) + 4] = p4
end

big = plot(panels...; layout=(ncases, 4), size=(1600, 320 * ncases),
    plot_title="Julia CSHORE.jl  vs  FORTRAN CSHORE_USACE — ERDC testbed",
    plot_titlefontsize=14, left_margin=4Plots.mm, bottom_margin=2Plots.mm,
    dpi=110)
savefig(big, outpath)
println("Saved → $outpath")
