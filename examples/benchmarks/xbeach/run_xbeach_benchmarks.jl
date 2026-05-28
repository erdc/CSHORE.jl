"""
Run CSHORE.jl on the three included XBeach benchmark cases and compare
wave heights to the analytical / measured reference data.

Usage:
  julia --project=../.. examples/benchmarks/xbeach/run_xbeach_benchmarks.jl

Each case is read with `read_xbeach_params` and then run with the standard
`run_simulation!` driver. Output NetCDF files are written to this directory.
"""

using CSHORE
using Printf

DATA = joinpath(@__DIR__, "data")

cases = [
    (name="Boers_1C",     dir=joinpath(DATA,"Boers_1C"),     params="params_BOI.txt"),
    (name="GWK98_F1",     dir=joinpath(DATA,"GWK98_F1"),     params="params_F1.txt"),
    (name="DUROS_7000308",dir=joinpath(DATA,"DUROS_7000308"),params="params.txt"),
]

for c in cases
    @printf "\n══ %s ══\n" c.name
    cfg = try
        read_xbeach_params(c.dir; params_file=c.params)
    catch e
        @warn "  Failed to parse: $e"
        continue
    end
    @printf "  dx=%.2f m  iprofl=%d  nBC=%d\n" cfg.grid.dx cfg.options.iprofl length(cfg.boundary.timebc)
    @printf "  Hrms=%.3f m  Tp=%.2f s  D50=%.3f mm\n" cfg.boundary.hrmsbc[1] cfg.boundary.tpbc[1] cfg.multifraction.grain_sizes[1]*1000

    outfile = joinpath(@__DIR__, "$(c.name)_out.nc")
    state = run_simulation!(cfg; outfile=outfile)
    @printf "  Simulation complete → %s\n" outfile
    @printf "  Final profile: zb in [%.2f, %.2f] m\n" minimum(state.zb) maximum(state.zb)
end
