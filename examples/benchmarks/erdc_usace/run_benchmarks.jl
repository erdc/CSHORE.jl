"""
Run the ERDC USACE testbed benchmarks and print a comparison table.

Usage:
  julia --project=../../.. examples/benchmarks/erdc_usace/run_benchmarks.jl

Each case is parsed with `read_infile` and run with `run_simulation!`, then
compared against the pre-computed FORTRAN CSHORE_USACE reference outputs stored
in cases/<name>/ref_OSETUP.txt, ref_OBPROF.txt, and ref_OXVELO.txt.

RMS errors are reported for:
  Hrms   — root-mean-square wave height (m)
  setup  — wave setup above SWL (m)
  Δzb    — bed-level change (m)  [0 for fixed-bed cases]
  umean  — mean cross-shore velocity (m/s)
"""

include(joinpath(@__DIR__, "compare_utils.jl"))
using Printf

println("\nERDC USACE testbed benchmarks — Julia CSHORE.jl vs FORTRAN CSHORE_USACE")
println("=" ^ 78)
@printf "  %-26s  %5s  %8s  %8s  %8s  %8s\n" "Case" "nodes" "Hrms(m)" "setup(m)" "Δzb(m)" "umean(m/s)"
println("  " * "─"^74)

results = BenchmarkResult[]
for case in CASES
    r = run_and_compare(case)
    push!(results, r)
    if r.ok
        @printf "  %-26s  %5d  %8.4f  %8.4f  %8.4f  %8.4f\n" \
            r.case.name r.n_wave r.rms_hrms r.rms_setup r.rms_dzb r.rms_umean
    else
        @printf "  %-26s  FAILED: %s\n" r.case.name r.err_msg
    end
end

n_ok   = count(r.ok for r in results)
n_fail = length(results) - n_ok
println("  " * "─"^74)
println("  $n_ok/$(length(results)) cases passed")
n_fail > 0 && println("  !! $n_fail case(s) failed — see messages above")
println()
