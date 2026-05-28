#==============================================================================
compare_utils.jl — Shared utilities for ERDC USACE benchmark comparisons.

Provides:
  - read_ref_block     : parse a single-block FORTRAN reference file
  - BenchmarkCase      : case descriptor (paths, SWL at final time)
  - CASES              : the 8 standard testbed cases
  - run_and_compare    : run Julia, diff against reference, return metrics
==============================================================================#

using CSHORE
using Statistics

# ── Single-block reference file parser ────────────────────────────────────────
"""
    read_ref_block(path; ncols) -> Matrix{Float64}

Read a single FORTRAN block-format reference file (header line + NJ data rows)
and return an NJ × ncols matrix.  The header line `L NJ TIME` is skipped.
"""
function read_ref_block(path::AbstractString; ncols::Int)
    lines = readlines(path)
    isempty(lines) && error("Empty reference file: $path")
    header = split(strip(lines[1]))
    length(header) == 3 || error("Expected header 'L NJ TIME' on line 1 of $path")
    nj = parse(Int, header[2])
    length(lines) == nj + 1 || error(
        "Expected $(nj+1) lines in $path, got $(length(lines))")
    data = Matrix{Float64}(undef, nj, ncols)
    for r in 1:nj
        toks = split(strip(lines[r + 1]))
        for c in 1:ncols
            data[r, c] = parse(Float64, toks[c])
        end
    end
    return data
end

# ── Case descriptors ───────────────────────────────────────────────────────────
const CASES_DIR = joinpath(@__DIR__, "cases")

struct BenchmarkCase
    name    ::String   # human-readable label
    dir     ::String   # path to cases/<slug>/
    swl     ::Float64  # still-water level at the final time step (m)
end

"""
The 8 standard ERDC USACE testbed benchmark cases with pre-computed
FORTRAN CSHORE_USACE reference outputs.
"""
const CASES = BenchmarkCase[
    BenchmarkCase("Agate Beach Sep-2013",    joinpath(CASES_DIR, "agate_sep"),      2.139),
    BenchmarkCase("Agate Beach Oct-2013",    joinpath(CASES_DIR, "agate_oct"),      3.260),
    BenchmarkCase("FRF BathyDuck runup #2",  joinpath(CASES_DIR, "frf_runup_2"),   -0.398),
    BenchmarkCase("FRF BathyDuck runup #3",  joinpath(CASES_DIR, "frf_runup_3"),    0.557),
    BenchmarkCase("FRF general morpho 070",  joinpath(CASES_DIR, "frf_morpho_070"), 0.0),
    BenchmarkCase("FRF general morpho 071",  joinpath(CASES_DIR, "frf_morpho_071"), 0.0),
    BenchmarkCase("GEE laboratory",          joinpath(CASES_DIR, "gee"),            0.0),
    BenchmarkCase("LSTF laboratory",         joinpath(CASES_DIR, "lstf"),           0.0),
]

# ── Per-case metrics struct ────────────────────────────────────────────────────
struct BenchmarkResult
    case      ::BenchmarkCase
    rms_hrms  ::Float64   # RMS |Hrms_julia - Hrms_fortran| (m)
    rms_setup ::Float64   # RMS |setup_julia - setup_fortran| (m)
    rms_dzb   ::Float64   # RMS |Δzb_julia - Δzb_fortran| (m)  (0 for fixed bed)
    rms_umean ::Float64   # RMS |umean_julia - umean_fortran| (m/s)
    n_wave    ::Int        # number of nodes used in wave comparison
    n_prof    ::Int        # number of nodes used in profile comparison
    ok        ::Bool       # did the simulation complete without error?
    err_msg   ::String
end

# ── Main comparison function ───────────────────────────────────────────────────
"""
    run_and_compare(case::BenchmarkCase; strict=false) -> BenchmarkResult

Parse `case.dir/infile`, run the Julia simulation, and compare against the
pre-computed reference blocks in `case.dir/ref_OSETUP.txt`, `ref_OBPROF.txt`,
and `ref_OXVELO.txt`.

Reference OSETUP columns: x, MWS (SWL+setup), h, sigma
Reference OBPROF columns: x, zb
Reference OXVELO columns: x, umean, ustd
"""
function run_and_compare(case::BenchmarkCase; strict::Bool=false)
    infile  = joinpath(case.dir, "infile")
    f_osetup = joinpath(case.dir, "ref_OSETUP.txt")
    f_obprof = joinpath(case.dir, "ref_OBPROF.txt")
    f_oxvelo = joinpath(case.dir, "ref_OXVELO.txt")

    # ── Run Julia ────────────────────────────────────────────────────────────
    cfg = try
        read_infile(infile; strict=strict)
    catch e
        return BenchmarkResult(case, NaN, NaN, NaN, NaN, 0, 0, false,
                               "parse error: $(sprint(showerror, e))")
    end

    # Capture initial bed for Δzb
    nb0   = cfg.bathymetry.nbinp[1]
    zb0   = cfg.bathymetry.zbinp[1:nb0, 1]

    state = try
        run_simulation!(cfg)
    catch e
        return BenchmarkResult(case, NaN, NaN, NaN, NaN, 0, 0, false,
                               "run error: $(sprint(showerror, e))")
    end

    jmax = state.jmax[1]
    jr   = state.jr
    x_jl    = state.xb[1:jmax]
    hrms_jl = state.hrms[1:jmax]
    setup_jl = state.wsetup[1:jmax]
    zb_jl   = state.zb[1:jmax, 1]
    umean_jl = state.umean[1:jmax]
    dzb_jl  = zb_jl .- zb0[1:jmax]

    # ── Load reference ───────────────────────────────────────────────────────
    ref_os = read_ref_block(f_osetup; ncols=4)
    ref_ob = read_ref_block(f_obprof; ncols=2)
    ref_ox = read_ref_block(f_oxvelo; ncols=3)

    ft_x_wave  = ref_os[:, 1]
    ft_hrms    = sqrt(8.0) .* ref_os[:, 4]
    ft_setup   = ref_os[:, 2] .- case.swl
    ft_x_prof  = ref_ob[:, 1]
    ft_zb      = ref_ob[:, 2]
    ft_dzb     = ft_zb .- zb0[1:length(ft_zb)]   # FORTRAN Δzb vs same initial
    ft_umean   = ref_ox[:, 2]

    # ── Compare ──────────────────────────────────────────────────────────────
    n_w = min(length(ft_x_wave), jr)
    n_p = min(length(ft_x_prof), jmax)

    Δhrms  = ft_hrms[1:n_w]  .- hrms_jl[1:n_w]
    Δsetup = ft_setup[1:n_w] .- setup_jl[1:n_w]
    Δumean = ft_umean[1:n_w] .- umean_jl[1:n_w]
    Δdzb   = ft_dzb[1:n_p]   .- dzb_jl[1:n_p]

    rms(v) = sqrt(mean(v .^ 2))

    return BenchmarkResult(case,
        rms(Δhrms), rms(Δsetup), rms(Δdzb), rms(Δumean),
        n_w, n_p, true, "")
end
