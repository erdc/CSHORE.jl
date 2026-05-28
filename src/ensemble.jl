# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
ensemble.jl — Parallel scenario / parameter sweep runner.

Runs a population of CSHORE simulations from a vector of CshoreConfigs (or a
base config + patch functions) using Julia threads. Designed for uncertainty
quantification, parameter sensitivity sweeps, and what-if scenario studies.

Each member of the ensemble is fully independent — no shared mutable state.
A failure in one member is captured as the corresponding result entry and
does not interrupt the rest of the batch.
==============================================================================#

"""
    EnsembleResult

Outcome of an ensemble sweep. `results[i]` is either the final `CshoreState`
of member `i`, or the `Exception` that caused it to fail. `configs[i]` is the
config that was actually run (after patches were applied), kept so that
post-hoc analysis can correlate parameters with outcomes. `succeeded` is a
BitVector for quick `count(r.succeeded)` style queries.
"""
struct EnsembleResult
    results::Vector{Any}                # CshoreState or Exception per member
    configs::Vector{CshoreConfig}
    succeeded::BitVector
    runtimes_s::Vector{Float64}
end

Base.length(r::EnsembleResult) = length(r.results)
Base.iterate(r::EnsembleResult, state=1) =
    state > length(r) ? nothing : ((r.results[state], r.configs[state]), state + 1)

"""
    successes(r::EnsembleResult) -> Vector{CshoreState}

Final states of members that ran to completion.
"""
successes(r::EnsembleResult) = [s for (s, ok) in zip(r.results, r.succeeded) if ok]

"""
    failures(r::EnsembleResult) -> Vector{Tuple{Int,Exception}}

`(member_index, exception)` pairs for members that errored.
"""
failures(r::EnsembleResult) =
    [(i, r.results[i]) for i in eachindex(r.results) if !r.succeeded[i]]

"""
    ensemble_run(configs::Vector{CshoreConfig};
                 outdir=".", filename_pattern=nothing,
                 threaded=true, max_concurrent=Threads.nthreads(),
                 run_kwargs...) -> EnsembleResult

Run each `configs[i]` as an independent CSHORE simulation. Returns an
`EnsembleResult` whose `results[i]` is either the final state or the
exception raised by that member.

# Arguments
- `configs`: vector of fully-built `CshoreConfig`s. Use the patch-based
  method below if you want to vary a base config along one or more axes.
- `outdir`: directory for per-member NetCDF outputs. Created if missing.
- `filename_pattern`: if non-`nothing`, must contain `{i}` (replaced with
  the 1-based member index). When `nothing`, no NetCDF is written.
- `threaded`: spawn members across `Threads.@spawn` tasks. Disable for
  reproducible serial runs or when debugging.
- `max_concurrent`: cap on simultaneously-running tasks. Default is
  `Threads.nthreads()`. Set to 1 to force serial execution.
- `run_kwargs...`: forwarded to `run_simulation!` (e.g. `provenance=…`,
  `output_interval_s=…`). The same kwargs apply to every member.

# Example
```julia
configs = [cfg_dxz(d50) for d50 in (0.15e-3, 0.20e-3, 0.30e-3, 0.50e-3)]
result  = ensemble_run(configs; outdir="ensembles/d50_sweep",
                       filename_pattern="d50_{i}.nc")
@info "succeeded: \$(count(result.succeeded))/\$(length(result))"
```
"""
function ensemble_run(configs::Vector{CshoreConfig};
                      outdir::AbstractString=".",
                      filename_pattern::Union{Nothing,AbstractString}=nothing,
                      threaded::Bool=true,
                      max_concurrent::Int=Threads.nthreads(),
                      run_kwargs...)
    n = length(configs)
    n == 0 && return EnsembleResult(Any[], CshoreConfig[], BitVector(), Float64[])

    if filename_pattern !== nothing
        occursin("{i}", filename_pattern) ||
            throw(ArgumentError("ensemble_run: filename_pattern must contain '{i}'"))
        isdir(outdir) || mkpath(outdir)
    end

    results   = Vector{Any}(undef, n)
    succeeded = falses(n)
    runtimes  = zeros(Float64, n)

    # Bound concurrency with a semaphore so very large ensembles don't
    # oversubscribe the BLAS threadpool that NetCDF writes rely on.
    sem = Base.Semaphore(max(1, max_concurrent))

    run_one = function (i)
        Base.acquire(sem)
        t0 = time()
        try
            outfile = filename_pattern === nothing ? nothing :
                      replace(filename_pattern, "{i}" => string(i))
            state = run_simulation!(configs[i]; outdir=outdir, outfile=outfile,
                                    run_kwargs...)
            results[i]   = state
            succeeded[i] = true
        catch err
            results[i]   = err
            succeeded[i] = false
        finally
            runtimes[i] = time() - t0
            Base.release(sem)
        end
    end

    if threaded && Threads.nthreads() > 1 && max_concurrent > 1
        tasks = [Threads.@spawn run_one(i) for i in 1:n]
        for t in tasks
            wait(t)
        end
    else
        for i in 1:n
            run_one(i)
        end
    end

    return EnsembleResult(results, configs, succeeded, runtimes)
end

"""
    ensemble_run(base::CshoreConfig, patches::Vector;
                 kwargs...) -> EnsembleResult

Convenience form: derive members from a `base` config by applying each
`patches[i]` to a deep copy. `patches[i]` may be either a function
`(cfg::CshoreConfig) -> CshoreConfig` (must return a config; in-place
mutation is fine but the returned value is what's actually run) or a
`NamedTuple` of field overrides applied via `Setfield.@set` semantics —
but to keep dependencies light we only support the function form for now.

# Example
```julia
sweep = ensemble_run(base_cfg, [cfg -> (cfg.sediment.d50 = d50; cfg) for d50 in d50s])
```
"""
function ensemble_run(base::CshoreConfig, patches::AbstractVector;
                      kwargs...)
    # Apply each patch defensively: a patch that throws becomes a captured
    # failure for that member, not a hard error that aborts the batch. The
    # placeholder config (a deepcopy of `base`) is kept in `configs` so the
    # rest of the batch still has a valid CshoreConfig vector to iterate.
    n = length(patches)
    configs = Vector{CshoreConfig}(undef, n)
    patch_errors = Dict{Int,Exception}()
    for (i, p) in enumerate(patches)
        if !(p isa Function)
            throw(ArgumentError(
                "ensemble_run: patches[$i] must be a Function (cfg -> cfg); " *
                "got $(typeof(p))"))
        end
        try
            configs[i] = p(deepcopy(base))::CshoreConfig
        catch err
            configs[i]      = deepcopy(base)
            patch_errors[i] = err
        end
    end
    r = ensemble_run(configs; kwargs...)
    # Overlay patch-time failures on the run-time result.
    for (i, err) in patch_errors
        r.results[i]   = err
        r.succeeded[i] = false
    end
    return r
end
