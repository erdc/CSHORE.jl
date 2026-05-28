# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
hydro.jl — Hydrodynamic helpers.
==============================================================================#

"""
    DerivedConstants

Scalars computed once by `compute_derived_constants!` from `CshoreConfig`.

This struct is a diagnostic container returned by `compute_derived_constants!`
for unit tests and debugging.
"""
Base.@kwdef struct DerivedConstants
    wkpo::Float64            # deep-water wavenumber = (2π/T)²/g
    sqrg1_pi::Float64        # √(2/π) — Gaussian-friction constant
    sqrg2_pi::Float64        # 2·√(2/π)
    # Porous parameters (zero unless iperm=1 or infilt=1)
    alpha::Float64 = 0.0
    beta1::Float64 = 0.0
    beta2::Float64 = 0.0
    alsta::Float64 = 0.0
    besta1::Float64 = 0.0
    besta2::Float64 = 0.0
    alsta2::Float64 = 0.0
    be2::Float64 = 0.0
    be4::Float64 = 0.0
    wpm::Float64 = 0.0
    # Swash parameters (zero unless iover=1)
    awd::Float64 = 0.0
    wdn::Float64 = 0.0
    ewd::Float64 = 0.0
    cwd::Float64 = 0.0
    aqwd::Float64 = 0.0
    bwd::Float64 = 0.0
    agwd::Float64 = 0.0
    auwd::Float64 = 0.0
end

"""
    compute_derived_constants!(state, config) -> DerivedConstants

- Sets `state.rbzero = 0.1` (roller slope).
- Computes the deep-water wavenumber `wkpo = (2π/Tp)²/g` using the first
  wave period in the boundary time series.
- If `iperm=1` or `infilt=1`: porous-flow constants (α, β1, β2, WPM).
- If `iover=1`: swash-zone constants (AWD, EWD, CWD, etc.).

The returned `DerivedConstants` struct is a diagnostic — the main data are
stored on `state` or computed inline as needed.
"""
function compute_derived_constants!(state::CshoreState, config::CshoreConfig)
    state.rbzero = 0.1
    tp0 = first(config.boundary.tpbc)
    tp0 > 0 || throw(ArgumentError("First boundary period must be > 0, got $tp0"))
    wkpo = (PI2 / tp0)^2 / GRAV

    dc_kwargs = Dict{Symbol,Float64}(
        :wkpo => wkpo,
        :sqrg1_pi => sqrt(2.0 / π),
        :sqrg2_pi => 2.0 * sqrt(2.0 / π),
    )

    if config.options.iperm == 1 || config.options.infilt == 1
        p = config.porous
        wnu = p === nothing ? 1.0e-6 : p.wnu
        A = 1000.0
        B = 5.0
        dump_p = if config.options.iperm == 1
            p === nothing ? 0.4 : p.snp
        else
            1.0 - config.sediment.sporo1
        end
        dumd = if config.options.iperm == 1
            p === nothing ? 0.05 : p.sdp
        else
            config.sediment.d50
        end
        C = 1.0 - dump_p
        α  = A * wnu * C^2 / (dump_p * dumd)^2
        β1 = B * C / dump_p^3 / dumd
        β2 = 7.5 * B * C / sqrt(2.0) / dump_p^2
        alsta  = α  / GRAV
        besta1 = β1 / GRAV
        besta2 = β2 / GRAV
        alsta2 = alsta * alsta
        be2    = 2.0 * besta1
        be4    = 2.0 * be2
        wpm    = (sqrt(alsta2 + be4) - alsta) / be2
        merge!(dc_kwargs, Dict(:alpha=>α, :beta1=>β1, :beta2=>β2, :alsta=>alsta,
                                :besta1=>besta1, :besta2=>besta2, :alsta2=>alsta2,
                                :be2=>be2, :be4=>be4, :wpm=>wpm))
    end

    if config.options.iover == 1
        awd = config.options.iprofl == 1 && config.options.iperm == 0 ? 1.6 : 2.0
        ewd = config.options.iperm == 1 ? 0.01 : 0.015
        cwd = 0.75 * sqrt(π)
        aqwd = cwd * awd
        agwd = awd * awd
        auwd = 0.5 * sqrt(π) * awd
        bwd = (2.0 - 9.0 * π / 16.0) * agwd + 1.0
        merge!(dc_kwargs, Dict(:awd=>awd, :ewd=>ewd, :cwd=>cwd, :aqwd=>aqwd,
                                :agwd=>agwd, :auwd=>auwd, :bwd=>bwd))
        # Wire slpot from sediment config into state — this is the overtopping-
        # driven onshore suspended transport coefficient. Was previously always
        # 0.0 when using build_config programmatically.
        state.slpot = config.sediment.slpot
    end

    return DerivedConstants(; dc_kwargs...)
end

"""
    compute_bed_slope!(state, config, l)

Populate `state.bslope[j, l] = -dzb/dx` (note the sign convention — positive
bslope means bed *rises* shoreward). Uses central differences for interior
nodes and one-sided at boundaries.

Called once per line before the landward march starts.
"""
function compute_bed_slope!(state::CshoreState, config::CshoreConfig, l::Int)
    dx = config.grid.dx
    n = state.jmax[l]
    zb = @view state.zb[:, l]
    @inbounds begin
        # Forward diff at seaward node
        state.bslope[1, l] = (zb[2] - zb[1]) / dx
        # Central diff in interior
        for j in 2:n-1
            state.bslope[j, l] = (zb[j+1] - zb[j-1]) / (2 * dx)
        end
        # Backward diff at landward node
        state.bslope[n, l] = (zb[n] - zb[n-1]) / dx
    end
    return state
end
