# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
fractions.jl — Per-fraction grain-size helpers.
==============================================================================#

"""
    fall_velocity(d, config) -> Float64

Soulsby (1997) formula for sediment fall velocity in water.

Inputs:
- `d`       : grain diameter (m)
- `config`  : CshoreConfig — used for `sediment.sg` (specific gravity) and
              a hard-coded kinematic viscosity ν = 1e-6 m²/s (20°C water).

Returns settling velocity (m/s).
- Stokes drag for d_star < 1
- Soulsby transitional formula for 1 ≤ d_star < 100
- Turbulent drag for d_star ≥ 100
"""
function fall_velocity(d::Float64, config::CshoreConfig)
    ν = 1.0e-6                     # kinematic viscosity of water at 20°C
    sgm1 = submerged_sgm1(config.sediment)   # ρ_s / ρ_w - 1 (rho_water-aware)
    dstar = d * (sgm1 * GRAV / ν^2)^(1 / 3)
    if dstar < 1.0
        return sgm1 * GRAV * d^2 / (18 * ν)                         # Stokes
    elseif dstar < 100.0
        return ν / d * (sqrt(10.36^2 + 1.049 * dstar^3) - 10.36)    # Soulsby transitional
    else
        return 1.05 * sqrt(sgm1 * GRAV * d)                         # Turbulent
    end
end

"""
    critical_shields(d, config) -> Float64

Soulsby (1997) critical Shields parameter θ_cr as a function of grain size.

Inputs:
- `d`       : grain diameter (m)
- `config`  : CshoreConfig

If `multifraction.use_size_dependent_shields` is `false`, returns the constant
`config.sediment.shield` (FORTRAN parity — the original code uses a single
θ_cr = 0.05 regardless of grain size).

Otherwise uses the Soulsby 1997 fit:
    θ_cr = 0.30 / (1 + 1.2·d_*) + 0.055 · (1 − exp(−0.020·d_*))

where d_* = d · ((s−1)g/ν²)^(1/3).
"""
function critical_shields(d::Float64, config::CshoreConfig)
    if !config.multifraction.use_size_dependent_shields
        return config.sediment.shield
    end
    sgm1 = submerged_sgm1(config.sediment)
    ν = 1.0e-6
    dstar = d * (sgm1 * GRAV / ν^2)^(1 / 3)
    return 0.30 / (1 + 1.2 * dstar) + 0.055 * (1 - exp(-0.020 * dstar))
end

"""
    hiding_exposure_factor(d_k, d_mean, config) -> Float64

Hiding-exposure correction for mixed-size
beds. Returns a multiplier on θ_cr for fraction `k`.

Disabled by default when `use_size_dependent_shields=true` to avoid
double-counting size dependence (the Shields curve already has it).

Fine grains (d_k/d_mean < 0.4) hide below coarse ones and are harder to move;
coarse grains are more exposed and easier to move — but the sign convention
here follows the Python port, where ξ_k > 1 means harder to move (increases
θ_cr), ξ_k < 1 means easier.
"""
function hiding_exposure_factor(d_k::Float64, d_mean::Float64, config::CshoreConfig)
    mf = config.multifraction
    if !mf.use_hiding_exposure || mf.use_size_dependent_shields
        return 1.0
    end
    ratio = d_k / d_mean
    if mf.hiding_method === :egiazaroff
        return ratio < 0.4 ? 0.843 * ratio^(-0.5) : ratio^(-1.0)
    elseif mf.hiding_method === :ashida_michiue
        return ratio^(-0.6)
    else
        throw(ArgumentError("Unknown hiding_method $(mf.hiding_method)"))
    end
end

"""
    grainsize_tadapt(d, h, config) -> Float64

Per-fraction suspended-sediment adaptation time.

    T_adapt_k = h / w_s(d_k)

clipped to [5, 500] seconds and scaled by `tadapt_multiplier`.
"""
function grainsize_tadapt(d::Float64, h::Float64, config::CshoreConfig)
    ws = fall_velocity(d, config)
    t = h / (ws + 1e-10)
    t *= config.multifraction.tadapt_multiplier
    return clamp(t, 5.0, 500.0)
end

"""
    surface_d50(state, j) -> Float64

Mean grain diameter of the surface bed layer at node `j`, weighted by mass
fraction. Used in hiding-exposure and as a diagnostic for grain sorting.
"""
function surface_d50(state::CshoreState, config::CshoreConfig, j::Int)
    nf = size(state.bed_mass, 3)
    nf == 0 && return config.sediment.d50
    total = 0.0
    weighted = 0.0
    @inbounds for k in 1:nf
        m = state.bed_mass[j, 1, k]
        total += m
        weighted += m * config.multifraction.grain_sizes[k]
    end
    return total > 0 ? weighted / total : config.sediment.d50
end
