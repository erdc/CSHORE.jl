# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
vegetation.jl — Vegetation drag and dissipation (IVEG=1,2,3).
==============================================================================#

"""
    veg_dissipation!(state, config, j, l)

Compute vegetation-induced wave energy dissipation rate `dvegsta[j]`.
Implements IDISS=1: Mendez & Losada (2004).
"""
function veg_dissipation!(state::CshoreState, config::CshoreConfig, j::Int, l::Int)
    veg = config.vegetation::VegetationInput
    d = state.h[j]

    # No dissipation if no vegetation or no water
    if veg.vegn[j, l] <= 0.0 || d <= 0.0 || state.hrms[j] <= 0.0
        state.dvegsta[j] = 0.0
        return nothing
    end

    # Effective vegetation height (submerged: full height; emergent: capped at depth).
    # Use dynamic vegh from state when uprooting is tracked (IVEG=1,3 + IPROFL=1);
    # otherwise fall back to the static config height.
    vh = size(state.vegh, 1) > 0 ? state.vegh[j, l] : veg.vegd[j, l]
    effvegh = min(vh, d)
    if effvegh <= 0.0
        state.dvegsta[j] = 0.0
        return nothing
    end

    k = state.wkp          # peak wavenumber (from lwave!)
    hrms = state.hrms[j]
    t = state.wt[j]        # local wave period

    # Mendez & Losada (2004) dissipation formula:
    #   Dv = (0.5/√π/g) · Cd · bv · Nv · (k·g·T/4π)³ ·
    #        (sinh³(k·αh) + 3·sinh(k·αh)) · Hrms³ / (3k·cosh³(k·h))
    kd = k * d
    kah = k * effvegh
    sinh_kah = sinh(kah)
    cosh_kd = cosh(kd)

    # Guard against very large arguments (cosh overflow)
    if cosh_kd > 1e100
        state.dvegsta[j] = 0.0
        return nothing
    end

    # Effective drag coefficient: optionally reduce Cd based on the
    # Keulegan-Carpenter number KC = u_max * T / bv, following the empirical
    # fit of He et al. (2019, HESS 25:4825):
    #   Cd_eff = Cd0 * exp(-0.043 * KC)  [+ baseline 1.2 absorbed into Cd0]
    # For KC < 5 (small waves relative to stem), Cd_eff ≈ Cd0.
    # For KC > 50 (large waves), Cd_eff drops to ~10-20% of Cd0.
    # This accounts for flow separation, vortex shedding regime changes,
    # and — for flexible vegetation — effective blade reconfiguration
    # (Luhar & Nepf 2016, JGR Oceans 121).
    #
    # When Cd0 is already low (≤0.5, typical for reef presets), the KC
    # correction is small and mostly a refinement. When Cd0 is high
    # (≥1.0, rigid stems), the correction provides significant reduction
    # at high wave heights.
    bv = veg.vegb[j, l]
    u_max = π * hrms / (t * sinh(kd))  # near-bed orbital velocity amplitude
    # Cap KC at 50 in the He et al. (2019) reduction. Beyond KC≈50 the
    # empirical fit flattens out at ~10-20% of Cd0; continuing to apply
    # exp(-0.043·KC) for KC = O(10³) (which happens for thin blades, e.g.
    # bv = 5 mm with O(1 m/s) orbital velocity) would unphysically drive
    # Cd_eff to zero and silently kill all wave-vegetation dissipation.
    KC_raw = bv > 0 ? u_max * t / bv : 0.0
    KC = min(KC_raw, 50.0)
    cd_eff = KC > 1.0 ? veg.vegcd * exp(-0.043 * KC) : veg.vegcd

    state.dvegsta[j] = (0.5 / sqrt(π) / GRAV) * cd_eff * bv * veg.vegn[j, l] *
                       (0.25 * k * GRAV * t / π)^3 *
                       (sinh_kah^3 + 3.0 * sinh_kah) *
                       hrms^3 /
                       (3.0 * k * cosh_kd^3)

    return nothing
end

"""
    apply_veg_friction!(state, config, j, l)

Apply vegetation friction enhancement for IVEG=1 or IVEG=2.
Multiplies `tbxsta[j]` and `dfsta[j]` by `VEGCV = 1 + min(vegd, h) * vegfb`.
"""
function apply_veg_friction!(state::CshoreState, config::CshoreConfig, j::Int, l::Int)
    veg = config.vegetation::VegetationInput
    if veg.vegn[j, l] <= 0.0
        return nothing
    end
    # Use dynamic vegh from state when uprooting is tracked, else config.
    vh = size(state.vegh, 1) > 0 ? state.vegh[j, l] : veg.vegh[j, l]
    dumh = min(vh, state.h[j])
    if dumh <= 0.0
        return nothing
    end
    vegcv = 1.0 + dumh * veg.vegfb[j, l]
    state.tbxsta[j] *= vegcv
    state.dfsta[j] *= vegcv
    return nothing
end

"""
    veg_momentum_stress(state, config, j, l) -> Float64

Compute the vegetation-enhanced stream stress multiplier for the momentum
equation when IVEG=3. Returns the multiplier for TBXSTA at node j:

    (1 + VEGCDM/VEGCD * min(VEGH, h) * VEGFB) * TBXSTA
"""
function veg_momentum_multiplier(state::CshoreState, config::CshoreConfig, j::Int, l::Int)
    veg = config.vegetation::VegetationInput
    if veg.vegn[j, l] <= 0.0
        return 1.0
    end
    # Use dynamic vegh from state when uprooting is tracked; else fallback to config.
    vh = size(state.vegh, 1) > 0 ? state.vegh[j, l] : veg.vegh[j, l]
    dumh = min(vh, state.h[j])
    return 1.0 + (veg.vegcdm / veg.vegcd) * dumh * veg.vegfb[j, l]
end

"""
    initialize_vegetation_bounds!(state, config)

Sets the fixed elevation bounds for vegetation at the start of a simulation.
The bounds do not change even as the bed evolves — they anchor the window in
which the vegetation can remain active. The dynamic canopy height `state.vegh`
is initialized from the config input and then updated by
`check_vegetation_uprooting!` after every morphodynamic sub-step.

Bounds:
- `vegzd = zb0 + vegh0` — upper elevation (top of canopy)
- `vegzr = zb0 - vegrh` — lower elevation (bottom of root zone)
- `uproot = 1.0` if vegetation present (vegfb > 0), else 0.0

Called once from `run_simulation!` / `step_bc_window!` initialization.
Only active when `iveg ∈ {1, 3}` and `iprofl == 1`.
"""
function initialize_vegetation_bounds!(state::CshoreState, config::CshoreConfig)
    iveg = config.options.iveg
    iprofl = config.options.iprofl
    (iveg in (1, 3) && iprofl == 1) || return nothing
    config.vegetation === nothing && return nothing

    veg = config.vegetation::VegetationInput
    nl = config.options.iline
    @inbounds for l in 1:nl
        jmax_l = state.jmax[l]
        for j in 1:jmax_l
            vegh0 = veg.vegh[j, l]
            vegrh = veg.vegrh[j, l]
            zb0 = state.zb[j, l]
            state.vegh[j, l] = vegh0          # dynamic height starts equal to static
            state.vegzd[j, l] = zb0 + vegh0    # top of canopy
            state.vegzr[j, l] = zb0 - vegrh    # bottom of roots
            state.uproot[j, l] = veg.vegfb[j, l] > 0.0 ? 1.0 : 0.0
        end
    end
    return nothing
end

"""
    check_vegetation_uprooting!(state, config, l)

After each morphodynamic sub-step, recompute the effective vegetation
height `state.vegh[j, l]` based on the current bed elevation `state.zb[j, l]`
relative to the fixed bounds `vegzd` and `vegzr`. Two failure modes:
"""
function check_vegetation_uprooting!(state::CshoreState, config::CshoreConfig, l::Int)
    iveg = config.options.iveg
    iprofl = config.options.iprofl
    (iveg in (1, 3) && iprofl == 1) || return nothing
    size(state.uproot, 1) == 0 && return nothing

    jmax_l = state.jmax[l]
    @inbounds for j in 1:jmax_l
        if state.uproot[j, l] == 1.0
            # Burial check: bed has risen above canopy → clamp height
            vh = state.vegzd[j, l] - state.zb[j, l]
            if vh < 0.0
                vh = 0.0
            end
            state.vegh[j, l] = vh
            # Uprooting check: erosion into the root zone
            if state.vegzr[j, l] >= state.zb[j, l]
                state.uproot[j, l] = 0.0
                state.vegh[j, l] = 0.0
            end
        end
    end
    return nothing
end
