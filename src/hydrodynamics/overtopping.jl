# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
overtopping.jl — Wave overtopping / runup / landward transmission.
==============================================================================#

"""
    overtopping_rate!(state, config, itime, l, iteqo) -> (converged::Bool)

Iterative convergence loop for the overtopping flux `QO(L)`:

1. Calls `wetdry!` to compute swash zone and overtopping flux `QOTF`.
2. Computes new combined rate `QONEW = QOTF`.
3. Checks convergence of `QO(L)`.
4. Blends old and new values with a damping factor.

Returns `true` if converged (no more iterations needed), `false` if not.
"""
function overtopping_rate!(state::CshoreState, config::CshoreConfig,
                            itime::Int, l::Int, iteqo::Int)
    wetdry!(state, config, itime, l, iteqo)

    qonew = state.qotf

    if qonew < 1e-5
        state.qo[l] = qonew
        return true
    end

    dum = abs(qonew - state.qo[l]) / qonew
    aer = 1e-4 / qonew
    if aer < 1e-2; aer = 1e-2; end
    if dum < aer
        state.qo[l] = qonew
        return true
    end

    # Blend with increasing weight on the new value as iterations progress
    fractn = 0.5 + 0.1 * iteqo
    if fractn > 0.9; fractn = 0.9; end
    if iteqo == 10; fractn = 0.5; end
    if iteqo == 20; fractn = 0.0; end
    state.qo[l] = fractn * state.qo[l] + (1.0 - fractn) * qonew

    return false
end

# Wave transmission to a landward water body (IWTRAN=1) is implemented in
# hydrodynamics/transmission.jl.
