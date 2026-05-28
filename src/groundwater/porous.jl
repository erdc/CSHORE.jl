# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
porous.jl — Porous-layer flow (IPERM=1).

Computes the mean and oscillatory discharge velocity, water flux, and
energy dissipation within a permeable layer (e.g., gravel beach,
rubble-mound revetment) using the Forchheimer resistance equation.

The porous layer acts as:
  - An energy sink (DPSTA added to total dissipation in the energy equation)
  - A water sink (QP reduces the effective cross-shore water flux)
  - Both effects reduce wave setup and modify swash excursion
==============================================================================#

"""
    porous_flow!(state, config, j, l, pkhsig, dedx)

Computes porous-layer velocity moments and dissipation at node `j`:
- `state.upmean[j]` — mean horizontal discharge velocity (m/s)
- `state.upstd[j]`  — std dev of discharge velocity (m/s)
- `state.qp[j]`     — water flux through porous layer (m²/s)
- `state.dpsta[j]`  — energy dissipation rate in porous layer

Arguments:
- `pkhsig`: dimensionless wave parameter `k·σ` (product of wavenumber and sigma)
- `dedx`: water-surface elevation gradient `∂η/∂x`
"""
function porous_flow!(state::CshoreState, config::CshoreConfig,
                      j::Int, l::Int, pkhsig::Float64, dedx::Float64)
    hp_j = state.hp[j, l]

    # No porous layer at this node — zero everything
    if hp_j == 0.0
        state.upmean[j] = 0.0
        state.upstd[j]  = 0.0
        state.dpsta[j]  = 0.0
        state.qp[j]     = 0.0
        return nothing
    end

    # Forchheimer resistance coefficients
    p = config.porous
    wnu = p === nothing ? 1.0e-6 : p.wnu
    snp = p === nothing ? 0.4 : p.snp
    sdp = p === nothing ? 0.05 : p.sdp
    C_por = 1.0 - snp
    α  = 1000.0 * wnu * C_por^2 / (snp * sdp)^2
    β1 = 5.0 * C_por / snp^3 / sdp
    β2 = 7.5 * 5.0 * C_por / sqrt(2.0) / snp^2
    alsta  = α  / GRAV
    besta1 = β1 / GRAV
    besta2 = β2 / GRAV

    # Standard deviation of porous velocity (Forchheimer resistance)
    a_coef = 1.9 * besta1
    b2_val = besta2 / state.wt[j]
    b_coef = alsta + 1.9 * b2_val
    state.upstd[j] = 0.5 * (sqrt(b_coef^2 + 4.0 * a_coef * pkhsig) - b_coef) / a_coef

    # Mean porous velocity from pressure gradient
    sqrg1 = sqrt(2.0 / π)   # √(2/π)
    a_fric = sqrg1 * (b2_val + besta1 * state.upstd[j])
    c_cos2 = state.ctheta[j]^2
    state.upmean[j] = -dedx / (alsta + a_fric * (1.0 + c_cos2))

    # Stability clipping: limit |UPMEAN/UPSTD| ≤ 0.5
    upstd_j = state.upstd[j]
    if upstd_j > 0.0
        ratio = state.upmean[j] / upstd_j
        if ratio > 0.5
            state.upmean[j] = 0.5 * upstd_j
        elseif ratio < -0.5
            state.upmean[j] = -0.5 * upstd_j
        end
    end

    # Water flux through porous layer
    state.qp[j] = state.upmean[j] * hp_j

    # Energy dissipation
    a2 = state.upmean[j]^2
    b2 = state.upstd[j]^2
    state.dpsta[j] = hp_j * (alsta * (a2 + b2) + a_fric * (2.0 * b2 + a2 * (1.0 + 2.0 * c_cos2)))

    return nothing
end
