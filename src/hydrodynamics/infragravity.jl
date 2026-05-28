# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
infragravity.jl — Infragravity (IG) wave energy in the time-averaged model.

Infragravity waves (periods ≈ 25–300 s) dominate swash and runup on
dissipative, low-sloping beaches (Iribarren ξ < 0.3). Because CSHORE
averages over BC windows of O(hours), IG waves cannot appear as resolved
oscillations. They enter as a second bulk energy component alongside Hrms_ss.

Entry point: compute_ig_field!(state, config, wave, l)

Called once per BC window from _run_line_step! in driver.jl, immediately
after the short-wave convergence loop. Populates:
  state.hrms_ig[j]   — RMS IG wave height at each node (m)
  state.ustd_ig[j]   — RMS IG orbital velocity at each node (m/s)

Both are zeroed when config.ig is nothing (backward-compatible no-op).

Physics layers
--------------
Layer 1 (always active when IgConfig present):
  Algebraic IG height from offshore BC, shoaled by Green's law (h^{-1/4}),
  capped by the saturation limit gamma_ig · h.

  Offshore IG height:
    hrms_ig_0 = hrms_ig_bc[itime]     if provided by user
              = kappa_ig · hrms0      otherwise

  Shoaling (Green's law): hrms_ig[j] = hrms_ig_0 · (h0/h[j])^0.25
  Saturation cap:          hrms_ig[j] = min(hrms_ig[j], gamma_ig · h[j])

  IG orbital velocity (shallow-water limit, kh → 0):
    ustd_ig[j] = hrms_ig[j] · sqrt(g / h[j]) / 2

  The /2 comes from the RMS orbital velocity of a shallow-water wave:
    U_rms = H_rms · c / (2h) = H_rms · sqrt(g/h) / 2

Layer 2 (when ig.ig_energy_balance = true):
  1D IG energy balance marched landward from the offshore boundary:

    d(F_ig)/dx = kappa_source · |dF_ss/dx| − f_bf_ig · E_ig · sqrt(g/h)

  where:
    F_ig = E_ig · cg_ig    (IG energy flux, cg_ig = sqrt(g·h) in shallow water)
    E_ig = (1/8) · rho · g · hrms_ig²
    |dF_ss/dx| = |d(E_ss · cg_ss)/dx|  (SW energy flux gradient, from qbreak/dbsta)

  Spatial integration uses a simple forward-Euler march from j=1 (offshore)
  to j=jmax. The algebraic Layer-1 profile is used as the initial condition
  at j=1; Layer 2 then replaces it landward.

Layer 3 (when ig.ig_swash_active = true):
  Augments the swash-front water depth by c_ig_swash · hrms_ig at jwd.
  Currently stored as a diagnostic — the actual wet/dry coupling is in
  wetdry.jl (future integration point).
==============================================================================#

const RHO_WATER = 1025.0   # kg/m³ (seawater density for IG energy flux)

"""
    compute_ig_field!(state, config, wave, l)

Populate `state.hrms_ig` and `state.ustd_ig` for cross-shore line `l`.

When `config.ig` is `nothing`, both arrays are zeroed and the function
returns immediately (backward-compatible no-op).

# Arguments
- `state::CshoreState` — mutable model state (modified in-place)
- `config::CshoreConfig` — immutable configuration
- `wave::WaveParams` — current BC window wave parameters (tp, hrms0, swl)
- `l::Int` — cross-shore line index (1-based)
"""
function compute_ig_field!(state::CshoreState, config::CshoreConfig,
                           wave::WaveParams, l::Int, itime::Int)
    nn = length(state.hrms_ig)

    # ── No-op path ────────────────────────────────────────────────────────────
    if config.ig === nothing
        @inbounds for j in 1:nn
            state.hrms_ig[j] = 0.0
            state.ustd_ig[j] = 0.0
        end
        return nothing
    end

    ig  = config.ig
    bc  = config.boundary
    jmax_l = state.jmax[l]

    # ── Offshore IG height (BC override or algebraic estimate) ────────────────
    hrms_ig_0 = if !isempty(bc.hrms_ig_bc) && itime <= length(bc.hrms_ig_bc)
        max(bc.hrms_ig_bc[itime], 0.0)
    else
        ig.kappa_ig * wave.hrms0
    end

    # Reference offshore depth at j=1 (used for Green's-law shoaling).
    # Use the total water depth (including setup) at the seaward node.
    h0 = max(state.h[1], 0.01)

    # Still-water shoreline node: first node above SWL.
    # state.h[j] = 0 for j ≥ jswl (dry nodes).  Beyond jswl, we use the BDJ
    # swash depth state.hwd[j] — which is the time-averaged water depth during
    # swash occupation — as the physically meaningful depth for the IG
    # saturation cap and orbital velocity.  This is the only reliable depth
    # estimate in the swash zone; state.h[j] = 0 there and would force the
    # saturation limit to gamma_ig × 0.01 m = 0.003 m (physically wrong).
    jswl_l = state.jswl[l]   # first dry node (zb > SWL)

    # ── Layer 1: algebraic Green's law shoaling + saturation cap ─────────────
    # Wet zone (j ≤ jswl): use state.h[j] (wave setup + SWL depth).
    # Swash zone (j > jswl): shoaling continues from the last wet-zone value
    # but the saturation cap and orbital velocity use state.hwd[j] (BDJ
    # swash depth), which is only available after overtopping_rate! has run.
    @inbounds for j in 1:jmax_l
        # Choose depth: wave-field depth in wet zone, swash depth in dry zone.
        # hwd[j] > 0 only in the active swash (jwd:jdry); outside it is 0.
        h_wave  = state.h[j]
        h_swash = state.hwd[j]
        hj = if j < jswl_l || h_wave > h_swash
            max(h_wave, 0.01)
        else
            # Swash zone: prefer hwd when it gives a physical depth estimate.
            # Floor at hwdmin so we don't divide by zero or saturate trivially.
            max(h_swash, 0.001)
        end

        # Green's law: E ∝ h^{-1/2} → H ∝ h^{-1/4}
        hrms_ig_j = hrms_ig_0 * (h0 / hj)^0.25

        # Saturation cap: IG cannot exceed gamma_ig × local depth
        hrms_ig_j = min(hrms_ig_j, ig.gamma_ig * hj)

        state.hrms_ig[j] = hrms_ig_j
    end

    # Zero beyond the profile
    @inbounds for j in (jmax_l + 1):nn
        state.hrms_ig[j] = 0.0
    end

    # ── Layer 2: 1D IG energy balance (optional) ──────────────────────────────
    if ig.ig_energy_balance
        _ig_energy_balance!(state, config, wave, l, hrms_ig_0, h0)
    end

    # ── Compute ustd_ig from hrms_ig ─────────────────────────────────────────
    # Shallow-water orbital velocity: u_rms = H_rms · sqrt(g/h) / 2
    # (derived from U_rms = H_rms · c / (2h) with c = sqrt(g·h))
    # Use the same depth selection as above so that orbital velocity is
    # consistent with the IG height (no depth mismatch between hrms_ig and hj).
    @inbounds for j in 1:jmax_l
        h_wave  = state.h[j]
        h_swash = state.hwd[j]
        hj = if j < jswl_l || h_wave > h_swash
            max(h_wave, 0.01)
        else
            max(h_swash, 0.001)
        end
        state.ustd_ig[j] = 0.5 * state.hrms_ig[j] * sqrt(GRAV / hj)
    end
    @inbounds for j in (jmax_l + 1):nn
        state.ustd_ig[j] = 0.0
    end

    return nothing
end

"""
    _ig_energy_balance!(state, config, wave, l, hrms_ig_0, h0)

Layer 2: forward-Euler march of the 1D IG energy balance landward from
the offshore boundary.  Updates `state.hrms_ig` in-place.

Energy balance per unit width (marched from j=1 to j=jmax_l):

    d(F_ig)/dx = S_ig − D_ig

where:
    F_ig = E_ig · cg_ig            IG energy flux
    cg_ig = sqrt(g·h)              shallow-water group velocity
    E_ig  = (ρg/8) · hrms_ig²     IG wave energy density
    S_ig  = kappa_source · |dF_ss/dx|   source from SW breaking gradient
    D_ig  = f_bf_ig · E_ig · sqrt(g/h) / h   IG bottom-friction dissipation

The source is computed from the already-solved short-wave field via
`state.dbsta` (wave-breaking dissipation) and `state.dfsta` (bottom-friction
dissipation).  Their sum is proportional to |dF_ss/dx|.
"""
function _ig_energy_balance!(state::CshoreState, config::CshoreConfig,
                             wave::WaveParams, l::Int,
                             hrms_ig_0::Float64, h0::Float64)
    ig     = config.ig
    dx     = config.grid.dx
    jmax_l = state.jmax[l]

    kappa_src = ig.kappa_source
    f_bf_ig   = ig.f_bf_ig

    # IG energy at the seaward boundary (from Layer 1 initial estimate)
    h0_safe = max(h0, 0.01)
    E_ig = (RHO_WATER * GRAV / 8.0) * hrms_ig_0^2
    cg_ig = sqrt(GRAV * h0_safe)
    F_ig  = E_ig * cg_ig

    # March landward from j=1
    @inbounds for j in 1:jmax_l
        hj = max(state.h[j], 0.01)

        # IG energy density from current flux
        cg_ig_j = sqrt(GRAV * hj)
        E_ig_j  = max(F_ig / cg_ig_j, 0.0)

        # Update hrms_ig from energy balance result
        hrms_ig_j = sqrt(max(8.0 * E_ig_j / (RHO_WATER * GRAV), 0.0))

        # Apply saturation cap
        hrms_ig_j = min(hrms_ig_j, ig.gamma_ig * hj)
        state.hrms_ig[j] = hrms_ig_j

        # Recompute consistent E_ig after saturation
        E_ig_sat = (RHO_WATER * GRAV / 8.0) * hrms_ig_j^2

        if j < jmax_l
            hj1 = max(state.h[j + 1], 0.01)
            cg_ig_j1 = sqrt(GRAV * hj1)
            dx_eff = dx

            # Source: short-wave breaking + friction dissipation at this node
            # dbsta = wave-breaking dissipation rate (W/m²); dfsta = friction diss.
            S_ig = kappa_src * (state.dbsta[j] + state.dfsta[j])

            # Dissipation: IG bottom friction
            D_ig = f_bf_ig * E_ig_sat * sqrt(GRAV / hj) / max(hj, 0.01)

            # Forward Euler: dF/dx = S - D  →  F_{j+1} = F_j + dx * (S - D)
            F_ig_next = F_ig + dx_eff * (S_ig - D_ig)
            F_ig = max(F_ig_next, 0.0)   # energy flux cannot go negative
        end
    end

    return nothing
end
