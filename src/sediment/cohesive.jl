# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
sediment/cohesive.jl — Minimal-viable cohesive (mud) sediment transport.

Adds a single suspended cohesive fraction tracked separately from the
existing multifraction capability.
==============================================================================#

"""
    initial_cohesive_bed_mass(cohesive::CohesiveSedimentConfig, nn::Int) -> Vector{Float64}

Build the per-node initial cohesive bed-mass vector (uniform fill).
Callers that need spatial variation should overwrite the returned vector
in place before stepping.
"""
initial_cohesive_bed_mass(cohesive::CohesiveSedimentConfig, nn::Int) =
    fill(cohesive.initial_bed_mass, nn)

"""
    cohesive_step!(state, config, l, dt) -> nothing

Advance the cohesive bed mass and suspended concentration on line `l` by
`dt` seconds using a one-step explicit Partheniades-Krone update.
"""
function cohesive_step!(state::CshoreState, config::CshoreConfig, l::Int, dt::Float64)
    coh = config.cohesive
    coh === nothing && return nothing
    dt > 0 || return nothing

    rho_w = coh.rho_water
    ws = coh.settling_velocity
    tau_ce = coh.tau_ce
    tau_cd = coh.tau_cd
    M = coh.M

    # Convert per-node depth-averaged concentration into a depth-times-C
    # so we can rewrite the explicit Euler step in mass-per-area units.
    jmax_l = state.jmax[l]
    @inbounds for j in 1:jmax_l
        # Skip dry nodes — no water column to support suspension.
        h_j = state.h[j]
        h_j > 1e-6 || continue

        # Bed shear in Pa from CSHORE's normalized tbxsta.
        tau_b = abs(state.tbxsta[j]) * rho_w * GRAV

        # Current pools (kg/m² and kg/m³ respectively).
        bed_j = state.cohesive_bed_mass[j]
        conc_j = state.cohesive_concentration[j]

        # Erosion flux E (kg/m²/s), only when τ_b exceeds the erosion threshold.
        E = tau_b > tau_ce ? M * (tau_b / tau_ce - 1.0) : 0.0
        # Cap by available bed mass over this step (mass-conservative).
        E_max = bed_j / dt
        E = min(E, E_max)

        # Deposition flux D (kg/m²/s), only when τ_b is below the deposition
        # threshold AND there is something suspended to fall out.
        D = (tau_b < tau_cd && conc_j > 0.0) ?
            ws * conc_j * (1.0 - tau_b / tau_cd) : 0.0
        # Cap by available suspended mass over this step (= conc·h).
        D_max = (conc_j * h_j) / dt
        D = min(D, D_max)

        # Update pools. ΔC over the column is (E - D)·dt / h.
        dC = (E - D) * dt / h_j
        new_conc = max(conc_j + dC, 0.0)
        new_bed = max(bed_j + (D - E) * dt, 0.0)

        state.cohesive_concentration[j] = new_conc
        state.cohesive_bed_mass[j] = new_bed
    end
    return nothing
end
