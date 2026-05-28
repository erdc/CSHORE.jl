# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
thermal.jl — 1D heat diffusion permafrost / active-layer model.
==============================================================================#

"""
    initialize_thermal_state(config, cstate) -> ThermalState

Allocate the `ThermalState` for a simulation. Distributes `n_rep`
representative columns uniformly across the active cross-shore domain
(j = 1..jmax), initializes `T ≡ T_init`, `ice_frac ≡ 1`, and sets
`ALT = 0` (winter start).
"""
function initialize_thermal_state(thconfig::ThermalConfig, state::CshoreState)
    dz = thconfig.dz
    nn = length(state.xb)

    nz_node = Vector{Int}(undef, nn)
    if isfinite(thconfig.z_bottom)
        @inbounds for j in 1:nn
            col_depth = state.zb[j, 1] - thconfig.z_bottom
            nz_node[j] = max(1, isfinite(col_depth) ? ceil(Int, col_depth / dz) : thconfig.nz)
        end
    else
        fill!(nz_node, thconfig.nz)
    end
    nz_max = maximum(nz_node)

    T = fill(thconfig.T_init, nz_max, nn)
    init_ice = thconfig.T_init < 0.0 ? 1.0 : 0.0
    ice_frac = fill(init_ice, nz_max, nn)

    prof = thconfig.T_init_profile
    if !isempty(prof)
        prof_depths = Float64[p[1] for p in prof]
        prof_temps = Float64[p[2] for p in prof]
        swl_init = 0.0

        @inbounds for j in 1:nn
            nz_j = nz_node[j]
            if state.zb[j, 1] > swl_init
                for iz in 1:nz_j
                    z = (iz - 0.5) * dz
                    if z <= prof_depths[1]
                        T[iz, j] = prof_temps[1]
                    elseif z >= prof_depths[end]
                        T[iz, j] = prof_temps[end]
                    else
                        for k in 1:length(prof_depths)-1
                            if prof_depths[k] <= z <= prof_depths[k+1]
                                frac = (z - prof_depths[k]) / (prof_depths[k+1] - prof_depths[k])
                                T[iz, j] = prof_temps[k] + frac * (prof_temps[k+1] - prof_temps[k])
                                break
                            end
                        end
                    end
                    ice_frac[iz, j] = T[iz, j] < 0.0 ? 1.0 : 0.0
                end
            end
        end
    end

    ALT = zeros(Float64, nn)
    z_top_last = copy(state.zb[:, 1])

    # Thermal intervention vectors: copy from config if provided, else zeros.
    R_insulation = if !isempty(thconfig.R_insulation)
        length(thconfig.R_insulation) >= nn ?
        thconfig.R_insulation[1:nn] : vcat(thconfig.R_insulation, zeros(nn - length(thconfig.R_insulation)))
    else
        zeros(Float64, nn)
    end
    Q_thermosyphon = if !isempty(thconfig.Q_thermosyphon)
        length(thconfig.Q_thermosyphon) >= nn ?
        thconfig.Q_thermosyphon[1:nn] : vcat(thconfig.Q_thermosyphon, zeros(nn - length(thconfig.Q_thermosyphon)))
    else
        zeros(Float64, nn)
    end

    snow_depth = zeros(Float64, nn)
    R_sod_base = copy(R_insulation)   # permanent sod component

    # Preserve the initial structural hardbottom so thermal updates
    # never erase non-erodible structures (breakwaters, revetments).
    zb_hard_init = copy(state.zb_hard[:, 1])

    zb_hard_scour_floor = state.zb[:, 1] .- (2.0 * thconfig.alt_max)

    return ThermalState(T, ice_frac, ALT, nz_node, z_top_last, R_insulation, Q_thermosyphon,
        snow_depth, R_sod_base, zb_hard_init, zb_hard_scour_floor)
end

"""
    _step_heat_column!(T, ice_frac, T_surface, thconfig, dt)

Advance one thermal column by wall-clock time `dt` seconds with surface
boundary condition `T_surface` (°C) and fixed lower BC `thconfig.T_lower`.
"""
function _step_heat_column!(T::AbstractVector{Float64},
    ice_frac::AbstractVector{Float64},
    T_surface::Float64,
    thconfig::ThermalConfig,
    dt::Float64,
    moisture_eff::Float64=thconfig.moisture,
    R_insul::Float64=0.0,
    Q_thermo::Float64=0.0,
    nz::Int=thconfig.nz)
    dz = thconfig.dz
    # Diffusivity per cell, evaluated from phase.
    @inline function α_of(cell::Int)
        if T[cell] < 0 && ice_frac[cell] > 0.5
            return thconfig.k_frozen / thconfig.C_frozen
        else
            return thconfig.k_thawed / thconfig.C_thawed
        end
    end
    @inline function k_of(cell::Int)
        if T[cell] < 0 && ice_frac[cell] > 0.5
            return thconfig.k_frozen
        else
            return thconfig.k_thawed
        end
    end
    @inline function C_of(cell::Int)
        if T[cell] < 0 && ice_frac[cell] > 0.5
            return thconfig.C_frozen
        else
            return thconfig.C_thawed
        end
    end

    # Sub-step count for CFL: dt_sub ≤ cfl * dz² / α_max
    α_max = max(thconfig.k_frozen / thconfig.C_frozen,
        thconfig.k_thawed / thconfig.C_thawed)
    dt_cfl = thconfig.cfl_safety * dz * dz / α_max
    n_sub = max(1, ceil(Int, dt / dt_cfl))
    dt_sub = dt / n_sub

    L_eff = thconfig.L * moisture_eff

    flux = Vector{Float64}(undef, nz + 1)   # face fluxes
    @inbounds for _ in 1:n_sub
        # --- Compute face fluxes q_face = -k * dT/dz ---
        # Top face: surface BC with optional insulation resistance.
        # Without insulation: flux = -k * (T[1] - T_surface) / (0.5*dz)
        # With insulation R (m²·K/W) in series with the half-cell:
        #   flux = -(T[1] - T_surface) / (0.5*dz/k + R)
        # When R=0 this reduces to the standard formula.
        k_top = k_of(1)
        flux[1] = -(T[1] - T_surface) / (0.5 * dz / k_top + R_insul)
        # Thermosyphon: passive one-way heat extraction (winter only).
        # Active when T_surface < T[1] (air colder than ground surface).
        # flux[1] sign convention: negative = upward (heat leaving cell 1).
        # Thermosyphon enhances upward flux (more cooling), so subtract.
        # Capped at the energy to bring T[1] to T_surface in one sub-step.
        if Q_thermo > 0.0 && T_surface < T[1]
            Q_max = C_of(1) * (T[1] - T_surface) * dz / dt_sub
            flux[1] -= min(Q_thermo, Q_max)
        end
        for i in 2:nz
            k_face = 0.5 * (k_of(i - 1) + k_of(i))
            flux[i] = -k_face * (T[i] - T[i-1]) / dz
        end
        # Bottom face: Dirichlet BC at T_lower (half-cell).
        k_bot = k_of(nz)
        flux[nz+1] = -k_bot * (thconfig.T_lower - T[nz]) / (0.5 * dz)

        # --- Update each cell ---
        # State machine:
        #   (A) 0 < ice_frac < 1  →  partially thawed, T pinned at 0, all
        #       heat goes to latent (melt/freeze).  If latent is exhausted
        #       (ice crosses 0 or 1), any remaining dQ is applied as
        #       single-phase heating/cooling in the new phase.
        #   (B) ice_frac == 1 && T ≤ 0 → fully frozen single-phase heating;
        #       if T would cross 0, the excess dQ converts to latent.
        #   (C) ice_frac == 0 && T ≥ 0 → fully thawed single-phase cooling;
        #       if T would cross 0, the deficit dQ converts to latent.
        #   Edge cases T==0 with ice_frac==1 or 0 are handled by direction
        #   of dQ: heat in at T=0,ice=1 → start melting (case A branch);
        #   heat out at T=0,ice=0 → start freezing (case A branch).
        for i in 1:nz
            C_i = C_of(i)
            dQ = (flux[i] - flux[i+1]) * dt_sub / dz    # J/m³ energy added
            ifr = ice_frac[i]
            Ti = T[i]

            # Decide which regime the cell is in BEFORE applying dQ.
            is_mid = 0.0 < ifr < 1.0
            is_frozen = ifr >= 1.0 && Ti <= 0.0
            is_thawed = ifr <= 0.0 && Ti >= 0.0
            # "T=0, ice=1, dQ>0" → start melting: treat as mid-phase.
            if !is_mid && Ti == 0.0 && ifr == 1.0 && dQ > 0.0
                is_mid = true
            end
            # "T=0, ice=0, dQ<0" → start freezing: treat as mid-phase.
            if !is_mid && Ti == 0.0 && ifr == 0.0 && dQ < 0.0
                is_mid = true
            end

            if is_mid
                # All heat goes into latent until phase is pushed to a
                # boundary, then any residual heats/cools the new phase.
                d_ice = -dQ / L_eff
                new_ice = ifr + d_ice
                T[i] = 0.0
                if new_ice < 0.0
                    # Finished melting; residual heat warms thawed cell.
                    residual = -new_ice * L_eff   # (J/m³, positive)
                    ice_frac[i] = 0.0
                    T[i] = residual / thconfig.C_thawed
                elseif new_ice > 1.0
                    # Finished freezing; residual removes heat from frozen.
                    residual = (new_ice - 1.0) * L_eff   # J/m³ still to remove
                    ice_frac[i] = 1.0
                    T[i] = -residual / thconfig.C_frozen
                else
                    ice_frac[i] = new_ice
                end
            elseif is_frozen
                # Single-phase frozen heating/cooling. If dQ > 0 and we
                # would cross 0, pin at 0 and send the remaining energy
                # to latent (start melting).
                dT_proposed = dQ / C_i
                if Ti + dT_proposed > 0.0 && dQ > 0.0
                    used = (-Ti) * C_i            # energy to reach 0
                    leftover = dQ - used              # remainder → latent
                    T[i] = 0.0
                    ice_frac[i] = clamp(1.0 - leftover / L_eff, 0.0, 1.0)
                else
                    T[i] = Ti + dT_proposed
                end
            elseif is_thawed
                # Single-phase thawed heating/cooling. If dQ < 0 and we
                # would cross 0, pin at 0 and send the remaining energy
                # to latent (start freezing).
                dT_proposed = dQ / C_i
                if Ti + dT_proposed < 0.0 && dQ < 0.0
                    used = (-Ti) * C_i            # energy to reach 0 (negative)
                    leftover = dQ - used              # remainder (negative) → latent
                    T[i] = 0.0
                    ice_frac[i] = clamp(-leftover / L_eff, 0.0, 1.0)
                else
                    T[i] = Ti + dT_proposed
                end
            else
                # Shouldn't happen given the above decision tree, but
                # leave a safe single-phase update as a fallback.
                T[i] = Ti + dQ / C_i
            end
        end
    end
    return nothing
end

"""
    _shift_column_up!(T, ice_frac, n_shift, thconfig)

Shift a thermal column upward by `n_shift` cells to account for bed
erosion of (approximately) `n_shift * dz` meters that has occurred
between thermal calls. After the shift:
  - cells 1..(nz - n_shift)  ← previous cells (n_shift+1..nz)
  - cells (nz - n_shift + 1..nz)  ← fully frozen at T_lower
"""
function _shift_column_up!(T::AbstractVector{Float64},
    ice_frac::AbstractVector{Float64},
    n_shift::Int, thconfig::ThermalConfig,
    nz::Int=thconfig.nz)
    n_shift = min(max(n_shift, 0), nz)
    n_shift == 0 && return nothing
    # Move cells upward: T[i] = T[i + n_shift] for i in 1..nz-n_shift
    @inbounds for i in 1:(nz-n_shift)
        T[i] = T[i+n_shift]
        ice_frac[i] = ice_frac[i+n_shift]
    end
    # Fill the newly-empty bottom with deep-permafrost state.
    @inbounds for i in (nz-n_shift+1):nz
        T[i] = thconfig.T_lower
        ice_frac[i] = 1.0
    end
    return nothing
end

"""
    _column_alt(T, ice_frac, thconfig) -> Float64

Compute the active-layer thickness (distance from the surface to the
shallowest 0 °C isotherm — i.e., the deepest fully-thawed cell) in the
given column. Returns 0 if no cell is thawed.
"""
function _column_alt(T::AbstractVector{Float64},
    ice_frac::AbstractVector{Float64},
    thconfig::ThermalConfig,
    nz::Int=thconfig.nz)
    alt = 0.0
    @inbounds for i in 1:nz
        if ice_frac[i] < 0.5     # mostly or fully thawed
            alt = i * thconfig.dz
        else
            break                 # below the thaw front — done
        end
    end
    return clamp(alt, thconfig.alt_min, thconfig.alt_max)
end

"""
    update_snow!(thstate, snow_config, state, config, T_air, dt, l=1;
                 prescribed_depth=NaN)

Update per-node snow depth and recompute R_insulation (sod + snow) before
the thermal solve. Two modes:

1. **Prescribed**: if `prescribed_depth` is not NaN, set all exposed nodes
   to that depth (user-supplied time series, interpolated by the caller).
2. **Degree-day model**: accumulate when T_air < 0, melt when T_air > 0,
   capped at `snow_config.max_depth`.

In both modes:
- Submerged nodes (zb < swl): snow_depth = 0
- Swash zone: snow reduced by wave action — `snow_depth *= (1 - pwet)`
- R_insulation = R_sod_base + snow_depth / k_snow
"""
function update_snow!(thstate::ThermalState, snow_config::SnowConfig,
    state::CshoreState, config::CshoreConfig,
    T_air::Float64, dt::Float64, l::Int=1;
    prescribed_depth::Float64=NaN)
    jmax = state.jmax[l]
    swl = config.boundary.swlbc[max(1, state.itime)]
    dt_hours = dt / 3600.0
    k_snow = snow_config.k_snow
    use_pwet = (config.options.iover == 1) && (state.jdry > 0) && (state.jwd > 0)
    smod = config.snow_modifier

    @inbounds for j in 1:jmax
        zb_now = state.zb[j, l]

        if zb_now < swl
            # Submerged: no snow
            thstate.snow_depth[j] = 0.0
        elseif !isnan(prescribed_depth)
            # Prescribed snow depth from user time series
            thstate.snow_depth[j] = clamp(prescribed_depth, 0.0, snow_config.max_depth)
        else
            # Degree-day accumulation/melt model
            sd = thstate.snow_depth[j]
            if T_air < 0.0
                sd += snow_config.accum_rate * dt_hours
            else
                sd -= snow_config.melt_rate * T_air * dt_hours
            end
            thstate.snow_depth[j] = clamp(sd, 0.0, snow_config.max_depth)
        end

        if smod !== nothing && zb_now >= swl
            sd = thstate.snow_depth[j]
            sd = max(sd, smod.depth_min[j])
            sd = min(sd, smod.depth_max[j])
            thstate.snow_depth[j] = clamp(sd, 0.0, snow_config.max_depth)
        end

        # Swash zone: wave action removes snow proportional to wetting
        if use_pwet && j >= state.jwd && j <= state.jdry
            pw = clamp(state.pwet[j], 0.0, 1.0)
            thstate.snow_depth[j] *= (1.0 - pw)
        end

        # Update total surface resistance: permanent sod + time-varying snow
        thstate.R_insulation[j] = thstate.R_sod_base[j] +
                                  thstate.snow_depth[j] / k_snow
    end
    return nothing
end

"""
    update_thermal!(state, config, thconfig, thstate, T_air, T_water, dt, l=1)

One thermal BC step:
1. For each representative column, pick the surface BC: `T_water` if the
   node is submerged (SWL above bed), else `T_air`.
2. Advance the column by `dt` seconds with `_step_heat_column!`.
3. Read the column ALT.
4. Interpolate ALT linearly along x onto every shore node.
5. Write `state.zb_hard[j, l] = state.zb[j, l] - ALT[j]`, which the
   existing hardbottom infrastructure (item 1) then enforces during the
   wave / sediment sub-steps.

`dt` is in seconds. Intended to be called once per BC window from
`step_bc_window!`.
"""
function update_thermal!(state::CshoreState, config::CshoreConfig,
    thconfig::ThermalConfig, thstate::ThermalState,
    T_air::Float64, T_water::Float64, dt::Float64,
    l::Int=1)
    jmax = state.jmax[l]
    swl = config.boundary.swlbc[max(1, state.itime)]
    dz = thconfig.dz

    use_pwet = (config.options.iover == 1) && (state.jdry > 0) && (state.jwd > 0)

    @inbounds for j in 1:jmax
        zb_now = state.zb[j, l]
        dz_eroded = thstate.z_top_last[j] - zb_now
        if dz_eroded > 0.0
            n_shift = isfinite(dz_eroded) ? floor(Int, dz_eroded / dz) : 0
            if n_shift > 0
                nz_j = thstate.nz_node[j]
                _shift_column_up!(view(thstate.T, :, j),
                    view(thstate.ice_frac, :, j),
                    n_shift, thconfig, nz_j)
                thstate.z_top_last[j] -= n_shift * dz
            end
            # else: keep accumulating sub-dz erosion in z_top_last
        else
            thstate.z_top_last[j] = zb_now
        end

        if use_pwet
            if j < state.jwd
                T_surface = T_water                         # seaward: always submerged
                moisture_eff = 1.0                          # saturated
            elseif j <= state.jdry
                pw = clamp(state.pwet[j], 0.0, 1.0)        # swash zone: intermittent
                T_surface = pw * T_water + (1.0 - pw) * T_air
                moisture_eff = pw * 1.0 + (1.0 - pw) * thconfig.moisture
            else
                T_surface = zb_now < swl ? T_water : T_air  # landward: sharp fallback
                moisture_eff = zb_now < swl ? 1.0 : thconfig.moisture
            end
        else
            T_surface = zb_now < swl ? T_water : T_air      # no swash data: sharp switch
            moisture_eff = zb_now < swl ? 1.0 : thconfig.moisture
        end

        Tcol = view(thstate.T, :, j)
        icol = view(thstate.ice_frac, :, j)
        R_j = thstate.R_insulation[j]
        Q_j = thstate.Q_thermosyphon[j]
        nz_j = thstate.nz_node[j]
        _step_heat_column!(Tcol, icol, T_surface, thconfig, dt, moisture_eff, R_j, Q_j, nz_j)
        alt_j = _column_alt(Tcol, icol, thconfig, nz_j)
        thstate.ALT[j] = alt_j
        thermal_hard = zb_now - alt_j
        structural_hard = thstate.zb_hard_init[j]
        scour_floor = thstate.zb_hard_scour_floor[j]
        state.zb_hard[j, l] = max(thermal_hard, structural_hard, scour_floor)
    end
    return nothing
end
