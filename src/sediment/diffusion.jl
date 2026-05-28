# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
diffusion.jl — Hillslope diffusion for non-wave gravity-driven mass wasting.
==============================================================================#

@inline function compute_diffusivity(slope::Float64, cfg)
    S = atan(abs(slope))
    tan_ratio = tan(S) / tan(cfg.critical_slope)
    D = cfg.D_base * (1.0 + tan_ratio^cfg.slope_exponent * cfg.critical_factor)
    return clamp(D, 0.0, cfg.max_diffusivity)
end

@inline function compute_diffusivity_local(slope::Float64, D_base::Float64,
    crit_slope::Float64,
    slope_exp::Float64,
    crit_factor::Float64,
    max_diff::Float64)
    S = atan(abs(slope))
    tan_ratio = tan(S) / tan(crit_slope)
    D = D_base * (1.0 + tan_ratio^slope_exp * crit_factor)
    return clamp(D, 0.0, max_diff)
end

function apply_hillslope_diffusion!(state::CshoreState, config::CshoreConfig,
    l::Int, dt_seconds::Float64)
    dcfg = config.diffusion
    dcfg === nothing && return nothing
    dcfg.D_base ≤ 0.0 && return nothing

    dx = config.grid.dx
    jmax = state.jmax[l]
    jmax < 3 && return nothing

    dt_days = dt_seconds / 86400.0
    dt_days ≤ 0.0 && return nothing

    has_wet_config = isfinite(dcfg.wet_critical_slope)
    wet_d_base = isfinite(dcfg.wet_d_base) ? dcfg.wet_d_base : dcfg.D_base
    wet_crit = isfinite(dcfg.wet_critical_slope) ? dcfg.wet_critical_slope : dcfg.critical_slope
    wet_max_diff = isfinite(dcfg.wet_max_diffusivity) ? dcfg.wet_max_diffusivity : dcfg.max_diffusivity

    is_wet_zone = falses(jmax)   # swash + buffer cells get wet diffusion
    is_full_submerged = falses(jmax) # fully submerged: zero diffusion
    if has_wet_config
        @inbounds for j in 1:jmax
            pw = state.pwet[j]
            if pw > dcfg.swash_pwet_max
                is_full_submerged[j] = true
            elseif pw > dcfg.swash_pwet_min
                is_wet_zone[j] = true
            end
        end
        # Add `swash_buffer_cells` landward of the most-landward swash cell.
        if dcfg.swash_buffer_cells > 0
            last_wet = findlast(is_wet_zone)
            if last_wet !== nothing
                last_buf = min(jmax, last_wet + dcfg.swash_buffer_cells)
                @inbounds for j in (last_wet+1):last_buf
                    is_wet_zone[j] = true
                end
            end
        end
    end

    D = Vector{Float64}(undef, jmax)
    @inbounds for j in 1:jmax
        if j == 1
            slope = (state.zb[2, l] - state.zb[1, l]) / dx
        elseif j == jmax
            slope = (state.zb[jmax, l] - state.zb[jmax-1, l]) / dx
        else
            slope = (state.zb[j+1, l] - state.zb[j-1, l]) / (2.0 * dx)
        end
        if is_wet_zone[j]
            D_raw = compute_diffusivity_local(slope, wet_d_base, wet_crit,
                dcfg.slope_exponent,
                dcfg.critical_factor, wet_max_diff)

            if dcfg.wet_pwet_scaling
                pw_scale = clamp(state.pwet[j], dcfg.wet_pwet_scale_floor, 1.0)
                D[j] = D_raw * pw_scale
            else
                D[j] = D_raw
            end
        else
            D[j] = compute_diffusivity(slope, dcfg)
        end
    end

    has_thermal = state.thermal !== nothing
    @inbounds for j in 1:jmax
        blocked = false

        if has_wet_config
            if is_full_submerged[j]
                blocked = true
            end
        elseif state.swldep[j, l] > 0.0
            blocked = true
        end

        # Thermal/hardbottom control (always, when thermal_control enabled)
        if dcfg.thermal_control && !blocked
            if has_thermal
                alt_j = state.thermal.ALT[j]
                if alt_j < dcfg.thaw_threshold
                    blocked = true
                end
            end
            zh = state.zb_hard[j, l]
            if zh > -1e20
                hp_j = state.zb[j, l] - zh
                if hp_j ≤ 0.0
                    blocked = true
                end
            end
        end

        if blocked
            D[j] = 0.0
        end
    end

    D_max = maximum(D)
    D_max ≤ 0.0 && return nothing  # nothing to diffuse

    dt_stable = dx^2 / (2.0 * D_max)  # CFL limit (days)
    n_sub = max(1, isfinite(dt_stable) && dt_stable > 0 ? ceil(Int, dt_days / dt_stable) : 1)
    dt_sub = dt_days / n_sub

    zb_work = Vector{Float64}(undef, jmax)
    @inbounds for j in 1:jmax
        zb_work[j] = state.zb[j, l]
    end
    flux = Vector{Float64}(undef, jmax - 1)  # face fluxes (j+½)

    for _isub in 1:n_sub
        @inbounds for f in 1:(jmax-1)
            D_face = min(D[f], D[f+1])
            dz = zb_work[f+1] - zb_work[f]
            flux[f] = -D_face * dz / dx
            if flux[f] > 0
                zh_src = state.zb_hard[f, l]
                if zh_src > -1e20
                    avail = max(0.0, zb_work[f] - zh_src)
                    max_flux = avail * dx / dt_sub
                    flux[f] = min(flux[f], max_flux)
                end
            elseif flux[f] < 0
                zh_src = state.zb_hard[f+1, l]
                if zh_src > -1e20
                    avail = max(0.0, zb_work[f+1] - zh_src)
                    max_flux = avail * dx / dt_sub
                    flux[f] = max(flux[f], -max_flux)
                end
            end
        end

        @inbounds for j in 2:(jmax-1)
            dzb = dt_sub * (flux[j-1] - flux[j]) / dx
            zb_work[j] += dzb
        end

        @inbounds for j in 1:jmax
            zh = state.zb_hard[j, l]
            if zh > -1e20 && zb_work[j] < zh
                zb_work[j] = zh
            end
        end

        @inbounds for j in 1:jmax
            if D[j] > 0  # only for active (non-blocked) cells
                if j == 1
                    slope = (zb_work[2] - zb_work[1]) / dx
                elseif j == jmax
                    slope = (zb_work[jmax] - zb_work[jmax-1]) / dx
                else
                    slope = (zb_work[j+1] - zb_work[j-1]) / (2.0 * dx)
                end
                if is_wet_zone[j]
                    D_raw = compute_diffusivity_local(slope, wet_d_base, wet_crit,
                        dcfg.slope_exponent,
                        dcfg.critical_factor, wet_max_diff)
                    if dcfg.wet_pwet_scaling
                        pw_scale = clamp(state.pwet[j], dcfg.wet_pwet_scale_floor, 1.0)
                        D[j] = D_raw * pw_scale
                    else
                        D[j] = D_raw
                    end
                else
                    D[j] = compute_diffusivity(slope, dcfg)
                end
            end
        end
    end

    nf = size(state.bed_mass, 3)
    nlayers = size(state.bed_mass, 2)
    ρs = config.sediment.sg * 1000.0
    one_minus_n = 1.0 - config.multifraction.porosity
    mass_per_m = ρs * one_minus_n

    @inbounds for j in 1:jmax
        dz = zb_work[j] - state.zb[j, l]
        if abs(dz) > 1e-15
            state.zb[j, l] = zb_work[j]

            total_surf = 0.0
            for k in 1:nf
                total_surf += max(0.0, state.bed_mass[j, 1, k])
            end
            dm_total = dz * mass_per_m  # positive = accretion, negative = erosion

            if total_surf > 0
                for k in 1:nf
                    frac_k = max(0.0, state.bed_mass[j, 1, k]) / total_surf
                    state.bed_mass[j, 1, k] += dm_total * frac_k
                    state.bed_mass[j, 1, k] = max(0.0, state.bed_mass[j, 1, k])
                end
            elseif dm_total > 0
                # Accretion onto empty surface: use initial fractions
                for k in 1:nf
                    state.bed_mass[j, 1, k] += dm_total * config.multifraction.initial_fractions[k]
                end
            end
        end
    end

    return nothing
end
