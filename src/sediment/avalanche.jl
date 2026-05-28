# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

@inline function _apply_avalanche_dz!(state::CshoreState, config::CshoreConfig,
    j::Int, l::Int, dz::Float64,
    mass_per_m::Float64, nf::Int)
    abs(dz) > 1e-15 || return
    state.zb[j, l] += dz
    total_surf = 0.0
    @inbounds for k in 1:nf
        total_surf += max(0.0, state.bed_mass[j, 1, k])
    end
    dm_total = dz * mass_per_m
    if total_surf > 0.0
        @inbounds for k in 1:nf
            frac_k = max(0.0, state.bed_mass[j, 1, k]) / total_surf
            state.bed_mass[j, 1, k] += dm_total * frac_k
            if state.bed_mass[j, 1, k] < 0.0
                state.bed_mass[j, 1, k] = 0.0
            end
        end
    elseif dm_total > 0.0
        @inbounds for k in 1:nf
            state.bed_mass[j, 1, k] +=
                dm_total * config.multifraction.initial_fractions[k]
        end
    end
    return
end

function apply_underwater_avalanche!(state::CshoreState, config::CshoreConfig,
    l::Int)
    tanphi = config.sediment.tanphi
    tanphi > 0.0 || return nothing
    jmax = state.jmax[l]
    jmax >= 2 || return nothing

    dx = config.grid.dx
    max_face_dz = tanphi * dx
    nf = size(state.bed_mass, 3)
    ρs = config.sediment.sg * 1000.0
    one_minus_n = 1.0 - config.multifraction.porosity
    mass_per_m = ρs * one_minus_n
    max_sweeps = 200

    for _sweep in 1:max_sweeps
        any_fix = false
        @inbounds for j in 1:(jmax-1)
            jp1 = j + 1
            # Underwater gate — both cells must be below SWL.
            (state.swldep[j, l] > 0.0 && state.swldep[jp1, l] > 0.0) || continue
            dz_face = state.zb[jp1, l] - state.zb[j, l]
            abs(dz_face) > max_face_dz || continue
            # Target: enforce slope = ±tanphi while preserving (zb[j] + zb[jp1])
            tot = state.zb[j, l] + state.zb[jp1, l]
            target_gap = sign(dz_face) * max_face_dz
            zb_new_j = 0.5 * (tot - target_gap)
            zb_new_jp1 = 0.5 * (tot + target_gap)
            # Hardbottom guards — skip face if relaxation would breach either bed.
            zh_j = state.zb_hard[j, l]
            zh_jp1 = state.zb_hard[jp1, l]
            if (zh_j > -1e20 && zb_new_j < zh_j) ||
               (zh_jp1 > -1e20 && zb_new_jp1 < zh_jp1)
                continue
            end
            _apply_avalanche_dz!(state, config, j, l, zb_new_j - state.zb[j, l], mass_per_m, nf)
            _apply_avalanche_dz!(state, config, jp1, l, zb_new_jp1 - state.zb[jp1, l], mass_per_m, nf)
            any_fix = true
        end
        any_fix || break
    end
    return nothing
end
