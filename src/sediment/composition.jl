# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
composition.jl — Per-fraction bed composition derived from AeoLiS
==============================================================================#
# ----------------------------------------------------------------------------
const MASS_DEBUG_ON = Ref(false)

mutable struct MassDebugCounters
    n_clamp_fired::Int           # Step 8: how many (j, layer, k) entered with neg
    m_clamp_unrecovered::Float64 # Step 8: total kg created by unrecovered deficits
    m_step5_active_removed::Float64  # Step 5: total kg removed from active layer
    m_step7_reservoir::Float64   # Step 7: total kg added by reservoir feed
    n_thin_cover_nodes::Int      # nodes where cover < layer_thickness this call
end
MassDebugCounters() = MassDebugCounters(0, 0.0, 0.0, 0.0, 0)
const MASS_DEBUG_COUNTERS = MassDebugCounters()

function enable_mass_debug!(on::Bool=true)
    MASS_DEBUG_ON[] = on
    return nothing
end

function reset_mass_debug!()
    c = MASS_DEBUG_COUNTERS
    c.n_clamp_fired = 0
    c.m_clamp_unrecovered = 0.0
    c.m_step5_active_removed = 0.0
    c.m_step7_reservoir = 0.0
    c.n_thin_cover_nodes = 0
    return nothing
end

function _normalize_frac(m::AbstractVector{Float64})
    s = sum(m)
    nf = length(m)
    out = similar(m)
    if s > 0
        @inbounds for k in 1:nf
            out[k] = m[k] / s
        end
    else
        fill!(out, 0.0)
    end
    return out
end

function prevent_negative_mass!(
    m::Array{Float64,3},
    dm::Vector{Float64},
    pickup::Matrix{Float64},
    jmax::Int,
    nlayers::Int,
    nf::Int,
)
    ero = zeros(nf)
    dep = zeros(nf)
    ddep = zeros(nf)

    @inbounds for j in 1:jmax
        if dm[j] > 0
            # Split pickup into erosional / depositional components (per fraction)
            erog = 0.0
            for k in 1:nf
                ero[k] = max(0.0, pickup[j, k])
                dep[k] = -min(0.0, pickup[j, k])    # depositional magnitude (positive)
                erog += ero[k]
            end
            if erog > 0
                # Normalized deposition distribution for swap
                dep_sum = 0.0
                for k in 1:nf
                    dep_sum += dep[k]
                end
                if dep_sum > 0
                    for k in 1:nf
                        ddep[k] = erog * (dep[k] / dep_sum)
                    end
                    # Rewrite pickup: eroding fractions cancelled, depositing
                    # fractions reduced by `ddep` (net deposition unchanged).
                    new_dm = 0.0
                    for k in 1:nf
                        pickup[j, k] = -dep[k] + ddep[k]
                        new_dm -= pickup[j, k]
                    end
                    dm[j] = new_dm
                    # Adjust layer-1 mass to reflect the swap: remove erosional
                    # fractions outright and add the swap-equivalent depositional
                    # mass (so the fraction composition is already updated before
                    # the pickup loop below).
                    for k in 1:nf
                        m[j, 1, k] -= (ero[k] - ddep[k])
                    end
                end
            end
        end

        # Deposition larger than one full layer mass. Fill top layers with
        # fresh deposit and shift existing sediment down whole-layers,
        # until remaining deposit is sub-layer-sized.
        if dm[j] > 0
            mx = 0.0
            for k in 1:nf
                mx += m[j, 1, k]
            end
            if mx > 0
                # Clamp before Int conversion: dm[j]/mx can be a large finite
                # number when mx is near-zero, causing floor(Int,...) overflow.
                n_full = floor(Int, clamp(dm[j] / mx, 0.0, Float64(nlayers)))
                if n_full > 0
                    # Fresh-deposit distribution (from the depositional part of pickup)
                    pk_sum = 0.0
                    for k in 1:nf
                        pk_sum += max(0.0, -pickup[j, k])   # total deposit magnitude
                    end
                    if pk_sum > 0
                        d_dep = similar(ero)
                        for k in 1:nf
                            d_dep[k] = max(0.0, -pickup[j, k]) / pk_sum
                        end
                        # Fill `n_full` top layers with fresh deposit, shifting
                        # existing sediment down one slot per filled layer.
                        n_fill = min(n_full, nlayers)
                        for i in 1:n_fill
                            # Shift: move layer (i..nlayers-1) down by one slot.
                            # This means bed_mass[j, end, :] is discarded
                            # (or rather, it overflows off the bottom of the stack).
                            for il in nlayers:-1:(i+1)
                                for k in 1:nf
                                    m[j, il, k] = m[j, il-1, k]
                                end
                            end
                            # Fill the now-empty layer i with fresh deposit
                            # (of mass mx, using the deposit composition).
                            for k in 1:nf
                                m[j, i, k] = mx * d_dep[k]
                            end
                            # Subtract the filled mass from pickup (pickup was
                            # negative; we're removing deposit magnitude so it
                            # becomes less negative). Subtract `mx * d_dep[k]`
                            # from the MAGNITUDE of each depositional pickup
                            # entry.
                            for k in 1:nf
                                if pickup[j, k] < 0
                                    pickup[j, k] += mx * d_dep[k]
                                end
                            end
                        end
                        # If all layers got filled, discard any leftover deposit
                        if n_full ≥ nlayers
                            for k in 1:nf
                                pickup[j, k] = 0.0
                            end
                        end
                        # Recompute dm after the shifts
                        new_dm = 0.0
                        for k in 1:nf
                            new_dm -= pickup[j, k]
                        end
                        dm[j] = new_dm
                    end
                end
            end
        end
    end
    return nothing
end

function arrange_layers!(
    m::Array{Float64,3},
    dm::Vector{Float64},
    d::Array{Float64,3},
    jmax::Int,
    nlayers::Int,
    nf::Int,
)
    @inbounds for j in 1:jmax
        dmj = dm[j]
        if dmj < 0
            # Erosion: source is the deeper layer. Mass flows upward.
            for i in 2:nlayers
                for k in 1:nf
                    movement = dmj * d[j, i, k]
                    m[j, i-1, k] -= movement   # gains (subtracting negative)
                    m[j, i, k] += movement   # loses (adding negative)
                end
            end
            # The deepest layer has lost `|dm| * d[nlayers, :]` mass via the
            # last iteration. `update_bed_composition!` refills it from the
            # initial-grain-dist reservoir afterward.
        elseif dmj > 0
            # Deposition: source is the shallower layer. Mass flows down.
            for i in 2:nlayers
                for k in 1:nf
                    movement = dmj * d[j, i-1, k]
                    m[j, i-1, k] -= movement   # loses
                    m[j, i, k] += movement   # gains
                end
            end
            # The deepest layer has gained `dm * d[nlayers-1, :]` on the
            # last iteration. Drain that excess by subtracting `dm * d[nlayers, :]`
            # from it — prevents the deepest layer from inflating without bound.
            for k in 1:nf
                m[j, nlayers, k] -= dmj * d[j, nlayers, k]
            end
        end
    end
    return nothing
end

function update_bed_composition!(state::CshoreState, config::CshoreConfig,
    l::Int, dt::Float64)
    jmax_l = state.jmax[l]
    nf = size(state.bed_mass, 3)
    nlayers = size(state.bed_mass, 2)

    pickup = Matrix{Float64}(undef, jmax_l, nf)
    @inbounds for k in 1:nf, j in 1:jmax_l
        pickup[j, k] = state.pickup_fractions[j, k]
    end

    dm = Vector{Float64}(undef, jmax_l)
    @inbounds for j in 1:jmax_l
        s = 0.0
        for k in 1:nf
            s += pickup[j, k]
        end
        dm[j] = -s
    end

    m = view(state.bed_mass, 1:jmax_l, :, :)

    m_arr = state.bed_mass    # full (nn, nlayers, nf)
    layer_h = config.multifraction.layer_thickness
    geom_guard_on = get(ENV, "CSHORE_DISABLE_GEOMETRIC_GUARD", "0") != "1"
    @inbounds for j in 1:jmax_l
        geom_guard_on || break
        zh = state.zb_hard[j]
        if zh > -1e29 && (state.zb[j, l] - zh) < layer_h
            for k in 1:nf
                if pickup[j, k] > 0.0
                    avail = 0.0
                    for il in 1:nlayers
                        avail += m_arr[j, il, k]
                    end
                    if pickup[j, k] > avail
                        excess = pickup[j, k] - avail
                        pickup[j, k] = avail
                        dm[j] += excess   # less negative → smaller zb drop
                    end
                end
            end
        end
    end

    prevent_negative_mass!(m_arr, dm, pickup, jmax_l, nlayers, nf)

    # d[j, layer, k] = bed_mass[j, layer, k] / sum_k(bed_mass[j, layer, :])
    d = zeros(jmax_l, nlayers, nf)
    @inbounds for j in 1:jmax_l, il in 1:nlayers
        s = 0.0
        for k in 1:nf
            s += m_arr[j, il, k]
        end
        if s > 0
            for k in 1:nf
                d[j, il, k] = m_arr[j, il, k] / s
            end
        end
    end

    @inbounds for k in 1:nf, j in 1:jmax_l
        m_arr[j, 1, k] -= pickup[j, k]
    end
    if MASS_DEBUG_ON[]
        s5 = 0.0
        @inbounds for k in 1:nf, j in 1:jmax_l
            s5 += pickup[j, k]
        end
        MASS_DEBUG_COUNTERS.m_step5_active_removed += s5
    end

    m_sub = Array{Float64,3}(undef, jmax_l, nlayers, nf)
    @inbounds for k in 1:nf, il in 1:nlayers, j in 1:jmax_l
        m_sub[j, il, k] = m_arr[j, il, k]
    end
    arrange_layers!(m_sub, dm, d, jmax_l, nlayers, nf)
    @inbounds for k in 1:nf, il in 1:nlayers, j in 1:jmax_l
        m_arr[j, il, k] = m_sub[j, il, k]
    end

    init_frac = config.multifraction.initial_fractions
    @inbounds for j in 1:jmax_l
        if dm[j] < 0   # erosion cell
            for k in 1:nf
                # `-dm * initial_frac` is positive → adds mass to deepest layer
                m_arr[j, nlayers, k] -= dm[j] * init_frac[k]
            end
            if MASS_DEBUG_ON[]
                MASS_DEBUG_COUNTERS.m_step7_reservoir += -dm[j]
            end
        end
    end

    if MASS_DEBUG_ON[]
        layer_h = config.multifraction.layer_thickness
        thin = 0
        @inbounds for j in 1:jmax_l
            zh = state.zb_hard[j]
            if zh > -1e29 && (state.zb[j, l] - zh) < layer_h
                thin += 1
            end
        end
        MASS_DEBUG_COUNTERS.n_thin_cover_nodes += thin
    end

    @inbounds for j in 1:jmax_l
        for k in 1:nf
            for il in 1:nlayers
                if m_arr[j, il, k] < 0.0
                    deficit = -m_arr[j, il, k]
                    if MASS_DEBUG_ON[]
                        MASS_DEBUG_COUNTERS.n_clamp_fired += 1
                    end
                    m_arr[j, il, k] = 0.0

                    for il2 in (il+1):nlayers
                        avail = m_arr[j, il2, k]
                        if avail > 0.0
                            take = min(avail, deficit)
                            m_arr[j, il2, k] -= take
                            deficit -= take
                            if deficit ≤ 0.0
                                break
                            end
                        end
                    end

                    if deficit > 0.0
                        for k2 in 1:nf
                            for il2 in 1:nlayers
                                if k2 == k && il2 == il
                                    continue
                                end
                                avail = m_arr[j, il2, k2]
                                if avail > 0.0
                                    take = min(avail, deficit)
                                    m_arr[j, il2, k2] -= take
                                    deficit -= take
                                    if deficit ≤ 0.0
                                        break
                                    end
                                end
                            end
                            if deficit ≤ 0.0
                                break
                            end
                        end
                    end

                    if deficit > 0.0
                        if get(ENV, "CSHORE_DISABLE_CLAMP_PASS3", "0") != "1"
                            dm[j] += deficit
                        end
                        if MASS_DEBUG_ON[]
                            MASS_DEBUG_COUNTERS.m_clamp_unrecovered += deficit
                        end
                    end
                end
            end
        end
    end
    _clip_tiny_negatives!(m_arr, jmax_l, nlayers, nf)

    ρs = config.sediment.sg * 1000.0
    one_minus_n = 1.0 - config.multifraction.porosity
    mass_per_m = ρs * one_minus_n
    max_dz = config.max_dzb_per_step
    @inbounds for j in 1:jmax_l
        dz = dm[j] / mass_per_m
        if isfinite(max_dz) && max_dz > 0.0 && abs(dz) > max_dz
            # Clamp dz to ±max_dz and back-propagate the clamp into dm[j]
            # so the bed_mass bookkeeping stays consistent with zb.
            dz_clamped = sign(dz) * max_dz
            dm_adjust = (dz - dz_clamped) * mass_per_m
            # Push the unrealizable mass back into the deepest layer
            # uniformly across fractions (keeps total node mass invariant).
            nf_total = size(state.bed_mass, 3)
            for k in 1:nf_total
                state.bed_mass[j, nlayers, k] += dm_adjust * init_frac[k]
            end
            dz = dz_clamped
        end
        state.delzb[j, l] = dz
        state.zb[j, l] += dz
    end

    if config.morph_diffusion > 0.0 && jmax_l >= 3
        dx_sq = config.grid.dx * config.grid.dx
        coef = config.morph_diffusion * dt / dx_sq
        coef = min(coef, 0.45)
        delz_diff = Vector{Float64}(undef, jmax_l)
        delz_diff[1] = 0.0
        delz_diff[jmax_l] = 0.0
        @inbounds for j in 2:(jmax_l-1)
            delz_diff[j] = coef * (state.zb[j-1, l] - 2.0 * state.zb[j, l] + state.zb[j+1, l])
        end
        @inbounds for j in 2:(jmax_l-1)
            state.zb[j, l] += delz_diff[j]
            state.delzb[j, l] += delz_diff[j]
        end
    end

    @inbounds for j in 1:jmax_l
        state.vbx[j, l] += dt * state.qbx[j, 1]
        state.vsx[j, l] += dt * state.qsx[j, 1]
    end

    return nothing
end

function _clip_tiny_negatives!(m::Array{Float64,3}, jmax_l::Int,
    nlayers::Int, nf::Int; tol_rel::Float64=1e-10)
    pos_sum = 0.0
    pos_cnt = 0
    @inbounds for k in 1:nf, il in 1:nlayers, j in 1:jmax_l
        v = m[j, il, k]
        if v > 0
            pos_sum += v
            pos_cnt += 1
        end
    end
    mean_pos = pos_cnt > 0 ? pos_sum / pos_cnt : 1.0
    thresh = tol_rel * mean_pos
    @inbounds for k in 1:nf, il in 1:nlayers, j in 1:jmax_l
        if m[j, il, k] < 0 && -m[j, il, k] < thresh
            m[j, il, k] = 0.0
        end
    end
    return nothing
end

function bed_elevation_from_mass(state::CshoreState, config::CshoreConfig,
    j::Int, l::Int=1)
    ρs = config.sediment.sg * 1000.0
    one_minus_n = 1.0 - config.multifraction.porosity
    total = 0.0
    nlayers = size(state.bed_mass, 2)
    nf = size(state.bed_mass, 3)
    @inbounds for ilay in 1:nlayers, k in 1:nf
        total += state.bed_mass[j, ilay, k]
    end
    return total / (ρs * one_minus_n)
end

function active_layer_fractions!(state::CshoreState, config::CshoreConfig, l::Int)
    mf = config.multifraction
    nf = size(state.bed_mass, 3)
    nlayers = size(state.bed_mass, 2)
    jmax_l = state.jmax[l]
    layer_h = mf.layer_thickness
    facDOD = mf.facDOD

    # Single-grain shortcut — degenerate to 1.0.
    if nf ≤ 1
        @inbounds for j in 1:jmax_l
            state.active_frac[j, 1] = 1.0
        end
        return nothing
    end

    @inbounds for j in 1:jmax_l
        DOD = facDOD * state.hrms[j]
        if DOD ≤ 0.0
            # Collapse to layer 1
            s = 0.0
            for k in 1:nf
                s += state.bed_mass[j, 1, k]
            end
            if s > 0
                for k in 1:nf
                    state.active_frac[j, k] = state.bed_mass[j, 1, k] / s
                end
            else
                for k in 1:nf
                    state.active_frac[j, k] = mf.initial_fractions[k]
                end
            end
            continue
        end

        n_full = isfinite(DOD) ? floor(Int, DOD / layer_h) : 0
        frac = (DOD - n_full * layer_h) / layer_h   # ∈ [0,1)
        # Clamp into storage.
        if n_full ≥ nlayers
            n_full = nlayers
            frac = 0.0
        end

        # Weighted mass over the virtual active layer.
        total_mass = 0.0
        for il in 1:n_full
            for k in 1:nf
                total_mass += state.bed_mass[j, il, k]
            end
        end
        if frac > 0 && (n_full + 1) ≤ nlayers
            for k in 1:nf
                total_mass += frac * state.bed_mass[j, n_full+1, k]
            end
        end

        if total_mass > 0
            for k in 1:nf
                s = 0.0
                for il in 1:n_full
                    s += state.bed_mass[j, il, k]
                end
                if frac > 0 && (n_full + 1) ≤ nlayers
                    s += frac * state.bed_mass[j, n_full+1, k]
                end
                state.active_frac[j, k] = s / total_mass
            end
        else
            # Empty stack — fall back to initial distribution.
            for k in 1:nf
                state.active_frac[j, k] = mf.initial_fractions[k]
            end
        end
    end
    return nothing
end

function smooth_active_composition!(state::CshoreState, config::CshoreConfig, l::Int)
    nf = size(state.active_frac, 2)
    nf ≤ 1 && return nothing
    jmax_l = state.jmax[l]
    jmax_l ≥ 3 || return nothing
    n_passes = config.multifraction.n_composition_smooth
    n_passes ≤ 0 && return nothing

    @inbounds for k in 1:nf
        for _ in 1:n_passes
            smooth_tridiagonal!(view(state.active_frac, 1:jmax_l, k))
        end
    end

    @inbounds for j in 1:jmax_l
        s = 0.0
        for k in 1:nf
            v = state.active_frac[j, k]
            if v < 0.0
                state.active_frac[j, k] = 0.0
            else
                s += v
            end
        end
        if s > 0
            inv_s = 1.0 / s
            for k in 1:nf
                state.active_frac[j, k] *= inv_s
            end
        end
    end
    return nothing
end

function mix_top_layer!(state::CshoreState, config::CshoreConfig, l::Int)
    mf = config.multifraction
    mf.process_mixtoplayer || return nothing
    nf = size(state.bed_mass, 3)
    nf ≤ 1 && return nothing   # single-grain — nothing to mix
    nlayers = size(state.bed_mass, 2)
    nlayers ≤ 1 && return nothing

    jmax_l = state.jmax[l]
    layer_h = mf.layer_thickness
    total_stack_h = nlayers * layer_h
    facDOD = mf.facDOD

    @inbounds for j in 1:jmax_l
        # Depth of disturbance = facDOD × Hrms at this node
        DOD = facDOD * state.hrms[j]
        DOD > 0 || continue

        n_full = isfinite(DOD) ? floor(Int, DOD / layer_h) : 0
        if n_full == 0
            n_full = 1   # at least mix the top layer
        end
        n_mix = min(n_full, nlayers)

        f = min(1.0, total_stack_h / DOD)

        # Compute per-fraction mean over top n_mix layers
        for k in 1:nf
            s = 0.0
            for il in 1:n_mix
                s += state.bed_mass[j, il, k]
            end
            mean_k = s / n_mix
            # Write the blend back to each mixed layer
            for il in 1:n_mix
                old = state.bed_mass[j, il, k]
                state.bed_mass[j, il, k] = mean_k * f + old * (1.0 - f)
            end
        end
    end
    return nothing
end
