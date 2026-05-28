# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
exner.jl — Exner equation / morphodynamic bed update.

Adaptive timestep uses `config.morph_courant` (default 0.3) as the
Courant-like safety factor in Δt = morph_courant · Δx / |cb_max|. Tune
lower for tighter stability at the cost of run time.
==============================================================================#

Base.@kwdef struct ExnerStepResult
    delt::Float64
    iend::Bool
end

function compute_timestep!(state::CshoreState, config::CshoreConfig, l::Int,
    t_window_end::Float64)
    state.supply_factor_applied ||
        error("supply factor was not applied before compute_timestep!. ")

    jmax_l = state.jmax[l]
    dx = config.grid.dx
    dx2 = 2.0 * dx

    q = view(state.q_total, 1:jmax_l)
    zb = view(state.zb, 1:jmax_l, l)

    cb = zeros(Float64, jmax_l)
    dzbmax = 0.1 * dx
    cbmax = 0.004
    @inbounds for j in 1:jmax_l
        if j == 1
            delq = q[2] - q[1]
            delzb1 = zb[2] - zb[1]
        elseif j == jmax_l
            delq = q[j] - q[j-1]
            delzb1 = zb[j] - zb[j-1]
        else
            delq = q[j+1] - q[j-1]
            delzb1 = zb[j+1] - zb[j-1]
        end
        cb[j] = abs(delzb1) > dzbmax ? delq / delzb1 : 0.0
        dumc = abs(cb[j])
        if dumc > cbmax
            cbmax = dumc
        end
    end

    # ---- Adaptive dt --------------------------------------------------------
    delt = config.morph_courant * dx / cbmax
    iend = false
    if state.time + delt ≥ t_window_end
        delt = t_window_end - state.time
        iend = true
    end
    state.delt = delt

    ρs = config.sediment.sg * 1000.0
    one_minus_n = 1.0 - config.multifraction.porosity
    mass_per_m = ρs * one_minus_n

    nf = size(state.qbx, 2)
    fill!(view(state.pickup_fractions, 1:jmax_l, :), 0.0)
    if nf ≥ 1
        sporo1 = config.sediment.sporo1

        n_faces = jmax_l - 1
        q_face = zeros(Float64, n_faces, nf)

        decayl_m = min(state.xb[state.jswl[l]] / 4.0,
            2.0 * state.wt[1] * state.cp[1])
        # Guard against NaN/Inf in cp[1] on fine grids (see transport.jl).
        jdecay = max(1, round(Int, (isfinite(decayl_m) ? decayl_m : 2.0) / dx))
        jr_l = state.jr

        f_lo = jdecay
        f_hi_base = jr_l - 1
        if state.jdry > jr_l && state.jwd > 0
            f_hi_base = state.jdry - 1
        end
        f_hi = min(n_faces, f_hi_base)
        @inbounds for k in 1:nf
            for f in f_lo:f_hi
                q_here = (state.qbx[f, k] + state.qsx[f, k]) / sporo1
                q_next = (state.qbx[f+1, k] + state.qsx[f+1, k]) / sporo1
                q_face[f, k] = 0.5 * (q_here + q_next)
            end
        end

        if f_lo > 1 && f_lo <= n_faces
            if config.options.ilab == 0
                # Field: open seaward BC
                @inbounds for k in 1:nf
                    q_edge = q_face[f_lo, k]
                    for f in 1:(f_lo-1)
                        q_face[f, k] = q_edge
                    end
                end
            end
            # Lab (ilab=1): faces 1..f_lo-1 remain zero — closed BC.
        end


        face_smoother_on = get(ENV, "CSHORE_DISABLE_FACE_SMOOTHER", "0") != "1"
        n_face_smooth = config.multifraction.n_face_flux_smooth
        if face_smoother_on && n_face_smooth > 0 && nf >= 1 && n_faces >= 3
            tmp = Vector{Float64}(undef, n_faces)
            pw_thresh = 0.01
            face_wet = Vector{Bool}(undef, n_faces)
            @inbounds for f in 1:n_faces
                pw_l = (f >= 1) ? state.pwet[f] : 0.0
                pw_r = (f + 1 <= jmax_l) ? state.pwet[f+1] : 0.0
                face_wet[f] = (max(pw_l, pw_r) > pw_thresh)
            end
            face_is_landward_bc = falses(n_faces)
            @inbounds for f in (f_hi+1):n_faces
                face_is_landward_bc[f] = true
            end

            @inbounds for k in 1:nf
                for pass in 1:n_face_smooth
                    for f in 1:n_faces
                        tmp[f] = q_face[f, k]
                    end
                    for f in 2:(n_faces-1)
                        face_wet[f] || continue
                        face_is_landward_bc[f] && continue

                        q_l = (face_wet[f-1] && !face_is_landward_bc[f-1]) ?
                              tmp[f-1] : tmp[f]
                        q_r = (face_wet[f+1] && !face_is_landward_bc[f+1]) ?
                              tmp[f+1] : tmp[f]
                        q_face[f, k] = 0.25 * q_l + 0.5 * tmp[f] + 0.25 * q_r
                    end
                end
            end
        end

        if config.options.isedav != 0
            # Accumulate outflow (in kg/m, per cell) over both faces.
            # The limit is on total bed mass above the hardbottom, not
            # per fraction; fractions are then scaled uniformly (composition preserved).
            cell_out = zeros(Float64, jmax_l)
            @inbounds for k in 1:nf
                for f in 1:n_faces
                    q = q_face[f, k]
                    if q > 0
                        # Outflow from cell f
                        cell_out[f] += q
                    elseif q < 0
                        # Outflow from cell f+1
                        cell_out[f+1] += -q
                    end
                end
            end
            # `cell_out[j]` is in m²/s (outflow transport rate).
            # Convert to mass leaving during `delt`: mass = q·dt·mass_per_m.
            # Available sand: hp[j] * mass_per_m.
            @inbounds for j in 1:jmax_l
                hp_j = state.hp[j, l]
                if hp_j == Inf
                    continue       # no constraint
                end
                mass_out = cell_out[j] * delt * mass_per_m
                mass_avail = hp_j * mass_per_m
                if mass_out > mass_avail && mass_out > 0
                    scale = mass_avail / mass_out
                    # Scale every OUTFLOWING face flux from this cell.
                    for k in 1:nf
                        if j <= n_faces
                            # Right face of cell j is face index j (f=j)
                            if q_face[j, k] > 0
                                q_face[j, k] *= scale
                            end
                        end
                        if j >= 2
                            # Left face of cell j is face index j-1
                            if q_face[j-1, k] < 0
                                q_face[j-1, k] *= scale
                            end
                        end
                    end
                end
            end
        end

        @inbounds for k in 1:nf
            # Boundary cells: Neumann (zero divergence)
            state.pickup_fractions[1, k] = 0.0
            state.pickup_fractions[jmax_l, k] = 0.0
            # Interior cells: standard two-face divergence
            for j in 2:(jmax_l-1)
                div_j = (q_face[j, k] - q_face[j-1, k]) / dx
                state.pickup_fractions[j, k] = div_j * delt * mass_per_m
            end
        end

        pickup_smoother_on = get(ENV, "CSHORE_DISABLE_PICKUP_SMOOTHER", "0") != "1"
        n_pickup_smooth = config.multifraction.n_pickup_smooth
        if pickup_smoother_on && n_pickup_smooth > 0 && jmax_l >= 3
            tmp2 = Vector{Float64}(undef, jmax_l)
            pw_thresh = 0.01

            hp_full = 1e-2          # 1 cm — full smoothing above this
            hp_none = 1e-3          # 1 mm — no smoothing below this
            w = Vector{Float64}(undef, jmax_l)

            @inbounds for j in 1:jmax_l
                if config.options.isedav == 0
                    w[j] = 1.0
                    continue
                end
                hp_j = state.hp[j, l]
                if !isfinite(hp_j) || hp_j >= hp_full
                    w[j] = 1.0
                elseif hp_j <= hp_none
                    w[j] = 0.0
                else
                    w[j] = (hp_j - hp_none) / (hp_full - hp_none)
                end
            end

            for pass in 1:n_pickup_smooth
                @inbounds for k in 1:nf
                    for j in 1:jmax_l
                        tmp2[j] = state.pickup_fractions[j, k]
                    end
                    for j in 2:(jmax_l-1)
                        # Skip entirely-no-smooth cells. For partial-w
                        # cells we still apply the kernel — just with
                        # face-shared weights.
                        w[j] > 0.0 || continue
                        # Face-shared kernel weights — guarantees that
                        # the mass cell j sends across its left face
                        # equals the mass cell j-1 receives. Σ pickup
                        # is preserved under varying w_j.
                        wL = min(w[j], w[j-1])
                        wR = min(w[j], w[j+1])
                        # Kernel: (1 − 0.25 wL − 0.25 wR) on self,
                        # 0.25 wL on left, 0.25 wR on right. Sums to 1.
                        state.pickup_fractions[j, k] =
                            (1.0 - 0.25 * wL - 0.25 * wR) * tmp2[j] +
                            0.25 * wL * tmp2[j-1] +
                            0.25 * wR * tmp2[j+1]
                    end
                end
                # Catch the residual after each pass.
                if config.options.isedav != 0 &&
                   get(ENV, "CSHORE_DISABLE_HARDBOTTOM_REDISTRIBUTE", "0") != "1"
                    _redistribute_pickup_for_hardbottom!(state, config, l, delt,
                        mass_per_m, jmax_l, nf)
                end
            end
        end

        eps_pk = 1e-30
        min_cells_for_adjust = 50
        @inbounds for k in 1:nf
            s = 0.0
            cnt = 0
            for j in 2:(jmax_l-1)
                if abs(state.pickup_fractions[j, k]) < eps_pk
                    continue
                end
                s += state.pickup_fractions[j, k]
                cnt += 1
            end
            if cnt >= min_cells_for_adjust
                mean_k = s / cnt
                if mean_k != 0.0
                    for j in 2:(jmax_l-1)
                        if abs(state.pickup_fractions[j, k]) < eps_pk
                            continue
                        end
                        state.pickup_fractions[j, k] -= mean_k
                    end
                end
            end
        end


        if config.options.isedav != 0 &&
           get(ENV, "CSHORE_DISABLE_HARDBOTTOM_REDISTRIBUTE", "0") != "1"
            _redistribute_pickup_for_hardbottom!(state, config, l, delt,
                mass_per_m, jmax_l, nf)

            @inbounds for k in 1:nf
                s2 = 0.0
                cnt2 = 0
                for j in 2:(jmax_l-1)
                    if abs(state.pickup_fractions[j, k]) < eps_pk
                        continue
                    end
                    s2 += state.pickup_fractions[j, k]
                    cnt2 += 1
                end
                if cnt2 >= min_cells_for_adjust && s2 != 0.0
                    mean2_k = s2 / cnt2
                    for j in 2:(jmax_l-1)
                        if abs(state.pickup_fractions[j, k]) < eps_pk
                            continue
                        end
                        state.pickup_fractions[j, k] -= mean2_k
                    end
                end
            end
        end
    end

    state.supply_factor_applied = false

    return ExnerStepResult(delt=delt, iend=iend)
end

function _redistribute_pickup_for_hardbottom!(state::CshoreState,
    config::CshoreConfig,
    l::Int, delt::Float64,
    mass_per_m::Float64,
    jmax_l::Int, nf::Int)
    pw_thresh = 0.01
    max_passes = 3
    eps_total = 1e-30   # numeric noise floor
    residual_total = 0.0
    n_residual_cells = 0

    for _pass in 1:max_passes
        any_violation = false

        @inbounds for j in 1:jmax_l
            hp_j = state.hp[j, l]
            isfinite(hp_j) || continue                  # hp = Inf → no constraint
            max_pickup_j = hp_j * mass_per_m / delt

            total_j = 0.0
            sum_pos = 0.0
            for k in 1:nf
                pjk = state.pickup_fractions[j, k]
                total_j += pjk
                if pjk > 0
                    sum_pos += pjk
                end
            end
            (total_j <= max_pickup_j) && continue       # under budget — fine
            (sum_pos <= eps_total) && continue       # all-deposition cell, can't reduce
            any_violation = true

            sum_neg = total_j - sum_pos
            new_sum_pos = max_pickup_j - sum_neg
            new_sum_pos = max(new_sum_pos, 0.0)
            scale_pos = new_sum_pos / sum_pos    # in [0, 1)
            spilled_total = sum_pos - new_sum_pos    # > 0

            original_pickup_j = zeros(Float64, nf)
            for k in 1:nf
                original_pickup_j[k] = state.pickup_fractions[j, k]
            end

            spilled = zeros(Float64, nf)
            for k in 1:nf
                pjk = state.pickup_fractions[j, k]
                if pjk > 0
                    new_pjk = pjk * scale_pos
                    spilled[k] = pjk - new_pjk
                    state.pickup_fractions[j, k] = new_pjk
                end
            end

            for offset in (-1, +1)
                jn = j + offset
                (1 <= jn <= jmax_l) || continue
                state.pwet[jn] > pw_thresh || continue

                hp_n = state.hp[jn, l]
                if isfinite(hp_n)
                    total_n = 0.0
                    for k in 1:nf
                        total_n += state.pickup_fractions[jn, k]
                    end
                    cap_n = (hp_n * mass_per_m / delt) - total_n
                    cap_n <= 0 && continue   # neighbour also full
                else
                    cap_n = Inf
                end

                spilled_remaining = sum(spilled)
                spilled_remaining > eps_total || break
                accepted = min(cap_n, spilled_remaining)
                share = accepted / spilled_remaining
                for k in 1:nf
                    if spilled[k] > 0
                        delta = spilled[k] * share
                        state.pickup_fractions[jn, k] += delta
                        spilled[k] -= delta
                    end
                end
            end

            spilled_left = 0.0
            for k in 1:nf
                spilled_left += spilled[k]
            end
            if spilled_left > eps_total
                for k in 1:nf
                    state.pickup_fractions[j, k] = original_pickup_j[k]
                end
                residual_total += spilled_left
                n_residual_cells += 1
            end
        end

        any_violation || break
    end

    if n_residual_cells > 0 && get(ENV, "CSHORE_VERBOSE_REDISTRIBUTE", "0") == "1"
        approx_mass_kg = residual_total * delt * config.grid.dx
        @info "redistribute: rescue path (original pickup restored — clamp will adjust zb)" n_cells = n_residual_cells approx_residual_mass_kg = approx_mass_kg
    end
    return state
end

function apply_hardbottom_clamp!(state::CshoreState, config::CshoreConfig, l::Int)
    config.options.isedav == 0 && return state
    jmax_l = state.jmax[l]
    warn_threshold = 1e-5     # 10 µm
    big_threshold = 1e-3
    n_clamped = 0
    max_d = 0.0
    sum_d = 0.0

    @inbounds for j in 1:jmax_l
        zh = state.zb_hard[j, l]
        if zh == -Inf
            state.hp[j, l] = Inf
            continue
        end
        d = state.zb[j, l] - zh
        if d < 0
            depth_below = -d
            n_clamped += 1
            sum_d += depth_below
            if depth_below > max_d
                max_d = depth_below
            end
            # Only emit a per-event warning for unusually large overshoots.
            if depth_below > big_threshold
                @warn "apply_hardbottom_clamp!: zb dropped $(depth_below) m below zb_hard at j=$j (limiter missed a case?)"
            end
            state.zb[j, l] = zh
            d = 0.0
        end
        state.hp[j, l] = d
    end

    if n_clamped > 0 && max_d > warn_threshold &&
       get(ENV, "CSHORE_VERBOSE_CLAMP", "0") == "1"
        ρs = config.sediment.sg * 1000.0
        sand_frac = 1.0 - config.multifraction.porosity
        mass_inj = sum_d * config.grid.dx * ρs * sand_frac
        @info "apply_hardbottom_clamp! summary" n_clamped max_d_m = max_d sum_d_m = sum_d approx_mass_injected_kg_per_m = mass_inj
    end
    return state
end
