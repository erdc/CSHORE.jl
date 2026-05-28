# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
wetdry.jl — Wet/dry swash zone dynamics.

Computes swash zone hydrodynamics: probability of wet (PWET), swash depth
(HWD), swash velocity (USWD/UMEAWD), and overtopping flux (QOTF).

Scope: IPERM=0,1; INFILT=0,1 (mutually exclusive); no ponding (IPOND=0);
no landward transmission (IWTRAN=0).  INFILT=1 models infiltration of
overtopping water into a sand bed landward of the dune crest, reducing
swash momentum and cumulative water flux.
==============================================================================#

"""
    probwd(pw, a, us, uc) -> Float64

Swash-zone probability of sediment motion. Uses an exponential decay scaled
by the wetting probability `pw`:

- `|us| ≤ uc`: P = pw · exp(-a·(uc-us)²)  (sub-critical: exponential decay)
- `us > uc`:   P = pw                       (super-critical: full mobilization)
- `us < -uc`:  P = pw · (1 - exp(-a·(uc+us)²) + exp(-a·(uc-us)²))

Arguments:
- `pw`: wetting probability (PWET[j])
- `a`:  PWET / (AGWD · g · H) — dimensionless swash intensity
- `us`: swash velocity (USWD[j])
- `uc`: critical velocity for motion (bedload: √(gsd50s/fb2), suspended: wf/fb2^⅓)
"""
@inline function probwd(pw::Float64, a::Float64, us::Float64, uc::Float64)
    if abs(us) <= uc
        return pw * exp(-a * (uc - us)^2)
    elseif us > uc
        return pw
    else
        return pw * (1.0 - exp(-a * (uc + us)^2) + exp(-a * (uc - us)^2))
    end
end

"""
    gbwd(r) -> Float64

Bottom shear stress function for wet-dry zone.
"""
function gbwd(r::Float64)
    if r ≥ 0.0
        return 1.0 + 1.77245 * r + r * r
    else
        r2 = r * r
        return 2.0 * exp(-r2) - r2 - 1.0 + 1.77245 * r * (3.0 - 2.0 * erfcc(r))
    end
end

"""
    gdwd(r, x) -> Float64

Vegetation drag function for wet-dry zone.
"""
function gdwd(r::Float64, x::Float64)
    ex = exp(-x)
    sx = sqrt(max(x, 0.0))
    fx = 1.0 - erfcc(sx)
    r2 = r * r
    if r ≥ 0.0
        return 2.0 - (x + 2.0) * ex + r * (1.77245 * x - 3.0 * sx * ex +
               1.77245 * (1.5 - x) * fx) + r2 * (1.0 - ex)
    else
        er2 = exp(-r2)
        fr = 1.0 - erfcc(r)
        c = (x + 2.0 + r2 + 3.0 * r * sx) * ex
        if x ≤ r2
            return 2.0 * x * er2 - 2.0 - r2 + 1.77245 * r * ((x - 1.5) * fx +
                   2.0 * x * fr + x) + c
        else
            return 4.0 * er2 - 2.0 - r2 + 5.317362 * r * fr + 1.77245 * r * (x +
                   (1.5 - x) * fx) - c
        end
    end
end

"""
    tranwd!(f1, jr, f2, js, je)

Connect wave-field vector `f1[1:jr]` with swash-field vector `f2[js:je]`.
Overwrites `f1[jr+1:je]` (or blends the overlap if `jr ≥ js`).
"""
function tranwd!(f1::AbstractVector, jr::Int, f2::AbstractVector, js::Int, je::Int)
    if jr ≥ js
        @inbounds for j in js:jr
            f1[j] = 0.5 * (f1[j] + f2[j])
        end
        @inbounds for j in (jr + 1):je
            f1[j] = f2[j]
        end
    else
        dsr = Float64(js - jr)
        @inbounds for j in (jr + 1):js
            dum = Float64(js - j) / dsr
            f1[j] = f1[jr] * dum + f2[js] * (1.0 - dum)
        end
        @inbounds for j in (js + 1):je
            f1[j] = f2[j]
        end
    end
    return nothing
end

"""
    wetdry!(state, config, itime, l, iteqo)

Compute swash-zone hydrodynamics from `JWD` to `JDRY`.

Sets: `state.pwet`, `state.hwd`, `state.uswd`, `state.umeawd`,
`state.ustdwd`, `state.sigwd`, `state.qotf`, `state.jwd`, `state.jdry`.
"""
function wetdry!(state::CshoreState, config::CshoreConfig, itime::Int, l::Int, iteqo::Int)
    dx     = config.grid.dx
    gamma  = config.gamma
    iveg   = config.options.iveg
    jmax_l = state.jmax[l]

    # Swash parameters (from SwashConfig or defaults)
    oc = config.swash
    awd  = oc.awd  > 0 ? oc.awd  : 1.6    # default for erodible profiles
    ewd  = oc.ewd  > 0 ? oc.ewd  : 0.01
    cwd  = 0.75 * sqrt(π) * awd
    aqwd = cwd * awd
    bwd  = (2.0 - 9.0 * π / 16.0) * awd^2 + 1.0
    agwd = awd^2
    auwd = 0.5 * sqrt(π) * awd
    # Maximum seepage velocity (WPM) for infiltration / porous flow.
    # When INFILT=1 we compute it from sand properties (d50, porosity);
    # when IPERM=1 the PorousInput controls it; otherwise fall back to any
    # value that may already be on the OvertoppingConfig.
    wpm = if config.options.infilt == 1
        # Forchheimer coefficients for sand
        d50_sand = config.sediment.d50
        poro = 1.0 - config.sediment.sporo1   # sand porosity (sporo1 = 1-n)
        C_s  = 1.0 - poro
        wnu  = 1.0e-6
        α  = 1000.0 * wnu * C_s^2 / (poro * d50_sand)^2
        β1 = 5.0 * C_s / poro^3 / d50_sand
        alsta  = α  / GRAV
        besta1 = β1 / GRAV
        alsta2 = alsta * alsta
        be2    = 2.0 * besta1
        be4    = 2.0 * be2
        (sqrt(alsta2 + be4) - alsta) / be2
    else
        oc.wpm > 0 ? oc.wpm : 0.0
    end

    # Expose the resolved wpm for use by step_groundwater! (infilt=1 coupling).
    # When infilt=0, this stays at 0.0 (no seepage).
    state.wpm_derived = wpm

    jcrest = state.jcrest[l]
    hwdmin = 0.001  # minimum swash depth (m) — used for state.h floors,
                    # NOT for the march-termination threshold below (which
                    # uses 1e-6 to match the OLD reference implementation
                    # and let swash march extend further up the bluff face).

    # ---- Initialize at wet-dry transition ----
    if iteqo ≤ 2
        state.jwd = state.jswl[l]
        if config.options.iprofl == 1 && state.jwd > state.jr
            state.jwd = state.jr
            state.jdry = state.jr
            state.h1 = state.h[state.jr]
            _finalize_qotf!(state, config, l, aqwd, jcrest)
            return nothing
        end
        if state.jwd > state.jr
            state.jwd = state.jr
        end
        state.h1 = state.h[state.jwd]
    end

    jwd = state.jwd
    h1  = state.h1
    h1 = max(h1, 1e-6)
    # IG Layer-3 coupling: augment swash-front depth by c_ig_swash · hrms_ig.
    # Reads state.hrms_ig from the previous substep (zero on cold start);
    # compute_ig_field! refreshes it at the end of each substep. The lag is
    # acceptable because IG amplitude changes slowly relative to the
    # adaptive morphodynamic Δt. No-op when config.ig === nothing.
    if config.ig !== nothing && config.ig.c_ig_swash > 0.0 && jwd <= length(state.hrms_ig)
        h1 += config.ig.c_ig_swash * state.hrms_ig[jwd]
    end
    state.hwd[jwd] = h1
    bgh3 = bwd * GRAV * h1^3

    # Surf similarity correction factor
    ssp_50 = 1.0
    a_corr = 1.0
    correct = gamma / SQR8
    correct = 0.5 * (1.0 + correct) + 0.5 * (1.0 - correct) * tanh(a_corr * (state.ssp - ssp_50))

    state.sigwd[jwd] = correct * h1
    pmg1 = awd / sqrt(GRAV)

    state.pwet[jwd] = 1.0
    qwx = state.qo[l]  # IPERM=0 → no QP subtraction
    # INFILT=1: QP starts accumulating past the crest; zero at/before jwd
    if config.options.infilt == 1
        state.qp[jwd] = 0.0
    end

    qs = qwx - aqwd * h1 * sqrt(GRAV * h1)
    if qs > 0.0; qs = 0.0; end
    state.uswd[jwd] = qs / h1
    state.umeawd[jwd] = auwd * sqrt(GRAV * h1) + state.uswd[jwd]
    dum = agwd * GRAV * h1 - (state.umeawd[jwd] - state.uswd[jwd])^2
    if dum < 0.0
        state.jdry = jwd
        _finalize_qotf!(state, config, l, aqwd, jcrest)
        return nothing
    end
    state.ustdwd[jwd] = sqrt(dum)

    a_val = qwx^2 / bgh3
    a1 = a_val

    # Empirical wet probability parameter
    wdn = 1.01 + 0.98 * tanh((state.qo[l]^2 / bgh3)^0.3)
    w1  = wdn - 1.0
    bnwd = bwd * (1.0 + a1) * (2.0 - wdn) / max(w1, 1e-12)

    # ---- Landward march ----
    jend = jmax_l - 1
    lstart = 1
    g_arr  = zeros(jmax_l)
    dg_arr = zeros(jmax_l)
    h2 = h1
    d_val = 0.0    # accumulation variable for upslope
    ah = agwd / max(h1, 1e-12)
    bn12 = 0.0

    # Crest and peak tracking
    jc_peak = jcrest
    hc_peak = h1
    pc_peak = 1.0
    qwc = qwx

    # Downslope variables (hoisted for scope)
    pci = 1.0
    qwc2 = qwc^2
    bg = bwd * GRAV
    cpc = 0.0
    ab_ds = 0.0

    if jwd > jend
        state.jdry = state.jr
        _finalize_qotf!(state, config, l, aqwd, jcrest)
        return nothing
    end

    for j in jwd:jend
        jp1 = j + 1

        # Determine if we're on the upslope (toward crest) or downslope
        iupslp = false
        jdum_crest = jcrest
        if j < jdum_crest && state.zb[jp1, l] ≥ state.zb[j, l]
            iupslp = true
        end
        if j == jwd
            iupslp = true
        end

        if iupslp
            # ---- UPSLOPE MARCH ----
            if lstart == 1
                h2 = state.hwd[j]
                bn12 = bnwd * (h1 / max(h2, 1e-12))^w1
                d_val = bn12 - state.zb[j, l] / h1
                ah = agwd / h1

                dum_qwx = qwx - qs
                r = abs(dum_qwx) < 1e-3 ? 0.0 : cwd * qs / dum_qwx
                dg_arr[j] = 0.0
                lstart = 0
            end

            cx = d_val + state.zb[jp1, l] / h1
            dgjp1 = dg_arr[j]

            converged_wd = false
            for _iteh in 1:20
                g_arr[jp1] = g_arr[j] + dgjp1
                c = cx + g_arr[jp1]
                if c ≤ 0.0
                    if get(ENV, "CSHORE_WETDRY_DEBUG", "0") == "1"
                        @info "wetdry march terminated (c≤0)" j jp1 c cx zb=state.zb[jp1, l] h1
                    end
                    state.jdry = j
                    _finalize_qotf!(state, config, l, aqwd, jcrest)
                    return nothing
                end
                y = (c / bn12)^(1.0 / w1)
                state.hwd[jp1] = h2 / y
                if state.hwd[jp1] > state.hwd[j]
                    state.hwd[jp1] = state.hwd[j]
                end
                y = h1 / state.hwd[jp1]
                dum_pw = (1.0 + a1) * y^wdn - a_val * y^3
                if dum_pw < 1.0
                    state.pwet[jp1] = state.pwet[j]
                else
                    state.pwet[jp1] = 1.0 / dum_pw
                    if state.pwet[jp1] > state.pwet[j]
                        state.pwet[jp1] = state.pwet[j]
                    end
                end
                qwave = aqwd * state.hwd[jp1] * sqrt(GRAV * state.hwd[jp1] / max(state.pwet[jp1], 1e-12))

                qs_new = qwx - qwave
                if qs_new > 0.0; qs_new = 0.0; end
                state.uswd[jp1] = qs_new / max(state.hwd[jp1], 1e-12)
                dum_qwx2 = qwx - qs_new
                r = abs(dum_qwx2) < 1e-3 ? 0.0 : cwd * qs_new / dum_qwx2
                dg_arr[jp1] = 0.0

                dum_iter = 0.5 * (dg_arr[j] + dg_arr[jp1])
                if abs(dum_iter - dgjp1) ≤ 1e-5
                    g_arr[jp1] = g_arr[j] + dum_iter
                    # Check if we've peaked (about to go downslope)
                    if jp1 < jmax_l && state.zb[jp1 + 1, l] < state.zb[jp1, l]
                        jc_peak = jp1
                        hc_peak = state.hwd[jp1]
                        pc_peak = state.pwet[jp1]
                        qwc = qwx
                        lstart = 2
                    else
                        if j == jwd; lstart = 1; end
                    end
                    converged_wd = true
                    break
                end
                dgjp1 = dum_iter
            end
            if !converged_wd
                # Accept last iteration values
            end
        else
            # ---- DOWNSLOPE MARCH ----
            if lstart == 2
                pci = 1.0 / max(pc_peak, 1e-12)
                qwc2 = qwc^2
                bg = bwd * GRAV
                cpc = 0.5 * pc_peak / bwd / max(hc_peak, 1e-12)
                ab_ds = 0.25 * pc_peak * qwc2 / (bg * max(hc_peak, 1e-12)^3)
                g_arr[j] = 0.0

                qs_ds = state.uswd[j] * state.hwd[j]
                dum_qwc = qwc - qs_ds
                r = abs(dum_qwc) < 1e-3 ? 0.0 : cwd * qs_ds / dum_qwc
                dg_arr[j] = agwd * dx * state.fb2[j, l] * gbwd(r)
                lstart = 0
            end

            dzb = state.zb[jc_peak, l] - state.zb[jp1, l]

            # INFILT=1, IPERM=0: infiltration reduces swash momentum via a
            # water-pressure gradient term WPGH added to dgjp1. Inside the
            # iteration loop, WPGH is recomputed with the latest HWD.
            infilt_active = (config.options.infilt == 1) &&
                            (config.options.iperm == 0) && (wpm > 0.0)
            wpgh_infilt = 0.0
            if infilt_active
                wpgh_infilt = pmg1 * wpm * state.pwet[j] * dx /
                              sqrt(max(state.hwd[j], hwdmin))
            end
            dgjp1 = dg_arr[j] + wpgh_infilt

            converged_ds = false
            for _iteh in 1:20
                g_arr[jp1] = g_arr[j] + dgjp1
                c_ds = cpc * (dzb - g_arr[jp1])
                if c_ds < 0.0; c_ds = 0.0; end
                if hc_peak < 1e-6
                    state.jdry = j
                    _finalize_qotf!(state, config, l, aqwd, jcrest)
                    return nothing
                end
                y = state.hwd[j] / hc_peak

                # Newton iteration for depth
                for _newton in 1:20
                    dum_inv = 1.0 / max(y * y, 1e-12)
                    f_val = y - 1.0 + ab_ds * (dum_inv - 1.0) - c_ds
                    df_val = 1.0 - 2.0 * ab_ds * dum_inv / max(y, 1e-12)
                    if abs(df_val) < 1e-6
                        state.jdry = j
                        _finalize_qotf!(state, config, l, aqwd, jcrest)
                        return nothing
                    end
                    ynew = y - f_val / df_val
                    if abs(ynew - y) ≤ 1e-6
                        y = ynew
                        break
                    end
                    y = ynew
                end

                state.hwd[jp1] = y * hc_peak
                if state.hwd[jp1] < 1e-6
                    state.jdry = j
                    _finalize_qotf!(state, config, l, aqwd, jcrest)
                    return nothing
                end
                if state.hwd[jp1] > state.hwd[j]
                    state.hwd[jp1] = state.hwd[j]
                end

                # Wet probability on downslope — constant for impermeable
                state.pwet[jp1] = pc_peak
                if state.pwet[jp1] > state.pwet[j]
                    state.pwet[jp1] = state.pwet[j]
                end

                qwave = aqwd * state.hwd[jp1] * sqrt(GRAV * state.hwd[jp1] / max(state.pwet[jp1], 1e-12))
                # INFILT=1, IPERM=0: recompute WPGH and accumulate QP with
                # the latest HWD(JP1).
                if infilt_active
                    if jp1 > jcrest
                        state.qp[jp1] = state.qp[j] +
                                        0.5 * dx * wpm * (state.pwet[j] + state.pwet[jp1])
                        h_face = max(0.5 * (state.hwd[j] + state.hwd[jp1]), hwdmin)
                        dum_infilt = wpm * dx * 0.5 * (state.pwet[j] + state.pwet[jp1])
                        wpgh_infilt = pmg1 * dum_infilt / sqrt(h_face)
                        qwx = state.qo[l] - state.qp[jp1]
                        a_val = qwx^2 / bgh3   # update Froude-like parameter
                    else
                        state.qp[jp1] = 0.0
                        wpgh_infilt = 0.0
                    end
                end

                qs_ds = qwx - qwave
                # On landward slope (past crest), qs must be ≥ 0
                if j ≥ jcrest && qs_ds < 0.0; qs_ds = 0.0; end
                if state.hwd[jp1] < 1e-3 && qs_ds > 1e-3; qs_ds = 1e-3; end
                state.uswd[jp1] = qs_ds / max(state.hwd[jp1], 1e-12)

                dum_qds = qwx - qs_ds
                r = abs(dum_qds) < 1e-3 ? 0.0 : cwd * qs_ds / dum_qds
                dg_arr[jp1] = agwd * dx * state.fb2[jp1, l] * gbwd(r)

                # Trapezoidal average of friction gradient + WPGH
                dum_iter = 0.5 * (dg_arr[j] + dg_arr[jp1]) + wpgh_infilt
                if abs(dum_iter - dgjp1) ≤ 1e-5
                    g_arr[jp1] = g_arr[j] + dum_iter
                    # Check if slope reverses (back to upslope before crest)
                    if jp1 < jcrest && jp1 < jmax_l && state.zb[jp1 + 1, l] ≥ state.zb[jp1, l]
                        lstart = 1
                    end
                    converged_ds = true
                    break
                end
                dgjp1 = dum_iter
            end
        end

        # ---- Velocity moments at jp1 ----
        state.umeawd[jp1] = auwd * sqrt(GRAV * state.pwet[jp1] * state.hwd[jp1]) +
                             state.pwet[jp1] * state.uswd[jp1]
        state.sigwd[jp1] = correct * state.hwd[jp1] * sqrt(max(2.0 / max(state.pwet[jp1], 1e-12) -
                            2.0 + state.pwet[jp1], 0.0))
        dum_u = state.umeawd[jp1] - state.uswd[jp1]
        dum1 = state.pwet[jp1] * dum_u^2 - 2.0 * dum_u *
               (state.umeawd[jp1] - state.pwet[jp1] * state.uswd[jp1])
        dum_ustd = agwd * GRAV * state.hwd[jp1] + dum1
        if dum_ustd > 0.0
            state.ustdwd[jp1] = sqrt(dum_ustd)
        else
            state.jdry = j
            _finalize_qotf!(state, config, l, aqwd, jcrest)
            return nothing
        end

        # Oblique (IANGLE=1) longshore velocity in swash
        if config.options.iangle == 1
            state.stheta[jp1] = state.stheta[jwd]
            state.vmeawd[jp1] = auwd * sqrt(GRAV * state.pwet[jp1] * state.hwd[jp1]) * state.stheta[jp1]
            dum_vs = 1.0 - 0.25 * π * state.pwet[jp1] * (2.0 - state.pwet[jp1])
            state.vstdwd[jp1] = awd * sqrt(GRAV * state.hwd[jp1] * max(dum_vs, 0.0)) * abs(state.stheta[jp1])
        end

        # Termination checks — relaxed to 1e-6 m (1 µm) to let the swash
        # march extend further up the bluff face. The OLD reference code
        # used this threshold; NEW had a 1000× more aggressive 1e-3 m
        # cutoff which prematurely terminated the march at the bluff toe,
        # leaving the upper face starved of wave-driven transport and
        # causing the "saturating bluff retreat" failure mode seen in
        # constant-storm stress tests.
        if state.hwd[jp1] < 1e-6
            if get(ENV, "CSHORE_WETDRY_DEBUG", "0") == "1"
                @info "wetdry march terminated (hwd small)" jp1 hwd=state.hwd[jp1]
            end
            state.jdry = jp1
            _finalize_qotf!(state, config, l, aqwd, jcrest)
            return nothing
        end

        # If depth barely changes on downslope, treat as dry
        if (jp1 - jwd > 6) && jp1 > 5
            dh_dx = (state.hwd[jp1] - state.hwd[jp1 - 5]) / (5.0 * dx)
            if dh_dx > -1e-6 && state.hwd[jp1] < 1e4 * hwdmin &&
               state.zb[jp1, l] ≤ state.zb[jp1 - 5, l]
                state.jdry = jp1
                _finalize_qotf!(state, config, l, aqwd, jcrest)
                return nothing
            end
        end

        if j == jend
            state.jdry = jp1
        end
    end

    _finalize_qotf!(state, config, l, aqwd, jcrest)
    return nothing
end

"""
    _finalize_qotf!(state, config, l, aqwd, jcrest)

Compute the combined overtopping flux `state.qotf` at the crest node.
"""
function _finalize_qotf!(state::CshoreState, config::CshoreConfig, l::Int,
                          aqwd::Float64, jcrest::Int)
    state.qotf  = 0.0
    state.sprate = 0.0
    jdam = state.jswl[l]
    if state.jdry ≥ jcrest && jdam < state.jmax[l]
        j = jcrest
        if state.jwd == state.jmax[l]
            j = state.jmax[l]
        end
        if j ≤ length(state.hwd) && state.hwd[j] > 0.0 && state.pwet[j] > 0.0
            state.qotf = aqwd * state.hwd[j] * sqrt(GRAV * state.hwd[j] / state.pwet[j])
        end
    end
    return nothing
end
