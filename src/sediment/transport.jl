# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
transport.jl — Sediment transport.
==============================================================================#

abstract type TransportMethod end

"""
    OriginalCshoreTransport

The FORTRAN SEDTRA formula: `QBX = bld·pb·gslope·ustd³`, with erfc-based
probabilities that depend on grain size only through the initiation
parameter `rb_k = √(gsd50s_k/fb2)/ustd`. This saturates to `pb ≈ 1` in
the surf zone, so it has essentially no grain-size discrimination when
forcing is strong
"""
struct OriginalCshoreTransport <: TransportMethod end

"""
    SvrTransport

Soulsby-Van Rijn (1997) combined bedload + suspended load formula.
Valid for all grain sizes from silt to gravel. Uses depth-averaged
velocity and wave orbital velocity to compute total transport
"""
struct SvrTransport <: TransportMethod end

"""
    SizeAdaptiveTransport(; gravel_threshold=2.0e-3)

Per-fraction adaptive formula: uses the original CSHORE formula for
sand-sized fractions (d < threshold) and Meyer-Peter-Müller (1948)
for gravel/cobble fractions (d ≥ threshold). MPM produces bedload
only — gravel does not suspend.

Default threshold: 2 mm (Wentworth sand-gravel boundary).
"""
struct SizeAdaptiveTransport <: TransportMethod
    gravel_threshold::Float64
end
SizeAdaptiveTransport() = SizeAdaptiveTransport(2.0e-3)

# Select transport method from config symbol
function _select_transport(formula::Symbol)
    if formula == :original
        return OriginalCshoreTransport()
    elseif formula == :size_adaptive
        return SizeAdaptiveTransport()
    elseif formula == :soulsby_vanrijn
        return SvrTransport()
    else
        error("Unknown transport_formula: $formula. Use :original, :size_adaptive, or :soulsby_vanrijn")
    end
end

default_transport(nf::Int) = OriginalCshoreTransport()

@inline function _extrapo_ramp!(F::AbstractMatrix{Float64},
    j_start::Int, j_end::Int,
    n_extrap::Int, nf::Int)
    j_start < j_end || return
    j_ramp_end = min(j_start + n_extrap, j_end)
    n_ramp = j_ramp_end - j_start
    @inbounds for k in 1:nf
        y = F[j_start, k]
        dely = y / (n_extrap + 1)
        for j in 1:n_ramp
            F[j_start+j, k] = y - dely * j
        end
        for j in (j_ramp_end+1):j_end
            F[j, k] = 0.0
        end
    end
    return
end

# ---------------------------------------------------------------------------
# Per-class pickup weighting w_k. The chosen scheme defines the mass
# reservoir that drives surface availability for each grain-size class.
# Selected by config.multifraction.pickup_weighting (see config.jl docstring).
# ---------------------------------------------------------------------------
@inline function _pickup_weight(state::CshoreState, config::CshoreConfig,
                                 j::Int, k::Int)
    mf = config.multifraction
    scheme = mf.pickup_weighting
    if scheme === :active_frac
        return state.active_frac[j, k]
    elseif scheme === :top_layer
        # Top stratigraphic layer mass fraction. Falls back to initial_fractions
        # when the top layer has been fully depleted.
        nf = size(state.bed_mass, 3)
        s = 0.0
        @inbounds for kk in 1:nf
            s += state.bed_mass[j, 1, kk]
        end
        if s ≤ 0.0
            return mf.initial_fractions[k]
        end
        return state.bed_mass[j, 1, k] / s
    elseif scheme === :full_bed
        # Full vertical inventory mass fraction at this node.
        nf       = size(state.bed_mass, 3)
        nlayers  = size(state.bed_mass, 2)
        s_k      = 0.0
        s_total  = 0.0
        @inbounds for il in 1:nlayers, kk in 1:nf
            v = state.bed_mass[j, il, kk]
            s_total += v
            if kk == k
                s_k += v
            end
        end
        if s_total ≤ 0.0
            return mf.initial_fractions[k]
        end
        return s_k / s_total
    else
        # Validated upstream; default fall-through preserves current behavior.
        return state.active_frac[j, k]
    end
end

"""
    sedtra!(state, config, l; method=OriginalCshoreTransport())

Compute per-fraction sediment transport rates along line `l`, restricted to
the normally-incident wet-zone branch.

Populates along line `l`:
- `state.gslope[j]` — smoothed bed-slope correction for bedload
- `state.aslope[j]` — smoothed bed-slope correction for suspended load
- `state.pb[j, k]`  — bedload probability per fraction
- `state.ps[j, k]`  — suspended-load probability per fraction
- `state.vs[j, k]`  — suspended sediment volume per unit area (m)
- `state.qbx[j, k]` — cross-shore bedload transport per fraction (m²/s)
- `state.qsx[j, k]` — cross-shore suspended transport per fraction (m²/s)
- `state.q[j, k]`   — total per-fraction transport `(qbx + qsx) / sporo1`
- `state.q_total[j]` — summed over k and smoothed; feeds `exner_step!`
- `state.pickup_fractions[j, k]` — zeroed here (filled in exner_step!)
"""
function sedtra!(state::CshoreState, config::CshoreConfig, l::Int;
    method::TransportMethod=_select_transport(config.multifraction.transport_formula))
    iangle = config.options.iangle
    iroll = config.options.iroll

    jmax_l = state.jmax[l]
    jr = state.jr
    sed = config.sediment
    tanphi = sed.tanphi
    sporo1 = sed.sporo1
    effb = sed.effb
    efff = sed.efff
    bld = sed.bld
    slp = sed.slp
    bedlm = sed.bedlm
    isedav = config.options.isedav

    mf = config.multifraction
    nf = nfractions(mf)
    nlayers = mf.nlayers

    gslmax = 10.0
    bslop1 = -tanphi * (gslmax - 1.0) / gslmax
    bslop2 = tanphi * (gslmax + 1.0) / (gslmax + 2.0)

    # ---- Slope corrections -----------------------------------------------
    gslraw = Vector{Float64}(undef, jmax_l)
    aslraw = Vector{Float64}(undef, jmax_l)
    @inbounds for j in 1:jmax_l
        bs = state.bslope[j, l]
        if bs < 0.0
            gslraw[j] = bs > bslop1 ? tanphi / (tanphi + bs) : gslmax
        else
            gslraw[j] = bs < bslop2 ? (tanphi - 2 * bs) / (tanphi - bs) : -gslmax
        end
        a = slp
        if bs > 0.0
            a += sqrt(bs / tanphi)
        end
        aslraw[j] = a
    end
    state.gslope[1:jmax_l] .= gslraw
    state.aslope[1:jmax_l] .= aslraw
    smooth_tridiagonal!(view(state.gslope, 1:jmax_l))
    smooth_tridiagonal!(view(state.aslope, 1:jmax_l))

    # Zero per-fraction rate arrays up to jmax.
    fill!(view(state.qbx, 1:jmax_l, :), 0.0)
    fill!(view(state.qsx, 1:jmax_l, :), 0.0)
    fill!(view(state.qby, 1:jmax_l, :), 0.0)
    fill!(view(state.qsy, 1:jmax_l, :), 0.0)
    fill!(view(state.pb, 1:jmax_l, :), 0.0)
    fill!(view(state.ps, 1:jmax_l, :), 0.0)
    fill!(view(state.vs, 1:jmax_l, :), 0.0)
    fill!(view(state.q, 1:jmax_l, :), 0.0)
    fill!(view(state.pickup_fractions, 1:jmax_l, :), 0.0)
    fill!(view(state.q_total, 1:jmax_l), 0.0)
    qraw_total = zeros(Float64, jmax_l)  # summed-over-k total, smoothed at end

    # Guard against NaN/Inf by implementing a decay length
    decayl = min(state.xb[state.jswl[l]] / 4.0, 2.0 * state.wt[1] * state.cp[1])
    jdecay = max(1, round(Int, (isfinite(decayl) ? decayl : 2.0) / config.grid.dx))

    # ---- Per-node, per-fraction loop ------------------------------------
    iv_transport = config.options.iv_transport
    @inbounds for j in 1:jr
        ustd_j = state.ustd[j]
        v_j = iv_transport == 1 ? state.vmean[j] : 0.0
        if ustd_j ≤ 1e-10 && abs(v_j) ≤ 1e-3
            continue
        end

        # Per-class pickup weight w_k = bed-mass fraction in the chosen
        # reservoir (active layer / top stratigraphic layer / full bed).
        # See `_pickup_weight` and MultifractionConfig.pickup_weighting.
        qraw_acc = 0.0
        for k in 1:nf
            comp_k = _pickup_weight(state, config, j, k)
            if comp_k ≤ 0 && sum(view(state.active_frac, j, :)) ≤ 0
                comp_k = mf.initial_fractions[k]
            end
            if iangle == 0 || iv_transport == 1
                qb_k, qs_k, pb_k, ps_k, vs_k = _transport_kernel(
                    method, state, config, j, l, k, comp_k, bld, effb, efff, bedlm, isedav,
                )
                state.qby[j, k] = 0.0
                state.qsy[j, k] = 0.0
            else
                qb_k, qs_k, qby_k, qsy_k, pb_k, ps_k, vs_k = _transport_kernel_oblique(
                    state, config, j, l, k, comp_k, bld, effb, efff, bedlm, isedav,
                )
                state.qby[j, k] = qby_k
                state.qsy[j, k] = qsy_k
            end
            state.pb[j, k] = pb_k
            state.ps[j, k] = ps_k
            state.vs[j, k] = vs_k
            state.qbx[j, k] = qb_k
            state.qsx[j, k] = qs_k
            qraw_k = (qb_k + qs_k) / sporo1
            state.q[j, k] = qraw_k
            qraw_acc += qraw_k
        end
        qraw_total[j] = qraw_acc
    end

    # ─── Phase-lag relaxation on suspended sediment volume v_s ────────────────
    # When PhaseLagConfig is supplied, replace local-equilibrium v_s,eq[j,k]
    # with a non-equilibrium v_s,lag[j,k]:
    #
    #     U(j) ∂v_s/∂x = (v_s,eq − v_s) / τ_lag
    #
    # Then recompute qsx[j,k]
    nl = nonlinearity(config)
    if nl.phase_lag_enabled && jr ≥ 3
        _apply_phase_lag!(state, config, jr, nf, sporo1, qraw_total)
    end

    # ─── Bailard (1981) velocity-moment correction ────────────────────────────
    # When bailard_enabled, add (or replace) qsx with the wave-current cross
    # term `K · u_rms³ · U` and the skewness term `K · u_rms⁴ · Sk` (and
    # optionally asymmetry term `K · u_rms⁴ · As`).
    if nl.bailard_enabled && jr ≥ 3
        _apply_bailard!(state, config, jr, nf, sporo1, qraw_total)
    end

    # Decay padding — flat-extend the seaward edge from jdecay to the
    # boundary so the Exner flux divergence is well-behaved at j=1.
    if jdecay ≥ 2 && jdecay ≤ jr
        @inbounds for k in 1:nf
            qb_d = state.qbx[jdecay, k]
            qs_d = state.qsx[jdecay, k]
            q_d = state.q[jdecay, k]
            for j in 1:(jdecay-1)
                state.qbx[j, k] = qb_d
                state.qsx[j, k] = qs_d
                state.q[j, k] = q_d
            end
        end
        q_d_total = qraw_total[jdecay]
        @inbounds for j in 1:(jdecay-1)
            qraw_total[j] = q_d_total
        end
    end

    # ---- Swash-zone transport -----------------------------------------------
    # When the swash zone extends past the wave runup limit, compute transport
    # from jwd to jdry
    jwd = state.jwd
    jdry = state.jdry
    swash_active = (jdry > jr) && (jwd > 0)

    if swash_active
        dx = config.grid.dx
        agwd_val = config.swash.awd^2
        slpot_val = state.slpot
        qo_l = state.qo[l]
        hwdmin = 0.001

        hdip = Vector{Float64}(undef, jmax_l)
        @inbounds for j in 1:jmax_l
            hdip[j] = max(state.h[j], hwdmin)
        end
        jcr = state.jcrest[l]
        # Forward scan from jwd: track peak bed elevation
        zbpeak = state.zb[jwd, l]
        @inbounds for j in (jwd+1):min(jcr, jdry - 1)
            if state.zb[j-1, l] < state.zb[j, l] && state.zb[j, l] >= state.zb[min(j + 1, jmax_l), l]
                zbpeak = state.zb[j, l]
            end
            dum = zbpeak - state.zb[j, l]
            hdip[j] = max(dum, state.h[j])
        end
        # Backward scan from jdry: track peak from the landward side
        if jdry <= jmax_l && jcr < jdry
            zbpeak = state.zb[min(jdry, jmax_l), l]
            @inbounds for j in (min(jdry, jmax_l)-1):-1:(jcr+1)
                if j >= 2 && state.zb[j-1, l] < state.zb[j, l] &&
                   state.zb[j, l] >= state.zb[min(j + 1, jmax_l), l]
                    zbpeak = state.zb[j, l]
                end
                dum = zbpeak - state.zb[j, l]
                hdip[j] = max(dum, state.h[j])
            end
        end

        # Per-fraction swash transport
        @inbounds for k in 1:nf
            d_k = mf.grain_sizes[k]
            ws_k = state.ws_fractions[k]
            gsd50s_k = state.gsd50s_fractions[k]

            csedia_k = 2.0 * d_k
            # Swash transport arrays (local scratch)
            pbwd = Vector{Float64}(undef, jmax_l)
            pswd = Vector{Float64}(undef, jmax_l)
            vswd = Vector{Float64}(undef, jmax_l)
            qbxwd = Vector{Float64}(undef, jmax_l)
            qsxwd = Vector{Float64}(undef, jmax_l)

            vbf = 1.0
            blds_k = 1.0

            for j in jwd:min(jdry, jmax_l)
                hj = max(state.h[j], hwdmin)
                pw_j = state.pwet[j]
                pw_j = max(pw_j, 1e-10)

                fb2_j = max(state.fb2[j, l], 1e-10)
                ucb = d_k < csedia_k ? sqrt(gsd50s_k / fb2_j) : gsd50s_k
                ucs = ws_k / fb2_j^0.333333

                pwagh = pw_j / max(agwd_val * GRAV * hj, 1e-20)
                pbwd[j] = probwd(pw_j, pwagh, state.uswd[j], ucb)
                pswd[j] = probwd(pw_j, pwagh, state.uswd[j], ucs)
                pswd[j] = min(pswd[j], pbwd[j])

                slope_fac = sqrt(1.0 + state.bslope[j, l]^2)
                vswd[j] = vbf * pswd[j] * slope_fac

                # Swash bedload orbital velocity
                ustd_sw_j = state.ustd[j]
                ustd_swash_j = if config.ig !== nothing && config.ig.ustd_ig_in_transport
                    sqrt(ustd_sw_j * ustd_sw_j + state.ustd_ig[j] * state.ustd_ig[j])
                else
                    ustd_sw_j
                end
                qbxwd[j] = blds_k * pbwd[j] * state.gslope[j] * ustd_swash_j^3

                # Force swash values to match wet-zone values at the jwd transition.
                if j == jwd
                    if vswd[j] > 1e-20
                        vbf = state.vs[j, k] / vswd[j]
                    else
                        vbf = 0.0
                    end
                    vswd[j] = state.vs[j, k]
                    if abs(qbxwd[j]) > 1e-20
                        blds_k = state.qbx[j, k] / qbxwd[j]
                    else
                        blds_k = bld
                    end
                    qbxwd[j] = state.qbx[j, k]
                end

                # Hardbottom supply factor
                if abs(isedav) >= 1
                    dum_hp = state.hp[j, l]
                    brf = dum_hp >= d_k ? 1.0 : (max(dum_hp, 0.0) / d_k)^bedlm
                    vswd[j] *= brf
                    qbxwd[j] *= brf
                end

                qsxwd[j] = state.aslope[j] * state.umean[j] * vswd[j]

                # SLPOT: overtopping-driven onshore suspended transport.
                if slpot_val > 0.0
                    hdip_j = max(hdip[j], hwdmin)
                    qsxwd[j] += slpot_val * vswd[j] * qo_l / hdip_j
                end
            end

            # TRANWD: blend wet-zone (1:jr) and swash (jwd:jdry) transport
            jd = min(jdry, jmax_l)
            if jd > jr
                tranwd!(view(state.pb, :, k), jr, view(pbwd, :), jwd, jd)
                tranwd!(view(state.ps, :, k), jr, view(pswd, :), jwd, jd)
                tranwd!(view(state.vs, :, k), jr, view(vswd, :), jwd, jd)
                tranwd!(view(state.qbx, :, k), jr, view(qbxwd, :), jwd, jd)
                tranwd!(view(state.qsx, :, k), jr, view(qsxwd, :), jwd, jd)
            end

            # Recompute q[j,k] and accumulate qraw_total for blended nodes
            for j in jwd:jd
                qraw_k = (state.qbx[j, k] + state.qsx[j, k]) / sporo1
                state.q[j, k] = qraw_k
                # Only add the increment (wet-zone already accumulated for j≤jr)
                if j > jr
                    qraw_total[j] += qraw_k
                else
                    # jwd..jr nodes were reblended — recalculate total
                    qraw_total[j] += qraw_k - (state.qbx[j, k] + state.qsx[j, k]) / sporo1 + qraw_k
                end
            end
        end

        jd = min(jdry, jmax_l)
        if jd < jmax_l
            hrms_max = maximum(config.boundary.hrmsbc)
            n_extrap = 1 + round(Int, hrms_max / (2.0 * config.grid.dx))
            _extrapo_ramp!(view(state.qbx, :, :), jd, jmax_l, n_extrap, nf)
            _extrapo_ramp!(view(state.qsx, :, :), jd, jmax_l, n_extrap, nf)
            if iangle == 1
                _extrapo_ramp!(view(state.qby, :, :), jd, jmax_l, n_extrap, nf)
                _extrapo_ramp!(view(state.qsy, :, :), jd, jmax_l, n_extrap, nf)
            end
            @inbounds for j in (jd+1):jmax_l
                acc = 0.0
                for k in 1:nf
                    qk = (state.qbx[j, k] + state.qsx[j, k]) / sporo1
                    state.q[j, k] = qk
                    acc += qk
                end
                qraw_total[j] = acc
            end
        end

        # Recompute qraw_total for the swash zone from scratch (clean sum)
        @inbounds for j in jwd:min(jdry, jmax_l)
            acc = 0.0
            for k in 1:nf
                acc += state.q[j, k]
            end
            qraw_total[j] = acc
        end

    else
        if jr < jmax_l
            hrms_max = maximum(config.boundary.hrmsbc)
            n_extrap = 1 + round(Int, hrms_max / (2.0 * config.grid.dx))
            _extrapo_ramp!(view(state.qbx, :, :), jr, jmax_l, n_extrap, nf)
            _extrapo_ramp!(view(state.qsx, :, :), jr, jmax_l, n_extrap, nf)
            if iangle == 1
                _extrapo_ramp!(view(state.qby, :, :), jr, jmax_l, n_extrap, nf)
                _extrapo_ramp!(view(state.qsy, :, :), jr, jmax_l, n_extrap, nf)
            end
            @inbounds for j in (jr+1):jmax_l
                acc = 0.0
                for k in 1:nf
                    qk = (state.qbx[j, k] + state.qsx[j, k]) / sporo1
                    state.q[j, k] = qk
                    acc += qk
                end
                qraw_total[j] = acc
            end
        end
    end

    # Landward-boundary consistency
    if jmax_l ≥ 2
        @inbounds for k in 1:nf
            state.qbx[jmax_l, k] = state.qbx[jmax_l-1, k]
            state.qsx[jmax_l, k] = state.qsx[jmax_l-1, k]
            state.q[jmax_l, k] = state.q[jmax_l-1, k]
            if iangle == 1
                state.qby[jmax_l, k] = state.qby[jmax_l-1, k]
                state.qsy[jmax_l, k] = state.qsy[jmax_l-1, k]
            end
        end
        qraw_total[jmax_l] = qraw_total[jmax_l-1]
    end

    state.q_total[1:jmax_l] .= qraw_total
    smooth_tridiagonal!(view(state.q_total, 1:jmax_l))

    # --- Root biomass transport reduction ---
    # When vegetation with root depth (vegrd > 0) is present, reduce
    # sediment transport by (1 - root_factor) where root_factor ∈ [0,1]
    # scales with root density × depth.  The reduction applies to total transport AND
    # per-fraction fluxes so that the Exner equation sees the reduced
    # rates.
    #
    #   root_factor_j = clamp(vegn_j * vegrd_j / root_ref, 0, max_reduction)
    #
    if config.vegetation !== nothing
        veg = config.vegetation
        if size(veg.vegrd, 1) >= jmax_l
            root_ref = 500.0
            max_reduction = 0.30
            @inbounds for j in 1:jmax_l
                rd = veg.vegrd[j, l]
                nn_veg = veg.vegn[j, l]
                if rd > 0.0 && nn_veg > 0.0
                    rf = clamp(nn_veg * rd / root_ref, 0.0, max_reduction)
                    scale = 1.0 - rf
                    state.q_total[j] *= scale
                    for k in 1:nf
                        state.qbx[j, k] *= scale
                        state.qsx[j, k] *= scale
                    end
                end
            end
        end
    end

    return state
end

#==============================================================================
Per-fraction transport kernels — dispatched on TransportMethod.

Return a 5-tuple `(qbx_k, qsx_k, pb_k, ps_k, vs_k)`. Each kernel assumes the
state has already been populated by `transform_waves!` and the slope
corrections by the outer `sedtra!` preamble.
==============================================================================#

# ---------------------------------------------------------------------------
# Skewness/asymmetry coupling helpers (ISKEW_SPATIAL)
#
# iskew_spatial=0 (default): return the global facSK, facAS scalars verbatim.
# iskew_spatial=1: weight by tanh(Ur / ur_sk_ref) so that coupling is full
#   strength in the nonlinear surf zone (high Ur) and tapers toward zero in
#   the linear shoaling region (low Ur). This prevents a surf-zone-calibrated
#   facSK from over-driving transport far offshore.
# ---------------------------------------------------------------------------
@inline function _skew_coupling(config::CshoreConfig, ur::Float64)
    nl = nonlinearity(config)
    nl.spatial_weighting || return nl.skewness, nl.asymmetry
    w = tanh(ur / nl.ur_reference)
    return nl.skewness * w, nl.asymmetry * w
end

# ---------------------------------------------------------------------------
# Undertow vertical-structure correction (UndertowConfig)
#
# The depth-averaged umean[j] underestimates the near-bed return current in
# the surf zone. _undertow_bed returns an amplified bed-velocity that drives
# suspended-load transport magnitude/direction inside the kernels. Outside
# this function umean is unchanged (drives momentum, IG, etc.).
#
#   :hrms_h        U_bed = U_mean · (1 + α · (Hrms/h)^p)         (capped)
#   :dissipation   U_bed = U_mean - β · D_b / (g · h · c_p)       (capped)
#                  ⤷ canonical CSHORE units: D_b carries m³/s²,
#                    so D_b / (g · h · c_p) has units m/s already.
#
# Returns U_mean unchanged when undertow is disabled or the kernel is in
# linear/dry conditions. The cap limits the amplification to prevent runaway
# in very shallow swash water (where Hrms/h ≈ γ_b but the kernel quickly
# blows up otherwise).
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Phase-lag relaxation on suspended sediment volume (Reniers et al. 2004
# style non-equilibrium closure).
#
# After the per-fraction kernels have populated state.vs[j,k] with the
# *equilibrium* suspended volume v_s,eq, this routine applies an upwind
# spatial relaxation:
#
#     U(j) ∂v_s/∂x = (v_s,eq − v_s) / τ_lag(j,k)
#
# in dimensionless form
#
#     v_s,lag[j] = v_s,lag[j_up] + α(j,k) · (v_s,eq[j] − v_s,lag[j_up])
#     α = dx / max(|U_bed(j)|·τ_lag(j,k),  L_min)
#
# τ_lag is either user-fixed (cfg.tau_lag), or derived per-cell as
# h(j) / w_s(k) when cfg.tau_lag is NaN.  α is capped at cfg.cap_alpha
# (1.0 → no lag, 0.0 → frozen at upstream value).
#
# The march is run twice — once forward (handles U_bed > 0 cells) and once
# backward (handles U_bed < 0 cells, the typical undertow regime).  Cells
# whose local U_bed has the wrong sign for the active sweep are left at
# their already-lagged value.  The boundary condition is v_s = v_s,eq
# at the upwind end.
#
# After the relaxation, qsx[j,k] is recomputed from the lagged volume:
#
#     qsx[j,k] = aslope[j] · U_bed[j] · v_s,lag[j,k]
# ---------------------------------------------------------------------------
function _apply_phase_lag!(state::CshoreState, config::CshoreConfig,
    jr::Int, nf::Int, sporo1::Float64,
    qraw_total::Vector{Float64})
    nl = nonlinearity(config)
    dx = config.grid.dx
    L_min = nl.phase_lag_L_min
    cap_α = nl.phase_lag_cap
    tau_cfg = nl.phase_lag_tau

    # Save equilibrium values (we'll overwrite state.vs with lagged values)
    vs_eq = copy(view(state.vs, 1:jr, 1:nf))

    # Per-cell undertow-corrected velocity
    Ubed = Vector{Float64}(undef, jr)
    @inbounds for j in 1:jr
        Ubed[j] = _undertow_bed(config, state, j, state.umean[j])
    end

    # Per-cell, per-fraction relaxation coefficient α(j,k)
    @inbounds for k in 1:nf
        ws_k = state.ws_fractions[k]
        # Boundary: at offshore (j=1), assume equilibrium (no upstream cloud)
        state.vs[1, k] = vs_eq[1, k]

        # ── Forward sweep (handles U_bed > 0 cells) ─────────────────────────
        for j in 2:jr
            U_j = Ubed[j]
            U_j > 0 || continue   # skip when undertow / no flow this direction
            τ = isnan(tau_cfg) ? max(state.h[j], 0.05) / max(ws_k, 1e-6) : tau_cfg
            L_lag = max(abs(U_j) * τ, L_min)
            α = clamp(dx / L_lag, 0.0, cap_α)
            state.vs[j, k] = state.vs[j-1, k] + α * (vs_eq[j, k] - state.vs[j-1, k])
        end

        # ── Backward sweep (handles U_bed < 0 cells — typical undertow) ─────
        # For these cells, the upstream value is at j+1 (sediment flowing
        # from landward toward seaward).  Run after forward sweep so any
        # forward-only nodes already have their relaxed values to feed in.
        for j in (jr-1):-1:1
            U_j = Ubed[j]
            U_j < 0 || continue
            τ = isnan(tau_cfg) ? max(state.h[j], 0.05) / max(ws_k, 1e-6) : tau_cfg
            L_lag = max(abs(U_j) * τ, L_min)
            α = clamp(dx / L_lag, 0.0, cap_α)
            # Upstream = j+1 (further landward = where the cloud originated)
            state.vs[j, k] = state.vs[j+1, k] + α * (vs_eq[j, k] - state.vs[j+1, k])
        end
    end

    # ── Recompute qsx[j,k] from the lagged volume ──────────────────────────────
    # Bedload qbx[j,k] is left untouched (no phase-lag for bedload).  The
    # kernels store state.vs[j,k] **already weighted by comp_k**, so the
    # recomputation does NOT multiply by comp_k a second time:
    #
    #     qsx[j,k] = aslope · U_bed · v_s,lag(comp-weighted)
    #
    # qraw_total must also be updated because the q_total smoothing pass
    # that follows uses it.
    @inbounds for j in 1:jr
        U_bed_j = Ubed[j]
        as_j = state.aslope[j]
        qraw_acc = 0.0
        for k in 1:nf
            qsx_new = as_j * U_bed_j * state.vs[j, k]
            state.qsx[j, k] = qsx_new
            state.q[j, k] = (state.qbx[j, k] + qsx_new) / sporo1
            qraw_acc += state.q[j, k]
        end
        qraw_total[j] = qraw_acc
    end
    return state
end

# ---------------------------------------------------------------------------
# Bailard (1981) velocity-moment suspended transport correction.
#
# Adds (or replaces) qsx with two physically-motivated terms from the
# decomposition of ⟨|u|³·u⟩ for u(t) = U + u_w(t):
#
#   q_xc(j,k) = aslope · K_s(k) · γ_xc · ustd³ · U_bed       (cross term)
#   q_sk(j,k) = aslope · K_s(k) · γ_sk · ustd⁴ · Sk          (skewness term)
#   q_as(j,k) = aslope · K_s(k) · γ_as · ustd⁴ · As          (asymmetry term)
#
# with K_s(k) = ε_s / ((s−1)·g·w_s,k)  [s³/m²].  Units check:
#   [s³/m²] · [m³/s³] · [m/s] = m²/s  ✓   (cross term)
#   [s³/m²] · [m⁴/s⁴] · [-]   = m²/s  ✓   (skewness/asymmetry)
#
# Sign convention: U_bed < 0 = undertow → q_xc < 0 = offshore;
# Sk > 0 (crest > trough) → q_sk > 0 = onshore.  Cross + skewness terms
# therefore sum to a profile that flips sign across the breakpoint —
# producing the sharp ∂q/∂x required for bar formation.
#
# Bedload Bailard term ⟨|u|²·u⟩ is *not* added here; the existing
# pb·gslope·ustd³ formulation is its instantaneous-equilibrium analog.
# ---------------------------------------------------------------------------
function _apply_bailard!(state::CshoreState, config::CshoreConfig,
    jr::Int, nf::Int, sporo1::Float64,
    qraw_total::Vector{Float64})
    nl = nonlinearity(config)
    sed = config.sediment
    sgm1 = submerged_sgm1(sed)
    # Note: nl.skewness / nl.asymmetry already absorb the old Bailard
    # gamma_sk / gamma_as factors (see fanout in nonlinearity()), so they
    # multiply u_rms⁴·Sk / u_rms⁴·As directly here.

    # Per-cell undertow-corrected velocity (consistent with phase-lag use)
    Ubed = Vector{Float64}(undef, jr)
    @inbounds for j in 1:jr
        Ubed[j] = _undertow_bed(config, state, j, state.umean[j])
    end

    @inbounds for j in 1:jr
        u_rms = state.ustd[j]
        u_rms <= 1e-6 && continue
        u_rms2 = u_rms * u_rms
        u_rms3 = u_rms2 * u_rms
        u_rms4 = u_rms2 * u_rms2
        U = Ubed[j]
        Sk_j = state.skewness[j]
        As_j = state.asymmetry[j]
        aslope = state.aslope[j]

        qraw_acc = 0.0
        for k in 1:nf
            ws_k = state.ws_fractions[k]
            ws_k <= 1e-6 && continue
            K_s = nl.bailard_eps_s / (sgm1 * GRAV * ws_k)

            comp_k = _pickup_weight(state, config, j, k)
            if comp_k ≤ 0 && sum(view(state.active_frac, j, :)) ≤ 0
                comp_k = config.multifraction.initial_fractions[k]
            end

            q_xc = nl.bailard_gamma_xc * u_rms3 * U
            q_sk = nl.skewness * u_rms4 * Sk_j
            q_as = nl.asymmetry * u_rms4 * As_j
            q_bailard = K_s * (q_xc + q_sk + q_as) * aslope * comp_k

            if nl.bailard_additive
                state.qsx[j, k] += q_bailard
            else
                state.qsx[j, k] = q_bailard
            end
            state.q[j, k] = (state.qbx[j, k] + state.qsx[j, k]) / sporo1
            qraw_acc += state.q[j, k]
        end
        qraw_total[j] = qraw_acc
    end
    return state
end

@inline function _undertow_bed(config::CshoreConfig, state::CshoreState,
    j::Int, U_mean::Float64)
    config.undertow === nothing && return U_mean
    und = config.undertow
    und.enabled || return U_mean
    h = max(state.h[j], und.h_min)
    h <= 0.0 && return U_mean

    if und.mode == :hrms_h
        ratio = state.hrms[j] / h
        # Only apply where wave breaking has nontrivial Hrms/h
        ratio <= 0.05 && return U_mean
        amp = 1.0 + und.alpha * ratio^und.exponent
        amp = min(amp, und.cap)
        return U_mean * amp
    elseif und.mode == :dissipation
        Db = state.dbsta[j]
        Db <= 0.0 && return U_mean
        cp_j = state.cp[j]
        cp_j <= 1e-3 && return U_mean
        # ΔU = β · D_b / (g · h · c_p) — always offshore-directed
        delta = und.alpha * Db / (GRAV * h * cp_j)
        # Sign convention: in CSHORE +x is onshore, undertow is negative.
        # We *strengthen* the existing return flow (sign-preserving on U_mean
        # if it's already negative; otherwise add an offshore component).
        U_bed = U_mean - delta
        # Cap relative amplification to prevent unphysical near-bed velocities
        if abs(U_mean) > 1e-6
            amp = abs(U_bed) / abs(U_mean)
            amp > und.cap && (U_bed = sign(U_bed) * und.cap * abs(U_mean))
        end
        return U_bed
    else
        return U_mean
    end
end

"""
    _transport_kernel(::OriginalCshoreTransport, ...)

FORTRAN SEDTRA bedload/suspended with `comp[k]` availability weighting.
Used for single-grain runs and as the FORTRAN-parity reference. In the
surf zone this formula produces little grain-size discrimination because
`pb_k` and `ps_k` saturate near 1 for all grain sizes.
"""
@inline function _transport_kernel(::OriginalCshoreTransport,
    state::CshoreState, config::CshoreConfig,
    j::Int, l::Int, k::Int, comp_k::Float64,
    bld::Float64, effb::Float64, efff::Float64,
    bedlm::Float64, isedav::Int)
    mf = config.multifraction
    sed = config.sediment
    ustd_j = state.ustd[j]
    us = state.usta[j]

    d_k = mf.grain_sizes[k]
    ws_k = state.ws_fractions[k]
    gsd50s_k = state.gsd50s_fractions[k]
    wfsgm1_k = ws_k * submerged_sgm1(sed)

    # IG wave orbital velocity contribution (IgConfig, Layer 1–2 output).
    # Added in quadrature to short-wave ustd: ustd_comb = sqrt(ustd_ss² + ustd_ig²).
    # Only applied when IgConfig is active AND ustd_ig_in_transport = true.
    ustd_ss = if config.ig !== nothing && config.ig.ustd_ig_in_transport
        ustd_ig_j = state.ustd_ig[j]
        sqrt(ustd_j * ustd_j + ustd_ig_j * ustd_ig_j)
    else
        ustd_j
    end

    # When IV_TRANSPORT=1, augment bed shear with the alongshore current vmean.
    # The current increases the combined bed orbital+current velocity used for
    # motion initiation (rb, rs) and bedload magnitude (ustd_eff³), while the
    # dissipation-driven suspended stirring (vs_k) and undertow-directed
    # suspended transport (qsx) remain unchanged.  This naturally increases
    # transport at bar depth where waves alone may be near-threshold, without
    # disrupting the wave-dissipation peak that tracks the outer bar.
    v_j = config.options.iv_transport == 1 ? state.vmean[j] : 0.0
    ustd_eff = sqrt(ustd_ss * ustd_ss + v_j * v_j)
    ustd_eff = max(ustd_eff, 1e-10)

    rb = sqrt(max(gsd50s_k / state.fb2[j, l], 0.0)) / ustd_eff
    rs = ws_k / ustd_eff / state.fb2[j, l]^0.3333

    pb_k = 0.5 * (erfcc((rb + us) / SQR2) + erfcc((rb - us) / SQR2))
    ps_k = 0.5 * (erfcc((rs + us) / SQR2) + erfcc((rs - us) / SQR2))
    ps_k = min(ps_k, pb_k)

    slope_fac = sqrt(1.0 + state.bslope[j, l]^2)
    # When roller is active, roller dissipation (rbeta*rq) replaces wave
    # breaking dissipation (dbsta) as the suspended sediment stirring term.
    db_eff = config.options.iroll == 1 ? state.rbeta[j] * state.rq[j] : state.dbsta[j]
    vs_k = ps_k * (efff * state.dfsta[j] + effb * db_eff) / wfsgm1_k
    vs_k *= slope_fac

    bq = bld * (0.5 + state.qbreak[j])
    qbx_k = bq * pb_k * state.gslope[j] * ustd_eff^3

    # Nonlinear wave asymmetry/skewness (Ruessink et al. 2012).
    # Adds an onshore-biased bedload term proportional to velocity
    # skewness Sk (captures crest-trough asymmetry) and asymmetry As
    # (acceleration skewness, sawtooth wave shape).
    # Sk drives steady onshore transport in the shoaling zone;
    # As captures the "boundary-layer streaming" acceleration effect.
    # Gated by WaveNonlinearityConfig.enabled (or legacy iasym>0 via fanout).
    if nonlinearity(config).enabled
        fsk_j, fas_j = _skew_coupling(config, state.ursell[j])
        qbx_k += bq * pb_k * state.gslope[j] * ustd_eff^3 *
                 (fsk_j * state.skewness[j] + fas_j * state.asymmetry[j])
    end

    if abs(isedav) ≥ 1
        dum = state.hp[j, l]
        brf = dum ≥ d_k ? 1.0 : (dum / d_k)^bedlm
        qbx_k *= brf
        vs_k *= brf
    end

    # Undertow vertical-structure correction (only inside the kernel; depth-
    # averaged umean elsewhere is unchanged). Falls back to umean when the
    # UndertowConfig is absent or disabled.
    U_bed = _undertow_bed(config, state, j, state.umean[j])
    qsx_k = state.aslope[j] * U_bed * vs_k

    # Availability weighting — avoids nf× over-transport for the original-
    # CSHORE formula (sorting is weak anyway, so this is the right trade).
    qbx_k *= comp_k
    qsx_k *= comp_k
    vs_k *= comp_k

    return qbx_k, qsx_k, pb_k, ps_k, vs_k
end

# ============================================================================
# Soulsby-Van Rijn (1997) kernel — valid for all grain sizes
# ============================================================================

"""
    _transport_kernel(::SvrTransport, ...) -> (qbx_k, qsx_k, pb_k, ps_k, vs_k)

Soulsby-Van Rijn (1997) combined bedload + suspended load formula.
Reference: Soulsby, R.L. (1997) "Dynamics of Marine Sands", Thomas Telford.

Key equations:
  Ucr = critical depth-averaged velocity (grain-size dependent)
  Asb = 0.005·h·(d/h)^1.2 / [(s-1)·g·d]^0.6  (bedload coefficient)
  Ass = 0.012·d·D*^(-0.6) / [(s-1)·g·d]^0.6   (suspended coefficient)
  qt  = (Asb + Ass) · [sqrt(U² + 0.018/Cd·Urms²) - Ucr]^2.4
"""
@inline function _transport_kernel(::SvrTransport,
    state::CshoreState, config::CshoreConfig,
    j::Int, l::Int, k::Int, comp_k::Float64,
    bld::Float64, effb::Float64, efff::Float64,
    bedlm::Float64, isedav::Int)
    mf = config.multifraction
    sed = config.sediment
    d_k = mf.grain_sizes[k]
    h_j = state.h[j]
    h_j = max(h_j, 0.01)  # minimum depth

    # Dimensionless grain size D* — uses ρ_s/ρ_w-1 (seawater-aware), not sg-1.
    ν = 1.0e-6
    sgm1 = submerged_sgm1(sed)
    dstar = d_k * (sgm1 * GRAV / ν^2)^(1.0 / 3.0)

    # Critical depth-averaged velocity (Soulsby 1997, eq. 75-76)
    Ucr_d = if d_k ≤ 0.5e-3
        # Fine: threshold from Shields via log-profile
        0.19 * d_k^0.1 * log10(4.0 * h_j / max(d_k, 1e-8))
    else
        # Coarse: modified formula
        8.5 * d_k^0.6 * log10(4.0 * h_j / max(d_k, 1e-8))
    end
    Ucr = max(Ucr_d, 0.01)

    # Drag coefficient from friction factor
    Cd = max(state.fb2[j, l], 1e-6)

    # Effective velocity. With IV_TRANSPORT=1, the alongshore mean current
    # `vmean` contributes to total bed-shear, increasing entrainment when
    # waves are weak. With IV_TRANSPORT=0 (default), behavior is FORTRAN-
    # parity (cross-shore mean only).
    # svr_wave_scale and svr_current_scale allow independent calibration of the
    # wave-orbital and mean-flow contributions to Ueff before the power law,
    # enabling separate tuning of wave-driven vs current-driven transport.
    U_x = _undertow_bed(config, state, j, state.umean[j])
    V_y = config.options.iv_transport == 1 ? state.vmean[j] : 0.0
    Umag = sqrt(U_x * U_x + V_y * V_y)
    # IG orbital velocity contribution (IgConfig). Add IG in quadrature to
    # short-wave Urms so the SVR Ueff integrates both energy components.
    Urms_ss = state.ustd[j]
    Urms = if config.ig !== nothing && config.ig.ustd_ig_in_transport
        sqrt(Urms_ss^2 + state.ustd_ig[j]^2)
    else
        Urms_ss
    end
    wave_sc = config.multifraction.svr_wave_scale
    curr_sc = config.multifraction.svr_current_scale
    Ueff = sqrt(curr_sc * Umag * Umag + wave_sc * 0.018 / Cd * Urms^2)

    if Ueff ≤ Ucr
        return 0.0, 0.0, 0.0, 0.0, 0.0
    end

    # Transport coefficients (Soulsby 1997)
    denom = (sgm1 * GRAV * d_k)^0.6
    Asb = 0.005 * h_j * (d_k / h_j)^1.2 / max(denom, 1e-12)
    Ass = 0.012 * d_k * max(dstar, 0.1)^(-0.6) / max(denom, 1e-12)

    # Total transport magnitude
    excess = (Ueff - Ucr)^2.4
    svr_scale = config.multifraction.svr_scale
    qb_mag = svr_scale * Asb * excess   # bedload component (scaled)
    qs_mag = svr_scale * Ass * excess   # suspended component (scaled)

    # Direction. With IV_TRANSPORT=1, the cross-shore direction smoothly
    # blends from the standard sign(U_x)·gslope formulation (when waves and
    # cross-shore current dominate) toward a slope-cascade direction
    # (sign(downhill)·bs_factor) when the alongshore current dominates.
    # The blend weight w = |V| / (wave_drive + |V|).
    #
    # Direction — split bedload vs suspended, matching the wave-averaged physics
    # of the original CSHORE formula:
    #
    #   Bedload  (qbx): wave-driven → net ONSHORE by default (gslope > 0).
    #     Waves are the primary driver of near-bed bedload; even in the surf
    #     zone the asymmetric wave orbital velocity produces a net onshore
    #     bedload tendency.  This mirrors original CSHORE: bld·pb·gslope·ustd³
    #     is always positive (onshore) because ustd > 0 and gslope ≈ +1.
    #
    #   Suspended load (qsx): undertow-driven → OFFSHORE in the surf zone.
    #     The mean cross-shore current (umean < 0 = undertow) advects suspended
    #     material offshore.  sign_u = sign(umean).
    #
    # With IV_TRANSPORT=1 and a significant V current, a blended downslope-
    # cascade term (weight w) represents the cross-shore consequence of
    # current-driven stirring settling under gravity.
    sign_u = U_x ≥ 0 ? 1.0 : -1.0   # used by suspended load only

    # ── Bedload: wave-driven, onshore by default ──────────────────────────────
    cross_factor_bed = state.gslope[j]   # +gslope → onshore (no sign_u)

    if config.options.iv_transport == 1 && abs(V_y) > 1e-6
        wave_drive = abs(U_x) + Urms * sqrt(0.018 / Cd)
        w = abs(V_y) / (wave_drive + abs(V_y) + 1e-6)
        bs = state.bslope[j]
        if abs(bs) > 1e-6
            sign_down = -sign(bs)                        # downhill in CSHORE +x convention
            bs_factor = min(abs(bs) / sed.tanphi, 1.0)   # 0 on flat, 1 at repose
            cross_factor_bed = (1.0 - w) * state.gslope[j] +
                               w * sign_down * bs_factor
        end
        # negligible slope: pure wave-driven onshore, no change
    end

    qbx_k = qb_mag * cross_factor_bed

    # ── Suspended load: undertow-driven, offshore in surf zone ────────────────
    asuspended = state.aslope[j]
    if config.options.iv_transport == 1 && abs(V_y) > 1e-6
        wave_drive = abs(U_x) + Urms * sqrt(0.018 / Cd)
        w = abs(V_y) / (wave_drive + abs(V_y) + 1e-6)
        bs = state.bslope[j]
        if abs(bs) > 1e-6
            sign_down = -sign(bs)
            bs_factor = min(abs(bs) / sed.tanphi, 1.0)
            qsx_k = (1.0 - w) * sign_u * qs_mag * asuspended +
                    w * sign_down * qs_mag * bs_factor
        else
            qsx_k = (1.0 - w) * sign_u * qs_mag * asuspended
        end
    else
        qsx_k = sign_u * qs_mag * asuspended
    end

    # Nonlinear wave asymmetry (Ruessink 2012): positive Sk biases bedload
    # onshore; applied as an additive term scaled by bedload magnitude.
    # Gated by WaveNonlinearityConfig.enabled (or legacy iasym>0 via fanout).
    if nonlinearity(config).enabled
        fsk_j, fas_j = _skew_coupling(config, state.ursell[j])
        qbx_k += qb_mag * state.gslope[j] *
                 (fsk_j * state.skewness[j] + fas_j * state.asymmetry[j])
    end

    # BRF supply limitation
    if abs(isedav) ≥ 1
        hp_j = state.hp[j, l]
        brf = hp_j ≥ d_k ? 1.0 : (hp_j / d_k)^bedlm
        qbx_k *= brf
        qsx_k *= brf
    end

    # Availability weighting
    qbx_k *= comp_k
    qsx_k *= comp_k

    # Probabilities (diagnostic — not used in formula but stored)
    pb_k = Ueff > Ucr ? 1.0 : 0.0
    ps_k = pb_k
    vs_k = abs(qsx_k) / max(abs(state.umean[j]), 1e-10)

    return qbx_k, qsx_k, pb_k, ps_k, vs_k
end

"""
    _transport_kernel_oblique(state, config, j, l, k, comp_k, bld, effb, efff, bedlm, isedav)

Oblique-incidence sediment transport kernel for IANGLE=1. Returns a 7-tuple:
  `(qbx_k, qsx_k, qby_k, qsy_k, pb_k, ps_k, vs_k)`
"""
@inline function _transport_kernel_oblique(
    state::CshoreState, config::CshoreConfig,
    j::Int, l::Int, k::Int, comp_k::Float64,
    bld::Float64, effb::Float64, efff::Float64,
    bedlm::Float64, isedav::Int,
)
    mf = config.multifraction
    sed = config.sediment

    d_k = mf.grain_sizes[k]
    ws_k = state.ws_fractions[k]
    gsd50s_k = state.gsd50s_fractions[k]
    wfsgm1_k = ws_k * submerged_sgm1(sed)

    ustd_j = state.ustd[j]
    cth = state.ctheta[j]
    sth = state.stheta[j]

    # Oblique velocity projection
    us_j = state.usta[j]
    vs_j = state.vsta[j]
    wsta = us_j * cth + vs_j * sth      # projected along wave dir
    vcus = vs_j * cth - us_j * sth      # cross-wave component
    vcus2 = vcus * vcus

    # Sigt for oblique = ustd / ctheta
    sigt = ustd_j / max(abs(cth), 1e-10)

    rb = sqrt(gsd50s_k / state.fb2[j, l]) / sigt
    rs = ws_k / sigt / state.fb2[j, l]^0.3333

    # Oblique threshold
    fs = rs * rs - vcus2
    if fs < 0.0
        ps_k = 1.0
    else
        fs = sqrt(fs)
        ps_k = 0.5 * (erfcc((fs + wsta) / SQR2) + erfcc((fs - wsta) / SQR2))
    end

    fb = rb * rb - vcus2
    if fb < 0.0
        pb_k = 1.0
    else
        fb = sqrt(fb)
        pb_k = 0.5 * (erfcc((fb + wsta) / SQR2) + erfcc((fb - wsta) / SQR2))
    end
    ps_k = min(ps_k, pb_k)

    # When roller is active, roller dissipation (rbeta*rq) replaces dbsta.
    slope_fac = sqrt(1.0 + state.bslope[j, l]^2)
    db_eff = config.options.iroll == 1 ? state.rbeta[j] * state.rq[j] : state.dbsta[j]
    vs_k = ps_k * (efff * state.dfsta[j] + effb * db_eff) / wfsgm1_k
    vs_k *= slope_fac

    bq = bld * (0.5 + state.qbreak[j])
    vsta2 = vs_j * vs_j
    twos = 2.0 * sth
    dum = bq * pb_k * (ustd_j^2 + state.vstd[j]^2)^1.5

    qbx_k = dum * state.gslope[j] * (1.0 + us_j * vsta2 + twos * vcus)
    qby_k = dum * (vs_j * (1.0 + us_j^2 + vsta2) + twos * wsta)

    # Nonlinear wave asymmetry (Ruessink 2012): onshore-biased bedload
    # from velocity skewness and acceleration asymmetry.
    # Gated by WaveNonlinearityConfig.enabled (or legacy iasym>0 via fanout).
    if nonlinearity(config).enabled
        fsk_j, fas_j = _skew_coupling(config, state.ursell[j])
        qbx_k += dum * state.gslope[j] *
                 (fsk_j * state.skewness[j] + fas_j * state.asymmetry[j])
    end

    # BRF supply limitation
    if abs(isedav) ≥ 1
        dum_hp = state.hp[j, l]
        brf = dum_hp ≥ d_k ? 1.0 : (dum_hp / d_k)^bedlm
        vs_k *= brf
        qbx_k *= brf
        qby_k *= brf
    end

    U_bed = _undertow_bed(config, state, j, state.umean[j])
    qsx_k = state.aslope[j] * U_bed * vs_k
    qsy_k = state.vmean[j] * vs_k

    # Availability weighting
    qbx_k *= comp_k
    qsx_k *= comp_k
    qby_k *= comp_k
    qsy_k *= comp_k
    vs_k *= comp_k

    return qbx_k, qsx_k, qby_k, qsy_k, pb_k, ps_k, vs_k
end

# ============================================================================
# Meyer-Peter-Müller (1948) kernel — gravel bedload only
# ============================================================================

"""
    _mpm_kernel(state, config, j, l, k, comp_k, bedlm, isedav)

Meyer-Peter-Müller (1948) bedload formula for gravel/cobble.
Reference: Meyer-Peter, E. & Müller, R. (1948), IAHR.

    τ* = fb2 · ustd² / [g · (s-1) · d_k]
    qb = 8 · sqrt(g · (s-1) · d³) · max(0, τ* - τ*_cr)^1.5

Returns `(qbx_k, 0.0, pb_k, 0.0, 0.0)` — **no suspended load** (gravel
does not suspend in the water column).
"""
@inline function _mpm_kernel(
    state::CshoreState, config::CshoreConfig,
    j::Int, l::Int, k::Int, comp_k::Float64,
    bedlm::Float64, isedav::Int,
)
    mf = config.multifraction
    sed = config.sediment
    d_k = mf.grain_sizes[k]

    sgm1 = submerged_sgm1(sed)   # ρ_s/ρ_w - 1 (seawater-aware)
    ustd_j = state.ustd[j]
    ustd_j ≤ 1e-10 && return (0.0, 0.0, 0.0, 0.0, 0.0)

    # Shields parameter: τ* = fb2 · u*² / [g(s-1)d]
    # where u* ≈ sqrt(fb2) · ustd (bottom shear velocity)
    fb2_j = state.fb2[j, l]
    tau_star = fb2_j * ustd_j^2 / (GRAV * sgm1 * d_k)

    # Critical Shields (Soulsby 1997 curve)
    tau_cr = state.theta_cr_fractions[k]

    excess = tau_star - tau_cr
    if excess ≤ 0.0
        return 0.0, 0.0, 0.0, 0.0, 0.0
    end

    # MPM formula: qb* = 8 · (τ* - τ*_cr)^1.5
    # Dimensional: qb = qb* · sqrt(g(s-1)d³)
    qb_dim = 8.0 * sqrt(GRAV * sgm1 * d_k^3) * excess^1.5

    # Direction: follows slope correction (downslope bedload)
    qbx_k = qb_dim * state.gslope[j]

    # BRF supply limitation
    if abs(isedav) ≥ 1
        hp_j = state.hp[j, l]
        brf = hp_j ≥ d_k ? 1.0 : (hp_j / d_k)^bedlm
        qbx_k *= brf
    end

    # Availability weighting
    qbx_k *= comp_k

    # Probability of motion (diagnostic)
    pb_k = excess > 0 ? min(1.0, excess / tau_cr) : 0.0

    # No suspended load for gravel
    return qbx_k, 0.0, pb_k, 0.0, 0.0
end

# ============================================================================
# SizeAdaptiveTransport — per-fraction dispatch by grain size
# ============================================================================

"""
    _transport_kernel(::SizeAdaptiveTransport, ...) -> 5-tuple

Dispatches per-fraction based on grain diameter:
- d < threshold (default 2 mm): Original CSHORE formula (sand) by default.
  When `iv_transport=1`, sand falls through to the Soulsby-Van Rijn kernel
  instead, because the original CSHORE bedload `~ ustd³` is purely wave-
  driven and zeroes out without waves; SVR carries a current-aware
  entrainment term that reads vmean.
- d ≥ threshold: Meyer-Peter-Müller (gravel bedload only).
"""
@inline function _transport_kernel(m::SizeAdaptiveTransport,
    state::CshoreState, config::CshoreConfig,
    j::Int, l::Int, k::Int, comp_k::Float64,
    bld::Float64, effb::Float64, efff::Float64,
    bedlm::Float64, isedav::Int)
    d_k = config.multifraction.grain_sizes[k]
    if d_k < m.gravel_threshold
        if config.options.iv_transport == 1
            # Sand + current-aware: Soulsby-Van Rijn (which honors iv_transport)
            return _transport_kernel(SvrTransport(),
                state, config, j, l, k, comp_k, bld, effb, efff, bedlm, isedav)
        else
            # Sand, FORTRAN parity: original CSHORE bedload
            return _transport_kernel(OriginalCshoreTransport(),
                state, config, j, l, k, comp_k, bld, effb, efff, bedlm, isedav)
        end
    else
        # Gravel: use MPM (bedload only)
        return _mpm_kernel(state, config, j, l, k, comp_k, bedlm, isedav)
    end
end
