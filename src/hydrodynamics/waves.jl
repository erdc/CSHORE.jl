# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
waves.jl — Wave transformation routines.
==============================================================================#

"""
    WaveParams

Offshore-boundary wave parameters interpolated from the BC time series at the
current model time.
"""
Base.@kwdef struct WaveParams
    tp::Float64        # period (s)
    hrms0::Float64     # offshore Hrms (m)
    angle::Float64     # offshore wave angle (rad from shorenormal)
    swl::Float64       # still water level (m)
    wkpo::Float64      # deep-water wavenumber ω²/g
end

"""
    lwave_dispersion(wd, wkpo, tp; qdisp=0.0, x0=nothing, tol=1e-7) -> (wkp, x, wn, cp, wt)

Solves the linear wave dispersion relation ω² = g·k·tanh(k·h) via Newton
iteration in the non-dimensional variable `x = k·h`, returning

    wkp  = k         wavenumber (1/m)
    x    = k·h
    wn   = Cg/Cp     group-velocity ratio
    cp   = ω/k       phase velocity (m/s)
    wt   = 2π/ω      local period (s; equal to `tp` with no current)

`qdisp` is a current-flux correction (zero for the common IWCINT=0 case).

`x0` is an optional initial guess — pass the previous node's `k·h` during the
landward march for faster convergence.
"""
function lwave_dispersion(wd::Float64, wkpo::Float64, tp::Float64;
                          qdisp::Float64=0.0, x0::Union{Nothing,Float64}=nothing,
                          tol::Float64=1e-7, max_iter::Int=200)
    wd > 0 || throw(ArgumentError("water depth must be > 0, got $wd"))
    D = wd * wkpo
    x = x0 === nothing ? D / sqrt(tanh(D)) : x0

    if qdisp == 0.0
        # No current — direct Newton on F(x) = x - D·coth(x)
        for _ in 1:max_iter
            coth = 1.0 / tanh(x)
            xnew = x - (x - D * coth) / (1.0 + D * (coth^2 - 1.0))
            if abs(xnew - x) ≤ tol
                x = xnew
                break
            end
            x = xnew
        end
        af = PI2 / tp
    else
        # With current — Newton on F(x) = x - D·C²·coth(x), C = 1 - B·x
        B = tp * qdisp / PI2 / wd / wd
        for _ in 1:max_iter
            coth = 1.0 / tanh(x)
            C = 1.0 - B * x
            F  = x - D * C * C * coth
            FD = 1.0 + D * C * (2.0 * B * coth + C * (coth*coth - 1.0))
            xnew = x - F / FD
            if abs(xnew - x) ≤ tol
                x = xnew
                break
            end
            x = xnew
        end
        af = sqrt(GRAV * x * tanh(x) / wd)
    end

    wkp = x / wd
    x2  = 2x
    wn  = 0.5 * (1.0 + x2 / sinh(x2))
    wt  = PI2 / af
    cp  = af / wkp
    return (wkp=wkp, x=x, wn=wn, cp=cp, wt=wt, af=af)
end

"""
    lwave!(state, config, j, l, wd, wave, qdisp=0.0)

Populates `state.cp[j]`, `state.wn[j]`, `state.stheta[j]`, `state.ctheta[j]`,
`state.fsx/fsy/fe` for node `j` along line `l`. Updates `state.wkp` (scalar)
and — on j=1 — also sets `state.wkpsin` (the Snell invariant `k·sin θ`).

Arguments:
- `wd` — mean water depth at the node (m).
- `wave` — current offshore `WaveParams` (holds `tp`, `angle`, `wkpo`).
- `qdisp` — water flux affecting wave period (0 unless IWCINT=1).
"""
function lwave!(state::CshoreState, config::CshoreConfig,
                j::Int, l::Int, wd::Float64, wave::WaveParams;
                qdisp::Float64=0.0, hrms_j::Float64=0.0)
    x0 = j == 1 ? nothing : state.wkp * wd
    disp = lwave_dispersion(wd, wave.wkpo, wave.tp; qdisp=qdisp, x0=x0)

    wkp_final = disp.wkp
    cp_final  = disp.cp
    wn_final  = disp.wn
    wt_final  = disp.wt
    af_final  = disp.af

    # ── Feature 1: Kirby & Dalrymple (1986) nonlinear dispersion correction ──
    # ω² = g·k·tanh(kh)·[1 + (k·a)²·f(kh)]
    # f(kh) = (8 + cosh(4kh) - 2·tanh²(kh)) / (8·sinh⁴(kh))
    # a = Hrms/√2 is the RMS amplitude of the representative wave.
    # Newton iteration from the linear solution. K&D always yields k_nl ≤ k_lin
    # (longer wavelength); if the correction is too large to converge we fall
    # back to the linear solution gracefully.
    if config.options.inl_dispersion == 1 && hrms_j > 1e-6
        ω   = PI2 / wave.tp
        ω2  = ω * ω
        a   = hrms_j / sqrt(2.0)      # RMS amplitude
        k   = disp.wkp                # start from linear solution
        for _ in 1:15                 # up to 15 Newton steps
            kh_loc = k * wd
            sh_kh  = sinh(kh_loc)
            sh_kh < 1e-10 && break    # extremely shallow water — keep current k
            # K&D f(kh), Kirby & Dalrymple (1986) Eq. 5
            fkh = (8.0 + cosh(4.0 * kh_loc) - 2.0 * tanh(kh_loc)^2) /
                  (8.0 * sh_kh^4)
            !isfinite(fkh) && (fkh = 0.0)   # deep water: cosh overflows but correction → 0
            ka2 = (k * a)^2
            F   = GRAV * k * tanh(kh_loc) * (1.0 + ka2 * fkh) - ω2
            # Numerical dF/dk
            dk   = k * 1e-5
            kh2  = (k + dk) * wd
            sh2  = sinh(kh2)
            fkh2 = (sh2 < 1e-10 || !isfinite(sh2)) ? fkh :
                   (8.0 + cosh(4.0*kh2) - 2.0*tanh(kh2)^2) / (8.0*sh2^4)
            !isfinite(fkh2) && (fkh2 = fkh)
            ka2b = ((k + dk) * a)^2
            Fdk  = GRAV * (k + dk) * tanh(kh2) * (1.0 + ka2b * fkh2) - ω2
            dFdk = (Fdk - F) / dk
            abs(dFdk) < 1e-20 && break
            k_new = k - F / dFdk
            # Physics: K&D always gives k_nl ≤ k_lin. Clamp strictly.
            # Lower bound: 0.3·k_lin (K&D not valid beyond ~3× wave-length change).
            k_new = clamp(k_new, disp.wkp * 0.3, disp.wkp)
            abs(k_new - k) < disp.wkp * 1e-7 && (k = k_new; break)
            k = k_new
        end
        # Recompute derived quantities from corrected k.
        # Period is conserved (ω fixed by offshore forcing); cp = ω/k_nl.
        kh_c     = k * wd
        x2_c     = 2.0 * kh_c
        wn_final  = 0.5 * (1.0 + x2_c / sinh(x2_c))
        cp_final  = ω / k               # correct K&D phase speed: ω/k_nl
        wt_final  = disp.wt             # period unchanged — set by offshore forcing
        wkp_final = k
    end

    state.wkp      = wkp_final
    state.wkp_arr[j] = wkp_final    # store local k for ursell/biphase
    state.cp[j]    = cp_final
    state.wn[j]    = wn_final
    state.wt[j]    = wt_final

    # Radiation-stress & energy-flux coefficients
    fsx = 2.0 * disp.wn - 0.5
    fsy = 0.0
    fe  = disp.wn * disp.cp * disp.wt

    if config.options.iangle == 0
        state.stheta[j] = 0.0
        state.ctheta[j] = 1.0
    else
        # Oblique incidence — Snell's law
        if j == 1
            line_angle = length(config.bathymetry.agline) ≥ l ? config.bathymetry.agline[l] : 0.0
            dum = wave.angle - line_angle
            dum = dum > 180.0 ? dum - 360.0 : dum
            dum = dum < -180.0 ? dum + 360.0 : dum
            dum = clamp(dum, -80.0, 80.0)
            θ = deg2rad(dum)
            state.stheta[1] = sin(θ)
            state.ctheta[1] = cos(θ)
            state.wkpsin = state.wkp * state.stheta[1]
        else
            # Clamp to [-1,1]: WCI can reduce k enough that k₀·sinθ₀/k_new > 1
            # (total internal reflection). Cap at grazing incidence rather than crash.
            state.stheta[j] = clamp(state.wkpsin / state.wkp, -1.0, 1.0)
            θ = asin(state.stheta[j])
            state.ctheta[j] = cos(θ)
        end
        fsx -= disp.wn * state.stheta[j]^2
        fsy  = disp.wn * state.stheta[j] * state.ctheta[j]
        fe  *= state.ctheta[j]
    end

    if config.options.iwcint == 1
        fe += disp.wt * state.qwx / wd
    end

    state.fsx = fsx
    state.fsy = fsy
    state.fe  = fe

    # Roller kinematics
    if config.options.iroll == 1
        if config.options.iangle == 0
            state.rx[j] = disp.cp / GRAV
            state.re[j] = state.rx[j] * disp.cp
        else
            dum = disp.cp * state.ctheta[j] / GRAV
            state.rx[j] = dum * state.ctheta[j]
            state.ry[j] = dum * state.stheta[j]
            state.re[j] = dum * disp.cp
        end
        state.rbeta[j] = state.rbzero
        if state.bslope[j, l] > 0.0
            state.rbeta[j] += state.bslope[j, l] * state.ctheta[j]
        end
    end

    return state
end

"""
    dbreak!(state, config, j, l, whrms, wd, wt_j)

Computes at node `j`:
- `state.abreak[j]` — steep-slope dissipation factor
- `state.qbreak[j]` — fraction of breaking waves (Newton iteration)
- `state.dbsta[j]`  — time-averaged normalized breaking dissipation

Arguments:
- `whrms` — local Hrms used for the breaking criterion
- `wd`    — local mean water depth
- `wt_j`  — local wave period at node `j` (equal to `tp` in the common IWCINT=0 case)
"""
function dbreak!(state::CshoreState, config::CshoreConfig,
                 j::Int, l::Int, whrms::Float64, wd::Float64, wt_j::Float64)
    ab = (PI2 / (state.wkp * wd)) * state.bslope[j, l] * state.ctheta[j] / 3.0
    ab = max(ab, 1.0)
    state.abreak[j] = ab

    # Max breaker height. The breaker index γ can be:
    #   :constant       — fixed `gamma` (FORTRAN parity)
    #   :ruessink2003   — γ(kh) = clamp(a·kh + b, γ_min, γ_max)
    #   :steepness_sr   — γ(Hs/L) = clamp(c·√(Hs/L), γ_min, γ_max)
    #                     where Hs = √2·Hrms and L = 2π/k.
    #                     Coefficient `c` (config.gamma_sr_slope) defaults to
    #                     3.9 — SR-trained on FRF Duck 2008–2026 breaking-only
    #                     observations (~29k tuples). Cuts surf-zone Hs RMSE
    #                     ~60 % vs constant γ=0.78 at FRF.
    kh = state.wkp * wd
    γ = if config.gamma_method === :ruessink2003
        clamp(config.gamma_a * kh + config.gamma_b, config.gamma_min, config.gamma_max)
    elseif config.gamma_method === :steepness_sr
        L_local = state.wkp > 0 ? (2π / state.wkp) : 0.0
        Hs_over_L = (L_local > 0 && whrms > 0) ?
                     (sqrt(2.0) * whrms / L_local) : 0.0
        clamp(config.gamma_sr_slope * sqrt(max(Hs_over_L, 0.0)),
              config.gamma_min, config.gamma_max)
    else
        config.gamma
    end
    hm = 0.88 / state.wkp * tanh(γ * kh / 0.88)

    # ── Feature 2: crest-height correction via local skewness ─────────────────
    # A skewed wave has a higher crest than Hrms/2 implies; it breaks at a
    # lower Hrms. Reduce hm proportionally to local Sk (Ruessink 2012):
    #   hm_eff = hm / (1 + alpha_sk · max(0, Sk))
    # Sk computed directly from the local Ursell (no extra state needed).
    nl_w = nonlinearity(config)
    if nl_w.alpha_crest > 0.0 && nl_w.enabled && whrms > 1e-6
        ur_loc = ursell_number(whrms, state.wkp, wd)
        sk_loc, _ = ruessink_skewness_asymmetry(ur_loc)
        hm = hm / (1.0 + nl_w.alpha_crest * max(0.0, sk_loc))
    end

    # Fraction of breaking waves: Newton iteration on
    #   Qb - (1 - Qb + B·ln Qb)/(B/Qb - 1) = 0, with B = (Hrms/Hm)²
    if whrms > 1e-10
        B = (whrms / hm)^2
        if B < 0.99999
            qbold = B / 2.0
            qb = qbold
            for _ in 1:200
                qb = qbold - (1.0 - qbold + B * log(qbold)) / (B / qbold - 1.0)
                if qb ≤ 0.0
                    qb = qbold / 2.0
                end
                if abs(qb - qbold) ≤ 1e-6
                    break
                end
                qbold = qb
            end
            state.qbreak[j] = qb
        else
            state.qbreak[j] = 1.0
            hm = whrms
        end
    else
        state.qbreak[j] = 0.0
        hm = whrms
    end

    # Apply breaker delay: smooth the breaking threshold Hm transition.
    # This spreads the onset of breaking over a transition zone rather than
    # a sharp threshold, reducing gradient discontinuities in dissipation.
    # B = (Hrms/Hm)² ranges from 0 (not breaking) to 1+ (fully broken).
    # Delayed threshold: H_delayed = H_base * (1 + delay*(1 - exp(-B)))
    # This expands Hm when B is small, delaying onset, then allows normal breaking.
    qbreak_delayed = state.qbreak[j]
    if config.breaker_delay > 0.0 && whrms > 1e-10
        B = (whrms / hm)^2
        # Smooth transition factor: 0 at B=0, approaches 1 as B increases
        smooth_factor = 1.0 - exp(-B / (config.breaker_delay + 0.1))
        # Recalculate qbreak with delayed/smoothed threshold
        # Re-solve with effective Hm expanded by delay factor
        hm_delayed = hm / (1.0 + config.breaker_delay * (1.0 - smooth_factor))
        B_delayed = (whrms / hm_delayed)^2
        if B_delayed < 0.99999
            qbold = B_delayed / 2.0
            qb = qbold
            for _ in 1:50
                qb = qbold - (1.0 - qbold + B_delayed * log(qbold)) / (B_delayed / qbold - 1.0)
                if qb ≤ 0.0
                    qb = qbold / 2.0
                end
                if abs(qb - qbold) ≤ 1e-6
                    break
                end
                qbold = qb
            end
            qbreak_delayed = qb
        else
            qbreak_delayed = 1.0
        end
    end

    # Normalized dissipation rate
    state.dbsta[j] = 0.25 * ab * qbreak_delayed * hm * hm / wt_j

    return state
end

"""
    friction_coefficients(ctheta, usigt, stheta, vsigt, iangle) -> (gbx, gf)

Pure function — computes bottom friction coefficients `gbx` and `gf` from
wave angle and normalized velocity amplitudes.

For `iangle = 0` (normal incidence) uses the closed-form expression involving
`erfc(usigt/√2)`. For `iangle = 1` (oblique) uses the approximate formulas.
"""
function friction_coefficients(ctheta::Float64, usigt::Float64,
                               stheta::Float64, vsigt::Float64, iangle::Int)
    if iangle == 1
        rm  = -usigt * ctheta - vsigt * stheta
        afm = abs(vsigt * ctheta - usigt * stheta)
        dum = usigt * usigt + vsigt * vsigt
        gbx = SQRG1_PI * (usigt - rm * ctheta) + usigt * afm
        gf  = SQRG2_PI + (1.0 + dum) * afm + SQRG1_PI * (dum + 2.0 * rm * rm)
        return (gbx, gf)
    else
        # Note: SQRG1_PI = √(2/π) — a Gaussian-friction constant, not √g.
        c1 = 1.0 - erfcc(usigt / SQR2)                     # = erf(usigt/√2)
        c2 = SQRG1_PI * exp(-usigt * usigt / 2.0)
        c3 = 1.0 + usigt * usigt
        gbx = c3 * c1 + c2 * usigt
        gf  = usigt * (c3 + 2.0) * c1 + (c3 + 1.0) * c2
        return (gbx, gf)
    end
end

# Gaussian-friction constants. SQRG1_PI = √(2/π), SQRG2_PI = 2√(2/π).
# Distinct from the physics-constant `SQRG1 = √g` in config.jl.
const SQRG1_PI = sqrt(2.0 / π)
const SQRG2_PI = 2.0 * SQRG1_PI

"""
    longshore_vstgby(ctheta, usigt, stheta, gby) -> vsigt

Pure function. Computes `vsigt = VMEAN/SIGT` from the longshore shear stress
factor `gby`, given the wave angle and cross-shore `usigt`.

For `gby == 0` returns 0 directly (the common case when `iangle == 0` and
hence no longshore current is driven).
"""
function longshore_vstgby(ctheta::Float64, usigt::Float64, stheta::Float64, gby::Float64)
    gby == 0.0 && return 0.0
    B = SQRG1_PI * (1.0 + stheta * stheta)
    C = gby
    if gby > 0.0
        D = B * B + 4.0 * ctheta * C
        vsigt = D ≥ 0.0 ? 0.5 * (sqrt(D) - B) / ctheta : 0.0
        return vsigt < 0.0 ? 0.0 : vsigt
    else
        D = B * B - 4.0 * ctheta * C
        vsigt = D ≥ 0.0 ? 0.5 * (B - sqrt(D)) / ctheta : 0.0
        return vsigt > 0.0 ? 0.0 : vsigt
    end
end

"""
    gby_from_vsigt(ctheta, stheta, vsigt) -> gby

Analytical inverse of `longshore_vstgby`. Given a desired normalized longshore
velocity `vsigt = vmean/sigt` at a node, returns the longshore shear-stress
coefficient `gby` that the forward solver would need to produce that velocity.

Derivation: from `longshore_vstgby`, with B = √(2/π)·(1+sin²θ),

  gby > 0, vsigt > 0:
    vsigt = 0.5·(√(B² + 4·c_θ·gby) − B)/c_θ
    ⇒ gby = c_θ·vsigt² + B·vsigt

  gby < 0, vsigt < 0:
    vsigt = 0.5·(B − √(B² − 4·c_θ·gby))/c_θ
    ⇒ gby = −c_θ·vsigt² + B·vsigt

Combined:  gby = c_θ·vsigt·|vsigt| + B·vsigt  =  vsigt · (c_θ·|vsigt| + B)

This is exact (no iteration). Used by ICURRENT=1 to back-solve the alongshore
water-surface gradient DETADY required to produce a user-prescribed
alongshore current at the offshore boundary.
"""
function gby_from_vsigt(ctheta::Float64, stheta::Float64, vsigt::Float64)
    vsigt == 0.0 && return 0.0
    B = SQRG1_PI * (1.0 + stheta * stheta)
    return vsigt * (ctheta * abs(vsigt) + B)
end

#==============================================================================
Wave nonlinearity — Ruessink et al. (2012) parameterization

Ruessink, B.G., G. Ramaekers, and L.C. van Rijn (2012), "On the parameteri-
zation of the free-stream non-linear wave orbital motion in nearshore
morphodynamic models", Coastal Engineering, 65, 56-63.

The parameterization estimates the velocity skewness (Sk) and asymmetry (As)
of shoaling nonlinear waves from the Ursell number alone. Widely used in
XBeach, Delft3D, and other process-based coastal models.
==============================================================================#

"""
    ursell_number(hrms, k, h) -> Float64

Ursell number Ur = (3/8) · Hs · k / (k·h)³  with Hs = √2·Hrms (Hm0).
Ur → 0 in deep water (linear waves); Ur → O(100) in the surf zone
(highly nonlinear, skewed/asymmetric waves).

Note: `k` is the local wavenumber (1/m) from `state.wkp_arr`, NOT Cg/Cp.
"""
@inline function ursell_number(hrms::Float64, k::Float64, h::Float64)
    (h <= 0.0 || k <= 0.0) && return 0.0
    kh = k * h
    hs = sqrt(2.0) * hrms   # Hm0 ≈ √2 Hrms
    return 0.375 * hs * k / kh^3
end

"""
    ruessink_skewness_asymmetry(ursell) -> (Sk, As)

Compute velocity skewness `Sk` (crest amplification) and asymmetry `As`
(sawtooth shape) from the Ursell number via the Ruessink 2012 Boltzmann
sigmoid. Returns `(0.0, 0.0)` for `ursell <= 0`.

Parameters (Ruessink et al. 2012, Table 2):
  p1=0.0, p2=0.857, p3=-0.471, p4=0.297, p5=0.815, p6=0.672

Sk is positive for shoaling waves (crest > trough); As is negative (front-
leaning sawtooth) before breaking. Both approach ~0.6 in the surf zone.
"""
@inline function ruessink_skewness_asymmetry(ursell::Float64)
    ursell <= 0.0 && return (0.0, 0.0)
    p1, p2, p3, p4, p5, p6 = 0.0, 0.857, -0.471, 0.297, 0.815, 0.672
    # Amplitude envelope B(Ur) — approaches p2 ≈ 0.86 for large Ur
    log_ur = log(ursell)
    B = p1 + (p2 - p1) / (1.0 + exp((p3 - log_ur) / p4))
    # Phase ψ(Ur) — rotates from π/2 (pure asymmetry) to 0 (pure skewness)
    psi = -π / 2.0 * (1.0 - tanh(p5 / ursell^p6))
    Sk = B * cos(psi)
    As = B * sin(psi)
    return (Sk, As)
end

"""
    stokes2_skewness(wn, h) -> Sk

Zero-calibration 2nd-order Stokes skewness for a regular wave of
wavenumber `wn` in depth `h`.  Derived from the 2nd-order velocity
amplitude ratio u₂/u₁ = 3/(4·sinh³(kh)) (Dean & Dalrymple 1991, Eq. 4.23
for infinite depth; full form below).

For a 2nd-order Stokes wave u(t) = u₁cos(ωt) + u₂cos(2ωt), the normalised
velocity skewness Sk = <u³> / σ_u³ has the closed-form expression

    Sk = 3 · (u₂/u₁) / (1 + (u₂/u₁)²)^(3/2)

which saturates at ≈0.7 as kh → 0 and → 0 as kh → ∞. No fitted coefficients.

2nd-order Stokes produces symmetric crest-trough skewness but no sawtooth
asymmetry (that requires 3rd-order + boundary-layer streaming); so this
function returns `As = 0`.

Validity: kh ≳ 0.3 (otherwise the 2nd-order expansion breaks down and
Ursell becomes large). For highly nonlinear waves in very shallow water,
Ruessink is a better option.

Returns `(Sk, As) = (0, 0)` for unphysical inputs.
"""
@inline function stokes2_skewness(k::Float64, h::Float64, hrms::Float64=0.0)
    (k <= 0.0 || h <= 0.0) && return (0.0, 0.0)
    kh = k * h
    sinh_kh = sinh(kh)
    sinh_kh < 1e-6 && return (0.0, 0.0)

    # Bed velocity amplitude ratio u₂/u₁ for 2nd-order Stokes
    # (Dean & Dalrymple 1991, eq. 4.43; see also Isobe & Horikawa 1982).
    # u₁ = a·ω / sinh(kh)
    # u₂ = (3/4) · a² · ω · k / sinh⁴(kh)
    # Ratio: u₂/u₁ = (3/4) · a · k / sinh³(kh) = (3/8) · H · k / sinh³(kh)
    # Use Hrms if provided (phase-averaged), else assume a=H/2 with H=Hm0≈√2·Hrms
    # which gives the "representative" wave amplitude in a random sea.
    a = hrms > 0.0 ? sqrt(2.0) * hrms / 2.0 : 0.0
    a > 0.0 || return (0.0, 0.0)
    r = 0.75 * a * k / sinh_kh^3

    # Cap in deep water: as kh → ∞ the ratio → 0 naturally; in very shallow
    # water the expansion breaks down (Stokes invalid for Ur >> 1).
    r = min(r, 0.5)

    # Skewness of u(t) = u₁cos(ωt) + u₂cos(2ωt):
    #   <u²>  = (u₁² + u₂²)/2
    #   <u³>  = (3/4) u₁² u₂
    #   Sk    = <u³> / <u²>^(3/2)
    #         = (3/4) u₂ / ((1 + r²)/2)^(3/2) · (u₁² scales out after norm)
    #         = (3/4) r / ((1 + r²)/2)^(3/2)
    # Let u₁ = 1; then u₂ = r.
    mean_u2 = 0.5 * (1.0 + r * r)
    mean_u3 = 0.75 * r
    Sk = mean_u3 / mean_u2^1.5
    return (Sk, 0.0)
end

# ---------------------------------------------------------------------------
# Spatially-varying friction (IFRICTION_SPATIAL=1)
#
# Updates state.fb2[j, l] based on the local Shields-parameter regime:
#
#   θ = f_min · ustd² / [g · (s−1) · d50]   ← computed from grain-roughness
#                                               friction (f_min) to avoid the
#                                               circular dependency on fb2.
#
# Three regimes:
#   θ < θ_cr           → grain roughness only  (fb2 = f_min)
#   θ_cr ≤ θ < θ_sheet → ripple regime         (fb2 = f_min·(θ/θ_cr)^exp)
#   θ ≥ θ_sheet        → sheet flow            (fb2 = f_sheet, capped)
#
# Called from the driver predictor step (before the corrector loop) so that
# the updated friction is used consistently in both tbxsta and dfsta.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Dynamic Manning's n friction (IFRICTION_SPATIAL=2)
#
# Recomputes the bottom friction factor each step from a depth-averaged
# Manning's n and the live water depth:
#
#     fb2 = g · n² / h^(1/3)
#
# Quadratic-drag interpretation: CSHORE's fb2 enters the bed stress as
# τ_bed = fb2 · ρ · U² (i.e. fb2 ≡ Cf), which is consistent with the
# Manning relation Cf = g·n²/h^(1/3).
#
# The depth is floored by `config.manning_h_min` to keep fb2 finite near the
# wet/dry front, and the resulting fb2 is capped at `config.manning_fb2_max`
# to prevent unphysical drag in very shallow flow.
# ---------------------------------------------------------------------------
@inline function _update_fb2_manning!(state::CshoreState, config::CshoreConfig,
                                      j::Int, l::Int)
    isempty(state.manning_n) && return
    n = state.manning_n[j, l]
    n <= 0.0 && return
    h_eff = max(state.h[j], config.manning_h_min)
    fb2_new = GRAV * n * n / cbrt(h_eff)
    if fb2_new > config.manning_fb2_max
        fb2_new = config.manning_fb2_max
    end
    state.fb2[j, l] = fb2_new
    return
end

@inline function _update_fb2_spatial!(state::CshoreState, config::CshoreConfig,
                                      j::Int, l::Int, ustd_j::Float64)
    ustd_j <= 1e-10 && return
    sed      = config.sediment
    f_min    = config.f_min
    gsgm1_d  = GRAV * submerged_sgm1(sed) * sed.d50   # g · (ρ_s/ρ_w - 1) · d50
    theta    = f_min * ustd_j^2 / gsgm1_d          # grain-roughness Shields
    theta_cr = sed.shield
    state.fb2[j, l] = if theta >= config.theta_sheet
        config.f_sheet
    elseif theta >= theta_cr
        f_min * (theta / theta_cr)^config.f_ripple_exp
    else
        f_min
    end
end
