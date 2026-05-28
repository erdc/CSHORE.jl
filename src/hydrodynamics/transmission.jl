# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
hydrodynamics/transmission.jl — Wave transmission across the crest into a
landward water body (FORTRAN IWTRAN=1).

Drives the back-side wave field and back-slope morphodynamic response for
submerged/low-crested breakwaters, dune-overtopping cases, and inlet/breach
scenarios. Three pieces:

1. `crest_transmission_coefficient(Hmo_c, Rc, B; ...)` — van der Meer &
   d'Angremond (1996) empirical KT for rubble/low-crested structures.
   Falls back to a clamped freeboard-ratio form for steep dune crests.

2. `_init_landward_wave_field!` — seeds Hrms / sigma / h / sxxsta at the
   crest's landward neighbor from the transmitted wave amplitude.

3. `transmission!` — orchestrates the landward wave march from `jcrest`
   to `jmax_l` using the existing `lwave!` and `dbreak!` infrastructure,
   referencing the landward SWL when supplied.

Out of scope for this MVP: explicit landward longshore current (treated as
zero on back side), nonlinear wave statistics on back side (Ursell/skewness
reset to zero), roller energy continuity across the crest (zeroed). These
are deferred for the Phase-2 follow-up.
==============================================================================#

"""
    crest_transmission_coefficient(Hmo_c, Rc, B; method=:dangremond_vandermeer,
                                                 alpha_slope=0.4) -> Float64

Empirical wave transmission coefficient `KT = Hmo_landward / Hmo_seaward` at
a crest of freeboard `Rc` (m, positive when crest above SWL), incident
significant wave height `Hmo_c` (m), and crest width `B` (m).

Methods:
- `:dangremond_vandermeer` (default) — d'Angremond, van der Meer, De Jong
  (1996) for low-crested breakwaters:
      KT = clamp(-0.4·Rc/Hmo + 0.64·(B/Hmo)^(-0.31) · (1 - exp(-0.5·ξ)),
                 0.075, 0.80)
  where ξ = surf-similarity parameter (taken from the local incident wave).
  Robust against negative freeboards (submerged crests).
- `:goda` — Goda (1969) for emergent vertical structures:
      KT = clamp(0.5·(1 - sin(π/(2·a)·(Rc/Hmo + b))), 0.0, 1.0)
  with a=2.6, b=0.8 (typical Japanese seawall values).
- `:freeboard_ratio` — simple linear fallback:
      KT = clamp(1 - 0.5·Rc/Hmo, 0.075, 0.85)

`alpha_slope` is the structure slope tan(α), used for the surf-similarity
parameter (default 0.4 ≈ 22° dune-face).
"""
function crest_transmission_coefficient(Hmo_c::Float64, Rc::Float64, B::Float64;
                                          method::Symbol=:dangremond_vandermeer,
                                          alpha_slope::Float64=0.4,
                                          tp::Float64=8.0)
    Hmo_c = max(Hmo_c, 1e-6)
    if method === :dangremond_vandermeer
        # Surf-similarity parameter ξ = tan(α) / sqrt(Hmo / L0); L0 = g·Tp²/(2π)
        L0 = GRAV * tp * tp / (2π)
        xi = alpha_slope / sqrt(Hmo_c / max(L0, 1e-6))
        kt = -0.4 * Rc / Hmo_c +
              0.64 * (B / Hmo_c)^(-0.31) * (1.0 - exp(-0.5 * xi))
        return clamp(kt, 0.075, 0.80)
    elseif method === :goda
        a = 2.6; b = 0.8
        arg = (π / (2a)) * (Rc / Hmo_c + b)
        return clamp(0.5 * (1.0 - sin(arg)), 0.0, 1.0)
    elseif method === :freeboard_ratio
        return clamp(1.0 - 0.5 * Rc / Hmo_c, 0.075, 0.85)
    else
        throw(ArgumentError("Unknown KT method: $method. Choose " *
            ":dangremond_vandermeer, :goda, or :freeboard_ratio."))
    end
end

"""
    transmission!(state, config, l)

Drive the wave field landward of the crest (IWTRAN=1). Computes the crest
transmission coefficient KT, seeds the wave field at jcrest+1 with
`Hrms = KT · Hrms[jcrest]`, then marches landward to `jmax_l` using the
existing `lwave!` machinery against the landward SWL (config.boundary.swl_landward
when provided, else the seaward swlbc).

State updated (landward of crest, j = jcrest+1 .. jmax_l):
- hrms, sigma, sigsta, h, wsetup, cp, wn, wkp_arr, wt
- sxxsta, sxysta, efsta
- dbsta, dfsta (zeroed for now — no breaking/friction integration on back side)

Idempotent within a single timestep; called by step_bc_window! after the
forward wave march completes when `options.iwtran == 1`.
"""
function transmission!(state::CshoreState, config::CshoreConfig, l::Int)
    opts = config.options
    opts.iwtran == 1 || return state
    jcrest = state.jcrest[l]
    jmax_l = state.jmax[l]
    jcrest >= jmax_l && return state   # nothing landward of crest

    # 1. Landward SWL — interpolate at current time (fallback to seaward)
    bc = config.boundary
    itime = max(1, min(state.itime, length(bc.timebc)))
    swl_l = if isempty(bc.swl_landward)
        bc.swlbc[itime]
    else
        @assert length(bc.swl_landward) == length(bc.timebc) "swl_landward must match timebc length"
        bc.swl_landward[itime]
    end

    # 2. Crest geometry — freeboard Rc and effective crest width B
    zb_crest = state.zb[jcrest, l]
    Rc = zb_crest - bc.swlbc[itime]   # seaward-side freeboard (m, +ve when emergent)
    B  = _effective_crest_width(state, l, jcrest, zb_crest)

    # 3. Transmission coefficient (d'Angremond/van der Meer)
    Hmo_c = state.hrms[jcrest] * sqrt(2.0)   # Hrms → Hmo conversion
    tp_c  = state.wt[jcrest]
    kt = crest_transmission_coefficient(Hmo_c, Rc, B;
                                         method=opts.iwtran_kt_method,
                                         tp=tp_c)

    # 4. Seed landward wave field at jcrest+1
    j0 = jcrest + 1
    hrms_trans = kt * state.hrms[jcrest]
    state.hrms[j0]    = hrms_trans
    state.sigma[j0]   = hrms_trans / SQR8
    h_l = max(swl_l - state.zb[j0, l], 0.0)
    state.h[j0]       = h_l
    state.wsetup[j0]  = 0.0   # no setup on landward water body in MVP
    state.swldep[j0, l] = h_l
    state.sigsta[j0]  = h_l > 1e-6 ? state.sigma[j0] / h_l : 0.0

    # Wave kinematics — re-solve dispersion at the new local depth
    if h_l > 1e-6
        wave = _current_wave_params(config, itime)
        lwave!(state, config, j0, l, h_l, wave; hrms_j=hrms_trans)
        sigma2 = state.sigma[j0]^2
        state.sxxsta[j0] = sigma2 * state.fsx
        state.efsta[j0]  = sigma2 * state.fe
        if opts.iangle == 1
            state.sxysta[j0] = sigma2 * state.fsy
        end
    end
    state.dbsta[j0] = 0.0
    state.dfsta[j0] = 0.0
    state.dvegsta[j0] = 0.0

    # 5. Landward march — j0+1 .. jmax_l. We use linear-amplitude propagation
    # with a friction-only attenuation rate (no breaking on the still landward
    # water body in the MVP; the transmitted wave is already amplitude-limited
    # by KT). The march stops when h ≤ eps1 (intersection with landward bed).
    eps1 = config.eps1
    @inbounds for j in (j0 + 1):jmax_l
        h_j = max(swl_l - state.zb[j, l], 0.0)
        state.h[j]      = h_j
        state.wsetup[j] = 0.0
        state.swldep[j, l] = h_j
        if h_j ≤ eps1
            # Bed surfaces — terminate landward march
            state.hrms[j]   = 0.0
            state.sigma[j]  = 0.0
            state.sigsta[j] = 0.0
            state.sxxsta[j] = 0.0
            state.efsta[j]  = 0.0
            state.dbsta[j]  = 0.0
            state.dfsta[j]  = 0.0
            state.dvegsta[j] = 0.0
            continue
        end

        # Shoaling from the previous landward node via constant energy flux
        # (Hrms²·cg conserved). Friction modifies this slightly via dfsta.
        h_prev = state.h[j - 1]
        if h_prev > 1e-6 && state.hrms[j - 1] > 1e-6
            # Update kinematics at j first to get cp[j], wn[j]
            wave = _current_wave_params(config, itime)
            lwave!(state, config, j, l, h_j, wave; hrms_j=state.hrms[j - 1])
            cg_prev = state.cp[j - 1] * state.wn[j - 1]
            cg_j    = state.cp[j]     * state.wn[j]
            ks_sq   = cg_prev / max(cg_j, 1e-9)
            hrms_j  = state.hrms[j - 1] * sqrt(max(ks_sq, 0.0))
            # Friction attenuation: simple Dean-style 1-h decay over one dx
            fb2_j   = state.fb2[j, l]
            dx      = config.grid.dx
            atten   = exp(-2.0 * fb2_j * dx / max(h_j, eps1))
            hrms_j *= atten
        else
            hrms_j = 0.0
        end

        state.hrms[j]   = hrms_j
        state.sigma[j]  = hrms_j / SQR8
        state.sigsta[j] = h_j > 1e-6 ? state.sigma[j] / h_j : 0.0
        # Refresh radiation-stress / energy-flux arrays from final state
        sigma2 = state.sigma[j]^2
        state.sxxsta[j] = sigma2 * state.fsx
        state.efsta[j]  = sigma2 * state.fe
        if opts.iangle == 1
            state.sxysta[j] = sigma2 * state.fsy
        end
        state.dbsta[j]  = 0.0
        state.dfsta[j]  = 0.0
        state.dvegsta[j] = 0.0
    end

    return state
end

# Estimate effective crest width — distance over which zb stays within
# `dz_tol` of the crest elevation. Falls back to one dx if the crest is sharp.
function _effective_crest_width(state::CshoreState, l::Int, jcrest::Int,
                                  zb_crest::Float64; dz_tol::Float64=0.10)
    nb = length(state.xb)
    j_lo = jcrest; j_hi = jcrest
    @inbounds while j_lo > 1 && (zb_crest - state.zb[j_lo - 1, l]) < dz_tol
        j_lo -= 1
    end
    @inbounds while j_hi < nb && (zb_crest - state.zb[j_hi + 1, l]) < dz_tol
        j_hi += 1
    end
    return max(state.xb[j_hi] - state.xb[j_lo], state.xb[2] - state.xb[1])
end
