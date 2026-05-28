# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
state.jl — Mutable simulation state.

Field names follow FORTRAN conventions for traceability.
Multi-fraction arrays (`bed_mass`, `ps`, `pb`, `qbx`, `qsx`, `concentration`)
have a trailing fraction dimension `k` even in single-grain mode (nf=1).
==============================================================================#

"""
    CshoreState

Mutable simulation state. Allocated once by `initialize_state(config)` and
mutated in-place by the time loop. Shape conventions:

- 1D node arrays: `Vector{Float64}` of length `nn` (wave/hydro state, per
  current line).
- 2D (nx, nlines): per-line fields that persist across the L loop
  (`zb`, `fb2`, `swldep`, `bslope`, `delzb`).
- 2D (nx, nf): per-fraction transport-rate / concentration state.
- 3D (nx, nlayers, nf): `bed_mass` — the multifraction bed composition.
"""
mutable struct CshoreState
    # ----- /BPROFL/ discretized bathymetry -----
    xb::Vector{Float64}          # (nn,) x-coordinates (shared across lines)
    zb::Matrix{Float64}          # (nn, nl) bed elevation — updated directly by
                                 # update_bed_composition! (mass-is-primary, but
                                 # zb is a concurrent tracker maintained to machine
                                 # precision via `zb += dm / (ρs(1-n))` per sub-step).
    fb2::Matrix{Float64}         # (nn, nl) friction factor (dynamic when
                                 # options.ifriction_spatial = 1 or 2)
    manning_n::Matrix{Float64}   # (nn, nl) Manning's n. Populated from
                                 # BathyInput.manning_n at init when supplied;
                                 # consumed when options.ifriction_spatial == 2
                                 # to drive fb2 = g·n²/h^(1/3) each step.
                                 # Empty matrix when Manning is not active.
    zb_hard::Matrix{Float64}     # (nn, nl) hardbottom elevation (erosion floor). Cells
                                 # with `zb_hard = -Inf` sentinel are unconstrained.
                                 # Used by apply_hardbottom_clamp! (ISEDAV≥1) and to
                                 # drive state.hp = max(0, zb - zb_hard) for the
                                 # per-fraction BRF supply limitation in transport.jl.
    swldep::Matrix{Float64}      # (nn, nl) still water depth
    bslope::Matrix{Float64}      # (nn, nl) bed slope dzb/dx
    jmax::Vector{Int}            # (nl,) max wet node index per line
    jswl::Vector{Int}            # (nl,) SWL intersection node

    # ----- /PREDIC/ wave variables -----
    hrms::Vector{Float64}        # (nn,)
    sigma::Vector{Float64}
    h::Vector{Float64}           # total depth (swl + setup - zb)
    wsetup::Vector{Float64}
    sigsta::Vector{Float64}

    # ----- /LINEAR/ linear wave + /PERIOD/ -----
    wkp::Float64                 # offshore wavenumber (scalar; updated each node by lwave!)
    cp::Vector{Float64}          # (nn,) phase speed
    wn::Vector{Float64}          # (nn,) Cg/Cp
    wkp_arr::Vector{Float64}     # (nn,) local wavenumber k at each node (stored in lwave!)
    wt::Vector{Float64}          # (nn,) local wave period per node (usually = tp)
    stheta::Vector{Float64}      # sin(theta)
    ctheta::Vector{Float64}      # cos(theta)
    wkpsin::Float64              # k*sin(theta) invariant for Snell
    fsx::Float64; fsy::Float64; fe::Float64
    qwx::Float64; qwy::Float64

    # ----- /FRICTN/ friction -----
    gbx::Vector{Float64}
    gby::Vector{Float64}
    gf::Vector{Float64}

    # ----- /WBREAK/ wave breaking -----
    qbreak::Vector{Float64}
    dbsta::Vector{Float64}
    abreak::Vector{Float64}

    # ----- Wave nonlinearity (Ruessink et al. 2012 parameterization) -----
    # Ursell number Ur = (3/8) Hs k / (kh)^3 quantifies shallow-water
    # nonlinearity. Skewness (crest>trough) drives onshore bedload;
    # asymmetry (sawtooth shape) drives acceleration-skewness transport.
    ursell::Vector{Float64}     # (nn,) dimensionless
    skewness::Vector{Float64}   # (nn,) Sk in [-0.86, 0.86]
    asymmetry::Vector{Float64}  # (nn,) As in [-0.86, 0.86]
    biphase::Vector{Float64}    # (nn,) evolved biphase β (rad); = static Ruessink when biphase_relax_L=0

    # ----- /CRSMOM/ + /LOGMOM/ + /ENERGY/ -----
    sxxsta::Vector{Float64}; tbxsta::Vector{Float64}
    sxysta::Vector{Float64}; tbysta::Vector{Float64}
    # ----- /WIND/ wind shear stress (per BC window) -----
    twxsta::Float64             # wind stress x-component (ρ_air/ρ_w · Cd · W10² · cos(θ) / g)
    twysta::Float64             # wind stress y-component (ρ_air/ρ_w · Cd · W10² · sin(θ) / g)
    efsta::Vector{Float64};  dfsta::Vector{Float64}
    dvegsta::Vector{Float64}  # /VEGDISS/

    # ----- /VELOCY/ horizontal velocities -----
    umean::Vector{Float64}; ustd::Vector{Float64}; usta::Vector{Float64}
    vmean::Vector{Float64}; vstd::Vector{Float64}; vsta::Vector{Float64}

    # ----- /IGWAVE/ infragravity wave state (IgConfig) -----
    # Populated by compute_ig_field! each BC window when config.ig !== nothing.
    # Zeroed when IgConfig is absent so all downstream code can read them
    # unconditionally (they simply contribute nothing when zero).
    hrms_ig::Vector{Float64}    # (nn,) RMS IG wave height (m)
    ustd_ig::Vector{Float64}    # (nn,) RMS IG orbital velocity (m/s) ≈ hrms_ig·√(g/h)/2

    # ----- Cohesive (mud) sediment (CohesiveSedimentConfig) -----
    # Tracked separately from the sand multifraction stack. Zero when
    # config.cohesive === nothing. Single line (l=1) for v1 — extend to
    # (nn, nl) when q2d cohesive coupling is needed.
    cohesive_bed_mass::Vector{Float64}      # (nn,) kg/m² of mud on the bed
    cohesive_concentration::Vector{Float64} # (nn,) kg/m³ depth-averaged suspended C

    # ----- /GWMODEL/ beach groundwater + surface moisture (GroundwaterConfig) -----
    # Populated by step_groundwater! each BC window when config.groundwater !== nothing.
    # gw_eta is initialized to SWL on the first call; theta is always ∈ [theta_res, theta_sat].
    # Both are zero when GroundwaterConfig is absent; step_aeolian! reads theta as the
    # physics-based moisture field, overriding the empirical ae_moisture decay.
    gw_eta::Matrix{Float64}     # (nn, nl) water table elevation above datum (m)
    theta::Matrix{Float64}      # (nn, nl) volumetric surface moisture (-), Van Genuchten
    # Forchheimer max seepage velocity (m/s); set by wetdry! each BC window for use
    # by step_groundwater! when infilt=1.  Zero when infilt=0.
    wpm_derived::Float64

    # ----- /RUNUP/ -----
    xr::Float64; zr::Float64; ssp::Float64; jr::Int

    # ----- /ROLLER/ -----
    rbzero::Float64
    rbeta::Vector{Float64}; rq::Vector{Float64}
    rx::Vector{Float64}; ry::Vector{Float64}; re::Vector{Float64}

    # ----- /SEDOUT/ sediment output — PER-FRACTION -----
    ps::Matrix{Float64}          # (nn, nf) suspension probability
    vs::Matrix{Float64}          # (nn, nf) suspended volume
    pb::Matrix{Float64}          # (nn, nf) bedload probability
    gslope::Vector{Float64}      # (nn,) bed slope correction (smoothed)
    aslope::Vector{Float64}      # (nn,) suspended-load slope correction (smoothed)
    qbx::Matrix{Float64}         # (nn, nf)
    qby::Matrix{Float64}
    qsx::Matrix{Float64}
    qsy::Matrix{Float64}
    q::Matrix{Float64}           # (nn, nf) per-fraction total transport (qbx+qsx)/sporo1
    q_total::Vector{Float64}     # (nn,) summed over fractions, smoothed — feeds CHANGE
    # Longshore volume flux summed over fractions: Σ_k(qby[:,k]+qsy[:,k])/sporo.
    # Populated by bmi_compute_qby_total! (called after each step_bc_window! when
    # the Q2D orchestrator is active).  Exposed via BMI "qby_total" for coupling.
    qby_total::Vector{Float64}   # (nn,) longshore transport for Q2D coupling

    # ----- /SEDVOL/ cumulative volumes per line -----
    vbx::Matrix{Float64}         # (nn, nl)
    vsx::Matrix{Float64}
    vby::Matrix{Float64}
    vsy::Matrix{Float64}
    vy::Matrix{Float64}
    dzx::Matrix{Float64}

    # ----- /PROCOM/ morph bookkeeping -----
    delt::Float64
    delzb::Matrix{Float64}       # (nn, nl)

    # ----- /POROUS/ porous flow state -----
    zp::Matrix{Float64}          # (nn, nl)
    hp::Matrix{Float64}
    upmean::Vector{Float64}
    upstd::Vector{Float64}
    dpsta::Vector{Float64}
    qp::Vector{Float64}
    upmwd::Vector{Float64}

    # ----- /DIKERO/ + /SOCLAY/ — EROSON state (IPROFL=2 or ICLAY=1) -----
    # Only allocated/used when iprofl=2 or iclay=1; otherwise zeros.
    edike::Matrix{Float64}       # (nn, nl) dike erosion depth (m)
    zb0::Matrix{Float64}         # (nn, nl) initial bed elevation
    dsta::Vector{Float64}        # (nn,) erosion forcing
    dsum::Vector{Float64}        # (nn,) cumulative forcing (IPROFL=2)
    bsf::Vector{Float64}         # (nn,) bed-slope amplification factor
    dfswd::Vector{Float64}       # (nn,) swash-zone dike forcing
    grsd::Matrix{Float64}        # (nn, nl) grass-layer thickness
    grsr::Matrix{Float64}        # (nn, nl) undamaged grass resistance
    grsrd::Matrix{Float64}       # (nn, nl) grass resistance at base
    grs1::Matrix{Float64}        # (nn, nl) derived dike constants
    grs2::Matrix{Float64}
    grs3::Matrix{Float64}
    grs4::Matrix{Float64}
    grs5::Matrix{Float64}
    fba3::Matrix{Float64}        # (nn, nl) friction constant (√g · AWD³ · fb2)
    # ICLAY=1:
    epclay::Matrix{Float64}      # (nn, nl) cumulative clay erosion depth
    zp0_clay::Matrix{Float64}    # (nn, nl) initial clay surface
    rclay::Matrix{Float64}       # (nn, nl) clay erosion rate coefficient
    fclay::Matrix{Float64}       # (nn, nl) effective sand fraction from clay
    eroson_initialized::Bool     # flag: time-invariant constants computed

    # ----- Vegetation uprooting state (IVEG=1,3 and IPROFL=1) -----
    # Dynamic canopy height that evolves with bed change. Initially copied
    # from VegetationInput.vegh; shrinks when bed buries vegetation,
    # zeroed when roots are exposed.
    vegh::Matrix{Float64}        # (nn, nl) effective canopy height above current bed
    vegzd::Matrix{Float64}       # (nn, nl) fixed upper elevation (zb0 + vegh_initial)
    vegzr::Matrix{Float64}       # (nn, nl) fixed lower elevation (zb0 - vegrh/vegrd)
    uproot::Matrix{Float64}      # (nn, nl) 1.0=active vegetation, 0.0=uprooted/none

    # ----- /TIDALC/ tidal forcing (ITIDE=1) -----
    # DETADY (dimensionless alongshore SWL gradient) and DSWLDT (m/s rate of
    # SWL change) are interpolated per BC window from the TidalInput time
    # series. QTIDE is the per-node cross-shore tidal volume flux used when
    # ILAB=0 to augment QWX = QO + QTIDE.
    detady_now::Float64          # alongshore SWL gradient at current time
    dswldt_now::Float64          # SWL rate of change at current time (m/s)
    qtide::Vector{Float64}       # (nn,) cross-shore tidal flux per node (m²/s)

    # ----- ICURRENT=1 imposed alongshore current at offshore boundary -----
    # When `options.icurrent == 1`, `vbc_now` is the user-prescribed alongshore
    # current speed (m/s) at the offshore boundary, interpolated from the
    # `CurrentInput.vbc` time series each BC window. The driver inverts the
    # longshore-momentum balance at the most-offshore valid cell to derive
    # `detady_now` so that `vmean` at that cell equals `vbc_now`.
    vbc_now::Float64             # imposed alongshore current at offshore boundary (m/s)
    icurrent_solved::Vector{Bool}  # (nl,) per-line latch: true once detady_now
                                   # has been back-solved this BC window

    # ----- Aeolian (wind-driven) sediment transport (IAEOLIAN=1) ------------
    # Per-cell cross-shore aeolian flux. Sign convention: positive = transport
    # in the +x direction (landward), negative = transport in -x (seaward).
    # Multifraction: per-fraction per-cell flux; total summed in `qae`.
    qae::Vector{Float64}                 # (nn,) total aeolian flux (m²/s, sand vol)
    qae_k::Matrix{Float64}               # (nn, nf) per-fraction aeolian flux (m²/s)
    ust::Vector{Float64}                 # (nn,) current friction velocity (m/s) — diag.
    ust_threshold::Matrix{Float64}       # (nn, nf) per-fraction threshold u*_t (m/s)
    ae_dz::Vector{Float64}               # (nn,) bed-elevation rate from aeolian (m/s)
    ae_loss_landward::Float64            # cumulative volume exiting landward (m³/m)
    ae_loss_seaward::Float64             # cumulative volume exiting seaward  (m³/m)
    # Per-cell time of most recent wetting (s, model time). Initialized to
    # -Inf so initially-dry cells start at full configured moisture (or
    # zero moisture if cfg.dry_time = 0). Updated each aeolian step:
    # cells currently below SWL or within the swash band have their
    # timestamp set to state.time; cells above the runup line keep their
    # previous timestamp. When IAEOLIAN=1 with `cfg.dry_time > 0`, the
    # per-cell moisture content M(j) decays from 1 (just wetted) toward
    # `cfg.moisture_dry` over `cfg.dry_time` seconds.
    ae_t_last_wet::Vector{Float64}       # (nn,) seconds (model time) of last wetting
    ae_moisture::Vector{Float64}         # (nn,) current surface moisture M ∈ [0,1] (diag.)

    # ----- Wind-flow-over-topography solver (IWINDSHEAR=1) ------------------
    # Per-cell Kroy-Sauermann-Herrmann shear-stress perturbation factor
    # τ/τ₀(x). With IWINDSHEAR=0, all entries stay at 1.0 (no perturbation).
    # The aeolian kernel multiplies u* by √(τ/τ₀) per cell when active.
    tau_perturbation::Vector{Float64}    # (nn,) τ/τ₀ at each cell
    lee_zone::Vector{Bool}               # (nn,) true inside lee separation bubbles

    # ----- Iterative wave-current interaction (IWCINT_ALONG=1) ---------------
    # Alongshore-current contribution to the Doppler-shifted dispersion
    # `(σ - k⃗·U⃗)² = gk·tanh(kh)` requires `vmean` to be known when LWAVE
    # solves at each node, but vmean is only available *after* transform_waves!
    # finishes. We resolve this with an outer Picard iteration: each pass
    # uses `vmean_prev` from the previous pass to construct an effective
    # qdisp = qwx + vmean_prev·sin(θ)·h. Converges in 2-4 iterations for
    # typical longshore/oblique-wave combinations.
    vmean_prev::Vector{Float64}    # (nn,) vmean from the previous outer
                                   # iteration of transform_waves!

    # ----- /OVERTF/ overtopping -----
    qo::Vector{Float64}          # (nl,) overtopping rate per line
    qotf::Float64; sprate::Float64; slpot::Float64
    jcrest::Vector{Int}          # (nl,) crest node per line

    # ----- /SWASHY/ wet/dry swash -----
    pwet::Vector{Float64}
    uswd::Vector{Float64}; hwd::Vector{Float64}; sigwd::Vector{Float64}
    umeawd::Vector{Float64}; ustdwd::Vector{Float64}
    vmeawd::Vector{Float64}; vstdwd::Vector{Float64}
    hewd::Vector{Float64}; uewd::Vector{Float64}; qewd::Vector{Float64}
    h1::Float64
    jwd::Int; jdry::Int

    # ----- Multifraction bed composition -----
    bed_mass::Array{Float64,3}   # (nn, nlayers, nf) kg/m²
    thlyr::Matrix{Float64}       # (nn, nlayers)
    concentration::Matrix{Float64}  # (nn, nf)
    pickup_fractions::Matrix{Float64}  # (nn, nf) per-fraction -dQ/dx, m/s — for grain sorting
    active_frac::Matrix{Float64}  # (nn, nf) virtual adaptive active-layer fractions (Hrms-weighted blend of top layers)

    # ----- Per-fraction scratch (grain-size-dependent coefficients) -----
    ws_fractions::Vector{Float64}       # (nf,) fall velocity per fraction
    gsd50s_fractions::Vector{Float64}   # (nf,) sqrt((s-1)*g*d_k)
    tadapt_fractions::Vector{Float64}   # (nf,) adaptation time
    theta_cr_fractions::Vector{Float64} # (nf,) critical Shields per fraction

    # ----- /SERIES/ time-series accumulators -----
    tsqo::Vector{Float64}; tsqbx::Vector{Float64}; tsqsx::Vector{Float64}

    # ----- Current step bookkeeping -----
    time::Float64
    itime::Int

    # Enforces that the supply factor is applied exactly once per step
    supply_factor_applied::Bool

    # ----- Thermal / permafrost submodel (optional) -----
    # Populated by `initialize_state` when `config.thermal !== nothing`.
    thermal::Union{Nothing, ThermalState}
end

"""
    initialize_state(config::CshoreConfig)

Allocate all mutable arrays for a simulation. Zeros-out state; `run_simulation!`
calls `apply_initial_conditions!` to populate bed_mass, zb, and fb2 from the
input bathymetry.
"""
function initialize_state(config::CshoreConfig)
    validate(config)
    nn = config.grid.nn
    nl = config.grid.nl
    nf = nfractions(config.multifraction)
    nlayers = config.multifraction.nlayers

    zeros_nn() = zeros(Float64, nn)
    zeros_nn_nl() = zeros(Float64, nn, nl)
    zeros_nn_nf() = zeros(Float64, nn, nf)

    state = CshoreState(
        # BPROFL
        zeros_nn(),                 # xb
        zeros_nn_nl(),              # zb
        zeros_nn_nl(),              # fb2
        zeros(Float64, 0, 0),       # manning_n (filled by apply_initial_conditions! when supplied)
        fill(-Inf, nn, nl),         # zb_hard — sentinel: no hardbottom constraint
        zeros_nn_nl(),              # swldep
        zeros_nn_nl(),              # bslope
        zeros(Int, nl),             # jmax
        zeros(Int, nl),             # jswl

        # PREDIC
        zeros_nn(), zeros_nn(), zeros_nn(), zeros_nn(), zeros_nn(),

        # LINEAR + PERIOD
        0.0, zeros_nn(), zeros_nn(), zeros_nn(), zeros_nn(),
        zeros_nn(), zeros_nn(),
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,

        # FRICTN
        zeros_nn(), zeros_nn(), zeros_nn(),

        # WBREAK
        zeros_nn(), zeros_nn(), zeros_nn(),

        # Wave nonlinearity (Ursell, skewness, asymmetry, biphase)
        zeros_nn(), zeros_nn(), zeros_nn(), zeros_nn(),

        # CRSMOM/LOGMOM/ENERGY/VEGDISS
        zeros_nn(), zeros_nn(), zeros_nn(), zeros_nn(),
        0.0, 0.0,    # twxsta, twysta (wind stress)
        zeros_nn(), zeros_nn(), zeros_nn(),

        # VELOCY
        zeros_nn(), zeros_nn(), zeros_nn(), zeros_nn(), zeros_nn(), zeros_nn(),

        # IGWAVE
        zeros_nn(), zeros_nn(),   # hrms_ig, ustd_ig

        # COHESIVE (mud) — bed_mass seeded by initialize_state below.
        zeros_nn(), zeros_nn(),   # cohesive_bed_mass, cohesive_concentration

        # GWMODEL
        zeros_nn_nl(),            # gw_eta  (seeded from SWL on first step_groundwater! call)
        zeros_nn_nl(),            # theta   (Van Genuchten moisture; 0 until first GW step)
        0.0,                      # wpm_derived

        # RUNUP
        0.0, 0.0, 0.0, 0,

        # ROLLER
        0.0, zeros_nn(), zeros_nn(), zeros_nn(), zeros_nn(), zeros_nn(),

        # SEDOUT (per-fraction)
        zeros_nn_nf(),              # ps
        zeros_nn_nf(),              # vs
        zeros_nn_nf(),              # pb
        zeros_nn(),                 # gslope
        zeros_nn(),                 # aslope
        zeros_nn_nf(),              # qbx
        zeros_nn_nf(),              # qby
        zeros_nn_nf(),              # qsx
        zeros_nn_nf(),              # qsy
        zeros_nn_nf(),              # q
        zeros_nn(),                 # q_total
        zeros_nn(),                 # qby_total

        # SEDVOL
        zeros_nn_nl(), zeros_nn_nl(), zeros_nn_nl(),
        zeros_nn_nl(), zeros_nn_nl(), zeros_nn_nl(),

        # PROCOM
        0.0, zeros_nn_nl(),

        # POROUS
        zeros_nn_nl(), zeros_nn_nl(),
        zeros_nn(), zeros_nn(), zeros_nn(), zeros_nn(), zeros_nn(),

        # DIKERO / SOCLAY — EROSON state
        zeros_nn_nl(),             # edike
        zeros_nn_nl(),             # zb0
        zeros_nn(),                # dsta
        zeros_nn(),                # dsum
        zeros_nn(),                # bsf
        zeros_nn(),                # dfswd
        zeros_nn_nl(),             # grsd
        zeros_nn_nl(),             # grsr
        zeros_nn_nl(),             # grsrd
        zeros_nn_nl(),             # grs1
        zeros_nn_nl(),             # grs2
        zeros_nn_nl(),             # grs3
        zeros_nn_nl(),             # grs4
        zeros_nn_nl(),             # grs5
        zeros_nn_nl(),             # fba3
        zeros_nn_nl(),             # epclay
        zeros_nn_nl(),             # zp0_clay
        zeros_nn_nl(),             # rclay
        zeros_nn_nl(),             # fclay
        false,                     # eroson_initialized

        # Vegetation uprooting state
        zeros_nn_nl(),             # vegh
        zeros_nn_nl(),             # vegzd
        zeros_nn_nl(),             # vegzr
        zeros_nn_nl(),             # uproot

        # TIDALC — tidal gradient forcing
        0.0,                       # detady_now
        0.0,                       # dswldt_now
        zeros_nn(),                # qtide

        # ICURRENT — imposed alongshore current at offshore boundary
        0.0,                       # vbc_now
        falses(nl),                # icurrent_solved (per-line latch)

        # IAEOLIAN — wind-driven sediment transport diagnostics & per-fraction flux
        zeros_nn(),                # qae   (total flux, m²/s)
        zeros_nn_nf(),             # qae_k (per-fraction flux, m²/s)
        zeros_nn(),                # ust   (friction velocity, m/s)
        zeros_nn_nf(),             # ust_threshold (per-fraction, m/s)
        zeros_nn(),                # ae_dz (bed change rate, m/s)
        0.0,                       # ae_loss_landward (cumulative, m³/m)
        0.0,                       # ae_loss_seaward  (cumulative, m³/m)
        fill(-Inf, nn),            # ae_t_last_wet (cells start "long-dry")
        zeros_nn(),                # ae_moisture (diagnostic)
        ones(Float64, nn),         # tau_perturbation (no-perturbation default)
        falses(nn),                # lee_zone

        # IWCINT_ALONG — iterative alongshore-current WCI scratch
        zeros_nn(),                # vmean_prev

        # OVERTF
        zeros(Float64, nl), 0.0, 0.0, 0.0, zeros(Int, nl),

        # SWASHY
        zeros_nn(), zeros_nn(), zeros_nn(), zeros_nn(),
        zeros_nn(), zeros_nn(), zeros_nn(), zeros_nn(),
        zeros_nn(), zeros_nn(), zeros_nn(),
        0.0, 0, 0,

        # Multifraction
        zeros(Float64, nn, nlayers, nf),
        fill(config.multifraction.layer_thickness, nn, nlayers),
        zeros_nn_nf(),
        zeros_nn_nf(),     # pickup_fractions
        zeros_nn_nf(),     # active_frac

        # Per-fraction coefficients
        zeros(Float64, nf), zeros(Float64, nf),
        zeros(Float64, nf), zeros(Float64, nf),

        # SERIES
        zeros(Float64, nl), zeros(Float64, nl), zeros(Float64, nl),

        # Bookkeeping
        0.0, 0,

        # Supply factor assertion flag
        false,

        # Thermal state (populated later if config.thermal !== nothing)
        nothing,
    )

    # Precompute per-fraction grain quantities.
    #
    # `theta_cr_fractions[k]` — per-fraction critical Shields parameter.
    # When `use_size_dependent_shields` is true (default, physically
    # correct for graded beds), this uses the Soulsby (1997) curve which
    # varies with d_k. When false, it falls back to the constant
    # `config.sediment.shield` (single-grain behavior with SHIELD=0.05
    # for all fractions).
    #
    # `gsd50s_fractions[k]` — velocity² representing the critical shear
    # stress for fraction k: `(s-1)*g*d_k*θ_cr_k`. For single-grain
    # parity, FORTRAN uses a constant 0.05; for multifraction physical
    # correctness we use the per-fraction Soulsby value. SEDTRA uses
    # `rb_k = sqrt(gsd50s_k/fb2)/ustd` as the grain-specific initiation
    # parameter, so per-fraction `θ_cr_k` gives per-fraction threshold
    # behavior — this is what produces grain-size discrimination in the
    # surf zone that the single-grain formula loses.
    sgm1 = submerged_sgm1(config.sediment)
    for k in 1:nf
        d_k = config.multifraction.grain_sizes[k]
        state.ws_fractions[k]       = fall_velocity(d_k, config)
        state.theta_cr_fractions[k] = critical_shields(d_k, config)
        state.gsd50s_fractions[k]   = sgm1 * GRAV * d_k * state.theta_cr_fractions[k]
        state.tadapt_fractions[k]   = 1.0  # placeholder; recomputed in driver with mean h
    end

    # Seed cohesive (mud) bed mass with uniform initial_bed_mass when configured.
    # Callers that want spatial variation should overwrite state.cohesive_bed_mass
    # in place between initialize_state and the first transport step.
    if config.cohesive !== nothing
        fill!(state.cohesive_bed_mass, config.cohesive.initial_bed_mass)
    end

    return state
end

"""
    apply_initial_bathymetry!(state, config)

Populate `state.zb`, `state.xb`, `state.fb2` from `config.bathymetry` by
interpolating raw input points onto the uniform grid. Mirrors FORTRAN
SUBROUTINE BOTTOM.
"""
function apply_initial_bathymetry!(state::CshoreState, config::CshoreConfig)
    dx = config.grid.dx
    nl = config.options.iline
    bi = config.bathymetry

    has_hardbottom = !isempty(bi.zhinp)
    has_manning    = !isempty(bi.manning_n)
    if has_manning && isempty(state.manning_n)
        state.manning_n = zeros(Float64, config.grid.nn, nl)
    end
    for l in 1:nl
        np_l = bi.nbinp[l]
        x_in = @view bi.xbinp[1:np_l, l]
        z_in = @view bi.zbinp[1:np_l, l]
        f_in = @view bi.fbinp[1:np_l, l]
        zh_in = has_hardbottom ? (@view bi.zhinp[1:np_l, l]) : nothing
        n_in  = has_manning    ? (@view bi.manning_n[1:np_l, l]) : nothing

        x0 = x_in[1]
        xN = x_in[end]
        nnodes = floor(Int, (xN - x0) / dx) + 1
        nnodes ≤ config.grid.nn ||
            error("Initial bathymetry requires $nnodes nodes but grid.nn=$(config.grid.nn)")

        for j in 1:nnodes
            xj = x0 + (j - 1) * dx
            state.xb[j] = xj
            state.zb[j, l] = interp1(x_in, z_in, xj)
            state.fb2[j, l] = interp1(x_in, f_in, xj)
            if has_manning
                state.manning_n[j, l] = interp1(x_in, n_in, xj)
            end
            if has_hardbottom
                zh = interp1(x_in, zh_in, xj)
                # Semantics:
                #   zh ≤ z0 → rock buried under a sand veneer of thickness (z0 - zh).
                #             Initial zb stays at the sand surface; erosion thins
                #             the veneer until hp → 0 and BRF shuts off transport.
                #   zh > z0 → rock protrudes above the input bathymetry (e.g. an
                #             exposed outcrop or submerged breakwater). The bed
                #             surface IS the rock top in that region, so we
                #             clamp zb up to zh. No sand is "created" — there was
                #             no sand above the rock to begin with. The Aeolis
                #             bed_mass stack is left at its default composition
                #             (it's latent — BRF=0 there so nothing touches it
                #             until deposition raises zb above zh).
                if state.zb[j, l] < zh
                    state.zb[j, l] = zh
                end
                state.zb_hard[j, l] = zh
            end
        end
        state.jmax[l] = nnodes
    end

    # Initialize bed_mass from initial fractions
    init_bed_mass!(state, config)

    # ZB0: initial bed surface (used by IPROFL=2 to compute cumulative erosion)
    if config.options.iprofl >= 1
        @inbounds for l in 1:config.options.iline, j in 1:state.jmax[l]
            state.zb0[j, l] = state.zb[j, l]
        end
    end
    # ZP0_CLAY + RCLAY + FCLAY (ICLAY=1): interpolate from ClayInput onto grid
    if config.options.iclay == 1 && config.clay !== nothing
        bi = config.bathymetry
        ci = config.clay
        sporo1 = config.sediment.sporo1
        for l in 1:config.options.iline
            np_l = bi.nbinp[l]
            x_in = @view bi.xbinp[1:np_l, l]
            # Porous (clay surface) bathymetry
            has_por = config.porous !== nothing && size(config.porous.zpinp, 1) > 0
            if has_por
                np_por = config.porous.npinp[l]
                xp_in = @view config.porous.xpinp[1:np_por, l]
                zp_in = @view config.porous.zpinp[1:np_por, l]
            end
            np_ci = size(ci.rcinp, 1)
            rc_in = @view ci.rcinp[1:np_ci, l]
            fc_in = @view ci.fcinp[1:np_ci, l]
            for j in 1:state.jmax[l]
                xj = state.xb[j]
                # Initial clay surface: from porous input if provided, else zb_hard
                zp_j = has_por ? interp1(xp_in, zp_in, xj) : state.zb_hard[j, l]
                state.zp0_clay[j, l] = zp_j
                state.zp[j, l] = zp_j
                # Clay resistance: FORTRAN computes RCLAY = g / RCINP
                rcinp_j = interp1(x_in, rc_in, xj)
                state.rclay[j, l] = rcinp_j > 0.0 ? GRAV / rcinp_j : 0.0
                # Sand-fraction-in-clay effective factor (FCLAY = 1 - FCINP/(1-n))
                fcinp_j = interp1(x_in, fc_in, xj)
                state.fclay[j, l] = 1.0 - fcinp_j / max(sporo1, 1e-6)
            end
        end
    end
    # GRSD/GRSR/GRSRD (IPROFL=2): interpolate DikeErosionInput onto grid
    if config.options.iprofl == 2 && config.dike !== nothing
        bi = config.bathymetry
        di = config.dike
        for l in 1:config.options.iline
            np_l = bi.nbinp[l]
            x_in = @view bi.xbinp[1:np_l, l]
            np_di = size(di.gdinp, 1)
            gd_in = @view di.gdinp[1:np_di, l]
            gr_in = @view di.grinp[1:np_di, l]
            grd_in = @view di.grdinp[1:np_di, l]
            for j in 1:state.jmax[l]
                xj = state.xb[j]
                state.grsd[j, l]  = interp1(x_in, gd_in, xj)
                state.grsr[j, l]  = interp1(x_in, gr_in, xj)
                state.grsrd[j, l] = interp1(x_in, grd_in, xj)
            end
        end
    end

    return state
end

"""
    init_bed_mass!(state, config)

Populate the (nn, nlayers, nf) bed_mass array from uniform initial fractions.
mass[i,j,k] = ρ_s * (1 - porosity) * layer_thickness * initial_fraction[k]
"""
function init_bed_mass!(state::CshoreState, config::CshoreConfig)
    mf = config.multifraction
    ρs = config.sediment.sg * 1000.0   # kg/m³
    nlayers = mf.nlayers
    nf = nfractions(mf)
    base_mass = ρs * (1 - mf.porosity) * mf.layer_thickness

    has_spatial = size(mf.initial_fractions_spatial, 1) > 0

    if has_spatial
        # Spatial fractions: (n_bathy_nodes, nf) at raw bathymetry resolution.
        # Interpolate onto the uniform grid via the bathymetry x-coordinates.
        bi = config.bathymetry
        for l in 1:config.options.iline
            np_l = bi.nbinp[l]
            x_in = @view bi.xbinp[1:np_l, l]
            n_spatial = size(mf.initial_fractions_spatial, 1)
            for j in 1:state.jmax[l]
                xj = state.xb[j]
                for ilay in 1:nlayers
                    # Interpolate each fraction from raw nodes to grid node j
                    frac_sum = 0.0
                    for k in 1:nf
                        fk = if n_spatial == np_l
                            interp1(x_in, view(mf.initial_fractions_spatial, :, k), xj)
                        else
                            # Fallback: direct indexing if sizes match grid
                            j ≤ n_spatial ? mf.initial_fractions_spatial[j, k] : mf.initial_fractions[k]
                        end
                        fk = max(fk, 0.0)
                        state.bed_mass[j, ilay, k] = base_mass * fk
                        frac_sum += fk
                    end
                    # Renormalize to ensure mass conservation
                    if frac_sum > 0 && abs(frac_sum - 1.0) > 1e-10
                        for k in 1:nf
                            state.bed_mass[j, ilay, k] *= 1.0 / frac_sum
                        end
                    end
                end
            end
        end
    else
        # Uniform fractions everywhere (original behavior)
        for l in 1:config.options.iline
            for j in 1:state.jmax[l]
                for ilay in 1:nlayers
                    for k in 1:nf
                        state.bed_mass[j, ilay, k] = base_mass * mf.initial_fractions[k]
                    end
                end
            end
        end
    end
    return state
end
