# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
config.jl — Immutable configuration structs.
==============================================================================#

"""
    OptionFlags

Runtime toggles set once at input time. Features marked (OUT-OF-SCOPE) are
retained as fields for input-file parity but must be 0 in the current Julia port.
"""
Base.@kwdef struct OptionFlags
    iprofl::Int = 1         # 0=fixed bed, 1=morphodynamic, 2=dike erosion (OUT-OF-SCOPE)
    iangle::Int = 0         # 0=normal incidence, 1=oblique
    iroll::Int  = 0         # 0=no roller, 1=roller model
    iwind::Int  = 0         # 0=no wind, 1=wind stress
    iperm::Int  = 0         # 0=impermeable, 1=porous layer (POFLOW)
    iover::Int  = 0         # 0=single-pass swash (qo fixed at 0), 1=iterate wave+swash flux to convergence
    iwcint::Int = 0         # wave-current interaction
    isedav::Int = 0         # 0=unlimited sand, 1=hardbottom (BEDLM=1), -1=hardbottom (BEDLM=0), 2=wire mesh (OUT-OF-SCOPE)
    iwtran::Int = 0         # landward water transmission: 0=off, 1=transmit waves across the
                            # crest into a landward water body. See `hydrodynamics/transmission.jl`.
                            # When 1, optionally supply `boundary.swl_landward` time series;
                            # otherwise the landward SWL falls back to seaward swlbc.
    iwtran_kt_method::Symbol = :dangremond_vandermeer  # transmission coefficient model.
                            # :dangremond_vandermeer (default) — d'Angremond/van der Meer 1996;
                            # :goda — Goda 1969; :freeboard_ratio — simple linear fallback.
    ivwall::Vector{Int} = Int[]   # per-line vertical wall flag (NL)
    ilab::Int   = 0         # lab vs field flag
    infilt::Int = 0         # infiltration
    ipond::Int  = 0         # ridge-runnel ponded water (OUT-OF-SCOPE)
    itide::Int  = 0         # 0=no tide, 1=tidal alongshore gradient + cross-shore SWL rate
    iline::Int  = 1         # number of cross-shore lines (NL)
    iqydy::Int  = 0         # alongshore gradient coupling
    iveg::Int   = 0         # 0=none, 1=phase-averaged drag, 2=phase-resolved, 3=full (DVEG)
    iclay::Int  = 0         # clay resistance (OUT-OF-SCOPE)
    ismooth::Int = 1        # smoothing filter flag
    idiss::Int  = 0         # dissipation formulation (0=default, 3=measured spectrum OUT-OF-SCOPE)
    ifv::Int    = 0         # flow velocity option
    iweibull::Int = 0       # Weibull wave stats (OUT-OF-SCOPE)
    iasym::Int  = 0         # Wave nonlinearity correction for bedload:
                            # 0 = linear waves (no correction)
                            # 1 = Ruessink 2012 (empirical, calibrated)
                            # 2 = Stokes 2nd-order (zero-calibration)
    iwcint_along::Int = 0   # 0 = WCI uses cross-shore flux only;
                            # 1 = adds the alongshore-current Doppler contribution
                            #     vmean·sin(θ)·h to qdisp via an outer Picard
                            #     iteration. Requires IWCINT=1 + IANGLE=1, and
                            #     is only meaningful when an alongshore current
                            #     is present (e.g. ICURRENT=1, ITIDE=1, or wave-
                            #     driven longshore current with oblique waves).
    iv_transport::Int = 0   # 0 = transport ignores alongshore current vmean;
                            # 1 = vmean contributes to bed-shear entrainment, and when V
                            #     dominates over U_cross + wave action the transport direction
                            #     blends toward the downhill bed-slope direction (slope cascade).
                            # Currently honored by SizeAdaptiveTransport and SvrTransport kernels.
                            # Use with ICURRENT=1 / ITIDE=1 / iangle=1 to get current-driven
                            # bank erosion in fluvial / tidal-channel cases.
    iwindshear::Int = 0     # 0 = uniform u* over the dry beach;
                            # 1 = solve the Kroy-Sauermann-Herrmann shear-stress
                            # perturbation over the local topography, mask out the
                            # wet zone (max(zb, swl)), and apply a lee-separation
                            # bubble downwind of brinks. Outputs τ/τ₀(x) per cell;
                            # the aeolian kernel multiplies u* by √(τ/τ₀). Free
                            # cost (O(N²) Hilbert convolution). Requires WindShearConfig.
    # ---- scaffolding flags (full implementation pending) ----
    iwave_aeolian_coupling::Int = 0   # 0 = off; 1 (planned) = couple swash/overwash
                                      # to aeolian: when waves overwash a dune, the
                                      # newly-exposed wet sand resets the moisture
                                      # clock and becomes aeolian-active at the
                                      # next dry tide. Currently emits a warning.
    iwind_threshold_extras::Int = 0   # 0 = Belly-Johnson moisture only (current);
                                      # 1 (planned) = AEOLIS-style threshold stack:
                                      # sheltering by non-erodible (coarse) fractions
                                      # via Raupach (1993), salt content (Nickling &
                                      # Ecclestone), bed slope (Dyer 1986). Currently
                                      # emits a warning.
    iveg_dynamics::Int = 0            # 0 = static vegetation density (current);
                                      # 1 (planned) = vegetation density evolves
                                      # with sand burial (Buckley/Durán) and
                                      # natural growth/dieback. Currently emits a
                                      # warning.
    iaeolian::Int = 0       # 0 = no aeolian transport; 1 = enable wind-driven cross-shore
                            # transport on the dry beach (DRT-style Kawamura) with optional
                            # contour-based vegetation capture. Requires AeolianConfig and a wind
                            # time series (BoundaryTimeSeries.w10, .wangle).
                            # Currently uses :size_adaptive multifraction integration —
                            # per-fraction Kawamura threshold gives armouring "for free".
    icurrent::Int = 0       # 0 = none; 1 = imposed alongshore current at offshore boundary.
                            # When 1, the alongshore water-surface gradient (DETADY) is
                            # back-calculated from the user-prescribed `vbc` time series via
                            # analytical inversion of the longshore-momentum balance at the
                            # most-offshore valid cell. Compatible with ITIDE=1 (which then
                            # only handles SWL-rate forcing). Requires IANGLE=1 and a
                            # CurrentInput passed via `current=`.
    inl_dispersion::Int = 0 # 0 = linear dispersion (default);
                            # 1 = Kirby & Dalrymple (1986) finite-amplitude
                            #     correction to k: ω²=g·k·tanh(kh)·[1+(ka)²·f(kh)]
                            #     where f(kh)=(8+cosh(4kh)-2tanh²(kh))/(8sinh⁴(kh)).
                            #     Gives longer wavelength (k_nl < k_lin) → larger Ursell,
                            #     stronger Sk/As, and shifts the steepness-based breaking
                            #     criterion toward the physically correct location.
    iskew_spatial::Int = 0  # 0 = global facSK / facAS scalars (default);
                            # 1 = Ursell-weighted local coupling: the skewness/asymmetry
                            #     contribution is modulated by tanh(Ur / ur_sk_ref), so
                            #     coupling is strong in the nonlinear surf zone and tapers
                            #     to zero in linear (offshore) wave conditions. Avoids
                            #     applying surf-zone-tuned facSK in the shoaling region.
    ifriction_spatial::Int = 0  # 0 = constant fb2 from input bathymetry (default);
                                # 1 = update fb2 at each grid node each step based on the
                                #     local Shields-parameter regime:
                                #     θ < θ_cr          → f_min  (grain roughness only)
                                #     θ_cr ≤ θ < θ_sheet → f_min·(θ/θ_cr)^f_ripple_exp
                                #     θ ≥ θ_sheet        → f_sheet (sheet-flow cap)
                                #     Grain-roughness θ is computed from f_min to break the
                                #     circular dependency on fb2. Tune f_min, f_sheet,
                                #     theta_sheet, and f_ripple_exp in CshoreConfig.
                                # 2 = dynamic Manning: recompute fb2 = g·n²/h^(1/3) each
                                #     step from state.h, using bathymetry.manning_n (set
                                #     via build_config(manning=...)). Quadratic-drag
                                #     interpretation (fb2 ≡ Cf where τ_bed = Cf·ρ·U²);
                                #     h is floored by config.manning_h_min to avoid
                                #     singularity at the wet/dry front.
end

"""
    GridConfig

Grid and discretization parameters.
"""
Base.@kwdef struct GridConfig
    dx::Float64              # cross-shore grid spacing (m)
    nn::Int = 20000          # max nodes
    nb::Int = 30000          # max BC time points
    nl::Int = 1              # number of cross-shore lines
end

"""
    SedimentConfig

Single representative grain parameters that feed the legacy single-grain
sediment transport formulas. When `MultifractionConfig.nf > 1`, `d50` is used
only as a fallback; per-fraction quantities take precedence.

Derived quantities (`wfsgm1`, `gsgm1`, `gsd50s`, `bld`, `csedia`) follow the
FORTRAN SUBROUTINE INPUT post-read calculations. Build one with the
`SedimentConfig(; ...)` keyword constructor — the derived fields default to
zero and are computed in the outer helper below.
"""
Base.@kwdef struct SedimentConfig
    wf::Float64      = 0.026       # settling velocity (m/s)
    sg::Float64      = 2.65        # sediment specific gravity (vs freshwater @ 4°C, 1000 kg/m³)
    rho_water::Float64 = 1025.0    # surrounding water density (kg/m³). Default 1025 (seawater);
                                   # set to 1000.0 for freshwater simulations. Drives the
                                   # submerged specific gravity used in Shields, bedload,
                                   # settling velocity, and suspended-load formulas via the
                                   # `submerged_sgm1` helper. Grain density = sg * 1000 always
                                   # (sg is by definition referenced to freshwater).
    sporo1::Float64  = 0.6         # 1 - porosity
    tanphi::Float64  = 0.63        # tan of friction angle (~32°)
    bslop1::Float64  = 0.0
    bslop2::Float64  = 0.0
    effb::Float64    = 0.005       # breaking dissipation efficiency
    efff::Float64    = 0.005       # friction dissipation efficiency
    d50::Float64     = 0.3e-3      # median grain size (m)
    shield::Float64  = 0.05        # critical Shields parameter
    blp::Float64     = 2e-3        # bedload parameter
    slp::Float64     = 0.2         # suspended load parameter
    slpot::Float64   = 0.1         # overtopping-driven onshore suspended transport coefficient (iover=1 only)
    bedlm::Float64   = 1.0         # hardbottom reduction mode (1=standard, 0=special)
    cstabn::Float64  = 0.4         # stability constant
    # Derived quantities --------
    # Derived "submerged" quantities. These use s_sub = ρ_s / ρ_w - 1 (NOT sg - 1),
    # so they're consistent with the surrounding `rho_water`. Default values
    # below assume sg=2.65, rho_water=1025 → s_sub ≈ 1.5854. Recompute via
    # `make_sediment` after any change to sg / rho_water.
    wfsgm1::Float64  = 0.026 * 1.5854      # wf * s_sub
    gsgm1::Float64   = 9.81 * 1.5854       # g  * s_sub
    gsd50s::Float64  = 9.81 * 1.5854 * 0.3e-3 * 0.05  # s_sub * g * d50 * shield
    bld::Float64     = 2e-3 / (9.81 * 1.5854)         # blp / gsgm1
    csedia::Float64  = 2.0 * 0.3e-3        # 2·d50
end

"""
    submerged_sgm1(sed::SedimentConfig) -> Float64

Submerged specific-gravity-minus-1, `ρ_s / ρ_w - 1`, used everywhere a
formula needs the underwater grain weight (Shields, settling velocity,
bedload, suspended load). Equivalent to the classic `sg - 1` when
`rho_water == 1000` (freshwater); ~3.9% lower for seawater (1025).
Replaces all hardcoded `(sg - 1)` / `sgm1 = sg - 1` formulations in
sediment and wave code paths.
"""
@inline submerged_sgm1(sed::SedimentConfig) = sed.sg * 1000.0 / sed.rho_water - 1.0

"""
    make_sediment(; d50=0.3e-3, sg=2.65, sporo=0.4, shield=0.05, blp=2e-3,
                    wf=nothing, kwargs...) -> SedimentConfig

Convenience constructor that derives `wfsgm1`, `gsgm1`, `gsd50s`, `bld`,
`csedia` from primary inputs. If `wf` is not given, the Soulsby 1997 fall
velocity is used (see `sediment/fractions.jl :: fall_velocity`).
"""
function make_sediment(; d50::Float64=0.3e-3, sg::Float64=2.65,
                         rho_water::Float64=1025.0,
                         sporo::Float64=0.4, shield::Float64=0.05,
                         blp::Float64=2e-3, tanphi::Float64=0.63,
                         effb::Float64=0.005, efff::Float64=0.005,
                         slp::Float64=0.2, slpot::Float64=0.1, bedlm::Float64=1.0,
                         cstabn::Float64=0.4, wf::Union{Nothing,Float64}=nothing)
    # Submerged specific-gravity minus 1: ρ_s / ρ_w - 1.
    # Replaces classic (sg - 1) for seawater consistency. See `submerged_sgm1`.
    sgm1 = sg * 1000.0 / rho_water - 1.0
    gsgm1 = GRAV * sgm1
    # Derive wf via Soulsby if not given. Avoid circular dep by inlining.
    wf_val = if wf === nothing
        ν = 1.0e-6
        dstar = d50 * (sgm1 * GRAV / ν^2)^(1/3)
        if dstar < 1.0
            sgm1 * GRAV * d50^2 / (18 * ν)
        elseif dstar < 100.0
            ν / d50 * (sqrt(10.36^2 + 1.049 * dstar^3) - 10.36)
        else
            1.05 * sqrt(sgm1 * GRAV * d50)
        end
    else
        wf
    end
    return SedimentConfig(
        wf=wf_val, sg=sg, rho_water=rho_water, sporo1=1.0 - sporo, tanphi=tanphi,
        effb=effb, efff=efff, d50=d50, shield=shield,
        blp=blp, slp=slp, slpot=slpot, bedlm=bedlm, cstabn=cstabn,
        wfsgm1=wf_val * sgm1,
        gsgm1=gsgm1,
        gsd50s=gsgm1 * d50 * shield,
        bld=blp / gsgm1,
        csedia=2.0 * d50,
    )
end

"""
    MultifractionConfig

Defines the multi-layer, multi-fraction bed composition.

When `nf == 1`, the port reduces to single-grain behavior — every transport
formula, the Exner solver, and the grain-sorting code still use the
`(nx, nlayers, nf)` arrays, but with the trivial `nf=1` slice.
"""
Base.@kwdef struct MultifractionConfig
    grain_sizes::Vector{Float64} = [0.3e-3]  # per fraction, meters (not mm!)
    nlayers::Int = 3                          # vertical bed layers
    layer_thickness::Float64 = 0.1            # per-layer thickness (m)
    porosity::Float64 = 0.4                   # bed porosity
    initial_fractions::Vector{Float64} = [1.0]   # sum to 1, length nf
    use_size_dependent_shields::Bool = true      # Soulsby 1997 curve per d_k
    use_hiding_exposure::Bool = false            # off by default to avoid double-counting
    use_grainsize_tadapt::Bool = true            # per-fraction T_adapt
    hiding_method::Symbol = :egiazaroff           # :egiazaroff or :ashida_michiue
    tadapt_multiplier::Float64 = 1.0
    # Top-layer mixing: averages the composition of the top layers above
    # the wave-induced "depth of disturbance" = facDOD × Hrms. The mixed
    # composition is written back to each of the top layers. Setting
    # process_mixtoplayer=false disables the mix entirely. Default
    # facDOD=0.2 (20% of Hrms) provides enough vertical mixing to damp
    # sub-step surface-composition swings while preserving the
    # stratigraphic record in deeper layers.
    process_mixtoplayer::Bool = true
    facDOD::Float64 = 0.2
    # ── Smoothing passes ──
    # These control the strength of spatial smoothing at different stages
    # of the morphodynamic update. Increase to reduce checkerboard noise;
    # decrease for sharper (noisier) gradients. Set to 0 to disable.
    n_face_flux_smooth::Int    = 3    # 1-2-1 passes on face fluxes (exner step 1d)
    n_pickup_smooth::Int       = 10   # 1-2-1 passes on per-fraction pickup (exner step 3b)
    n_composition_smooth::Int  = 3    # 1-2-1 passes on active_frac (before sedtra!)
    # ── Per-class pickup weighting (sets how the per-class potential
    # transport rate is allocated to each grain-size class):
    #     q_k = q_k_pot × w_k
    # with w_k determined by the scheme below. All three schemes preserve
    # the relative AeoLiS-equivalent "surface-mass-fraction × class-physics"
    # structure; they differ only in WHICH bed reservoir defines the
    # surface mass fraction. Default reproduces existing behavior.
    #
    # :active_frac (DEFAULT) — w_k = state.active_frac[j,k], the
    #     composition of the wave-mixed surface layer of thickness
    #     facDOD·Hrms. Most physical for wave-driven mixing; current
    #     CSHORE.jl behavior.
    # :top_layer — w_k = bed_mass[j,1,k] / Σ_j bed_mass[j,1,j], the
    #     top stratigraphic layer alone, ignoring DOD wave mixing.
    #     Closer to "instantaneous bed surface composition." Useful when
    #     wave mixing is weak relative to bed evolution, or when the
    #     bed scheme is the more reliable indicator of availability.
    # :full_bed — w_k = Σ_l bed_mass[j,l,k] / Σ_{l,m} bed_mass[j,l,m],
    #     the full vertical bed inventory at this node. Treats the bed
    #     as a single well-mixed reservoir; appropriate for thin bed
    #     stacks or when stratigraphy is not resolved.
    pickup_weighting::Symbol = :active_frac
    # ── Transport formula selection ──
    # :original       — CSHORE erfc-based formula for all sizes (default)
    # :size_adaptive  — CSHORE for d < 2mm (sand), MPM bedload-only for d ≥ 2mm (gravel)
    # :soulsby_vanrijn — Soulsby-Van Rijn (1997) for all sizes
    transport_formula::Symbol = :original
    # ── SVR transport scale factors ──
    # svr_scale     : overall multiplier on final qb_mag / qs_mag (both components).
    # svr_wave_scale: multiplier on the wave-orbital contribution to Ueff²
    #                 (the 0.018/Cd·Urms² term). Boosts wave-driven entrainment
    #                 independently of the current.
    # svr_current_scale: multiplier on the mean-flow contribution to Ueff²
    #                 (Umag² = Ux²+Vy²). Boosts current-driven entrainment
    #                 independently of the waves.
    # Physical motivation: SVR Asb/Ass coefficients are ~10–50× smaller than
    # calibrated CSHORE effb/efff, so scaling is needed to match field data.
    # Separate wave/current scales allow different calibration for each driver.
    # Defaults = 1.0 (unmodified Soulsby 1997). Typical calibrated range: 5–50.
    svr_scale::Float64         = 1.0
    svr_wave_scale::Float64    = 1.0
    svr_current_scale::Float64 = 1.0
    # ── Spatial grain-size distribution ──
    # When non-empty, overrides `initial_fractions` with per-node values.
    # Shape: (n_bathy_nodes, nf). Each row must sum to 1.0.
    # Use `spatially_varying_fractions()` from presets.jl to generate.
    initial_fractions_spatial::Matrix{Float64} = zeros(0, 0)
end

nfractions(m::MultifractionConfig) = length(m.grain_sizes)

function validate(m::MultifractionConfig)
    nf = nfractions(m)
    length(m.initial_fractions) == nf ||
        throw(ArgumentError("initial_fractions length $(length(m.initial_fractions)) ≠ nf=$nf"))
    sum_f = sum(m.initial_fractions)
    isapprox(sum_f, 1.0; atol=1e-9) ||
        throw(ArgumentError("initial_fractions must sum to 1.0 (got $sum_f)"))
    all(d -> d > 0, m.grain_sizes) ||
        throw(ArgumentError("grain_sizes must be positive"))
    m.nlayers >= 1 || throw(ArgumentError("nlayers must be ≥ 1"))
    0 ≤ m.porosity < 1 || throw(ArgumentError("porosity must be in [0,1)"))
    m.pickup_weighting in (:active_frac, :top_layer, :full_bed) ||
        throw(ArgumentError("pickup_weighting must be :active_frac, :top_layer, or :full_bed " *
                            "(got :$(m.pickup_weighting))"))
    return nothing
end

"""
    BoundaryTimeSeries

Time-ordered arrays of offshore boundary conditions. The model interpolates
these at each model time via TSINTP.
"""
Base.@kwdef struct BoundaryTimeSeries
    # Wave boundary conditions
    timebc::Vector{Float64}       # time values (s)
    tpbc::Vector{Float64}         # period (s)
    hrmsbc::Vector{Float64}       # offshore Hrms (m)
    wsetbc::Vector{Float64}       # offshore setup
    swlbc::Vector{Float64}        # still water level (m)
    wangbc::Vector{Float64}       # wave angle (rad or deg — FORTRAN uses deg in input, rad internal)
    # Infragravity (IG) wave offshore boundary condition (optional).
    # When non-empty and IgConfig is active, hrms_ig_bc[itime] is used as
    # the seaward IG Hrms rather than the internally-derived kappa_ig · hrms0.
    # Length must equal length(timebc) if provided.
    hrms_ig_bc::Vector{Float64} = Float64[]
    # Landward water-body still-water level (IWTRAN=1). When non-empty, used
    # for the landward wave march and back-side h computation. When empty,
    # the landward SWL falls back to swlbc (seaward water body extends behind).
    # Length must equal length(timebc) if provided.
    swl_landward::Vector{Float64} = Float64[]
    # Wind forcing
    w10::Vector{Float64}      = Float64[]   # 10 m wind speed
    wangle::Vector{Float64}   = Float64[]   # wind direction
    windcd::Vector{Float64}   = Float64[]   # drag coefficient
end

ntime(bc::BoundaryTimeSeries) = length(bc.timebc)

"""
    BathyInput

Raw (unsampled) input bathymetry for each of `iline` cross-shore lines. The
model later interpolates onto a uniform grid via SUBROUTINE BOTTOM.
"""
Base.@kwdef struct BathyInput
    xbinp::Matrix{Float64}        # (max_nbinp, iline)
    zbinp::Matrix{Float64}
    fbinp::Matrix{Float64}        # bottom friction factor
    nbinp::Vector{Int}            # (iline,) — number of input points per line
    xs::Vector{Float64}           # (iline,) x-location of line
    yline::Vector{Float64} = Float64[]
    dyline::Vector{Float64} = Float64[]
    agline::Vector{Float64} = Float64[]   # angle of each line (deg)
    # Hardbottom elevation (ISEDAV≥1). Empty matrix means "no hardbottom" — the
    # initializer fills state.zb_hard with the -Inf sentinel. Otherwise this
    # must be `(max_nbinp, iline)` in the same layout as `zbinp`.
    zhinp::Matrix{Float64} = zeros(0, 0)
    # Manning's n per (input-row, iline). Empty means "not using Manning";
    # bottom friction stays as fbinp. Consumed when
    # `options.ifriction_spatial == 2` to recompute fb2 dynamically each step
    # from the live water depth (fb2 = g·n²/h^(1/3)).
    manning_n::Matrix{Float64} = zeros(0, 0)
end

"""
    VegetatedManningField

A per-node Manning's n field whose values already encode a vegetation
category's frictional effect (via `nbs_vegetation_manning_field`).
`build_config` refuses to combine this with a `VegetationInput` to prevent
double-counting (vegetation drag through bed friction AND through Cd-based
stem drag).

Plain `Vector{Float64}` Manning fields are assumed to be bed-only and may
freely combine with `VegetationInput`.
"""
struct VegetatedManningField
    values::Vector{Float64}
    categories::Vector{Tuple{Symbol,Symbol}}  # (category, region) entries composed in
end

Base.length(f::VegetatedManningField) = length(f.values)
Base.getindex(f::VegetatedManningField, i) = f.values[i]

"""
    PorousInput

Present only when `options.iperm == 1`.
"""
Base.@kwdef struct PorousInput
    xpinp::Matrix{Float64}
    zpinp::Matrix{Float64}
    npinp::Vector{Int}
    wnu::Float64    = 1.0e-6      # kinematic viscosity (m²/s)
    snp::Float64    = 0.4         # porosity of permeable layer (gravel/stone)
    sdp::Float64    = 0.02        # nominal stone diameter D_n50 (m) — 20mm default for gravel
    alpha::Float64  = 1000.0      # laminar resistance coefficient α₀ (Eq. 68)
    beta1::Float64  = 5.0         # turbulent resistance coefficient β₀ (Eq. 69)
    beta2::Float64  = 5.0         # (unused — computed internally from β₀, snp, sdp)
end

"""
    PorousInput(x, zp; porosity=0.4, stone_diameter=0.02)

Convenience constructor. `x` and `zp` are 1D vectors defining the porous
layer bottom elevation along the profile. Wraps them into the matrix
format required internally.

# Example
```julia
por = PorousInput(x, z_floor; porosity=0.4, stone_diameter=0.05)
```
"""
function PorousInput(x::Vector{Float64}, zp::Vector{Float64};
                     porosity::Float64=0.4, stone_diameter::Float64=0.02,
                     viscosity::Float64=1.0e-6)
    np = length(x)
    length(zp) == np || throw(DimensionMismatch("x and zp must have the same length"))
    PorousInput(
        xpinp = reshape(copy(x), np, 1),
        zpinp = reshape(copy(zp), np, 1),
        npinp = [np],
        snp   = porosity,
        sdp   = stone_diameter,
        wnu   = viscosity,
    )
end

"""
    DikeErosionInput

Grassed-dike erosion parameters (IPROFL=2). `gdinp`, `grinp`, `grdinp` are
spatial input arrays (np_raw, nl) that get interpolated onto the uniform
grid by `build_config`. After initialization these become `grsd`, `grsr`,
`grsrd` per-node state fields.

- `deeb`: wave-breaking dissipation efficiency → dike erosion (default 0.005)
- `deef`: bed-friction dissipation efficiency → dike erosion (default 0.005)
- `gdinp`: grass layer thickness profile (m)
- `grinp`: undamaged grass resistance (J/m²)
- `grdinp`: damaged grass resistance at layer base (J/m²)
"""
Base.@kwdef struct DikeErosionInput
    deeb::Float64 = 0.005
    deef::Float64 = 0.005
    gdinp::Matrix{Float64}  = zeros(0, 0)
    grinp::Matrix{Float64}  = zeros(0, 0)
    grdinp::Matrix{Float64} = zeros(0, 0)
end

"""
    ClayInput

Sand-over-clay erosion parameters (ICLAY=1). When the sand layer is
thinner than `d50` at a given node, the underlying clay is exposed and
erodes at rate `rclay · dsta`, where `rclay` is derived from the input
erosion rate `rcinp` as `rclay = g / rcinp`.

- `deeb`: wave-breaking dissipation efficiency → clay erosion (default 0.005)
- `deef`: bed-friction dissipation efficiency → clay erosion (default 0.005)
- `rcinp`: clay erosion-rate input (W·s/m³ equivalent); raw input (np_raw, nl)
- `fcinp`: volumetric sand fraction in clay; raw input (np_raw, nl)
  (FCLAY = 1 - FCINP/(1-n) applied inside build_config to give the fraction
   of sand liberated per unit clay erosion; must satisfy 0 ≤ fcinp ≤ 1-n)
"""
Base.@kwdef struct ClayInput
    deeb::Float64 = 0.005
    deef::Float64 = 0.005
    rcinp::Matrix{Float64} = zeros(0, 0)
    fcinp::Matrix{Float64} = zeros(0, 0)
end

"""
    VegetationInput

Present only when `options.iveg > 0`.

Fields:
- `vegcd`: drag coefficient (dimensionless, ~1.0 for rigid stems)
- `vegcdm`: inertia coefficient (for IVEG=3 momentum; ratio VEGCDM/VEGCD)
- `vegn`: stem density (stems/m²) — spatially varying (nn, nl)
- `vegb`: blade/stem width (m) — spatially varying (nn, nl)
- `vegd`: vegetation height above bed (m) — spatially varying (nn, nl)
- `vegh`: effective canopy height (m) = min(vegd, water depth) — used for IVEG=1,2 friction
- `vegrd`: root depth (m) — for uprooting check
- `vegrh`: root height (m) — for uprooting check
- `vegfb`: precomputed friction modifier = vegcd * vegn * vegb / fb2 (computed at init)
"""
Base.@kwdef struct VegetationInput
    vegcd::Float64 = 1.0
    vegcdm::Float64 = 1.0
    vegn::Matrix{Float64}   = zeros(0, 0)   # density (stems/m²)
    vegb::Matrix{Float64}   = zeros(0, 0)   # blade width (m)
    vegd::Matrix{Float64}   = zeros(0, 0)   # vegetation height above bed (m)
    vegh::Matrix{Float64}   = zeros(0, 0)   # effective canopy height (m) — same as vegd at input
    vegrd::Matrix{Float64}  = zeros(0, 0)   # root depth (m)
    vegrh::Matrix{Float64}  = zeros(0, 0)   # root height (m)
    vegfb::Matrix{Float64}  = zeros(0, 0)   # friction modifier (precomputed)
end

"""
    SwashConfig

Configuration for the BDJ swash-zone model (`wetdry!`) and overtopping flux.

These parameters control the time-averaged swash hydrodynamics that run
regardless of whether `iover=0` or `iover=1`.  With `iover=0` the overtopping
flux `qo` is held fixed at zero (single-pass, no wave-solver feedback); with
`iover=1` the solver iterates until `qo` converges.

Fields
------
- `awd`: BDJ swash amplitude parameter (default 1.6 from empirical fit)
- `ewd`, `cwd`, `aqwd`, `bwd`, `agwd`, `auwd`: derived BDJ parameters
  (computed internally in `wetdry!` when set to zero)
- `wpm`: maximum seepage velocity for infiltration (m/s); computed from
  sediment properties when `infilt=1`
- `rcrest`: optional crest elevation override per line (m); auto-detected
  from the profile when empty
- `rwh`, `diketoe`, `runup_kappa`, `runup_phi`: dike/runup parameters
"""
Base.@kwdef struct SwashConfig
    rwh::Float64    = 0.0         # runup wire height
    rcrest::Vector{Float64} = Float64[]   # crest elevation per line
    diketoe::Float64 = 0.0
    runup_kappa::Float64 = 1.0
    runup_phi::Float64   = 1.0
    # Swash parameters
    awd::Float64 = 0.0; wdn::Float64 = 0.0; ewd::Float64 = 0.0; cwd::Float64 = 0.0
    aqwd::Float64 = 0.0; bwd::Float64 = 0.0; agwd::Float64 = 0.0; auwd::Float64 = 0.0
    wpm::Float64 = 0.0
end

# Backward-compatibility alias — existing code using OvertoppingConfig continues to work.
const OvertoppingConfig = SwashConfig

"""
    GroundwaterConfig

Configuration for the coupled beach groundwater and surface-moisture model.

Implements the 1D nonlinear Boussinesq equation for water-table dynamics
in the beach aquifer, forced by CSHORE's native wave-setup and BDJ swash
fields (replacing the Stockdon parameterisation used in Psamathe).
The Van Genuchten retention curve converts the water-table depth to
volumetric surface moisture, which feeds the aeolian transport threshold.

# Groundwater physics
- `K`   : hydraulic conductivity (m/s); ≈ 1×10⁻³ for medium beach sand,
          ≈ 5×10⁻³ for coarse sand, ≈ 1×10⁻⁴ for fine sand.
- `ne`  : effective (drainable) porosity (-); typically 0.20–0.35.
- `D`   : aquifer thickness below datum (m); distance from MSL to impermeable base.

# Van Genuchten moisture retention (Carsel & Parrish 1988 medium-sand defaults)
- `vg_alpha`  : α (1/m); inverse of air-entry capillary head.
- `vg_n`      : n (-); controls curve steepness; m = 1 − 1/n computed internally.
- `theta_res` : residual volumetric moisture (-); θ as depth-to-table → ∞.
- `theta_sat` : saturated volumetric moisture (-); θ when table at surface.

# Infiltration source coupling
- `infiltration_rate` : when `nothing` (default) and `infilt=1`, the rate is
  automatically set to `wpm` (from `wetdry!`'s Forchheimer calculation) so the
  swash momentum drain and aquifer recharge use the same flux — mass is conserved.
  When `nothing` and `infilt=0`, there is no swash source term (seaward BC only).
  Set explicitly only to override both (raises a validation warning with `infilt=1`).

# Rainfall recharge
- `rainfall_rate` : rainfall intensity (m/s) applied over the dry backshore
  (landward of the swash zone).  Can be:
    - `nothing` (default) — no rainfall.
    - `Float64` — constant rate for all BC windows.
    - `Vector{Float64}` of length `length(boundary.timebc)` — time-varying rate
      that is linearly interpolated at each BC window, e.g. to simulate a storm
      hyetograph, diurnal convective pattern, or seasonal rainfall cycle.
  Positive values add recharge to the Boussinesq equation at every dry-backshore
  node (`j > jdry`).  When active, surface moisture is forced to `theta_sat`
  (or `rainfall_wet_theta`) at rain-receiving nodes, overriding Van Genuchten —
  rainfall suppresses aeolian transport regardless of water-table depth.
  Typical values: light rain ≈ 5e-7 m/s (2 mm/h); heavy storm ≈ 5e-6 m/s
  (18 mm/h); extreme (1-in-10 yr) ≈ 2e-5 m/s (72 mm/h).
- `rainfall_wet_theta` : surface moisture during rain events.  Scalar or
  time-series (same length rules as `rainfall_rate`).  Defaults to `theta_sat`.
  Set lower (e.g. 0.25) to represent partial wetting from light rain, or vary
  in time to simulate drying between events.

# Boundary conditions
- `gw_eta_landward` : fixed water-table elevation (m) at the landward profile
  boundary. `NaN` (default) applies a no-flux (closed) Neumann condition.

# Solver
- `dt_gw` : internal Boussinesq sub-step (s). Must satisfy the CFL stability
  limit `dt ≤ nₑ dx² / (2 K D)`.  For typical beach parameters (K=10⁻³ m/s,
  dx=1 m, D=3 m, nₑ=0.30) the CFL limit is ≈ 50 s; default 30 s is safe.
"""
Base.@kwdef struct GroundwaterConfig
    # Aquifer hydraulic parameters
    K::Float64          = 1e-3     # hydraulic conductivity (m/s)
    ne::Float64         = 0.30     # effective (drainable) porosity (-)
    D::Float64          = 3.0      # aquifer thickness below datum (m)
    # Van Genuchten retention curve — medium beach sand defaults
    vg_alpha::Float64   = 14.5     # α (1/m)
    vg_n::Float64       = 2.68     # n (-);  m = 1 − 1/n computed in groundwater.jl
    theta_res::Float64  = 0.045    # residual volumetric moisture (-)
    theta_sat::Float64  = 0.38     # saturated volumetric moisture (-)
    # Swash infiltration source coupling (see docstring)
    infiltration_rate::Union{Nothing,Float64} = nothing
    # Rainfall recharge (see docstring)
    # Scalar: constant rate for all BC windows.
    # Vector: time series of length == length(boundary.timebc); interpolated each BC window.
    rainfall_rate::Union{Nothing,Float64,Vector{Float64}} = nothing
    # Surface θ during rain: scalar constant or time series (same length rules as rainfall_rate).
    rainfall_wet_theta::Union{Nothing,Float64,Vector{Float64}} = nothing
    # Landward boundary: NaN → no-flux; otherwise fixed WT elevation (m)
    gw_eta_landward::Float64 = NaN
    # Solver
    dt_gw::Float64      = 30.0     # internal Boussinesq sub-step (s)
end

"""
    IgConfig

Configuration for infragravity (IG) wave energy in the time-averaged
framework.  Infragravity waves (periods ≈ 25–300 s) are generated by
nonlinear interactions between short-wave groups and dominate swash
dynamics on dissipative, low-sloping beaches (Iribarren ξ < 0.3).

Because CSHORE averages over both short-wave (~8-12 s) and IG (~25-300 s)
periods (BC windows are hourly), IG energy cannot appear as a resolved
oscillation.  Instead it enters as a second bulk energy component that
propagates alongside Hrms_ss, contributes to bottom orbital velocity, and
— when `ig_swash_active = true` — augments the swash flux at the runup front.

# Layers of physics

**Layer 1 — IG height profile (always computed when `IgConfig` is provided):**

The offshore IG height is estimated as

    hrms_ig_0 = kappa_ig · hrms0

where `kappa_ig` is a dimensionless coupling coefficient (~0.3–0.5 for
open-ocean swell, ~0.1–0.2 for wind-sea dominated conditions, and ~0.5–0.8
on very dissipative beaches).  If `hrms_ig_bc` is provided in
`BoundaryTimeSeries` it overrides this estimate.

IG height shoals landward via depth-averaged Green's law (h^{-1/4})
until it reaches either the saturation cap (γ_ig · h) or the dry boundary.
The IG wave period is `Tig = Tig_factor · Tp` (default 3.0).

**Layer 2 — IG energy balance (when `ig_energy_balance = true`):**

Propagates a 1D IG energy budget from offshore to shore:

    d(F_ig)/dx = kappa_source · |dF_ss/dx| − f_bf_ig · E_ig · sqrt(g/h)

where the source term tracks short-wave breaking as the IG generator and
the dissipation term accounts for IG bottom friction.

**Layer 3 — IG swash augmentation (when `ig_swash_active = true`):**

Adds `c_ig_swash · hrms_ig` to the swash-front water depth used in the
wet/dry flux calculation, increasing the effective swash excursion on
dissipative profiles where IG dominates runup.

# Transport coupling

`ustd_ig[j]` (shallow-water orbital velocity ≈ hrms_ig · √(g/h) / 2) is
added in quadrature to the short-wave `ustd[j]` when computing transport
probabilities and bedload magnitudes in `sedtra!`:

    ustd_combined = sqrt(ustd_ss² + ustd_ig²)

This increases IG-driven sediment mobility in the inner surf and swash zones
without affecting the wave setup, rollers, or radiation stress — only the
orbital velocity in the transport formula is augmented.

# Fields

- `kappa_ig`         : offshore IG / SS height ratio (seaward BC; default 0.35).
- `gamma_ig`         : IG saturation limit H_ig/h (default 0.3).
- `Tig_factor`       : T_IG = Tig_factor × Tp (default 3.0).
- `ig_energy_balance`: enable Layer-2 energy-balance propagation (default false).
- `f_bf_ig`          : IG bottom-friction coefficient for Layer-2 dissipation
                       (default 0.01 m²/s³ — weak friction for long waves).
- `kappa_source`     : coupling coefficient from short-wave breaking to IG
                       source term in Layer 2 (default 0.3).
- `ig_swash_active`  : enable Layer-3 swash augmentation (default false).
- `c_ig_swash`       : Layer-3 IG contribution to swash depth (default 0.5).
- `ustd_ig_in_transport`: couple `ustd_ig` into transport probabilities / bedload
                       (default true when IgConfig is provided; set false for
                       diagnostic-only runs where you want the IG field computed
                       and output but NOT coupled to sediment transport).
"""
Base.@kwdef struct IgConfig
    kappa_ig::Float64          = 0.35   # offshore IG/SS height ratio
    gamma_ig::Float64          = 0.30   # IG saturation cap H_ig/h
    Tig_factor::Float64        = 3.0    # T_IG = Tig_factor * Tp
    ig_energy_balance::Bool    = false  # Layer 2: 1D energy balance
    f_bf_ig::Float64           = 0.01   # IG bottom-friction coeff (Layer 2)
    kappa_source::Float64      = 0.30   # SW-breaking → IG source coupling (Layer 2)
    ig_swash_active::Bool      = false  # Layer 3: swash augmentation
    c_ig_swash::Float64        = 0.5    # IG contribution to swash front depth
    ustd_ig_in_transport::Bool = true   # couple ustd_ig into sedtra!
end

#=============================================================================
UndertowConfig — vertical-structure correction for the cross-shore mean
                 (return) current used in sediment transport.

Background
──────────
state.umean[j] is the depth-averaged cross-shore mean current obtained from
the radiation-stress balance (and the roller mass-flux multiplier). In the
surf zone, however, the offshore return flow is concentrated near the bed —
the depth-averaged value typically underestimates the *near-bed* undertow
that actually drives suspended-load transport at the bar crest.

References: Reniers et al. (2004 JGR-Oceans 109 C01030),
            Hoefel & Elgar (2003 Science 299:1885),
            Garcez Faria et al. (2000 JGR-Oceans 105 C7).

Two amplification modes
───────────────────────
:hrms_h        U_bed = U_mean · (1 + α · (Hrms/h)^p)
               • α  = `alpha`     (typical 0.4 – 1.0)
               • p  = `exponent`  (typical 2 – 4)
               • Hrms/h grows from ≈ 0.1 offshore to ≈ γ_b at break point,
                 so amplification grows naturally with breaking strength.

:dissipation   U_bed = U_mean - β · D_b / (ρ · g · h · c_p)
               • β  = `alpha`   (typical 1.0 – 2.0)
               • D_b is the breaker-dissipation (state.dbsta in m³/s²)
               • c_p is local phase celerity
               • Density ρ is implicit (canonical CSHORE units already
                 normalize by ρg; see _undertow_bed in transport.jl).

The bed velocity is used **only** inside sedtra! for the suspended-load
direction/magnitude — the depth-averaged umean is unchanged everywhere
else (longshore momentum, IG forcing, runup, etc.).
=============================================================================#
Base.@kwdef struct UndertowConfig
    enabled::Bool       = false
    mode::Symbol        = :hrms_h     # :hrms_h | :dissipation
    alpha::Float64      = 0.7         # amplification coefficient (or β for :dissipation)
    exponent::Float64   = 3.0         # power on Hrms/h  (only :hrms_h)
    h_min::Float64      = 0.05        # lower clamp on h to avoid singularity
    cap::Float64        = 3.0         # max amplification factor (1 + α·… ≤ cap)
end

#=============================================================================
AsymmetryConfig — explicit physically-grounded wave-shape closure that
                  replaces the empirical Ruessink Ur→As mapping near
                  the breakpoint.

Background
──────────
Ruessink et al. (2012) ties Sk and As to the Ursell number alone — works
well in shoaling but produces non-zero As far offshore of breaking and
misses the sharp acceleration-skewness signature at the breakpoint.

The Doering & Bowen (1995) / Reniers et al. (2004) closure derives the
asymmetry from the *gradient of breaker dissipation*:

    As(j) ≈ -K_as · (Db[j+1] - Db[j-1]) / (2·dx)
                       · h[j] / (ρ·g·max(Hrms²,ε)·ω)

The dissipation gradient naturally:
    • vanishes outside the surf zone        → As ≈ 0
    • is negative just inside breaking      → As < 0  (front-leaning)
    • can flip positive across the bar crest → As > 0  (back-leaning)

`K_as` is dimensionless, with literature values 0.2 – 0.6 (Walstra+ 2012).
Sk is left to the chosen Sk method (Stokes2 or Ruessink).

Field selection (`method`):
    :ruessink         — current default (Ur → Sk & As, Ruessink 2012)
    :stokes2          — Stokes 2nd-order Sk; As = 0
    :boussinesq_diss  — Sk from Stokes2; As from breaker-dissipation gradient

The legacy `OptionFlags.iasym` integer is retained for backwards
compatibility; when an `AsymmetryConfig` is supplied it takes precedence.
=============================================================================#
Base.@kwdef struct AsymmetryConfig
    method::Symbol           = :ruessink     # :ruessink | :stokes2 | :boussinesq_diss
    K_as::Float64            = 0.4           # gradient coefficient (only :boussinesq_diss)
    smooth_diss_window::Int  = 3             # boxcar window over Db before differencing
    diss_grad_cap::Float64   = 0.8           # |As| cap from this term alone
end

#=============================================================================
PhaseLagConfig — non-equilibrium suspended sediment with cross-shore lag.

Background
──────────
The default CSHORE transport closure assumes *local equilibrium*: the
suspended sediment volume v_s at node j is determined entirely by the
local stirring rate (ε_b·D_b + ε_f·D_f) divided by w_s.  In reality,
sediment stirred at the break point requires time τ ≈ h/w_s (~10–60 s
in the surf zone) to settle — and during that time it is advected by
the mean cross-shore current (undertow), forming a "suspended cloud"
that deposits some distance offshore of where it was picked up.

This lag is the primary mechanism by which subtidal bars form and
migrate (Reniers et al. 2004 JGR-Oceans 109:C01030; Houwman & Ruessink
1996; Camenen & Larson 2007/2008).  Equilibrium closures cannot
reproduce this cloud advection — they predict deposition wherever
stirring is reduced — and therefore tend to flatten bar–trough
morphology.

Implementation
──────────────
After each kernel populates state.vs[j,k] = v_s,eq(j,k) (local-
equilibrium suspended volume), apply a 1-D upwind relaxation along the
profile:

    U(j)·∂v_s/∂x = (v_s,eq − v_s) / τ_lag

with τ_lag the user-supplied settling timescale.  In practice we treat
this as `α(j) = dx / max(|U|·τ_lag, L_min)` and march upwind so that
the lagged volume on each cell relaxes toward the local equilibrium
value over an effective length scale L = |U|·τ_lag (typically 2–20 m
in the surf zone for fine-medium sand).

Once v_s,lag is populated, qsx = aslope · U_bed · v_s,lag · comp_k
replaces the equilibrium qsx everywhere the relaxation is active.

Parameters
──────────
  enabled : master toggle
  tau_lag : Eulerian settling timescale (s).  Set explicitly if you
            want a constant lag, or leave NaN to derive per-cell from
            h(j)/w_s,k(k) automatically.
  L_min   : floor on |U|·τ_lag to avoid α > 1 when undertow is weak (m).
            With L_min ≤ dx, α saturates at 1 → recovers equilibrium.
  cap_alpha : max relaxation per cell (≤ 1; 1.0 = full equilibrium).

Notes
─────
  • A single backward sweep handles the undertow-dominated surf-zone
    case (U_mean < 0 throughout the breaking region).  A second forward
    sweep handles any U_mean > 0 cells (onshore-flowing currents,
    swash uprush) so the routine is sign-aware.
  • Bedload is unchanged — bedload responds in fractions of a wave
    period and is well-described by a local equilibrium.
=============================================================================#
Base.@kwdef struct PhaseLagConfig
    enabled::Bool       = false
    tau_lag::Float64    = NaN     # settling time [s]; NaN → derive h/w_s per cell
    L_min::Float64      = 0.5     # min lag length [m] — floors α below 1
    cap_alpha::Float64  = 1.0     # max relaxation per cell (1.0 = equilibrium)
end

#=============================================================================
BailardConfig — Bailard (1981) velocity-moment transport correction.

Background
──────────
The default CSHORE suspended-load closure approximates the time-averaged
moment ⟨|u|³·u⟩ as `aslope · U · v_s`, where v_s is a stirring-rate proxy.
This collapses two distinct physical mechanisms into a single term:

  ⟨|u|³·u⟩  ≈  K_xc · u_rms³ · U      (wave-current cross term)
              + K_sk · u_rms⁴ · Sk    (skewness term)
              + K_as · u_rms⁴ · As    (asymmetry term)
              + |U|³·U                (current-only — small in surf zone)

The cross term `u_rms³·U` and the skewness term `u_rms⁴·Sk` have *opposite
signs* across the breakpoint:

  • Shoaling (U ≈ 0, Sk > 0)  → onshore transport (q > 0)
  • Surf zone (U < 0, Sk ≈ 0.4) → offshore transport (q < 0)
  • Crossover at the breakpoint → sharp ∂q/∂x

This is the spatial gradient signature absent from the equilibrium kernel
and identified as missing in the Tier-1 / Tier-2 residual analysis. Adding
these terms additively over the existing qsx tests whether the missing
physics is the velocity-moment representation or something else.

References: Bailard (1981) JGR-Oceans 86 C11; Stive (1986); Roelvink & Stive
(1989); Hoefel & Elgar (2003).

Parameters
──────────
  enabled  : master toggle
  eps_s    : suspended-load Bailard efficiency (0.01 – 0.04 typical)
  gamma_xc : multiplier on the wave-current cross term `u_rms³ · U`
  gamma_sk : multiplier on the skewness term       `u_rms⁴ · Sk`
  gamma_as : multiplier on the asymmetry term      `u_rms⁴ · As`
  additive : true → add to existing qsx (perturbation test)
             false → REPLACE qsx (full Bailard formulation; needs re-cal)

The sign convention follows CSHORE: +x = onshore, U_bed < 0 = undertow,
Sk > 0 = onshore-biased (crest > trough), As < 0 = front-leaning.

Notes
─────
  • Per-fraction: q_bailard ÷ ((s−1)·g·w_s,k) gives proper m²/s units.
  • Bedload Bailard term ⟨|u|²·u⟩ is *not* added here — bedload responds
    in fractions of a wave period and is well-captured by the existing
    pb·gslope·ustd³ formulation. Only the suspended-load cross + skewness
    terms are added.
=============================================================================#
Base.@kwdef struct BailardConfig
    enabled::Bool      = false
    eps_s::Float64     = 0.025    # Bailard suspended-load efficiency
    gamma_xc::Float64  = 1.0      # wave-current cross term weight
    gamma_sk::Float64  = 1.0      # skewness term weight
    gamma_as::Float64  = 0.0      # asymmetry term weight (default off)
    additive::Bool     = true     # add to existing qsx (true) or replace (false)
end

#=============================================================================
WaveNonlinearityConfig — unified configuration for every wave-nonlinearity
mechanism in the model. Owns the closure family (Sk/As), coupling strengths,
spatial weighting, crest-height correction, biphase relaxation, phase-lag
on suspended sediment, and Bailard velocity-moment cross/skew terms.

Replaces the older scattered surfaces:
  - `OptionFlags.iasym` / `OptionFlags.iskew_spatial`
  - top-level `facSK`, `facAS`, `ur_sk_ref`, `alpha_sk`, `biphase_relax_L`
  - the standalone `AsymmetryConfig`, `PhaseLagConfig`, `BailardConfig`
    structs (kept as deprecated aliases — emit a warning when passed to
    `build_config`; their fields are mirrored as flat fields below).

Recommended usage:

    cfg = build_config(...; wave_nonlinearity = WaveNonlinearityConfig(
        enabled = true,
        closure = :ruessink,
        skewness = 1.0,
        asymmetry = 0.0,
        spatial_weighting = false,
    ))

For Bailard velocity-moment cross/skew terms, set `bailard_enabled = true`
and pick `bailard_eps_s`, `bailard_gamma_xc`, etc.  For phase-lag on
suspended sediment, set `phase_lag_enabled = true`.

All fields default to "wave-nonlinearity disabled" — equivalent to the
pre-commit behavior with `iasym=0`.
=============================================================================#
Base.@kwdef struct WaveNonlinearityConfig
    # ── Master enable + closure family ──────────────────────────────────────
    enabled::Bool      = false               # master toggle for the whole module
    # Closure family for skewness Sk and asymmetry As:
    #   :linear            — no nonlinearity (Sk=As=0; same as enabled=false)
    #   :ruessink          — Ruessink et al. 2012 (Ursell-based)
    #   :stokes2           — Stokes 2nd-order
    #   :boussinesq_diss   — Stokes2 Sk + dissipation-gradient As
    closure::Symbol    = :ruessink

    # ── Coupling strengths ──────────────────────────────────────────────────
    skewness::Float64  = 1.0                 # was facSK / Bailard gamma_sk
    asymmetry::Float64 = 0.0                 # was facAS / Bailard gamma_as

    # ── Spatial weighting (Ursell-based shoaling-zone taper) ────────────────
    spatial_weighting::Bool = false          # was OptionFlags.iskew_spatial > 0
    ur_reference::Float64   = 0.20           # was ur_sk_ref

    # ── Wave-form corrections ───────────────────────────────────────────────
    alpha_crest::Float64           = 0.0     # was alpha_sk (crest-height taper)
    biphase_relax_length::Float64  = 0.0     # was biphase_relax_L (0 = static)

    # ── :boussinesq_diss closure parameters (formerly AsymmetryConfig) ──────
    K_as::Float64           = 0.4            # gradient coefficient
    smooth_diss_window::Int = 3              # boxcar window over Db
    diss_grad_cap::Float64  = 0.8            # |As| cap from gradient term

    # ── Phase-lag on suspended sediment (formerly PhaseLagConfig) ───────────
    phase_lag_enabled::Bool   = false
    phase_lag_tau::Float64    = NaN          # NaN → derive from h/w_s per cell
    phase_lag_L_min::Float64  = 0.5
    phase_lag_cap::Float64    = 1.0

    # ── Bailard velocity-moment closure (formerly BailardConfig) ────────────
    bailard_enabled::Bool     = false
    bailard_eps_s::Float64    = 0.025
    bailard_gamma_xc::Float64 = 1.0          # wave-current cross term weight
    bailard_additive::Bool    = true         # add to qsx (true) or replace
end

# Convenient constructor that accepts a single preset symbol.
#   WaveNonlinearityConfig(:off)       — equivalent to enabled=false
#   WaveNonlinearityConfig(:ruessink)  — typical default-on configuration
#   WaveNonlinearityConfig(:bailard)   — Bailard cross + skewness terms
function WaveNonlinearityConfig(preset::Symbol; kwargs...)
    if preset === :off || preset === :linear
        return WaveNonlinearityConfig(; enabled = false, closure = :linear, kwargs...)
    elseif preset === :ruessink
        return WaveNonlinearityConfig(; enabled = true,  closure = :ruessink, kwargs...)
    elseif preset === :ruessink_asym
        return WaveNonlinearityConfig(; enabled = true,  closure = :ruessink,
                                        asymmetry = 1.0, kwargs...)
    elseif preset === :stokes2
        return WaveNonlinearityConfig(; enabled = true,  closure = :stokes2, kwargs...)
    elseif preset === :bailard
        return WaveNonlinearityConfig(; enabled = true,  closure = :ruessink,
                                        bailard_enabled = true, kwargs...)
    else
        throw(ArgumentError("Unknown WaveNonlinearityConfig preset: $preset. " *
                            "Try :off, :ruessink, :ruessink_asym, :stokes2, :bailard."))
    end
end

"""
    ThermalConfig

Parameters for the minimal 1D thermal permafrost model. Passed via
`build_config(..., thermal=ThermalConfig(...))`. Defaults approximate a
high-latitude organic-rich silt with ~35% moisture.
"""
Base.@kwdef struct ThermalConfig
    nz::Int           = 30          # number of vertical cells per column
    dz::Float64       = 0.1         # vertical cell thickness (m)
    n_rep::Int        = 0           # DEPRECATED: thermal columns are now tracked per shore node; kept for backwards-compat with old scripts
    k_frozen::Float64 = 1.5         # W/m/K
    k_thawed::Float64 = 0.8         # W/m/K
    C_frozen::Float64 = 1.8e6       # J/m³/K
    C_thawed::Float64 = 3.0e6       # J/m³/K
    L::Float64        = 3.34e8      # J/m³ latent heat of fusion
    moisture::Float64 = 0.35        # volumetric soil moisture (0–1)
    T_init::Float64   = -5.0        # initial column temperature (°C) — used as uniform fallback
    T_lower::Float64  = -5.0        # fixed bottom-of-column temperature (°C)
    alt_min::Float64  = 0.0         # minimum reported ALT (m)
    alt_max::Float64  = 10.0        # maximum reported ALT (m)
    cfl_safety::Float64 = 0.4       # CFL safety factor for heat-equation sub-stepping
    # Absolute elevation of the thermal column bottom (m). When finite,
    # each node's column depth = zb[j] - z_bottom, giving variable nz per
    # node so tall bluffs have deep columns while shallow nodes have short
    # ones. The deep BC T_lower is applied at this elevation everywhere.
    # Default -Inf = legacy mode: every column is nz × dz below the local
    # bed surface (backward compatible).
    z_bottom::Float64 = -Inf
    # Optional depth-dependent initial temperature profile.
    # Vector of (depth_m, T_degC) tuples ordered by increasing depth.
    # When provided, each thermal column above the SWL is initialized by
    # linearly interpolating this profile onto the column's cell centers.
    # Submerged columns still get T_init (uniform).  If empty, T_init is
    # used everywhere (backward-compatible default).
    T_init_profile::Vector{Tuple{Float64,Float64}} = Tuple{Float64,Float64}[]
    # --- Thermal intervention fields (optional, per-node) ---
    # Surface thermal resistance from sod/moss insulation (m²·K/W).
    # Length must equal the number of cross-shore nodes. Empty = no insulation.
    # Populated by nbs_sod_insulation() preset.
    R_insulation::Vector{Float64} = Float64[]
    # Thermosyphon heat extraction rate (W/m²), Gaussian-blended from pipe
    # locations. Empty = no thermosyphons. Populated by nbs_thermosyphon().
    Q_thermosyphon::Vector{Float64} = Float64[]
end

"""
    DiffusionConfig

Hillslope diffusion parameters for non-wave gravity-driven mass wasting.

The diffusivity D increases nonlinearly with slope:
  D = D_base × [1 + (tan(S)/tan(S_c))^n × f_c]

When `thermal_control=true`, diffusion is blocked where the active layer
thickness (ALT) is below `thaw_threshold` — frozen permafrost acts as a
rigid floor. Face fluxes use min(D_left, D_right) so frozen cells also
block flux from adjacent thawed cells.
"""
Base.@kwdef struct DiffusionConfig
    D_base::Float64          = 0.01     # base diffusivity for gentle slopes (m²/day)
    critical_slope::Float64  = 0.7      # critical slope angle (rad, ~40°)
    slope_exponent::Float64  = 2.0      # nonlinearity: 1=linear, 2=quadratic
    critical_factor::Float64 = 10.0     # enhancement multiplier at critical slope
    max_diffusivity::Float64 = 1000.0   # ceiling on D (m²/day) for CFL stability
    thermal_control::Bool    = true     # block diffusion where frozen
    thaw_threshold::Float64  = 0.01     # minimum ALT (m) to allow diffusion
    # ── Optional swash-zone enhancement ──
    # When `wet_critical_slope` is finite, cells where the wave-induced wet
    # probability `pwet` is between `swash_pwet_min` and `swash_pwet_max`
    # (i.e. the swash band — alternately wet & dry under wave action), plus
    # `swash_buffer_cells` landward, use the wet_* params instead of the dry
    # defaults — typically a much shallower critical slope so wave-loaded
    # sediment slumps while the dry bluff above keeps its steep-scarp
    # behaviour. Fully submerged cells (pwet > swash_pwet_max) are handled
    # by the underwater sediment-transport solver and are NOT diffused
    # here. Thermal control still applies (frozen → blocked).
    # Defaults = NaN preserve the legacy "block submerged cells" behaviour.
    wet_critical_slope::Float64  = NaN
    wet_d_base::Float64          = NaN
    wet_max_diffusivity::Float64 = NaN
    swash_buffer_cells::Int      = 0
    swash_pwet_min::Float64      = 0.01    # below this → "dry" cell
    swash_pwet_max::Float64      = 0.99    # above this → fully submerged
    # ── Optional pwet scaling on the wet-zone diffusivity ──
    # When `true`, multiply the wet-zone D by `clamp(pwet, pwet_scale_floor, 1.0)`
    # before use. Creates a gradient-driven flux divergence on a straight
    # slope, enabling continuous swash-face erosion. Originally added in
    # commit 441175f for the Tuk hindcast; calibrations such as
    # `wet_max_diffusivity=50` paired with `pwet=0.3 → D≈15` expect this.
    # Default `false` preserves the safer uniform-D behavior; set to `true`
    # to recover the Tuk-style continuous swash-face erosion. Floor
    # default 0.05 keeps the upper-swash edge at 5% D rather than 2% so
    # bluff-feed transport doesn't collapse entirely.
    wet_pwet_scaling::Bool         = false
    wet_pwet_scale_floor::Float64  = 0.05
end

"""
    CohesiveSedimentConfig

Minimal cohesive (mud) sediment configuration. A single suspended fraction
with Partheniades-Krone erosion / deposition. Kernel implementation lives
in `src/sediment/cohesive.jl`.

# Fields
- `settling_velocity::Float64`   — w_s (m/s). Constant floc settling.
- `density::Float64`             — dry sediment density (kg/m³, default 2650).
- `tau_ce::Float64`              — critical shear for erosion (Pa).
- `tau_cd::Float64`              — critical shear for deposition (Pa). Must
                                   satisfy `0 < tau_cd ≤ tau_ce` (the
                                   no-net-change buffer band).
- `M::Float64`                   — erosion rate constant (kg/m²/s).
- `initial_bed_mass::Float64`    — kg/m² of mud per node at t=0. Uniform
                                   fill; overwrite the state field for
                                   spatial variation.
- `rho_water::Float64`           — overlying water density (kg/m³, default
                                   1025). Used to convert CSHORE's
                                   normalized `tbxsta` (length units) to Pa.

# Intentional v1 limitations
- Single cohesive fraction (no floc-size dynamics)
- Constant w_s (no hindered settling)
- No bed consolidation (constant τ_ce)
- No advection of suspended concentration
- No feedback into `state.zb` (mass tracked in dedicated state vectors)

Suitable for screening-level mud deposition / resuspension studies;
not production estuarine modeling.
"""
Base.@kwdef struct CohesiveSedimentConfig
    settling_velocity::Float64 = 0.5e-3
    density::Float64           = 2650.0
    tau_ce::Float64            = 0.15
    tau_cd::Float64            = 0.06
    M::Float64                 = 1e-4
    initial_bed_mass::Float64  = 100.0
    rho_water::Float64         = 1025.0
end

"""
    SnowConfig

Simple snow insulation model. Snow accumulates when T_air < 0 and melts
when T_air > 0 (degree-day model), capped at `max_depth`. The snow layer
adds a thermal resistance `R = depth / k_snow` in series with the surface
boundary condition, reducing winter heat loss and (slightly) slowing summer
warming. Wave action in the swash zone removes snow via `pwet`.

When a prescribed `snow_depth` time series is provided in the boundary
forcing, it overrides the degree-day model entirely.
"""
Base.@kwdef struct SnowConfig
    k_snow::Float64      = 0.15    # thermal conductivity (W/m/K), ~fresh snow
    max_depth::Float64   = 0.3     # maximum snow depth cap (m)
    accum_rate::Float64  = 0.002   # accumulation when T_air < 0 (m/hour)
    melt_rate::Float64   = 0.005   # melt when T_air > 0 (m/hour per °C above 0)
end

"""
    SnowSpatialModifier

Per-node spatial overrides for snow depth, used to model snow-management NbS
interventions (snow fences, engineered drifts, snow clearing). Applied
inside `update_snow!` after the degree-day or prescribed depth is computed
and before the swash-pwet reduction:

    sd := clamp(sd, depth_min[j], depth_max[j])

- `depth_min[j]` — per-node forced minimum snow depth (m). Used by snow
  fences and insulating drifts to model engineered accumulation.
- `depth_max[j]` — per-node maximum cap (m). Use `0.0` for clearing zones
  (snow removal) or `Inf` for "no cap".

Construct via [`nbs_snow_fence`](@ref), [`nbs_snow_clearing`](@ref), or
[`nbs_insulating_drift`](@ref); compose by passing `base=` to a subsequent
call. Final cap of `snow_config.max_depth` still applies on top.
"""
struct SnowSpatialModifier
    depth_min::Vector{Float64}
    depth_max::Vector{Float64}
end

Base.length(m::SnowSpatialModifier) = length(m.depth_min)

"""
    ThermalBoundaryTimeSeries

Surface thermal forcing over time. `time` is in seconds since simulation
start (same clock as `BoundaryTimeSeries.timebc`). `T_air` (°C) drives
exposed nodes; `T_water` (°C) drives submerged nodes. Optional `snow_depth`
(m) prescribes snow cover, overriding the degree-day model when provided.
"""
Base.@kwdef struct ThermalBoundaryTimeSeries
    time::Vector{Float64}
    T_air::Vector{Float64}
    T_water::Vector{Float64}
    snow_depth::Vector{Float64} = Float64[]   # optional prescribed snow depth (m)
end

"""
    TidalInput

Tidal forcing time series for ITIDE=1. Applied to the longshore momentum
equation (via the `DETADY` alongshore water-surface gradient) and, when
`ILAB=0`, to the cross-shore volume flux via `DSWLDT` (rate of SWL change
→ `QTIDE` flux across each column).

- `time` (s): time grid (must span the wave BC time range)
- `detady` (dimensionless): alongshore water surface gradient `∂η/∂y`,
  drives longshore currents independent of wave action
- `dswldt` (m/s): rate of change of still-water level. Optional — if empty
  and `ILAB=0`, it's auto-computed from the SWL time series via finite
  difference. Only used when `ILAB=0` (lab experiments have no real tide).
"""
Base.@kwdef struct TidalInput
    time::Vector{Float64}
    detady::Vector{Float64}
    dswldt::Vector{Float64} = Float64[]
end

"""
    CurrentInput

Imposed alongshore-current forcing for ICURRENT=1. The user prescribes a target
longshore current speed `vbc` (m/s) at the offshore boundary as a time series,
and the model back-calculates the alongshore water-surface gradient `DETADY`
that produces it via analytical inversion of the longshore-momentum balance
(see `gby_from_vsigt` in `hydrodynamics/waves.jl`).

This is the inverse of ITIDE=1: instead of "give me detady, get back V", you
give it "I want V, what detady is required?". Useful when the ambient ocean
or estuarine current is known (from ADCPs, regional models, or river
discharge) but the local surface slope is not.

- `time` (s): time grid (must span the wave BC time range)
- `vbc` (m/s): alongshore current speed at the offshore boundary. Sign
  follows the same convention as `vmean`/`vsta` — positive = +y direction.
"""
Base.@kwdef struct CurrentInput
    time::Vector{Float64}
    vbc::Vector{Float64}
end

#===============================================================================
Aeolian (wind-driven) sediment transport

Adds wind-driven cross-shore sediment transport on the dry beach above the
runup line. Loosely modeled on the Dune Response Tool
(https://github.com/erdc/dune-response-tool): Kawamura (1951) transport
formula with Delgado-Fernández fetch limitation. Two extensions over DRT:

  • Per-fraction operation. With `MultifractionConfig.nf > 1`, each grain
    size has its own threshold velocity and contributes to flux based on
    its surface availability (active layer composition). This gives
    armouring "for free" — fines move first, coarse stays put.

  • Bypass mode. If no vegetation contour is present (or the wind blows a
    direction where it never crosses the contour), aeolian flux exits
    the model domain at the appropriate boundary and is removed from the
    sediment budget. Tracked as a diagnostic.

Vegetation-mediated capture is modeled via an abstract
`AeolianVegetationModel` so later phases can add density-based vegetation,
sheltering, etc., without changing the kernel.
===============================================================================#

"""
Abstract base type for the aeolian vegetation-capture model. Different
concrete types encode different physical assumptions about where
wind-blown sand stops:

- `ContourVegetation`          : Heaviside capture above an elevation contour (DRT-like).
- `DensityVegetation`          : continuous capture based on cover fraction.
- `MultiSpeciesVegetation`     : per-species biomass + aggregated capture.

The kernel calls `capture_efficiency(model, j, state, l)` to get a per-cell
capture fraction `η ∈ [0, 1]`: the fraction of incoming flux that is
deposited at this cell. Anything not captured passes downwind.
"""
abstract type AeolianVegetationModel end

"""
    ContourVegetation(z_contour=NaN; decay_length=0.0)

DRT-style capture. Above `z_contour` (m, in the same datum as `bed.dep`),
sand capture is active; below, flux passes through unchanged.

- `z_contour` — elevation marking the seaward edge of the capture zone (m).
  In nature this is typically the dune-toe contour, identified by the
  seaward limit of perennial vegetation. Set to `NaN` to disable capture
  entirely (wind blows sand right out of the domain).
- `decay_length` — exponential ramp-up of capture efficiency above the
  contour (m). 0 = sharp Heaviside (all flux deposits in the first
  capture-zone cell). >0 = `η = 1 - exp(-(z - z_contour)/decay_length)`,
  spreading deposition over a transition band — a common engineering
  smoothing.
"""
Base.@kwdef struct ContourVegetation <: AeolianVegetationModel
    z_contour::Float64 = NaN
    decay_length::Float64 = 0.0
end

"""
    DensityVegetation(; vegrho, alpha=1.0, ...)

Continuous, density-based capture. Each cell has a fractional vegetation
cover `ρ_veg ∈ [0, 1]` and the per-cell capture efficiency follows a Raupach-
style saturating exponential:

    η(j) = 1 - exp(-α · ρ_veg(j))

so that low-density cells capture a small fraction of incoming flux,
high-density cells capture nearly all of it, and bare cells (ρ_veg = 0)
let flux pass through. `alpha` controls the steepness — large `alpha`
makes the response close to a Heaviside (similar to `ContourVegetation`),
small `alpha` (~1) gives a smooth ramp.

Optional dynamic evolution of `vegrho` (when `dynamics_enabled = true`),
vegetation density evolves between BC windows according to:

    dρ/dt = (ρ_eq - ρ) / τ_growth    -    max(0, dz_aeolian/dt) · k_burial
                  (regrowth                 (burial: sand suppresses cover)
                  toward equilibrium)

Burial-driven loss kicks in when the local aeolian deposition rate
exceeds vegetation's ability to keep up. Bare beach (ρ_eq = 0) stays
bare; the back-beach / dune (ρ_eq > 0) regenerates after burial.

- `vegrho::Vector{Float64}` — required, length = number of grid nodes.
  Initial vegetation cover. Set ρ_veg = 0 on the seaward beach and > 0
  on the back-beach / dune.
- `alpha::Float64 = 1.0` — Raupach steepness for capture efficiency.
- `dynamics_enabled::Bool = false` — turn dynamic ρ evolution on.
- `vegrho_eq::Vector{Float64} = vegrho` — equilibrium cover; if not
  supplied, the initial `vegrho` is used as the equilibrium target.
- `tau_growth::Float64 = 31536000.0` — regrowth time constant (s).
  Default = 1 year; vegetation typically reaches equilibrium over a
  growing season to a few seasons.
- `k_burial::Float64 = 1.0` — burial sensitivity (per metre of bed-
  elevation gain). Larger values make vegetation more sensitive to
  burial (e.g. weakly rooted grasses).
"""
Base.@kwdef struct DensityVegetation <: AeolianVegetationModel
    vegrho::Vector{Float64}
    alpha::Float64 = 1.0
    dynamics_enabled::Bool = false
    vegrho_eq::Vector{Float64} = Float64[]
    tau_growth::Float64 = 31536000.0
    k_burial::Float64 = 1.0
end

"""
    VegetationSpecies(name; biomass_eq, ...)

Per-species parameter set used by `MultiSpeciesVegetation`. Separates state
(biomass, cover) from per-species response coefficients. Three default
presets are exposed below: `species_dune_grass()`, `species_shrub()`,
`species_forb()`.

State is **biomass** (kg/m² above-ground); **cover fraction** is a
diagnostic computed via `cover = 1 − exp(−c · B / B_eq)`. Capture
efficiency is then `α · cover`, summed across species.

- `name::String` — human-readable label, e.g. "dune_grass".
- `biomass_eq::Float64` — equilibrium biomass (kg/m²) at the species'
  preferred growing conditions. Coastal dune-grass: ~0.5; shrub:
  ~3.0; forb: ~0.2.
- `canopy_h_at_eq::Float64` — canopy height (m) at equilibrium biomass.
  Sets the roughness scale (z₀ ≈ 0.1·h) and the frontal-area-index
  contribution.
- `alpha::Float64` — Raupach capture-efficiency steepness for THIS
  species (replaces the single-species scalar `alpha` of
  `DensityVegetation`).
- `cover_per_biomass::Float64` — exponential coefficient `c`
  controlling how quickly cover saturates with biomass. Higher c =
  cover hits 1.0 sooner. Coastal grasses ~3-5; sparse shrubs ~1-2.
- `tau_growth::Float64` — biomass-recovery time constant (s). Coastal
  grasses recover in ~0.5-1 yr; shrubs take 3-10 yr.
- `k_burial::Float64` — burial sensitivity per metre of bed-elevation
  gain. Brittle forbs have high k (≥5); grasses tolerate burial
  (k ≤ 1).
- `frontal_area_factor::Float64` — multiplier on `cover · canopy_h`
  when computing the species' contribution to the aggregated capture
  Raupach term. Default 1.0; tune for species with unusual canopy
  architecture.
"""
Base.@kwdef struct VegetationSpecies
    name::String
    biomass_eq::Float64                              # kg/m²
    canopy_h_at_eq::Float64                          # m
    alpha::Float64 = 4.0
    cover_per_biomass::Float64 = 4.0                 # 1/(B/B_eq); cover→1-1/e at B/B_eq=0.25
    tau_growth::Float64 = 31536000.0                 # 1 yr
    k_burial::Float64 = 1.0
    frontal_area_factor::Float64 = 1.0
end

"""
    species_dune_grass(; kwargs...) -> VegetationSpecies

Preset matching coastal pioneer dune-grass species (e.g. *Ammophila*
on temperate coasts, *Leymus* on Arctic coasts). Burial-tolerant,
quick recovery, modest canopy.
"""
species_dune_grass(; kwargs...) = VegetationSpecies(;
    name = "dune_grass",
    biomass_eq = 0.5, canopy_h_at_eq = 0.4,
    alpha = 4.0, cover_per_biomass = 4.0,
    tau_growth = 1.0 * 365 * 86400.0,                # ~1 yr
    k_burial = 1.0,                                   # burial-tolerant
    frontal_area_factor = 1.0,
    kwargs...)

"""
    species_shrub(; kwargs...) -> VegetationSpecies

Preset for back-dune shrub vegetation (e.g. *Hudsonia*, beach plum).
High biomass, slow recovery, taller canopy than grass — strong
aerodynamic shelter when established but slow to re-establish after
disturbance.
"""
species_shrub(; kwargs...) = VegetationSpecies(;
    name = "shrub",
    biomass_eq = 3.0, canopy_h_at_eq = 1.2,
    alpha = 6.0, cover_per_biomass = 2.0,
    tau_growth = 5.0 * 365 * 86400.0,                # ~5 yr
    k_burial = 0.5,                                   # somewhat burial-tolerant
    frontal_area_factor = 1.5,                        # taller canopy → more shelter
    kwargs...)

"""
    species_forb(; kwargs...) -> VegetationSpecies

Preset for short-lived forbs / annual flowering plants. Low biomass,
fast turnover, brittle — easily killed by burial and slow to provide
substantial shelter.
"""
species_forb(; kwargs...) = VegetationSpecies(;
    name = "forb",
    biomass_eq = 0.2, canopy_h_at_eq = 0.15,
    alpha = 2.0, cover_per_biomass = 3.0,
    tau_growth = 0.5 * 365 * 86400.0,                # ~6 mo
    k_burial = 5.0,                                   # brittle
    frontal_area_factor = 0.8,
    kwargs...)

"""
    MultiSpeciesVegetation(; species, biomass, biomass_eq, dynamics_enabled=false)

Multi-species vegetation model. Each cell tracks per-species biomass (kg/m²);
cover is computed diagnostically and the aggregated capture efficiency adds
Raupach contributions across species:

    cover_s(j) = 1 − exp(−c_s · B[j,s] / B_eq,s)
    λ(j)       = Σ_s f_s · α_s · cover_s(j)
    η(j)       = 1 − exp(−λ(j))

This generalises `DensityVegetation`. With one species and tuned
parameters it can reproduce the single-density behaviour exactly,
but the natural multi-species form lets users distinguish e.g. a
burial-tolerant grass from a brittle forb at the same surface
density.

Dynamics (when `dynamics_enabled = true`): each species evolves
independently with first-order regrowth and burial-driven loss:

    dB_s/dt = (B_eq,s − B_s) / τ_growth,s   −   max(0, dz_aeolian/dt) · k_burial,s

For the simplest multi-species use case, build with two arrays:
species (length nspecies) and biomass (n_grid × n_species).
"""
Base.@kwdef struct MultiSpeciesVegetation <: AeolianVegetationModel
    species::Vector{VegetationSpecies}
    biomass::Matrix{Float64}        # (nn, nspecies) kg/m²
    biomass_eq::Matrix{Float64}     # (nn, nspecies) kg/m² (per-cell equilibrium target)
    dynamics_enabled::Bool = false
end

"""
    AeolianConfig(; vegetation, kwargs...)

Top-level configuration for IAEOLIAN=1.

- `vegetation::AeolianVegetationModel` — required. Specify `ContourVegetation(...)`.
- `Ck::Float64 = 1.8` — Kawamura (1951) coefficient (dimensionless). DRT
  exposes this as the only main tunable.
- `z_meas::Float64 = 10.0` — anemometer height for the wind-speed time
  series (m). Used in the log-law to compute `u*`.
- `karman::Float64 = 0.4` — von Kármán constant.
- `z0_factor::Float64 = 2.0/30.0` — `z0 = z0_factor · D50` (Nikuradse
  roughness parameterization).
- `rho_air::Float64 = 1.225` — air density (kg/m³).
- `moisture::Float64 = 0.0` — dimensionless surface moisture content
  (0..1). When `dry_time = 0`, this is a uniform static value applied
  to every cell. When `dry_time > 0`, `moisture` is the **wet-cell**
  value (used while a cell is at or below SWL / inside the swash zone)
  — see `dry_time` and `moisture_dry`. Modifies the threshold velocity
  via Belly-Johnson: `u*_t' = u*_t · √(1 + C·M)` with `C = 1.87`.
- `dry_time::Float64 = 0.0` — characteristic drying time (s) over which
  freshly wetted cells return to the dry-cell threshold. Set > 0 to
  enable per-cell dynamic moisture tracking: each aeolian step
  identifies cells currently below SWL or within the swash band and
  resets their wet-time clock; for cells above the runup line, M(j)
  decays linearly from `moisture` (at the moment of wetting) toward
  `moisture_dry` over `dry_time` seconds. Default 0 = static-moisture
  behaviour.
- `moisture_dry::Float64 = 0.0` — asymptotic moisture content of fully-
  dry cells. Reached `dry_time` seconds after the last wetting event.
  Typical: 0 for sand (truly dry), > 0 for cohesive bed types.
- `iuth_sheltering::Bool = false` — when `true`, raise the per-fraction
  threshold u*_t,k by Raupach (1993)-style sheltering from coarser
  non-erodible fractions in the same cell. Adds a multiplicative
  factor `√(1 + m·σ·λ)` where λ is the area-fraction of non-erodible
  grains (those with u*_t,k > current u*) and `m·σ ≈ 1` for sand-on-
  sand mixes. Only matters when multifraction (nf > 1). When `false`
  each fraction's threshold is fixed by D50,k alone (and moisture).
- `iuth_bedslope::Bool = false` — when `true`, modify the per-fraction
  threshold by the bed slope using the Dyer (1986) form:

      u*_t' = u*_t · √(cos α · (1 ± tan α / tan φ_r))

  with `+` on slopes that *oppose* the wind direction (uphill, harder
  to entrain) and `-` on slopes that *support* it (downhill, easier).
  `tan φ_r` is the angle-of-repose tangent from `SedimentConfig.tanphi`.
  Off by default — bed slope effects are typically modest on coastal
  beaches but become important on dune slip faces.
- `sheltering_msigma::Float64 = 1.0` — Raupach `m·σ` parameter, when
  `iuth_sheltering = true`. Larger values → stronger sheltering for
  the same coarse-fraction coverage.
- `fetch_critical_a::Float64 = 4.38`, `fetch_critical_b::Float64 = -8.23` —
  Delgado-Fernández (2011): `Fc = a·u_w + b`. Below `Fc`, flux is
  reduced by `sin(π·F/2Fc)`.
- `runup_buffer::Float64 = 0.0` — optional elevation offset (m) added
  above CSHORE's own runup tracker when defining the seaward edge of
  the dry beach. The base runup node is taken directly from CSHORE:
  `state.jdry` when `IOVER=1` (post-swash) or `state.jswl` when
  `IOVER=0`. No parametric runup formula is applied on top — the
  wave-model output is taken as fact. Leave at 0 unless you want a
  safety margin above the wet/dry boundary.
- `saturation_length::Float64 = 5.0` — Q(x) approaches local capacity as
  `Q = Q_cap · (1 - exp(-x_active/L_sat))`. Wind blowing over a long
  uniform sand surface saturates within a few meters; downwind of that
  ramp, the bed is in **bypass** (entrainment ≡ deposition, no net Δz).
  Only the upwind portion (within `L_sat`) shows actual erosion.
- `dune_decay_length::Float64 = 5.0` — characteristic length scale of
  the dune-toe deposit. Interpretation depends on `deposition_shape`:
  * `:exponential` — exponential e-folding length: `Q(x) = Q_in · exp(-η·Δx/L)`
  * `:triangular`  — total base length L of the triangular footprint
                     (DRT MorphologyResolving Style 1: linear rise from
                     0 at the seaward edge to peak at L/3, linear fall
                     to 0 at L). Mass-conserving — total area equals the
                     incoming saturated flux × dt.
  * `:right_triangle` — base length L, peak at the seaward edge, linear
                        falling slope to 0 at L (DRT Style 2).
  * `:gaussian`    — DRT MorphologyResolving Style 3: a Style-1 triangle
                     followed by a Gaussian smoothing pass that softens
                     the linear edges into a smoother dune-toe profile.
                     The smoothing is mass-conserving — the deposited
                     volume is rescaled after the convolution so the
                     total ∫Δz dx exactly matches the incoming flux × dt.
- `deposition_shape::Symbol = :exponential` — selects the kernel above.
  `:exponential` is the simplest and matches DRT-main behaviour;
  `:triangular` mirrors the DRT MorphologyResolving (Style 1) branch
  and gives a more realistic dune profile with the peak set back from
  the seaward edge by L/3 (the natural slip-face geometry);
  `:gaussian` is DRT Style 3 (Style 1 + Gaussian smoothing).
- `gaussian_smooth_sigma::Float64 = 2.0` — Gaussian smoothing standard
  deviation in metres, used only when `deposition_shape == :gaussian`.
  The kernel is a Gaussian with this σ, truncated at ±3σ, normalized
  so its discrete sum equals 1, and applied to the triangular profile.
  The post-smoothing deposit is then rescaled to preserve the exact
  incoming-flux mass.
- `veg_deposition_center::Float64 = 0.0` — cross-shore offset (m) shifting
  the deposition zone landward of the vegetation contour. Mirrors the
  `VegDepositionCenter` parameter from the DRT MorphologyResolving branch
  (https://github.com/erdc/dune-response-tool/tree/MorphologyResolving):
  deposition begins at `x_contour + veg_deposition_center` instead of
  right at the contour. Cells between the contour and the offset
  position act as a vegetation-protected bypass zone (no entrainment,
  no deposition). Useful when the vegetation line marks where wind
  starts to slow but actual sand capture (e.g. by a denser back-dune)
  occurs further inland. Set to 0 for DRT-main behaviour. Negative
  values are not yet supported.
- `dt_aeolian_max::Float64 = Inf` — sub-step cap (s). When the BC window
  is long, the kernel sub-steps to keep bed change per step bounded.
  `Inf` = single integration over the full BC window.
"""
Base.@kwdef struct AeolianConfig
    vegetation::AeolianVegetationModel
    Ck::Float64 = 1.8
    z_meas::Float64 = 10.0
    karman::Float64 = 0.4
    z0_factor::Float64 = 2.0 / 30.0
    rho_air::Float64 = 1.225
    moisture::Float64 = 0.0
    dry_time::Float64 = 0.0
    moisture_dry::Float64 = 0.0
    iuth_sheltering::Bool = false
    iuth_bedslope::Bool = false
    sheltering_msigma::Float64 = 1.0
    fetch_critical_a::Float64 = 4.38
    fetch_critical_b::Float64 = -8.23
    runup_buffer::Float64 = 0.0
    saturation_length::Float64 = 5.0
    dune_decay_length::Float64 = 5.0
    deposition_shape::Symbol = :exponential
    gaussian_smooth_sigma::Float64 = 2.0
    veg_deposition_center::Float64 = 0.0
    dt_aeolian_max::Float64 = Inf
end

"""
    WindShearConfig(; ...)

Configuration for the wind-flow-over-topography solver (`IWINDSHEAR=1`).
Uses the Kroy-Sauermann-Herrmann (2002) "minimal model" perturbation for
1D, which is what DUNA also uses.

- `kroy_alpha::Float64 = 3.0` — overall amplitude. Tunes how strongly
  the local shear responds to bed-slope features. DUNA default 3.0.
- `kroy_beta::Float64 = 0.2` — local-slope contribution weight. Controls
  stoss-vs-lee asymmetry of the perturbation. DUNA default 0.2.
- `tau_clamp_floor::Float64 = 0.1` — minimum allowed `τ/τ₀` (default 0.1).
  Prevents unphysical negative shear when the Hilbert transform produces
  large negative values in deep lee zones; the lee separation mask is the
  proper handler for such cells.
- `mask_bathymetry::Bool = true` — pass `max(zb, swl)` to the shear solver
  rather than the raw bed. Submerged bars / scour holes don't generate
  spurious shear bumps. Matches DUNA's `duna.m:113` convention.
- `lee_slope::Float64 = 0.4` — recirculation-bubble slope (~22°). Cells
  under a downwind plane projected from each landward-facing brink at
  this slope have `τ` set to 0 (no transport).
- `min_brink_drop::Float64 = 0.0` — only consider brinks where the bed
  drops at least this much within ~5 cells downwind. Set to a small value
  (e.g. 0.05 m) to filter out ripple-scale brinks.
- `tau_smooth_sigma::Float64 = 0.0` — Gaussian smoothing standard
  deviation in metres applied to τ/τ₀ after the Kroy convolution. The
  discrete Hilbert transform tends to amplify grid-scale noise in
  `∂z/∂x`, producing spurious sub-meter spikes in τ/τ₀ even on smooth
  profiles. Set σ ≈ 1–3 m to filter those out without broadening the
  macroscopic crest speed-up. Mass-conserving rescale (the smoothed
  area equals the raw area). Default 0 disables smoothing.
- `shear_method::Symbol = :duna_kroy` — selects the wind-shear closure:
  * `:duna_kroy`  — fixed `(α, β)` from `kroy_alpha` / `kroy_beta`.
  * `:aeolis_kroy` — Kroy 2002 self-consistent `(α, β)` from feature
    length scale `aeolis_length_scale` and local roughness `z0`.
    Same Hilbert kernel; the differentiator is principled α, β.
- `aeolis_length_scale::Float64 = 20.0` — characteristic feature half-
  length L (m) used by `:aeolis_kroy`. For coastal foredunes ~10–30 m.
- `kroy_solver::Symbol = :direct` — Hilbert-transform implementation:
  * `:direct` — O(N²) discrete real-space convolution (DUNA-style).
                No boundary artifacts. Default for backward compat.
  * `:fft`    — O(N log N) spectral implementation via FFTW. Uses
                even-symmetric reflection of the slope field around
                the boundaries to suppress Gibbs ringing. Materially
                faster for N ≳ 500 cells.
"""
Base.@kwdef struct WindShearConfig
    kroy_alpha::Float64 = 3.0
    kroy_beta::Float64 = 0.2
    tau_clamp_floor::Float64 = 0.1
    mask_bathymetry::Bool = true
    lee_slope::Float64 = 0.4
    min_brink_drop::Float64 = 0.0
    tau_smooth_sigma::Float64 = 0.0
    shear_method::Symbol = :duna_kroy
    aeolis_length_scale::Float64 = 20.0
    kroy_solver::Symbol = :direct
end

"""
    CshoreConfig

Top-level immutable configuration. Constructed by `read_infile` from a FORTRAN
`.infile`, or programmatically for tests.
"""
Base.@kwdef struct CshoreConfig
    options::OptionFlags
    grid::GridConfig
    sediment::SedimentConfig
    multifraction::MultifractionConfig = MultifractionConfig()
    boundary::BoundaryTimeSeries
    bathymetry::BathyInput
    porous::Union{Nothing, PorousInput} = nothing
    vegetation::Union{Nothing, VegetationInput} = nothing
    dike::Union{Nothing, DikeErosionInput} = nothing
    clay::Union{Nothing, ClayInput} = nothing
    tidal::Union{Nothing, TidalInput} = nothing
    current::Union{Nothing, CurrentInput} = nothing
    aeolian::Union{Nothing, AeolianConfig} = nothing
    windshear::Union{Nothing, WindShearConfig} = nothing
    ig::Union{Nothing, IgConfig} = nothing    # infragravity wave energy (Layer 1/2/3)
    undertow::Union{Nothing, UndertowConfig}  = nothing  # near-bed return-flow amplification
    asymm::Union{Nothing, AsymmetryConfig}    = nothing  # DEPRECATED — use wave_nonlinearity
    phase_lag::Union{Nothing, PhaseLagConfig} = nothing  # DEPRECATED — use wave_nonlinearity
    bailard::Union{Nothing, BailardConfig}    = nothing  # DEPRECATED — use wave_nonlinearity
    # Unified wave-nonlinearity surface (recommended). When set, this
    # overrides every legacy field — `facSK`, `facAS`, `ur_sk_ref`,
    # `alpha_sk`, `biphase_relax_L`, `OptionFlags.iasym`,
    # `OptionFlags.iskew_spatial`, and the `asymm` / `phase_lag` /
    # `bailard` sub-structs above. Access via `nonlinearity(config)` from
    # consumer code so the fallback path stays compatible with old infiles.
    wave_nonlinearity::Union{Nothing, WaveNonlinearityConfig} = nothing
    groundwater::Union{Nothing, GroundwaterConfig} = nothing  # beach aquifer + surface moisture
    swash::SwashConfig = SwashConfig()
    thermal::Union{Nothing, ThermalConfig} = nothing
    thermal_bc::Union{Nothing, ThermalBoundaryTimeSeries} = nothing
    snow::Union{Nothing, SnowConfig} = nothing
    snow_modifier::Union{Nothing, SnowSpatialModifier} = nothing
    diffusion::Union{Nothing, DiffusionConfig} = nothing
    cohesive::Union{Nothing, CohesiveSedimentConfig} = nothing  # mud / cohesive sediment
    # Wave breaking / stability
    gamma::Float64 = 0.78         # breaker index (used as fallback when ruessink is off)
    # `gamma_method` selects how the local breaker index γ is computed:
    #   :constant       — fixed `gamma` (FORTRAN parity)
    #   :ruessink2003   — Ruessink, Walstra & Southgate (2003) kh-dependent form
    #                       γ(kh) = clamp(γ_a · kh + γ_b, γ_min, γ_max)
    #                     defaults match RWS2003 Eq. 7.
    #   :steepness_sr   — Symbolic-regression formula derived from a 13.4-year
    #                     FRF Duck NC inverse pipeline (n≈29k breaking-only obs):
    #                       γ(Hs/L) = clamp(γ_a · √(Hs/L), γ_min, γ_max)
    #                     with γ_a ≈ 3.90 (Cohn et al., 2026, in prep).
    #                     Hs is approximated from local Hrms (Hs = √2·Hrms);
    #                     L is the linear-theory wavelength 2π/k. Provides
    #                     ~–60 % surf-zone Hs RMSE vs constant γ=0.78 at FRF.
    gamma_method::Symbol = :constant
    gamma_a::Float64 = 0.76        # RWS2003 slope on kh
    gamma_b::Float64 = 0.29        # RWS2003 intercept
    gamma_min::Float64 = 0.35
    gamma_max::Float64 = 0.90
    # Coefficient for `:steepness_sr`: γ = gamma_sr_slope · √(Hs/L).
    # Default 3.9 is the SR-trained value on FRF 2008–2026 breaking-only data.
    gamma_sr_slope::Float64 = 3.9
    breaker_delay::Float64 = 0.0     # [0.0, 1.0] — ramping factor for dissipation onset
                                      # 0.0 = current behavior (sharp onset, no delay)
                                      # 1.0 = full smooth ramp over transition zone
    sismax::Float64 = 1.0
    # Morphology stability
    # ─ morph_courant : Courant-like safety factor in the adaptive
    #   morphodynamic timestep: Δt = morph_courant · Δx / |cb_max|.
    #   Lower = smaller Δt = more stable (suppresses 2Δx checkerboard
    #   and Exner ringing) at the cost of run time. Range typical 0.2-0.5;
    #   default 0.3 errs conservative (was 0.5 historically).
    # ─ max_dzb_per_step : hard guard against pathological per-step bed
    #   excursions (m). Not a CFL knob — it's a safety cap to catch
    #   unphysical transport blow-ups from numerical glitches before they
    #   propagate. Should rarely bind in well-resolved runs. Default 0.1 m.
    # ─ morph_diffusion : optional explicit smoothing applied each step
    #   (m²/s, default 0 = off). Effective against 2Δx waves; introduces
    #   small diffusion error.
    # ─ min_depth_morph : minimum water depth (m) for morphology to engage.
    morph_courant::Float64     = 0.3
    max_dzb_per_step::Float64  = 0.1
    morph_diffusion::Float64   = 0.0
    min_depth_morph::Float64   = 0.05
    # Minimum water depth for wave-current interaction (IWCINT=1).
    # Below this depth, qdisp is zeroed to prevent spurious WCI
    # amplification in very shallow water.
    min_depth_wcint::Float64   = 0.10
    # Wave asymmetry / skewness scaling (IASYM=1). The bedload gains an
    # additional term proportional to skewness × ustd³. facSK=1.0 applies
    # the Ruessink 2012 parameterization at full strength; lower values
    # scale the effect down. facAS scales the asymmetry (acceleration-
    # driven) contribution similarly.
    facSK::Float64 = 1.0           # skewness coupling strength [0-1]
    facAS::Float64 = 0.0           # asymmetry coupling strength [0-1] (default off)
    # Ursell-weighted skewness coupling (ISKEW_SPATIAL=1).
    # The effective coupling at each node is facSK·tanh(Ur / ur_sk_ref).
    # ur_sk_ref ≈ 0.2 gives half-weight at moderate shoaling nonlinearity
    # (kh ~ 0.8); increase to extend coupling further offshore.
    ur_sk_ref::Float64 = 0.20
    # ── Wave nonlinearity — breaking crest-height correction ──────────────────
    # When alpha_sk > 0, the maximum breaker height hm is reduced by the local
    # wave skewness before computing the breaking fraction:
    #   hm_eff = hm / (1 + alpha_sk · max(0, Sk))
    # A skewed wave has a higher crest than Hrms/2 implies, so it breaks at
    # lower Hrms. Typical literature values: 0.0 (off) to ~0.5.
    # Only active when iasym ≥ 1 (Sk computed from Ruessink or Stokes).
    alpha_sk::Float64 = 0.0
    # ── Wave nonlinearity — biphase relaxation along profile ──────────────────
    # When biphase_relax_L > 0, the Ruessink (2012) biphase β is not applied
    # instantaneously from the local Ursell number, but instead evolves along
    # the cross-shore profile with a spatial relaxation length L_relax (m):
    #   dβ/dx = (β_eq(Ur) - β) / L_relax
    # This captures the spatial lag between the nonlinear wave shape and the
    # local depth — waves stay asymmetric past the breakpoint and only relax
    # to the equilibrium skewed state over O(L_relax). Typical: 50–200 m.
    # 0.0 = static Ruessink (default, original behaviour).
    biphase_relax_L::Float64 = 0.0
    # Shields-regime-adaptive friction (IFRICTION_SPATIAL=1).
    # f_min     — grain-roughness friction factor (no bedforms); also used as
    #             the reference for Shields computation (breaks the fb2 circularity).
    # f_sheet   — friction factor in sheet-flow regime (high Shields).
    # theta_sheet — Shields number above which sheet flow is assumed.
    # f_ripple_exp — power-law exponent for ripple-regime friction:
    #             f = f_min · (θ / θ_cr)^f_ripple_exp   (θ_cr ≤ θ < theta_sheet)
    f_min::Float64        = 0.002   # grain-roughness baseline
    f_sheet::Float64      = 0.015   # sheet-flow friction cap
    theta_sheet::Float64  = 1.0     # Shields threshold for sheet flow
    f_ripple_exp::Float64 = 0.5     # ripple-regime power-law exponent
    # Dynamic Manning (ifriction_spatial=2):
    # fb2 = g·n²/h^(1/3) recomputed each step from state.h. h is floored by
    # manning_h_min (m) to avoid singularity near the wet/dry front; fb2 is
    # capped at manning_fb2_max to prevent runaway in very shallow flow.
    manning_h_min::Float64   = 0.05
    manning_fb2_max::Float64 = 0.1
    # Convergence tolerances
    eps1::Float64 = 1e-3
    eps2::Float64 = 1e-6
    maxite::Int   = 20
    # Iterative WCI (IWCINT_ALONG=1). Outer Picard iteration on the alongshore-
    # current Doppler term in the wave dispersion. Default tol of 5 mm/s and
    # 5 iterations typically converges in 2-3 sweeps for usual conditions.
    iwcint_along_tol::Float64    = 5e-3   # m/s; max |Δvmean| across nodes
    iwcint_along_maxite::Int     = 5
end

"""
    nonlinearity(config::CshoreConfig) -> WaveNonlinearityConfig

Return the effective wave-nonlinearity configuration for a `CshoreConfig`.

If `config.wave_nonlinearity` is set, returns it unchanged. Otherwise, fans
the legacy scattered fields (`OptionFlags.iasym`, `OptionFlags.iskew_spatial`,
top-level `facSK` / `facAS` / `ur_sk_ref` / `alpha_sk` / `biphase_relax_L`,
and the `asymm` / `phase_lag` / `bailard` sub-structs) into a single
`WaveNonlinearityConfig`. The fanout preserves the original behavior:
old infiles and scripts continue to work without modification.

All consumers (transport.jl, waves.jl, driver.jl) should read through
this helper rather than touching the legacy fields directly.
"""
function nonlinearity(config::CshoreConfig)
    config.wave_nonlinearity === nothing || return config.wave_nonlinearity

    # ── Fanout from legacy fields ──
    iasym = config.options.iasym
    closure_legacy = if config.asymm !== nothing
        config.asymm.method
    elseif iasym == 2
        :stokes2
    elseif iasym == 1
        :ruessink
    else
        :linear
    end
    enabled_legacy = (iasym >= 1) || (config.asymm !== nothing) ||
                     (config.bailard !== nothing && config.bailard.enabled) ||
                     (config.phase_lag !== nothing && config.phase_lag.enabled)

    # Sub-config field passthrough (use sub-struct values when present,
    # otherwise the defaults from WaveNonlinearityConfig).
    K_as = config.asymm !== nothing ? config.asymm.K_as : 0.4
    smooth_diss_window = config.asymm !== nothing ? config.asymm.smooth_diss_window : 3
    diss_grad_cap = config.asymm !== nothing ? config.asymm.diss_grad_cap : 0.8

    pl_en  = config.phase_lag !== nothing && config.phase_lag.enabled
    pl_tau = config.phase_lag !== nothing ? config.phase_lag.tau_lag   : NaN
    pl_lm  = config.phase_lag !== nothing ? config.phase_lag.L_min     : 0.5
    pl_cap = config.phase_lag !== nothing ? config.phase_lag.cap_alpha : 1.0

    ba_en  = config.bailard !== nothing && config.bailard.enabled
    ba_eps = config.bailard !== nothing ? config.bailard.eps_s    : 0.025
    ba_xc  = config.bailard !== nothing ? config.bailard.gamma_xc : 1.0
    ba_add = config.bailard !== nothing ? config.bailard.additive : true
    # Bailard's gamma_sk / gamma_as carry their own multipliers separate from
    # the global skewness/asymmetry strengths. When a BailardConfig is present
    # AND a legacy facSK/facAS are also set, the consumer code currently
    # multiplies them together; the fanout preserves that by leaving global
    # skewness/asymmetry as the configured top-level values.
    skew_strength  = config.facSK
    asym_strength  = config.facAS
    if config.bailard !== nothing
        skew_strength *= config.bailard.gamma_sk
        asym_strength *= config.bailard.gamma_as
    end

    return WaveNonlinearityConfig(;
        enabled              = enabled_legacy,
        closure              = closure_legacy,
        skewness             = skew_strength,
        asymmetry            = asym_strength,
        spatial_weighting    = config.options.iskew_spatial >= 1,
        ur_reference         = config.ur_sk_ref,
        alpha_crest          = config.alpha_sk,
        biphase_relax_length = config.biphase_relax_L,
        K_as                 = K_as,
        smooth_diss_window   = smooth_diss_window,
        diss_grad_cap        = diss_grad_cap,
        phase_lag_enabled    = pl_en,
        phase_lag_tau        = pl_tau,
        phase_lag_L_min      = pl_lm,
        phase_lag_cap        = pl_cap,
        bailard_enabled      = ba_en,
        bailard_eps_s        = ba_eps,
        bailard_gamma_xc     = ba_xc,
        bailard_additive     = ba_add,
    )
end

function validate(c::CshoreConfig)
    validate(c.multifraction)
    c.grid.dx > 0 || throw(ArgumentError("dx must be > 0"))
    0.05 ≤ c.morph_courant ≤ 1.0 || throw(ArgumentError(
        "morph_courant must be in [0.05, 1.0] (got $(c.morph_courant)); " *
        "default 0.3, lower = more stable / slower"))
    c.max_dzb_per_step > 0 || throw(ArgumentError(
        "max_dzb_per_step must be > 0 (got $(c.max_dzb_per_step))"))
    c.options.iline == size(c.bathymetry.xbinp, 2) ||
        throw(ArgumentError("iline mismatch with bathymetry.xbinp columns"))

    c.options.iover in (0, 1) || throw(ArgumentError(
        "iover must be 0 or 1 (got $(c.options.iover))"))
    if c.options.iover == 1
        # Crest auto-detected when rcrest is empty; warn so the user is aware.
        any(c.swash.rcrest .> -1e10) || @warn "IOVER=1 but no RCREST specified — crest will be auto-detected from bathymetry"
    end
    c.options.iperm in (0, 1) || throw(ArgumentError(
        "iperm must be 0 or 1 (got $(c.options.iperm))"))
    if c.options.iperm == 1
        c.porous !== nothing || @warn "iperm=1 but no PorousInput provided — using default parameters (porosity=0.4, Dn50=0.02m). Pass porous= or porous_z= to build_config for control."
    end
    if c.porous !== nothing && c.options.iperm == 0
        @warn "PorousInput provided but iperm=0 — porous flow will NOT be computed. Set iperm=1 in OptionFlags to activate."
    end
    c.options.iveg in (0, 1, 2, 3) || throw(ArgumentError(
        "iveg must be 0, 1, 2, or 3 (got $(c.options.iveg))"))
    if c.options.iveg > 0
        c.vegetation !== nothing || throw(ArgumentError(
            "iveg=$(c.options.iveg) requires a VegetationInput — pass vegetation= to build_config"))
        c.options.iveg == 3 && c.options.idiss ∉ (0, 1) && throw(ArgumentError(
            "iveg=3 currently supports idiss=1 only (Mendez-Losada); got idiss=$(c.options.idiss)"))
    end
    c.options.iwtran in (0, 1) || throw(ArgumentError(
        "iwtran must be 0 or 1 (got $(c.options.iwtran))"))
    if c.options.iwtran == 1
        c.options.iwtran_kt_method in (:dangremond_vandermeer, :goda, :freeboard_ratio) ||
            throw(ArgumentError("iwtran_kt_method must be :dangremond_vandermeer, :goda, " *
                                "or :freeboard_ratio (got :$(c.options.iwtran_kt_method))"))
        if !isempty(c.boundary.swl_landward) &&
           length(c.boundary.swl_landward) != length(c.boundary.timebc)
            throw(ArgumentError("boundary.swl_landward (length $(length(c.boundary.swl_landward))) " *
                                "must match timebc length ($(length(c.boundary.timebc))) when iwtran=1"))
        end
    end
    c.options.iroll in (0, 1) || throw(ArgumentError(
        "iroll must be 0 or 1 (got $(c.options.iroll))"))
    c.options.iwcint in (0, 1) || throw(ArgumentError(
        "iwcint must be 0 or 1 (got $(c.options.iwcint))"))
    c.options.iwcint_along in (0, 1) || throw(ArgumentError(
        "iwcint_along must be 0 or 1 (got $(c.options.iwcint_along))"))
    if c.options.iwcint_along == 1
        c.options.iwcint == 1 || throw(ArgumentError(
            "iwcint_along=1 requires iwcint=1 (cross-shore WCI must be enabled too)"))
        c.options.iangle == 1 || @warn(
            "iwcint_along=1 has no effect with iangle=0 (no oblique-wave projection of the alongshore current onto the wave-propagation direction)")
    end
    c.options.iwind in (0, 1) || throw(ArgumentError(
        "iwind must be 0 or 1 (got $(c.options.iwind))"))
    c.options.iasym in (0, 1, 2) || throw(ArgumentError(
        "iasym must be 0, 1 (Ruessink 2012), or 2 (Stokes 2nd-order); got $(c.options.iasym)"))
    if c.asymm !== nothing
        c.asymm.method in (:ruessink, :stokes2, :boussinesq_diss) || throw(ArgumentError(
            "asymm.method must be :ruessink, :stokes2, or :boussinesq_diss; got $(c.asymm.method)"))
        c.asymm.K_as >= 0.0 || throw(ArgumentError("asymm.K_as must be ≥ 0"))
        c.asymm.smooth_diss_window >= 0 || throw(ArgumentError("asymm.smooth_diss_window must be ≥ 0"))
        c.asymm.diss_grad_cap > 0.0 || throw(ArgumentError("asymm.diss_grad_cap must be > 0"))
    end
    if c.undertow !== nothing
        c.undertow.mode in (:hrms_h, :dissipation) || throw(ArgumentError(
            "undertow.mode must be :hrms_h or :dissipation; got $(c.undertow.mode)"))
        c.undertow.alpha >= 0.0 || throw(ArgumentError("undertow.alpha must be ≥ 0"))
        c.undertow.exponent >= 0.0 || throw(ArgumentError("undertow.exponent must be ≥ 0"))
        c.undertow.h_min > 0.0 || throw(ArgumentError("undertow.h_min must be > 0"))
        c.undertow.cap >= 1.0 || throw(ArgumentError("undertow.cap must be ≥ 1"))
    end
    if c.phase_lag !== nothing
        pl = c.phase_lag
        (isnan(pl.tau_lag) || pl.tau_lag > 0.0) ||
            throw(ArgumentError("phase_lag.tau_lag must be > 0 or NaN; got $(pl.tau_lag)"))
        pl.L_min > 0.0 || throw(ArgumentError("phase_lag.L_min must be > 0"))
        0.0 < pl.cap_alpha <= 1.0 || throw(ArgumentError(
            "phase_lag.cap_alpha must be in (0, 1]"))
    end
    if c.bailard !== nothing
        b = c.bailard
        b.eps_s >= 0.0 || throw(ArgumentError("bailard.eps_s must be ≥ 0"))
        # gamma_xc, gamma_sk, gamma_as can be any sign (allow negative for tests)
    end
    c.options.inl_dispersion in (0, 1) || throw(ArgumentError(
        "inl_dispersion must be 0 (linear) or 1 (Kirby-Dalrymple nonlinear); got $(c.options.inl_dispersion)"))
    c.alpha_sk >= 0.0 || throw(ArgumentError(
        "alpha_sk must be ≥ 0; got $(c.alpha_sk)"))
    c.biphase_relax_L >= 0.0 || throw(ArgumentError(
        "biphase_relax_L must be ≥ 0 (0 = static Ruessink); got $(c.biphase_relax_L)"))
    if c.alpha_sk > 0.0 && c.options.iasym == 0
        @warn "alpha_sk=$(c.alpha_sk) has no effect when iasym=0 (skewness not computed)"
    end
    0.0 ≤ c.breaker_delay ≤ 1.0 || throw(ArgumentError(
        "breaker_delay must be in [0.0, 1.0], got $(c.breaker_delay)"))
    c.options.infilt in (0, 1) || throw(ArgumentError(
        "infilt must be 0 or 1 (got $(c.options.infilt))"))
    if c.options.infilt == 1
        c.options.iperm == 0 || throw(ArgumentError(
            "infilt=1 requires iperm=0 (infiltration and porous layer are mutually exclusive)"))
        c.options.iover == 1 || throw(ArgumentError(
            "infilt=1 requires iover=1 (infiltration uses the iterated overtopping flux qo)"))
    end
    if c.groundwater !== nothing
        gw = c.groundwater
        gw.K  > 0 || throw(ArgumentError("groundwater.K must be > 0"))
        gw.ne > 0 || throw(ArgumentError("groundwater.ne must be > 0"))
        gw.D  > 0 || throw(ArgumentError("groundwater.D must be > 0"))
        (0.0 < gw.theta_res < gw.theta_sat <= 1.0) ||
            throw(ArgumentError("groundwater: require 0 < theta_res < theta_sat ≤ 1"))
        gw.vg_n > 1.0 || throw(ArgumentError("groundwater.vg_n must be > 1"))
        gw.dt_gw > 0  || throw(ArgumentError("groundwater.dt_gw must be > 0"))
        if gw.infiltration_rate !== nothing && c.options.infilt == 1
            @warn "groundwater.infiltration_rate is set explicitly while infilt=1. " *
                  "This may break mass conservation with the INFILT swash drain — the " *
                  "swash drains at wpm but the aquifer recharges at infiltration_rate. " *
                  "Leave infiltration_rate=nothing to auto-derive from wpm."
        end
        if gw.infiltration_rate !== nothing && gw.infiltration_rate < 0
            throw(ArgumentError("groundwater.infiltration_rate must be ≥ 0"))
        end
        if gw.rainfall_rate !== nothing
            nbc = length(c.boundary.timebc)
            if gw.rainfall_rate isa Float64
                gw.rainfall_rate < 0 &&
                    throw(ArgumentError("groundwater.rainfall_rate must be ≥ 0 (m/s)"))
            else  # Vector
                length(gw.rainfall_rate) == nbc ||
                    throw(ArgumentError("groundwater.rainfall_rate vector length " *
                        "($(length(gw.rainfall_rate))) must equal length(boundary.timebc) ($nbc)"))
                any(<(0), gw.rainfall_rate) &&
                    throw(ArgumentError("all groundwater.rainfall_rate values must be ≥ 0"))
            end
        end
        if gw.rainfall_wet_theta !== nothing
            nbc = length(c.boundary.timebc)
            _check_theta = (θ, tag) -> begin
                (gw.theta_res ≤ θ ≤ gw.theta_sat) ||
                    throw(ArgumentError("groundwater.rainfall_wet_theta $tag must be " *
                        "in [theta_res, theta_sat] = [$(gw.theta_res), $(gw.theta_sat)]"))
            end
            if gw.rainfall_wet_theta isa Float64
                _check_theta(gw.rainfall_wet_theta, "")
            else  # Vector
                length(gw.rainfall_wet_theta) == nbc ||
                    throw(ArgumentError("groundwater.rainfall_wet_theta vector length " *
                        "($(length(gw.rainfall_wet_theta))) must equal length(boundary.timebc) ($nbc)"))
                for (i, θ) in enumerate(gw.rainfall_wet_theta)
                    _check_theta(θ, "[$i]")
                end
            end
        end
        # Warn if aeolian is not active — the main consumer of theta
        if c.aeolian === nothing
            @warn "GroundwaterConfig is set but no AeolianConfig is present. " *
                  "The surface moisture field (state.theta) will be computed but " *
                  "not used for transport. Add aeolian= to build_config to couple."
        end
    end
    if c.options.infilt == 1
        c.options.iprofl == 1 || throw(ArgumentError(
            "infilt=1 requires iprofl=1 (morphodynamic mode)"))
    end
    if c.options.iwind == 1
        !isempty(c.boundary.w10) || throw(ArgumentError(
            "iwind=1 requires non-empty w10 wind speed array"))
    end

    # IPROFL=2 (grassed dike erosion / EROSON): supported with DikeErosionInput
    c.options.iprofl in (0, 1, 2) || throw(ArgumentError(
        "iprofl must be 0, 1, or 2 (got $(c.options.iprofl))"))
    if c.options.iprofl == 2
        c.dike !== nothing || throw(ArgumentError(
            "iprofl=2 (dike erosion) requires a DikeErosionInput — pass dike= to build_config"))
        c.options.iclay == 0 || throw(ArgumentError(
            "iprofl=2 is incompatible with iclay=1"))
    end
    # ICLAY=1 (sand over clay / EROSON clay branch)
    c.options.iclay in (0, 1) || throw(ArgumentError(
        "iclay must be 0 or 1 (got $(c.options.iclay))"))
    if c.options.iclay == 1
        c.clay !== nothing || throw(ArgumentError(
            "iclay=1 requires a ClayInput — pass clay= to build_config"))
        c.options.iperm == 0 || throw(ArgumentError(
            "iclay=1 requires iperm=0 (cannot combine with porous layer)"))
        c.options.iprofl == 1 || throw(ArgumentError(
            "iclay=1 requires iprofl=1 (morphodynamic sand bed)"))
        abs(c.options.isedav) == 1 || throw(ArgumentError(
            "iclay=1 requires isedav=±1 (sand-over-hardbottom framework)"))
    end
    c.options.isedav == 2 && throw(ArgumentError("isedav=2 (wire mesh) not supported"))
    c.options.ipond  == 1 && throw(ArgumentError("ipond=1 (ridge-runnel) not supported"))
    c.options.itide in (0, 1) || throw(ArgumentError(
        "itide must be 0 or 1 (got $(c.options.itide))"))
    if c.options.itide == 1
        c.tidal !== nothing || throw(ArgumentError(
            "itide=1 requires a TidalInput — pass tidal= to build_config"))
        c.options.iangle == 1 || @warn(
            "itide=1 has no longshore effect with iangle=0 (normal incidence)")
        c.options.ilab == 1 && @warn(
            "itide=1 with ilab=1 (lab mode) skips the QTIDE cross-shore term")
    end
    c.options.itide  == 1 && c.tidal === nothing && throw(ArgumentError(
        "itide=1 requires tidal=TidalInput(...)"))
    c.options.icurrent in (0, 1) || throw(ArgumentError(
        "icurrent must be 0 or 1 (got $(c.options.icurrent))"))
    c.options.iaeolian in (0, 1) || throw(ArgumentError(
        "iaeolian must be 0 or 1 (got $(c.options.iaeolian))"))
    if c.options.iaeolian == 1
        c.aeolian !== nothing || throw(ArgumentError(
            "iaeolian=1 requires an AeolianConfig — pass aeolian=AeolianConfig(vegetation=...)"))
        !isempty(c.boundary.w10) || throw(ArgumentError(
            "iaeolian=1 requires non-empty w10 wind speed in BoundaryTimeSeries"))
        c.options.iprofl == 1 || @warn(
            "iaeolian=1 with iprofl=$(c.options.iprofl) — bed will not evolve from aeolian flux unless iprofl=1")
        c.aeolian.deposition_shape in (:exponential, :triangular, :right_triangle, :gaussian) ||
            throw(ArgumentError(
                "AeolianConfig.deposition_shape must be one of (:exponential, :triangular, :right_triangle, :gaussian); got :$(c.aeolian.deposition_shape)"))
        c.aeolian.deposition_shape == :gaussian && c.aeolian.gaussian_smooth_sigma > 0 ||
            c.aeolian.deposition_shape != :gaussian || throw(ArgumentError(
                "AeolianConfig.gaussian_smooth_sigma must be > 0 when deposition_shape=:gaussian"))
        c.aeolian.veg_deposition_center >= 0.0 || throw(ArgumentError(
            "AeolianConfig.veg_deposition_center must be >= 0 (got $(c.aeolian.veg_deposition_center)); negative offsets aren't supported yet"))
    end
    if c.aeolian !== nothing && c.options.iaeolian == 0
        @warn "AeolianConfig provided but iaeolian=0 — wind transport will NOT be applied. Set iaeolian=1 in OptionFlags to activate."
    end
    c.options.iwindshear in (0, 1) || throw(ArgumentError(
        "iwindshear must be 0 or 1 (got $(c.options.iwindshear))"))
    if c.options.iwindshear == 1
        c.windshear !== nothing || throw(ArgumentError(
            "iwindshear=1 requires a WindShearConfig — pass windshear=WindShearConfig()"))
        c.options.iaeolian == 1 || @warn(
            "iwindshear=1 with iaeolian=0 — shear perturbation will be computed but never used (no aeolian transport active)")
        c.windshear.shear_method in (:duna_kroy, :aeolis_kroy) || throw(ArgumentError(
            "WindShearConfig.shear_method must be :duna_kroy or :aeolis_kroy (got :$(c.windshear.shear_method))"))
        if c.windshear.shear_method == :aeolis_kroy
            c.windshear.aeolis_length_scale > 0 || throw(ArgumentError(
                "aeolis_length_scale must be > 0 when shear_method=:aeolis_kroy"))
        end
        c.windshear.kroy_solver in (:direct, :fft) || throw(ArgumentError(
            "WindShearConfig.kroy_solver must be :direct or :fft (got :$(c.windshear.kroy_solver))"))
        c.windshear.tau_smooth_sigma >= 0 || throw(ArgumentError(
            "tau_smooth_sigma must be ≥ 0 (got $(c.windshear.tau_smooth_sigma))"))
    end
    # ---- scaffolding flags (full implementation pending) ----
    c.options.iwave_aeolian_coupling in (0, 1) || throw(ArgumentError(
        "iwave_aeolian_coupling must be 0 or 1 (got $(c.options.iwave_aeolian_coupling))"))
    c.options.iwave_aeolian_coupling == 1 && @warn(
        "iwave_aeolian_coupling=1 is a scaffolding flag — full implementation pending. " *
        "Set to 0 for current behaviour.")
    c.options.iwind_threshold_extras in (0, 1) || throw(ArgumentError(
        "iwind_threshold_extras must be 0 or 1 (got $(c.options.iwind_threshold_extras))"))
    c.options.iwind_threshold_extras == 1 && @warn(
        "iwind_threshold_extras=1 is a scaffolding flag — full implementation pending. " *
        "Set to 0 for current behaviour.")
    c.options.iveg_dynamics in (0, 1) || throw(ArgumentError(
        "iveg_dynamics must be 0 or 1 (got $(c.options.iveg_dynamics))"))
    c.options.iveg_dynamics == 1 && @warn(
        "iveg_dynamics=1 is a scaffolding flag — full implementation pending. " *
        "Set to 0 for current behaviour.")
    c.options.iv_transport in (0, 1) || throw(ArgumentError(
        "iv_transport must be 0 or 1 (got $(c.options.iv_transport))"))
    if c.options.iv_transport == 1
        # With icurrent=1 (imposed current), vmean is propagated even at iangle=0,
        # so iv_transport=1 DOES have an effect in that case.
        if c.options.iangle == 0 && c.options.icurrent == 0
            @warn("iv_transport=1 has no effect with iangle=0 and icurrent=0 " *
                  "(vmean is not solved without oblique waves or imposed current). " *
                  "Set iangle=1 for wave-driven currents, or icurrent=1 for imposed current.")
        end
        # All transport formulas honor iv_transport=1:
        #   :original        — ustd_eff = sqrt(ustd²+vmean²) augments bed shear
        #   :size_adaptive   — routes sand to SVR kernel which reads vmean
        #   :soulsby_vanrijn — SVR kernel reads vmean directly
    end
    if c.options.icurrent == 1
        c.current !== nothing || throw(ArgumentError(
            "icurrent=1 requires a CurrentInput — pass current=CurrentInput(time=..., vbc=...)"))
        length(c.current.time) == length(c.current.vbc) || throw(ArgumentError(
            "current.time and current.vbc must have equal length"))
        length(c.current.time) >= 1 || throw(ArgumentError(
            "current.time must have at least 1 element"))
        if c.options.iangle == 0
            @warn("icurrent=1 with iangle=0: imposed current will affect sediment " *
                  "transport (via vmean when iv_transport=1) but NOT the longshore " *
                  "momentum balance (Sxy radiation stress requires oblique waves). " *
                  "Set iangle=1 for full wave-driven current dynamics.")
        end
    end
    if c.current !== nothing && c.options.icurrent == 0
        @warn "CurrentInput provided but icurrent=0 — imposed alongshore current will NOT be applied. Set icurrent=1 in OptionFlags to activate."
    end
    c.options.iweibull == 1 && throw(ArgumentError("iweibull=1 not supported"))
    c.options.idiss  == 3 && throw(ArgumentError("idiss=3 (measured spectrum) not supported"))
    c.options.iqydy == 1 && throw(ArgumentError(
        "iqydy=1 (alongshore transport coupling) not supported in this port"))

    return nothing
end

# Physical constants
const GRAV  = 9.81
const SQR2  = sqrt(2.0)
const SQR8  = sqrt(8.0)
const PI2   = 2π
const SQRG1 = sqrt(GRAV)
const SQRG2 = GRAV * sqrt(GRAV)
