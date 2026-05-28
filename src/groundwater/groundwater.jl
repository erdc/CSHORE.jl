# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
groundwater.jl — Coupled beach groundwater and surface-moisture model.

Translates the core physics of Psamathe (Utrecht Coastal Group; Brakenhoff
et al. 2019, Coastal Engineering) into CSHORE.jl's time-averaged framework,
replacing Psamathe's Stockdon (2006) runup parameterisation with CSHORE's
native wave-setup and BDJ swash fields.

Entry point: step_groundwater!(state, config, wave, l, dt_window)

Called once per BC window from driver.jl immediately after the wave+swash
solve.  Populates:
  state.gw_eta[j, l]  — water table elevation at each node (m above datum)
  state.theta[j, l]   — volumetric surface moisture (-), from Van Genuchten

Physics layers
--------------
1. Nonlinear 1D Boussinesq groundwater equation (cross-shore, transient):

     nₑ ∂η/∂t = K ∂/∂x [(D + η) ∂η/∂x] + q_src(x)

   where:
     η(x,t)   = water table elevation above datum (m)
     K         = hydraulic conductivity (m/s)
     nₑ        = effective (drainable) porosity (-)
     D         = aquifer thickness below datum (m)
     q_src(x)  = swash infiltration source (m/s); non-zero over jwd:jdry

   Boundary conditions:
     Seaward (j ≤ jswl):   η[j] = swl + wsetup[j]   (ocean-controlled; not marched)
     Seaward face (j=jbc):  Dirichlet: η = swl + wsetup[jbc]
     Landward (j = jmax):   dη/dx = 0 (no-flux, default)
                          or η = gw_eta_landward (fixed WT override)

   Solved with explicit forward-Euler finite differences, sub-stepping at
   dt_gw (default 30 s) within each BC window.  CFL stability requires:
     dt_gw ≤ nₑ dx² / (2 K (D + max(η)))
   where max(η) is the current maximum water-table elevation in the domain.
   The solver computes this bound each BC window and auto-corrects dt_sub
   when gw.dt_gw would exceed it.

2. Van Genuchten soil-water retention curve (unsaturated zone):

     θ(h) = θ_res + (θ_sat − θ_res) / [1 + (α h)ⁿ]ᵐ,   m = 1 − 1/n

   where h = max(zb − η, 0) is the depth from surface to water table.
   When η ≥ zb (table at or above surface), θ = θ_sat (fully saturated).
   Defaults are calibrated for medium beach sand (Carsel & Parrish 1988).

3. Infiltration source / INFILT coupling:
   When infilt=1, the per-node swash infiltration source rate uses wpm
   (Forchheimer-derived max seepage velocity, stored in state.wpm_derived
   by wetdry!).  This equates the swash momentum drain (INFILT side) with
   the aquifer recharge source (groundwater side), conserving mass.

   When GroundwaterConfig.rainfall_rate is set (m/s), a uniform vertical
   recharge flux is added at every dry-backshore node (landward of the swash
   runup limit).  During active rain, state.theta[j, l] is set to theta_sat
   (or rainfall_wet_theta if specified), overriding Van Genuchten — rainfall
   fully wets the surface and suppresses aeolian transport.

4. Aeolian coupling:
   state.theta[j, l] is consumed by step_aeolian! as the local volumetric
   moisture M for the Belly-Johnson / Kawamura threshold modification.
   It overrides the empirical ae_moisture decay when GroundwaterConfig is
   present.

References
----------
  Brakenhoff et al. (2019) Simulating surface soil moisture on sandy beaches
    using a 1D groundwater model, Coastal Engineering 152.
  Van Genuchten (1980) A closed-form equation for predicting the hydraulic
    conductivity of unsaturated soils, SSSA 44(5).
  Carsel & Parrish (1988) Developing joint probability distributions of soil
    water retention characteristics, Water Resources Research 24(5).
==============================================================================#

"""
    step_groundwater!(state, config, wave, l, dt_window)

Advance the beach groundwater model by `dt_window` seconds for cross-shore
line `l`.  Updates `state.gw_eta[:, l]` and `state.theta[:, l]` in-place.

The function is a no-op when `config.groundwater === nothing`.

# Arguments
- `state`      — mutable model state
- `config`     — immutable configuration (reads `config.groundwater`)
- `wave`       — current BC-window wave parameters (swl, hrms0, tp)
- `l`          — cross-shore line index (1-based)
- `dt_window`  — BC window duration (s); the solver sub-steps internally
"""
function step_groundwater!(state::CshoreState, config::CshoreConfig,
                           wave::WaveParams, l::Int, dt_window::Float64)
    config.groundwater === nothing && return nothing

    gw     = config.groundwater
    dx     = config.grid.dx
    jmax_l = state.jmax[l]
    jswl   = state.jswl[l]
    jr     = state.jr
    jwd    = state.jwd
    jdry   = state.jdry

    # ── Seaward Dirichlet BC: setup-elevated tidal water level ───────────────
    # Use the last valid setup node (min of jswl and jr to avoid reading
    # zero-setup values in the dry zone).  wsetup is the wave-averaged
    # setup above SWL, already solved by transform_waves!.
    j_bc = max(1, min(jswl, jr))
    eta_bc = wave.swl + state.wsetup[j_bc]

    # ── Initialise gw_eta on cold start ──────────────────────────────────────
    # On the first call (gw_eta is all zeros), seed with ocean water level
    # so the model doesn't start from an unrealistic dry state.
    if all(iszero, view(state.gw_eta, 1:jmax_l, l))
        @inbounds for j in 1:jmax_l
            state.gw_eta[j, l] = wave.swl + state.wsetup[min(j, jr)]
        end
    end

    # ── Swash infiltration source (m/s) ──────────────────────────────────────
    # Resolve the per-node source rate.  Three cases:
    #   (a) explicit override in GroundwaterConfig — user is responsible
    #   (b) infilt=1 — use wpm (Forchheimer rate from wetdry!), same flux
    #       that drains the swash momentum balance → mass conservation
    #   (c) neither — no swash source term; forcing comes via seaward BC only
    q_infilt_rate = if gw.infiltration_rate !== nothing
        gw.infiltration_rate
    elseif config.options.infilt == 1
        state.wpm_derived          # set by wetdry! each BC window
    else
        0.0
    end

    # ── Rainfall recharge (m/s) ───────────────────────────────────────────────
    # Uniform rainfall applied over the dry backshore (j > jdry or j > jswl
    # when no swash zone exists).  Positive recharge → raises the water table.
    # The rain rate is independent of the swash infiltration pathway: it
    # operates at the surface and enters the aquifer as vertical infiltration.
    # CFL: rainfall rate is typically O(10⁻⁷–10⁻⁵) m/s, far below the
    # hydraulic-diffusivity driven CFL, so it does not constrain dt_gw.
    #
    # Scalar: constant for all BC windows.
    # Vector: linearly interpolated from the boundary time series at state.time.
    t_now = state.time
    q_rain = if gw.rainfall_rate === nothing
        0.0
    elseif gw.rainfall_rate isa Float64
        gw.rainfall_rate
    else  # Vector — interpolate at current simulation time
        interp1(config.boundary.timebc, gw.rainfall_rate, t_now)
    end

    # Resolve the surface-moisture value to assign during rain (θ_rain).
    # Scalar or time-interpolated vector; defaults to theta_sat when not set.
    theta_rain_cfg = if gw.rainfall_wet_theta === nothing
        gw.theta_sat
    elseif gw.rainfall_wet_theta isa Float64
        gw.rainfall_wet_theta
    else
        interp1(config.boundary.timebc, gw.rainfall_wet_theta, t_now)
    end

    # ── Build per-node source vector ──────────────────────────────────────────
    # gw_source[j] (m/s): vertical recharge flux entering node j.
    # The Boussinesq solver divides by nₑ to convert flux → rate-of-table-rise.
    # Two source types share the vector (additive):
    #   (1) Swash zone infiltration (jwd:jdry)
    #   (2) Rainfall recharge      (jdry+1:jmax_l, i.e. dry backshore)
    gw_source = zeros(Float64, jmax_l)   # allocation is cheap; jmax_l ~ O(100–1000)

    # Swash infiltration over the wet-dry front
    if q_infilt_rate > 0.0 && jdry > jr && jwd > 0
        @inbounds for j in jwd:min(jdry, jmax_l)
            # pwet[j] is the time-averaged wet probability at each swash node.
            # Clamp to avoid negative source from numerical noise.
            gw_source[j] = q_infilt_rate * max(state.pwet[j], 0.0)
        end
    end

    # Rainfall recharge over the dry backshore (landward of swash runup limit)
    if q_rain > 0.0
        j_rain_start = max(j_bc + 1, jdry + 1)   # start landward of swash
        @inbounds for j in j_rain_start:jmax_l
            gw_source[j] += q_rain   # additive in case swash zone overlaps
        end
    end

    eta = view(state.gw_eta, :, l)   # working alias — modified in-place

    # ── Boussinesq sub-stepping ───────────────────────────────────────────────
    # Explicit forward-Euler CFL stability:
    #   dt ≤ nₑ dx² / (2 K (D + max(η)))
    #
    # IMPORTANT: the transmissivity is T = K·(D + η), NOT K·D.  When the water
    # table is elevated (e.g. η = 2.5 m at IC) the actual CFL limit is
    # significantly tighter than the D-only estimate.  Using the D-only bound
    # under-counts sub-steps and produces a numerically unstable (oscillating,
    # drifting-negative) solution.  We therefore compute the actual CFL from the
    # current maximum WT and use that to set n_sub, overriding gw.dt_gw when
    # needed to guarantee stability.
    eta_max_domain = j_bc < jmax_l ?
        maximum(view(eta, j_bc:jmax_l)) : eta[j_bc]
    T_max    = gw.K * max(gw.D + eta_max_domain, gw.D)   # K·max(D+η)
    dt_cfl   = gw.ne * dx^2 / (2.0 * T_max)              # actual CFL limit (s)
    dt_safe  = min(gw.dt_gw, dt_cfl * 0.95)              # 5% safety margin
    n_sub    = max(1, ceil(Int, dt_window / dt_safe))
    dt_sub   = dt_window / n_sub

    # Warn if user-specified dt_gw exceeded the actual CFL (once per session is
    # sufficient; the code corrects automatically so this is advisory only).
    if gw.dt_gw > dt_cfl * 1.05
        @warn "Groundwater: dt_gw=$(gw.dt_gw) s exceeds actual Boussinesq CFL " *
              "($(round(dt_cfl, sigdigits=3)) s for K=$(gw.K), ne=$(gw.ne), " *
              "D=$(gw.D), dx=$(dx), max_eta=$(round(eta_max_domain,digits=2)) m). " *
              "Auto-corrected to $(round(dt_sub,sigdigits=3)) s ($(n_sub) sub-steps)." maxlog=3
    end

    for _ in 1:n_sub
        _boussinesq_step!(eta, state.zb, l, j_bc, jmax_l, eta_bc, gw_source, gw, dx, dt_sub)
    end

    # ── Van Genuchten surface moisture ────────────────────────────────────────
    # theta_rain_cfg was resolved above (scalar or time-interpolated).
    # During active rain, override Van Genuchten with the wet-surface value.
    j_rain_start = max(j_bc + 1, jdry + 1)   # same threshold as recharge loop
    rain_active  = q_rain > 0.0

    @inbounds for j in 1:jmax_l
        if rain_active && j >= j_rain_start
            # Active rainfall: surface is wet regardless of water-table depth.
            # Suppresses aeolian transport; overrides Van Genuchten.
            state.theta[j, l] = theta_rain_cfg
        else
            h_to_table = state.zb[j, l] - eta[j]   # depth below surface to WT (m)
            state.theta[j, l] = _van_genuchten(h_to_table, gw)
        end
    end
    # Zero beyond the active profile
    @inbounds for j in (jmax_l + 1):size(state.theta, 1)
        state.theta[j, l] = 0.0
    end

    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal: one explicit Boussinesq finite-difference step
# ─────────────────────────────────────────────────────────────────────────────

"""
    _boussinesq_step!(eta, zb, l, j_bc, jmax_l, eta_bc, gw_source, gw, dx, dt)

Advance the water-table elevation vector `eta` by one explicit time step `dt`.

Domain: j = j_bc+1 … jmax_l (the dry backshore).
  j ≤ j_bc  — ocean-controlled; forced to `eta_bc` each sub-step.
  j = jmax_l — no-flux landward BC (or fixed elevation from gw.gw_eta_landward).

Transmissivity at the interface between nodes j and j+1:
  T_{j+½} = K · (D + ½(η_j + η_{j+1}))     [clamped ≥ 0]

Flux at interface:
  F_{j+½} = T_{j+½} · (η_{j+1} − η_j) / dx

Update:
  η_j^{n+1} = η_j^n + (dt/(nₑ·dx)) · (F_{j+½} − F_{j−½}) + dt·q_src_j/nₑ
"""
function _boussinesq_step!(eta::AbstractVector{Float64},
                           zb::AbstractMatrix{Float64}, l::Int,
                           j_bc::Int, jmax_l::Int,
                           eta_bc::Float64,
                           gw_source::Vector{Float64},
                           gw::GroundwaterConfig,
                           dx::Float64, dt::Float64)
    K  = gw.K
    ne = gw.ne
    D  = gw.D
    ne_inv = 1.0 / ne
    dt_ne   = dt * ne_inv
    dt_ne_dx = dt_ne / dx

    # Apply seaward Dirichlet BC: ocean nodes stay at eta_bc
    @inbounds for j in 1:j_bc
        eta[j] = eta_bc
    end

    # Apply landward BC: fixed WT override or leave (no-flux handled implicitly)
    if isfinite(gw.gw_eta_landward)
        eta[jmax_l] = gw.gw_eta_landward
    end

    # Build flux vector: F[j] = flux at the j+½ interface (between j and j+1)
    # F[j_bc] = 0 (left wall) is handled by not updating j_bc.
    # Reuse eta scratch to avoid allocation: compute fluxes in a temp buffer.
    # For in-place efficiency we update eta directly from j_bc+1 to jmax_l-1
    # using staggered fluxes; store F in a local vector.
    n_domain = jmax_l - j_bc   # number of interior nodes to update
    n_domain <= 0 && return nothing

    # Compute interface transmissivities and fluxes
    # F[i] is the flux at the right face of node (j_bc + i), i ∈ 1:n_domain
    # We need F at i-1 and F at i for each node update.
    # Compute inline to avoid extra allocation: loop once, updating in a second pass.

    # Step 1: compute fluxes F_{j+½} for j = j_bc … jmax_l-1
    # (stored as F_right[j] for node j, i.e., the interface on the right of j)
    F_right = Vector{Float64}(undef, jmax_l)

    @inbounds for j in j_bc:(jmax_l - 1)
        T_half = K * max(D + 0.5 * (eta[j] + eta[j + 1]), 0.0)
        F_right[j] = T_half * (eta[j + 1] - eta[j]) / dx
    end
    # Zero flux at the left boundary (no incoming flux from seaward ocean column)
    F_right[j_bc - 1 < 1 ? 1 : j_bc - 1] = 0.0   # left wall
    # Landward BC: no-flux → F at j = jmax_l-½ is already computed; no adjustment
    # needed unless fixed WT is set (then eta[jmax_l] is fixed, flux naturally adjusts).

    # Step 2: update eta for j = j_bc+1 … jmax_l-1 (interior)
    @inbounds for j in (j_bc + 1):(jmax_l - 1)
        F_left  = F_right[j - 1]
        F_right_j = F_right[j]
        src_j   = j <= length(gw_source) ? gw_source[j] : 0.0
        eta[j] += dt_ne_dx * (F_right_j - F_left) + dt * src_j * ne_inv
        # Clamp: water table cannot rise above bed surface in the dry backshore
        # (ponding is handled by the wetdry!/overtopping modules).
        # GUARD: only apply the cap where the bed is still elevated above the
        # seaward BC.  If morphological erosion has lowered zb[j] below eta_bc,
        # the node is effectively wet/submerged; forcing eta = zb would create
        # an artificial WT sink that the diffusion cannot recover from.
        # For those eroded nodes we let Boussinesq diffusion equilibrate the WT
        # naturally between the seaward BC and the landward reservoir.
        zb_j = zb[j, l]
        if zb_j >= eta_bc
            eta[j] = min(eta[j], zb_j)
        end
    end

    # Landward boundary node (jmax_l): no-flux → F_{jmax_l+½} = F_{jmax_l-½}
    # This makes the net flux zero at jmax_l, equivalent to a mirror condition.
    if !isfinite(gw.gw_eta_landward)
        j = jmax_l
        F_left = F_right[j - 1]
        # No-flux: F_right[j] = 0  (closed boundary)
        src_j  = j <= length(gw_source) ? gw_source[j] : 0.0
        eta[j] += dt_ne_dx * (0.0 - F_left) + dt * src_j * ne_inv
        zb_j = zb[j, l]
        if zb_j >= eta_bc
            eta[j] = min(eta[j], zb_j)
        end
    end

    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Van Genuchten soil-water retention curve
# ─────────────────────────────────────────────────────────────────────────────

"""
    _van_genuchten(h, gw) -> θ

Volumetric soil moisture at depth `h` (m) above the water table (positive
upward, so h = zb − η).

- h ≤ 0: water table at/above surface → fully saturated: θ = θ_sat
- h > 0: unsaturated zone; moisture decreases with depth-to-table

Uses the Van Genuchten (1980) retention curve:
  θ = θ_res + (θ_sat − θ_res) / [1 + (α h)ⁿ]ᵐ,   m = 1 − 1/n
"""
@inline function _van_genuchten(h::Float64, gw::GroundwaterConfig)::Float64
    h <= 0.0 && return gw.theta_sat     # saturated (table at/above surface)
    m = 1.0 - 1.0 / gw.vg_n
    denom = (1.0 + (gw.vg_alpha * h)^gw.vg_n)^m
    return gw.theta_res + (gw.theta_sat - gw.theta_res) / denom
end
