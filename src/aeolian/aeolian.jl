# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
aeolian.jl — Wind-driven sediment transport (IAEOLIAN=1).

Kawamura (1951) transport with Delgado-Fernández (2011) fetch limitation,
extended for multifraction sand and contour-based vegetation capture.

- Operates per cross-shore line, on the dry beach above SWL + runup_buffer.
- Computes the friction velocity u* from a 10 m wind speed via a log-law
  with z0 = z0_factor · D50 (Nikuradse-style roughness).
- Applies a per-fraction Bagnold/Kawamura threshold u*_t,k from each grain
  diameter d_k (so fine sand mobilises before coarse, giving armouring
  for free without any explicit armouring term).
- Sums per-fraction Kawamura flux  Q_k = Ck · ρ_a/g · (u* - u*_t,k) · (u* + u*_t,k)²
  weighted by surface availability comp_k (active-layer fractions).
- Applies fetch limitation: F = beach_width / |cos θ_w|, multiplied by
  sin(π F / 2 Fc) when F < Fc = a · u_w + b. Cross-shore projection by cos.
- Vegetation capture is delegated to an `AeolianVegetationModel`. `ContourVegetation`
  deposits 100% above z_contour with optional exponential ramp-up.

Bed change is applied to `state.zb` directly. Losses exiting the domain
landward / seaward are accumulated in `state.ae_loss_landward / _seaward`
for mass-balance bookkeeping.
==============================================================================#

"""
    capture_efficiency(model::AeolianVegetationModel, j::Int, state, l::Int) -> η

Return the per-cell aeolian capture efficiency `η ∈ [0, 1]`. Cell `j` retains
fraction `η` of the incoming wind-blown flux; the remainder passes downwind.
"""
function capture_efficiency end

@inline function capture_efficiency(t::ContourVegetation, j::Int, state::CshoreState, l::Int)
    isnan(t.z_contour) && return 0.0          # no contour configured → flux passes through
    z = state.zb[j, l]
    if z < t.z_contour
        return 0.0
    end
    if t.decay_length <= 0.0
        return 1.0                             # sharp Heaviside
    end
    return 1.0 - exp(-(z - t.z_contour) / t.decay_length)
end

@inline function capture_efficiency(t::DensityVegetation, j::Int,
                                 state::CshoreState, l::Int)
    rho = j <= length(t.vegrho) ? t.vegrho[j] : 0.0
    rho <= 0.0 && return 0.0
    return 1.0 - exp(-t.alpha * rho)
end

@inline function capture_efficiency(t::MultiSpeciesVegetation, j::Int,
                                  state::CshoreState, l::Int)
    nrows = size(t.biomass, 1)
    j <= nrows || return 0.0
    nspecies = size(t.biomass, 2)
    nspecies == length(t.species) || return 0.0
    # Aggregate Raupach term across species: λ = Σ_s f_s · α_s · cover_s
    lambda = 0.0
    @inbounds for s in 1:nspecies
        sp = t.species[s]
        B = t.biomass[j, s]
        B_eq = t.biomass_eq[j, s]
        (B <= 0.0 || B_eq <= 0.0) && continue
        cover_s = 1.0 - exp(-sp.cover_per_biomass * B / B_eq)
        lambda += sp.frontal_area_factor * sp.alpha * cover_s
    end
    return lambda <= 0.0 ? 0.0 : (1.0 - exp(-lambda))
end

"""
    _ust_threshold_kawamura(d_k, moisture)

Belly-Johnson moisture-modified Bagnold threshold. Uses the form

    u*_t = 0.1 · √( g · D · (ρ_s/ρ_a) · (1 + C·M) )

with C = 1.87 and dry-sand specific gravity 2.65 (so ρ_s/ρ_a ≈ 2160).
M ∈ [0, 1] is the volumetric surface-moisture content. M = 0 (default)
recovers the canonical dry threshold; higher M steepens the threshold,
suppressing transport on damp sand.
"""
@inline function _ust_threshold_kawamura(d_k::Float64, moisture::Float64)
    rho_s = 2650.0
    rho_a = 1.225
    return 0.1 * sqrt(GRAV * d_k * (rho_s / rho_a) * (1.0 + 1.87 * moisture))
end

"""
    _surface_d50(state, j, k_range, mf)

Surface-fraction-weighted D50 used to set the aerodynamic roughness z0.
With single-grain sand this is just the grain diameter; with multifraction
it reflects the active-layer composition.
"""
@inline function _surface_d50(state::CshoreState, j::Int, mf::MultifractionConfig)
    nf = nfractions(mf)
    d_k_sum = 0.0
    w_sum   = 0.0
    @inbounds for k in 1:nf
        w = state.active_frac[j, k]
        d_k_sum += w * mf.grain_sizes[k]
        w_sum   += w
    end
    return w_sum > 0 ? d_k_sum / w_sum : mf.grain_sizes[1]
end

"""
    step_aeolian!(state, config, itime; dt::Float64)

Apply one aeolian sub-step over duration `dt` (seconds). Reads the wind
record (BoundaryTimeSeries.w10, .wangle) at the current `state.time`,
computes per-fraction Kawamura flux on the dry beach, applies vegetation
capture according to `config.aeolian.vegetation`, and updates `state.zb`. Tracks
mass loss to the domain edges.

Wind-direction convention:
  `wangle` (rad) is the wind direction in the same sense as `wangbc`
  (used for wind stress on water surface): the angle measured from
  shore-normal toward +y, in the direction the wind is blowing TO.
  cos(wangle) > 0 ⇒ onshore (toward landward), cos(wangle) < 0 ⇒ offshore.

Skips silently if IAEOLIAN=0 or no AeolianConfig is provided.
"""
function step_aeolian!(state::CshoreState, config::CshoreConfig,
                       itime::Int; dt::Float64)
    config.options.iaeolian == 1 || return nothing
    cfg = config.aeolian
    cfg === nothing && return nothing
    dt > 0.0 || return nothing

    bc      = config.boundary
    mf      = config.multifraction
    nf      = nfractions(mf)
    dx      = config.grid.dx
    rho_a   = cfg.rho_air
    rho_s   = config.sediment.sg * 1000.0      # kg/m³
    sporo1  = config.sediment.sporo1           # 1 - porosity
    g       = GRAV

    # Wind at current time
    u10  = _aeolian_interp(bc.w10,    bc.timebc, state.time)
    wdir = _aeolian_interp(bc.wangle, bc.timebc, state.time)

    cos_w = cos(wdir)
    abs_cos_w = abs(cos_w)
    if abs_cos_w < 1e-3 || u10 <= 0
        # Pure alongshore or calm — no cross-shore aeolian transport
        return nothing
    end
    onshore = cos_w > 0     # wind blowing in +x direction (toward landward)

    for l in 1:config.options.iline
        jmax = state.jmax[l]

        # Fill zero by default (so diagnostic reads make sense)
        @inbounds for j in 1:jmax
            state.qae[j] = 0.0
            state.ae_dz[j] = 0.0
            for k in 1:nf
                state.qae_k[j, k] = 0.0
            end
        end

        # ------------------------------------------------------------
        # Surface moisture book-keeping (when cfg.dry_time > 0).
        # Mark cells currently wet (below SWL or inside the swash band)
        # by stamping their last-wet timestamp with state.time. Cells
        # above the runup line keep their previous timestamp; their
        # diagnostic moisture decays linearly from cfg.moisture toward
        # cfg.moisture_dry over cfg.dry_time seconds. The actual
        # threshold-velocity modifier is read off `state.ae_moisture`.
        swl_now = _aeolian_interp(bc.swlbc, bc.timebc, state.time)
        if cfg.dry_time > 0.0
            j_wet_landward = config.options.iover == 1 && state.jdry > 0 ?
                              state.jdry : (state.jswl[l] > 0 ? state.jswl[l] : 0)
            @inbounds for j in 1:jmax
                if state.zb[j, l] <= swl_now || (j_wet_landward > 0 && j <= j_wet_landward)
                    # Currently wet — reset the clock
                    state.ae_t_last_wet[j] = state.time
                    state.ae_moisture[j] = cfg.moisture
                else
                    # Dry — linear interpolation between the wet-time and
                    # `dry_time` later. Cells with `t_last_wet = -Inf`
                    # (never wetted) jump straight to moisture_dry.
                    if isfinite(state.ae_t_last_wet[j])
                        f = clamp((state.time - state.ae_t_last_wet[j]) / cfg.dry_time,
                                  0.0, 1.0)
                        state.ae_moisture[j] = cfg.moisture +
                                               f * (cfg.moisture_dry - cfg.moisture)
                    else
                        state.ae_moisture[j] = cfg.moisture_dry
                    end
                end
            end
        end

        # ------------------------------------------------------------
        # 1. Identify the dry-beach (active) zone vs the deposition zone
        # ------------------------------------------------------------
        # Seaward edge of the truly-dry beach is taken directly from
        # CSHORE's own tracker — `state.jdry` when IOVER=1 (post-swash
        # march), `state.jswl` otherwise. See `_aeolian_runup_node` for
        # full source priority.
        swl = _aeolian_interp(bc.swlbc, bc.timebc, state.time)
        i_runup = _aeolian_runup_node(state, config, l, swl)
        i_runup == 0 && continue   # no usable dry-beach line

        # Vegetation contour: first cell where the vegetation model gives η > 0.
        # For ContourVegetation this is z >= z_contour. Marks the seaward edge of
        # the vegetated zone (and therefore the LANDWARD edge of the bare
        # erodible beach).
        i_contour_start = jmax + 1   # default: no contour → no vegetation
        @inbounds for j in i_runup:jmax
            if capture_efficiency(cfg.vegetation, j, state, l) > 0.0
                i_contour_start = j; break
            end
        end

        # Deposition zone start. Equals i_contour_start by default. With
        # `veg_deposition_center > 0` the deposition footprint is shifted
        # landward of the vegetation contour by that many meters: cells in
        # [i_contour_start, i_dep_start) become a vegetation-protected
        # bypass zone (no entrainment, no deposition). Useful when the
        # seaward vegetation edge marks where wind starts to slow but the
        # actual sand capture sits further inland.
        i_dep_start = i_contour_start
        if i_contour_start <= jmax && cfg.veg_deposition_center > 0.0
            x_dep_target = state.xb[i_contour_start] + cfg.veg_deposition_center
            found = jmax + 1
            @inbounds for j in i_contour_start:jmax
                if state.xb[j] >= x_dep_target
                    found = j; break
                end
            end
            i_dep_start = found
        end

        active_lo = i_runup
        active_hi = min(i_contour_start - 1, jmax)
        if active_lo > active_hi
            # Vegetation contour reaches to the runup line — no fetch / bare beach
            continue
        end

        # ------------------------------------------------------------
        # 1b. Optional wind-flow-over-topography update (IWINDSHEAR=1)
        # ------------------------------------------------------------
        # Builds state.tau_perturbation[j] (= τ/τ₀ at each cell) and
        # state.lee_zone[j] (true inside lee separation bubbles). The
        # per-fraction Kawamura entrainment below multiplies the local
        # u* by √(τ/τ₀) so cells over the dune crest see speed-up and
        # lee-shadowed cells see no transport. With IWINDSHEAR=0 the
        # tau_perturbation array stays at 1.0 → uniform u*.
        compute_wind_shear!(state, config, l, swl;
                             wind_direction = onshore ? :onshore : :offshore)

        # ------------------------------------------------------------
        # 2. Fetch & cross-shore projection
        # ------------------------------------------------------------
        beach_width = (active_hi - active_lo + 1) * dx
        # Effective fetch over erodible bed (oblique winds traverse a longer
        # path before crossing the active zone).
        fetch = beach_width / abs_cos_w
        Fc = cfg.fetch_critical_a * u10 + cfg.fetch_critical_b
        fetch_factor = if Fc <= 0
            1.0
        elseif fetch >= Fc
            1.0
        else
            sin(0.5 * π * fetch / Fc)
        end
        # Cosine projection of the flux onto the cross-shore direction
        proj = abs_cos_w

        # ------------------------------------------------------------
        # 3. Per-fraction Kawamura flux at every active-beach cell
        # ------------------------------------------------------------
        # Per-cell roughness from surface D50; per-fraction threshold from d_k.
        # Q_k(j) is the equilibrium-capacity sand-mass flux at cell j (kg/m/s).
        @inbounds for j in active_lo:active_hi
            d50_local = _surface_d50(state, j, mf)
            z0 = max(cfg.z0_factor * d50_local, 1e-6)
            ust_uniform = u10 * cfg.karman / log(cfg.z_meas / z0)
            # Apply per-cell shear-stress perturbation (IWINDSHEAR=1).
            # τ_pert is 1.0 everywhere when IWINDSHEAR=0. In lee-separation
            # cells the perturbation is 0 → ust = 0 → no transport.
            ust = ust_uniform * sqrt(state.tau_perturbation[j])
            state.ust[j] = ust
            # Per-cell moisture: three-way priority —
            #   1. Physics-based Van Genuchten θ from groundwater model (authoritative)
            #   2. Empirical time-decay ae_moisture (fallback when no groundwater)
            #   3. Static cfg.moisture scalar (simplest fallback)
            M_local = if config.groundwater !== nothing
                state.theta[j, l]                       # physics-based (groundwater)
            elseif cfg.dry_time > 0.0
                state.ae_moisture[j]                    # empirical decay
            else
                cfg.moisture                            # static scalar
            end
            # Hardbottom BRF (sediment-thickness limiter). When ISEDAV ≠ 0
            # and the cell has a finite hardbottom (state.zb_hard > -Inf),
            # state.hp[j,l] gives the available sand depth above hardbottom.
            # We damp the per-fraction capacity by `(hp/d_k)^bedlm` so the
            # aeolian kernel cannot erode below the hardbottom. When hp >= d_k
            # the BRF is 1 (no reduction); when hp = 0 the cell is fully armored.
            isedav = config.options.isedav
            apply_brf = abs(isedav) ≥ 1
            hp_j = apply_brf ? state.hp[j, l] : Inf
            bedlm = config.sediment.bedlm

            # ---- Sheltering coverage (Raupach 1993) ----
            # Any fraction with u*_t,k > ust is "non-erodible" at this u*.
            # The sum of those fractions' surface coverage is λ — the area
            # density of non-erodible roughness elements that shelter the
            # erodible fractions. Only meaningful for nf > 1.
            shelter_lambda = 0.0
            if cfg.iuth_sheltering && nf > 1
                for k in 1:nf
                    d_k = mf.grain_sizes[k]
                    ust_t_k = _ust_threshold_kawamura(d_k, M_local)
                    if ust_t_k > ust
                        shelter_lambda += state.active_frac[j, k]
                    end
                end
            end
            shelter_factor = cfg.iuth_sheltering ?
                sqrt(1.0 + cfg.sheltering_msigma * shelter_lambda) : 1.0

            # ---- Bed-slope correction (Dyer 1986) ----
            # state.bslope[j, l] is dz/dx (positive = slope rising in +x =
            # landward). With onshore wind (going landward), positive slope
            # is uphill (harder to entrain → factor > 1); negative slope is
            # downhill (easier → factor < 1). For offshore wind, signs flip.
            slope_factor = 1.0
            if cfg.iuth_bedslope
                bs = state.bslope[j, l]
                tan_phi = config.sediment.tanphi
                # Slope angle (radians)
                α = atan(bs)
                cosα = cos(α)
                tanα = tan(α)
                # Sign convention: onshore + uphill → harder; offshore + uphill → easier
                slope_sign = onshore ? 1.0 : -1.0
                inner = cosα * (1.0 + slope_sign * tanα / tan_phi)
                # Clamp to avoid imaginary or unphysical values for very
                # steep downslopes (where the bed is at or beyond the
                # angle of repose).
                slope_factor = inner > 0.05 ? sqrt(inner) : 0.05
            end

            q_total = 0.0
            for k in 1:nf
                d_k = mf.grain_sizes[k]
                ust_t = _ust_threshold_kawamura(d_k, M_local) *
                        shelter_factor * slope_factor
                state.ust_threshold[j, k] = ust_t
                comp_k = state.active_frac[j, k]
                if comp_k <= 0 || ust <= ust_t
                    state.qae_k[j, k] = 0.0
                    continue
                end
                q_k_kg = cfg.Ck * (rho_a / g) *
                         (ust - ust_t) * (ust + ust_t)^2 *
                         comp_k * fetch_factor * proj      # kg/m/s
                # Convert to volume flux (m²/s sand grains, no pore space)
                q_k_vol = q_k_kg / rho_s
                # Hardbottom limiter
                if apply_brf && isfinite(hp_j)
                    brf = hp_j >= d_k ? 1.0 : (hp_j / d_k)^bedlm
                    q_k_vol *= brf
                end
                state.qae_k[j, k] = q_k_vol
                q_total += q_k_vol
            end
            state.qae[j] = q_total
        end

        # ------------------------------------------------------------
        # 4. Saturated-flux Exner march
        # ------------------------------------------------------------
        # Build a 1D advection-deposition transport equation. Track per-
        # fraction load q_k (m²/s sand-vol) along the wind direction,
        # cell by cell. At each cell:
        #
        #   • Active beach (bare, between runup and deposition zone):
        #         Q approaches local capacity Q_cap with saturation length
        #         L_sat:  Q' = Q + (Q_cap - Q) · (1 - exp(-Δx/L_sat))
        #     If Q < Q_cap: net erosion (the bed supplies the deficit).
        #     If Q ≥ Q_cap: net zero — the bypass condition.
        #
        #   • Deposition zone (above contour, with onshore wind):
        #         Q decays as  Q' = Q · exp(-η · Δx / L_dune)
        #     Mass lost from Q is deposited in this cell.
        #
        #   • Wet zone (below runup, with offshore wind):
        #         Q decays as  Q' = Q · exp(-Δx / L_dune)
        #     (the wave bed eats the deposit on the next BC window.)
        #
        # Bed change per cell:  Δz = -(Q_out - Q_in) · dt / (dx · sporo1).
        # Q_out > Q_in → erosion. Q_out < Q_in → deposition. Q_out = Q_in →
        # bypass (no net change).
        L_sat  = cfg.saturation_length
        L_dune = cfg.dune_decay_length
        inv_dx_sporo = 1.0 / (dx * sporo1)

        # Determine march direction & starting load.
        # Onshore wind: march from `active_lo` landward to `jmax`. Anything
        # leaving at j = jmax is recorded as landward loss.
        # Offshore wind: start at the most landward dry cell (or active_hi
        # if there's a deposition zone). March seaward, crossing the active
        # beach (entrainment), then into the wet zone (deposition), and
        # finally exit at j = 1 if anything is still airborne.
        if onshore
            j_first = max(1, active_lo)
            march_step = +1
            j_last = jmax
        else
            # Start at the seaward edge of the vegetated zone if it exists,
            # else at the landward boundary; this is the upwind boundary
            # for offshore wind. The veg-protected band between the contour
            # and the deposition zone is bypass for offshore wind too.
            j_first = i_contour_start <= jmax ? i_contour_start - 1 : jmax
            j_first = clamp(j_first, 1, jmax)
            march_step = -1
            j_last = 1
        end

        # Initialize per-fraction load
        @inbounds q_k = zeros(nf)
        # Latch for triangular kernels: applied once when the march first
        # enters the deposition zone (onshore wind only).
        triangle_applied = false

        @inbounds j = j_first
        while true
            in_active = active_lo <= j <= active_hi
            in_dep   = !in_active && (j >= i_dep_start)
            below_runup = j < active_lo
            # Cells in [i_contour_start, i_dep_start) are vegetation-protected
            # bypass — no entrainment, no deposition.

            # Triangular / right-triangle deposition kernels: deposit the entire
            # incoming saturated flux in a fixed-length footprint at the
            # seaward edge of the deposition zone. This is a one-shot
            # pre-deposition done before the per-fraction loop so the
            # per-fraction code can be a no-op in the deposition zone for
            # these shapes.
            if onshore && in_dep && !triangle_applied &&
               cfg.deposition_shape in (:triangular, :right_triangle, :gaussian)
                _apply_triangular_deposit!(state, cfg, l, j, dx, dt, sporo1,
                                            q_k, jmax; rho_s = rho_s)
                triangle_applied = true
            end

            # Mass-units factor for per-fraction bed_mass updates
            # (kg/m² per metre of bed elevation change).
            mass_per_m = sporo1 * rho_s

            for k in 1:nf
                if in_active
                    qcap = state.qae_k[j, k]
                    # Saturating ramp toward local capacity
                    q_new = q_k[k] + (qcap - q_k[k]) * (1.0 - exp(-dx / L_sat))
                    dz = -(q_new - q_k[k]) * dt * inv_dx_sporo  # neg = erosion
                    state.zb[j, l] += dz
                    state.ae_dz[j] += dz / dt
                    _update_bed_mass!(state, j, l, k, dz, mass_per_m)
                    q_k[k] = q_new
                elseif in_dep && onshore
                    if cfg.deposition_shape == :exponential
                        eta = capture_efficiency(cfg.vegetation, j, state, l)
                        q_new = q_k[k] * exp(-eta * dx / L_dune)
                        dz = -(q_new - q_k[k]) * dt * inv_dx_sporo  # pos = deposition
                        state.zb[j, l] += dz
                        state.ae_dz[j] += dz / dt
                        _update_bed_mass!(state, j, l, k, dz, mass_per_m)
                        q_k[k] = q_new
                    end
                    # :triangular / :right_triangle / :gaussian handled above
                    # by the one-shot _apply_triangular_deposit!; q_k is now
                    # zero, so this branch is a no-op for those shapes.
                elseif below_runup && !onshore
                    q_new = q_k[k] * exp(-dx / L_dune)
                    dz = -(q_new - q_k[k]) * dt * inv_dx_sporo
                    state.zb[j, l] += dz
                    state.ae_dz[j] += dz / dt
                    _update_bed_mass!(state, j, l, k, dz, mass_per_m)
                    q_k[k] = q_new
                else
                    # Bypass cell: no change to load, no bed change.
                end
            end

            # Total flux at this cell after the update (diagnostic)
            qtot = 0.0
            for k in 1:nf; qtot += q_k[k]; end
            state.qae[j] = qtot

            # Step
            if (march_step == +1 && j == j_last) ||
               (march_step == -1 && j == j_last)
                break
            end
            j += march_step
        end

        # Anything still airborne at the downwind edge exits the domain
        q_residual = 0.0
        @inbounds for k in 1:nf
            q_residual += q_k[k]
        end
        # Convert flux (m²/s) × dt → m³/m alongshore
        residual_vol = q_residual * dt
        if onshore
            state.ae_loss_landward += residual_vol
        else
            state.ae_loss_seaward  += residual_vol
        end

        # ---- Hardbottom clamp (revetment / bedrock) ----
        # Belt-and-suspenders: even with the BRF damper above, numerical
        # accumulation could nudge zb fractionally below zb_hard. Clamp
        # back to zb_hard and zero out any per-fraction surface mass for
        # cells that bottomed out. `zb_hard = -Inf` (the "no hardbottom"
        # sentinel) is a no-op.
        if abs(config.options.isedav) >= 1
            @inbounds for j in 1:jmax
                zh = state.zb_hard[j, l]
                isfinite(zh) || continue
                if state.zb[j, l] < zh
                    state.zb[j, l] = zh
                    for k in 1:nf
                        state.bed_mass[j, 1, k] = 0.0
                    end
                end
            end
        end
    end

    # Vegetation dynamics: regrowth toward equilibrium minus burial loss
    # proportional to the local aeolian deposition rate. Applied last so
    # the next BC window sees the updated cover.
    update_vegetation_density!(state, config, dt)

    return nothing
end

"""
    update_vegetation_density!(state, config, dt)

Evolve `DensityVegetation.vegrho` per cell over duration `dt` (s)
when `dynamics_enabled = true`. No-op for any other vegetation model or when
dynamics are disabled.

Per cell:

    dρ/dt = (ρ_eq - ρ) / τ_growth   −   max(0, dz_aeolian/dt) · k_burial

Forward-Euler step (the BC window is hours-to-days, much shorter than
τ_growth which is months-to-years, so explicit integration is stable
without sub-stepping). Result clamped to [0, 1].

The burial term uses `state.ae_dz[j]` — the aeolian deposition RATE
(m/s) at this cell, computed by `step_aeolian!`. Erosion (negative
ae_dz) doesn't directly remove vegetation; vegetation just regrows
toward `ρ_eq` once buried sand is gone.
"""
function update_vegetation_density!(state::CshoreState, config::CshoreConfig,
                                     dt::Float64)
    cfg = config.aeolian
    cfg === nothing && return nothing
    veg = cfg.vegetation
    if veg isa DensityVegetation
        _update_density_vegetation!(state, veg, dt)
    elseif veg isa MultiSpeciesVegetation
        _update_multispecies_vegetation!(state, veg, dt)
    end
    return nothing
end

@inline function _update_density_vegetation!(state::CshoreState,
                                             veg::DensityVegetation, dt::Float64)
    veg.dynamics_enabled || return nothing
    inv_tau = 1.0 / max(veg.tau_growth, 1.0)
    ρ_eq_arr = isempty(veg.vegrho_eq) ? veg.vegrho : veg.vegrho_eq
    n = length(veg.vegrho)
    @inbounds for j in 1:n
        ρ = veg.vegrho[j]
        ρeq = j <= length(ρ_eq_arr) ? ρ_eq_arr[j] : ρ
        burial_rate = j <= length(state.ae_dz) ? max(0.0, state.ae_dz[j]) : 0.0
        dρdt = (ρeq - ρ) * inv_tau - burial_rate * veg.k_burial
        ρ_new = ρ + dρdt * dt
        veg.vegrho[j] = clamp(ρ_new, 0.0, 1.0)
    end
    return nothing
end

@inline function _update_multispecies_vegetation!(state::CshoreState,
                                                  veg::MultiSpeciesVegetation,
                                                  dt::Float64)
    veg.dynamics_enabled || return nothing
    nrows = size(veg.biomass, 1)
    nspecies = size(veg.biomass, 2)
    nspecies == length(veg.species) || return nothing
    @inbounds for s in 1:nspecies
        sp = veg.species[s]
        inv_tau_s = 1.0 / max(sp.tau_growth, 1.0)
        for j in 1:nrows
            B = veg.biomass[j, s]
            B_eq = veg.biomass_eq[j, s]
            burial_rate = j <= length(state.ae_dz) ? max(0.0, state.ae_dz[j]) : 0.0
            # First-order regrowth toward B_eq − burial-driven loss.
            # Regrowth shuts off if B_eq is zero (cell is bare).
            dBdt = (B_eq - B) * inv_tau_s - burial_rate * sp.k_burial * B
            B_new = B + dBdt * dt
            # Clamp to [0, 2·B_eq] so a regrowth overshoot can't run away
            # if dt is large relative to τ.
            B_max = 2.0 * max(B_eq, 0.0)
            veg.biomass[j, s] = clamp(B_new, 0.0, B_max)
        end
    end
    return nothing
end

"""
    _apply_triangular_deposit!(state, cfg, l, j_seaward, dx, dt, sporo1,
                                q_k, jmax)

Applies the entire current per-fraction load `q_k` (m²/s sand-vol) as a
triangular bed-elevation footprint of total base length `L = cfg.dune_decay_length`
starting at cell `j_seaward`. `:triangular` places the peak at L/3 from the
seaward edge; `:right_triangle` places the peak at the seaward edge with
linear fall-off to L.

Mass-conserving: ∫ Δz(x) dx over [0, L] = (Σₖ q_k) · dt / sporo1, so the
total deposited bed-elevation-volume per unit alongshore length matches
the saturated incoming flux exactly. Per-fraction tracking: each fraction
`k` is deposited proportionally — same shape, scaled by q_k(in). After
the call, `q_k` is reset to zero (everything was deposited).

Cells beyond `jmax` are clamped (any residual mass that would have
landed past the landward boundary is added to `state.ae_loss_landward`).
"""
function _apply_triangular_deposit!(state::CshoreState, cfg::AeolianConfig,
                                    l::Int, j_seaward::Int,
                                    dx::Float64, dt::Float64, sporo1::Float64,
                                    q_k::Vector{Float64}, jmax::Int;
                                    rho_s::Float64 = 2650.0)
    nf = length(q_k)
    Q_in_total = 0.0
    @inbounds for k in 1:nf
        Q_in_total += q_k[k]
    end
    Q_in_total > 0 || (return nothing)
    # Per-fraction shares of the incoming load (preserved through the
    # triangle so the deposited material's composition matches what the
    # wind brought in).
    frac_share = Vector{Float64}(undef, nf)
    @inbounds for k in 1:nf
        frac_share[k] = q_k[k] / Q_in_total
    end
    mass_per_m = sporo1 * rho_s

    L_dep = cfg.dune_decay_length
    n_dep = max(1, ceil(Int, L_dep / dx))
    # :triangular places peak at L/3 from seaward edge. :right_triangle
    # places peak at the seaward edge (xpk = 0). :gaussian uses the
    # triangular shape then smooths.
    xpk = cfg.deposition_shape == :right_triangle ? 0.0 : (L_dep / 3.0)

    # Total deposited bed-elevation × cross-shore-length (m²) per unit alongshore
    # = (Σ q_k) · dt / sporo1 (sand vol with pores → bed elevation × length).
    # Triangle area: ∫ h(x) dx = h_peak · L / 2  (peak at xpk; same total area
    # regardless of where the peak sits along [0, L]).
    bed_vol_per_m = Q_in_total * dt / sporo1
    h_peak = 2.0 * bed_vol_per_m / L_dep

    # Build the triangle into a scratch array first; this lets us apply
    # Gaussian smoothing before depositing.
    h_profile = zeros(Float64, n_dep)
    for kshift in 0:(n_dep - 1)
        x_local = (kshift + 0.5) * dx
        if x_local >= L_dep
            h_profile[kshift + 1] = 0.0
        elseif x_local <= xpk
            h_profile[kshift + 1] = xpk > 0 ? (h_peak / xpk) * x_local : h_peak
        else
            h_profile[kshift + 1] = h_peak * (L_dep - x_local) / (L_dep - xpk)
        end
    end

    # Optional Gaussian smoothing.
    # Convolve with a discrete Gaussian truncated at ±3σ, then RESCALE
    # the result so the total area exactly matches the pre-smoothing
    # area. This guarantees mass conservation regardless of how the
    # convolution treats the boundaries.
    if cfg.deposition_shape == :gaussian
        h_profile = _gaussian_smooth_mass_conserving(h_profile, cfg.gaussian_smooth_sigma, dx)
    end

    # Apply to the bed. Anything overflowing past the landward boundary
    # is recorded as a domain loss for mass-balance bookkeeping.
    residual_overflow = 0.0
    for kshift in 0:(n_dep - 1)
        h = h_profile[kshift + 1]
        h <= 0 && continue
        j_dep = j_seaward + kshift

        if j_dep > jmax
            residual_overflow += h * dx * sporo1
            continue
        end

        @inbounds state.zb[j_dep, l] += h
        @inbounds state.ae_dz[j_dep] += h / dt
        # Per-fraction bed_mass update — split the local deposition `h`
        # between fractions in the same proportion they arrive in q_k.
        @inbounds for k in 1:nf
            _update_bed_mass!(state, j_dep, l, k,
                              h * frac_share[k], mass_per_m)
        end
    end

    # Drain the per-fraction load.
    @inbounds for k in 1:nf
        q_k[k] = 0.0
    end

    # Mass conservation: anything that overflowed past the landward edge
    # is recorded as loss.
    if residual_overflow > 0
        state.ae_loss_landward += residual_overflow
    end
    return nothing
end

"""
    _gaussian_smooth_mass_conserving(h, sigma, dx) -> Vector

Convolve `h` with a Gaussian of standard deviation `sigma` (in metres),
discretized at grid spacing `dx`, truncated at ±3σ. Boundary mode:
zero-pad. After the convolution, rescale the smoothed profile so its
total area `sum(h_smooth) · dx` exactly matches the pre-smoothing area
`sum(h_in) · dx`. The rescale is necessary because zero-padding leaks
mass at the boundaries; the rescale recovers it perfectly without
changing the smoothed *shape*.
"""
function _gaussian_smooth_mass_conserving(h::Vector{Float64},
                                          sigma::Float64, dx::Float64)
    n = length(h)
    sigma > 0 || return h
    # Build truncated Gaussian kernel (discrete, normalized to sum=1)
    half_w = max(1, ceil(Int, 3 * sigma / dx))
    inv2sig2 = 1.0 / (2 * sigma * sigma)
    kernel = Float64[exp(-((i * dx)^2) * inv2sig2) for i in -half_w:half_w]
    ksum = sum(kernel)
    kernel ./= ksum

    # Pre-smoothing area
    area_in = sum(h)

    # Zero-padded convolution
    out = zeros(Float64, n)
    @inbounds for i in 1:n
        s = 0.0
        for k in -half_w:half_w
            ji = i + k
            if 1 <= ji <= n
                s += kernel[k + half_w + 1] * h[ji]
            end
        end
        out[i] = s
    end

    # Mass-conserving rescale
    area_out = sum(out)
    if area_out > 1e-15
        out .*= area_in / area_out
    end
    return out
end

"""
    _aeolian_runup_node(state, config, l, swl) -> Int

Return the seaward index of the truly-dry beach (the first cell at or
landward of CSHORE's runup line, plus an optional `runup_buffer`
elevation offset). Returns 0 when no dry cells exist.

Source priority (CSHORE's own wet/dry tracker, taken as fact):

  • If `IOVER=1`: `state.jdry` — set by the wet/dry / overtopping march
    (wave-by-wave swash), gives the seaward edge of the truly dry zone
    after accounting for swash uprush.
  • If `IOVER=0`: `state.jswl[l]` — the still-water-level crossing node.
    The wet/dry march never runs in this mode, so SWL is the best
    available runup proxy.

`runup_buffer` (default 0) is an extra elevation buffer above the cell
returned by the source above.
"""
function _aeolian_runup_node(state::CshoreState, config::CshoreConfig,
                              l::Int, swl::Float64)
    cfg = config.aeolian
    cfg === nothing && return 0
    jmax = state.jmax[l]

    # Pick the CSHORE-supplied node: jdry when the swash march ran, jswl
    # otherwise.
    j_src = if config.options.iover == 1 && state.jdry > 0 && state.jdry <= jmax
        state.jdry
    elseif state.jswl[l] > 0 && state.jswl[l] <= jmax
        state.jswl[l]
    else
        return 0
    end

    # Optional elevation buffer above the source cell.
    if cfg.runup_buffer <= 0
        return j_src
    end
    z_target = state.zb[j_src, l] + cfg.runup_buffer
    @inbounds for j in j_src:jmax
        if state.zb[j, l] >= z_target
            return j
        end
    end
    return jmax
end

"""
    _update_bed_mass!(state, j, l, k, dz_k, mass_per_m)

Apply a per-fraction bed-elevation change `dz_k` (m, signed:
positive = deposition, negative = erosion) to the surface layer of
cell `j` for grain-size fraction `k`. Multiplies by `mass_per_m =
sporo1 · ρ_s` (kg/m³) to convert metres-of-bed to kg/m² and updates
`state.bed_mass[j, 1, k]` in place. Floors mass at 0 to avoid
unphysical negative composition (any drawdown beyond the surface
layer's content for that fraction is silently capped — the
active-layer recompute at the next BC window absorbs the bookkeeping).

This is the link between the aeolian zb update and the multifraction
composition machinery: without it, the surface fractions stay at their
initial values and armouring never appears.
"""
@inline function _update_bed_mass!(state::CshoreState, j::Int, l::Int, k::Int,
                                   dz_k::Float64, mass_per_m::Float64)
    dm = dz_k * mass_per_m
    new_mass = state.bed_mass[j, 1, k] + dm
    state.bed_mass[j, 1, k] = max(new_mass, 0.0)
    return nothing
end

"""
Linear interpolation helper local to this module so we don't depend on
`driver.jl` being loaded first.
"""
@inline function _aeolian_interp(values::Vector{Float64},
                                 times::Vector{Float64}, t::Float64)
    n = length(times)
    n == 0 && return 0.0
    t <= times[1]   && return values[1]
    t >= times[n]   && return values[n]
    i = searchsortedlast(times, t)
    i = clamp(i, 1, n - 1)
    w = (t - times[i]) / (times[i + 1] - times[i])
    return (1 - w) * values[i] + w * values[i + 1]
end
