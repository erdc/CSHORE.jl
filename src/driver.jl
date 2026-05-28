# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
driver.jl — Top-level time-stepping orchestrator.

The nested loop structure is:

    for ITIME in 1:NTIME          # boundary-condition time steps
      for L in 1:ILINE            # cross-shore lines
        for ITEQO in 1:ITEMAX     # overtopping iteration (IOVER=1 only)
          # Landward march:
          for J in 1:JMAX-1
            for ITE in 1:MAXITE   # wave-solver convergence
              call LWAVE, DBREAK, DVEG, POFLOW, GBXAGF, VSTGBY, ...
            end
          end
          if IOVER==1: CALL QORATE
        end
        if IPROFL==1:
          CALL SEDTRA; CALL CHANGE; CALL WETDRY
        end
        CALL OUTPUT
      end
    end
==============================================================================#

"""
    run_simulation!(config::CshoreConfig; outdir=".", outfile=nothing,
                    ascii_outfile=nothing, ...) -> CshoreState

Entry point. Allocates state, applies initial bathymetry, and runs the full
CSHORE time loop.

## Output writers — independently opt-in

The two output writers are entirely independent. Either, both, or neither
can be enabled per run:

| Combination | NetCDF | ASCII |
|---|---|---|
| `run_simulation!(cfg)` | off | off |
| `run_simulation!(cfg; outfile="run.nc")` | on | off |
| `run_simulation!(cfg; ascii_outfile="run")` | off | on |
| `run_simulation!(cfg; outfile="run.nc", ascii_outfile="run")` | on | on |

- `outfile::Union{Nothing,String}` — when set, opens a CF-1.10 NetCDF file
  at `joinpath(outdir, outfile)` (writer in `io/netcdf.jl`).
- `ascii_outfile::Union{Nothing,String}` — when set, opens the standard
  CSHORE_USACE ASCII files (`OBPROF.OUT`, `OENERG.OUT`, …) in `outdir`
  using the given string as a filename prefix (`""` means no prefix). Set
  to `"run"` to get `runOBPROF.OUT` etc. (writer in `io/output.jl`).

Other kwargs:
- `output_interval_s` — NetCDF write cadence (0 = every step).
- `write_composition`, `write_T_profile`, `write_transport` — NetCDF variable groups.
- `provenance` — optional provenance tracking, both writers honor it.
- `progress_callback` — GUI hook called after each BC window.
"""
function run_simulation!(config::CshoreConfig;
                         outdir::AbstractString=".",
                         outfile::Union{Nothing,AbstractString}=nothing,
                         output_interval_s::Real=0.0,
                         write_composition::Bool=true,
                         write_T_profile::Bool=true,
                         write_transport::Bool=true,
                         ascii_outfile::Union{Nothing,AbstractString}=nothing,
                         provenance::Union{ProvenanceConfig,Nothing}=nothing,
                         progress_callback::Union{Function,Nothing}=nothing)
    state = initialize_state(config)
    apply_initial_bathymetry!(state, config)
    compute_derived_constants!(state, config)

    # Bed slope per line (needed by DBREAK's steep-slope factor)
    for l in 1:config.options.iline
        compute_bed_slope!(state, config, l)
    end

    # Vegetation uprooting: initialize fixed elevation bounds (vegzd, vegzr)
    # and the uproot flag. `state.vegh` starts equal to the config canopy
    # height and then evolves with the bed via `check_vegetation_uprooting!`.
    initialize_vegetation_bounds!(state, config)

    # Optional thermal / permafrost submodel. When configured, we allocate
    # the thermal state, seed the initial ALT (zero — winter start), and
    # prime `state.zb_hard` from `zb - ALT` so the very first wave sub-step
    # already sees a consistent hardbottom. The per-BC thermal advance
    # happens inside `step_bc_window!` below.
    if config.thermal !== nothing
        state.thermal = initialize_thermal_state(config.thermal, state)
        for l in 1:config.options.iline
            jmax = state.jmax[l]
            @inbounds for j in 1:jmax
                state.zb_hard[j, l] = state.zb[j, l] - state.thermal.ALT[j]
            end
        end
    end

    # Initialize provenance state if requested
    prov_state = provenance === nothing ? nothing :
                 init_provenance(provenance, state, config)

    # Optional NetCDF recorder (opt-in via `outfile` kwarg)
    writer = outfile === nothing ? nothing :
             open_netcdf(joinpath(outdir, outfile), config, state;
                         output_interval_s=output_interval_s,
                         write_composition=write_composition,
                         write_T_profile=write_T_profile,
                         write_transport=write_transport,
                         provenance=provenance,
                         prov_state=prov_state)

    # Optional FORTRAN-compatible ASCII recorder (opt-in via `ascii_outfile`).
    # The string is used as a filename prefix in `outdir`; pass "" for no prefix.
    ascii_writer = ascii_outfile === nothing ? nothing :
        open_ascii_outputs(config; outdir=outdir, prefix=String(ascii_outfile))

    try
        # Record the initial state at t=0 (before any BC window advances)
        if writer !== nothing
            write_step!(writer, state, state.time; prov_state=prov_state)
        end
        if ascii_writer !== nothing
            write_ascii_step!(ascii_writer, state, config, state.time)
        end

        ntimes = ntime(config.boundary)
        # ILAB=0 (field mode): timebc has NWAVE+1 entries; NWAVE windows between them.
        # ILAB=1 (lab mode):   timebc has NWAVE entries;   each entry is its own run.
        n_steps = config.options.ilab == 1 ? ntimes : ntimes - 1
        n_steps_total = max(1, n_steps)
        for itime in 1:n_steps_total
            step_bc_window!(state, config, itime; prov_state=prov_state,
                            prov_cfg=provenance)
            if writer !== nothing
                maybe_write_step!(writer, state, state.time; prov_state=prov_state)
            end
            if ascii_writer !== nothing
                write_ascii_step!(ascii_writer, state, config, state.time)
            end
            # Optional progress hook — used by the GUIs to update a progress
            # bar between BC windows. Wrapped in a try/catch so a misbehaving
            # callback can't take down the simulation.
            if progress_callback !== nothing
                try
                    progress_callback(itime, n_steps_total, state)
                catch err
                    @warn "run_simulation!: progress_callback raised — disabling" exception=err
                    progress_callback = nothing
                end
            end
        end
    finally
        writer !== nothing && close_netcdf!(writer)
        ascii_writer !== nothing && close_ascii_outputs!(ascii_writer)
    end

    return state
end

"""
    step_bc_window!(state, config, itime) -> Nothing

Advance the simulation by one boundary-condition window — from
`config.boundary.timebc[itime]` to `config.boundary.timebc[itime+1]`. All
cross-shore lines are advanced. On return `state.time == timebc[itime+1]`
(modulo floating-point) and `state.itime == itime`.

This is the natural "one step" unit for the Basic Model Interface
(`BMI.update`). It's also the inner unit that `run_simulation!` calls
in a loop.
"""
function step_bc_window!(state::CshoreState, config::CshoreConfig, itime::Int;
                         prov_state::Union{ProvenanceState,Nothing}=nothing,
                         prov_cfg::Union{ProvenanceConfig,Nothing}=nothing)
    state.itime = itime
    state.time  = config.boundary.timebc[itime]
    wave = _current_wave_params(config, itime)
    ntimes = ntime(config.boundary)
    t_window_end = if ntimes ≥ itime + 1
        config.boundary.timebc[itime + 1]
    elseif ntimes ≥ 2
        # Last ILAB=1 step: extrapolate duration from the preceding interval.
        state.time + (config.boundary.timebc[ntimes] - config.boundary.timebc[ntimes - 1])
    else
        state.time + 3600.0
    end

    # Advance the thermal submodel BEFORE the wave / sediment sub-steps so
    # the hardbottom (state.zb_hard) reflects the latest active-layer
    # thickness for this BC window. We use the mean `(T_air, T_water)` over
    # the BC window as a simple time-average forcing.
    if config.thermal !== nothing && config.thermal_bc !== nothing
        dt_bc = t_window_end - state.time
        if dt_bc > 0
            Ta = _interp_thermal(config.thermal_bc.T_air,
                                 config.thermal_bc.time, state.time, t_window_end)
            Tw = _interp_thermal(config.thermal_bc.T_water,
                                 config.thermal_bc.time, state.time, t_window_end)
            # Snow insulation: update snow depth and R_insulation before
            # the thermal solve so the heat equation sees current snow cover.
            if config.snow !== nothing
                prescribed = if !isempty(config.thermal_bc.snow_depth)
                    _interp_thermal(config.thermal_bc.snow_depth,
                                    config.thermal_bc.time, state.time, t_window_end)
                else
                    NaN   # use degree-day model
                end
                for l in 1:config.options.iline
                    update_snow!(state.thermal, config.snow, state, config,
                                 Ta, dt_bc, l; prescribed_depth=prescribed)
                end
            end
            for l in 1:config.options.iline
                update_thermal!(state, config, config.thermal, state.thermal,
                                Ta, Tw, dt_bc, l)
            end
        end
    end

    # Wind stress: interpolate W10, angle, Cd at the current BC time and
    # convert to the dimensionless stress used in the momentum equation.
    # tau_w = rho_air * Cd * W10^2;  twxsta = tau_w * cos(angle) / (rho_w * g)
    if config.options.iwind == 1 && !isempty(config.boundary.w10)
        bc = config.boundary
        t_now = state.time
        w10   = _interp_at(bc.w10, bc.timebc, t_now)
        wangl = _interp_at(bc.wangle, bc.timebc, t_now)
        wcd   = !isempty(bc.windcd) ? _interp_at(bc.windcd, bc.timebc, t_now) : 0.0015
        rho_ratio = 1.225 / 1025.0   # rho_air / rho_water
        tau_norm  = rho_ratio * wcd * w10^2 / GRAV
        ang_rad   = deg2rad(wangl)
        state.twxsta = tau_norm * cos(ang_rad)
        state.twysta = tau_norm * sin(ang_rad)
    end

    # Tidal forcing (ITIDE=1): interpolate alongshore water-surface gradient
    # DETADY and (for ILAB=0) SWL rate of change DSWLDT at current time, then
    # build the per-node cross-shore tidal volume flux QTIDE from the
    # instantaneous tidal prism geometry.
    if config.options.itide == 1 && config.tidal !== nothing
        ti = config.tidal
        state.detady_now = _interp_at(ti.detady, ti.time, state.time)
        if config.options.ilab == 0
            state.dswldt_now = if !isempty(ti.dswldt)
                _interp_at(ti.dswldt, ti.time, state.time)
            else
                _dswldt_from_swlbc(config.boundary, state.time)
            end
            # Build QTIDE[j] = (xb[jswl] - xb[j]) * dswldt for j < jswl
            fill!(state.qtide, 0.0)
            for l_i in 1:config.options.iline
                jswl = state.jswl[l_i]
                jmax = state.jmax[l_i]
                if jswl > 0 && jswl <= jmax
                    xswl = state.xb[jswl]
                    @inbounds for j in 1:jswl
                        state.qtide[j] = (xswl - state.xb[j]) * state.dswldt_now
                    end
                end
            end
        else
            state.dswldt_now = 0.0
        end
    else
        state.detady_now = 0.0
        state.dswldt_now = 0.0
    end

    # Imposed alongshore current at offshore boundary (ICURRENT=1).
    # Interpolate the user-supplied vbc time series here; the actual
    # back-solve for DETADY is done inside transform_waves! at the most-
    # offshore valid cell (after wave kinematics are known).
    if config.options.icurrent == 1 && config.current !== nothing
        ci = config.current
        state.vbc_now = _interp_at(ci.vbc, ci.time, state.time)
        # Reset per-line latch so the inverse solve runs once per BC window
        # at the most-offshore valid cell of each cross-shore line.
        fill!(state.icurrent_solved, false)
    else
        state.vbc_now = 0.0
    end

    for l in 1:config.options.iline
        _run_line_step!(state, config, wave, itime, l, t_window_end;
                        prov_state=prov_state, prov_cfg=prov_cfg)
    end

    # Aeolian (wind-driven) sediment transport — IAEOLIAN=1.
    # Runs once per BC window AFTER the wave-driven hydrodynamics +
    # morphodynamics so that the active layer (which sets the surface
    # composition for per-fraction Kawamura thresholds) reflects the
    # post-wave bed. Uses the same BC cadence as waves.
    if config.options.iaeolian == 1 && config.aeolian !== nothing
        dt_ae = t_window_end - config.boundary.timebc[itime]
        if dt_ae > 0
            # Refresh hp = max(0, zb - zb_hard) per line so the aeolian
            # kernel can see how much sediment is available above the
            # hardbottom (revetment, bedrock, permafrost active-layer
            # base, etc.) before entraining.
            for l in 1:config.options.iline
                _update_hp_from_hardbottom!(state, l)
            end
            step_aeolian!(state, config, itime; dt = dt_ae)
        end
    end

    # Hillslope diffusion: applied AFTER wave/morphodynamic sub-steps so that
    # wave-oversteepened slopes get relaxed before the next BC window. Runs
    # once per BC window (same cadence as the thermal update).
    if config.diffusion !== nothing
        dt_bc = t_window_end - config.boundary.timebc[itime]
        if dt_bc > 0
            for l in 1:config.options.iline
                apply_hillslope_diffusion!(state, config, l, dt_bc)
            end
        end
    end

    return nothing
end

"""
    _interp_at(values, times, t) -> Float64

Linear interpolation of a piecewise-linear time series at a single time.
Clamps to endpoints for out-of-range t.
"""
function _interp_at(values::Vector{Float64}, times::Vector{Float64}, t::Float64)
    n = length(times)
    n == 0 && return 0.0
    t ≤ times[1] && return values[1]
    t ≥ times[n] && return values[n]
    i = searchsortedlast(times, t)
    i = clamp(i, 1, n - 1)
    w = (t - times[i]) / (times[i + 1] - times[i])
    return (1 - w) * values[i] + w * values[i + 1]
end

"""
    _dswldt_from_swlbc(boundary, t) -> Float64

Finite-difference rate of change of SWL at time `t`, computed from the
piecewise-linear `boundary.swlbc` time series. Returns 0 at the endpoints
or when only one time point exists. Used for ITIDE=1 + ILAB=0 to derive
DSWLDT when it isn't explicitly supplied.
"""
function _dswldt_from_swlbc(boundary, t::Float64)
    tb = boundary.timebc
    sb = boundary.swlbc
    n = length(tb)
    (n < 2) && return 0.0
    t ≤ tb[1] && return (sb[2] - sb[1]) / (tb[2] - tb[1])
    t ≥ tb[n] && return (sb[n] - sb[n - 1]) / (tb[n] - tb[n - 1])
    i = searchsortedlast(tb, t)
    i = clamp(i, 1, n - 1)
    return (sb[i + 1] - sb[i]) / (tb[i + 1] - tb[i])
end

"""
    _interp_thermal(values, times, t0, t1) -> Float64

Trapezoid-averaged value of a piecewise-linear time series `(times, values)`
over the interval `[t0, t1]`. Used to get a mean air or water temperature
for a BC window. Falls back to the nearest endpoint for out-of-range t0/t1.
"""
function _interp_thermal(values::Vector{Float64}, times::Vector{Float64},
                         t0::Float64, t1::Float64)
    n = length(times)
    n == 0 && return 0.0
    if t1 ≤ times[1];  return values[1];   end
    if t0 ≥ times[n];  return values[n];   end
    # Linear interpolation at the two endpoints, then average.
    @inline function vat(t)
        if t ≤ times[1]; return values[1]; end
        if t ≥ times[n]; return values[n]; end
        i = searchsortedlast(times, t)
        i = clamp(i, 1, n - 1)
        w = (t - times[i]) / (times[i + 1] - times[i])
        return (1 - w) * values[i] + w * values[i + 1]
    end
    return 0.5 * (vat(t0) + vat(t1))
end

"""
    _current_wave_params(config, itime) -> WaveParams

Read the offshore boundary conditions at time index `itime` and pack into a
`WaveParams` struct. Handles the fact that `wangbc` is stored in degrees in
the input file but used internally as the raw angle value (LWAVE converts
to radians via its own degree-based Snell math).
"""
function _current_wave_params(config::CshoreConfig, itime::Int)
    bc = config.boundary
    tp    = bc.tpbc[itime]
    hrms0 = bc.hrmsbc[itime]
    angle = bc.wangbc[itime]
    swl   = bc.swlbc[itime]
    wkpo  = (PI2 / tp)^2 / GRAV
    return WaveParams(; tp=tp, hrms0=hrms0, angle=angle, swl=swl, wkpo=wkpo)
end

"""
    _run_line_step!(state, config, wave, itime, l, t_window_end)

One (itime, l) cell of the nested time loop:

1. Initialize SWL depth and JSWL shoreline node.
2. Sub-step loop until `state.time` reaches `t_window_end`:
   a. Solve the wave field via `transform_waves!`.
   b. If `iprofl == 0` (fixed bed), set `delt = t_window_end - time` and
      break (one wave solve covers the whole BC window).
   c. If `iprofl == 1`, run `sedtra!` + `exner_step!` + `apply_bed_update!`,
      advance `state.time += delt`, and repeat unless `iend == true`.
3. Recompute bed slope after bed updates so the next wave solve sees the
   current bathymetry.
"""
function _locate_crest!(state::CshoreState, config::CshoreConfig, l::Int)
    jmax_l = state.jmax[l]
    rcrest = config.swash.rcrest
    if !isempty(rcrest) && rcrest[min(l, length(rcrest))] > -1e10
        # User-supplied crest elevation: first node at or above it.
        rc = rcrest[min(l, length(rcrest))]
        state.jcrest[l] = 1
        @inbounds for j in 1:jmax_l
            if state.zb[j, l] ≥ rc
                state.jcrest[l] = j
                break
            end
        end
    else
        # Auto-detect: highest point on the profile.
        zmax = -Inf
        state.jcrest[l] = jmax_l
        @inbounds for j in 1:jmax_l
            if state.zb[j, l] > zmax
                zmax = state.zb[j, l]
                state.jcrest[l] = j
            end
        end
    end
    return nothing
end

"""
    _morphodynamic_substep!(state, config, l, t_window_end;
                            prov_state=nothing, prov_cfg=nothing) -> iend::Bool

One morphodynamic sub-step on line `l`: sediment transport → adaptive Δt →
per-fraction mass balance + Hirano rebalance → hardbottom clamp →
underwater avalanche → optional ICLAY/IVEG hooks. Returns `iend=true` when
`compute_timestep!` signals the BC window has been reached.

Assumes the hydrodynamic sub-step (wave transform + swash + IG + groundwater)
has already run for the same `(state, l)`.
"""
function _morphodynamic_substep!(state::CshoreState, config::CshoreConfig,
                                 l::Int, t_window_end::Float64;
                                 prov_state::Union{ProvenanceState,Nothing}=nothing,
                                 prov_cfg::Union{ProvenanceConfig,Nothing}=nothing)
    # Aeolis-style virtual adaptive active layer + smoothing (replaces the
    # earlier mix_top_layer! / per-fraction flux smoothing in sedtra!).
    active_layer_fractions!(state, config, l)
    smooth_active_composition!(state, config, l)

    # Refresh per-fraction BRF supply tracker (sand-thickness-above-hardbottom).
    _update_hp_from_hardbottom!(state, l)

    sedtra!(state, config, l)

    # Supply-factor (hardbottom / permafrost availability) applied to fluxes
    # before the bed update; sets the flag asserted by compute_timestep!.
    _apply_supply_factor!(state, config, l)
    state.supply_factor_applied = true

    # 1. Adaptive Δt + per-fraction pickup divergence
    result = compute_timestep!(state, config, l, t_window_end)

    # Provenance: snapshot surface-layer mass BEFORE the bed update.
    jmax_l = state.jmax[l]
    bed_mass_before = if prov_state !== nothing
        copy(state.bed_mass[1:jmax_l, 1, :])
    else
        nothing
    end

    # 2. Per-fraction mass balance + Hirano active-layer rebalance + zb recompute.
    # This is the ONLY function that mutates bed_mass / zb / delzb on a morpho step.
    update_bed_composition!(state, config, l, result.delt)

    if prov_state !== nothing && prov_cfg !== nothing && bed_mass_before !== nothing
        step_provenance!(prov_state, prov_cfg, bed_mass_before, state, config, l)
    end

    # 3. Hardbottom clamp + underwater angle-of-repose avalanche.
    apply_hardbottom_clamp!(state, config, l)
    apply_underwater_avalanche!(state, config, l)

    # 4. ICLAY=1: additive clay erosion below the sand layer when HP<D50.
    if config.options.iclay == 1
        eroson!(state, config, l; t_window_end=t_window_end)
    end

    # 5. Vegetation uprooting check (IVEG=1 or 3, IPROFL=1).
    if config.options.iveg in (1, 3) && config.options.iprofl == 1
        check_vegetation_uprooting!(state, config, l)
    end

    # 6. Cohesive (mud) sediment Partheniades-Krone step. Operates on a
    # separate bed-mass / concentration pool from the sand multifraction
    # stack and does not feed back to zb in v1. No-op when
    # config.cohesive === nothing.
    cohesive_step!(state, config, l, result.delt)

    state.time += result.delt
    return result.iend
end

function _run_line_step!(state::CshoreState, config::CshoreConfig,
                         wave::WaveParams, itime::Int, l::Int,
                         t_window_end::Float64;
                         prov_state::Union{ProvenanceState,Nothing}=nothing,
                         prov_cfg::Union{ProvenanceConfig,Nothing}=nothing)
    _reject_out_of_scope_flags(config)

    # Seaward boundary conditions (set once per BC window).
    # Note: wsetup[1] is initialized to 0.0 here but will be recalculated from
    # the radiation stress at the seaward boundary during each substep's transform_waves!
    state.wsetup[1] = 0.0
    state.hrms[1]   = wave.hrms0
    state.qo[l]     = 0.0

    iover = config.options.iover

    # Locate crest node on the current bathymetry. Always needed by wetdry!.
    _locate_crest!(state, config, l)

    max_substeps = 100_000
    for _substep in 1:max_substeps
        # Recompute swldep / jswl against current zb
        swl = wave.swl
        jmax_l = state.jmax[l]
        ichk = false
        state.jswl[l] = jmax_l
        @inbounds for j in 1:jmax_l
            d = swl - state.zb[j, l]
            state.swldep[j, l] = d
            if !ichk && d < 0.0
                state.jswl[l] = j
                ichk = true
            end
        end
        # Bed slope (recomputed each sub-step since zb evolves)
        compute_bed_slope!(state, config, l)

        # Reset wsetup at boundary each sub-step
        state.wsetup[1] = 0.0

        # ---- Outer iterative WCI loop (Picard) when iwcint_along=1 ----
        # Each outer pass uses the previous-pass vmean profile to compute the
        # alongshore-current Doppler contribution to qdisp. With iwcint_along=0
        # the loop runs exactly once and degenerates to the FORTRAN-parity path.
        # First pass uses vmean_prev = current vmean (zeros on cold start, or
        # last-substep values on warm restart).
        wci_along_active = config.options.iwcint_along == 1 &&
                           config.options.iwcint == 1 &&
                           config.options.iangle == 1
        wci_along_maxite = wci_along_active ? config.iwcint_along_maxite : 1
        wci_along_tol    = config.iwcint_along_tol

        # Seed vmean_prev with the current vmean (so first pass = baseline).
        if wci_along_active
            jmax_seed = state.jmax[l]
            @inbounds for j in 1:jmax_seed
                state.vmean_prev[j] = state.vmean[j]
            end
        end

        wci_residual = NaN
        wci_iters = 0
        for iwci_pass in 1:wci_along_maxite
            wci_iters = iwci_pass

            # ---- Wave transform + swash hydrodynamics ----
            # iover=1: iterate transform_waves! + wetdry! until qo converges.
            # iover=0: single-pass wave transform, then wetdry! once with qo=0.
            # In both cases tranwd! blends wave and swash fields when a swash
            # zone exists (jdry > jr), so that sedtra! always sees consistent
            # h, ustd, umean across the wet→swash transition.
            if iover == 1
                itemax = config.options.iprofl ≥ 1 ? 4 : 20
                for iteqo in 1:itemax
                    transform_waves!(state, config, wave, l)
                    converged = overtopping_rate!(state, config, itime, l, iteqo)
                    if converged; break; end
                end
            else
                transform_waves!(state, config, wave, l)
                # Single-pass swash: qo[l] stays 0, no iteration needed.
                overtopping_rate!(state, config, itime, l, 1)
            end

            # Blend wave field (1:jr) with swash field (jwd:jdry).
            # Runs whenever a swash zone exists, regardless of iover.
            jr   = state.jr
            jwd  = state.jwd
            jdry = state.jdry
            if jdry > jr
                tranwd!(state.h,    jr, state.hwd,    jwd, jdry)
                tranwd!(state.sigma, jr, state.sigwd,  jwd, jdry)
                tranwd!(state.umean, jr, state.umeawd, jwd, jdry)
                tranwd!(state.ustd,  jr, state.ustdwd, jwd, jdry)
            end
            # Mark the wet zone as fully wet (pwet = 1 seaward of jwd)
            @inbounds for j in 1:min(jwd, jmax_l)
                state.pwet[j] = 1.0
            end

            # ---- Landward wave transmission (IWTRAN=1) -----------------------
            # Drives the back-side wave field for low-crested breakwaters and
            # dune-overtopping cases. No-op when iwtran=0.
            if config.options.iwtran == 1
                transmission!(state, config, l)
            end

            # ---- Infragravity wave energy (IgConfig) -------------------------
            # Runs after the short-wave field is converged so that hrms, h,
            # dbsta, and dfsta are all up to date. Populates state.hrms_ig and
            # state.ustd_ig; is a no-op when config.ig === nothing.
            compute_ig_field!(state, config, wave, l, itime)

            # ---- Beach groundwater + surface moisture (GroundwaterConfig) ----
            # Runs after wave+swash so that wsetup, jswl, pwet, hwd, and
            # wpm_derived are all current.  Advances the Boussinesq water-table
            # field over the BC window and updates state.theta (Van Genuchten).
            # No-op when config.groundwater === nothing.
            if config.groundwater !== nothing
                step_groundwater!(state, config, wave, l,
                                  t_window_end - config.boundary.timebc[itime])
            end

            # Convergence check on max change in vmean across the line
            if !wci_along_active
                break
            end
            jmax_cur = state.jmax[l]
            local res = 0.0
            @inbounds for j in 1:jmax_cur
                d = abs(state.vmean[j] - state.vmean_prev[j])
                d > res && (res = d)
            end
            wci_residual = res
            if res < wci_along_tol
                break
            end
            # Update vmean_prev for next outer pass
            @inbounds for j in 1:jmax_cur
                state.vmean_prev[j] = state.vmean[j]
            end
        end

        if config.options.iprofl == 0
            # Fixed bed: no sub-stepping, one wave solve covers the whole BC window.
            state.delt = t_window_end - state.time
            state.time = t_window_end
            return state
        end

        # IPROFL=2: grassed-dike erosion entirely replaces the sediment
        # transport + Exner chain. EROSON computes its own DELT from the
        # DELEM/DMAX stability criterion and updates ZB directly from the
        # piecewise-grass erosion integral.
        if config.options.iprofl == 2
            # EROSON needs ZB0 (initial profile); capture on first call.
            if !state.eroson_initialized
                @inbounds for j in 1:state.jmax[l]
                    state.zb0[j, l] = state.zb[j, l]
                end
            end
            iend, dt_used = eroson!(state, config, l; t_window_end=t_window_end)
            state.time += dt_used
            if iend || state.time >= t_window_end
                return state
            end
            continue
        end

        # --- Morphodynamic sub-step (Aeolis-style per-fraction update) ---
        if _morphodynamic_substep!(state, config, l, t_window_end;
                                   prov_state=prov_state, prov_cfg=prov_cfg)
            return state
        end
    end
    error("_run_line_step!: exceeded $max_substeps sub-steps in BC window " *
          "($(t_window_end - state.time) s remaining) — check CFL / adaptive dt")
end

"""
    _apply_supply_factor!(state, config, l)

Multiply the per-fraction transport fluxes `state.qbx`, `state.qsx` by a
node-dependent supply factor `fac_supply(j) ∈ [0, 1]` that represents
hardbottom or permafrost sand-availability limitation. This must be done
exactly once per sub-step; `state.supply_factor_applied` is asserted by
`exner_step!`.

Per-fraction BRF supply limitation is applied inside `_transport_kernel`
when `abs(isedav) ≥ 1`, using the `state.hp[j, l]` tracker that
`_update_hp_from_hardbottom!` refreshed. This function is the enforcement
point for the supply_factor_applied flag, keeping any additional global
limiters (e.g. permafrost thaw rate) localized here.
"""
function _apply_supply_factor!(state::CshoreState, config::CshoreConfig, l::Int)
    return state
end

"""
    _update_hp_from_hardbottom!(state, l)

Refresh `state.hp[j, l] = max(0, state.zb[j, l] - state.zb_hard[j, l])`
for every active node on line `l`. Cells with `zb_hard = -Inf` (the
sentinel for "no hardbottom") get `hp = +Inf`, which makes the BRF
formula in `_transport_kernel` collapse to 1.0 (no reduction).

This must run BEFORE `sedtra!` so the transport kernel sees the current
sand-thickness-over-hardbottom when computing the per-fraction BRF.
"""
function _update_hp_from_hardbottom!(state::CshoreState, l::Int)
    jmax_l = state.jmax[l]
    @inbounds for j in 1:jmax_l
        zh = state.zb_hard[j, l]
        if zh == -Inf
            state.hp[j, l] = Inf
        else
            d = state.zb[j, l] - zh
            state.hp[j, l] = d > 0 ? d : 0.0
        end
    end
    return state
end

"""
    _reject_out_of_scope_flags(config)

Guard against unsupported flag combinations. Errors with a clear message
pointing at the unsupported subroutine.
"""
function _reject_out_of_scope_flags(config::CshoreConfig)
    opt = config.options
    # IOVER=0,1 now supported (wave overtopping via WETDRY + QORATE)
    # IPERM=0,1 now supported (porous layer flow via POFLOW)
    # IVEG=1,2,3 now supported (vegetation dissipation + friction enhancement)
    # IWCINT=0,1 now supported (wave-current interaction in dispersion)
    # IROLL=0,1 now supported (wave roller energy equation)
    # IWTRAN=0,1 now supported (landward wave transmission via transmission!)
    return nothing
end

"""
    _transform_waves_seaward_bc!(state, config, wave, l) -> Bool

Set the j=1 seaward boundary state for the wave-energy-flux march. Returns
`true` if the march should proceed, `false` if the seaward node is dry
(in which case `state.jr/xr/zr` are pinned to j=1 and the caller should
return early without marching).

Side effects on `state` at j=1: hrms (floored to 1e-5), sigma, h, sigsta,
qwx, qwy, sxxsta, efsta, sxysta, gbx, gf, tbxsta, dfsta, rq, wsetup, plus
the lwave / dbreak / porous / veg dissipation kernels.
"""
function _transform_waves_seaward_bc!(state::CshoreState, config::CshoreConfig,
                                      wave::WaveParams, l::Int)
    eps1   = config.eps1
    iangle = config.options.iangle
    iveg   = config.options.iveg
    iroll  = config.options.iroll
    iwcint = config.options.iwcint

    # Wave-energy-flux march termination requires a non-zero boundary energy
    # to step at all (the predictor aborts when efsta[1] ≤ 0). When the user
    # prescribes Hrms=0 (e.g. fluvial cases forced only by ICURRENT/ITIDE),
    # floor Hrms[1] to keep the march alive. The resulting sigt stays well
    # below SIGT_NOWAVE so the longshore-momentum solver falls back to the
    # wave-bypass quadratic-friction branch and transport stays current-only.
    if state.hrms[1] < 1e-5
        state.hrms[1] = 1e-5
    end
    state.sigma[1]  = state.hrms[1] / SQR8
    state.h[1]      = state.wsetup[1] + state.swldep[1, l]
    if state.h[1] ≤ eps1
        state.jr = 1
        state.xr = state.xb[1]
        state.zr = state.zb[1, l]
        return false
    end
    state.sigsta[1] = state.sigma[1] / state.h[1]
    state.qwx = 0.0; state.qwy = 0.0

    # Dynamic Manning's n: refresh fb2[1,l] from state.h[1] before it is read.
    if config.options.ifriction_spatial == 2
        _update_fb2_manning!(state, config, 1, l)
    end

    lwave!(state, config, 1, l, state.h[1], wave; hrms_j=state.hrms[1])

    sigma2 = state.sigma[1]^2
    state.sxxsta[1] = sigma2 * state.fsx
    state.efsta[1]  = sigma2 * state.fe
    if iangle == 1
        state.sxysta[1] = sigma2 * state.fsy
    end

    if state.fb2[1, l] > 0.0
        usigt = -state.sigsta[1] * GRAV * state.h[1] / state.cp[1]^2
        if iangle == 1
            usigt *= state.ctheta[1]
        end
        dum = state.sigsta[1] * state.cp[1]   # = sigt
        gbx, gf = friction_coefficients(state.ctheta[1], usigt, state.stheta[1], 0.0, iangle)
        state.gbx[1] = gbx
        state.gf[1]  = gf
        state.tbxsta[1] = state.fb2[1, l] * gbx * dum^2 / GRAV
        state.dfsta[1]  = state.fb2[1, l] * gf  * dum^3 / GRAV
        if iveg in (1, 2)
            apply_veg_friction!(state, config, 1, l)
        end
    else
        state.tbxsta[1] = 0.0
        state.dfsta[1]  = 0.0
    end

    dbreak!(state, config, 1, l, state.hrms[1], state.h[1], state.wt[1])
    if config.options.iperm == 1
        porous_flow!(state, config, 1, l, state.wn[1] * state.sigma[1], 0.0)
    end
    if iveg == 3
        veg_dissipation!(state, config, 1, l)
    end

    # Roller BC: zero at seaward boundary
    if iroll == 1
        state.rq[1] = 0.0
        state.sxxsta[1] += state.rx[1] * state.rq[1]   # = 0, for consistency
    end

    # Seaward setup from radiation-stress momentum balance:
    # η = Sxx / (ρ g h). Replaces the default zero BC so the boundary is
    # self-consistent with the nearshore momentum balance landward.
    h_1 = state.h[1]
    state.wsetup[1] = h_1 > eps1 ? state.sxxsta[1] / (GRAV * h_1) : 0.0

    # Wave-current interaction: cross-shore water flux at j=1, QDISP = UMEAN·H.
    if iwcint == 1
        if state.h[1] > config.min_depth_wcint
            usigt_1 = state.fb2[1, l] > 0.0 ?
                -state.sigsta[1] * GRAV * state.h[1] / state.cp[1]^2 : 0.0
            if iangle == 1; usigt_1 *= state.ctheta[1]; end
            sigt_1 = state.sigsta[1] * state.cp[1]
            state.qwx = usigt_1 * sigt_1 * state.h[1]
        else
            state.qwx = 0.0
        end
        if config.options.itide == 1 && config.options.ilab == 0
            state.qwx += state.qtide[1]
        end
    end
    return true
end

"""
    transform_waves!(state, config, wave, l)

Landward-marching wave solver. Solves, node by node from J=1 to J=JMAX-1,
the coupled:
- wave-action (energy-flux) equation for σ²
- cross-shore momentum equation for setup η̄

via a **predictor-corrector** iteration (outer loop is the J march; inner
loop is the MAXITE convergence sweep that uses trapezoidal averaging).

Boundary condition: at J=1, `HRMS = hrms0`, `WSETUP = 0`.
Termination: either J reaches JMAX-1 or the water depth drops below EPS1.
"""
function transform_waves!(state::CshoreState, config::CshoreConfig,
                          wave::WaveParams, l::Int)
    eps1   = config.eps1
    eps2   = config.eps2
    maxite = config.maxite
    dx     = config.grid.dx
    dxd2   = 0.5 * dx
    iangle = config.options.iangle
    iveg   = config.options.iveg
    iroll  = config.options.iroll
    iwcint = config.options.iwcint
    iwcint_along = config.options.iwcint_along

    # Helper: compose the effective qdisp at node `j`. With IWCINT_ALONG=1,
    # adds the alongshore-current Doppler contribution v·sin(θ)·h, where
    # `vmean_prev` carries the alongshore current from the previous outer
    # iteration (zeros on the first pass, matching FORTRAN parity).
    @inline qdisp_at(j::Int) =
        if iwcint == 0
            0.0
        elseif iwcint_along == 1 && iangle == 1
            state.qwx + state.vmean_prev[j] * state.stheta[j] * state.h[j]
        else
            state.qwx
        end

    # J=1 boundary: set hrms/sigma/h/sigsta, lwave, friction, dbreak, porous,
    # veg, roller, setup, qwx. Returns false if the seaward node is dry, in
    # which case the entire landward march is skipped.
    _transform_waves_seaward_bc!(state, config, wave, l) || return state

    # ---- Landward march J = 1, 2, ... -----------------------------------
    jdum = state.jmax[l]
    jp1 = 1
    march_aborted = false
    for j in 1:(state.jmax[l] - 1)
        jp1 = j + 1

        # --- Predictor: explicit Euler across node j ---
        diss_j = state.dfsta[j] + state.dbsta[j]
        if config.options.iperm == 1; diss_j += state.dpsta[j]; end
        if iveg == 3; diss_j += state.dvegsta[j]; end
        dum = diss_j * state.wt[j]
        dum = (state.efsta[j] - dx * dum) / state.fe
        if dum ≤ 0.0
            jp1 = j   # "accept up to node j"
            march_aborted = true
            break
        end
        sigite = sqrt(dum)

        # Roller predictor: forward Euler
        rqite = 0.0
        if iroll == 1
            rqite = state.rq[j] + dx * (state.dbsta[j] - state.rbeta[j] * state.rq[j]) / state.re[j]
            rqite = max(0.0, rqite)
        end

        state.sxxsta[jp1] = state.fsx * sigite^2
        if iangle == 1
            state.sxysta[jp1] = state.fsy * sigite^2
        end
        # Roller adds to radiation stress
        if iroll == 1
            state.sxxsta[jp1] += state.rx[jp1] * rqite
            if iangle == 1
                state.sxysta[jp1] += state.ry[jp1] * rqite
            end
        end
        # Momentum equation: WSETUP at JP1
        # Wind stress (twxsta) adds to the x-momentum balance.
        # Clamp h[j] denominator to eps1 to prevent setup blow-up at near-zero depth.
        state.wsetup[jp1] = state.wsetup[j] -
            (state.sxxsta[jp1] - state.sxxsta[j] +
             (state.tbxsta[j] + state.twxsta) * dx) / max(state.h[j], eps1)
        hite = state.wsetup[jp1] + state.swldep[jp1, l]

        if hite < eps1
            jp1 = j
            march_aborted = true
            break
        end

        qdisp_jp1 = qdisp_at(jp1)
        lwave!(state, config, jp1, l, hite, wave; qdisp=qdisp_jp1, hrms_j=sigite*SQR8)

        # Spatially-varying friction (IFRICTION_SPATIAL=1): update fb2 at jp1
        # using the predictor orbital velocity before the corrector loop so that
        # both tbxsta and dfsta see a consistent Shields-regime friction value.
        # sigt_pre = cp · (sigite/hite) is the wave orbital velocity scale;
        # multiply by ctheta for the cross-shore component (normal incidence: =1).
        if config.options.ifriction_spatial == 1
            sigt_pre = hite > eps1 ? state.cp[jp1] * sigite / hite : 0.0
            ustd_pre = sigt_pre * state.ctheta[jp1]
            _update_fb2_spatial!(state, config, jp1, l, ustd_pre)
        elseif config.options.ifriction_spatial == 2
            _update_fb2_manning!(state, config, jp1, l)
        end

        # --- Corrector iteration ---
        converged = false
        ite_aborted = false
        for ite in 1:maxite
            hrmite = sigite * SQR8
            dbreak!(state, config, jp1, l, hrmite, hite, state.wt[jp1])
            # Porous flow in corrector
            if config.options.iperm == 1
                pkhsig = state.wn[jp1] * sigite
                dedx = (state.wsetup[jp1] - state.wsetup[j]) / dx
                porous_flow!(state, config, jp1, l, pkhsig, dedx)
            end
            state.sigsta[jp1] = sigite / hite
            if state.sigsta[jp1] > config.sismax
                state.sigsta[jp1] = config.sismax
            end
            sigt = state.cp[jp1] * state.sigsta[jp1]

            if state.fb2[jp1, l] > 0.0
                dum2 = GRAV * hite / state.cp[jp1]^2
                usigt = -state.ctheta[jp1] * state.sigsta[jp1] * dum2
                gbx, gf = friction_coefficients(state.ctheta[jp1], usigt,
                                                state.stheta[jp1], 0.0, iangle)
                state.gbx[jp1] = gbx
                state.gf[jp1]  = gf
                state.tbxsta[jp1] = state.fb2[jp1, l] * gbx * sigt^2 / GRAV
                state.dfsta[jp1]  = state.fb2[jp1, l] * gf  * sigt^3 / GRAV
                # Vegetation friction enhancement (IVEG=1,2)
                if iveg in (1, 2)
                    apply_veg_friction!(state, config, jp1, l)
                end
            else
                usigt = 0.0
                state.tbxsta[jp1] = 0.0
                state.dfsta[jp1]  = 0.0
            end

            # Wave-current: update qwx from latest usigt.
            # Disable in very shallow water to prevent spurious WCI amplification.
            if iwcint == 1
                state.qwx = hite > config.min_depth_wcint ? usigt * sigt * hite : 0.0
                # ITIDE=1 + ILAB=0: add tidal cross-shore flux
                if config.options.itide == 1 && config.options.ilab == 0
                    state.qwx += state.qtide[jp1]
                end
            end

            # Vegetation dissipation (IVEG=3)
            if iveg == 3
                veg_dissipation!(state, config, jp1, l)
            end

            # Trapezoidal energy flux
            # IVEG=3: add vegetation dissipation to total
            dumd = state.dfsta[jp1] + state.dfsta[j] + state.dbsta[jp1] + state.dbsta[j]
            if config.options.iperm == 1
                dumd += state.dpsta[jp1] + state.dpsta[j]
            end
            if iveg == 3
                dumd += state.dvegsta[jp1] + state.dvegsta[j]
            end
            dumd *= 0.5 * (state.wt[j] + state.wt[jp1])
            dum  = (state.efsta[j] - dxd2 * dumd) / state.fe
            if dum ≤ 0.0
                jp1 = j
                ite_aborted = true
                break
            end
            state.sigma[jp1] = sqrt(dum)

            # Roller corrector: trapezoidal/implicit
            if iroll == 1
                dum1 = state.re[jp1] + dxd2 * state.rbeta[jp1]
                dum2 = (state.re[j] - dxd2 * state.rbeta[j]) * state.rq[j] +
                       dxd2 * (state.dbsta[jp1] + state.dbsta[j])
                state.rq[jp1] = max(0.0, dum2 / dum1)
            end

            # Stream stress + momentum equation
            state.sxxsta[jp1] = state.fsx * state.sigma[jp1]^2
            if iangle == 1
                state.sxysta[jp1] = state.fsy * state.sigma[jp1]^2
            end
            # Roller adds to radiation stress
            if iroll == 1
                state.sxxsta[jp1] += state.rx[jp1] * state.rq[jp1]
                if iangle == 1
                    state.sxysta[jp1] += state.ry[jp1] * state.rq[jp1]
                end
            end
            if iveg == 3
                # IVEG=3: vegetation-enhanced momentum.
                # Stream stress includes vegetation dissipation (IFV=1).
                stream_stress = state.fsx * (state.dfsta[jp1] + state.dfsta[j] +
                                             state.dvegsta[jp1] + state.dvegsta[j]) /
                                (state.wn[jp1] * state.cp[jp1] + state.wn[j] * state.cp[j])
                # Veg-modified bottom stress: (1 + vegcdm/vegcd * min(vegh,h) * vegfb) * tbxsta
                tb_jp1 = veg_momentum_multiplier(state, config, jp1, l) * state.tbxsta[jp1]
                tb_j   = veg_momentum_multiplier(state, config, j, l)   * state.tbxsta[j]
                state.wsetup[jp1] = state.wsetup[j] -
                    (2.0 * (state.sxxsta[jp1] - state.sxxsta[j]) +
                     dx * (tb_jp1 + tb_j + 2.0 * (stream_stress + state.twxsta))) /
                    max(hite + state.h[j], eps1)
            else
                # Standard (IVEG=0,1,2) momentum
                stream_stress = state.fsx * (state.dfsta[jp1] + state.dfsta[j]) /
                                (state.wn[jp1] * state.cp[jp1] + state.wn[j] * state.cp[j])
                state.wsetup[jp1] = state.wsetup[j] -
                    (2.0 * (state.sxxsta[jp1] - state.sxxsta[j]) +
                     dx * (state.tbxsta[jp1] + state.tbxsta[j] +
                           2.0 * (stream_stress + state.twxsta))) /
                    max(hite + state.h[j], eps1)
            end

            state.h[jp1] = state.wsetup[jp1] + state.swldep[jp1, l]
            state.sigsta[jp1] = state.sigma[jp1] / state.h[jp1]
            if state.sigsta[jp1] > config.sismax
                state.sigsta[jp1] = config.sismax
            end
            if state.h[jp1] ≤ eps1
                jp1 = j
                ite_aborted = true
                break
            end

            qdisp_jp1 = qdisp_at(jp1)
            lwave!(state, config, jp1, l, state.h[jp1], wave; qdisp=qdisp_jp1, hrms_j=sigite*SQR8)

            # Convergence check
            esigma = abs(state.sigma[jp1] - sigite)
            eh     = abs(state.h[jp1] - hite)
            erq    = iroll == 1 ? abs(state.rq[jp1] - rqite) : 0.0
            if esigma < eps1 && eh < eps1 && erq < eps2
                converged = true
                break
            end
            # Average to accelerate convergence
            sigite = 0.5 * (state.sigma[jp1] + sigite)
            hite   = 0.5 * (state.h[jp1] + hite)
            if iroll == 1
                rqite = 0.5 * (state.rq[jp1] + rqite)
            end
        end
        # Non-convergence: adopt last-iteration values and terminate the march.
        if ite_aborted
            march_aborted = true
            break
        end

        state.hrms[jp1]  = SQR8 * state.sigma[jp1]
        state.wsetup[jp1] = state.h[jp1] - state.swldep[jp1, l]

        # Final LWAVE + DBREAK at the converged state
        qdisp_jp1 = qdisp_at(jp1)
        lwave!(state, config, jp1, l, state.h[jp1], wave; qdisp=qdisp_jp1, hrms_j=state.hrms[jp1])
        dbreak!(state, config, jp1, l, state.hrms[jp1], state.h[jp1], state.wt[jp1])
        # Final vegetation dissipation (IVEG=3) at the converged state
        if iveg == 3
            veg_dissipation!(state, config, jp1, l)
        end
        state.sigsta[jp1] = state.sigma[jp1] / state.h[jp1]
        if state.sigsta[jp1] > config.sismax
            state.sigsta[jp1] = config.sismax
        end
        sigt = state.sigsta[jp1] * state.cp[jp1]

        # Refresh dynamic Manning fb2 from the corrector-converged depth before
        # the final tbxsta/dfsta write.
        if config.options.ifriction_spatial == 2
            _update_fb2_manning!(state, config, jp1, l)
        end

        # Re-store radiation-stress / energy-flux arrays at jp1 from final state
        sigma2 = state.sigma[jp1]^2
        state.sxxsta[jp1] = sigma2 * state.fsx
        state.efsta[jp1]  = sigma2 * state.fe
        if iangle == 1
            state.sxysta[jp1] = sigma2 * state.fsy
        end
        # Roller: final Sxx/Sxy update at converged state
        if iroll == 1
            state.sxxsta[jp1] += state.rx[jp1] * state.rq[jp1]
            if iangle == 1
                state.sxysta[jp1] += state.ry[jp1] * state.rq[jp1]
            end
        end
        if state.fb2[jp1, l] > 0.0
            usigt = -state.ctheta[jp1] * state.sigsta[jp1] * state.h[jp1] /
                    (state.cp[jp1]^2 / GRAV)
            gbx, gf = friction_coefficients(state.ctheta[jp1], usigt,
                                            state.stheta[jp1], 0.0, iangle)
            state.gbx[jp1] = gbx
            state.gf[jp1]  = gf
            state.tbxsta[jp1] = state.fb2[jp1, l] * gbx * sigt^2 / GRAV
            state.dfsta[jp1]  = state.fb2[jp1, l] * gf  * sigt^3 / GRAV
            # Vegetation friction enhancement (IVEG=1,2)
            if iveg in (1, 2)
                apply_veg_friction!(state, config, jp1, l)
            end
        else
            state.tbxsta[jp1] = 0.0
            state.dfsta[jp1]  = 0.0
        end

        if state.h[jp1] < eps1 || jp1 == jdum
            break
        end
    end

    # Runup limit
    state.jr = jp1
    state.xr = state.xb[jp1]
    state.zr = state.zb[jp1, l]

    # Velocity moments
    @inbounds for i in 1:jp1
        sigt = state.cp[i] * state.sigsta[i]
        state.ustd[i]  = sigt * state.ctheta[i]
        state.umean[i] = -state.ustd[i] * state.sigsta[i] * GRAV * state.h[i] /
                         state.cp[i]^2
        # Roller modifies mean current
        if iroll == 1 && state.sigma[i] > 1e-10
            state.umean[i] *= (1.0 + (state.cp[i] / GRAV) * state.rq[i] / state.sigma[i]^2)
        end
        state.usta[i]  = sigt > 1e-10 ? min(state.umean[i] / sigt, 1.0) : 0.0

        # Wave nonlinearity: always store Ursell; compute Sk/As by chosen
        # closure family. Method comes from WaveNonlinearityConfig.closure;
        # fanned out from legacy AsymmetryConfig / iasym int when needed.
        # Note: use wkp_arr[i] (local wavenumber) — NOT wn[i] (Cg/Cp ratio).
        state.ursell[i] = ursell_number(state.hrms[i], state.wkp_arr[i], state.h[i])
        nl_d = nonlinearity(config)
        sk_method = nl_d.closure
        sk, as_ = if sk_method == :stokes2 || sk_method == :boussinesq_diss
            # Sk from Stokes2 (zero-calibration); As filled by Ruessink as a
            # fallback so :stokes2 still has a non-zero As value if the user
            # has asymmetry > 0. For :boussinesq_diss As gets overwritten in
            # a second pass below using the dissipation gradient.
            sk_s2, _ = stokes2_skewness(state.wkp_arr[i], state.h[i], state.hrms[i])
            _, as_r = ruessink_skewness_asymmetry(state.ursell[i])
            (sk_s2, sk_method == :boussinesq_diss ? 0.0 : as_r)
        elseif sk_method == :linear
            (0.0, 0.0)
        else
            # :ruessink (default) — Ursell-based Sk and As
            ruessink_skewness_asymmetry(state.ursell[i])
        end

        # ── Biphase relaxation along the profile ──────────────────────────
        # When biphase_relax_length > 0 (and the closure is :ruessink), the
        # equilibrium biphase β_eq(Ur) is not applied instantaneously;
        # instead it relaxes spatially:
        #   β[i] = β[i-1] + (β_eq[i] - β[i-1]) * dx / L_relax
        # The amplitude B is taken from the equilibrium value at each node.
        if nl_d.biphase_relax_length > 0.0 && sk_method == :ruessink
            ur_i   = state.ursell[i]
            if ur_i > 0.0
                # Extract equilibrium B and β from Ruessink
                p1, p2, p3, p4, p5, p6 = 0.0, 0.857, -0.471, 0.297, 0.815, 0.672
                B_eq   = p1 + (p2 - p1) / (1.0 + exp((p3 - log(ur_i)) / p4))
                β_eq   = -π / 2.0 * (1.0 - tanh(p5 / ur_i^p6))
                β_prev = i == 1 ? β_eq : state.biphase[i - 1]
                α      = min(1.0, dx / nl_d.biphase_relax_length)
                β_i    = β_prev + α * (β_eq - β_prev)
                state.biphase[i] = β_i
                sk  = B_eq * cos(β_i)
                as_ = B_eq * sin(β_i)
            else
                state.biphase[i] = 0.0
            end
        else
            state.biphase[i] = sk_method == :ruessink ?
                atan(as_, sk) : 0.0   # store static biphase for diagnostics
        end

        state.skewness[i]  = sk
        state.asymmetry[i] = as_
        if iangle == 1
            state.vstd[i] = sigt * abs(state.stheta[i])
            # Longshore momentum balance:
            #   tbysta = -dSxy/dx + twysta - h*detady  (+ veg attenuation)
            #
            # Two branches based on whether wave action is significant:
            # ─── wave-driven (sigt > SIGT_NOWAVE) ──────────────────────
            #   Use the FORTRAN-parity Gaussian-friction closure (vstgby).
            #   gby = tbysta / (fb2 * cp²/g * sigsta²); vmean = vstgby·sigt.
            # ─── current-only (sigt ≤ SIGT_NOWAVE) ─────────────────────
            #   Wave bottom shear → 0, so the Gaussian closure degenerates
            #   (vmean = vsigt·sigt → 0). Fall back to plain quadratic-
            #   friction open-channel closure:  tbysta = (fb2/2)·V·|V|·vegcv
            #   ⇒ V = sign(tbysta)·√(2·|tbysta| / (fb2·vegcv)).
            #   This makes ICURRENT=1 / ITIDE=1 work in fluvial / tidal-
            #   channel cases with no waves at the boundary.
            #
            # Threshold: 1e-3 m/s of wave-RMS orbital velocity. Below this,
            # the wave contribution to bed shear is negligible compared to
            # any imposed/longshore current.
            SIGT_NOWAVE = 1e-3
            wave_driven = sigt > SIGT_NOWAVE

            if i > 1 && state.fb2[i, l] > 0.0
                # Common forcing terms (wave radstress + wind)
                dsxy = state.sxysta[i] - state.sxysta[i - 1]
                if state.stheta[i] * dsxy > 0.0; dsxy = 0.0; end
                tbysta_wave_wind = -dsxy / dx + state.twysta

                # Vegetation attenuation factor (1.0 when no vegetation)
                vegcv = 1.0
                if config.options.iveg >= 1 && config.vegetation !== nothing
                    veg = config.vegetation
                    vh = size(state.vegh, 1) > 0 ? state.vegh[i, l] : veg.vegh[i, l]
                    dumh = min(vh, state.h[i])
                    vegcv = 1.0 + dumh * veg.vegfb[i, l]
                end

                # ICURRENT=1: back-solve detady_now at the first valid cell
                # using whichever closure is appropriate for this node.
                if config.options.icurrent == 1 && !state.icurrent_solved[l]
                    tbysta_target = if wave_driven
                        # Wave-driven Gaussian closure inversion
                        vsigt_target = state.vbc_now / sigt
                        gby_target = gby_from_vsigt(state.ctheta[i],
                                                    state.stheta[i],
                                                    vsigt_target)
                        dum3_ref = state.cp[i]^2 / GRAV
                        gby_target * vegcv *
                            state.fb2[i, l] * dum3_ref * state.sigsta[i]^2
                    else
                        # Current-only quadratic-friction inversion:
                        # tbysta = (fb2/2) · V·|V| · vegcv
                        0.5 * state.fb2[i, l] * vegcv *
                            state.vbc_now * abs(state.vbc_now)
                    end
                    h_ref = state.h[i]
                    if h_ref > 1e-6
                        state.detady_now = (tbysta_wave_wind - tbysta_target) / h_ref
                        state.icurrent_solved[l] = true
                    end
                end

                # Forward solve at this node
                tbysta_i = tbysta_wave_wind
                if config.options.itide == 1 || config.options.icurrent == 1
                    tbysta_i -= state.h[i] * state.detady_now
                end
                state.tbysta[i] = tbysta_i

                if wave_driven
                    # Existing FORTRAN-parity closure
                    dum3 = state.cp[i]^2 / GRAV
                    gby_raw = tbysta_i / (state.fb2[i, l] * dum3 * state.sigsta[i]^2)
                    if vegcv != 1.0
                        gby_raw /= vegcv
                    end
                    state.gby[i] = gby_raw
                    usigt = state.usta[i]
                    vsigt = longshore_vstgby(state.ctheta[i], usigt,
                                             state.stheta[i], state.gby[i])
                    state.vmean[i] = vsigt * sigt
                    state.vsta[i] = vsigt
                else
                    # Wave-bypass: open-channel quadratic-friction closure.
                    # tbysta = (fb2/2)·V·|V|·vegcv ⇒ V = sign(tbysta)·√(2·|tbysta|/(fb2·vegcv))
                    cf2 = state.fb2[i, l] * vegcv
                    if cf2 > 1e-12 && abs(tbysta_i) > 1e-15
                        state.vmean[i] = sign(tbysta_i) *
                                         sqrt(2.0 * abs(tbysta_i) / cf2)
                    else
                        state.vmean[i] = 0.0
                    end
                    state.gby[i]  = 0.0   # not meaningful here
                    state.vsta[i] = 0.0   # wave-normalized velocity n/a
                end
            else
                state.gby[i] = 0.0
                state.tbysta[i] = 0.0
                state.vmean[i] = 0.0
                state.vsta[i] = 0.0
            end

            # Longshore water flux: QWY += RQ*sin(theta)
            if iroll == 1
                state.qwy += state.rq[i] * state.stheta[i]
            end
        end

        # ICURRENT=1 with iangle=0: propagate imposed current as spatially-uniform
        # vmean so iv_transport=1 can use it for enhanced bed-shear entrainment.
        # The longshore momentum balance (Sxy-driven) is not solved at normal
        # incidence, but the imposed vbc still enhances threshold exceedance and
        # suspension probability in the transport kernels.
        if iangle == 0 && config.options.icurrent == 1
            state.vmean[i] = state.vbc_now
        end
    end

    # ── Boussinesq asymmetry from breaker-dissipation gradient ──────────────
    # When closure == :boussinesq_diss, override state.asymmetry with the
    # Doering & Bowen (1995) / Reniers et al. (2004) closure:
    #
    #   As(j) = -K_as · (∂Db/∂x · h) / (g · max(Hrms², ε) · ω)
    #
    # The breaker-dissipation gradient is naturally peaked across the breaking
    # transition zone and changes sign across the bar crest, capturing the
    # acceleration-skewness signature that the Ursell-based Ruessink form
    # misses. Sk has already been set above (Stokes2). The result is clamped
    # to ±diss_grad_cap to prevent runaway in numerically rough zones.
    nl_bd = nonlinearity(config)
    if nl_bd.closure == :boussinesq_diss && jp1 ≥ 3
        K_as = nl_bd.K_as
        cap  = nl_bd.diss_grad_cap
        win  = max(1, nl_bd.smooth_diss_window)

        # Optional pre-smooth Db with a centered boxcar (helps numerical noise)
        Db_smooth = similar(state.dbsta, jp1)
        @inbounds for i in 1:jp1
            i0 = max(1, i - win); i1 = min(jp1, i + win)
            s = 0.0; n = 0
            for k in i0:i1
                s += state.dbsta[k]; n += 1
            end
            Db_smooth[i] = n > 0 ? s / n : 0.0
        end

        @inbounds for i in 2:(jp1-1)
            # Centered gradient (∂Db/∂x). Positive on the seaward face of
            # the bar (Db rising toward break point), negative beyond.
            dDbdx = (Db_smooth[i+1] - Db_smooth[i-1]) / (2.0 * dx)
            # ω = 2π / Tp from current BC window (Tp is window-constant).
            ωj = wave.tp > 0.0 ? (2.0 * π / wave.tp) : 0.0
            ωj <= 0.0 && continue
            Hrms2 = max(state.hrms[i] * state.hrms[i], 1e-6)
            h_j   = max(state.h[i], 0.05)
            # Dimensionless As proxy. Sign convention: positive ∂Db/∂x
            # (waves breaking → dissipation rising) → As < 0 (front-leaning,
            # acceleration-skewed) which drives onshore bedload via the
            # Ruessink-style coupling. Hence the leading minus.
            As_j = -K_as * dDbdx * h_j / (GRAV * Hrms2 * ωj)
            state.asymmetry[i] = clamp(As_j, -cap, cap)
        end
        # Boundary nodes: copy adjacent value to avoid spurious endpoints.
        state.asymmetry[1]   = state.asymmetry[2]
        state.asymmetry[jp1] = state.asymmetry[jp1-1]
    end

    # Smooth h, sigma, ustd, umean, usta, dfsta, dbsta
    if config.options.ismooth == 1 && jp1 ≥ 3
        smooth_tridiagonal!(view(state.h,      1:jp1))
        smooth_tridiagonal!(view(state.sigma,  1:jp1))
        smooth_tridiagonal!(view(state.ustd,   1:jp1))
        smooth_tridiagonal!(view(state.umean,  1:jp1))
        smooth_tridiagonal!(view(state.usta,   1:jp1))
        smooth_tridiagonal!(view(state.dfsta,  1:jp1))
        smooth_tridiagonal!(view(state.dbsta,  1:jp1))
        # Smooth the Boussinesq-derived As alongside dbsta (consistent with
        # the smoothing applied to its source field).
        if nl_bd.closure == :boussinesq_diss
            smooth_tridiagonal!(view(state.asymmetry, 1:jp1))
        end
        if iveg == 3
            smooth_tridiagonal!(view(state.dvegsta, 1:jp1))
        end
    end

    return state
end
