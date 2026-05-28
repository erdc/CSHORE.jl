# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
clay_dike_erosion.jl — Direct port of the FORTRAN CSHORE EROSON subroutine.

Two modes:

- `IPROFL=2`: grassed-dike erosion. Uses its own time-stepping (DELT from
  DELEM/DMAX stability criterion) and directly updates the bed elevation
  from cumulative forcing DSUM.

- `ICLAY=1`: sand layer over erodible clay. Uses DELT from the Exner
  solver. When the sand layer is thinner than D50, the clay below erodes
  at a rate set by the spatially-varying erodibility `rclay` and the
  combined wave+friction forcing DSTA.

The forcing DSTA = DEEB·DBSTA + DEEF·DFSTA (wave breaking + friction
dissipation, both with efficiency factors) amplified by a bed-slope
factor BSF that tracks local steepness. In the swash zone the forcing
comes from DFSWD with a separate calibration scheme.

Current scope: everything inside EROSON (both ICLAY=1 and IPROFL=2)
except IPOND=1 (ponded-water interactions) and IWTRAN=1 (landward
water-body transmission).
==============================================================================#

"""
    gfdwd(r) -> Float64

Dike-erosion shear-stress function for the wet-dry swash zone.
"""
@inline function gfdwd(r::Float64)
    tr = 3.0 * r
    r2 = r * r
    r3 = r2 * r
    if r >= 0.0
        return 1.32934 + tr + 2.658681 * r2 + r3
    else
        return 1.32934 * (1.0 + 2.0 * r2) * (2.0 * erfcc(r) - 1.0) -
               tr - r3 + (16.0 * r3 + 9.0 * r) * exp(-r2)
    end
end

"""
    initialize_erosion!(state, config, l)

Compute time-invariant EROSON constants at `time == 0` for line `l`. Sets:
- `fba3[j,l] = fb2[j,l] · √g · AWD³`
- For ICLAY=1: `epclay = 0`, `bsf = 1`
- For IPROFL=2: `grs1..grs5` from grass resistance parameters; `dsum = 0`
"""
function initialize_erosion!(state::CshoreState, config::CshoreConfig, l::Int)
    jmax = state.jmax[l]
    awd = config.swash.awd > 0 ? config.swash.awd : 1.6
    sqrg = sqrt(GRAV)
    @inbounds for j in 1:jmax
        state.fba3[j, l] = state.fb2[j, l] * sqrg * awd^3
    end

    if config.options.iclay == 1
        @inbounds for j in 1:jmax
            state.epclay[j, l] = 0.0
            state.bsf[j] = 1.0
        end
    else
        # IPROFL=2 dike mode: derived constants from grass resistance
        @inbounds for j in 1:jmax
            grsd_j  = state.grsd[j, l]
            grsr_j  = state.grsr[j, l]
            grsrd_j = state.grsrd[j, l]
            if grsrd_j > 0.0
                state.grs3[j, l] = GRAV / grsrd_j
            end
            state.grs4[j, l] = 0.5 * grsd_j * (grsr_j - grsrd_j) / max(grsrd_j, 1e-12)
            state.grs5[j, l] = 0.5 * grsd_j * (grsr_j + grsrd_j) / GRAV
            if grsd_j <= 0.0
                state.grs1[j, l] = 0.0
                state.grs2[j, l] = 0.0
            else
                dum = grsr_j - grsrd_j
                if abs(dum) > 1e-12
                    state.grs1[j, l] = grsd_j * grsr_j / dum
                    state.grs2[j, l] = 2.0 * GRAV * dum / grsd_j / grsr_j^2
                end
            end
            state.dsum[j] = 0.0
        end
    end
    state.eroson_initialized = true
    return state
end

"""
    eroson!(state, config, l) -> (iend, dt_used)

Computes:
1. BSF (bed-slope amplification) via smoothed inverse-slope function
2. DSTA (wet zone) and DFSWD (swash) with ED calibration at JWD
3. Connect DSTA/DFSWD via `tranwd!` and smooth
4. Branch:
   - IPROFL=2: compute DELT from DELEM/DMAX, integrate DSUM, update ZB from EDIKE
   - ICLAY=1: use existing `state.delt` from Exner, apply HP<D50 gate to
     update EPCLAY, ZP, and reduce ZB by sand fraction in clay

Returns `(iend, dt_used)`:
- `iend=true` when the computed DELT hits the BC window boundary (IPROFL=2 only)
- `dt_used` is the time step actually applied (identical to state.delt for ICLAY=1)
"""
function eroson!(state::CshoreState, config::CshoreConfig, l::Int;
                 t_window_end::Float64=0.0)
    if !state.eroson_initialized
        initialize_erosion!(state, config, l)
    end

    jmax = state.jmax[l]
    jr   = state.jr
    jwd  = state.jwd
    jdry = state.jdry == 0 ? jr : state.jdry
    iroll = config.options.iroll

    # Erosion efficiencies (different sources for dike vs clay)
    deeb, deef = if config.options.iclay == 1
        config.clay.deeb, config.clay.deef
    elseif config.options.iprofl == 2
        config.dike.deeb, config.dike.deef
    else
        0.005, 0.005
    end

    # ---- 1. Bed-slope amplification factor BSF (IPROFL=2 only) ----
    # For ICLAY=1 BSF stays at 1.0 as set in initialize_erosion!
    if config.options.iprofl == 2
        scp = 1.2       # max clay slope (FORTRAN DATA)
        dumvec = zeros(jdry)
        @inbounds for j in 1:jdry
            asb = abs(state.bslope[j, l])
            dum = asb / scp
            dumvec[j] = dum >= 0.9 ? 10.0 : 1.0 / (1.0 - dum)
        end
        bsf_slice = @view state.bsf[1:jdry]
        bsf_slice .= dumvec
        smooth_tridiagonal!(bsf_slice)
    end

    # ---- 2. DSTA in wet zone ----
    @inbounds for j in 1:jr
        d = if iroll == 1
            deeb * state.rbeta[j] * state.rq[j] + deef * state.dfsta[j]
        else
            deeb * state.dbsta[j] + deef * state.dfsta[j]
        end
        state.dsta[j] = d * state.bsf[j]
    end

    # ---- 3. DFSWD in swash zone with ED calibration ----
    aqwd = config.swash.aqwd > 0 ? config.swash.aqwd :
           (0.75 * sqrt(π) * (config.swash.awd > 0 ? config.swash.awd : 1.6) *
            (config.swash.awd > 0 ? config.swash.awd : 1.6))
    cwd = config.swash.cwd > 0 ? config.swash.cwd :
          (0.75 * sqrt(π) * (config.swash.awd > 0 ? config.swash.awd : 1.6))
    ed = 1.0
    @inbounds for j in jwd:jdry
        hj = max(state.h[j], 1e-12)
        pw = max(state.pwet[j], 1e-12)
        dum = aqwd * hj * sqrt(GRAV * hj / pw)
        rs = dum < 1e-6 ? 0.0 : cwd * (state.qo[l] - dum) / dum
        state.dfswd[j] = ed * state.fba3[j, l] * hj * sqrt(hj / pw) * gfdwd(rs)
        state.dfswd[j] *= state.bsf[j]
        if j == jwd && state.dfswd[j] > 1e-12
            ed = state.dsta[j] / state.dfswd[j]
            state.dfswd[j] = state.dsta[j]
        end
    end

    # ---- 4. Connect and smooth ----
    if jdry > jr
        tranwd!(view(state.dsta, :), jr, view(state.dfswd, :), jwd, jdry)
    else
        jdry = jr
    end
    smooth_tridiagonal!(view(state.dsta, 1:jdry))
    if jdry < jmax
        extrapolate_boundary!(view(state.dsta, jdry:jmax))
    end

    # ---- 5. Branch on mode ----
    iend = false
    dt_used = state.delt

    if config.options.iprofl == 2
        # Compute DELT internally
        delem = 0.05    # max dike-erosion increment (m) per step
        dmax = 1e-6
        @inbounds for j in 1:jmax
            dum = state.dsta[j] * state.grs3[j, l]
            if dum > dmax; dmax = dum; end
        end
        delt_local = delem / dmax
        dt_bc = t_window_end - state.time
        if delt_local > 0.5 * dt_bc; delt_local = 0.5 * dt_bc; end
        if state.time + delt_local >= t_window_end
            delt_local = t_window_end - state.time
            iend = true
        end
        state.delt = delt_local
        dt_used = delt_local

        # Integrate cumulative erosion
        @inbounds for j in 1:jmax
            state.dsum[j] += delt_local * state.dsta[j]
            grsd_j = state.grsd[j, l]
            edike_j = 0.0
            if grsd_j > 0.0
                if state.dsum[j] < state.grs5[j, l]
                    disc = 1.0 - state.grs2[j, l] * state.dsum[j]
                    edike_j = state.grs1[j, l] * (1.0 - sqrt(max(disc, 0.0)))
                else
                    edike_j = state.grs3[j, l] * state.dsum[j] - state.grs4[j, l]
                end
            else
                edike_j = state.grs3[j, l] * state.dsum[j]
            end
            state.edike[j, l] = edike_j
            state.zb[j, l] = state.zb0[j, l] - edike_j
        end
    end

    if config.options.iclay == 1
        # Clay erosion below sand layer — uses state.delt from the Exner solver
        d50 = config.sediment.d50
        @inbounds for j in 1:jmax
            if state.hp[j, l] < d50
                dum = state.delt * state.rclay[j, l] * state.dsta[j]
                state.epclay[j, l] += dum
                state.zp[j, l] = state.zp0_clay[j, l] - state.epclay[j, l]
                state.zb[j, l] -= dum * state.fclay[j, l]
            end
        end
    end

    return iend, dt_used
end
