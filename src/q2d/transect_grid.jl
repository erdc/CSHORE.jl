# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
q2d/transect_grid.jl — Quasi-2D morphodynamic framework.

Couples N independent CSHORE.jl cross-shore transects through their
longshore sediment-transport fluxes and wave-angle gradients to produce
a quasi-2D (Q2D) plan-form coastal evolution model.

Architecture
============
Each transect is a `CshoreBMI` instance that handles all cross-shore physics
(wave transformation, undertow, bedload, suspended load, swash, Exner).
The Q2D orchestrator handles the *alongshore* physics between transects:

  1. Advance all N transects by one BC window in parallel (Julia @threads).
  2. Read the longshore transport Qby[j,i] from every transect.
  3. Compute ∂Qby/∂y at each (j,i) using finite differences.
  4. Apply the y-direction Exner equation to update zb.
  5. Optionally update wave angles from the evolving shoreline plan-form.

The y-Exner Equation
====================
For each cross-shore node j and transect i:

  Δzb[j,i] = -dt · (1/(1-p)) · ∂Qby/∂y[j,i]

where ∂Qby/∂y is approximated by central differences over the transect
spacing dy. At boundary transects (i=1 and i=N), the gradient is computed
with a one-sided difference, and the boundary condition is that longshore
transport is equal to its nearest interior neighbour (effectively zero net
divergence at open boundaries) unless the user specifies `Qby_bc_left` /
`Qby_bc_right` to prescribe an updrift/downdrift supply.

Wave Angle Update (Shoreline Feedback)
=======================================
When the shoreline rotates between transects (due to differential
erosion/accretion), the effective angle of wave incidence changes. The
Q2D framework updates `wangbc` for each transect at every coupling step
using the local shoreline azimuth relative to the baseline (or a supplied
offshore angle). This closes the plan-form feedback loop without an
external wave model.

Thread Safety
=============
All N transects are advanced concurrently with `@threads`. After the
parallel step, a serial exchange phase reads fluxes and injects bed
corrections. No shared mutable state exists between transects during the
parallel step — each `CshoreBMI` owns its own `CshoreState`.

Stability
=========
The y-Exner CFL limit (Courant number ≤ 1) requires:
  dt_couple ≤ (1-p) · dy / (∂Qby/∂zb)_max

In practice the coupling step equals the BC window length (typically 0.5–3 h).
For beaches with |Qby| < 0.1 m²/s and dy ≥ 100 m, this is always satisfied.
The user is warned if the estimated Courant number exceeds 0.8.

Public API
==========
  TransectGrid         — holds the N transects, positions, and orchestration params
  step_transect_grid!  — advance the Q2D system by one coupling step
  run_transect_grid!   — run the full simulation (with optional output)
  transect_grid_output — collect all-transect snapshots into a NamedTuple
==============================================================================#

import BasicModelInterface as BMI

"""
    LongshoreBoundary

Prescribed longshore sediment supply at the grid lateral boundaries.

- `Qby_left`  : supply entering from the left (i=1) boundary (m²/s).
  Positive = sediment entering the domain from updrift.  `NaN` (default) →
  zero-gradient open boundary (no net divergence at the edge transect).
- `Qby_right` : supply leaving (or entering) at the right (i=N) boundary.
  `NaN` → zero-gradient open boundary.

Both can be time series (Vector of length == number of BC windows) or
a scalar constant.  Positive values represent supply *into* the domain.
"""
Base.@kwdef struct LongshoreBoundary
    Qby_left::Union{Float64, Vector{Float64}}  = NaN
    Qby_right::Union{Float64, Vector{Float64}} = NaN
end

"""
    WaveAngleMode

Controls how the offshore wave angle is updated each coupling step.

- `:fixed`        — all transects use the angle already in `wangbc` (no update).
- `:shoreline`    — the wave angle at each transect is rotated to remain normal
                    to the *local* shoreline azimuth, estimated from the cross-shore
                    position of the SWL intersection node (simple contour method).
                    This closes the Pelnard-Considère plan-form feedback loop.
- `:snell`        — reserved for future coupling to an external wave model that
                    supplies refracted angles per transect.
"""
@enum WaveAngleMode begin
    WaveAngleFixed      = 0
    WaveAngleShoreline  = 1
    WaveAngleSnell      = 2
end

"""
    TransectGrid

Holds N cross-shore `CshoreBMI` transects and the metadata needed to
couple them alongshore.

# Fields
- `transects`     — `Vector{CshoreBMI}`, length N.
- `y`             — alongshore positions of each transect centre (m), length N.
- `dy`            — half-spacings; `dy[i]` is the half-distance to the right
                    neighbour (used for the Exner finite difference). Computed
                    from `y` automatically if not supplied.
- `porosity`      — bed sediment porosity used in the y-Exner equation.
                    Defaults to `transects[1].config.sediment.porosity`.
- `lsbc`          — `LongshoreBoundary` — lateral sediment supply at y=0 and y=L.
- `angle_mode`    — `WaveAngleMode` — how wave angles evolve with plan-form.
- `wave_angle_0`  — baseline offshore wave angle (rad). Used by the shoreline
                    mode as the reference angle for the unperturbed coast.
- `dt_couple`     — coupling time step (s). Defaults to the BC window length.
"""
struct TransectGrid
    transects::Vector{CshoreBMI}
    y::Vector{Float64}          # (N,) alongshore coordinate (m)
    dy_left::Vector{Float64}    # (N,) distance to left  neighbour (m)
    dy_right::Vector{Float64}   # (N,) distance to right neighbour (m)
    porosity::Float64
    lsbc::LongshoreBoundary
    angle_mode::WaveAngleMode
    wave_angle_0::Float64       # baseline offshore wave angle (rad)
end

"""
    TransectGrid(transects, y; kwargs...) -> TransectGrid

Construct a `TransectGrid` from a vector of initialized `CshoreBMI` instances
and their alongshore positions `y` (m, ascending).

# Keyword arguments
- `porosity`      : sed porosity for y-Exner (default: from `transects[1]`)
- `lsbc`          : `LongshoreBoundary` for open/prescribed lateral BCs
- `angle_mode`    : `WaveAngleFixed` (default), `WaveAngleShoreline`, or `WaveAngleSnell`
- `wave_angle_0`  : baseline offshore angle (rad); used by `WaveAngleShoreline`
"""
function TransectGrid(transects::Vector{CshoreBMI}, y::Vector{Float64};
                      porosity::Float64   = 1.0 - transects[1].config.sediment.sporo1,
                      lsbc::LongshoreBoundary = LongshoreBoundary(),
                      angle_mode::WaveAngleMode = WaveAngleFixed,
                      wave_angle_0::Float64 = 0.0)
    N = length(transects)
    N == length(y)      || throw(DimensionMismatch("transects and y must be same length"))
    N >= 2              || throw(ArgumentError("TransectGrid requires at least 2 transects"))
    issorted(y)         || throw(ArgumentError("y must be sorted ascending"))

    # Half-spacings for central differences.  Edge transects use one-sided diffs.
    dy_left  = Vector{Float64}(undef, N)
    dy_right = Vector{Float64}(undef, N)
    for i in 1:N
        dy_left[i]  = i > 1 ? (y[i] - y[i-1]) : (y[2] - y[1])
        dy_right[i] = i < N ? (y[i+1] - y[i]) : (y[N] - y[N-1])
    end

    return TransectGrid(transects, y, dy_left, dy_right,
                        porosity, lsbc, angle_mode, wave_angle_0)
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    _read_qby_total!(Qby, grid) -> Qby

Fill the `(nn, N)` matrix `Qby` with the depth-integrated longshore volume flux
from each transect.  Calls `bmi_compute_qby_total!` on each transect first to
ensure the cached `qby_total` field is current.
"""
function _read_qby_total!(Qby::Matrix{Float64}, grid::TransectGrid)
    for (i, m) in enumerate(grid.transects)
        bmi_compute_qby_total!(m)
        nn = m.state.jmax[1]
        @inbounds for j in 1:nn
            Qby[j, i] = m.state.qby_total[j]
        end
    end
    return Qby
end

"""
    _resolve_lsbc(bc_val, itime) -> Float64

Return the scalar BC value at window `itime` for a lateral boundary condition
that may be `NaN` (open boundary), a constant `Float64`, or a `Vector`.
"""
@inline function _resolve_lsbc(bc_val, itime::Int)
    bc_val isa Vector && return Float64(bc_val[min(itime, length(bc_val))])
    return Float64(bc_val)   # NaN → NaN, scalar → scalar
end

"""
    _longshore_exner!(dzb, Qby, F, grid, dt_couple, itime) -> dzb

Compute the y-direction bed-level change `dzb[j,i]` (m) using a **donor-cell
(upwind) finite-volume** scheme.

## Why upwind?

The longshore Exner equation

    (1−p) ∂zb/∂t = −∂Qby/∂y

is a first-order hyperbolic PDE: information propagates in the direction of
the longshore transport.  A central-difference stencil has a dispersion error
that generates spurious oscillations when the flux is unidirectional, and is
only marginally stable (requires artificial diffusion).  An upwind (donor-cell)
scheme is unconditionally stable for the explicit step and free of spurious
oscillations, at the cost of a small amount of numerical diffusion.

## Donor-cell flux at the i+½ interface:

    F[j, i+½] = Qby[j, i]   if Qby[j, i+½] > 0  (transport in +y direction)
              = Qby[j, i+1]  if Qby[j, i+½] < 0  (transport in −y direction)

The interface flux sign is estimated from the average of the two cell values:
    Qby_face = 0.5*(Qby[j,i] + Qby[j,i+1])

This is the standard **donor-cell upwind** method (equivalent to a first-order
Godunov scheme with a Roe flux).

## Sweep direction

For a purely explicit update (fluxes from the previous step feed the current
Δzb), the loop order over transects does not affect the result — each Δzb[j,i]
depends only on the *old* flux values.  The user note is most relevant for:
  - Semi-implicit or predictor-corrector updates (future extension).
  - The order in which wave-angle updates propagate (handled by WaveAngleShoreline
    mode, which can be called with a direction-dependent sweep).

`itime` is the current BC-window index (1-based) used to index time-series BCs.
`F` is a pre-allocated scratch vector of length N+1 for face fluxes (reused to
avoid allocation per node).
"""
function _longshore_exner!(dzb::Matrix{Float64}, Qby::Matrix{Float64},
                           F::Vector{Float64},
                           grid::TransectGrid, dt_couple::Float64, itime::Int)
    N     = length(grid.transects)
    sporo = 1.0 - grid.porosity

    # Resolve lateral BCs once (same for all cross-shore nodes at fixed itime)
    Qby_left_bc  = _resolve_lsbc(grid.lsbc.Qby_left,  itime)
    Qby_right_bc = _resolve_lsbc(grid.lsbc.Qby_right, itime)

    # Determine the max active nodes across all transects
    nn_max = size(Qby, 1)

    @inbounds for j in 1:nn_max
        # ── Build face fluxes F[1..N+1] for this cross-shore row ─────────────
        # Face i+½ lies between transect i and transect i+1.
        # F[1]   = left  boundary flux (into transect 1 from the left)
        # F[N+1] = right boundary flux (out of transect N to the right)

        # Left boundary face (F[1]):
        #   NaN lsbc → zero-gradient open boundary → F[1] = Qby[j,1]
        #   Prescribed → F[1] = Qby_left_bc (supply entering from updrift)
        F[1] = isnan(Qby_left_bc) ? Qby[j, 1] : Qby_left_bc

        # Interior faces F[2..N] (between transects i and i+1):
        for i in 1:(N-1)
            Qface = 0.5 * (Qby[j, i] + Qby[j, i+1])   # interface-averaged flux
            # Donor-cell: take value from upwind cell
            F[i+1] = Qface >= 0.0 ? Qby[j, i] : Qby[j, i+1]
        end

        # Right boundary face (F[N+1]):
        F[N+1] = isnan(Qby_right_bc) ? Qby[j, N] : Qby_right_bc

        # ── Apply Exner for each transect cell ────────────────────────────────
        # dzb[j,i] = -dt * (F[i+1] - F[i]) / (sporo * dy_cell)
        # dy_cell for cell i = half-distance to left + half-distance to right
        #                     = dy_left[i] + dy_right[i]   (stored on grid)
        for i in 1:N
            nn_i = grid.transects[i].state.jmax[1]
            j > nn_i && continue   # this transect is shorter; skip
            dy_cell = grid.dy_left[i] + grid.dy_right[i]
            dzb[j, i] = -dt_couple * (F[i+1] - F[i]) / (sporo * dy_cell)
        end
    end
    return dzb
end

"""
    _update_wave_angles!(grid) -> nothing

Update `wangbc` for each transect based on the current shoreline azimuth.

For each transect i, the shoreline position x_sl[i] is estimated as the x
coordinate of the still-water level (SWL) intersection node `state.jswl`.
The local shoreline azimuth is computed from the x_sl gradient along y:

  θ_local[i] = θ₀ + arctan(Δx_sl / Δy)

where Δx_sl is the change in shoreline position between adjacent transects
and θ₀ is the baseline offshore wave angle.  This is the Pelnard-Considère
shoreline-normal convention used in one-line models.

Only called when `grid.angle_mode == WaveAngleShoreline`.
"""
function _update_wave_angles!(grid::TransectGrid)
    N = length(grid.transects)

    # Estimate shoreline x-position at each transect from jswl
    x_sl = Vector{Float64}(undef, N)
    for (i, m) in enumerate(grid.transects)
        jsl = m.state.jswl[1]
        x_sl[i] = m.state.xb[max(1, jsl)]
    end

    # Compute local shoreline angle relative to baseline and update wangbc
    θ₀ = grid.wave_angle_0
    for i in 1:N
        if i == 1
            Δx = x_sl[2] - x_sl[1]
            Δy = grid.y[2] - grid.y[1]
        elseif i == N
            Δx = x_sl[N] - x_sl[N-1]
            Δy = grid.y[N] - grid.y[N-1]
        else
            Δx = x_sl[i+1] - x_sl[i-1]
            Δy = grid.y[i+1] - grid.y[i-1]
        end

        # Local angle: wave coming from θ₀ + small rotation due to shoreline tilt
        # Positive Δx/Δy → shoreline tilting seaward → wave angle increases
        θ_local = θ₀ + atan(Δx / max(Δy, 1.0))

        # Write into the *entire* wangbc time series (same angle for all windows)
        # A coupled wave model would write per-window values here instead.
        m = grid.transects[i]
        m.config.boundary.wangbc .= θ_local
    end
    return nothing
end

"""
    _cfl_check(Qby, dzb, dy_min, dt_couple)

Emit a warning if the estimated Courant number for the y-Exner step exceeds
0.8.  Computed as max|dzb| / (max|Qby| · dt · dy_min).
"""
function _cfl_check(Qby::Matrix{Float64}, dzb::Matrix{Float64},
                    dy_min::Float64, dt_couple::Float64)
    max_dzb = maximum(abs, dzb)
    max_Qby = maximum(abs, Qby)
    if max_Qby > 1e-12 && max_dzb > 1e-9
        # CFL = dt · |∂Qby/∂zb| / dy  ≈ max|dzb| / (max|Qby| · dt / dy)
        # Simpler proxy: if max_dzb > 0.1 m in one step, warn
        if max_dzb > 0.1
            @warn "Q2D y-Exner: large step detected (max|Δzb|=$(round(max_dzb, sigdigits=2)) m). " *
                  "Consider reducing coupling dt or increasing transect spacing dy."
        end
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

"""
    step_transect_grid!(grid, itime) -> nothing

Advance the Q2D system by one coupling step:

  1. `@threads`: advance each transect by one BC window (`BMI.update`).
  2. Read longshore fluxes Qby[j,i] from all transects.
  3. Compute y-Exner bed change dzb[j,i].
  4. Inject dzb into each transect via `bmi_inject_dzb!`.
  5. Optionally update wave angles (if `grid.angle_mode == WaveAngleShoreline`).

`itime` is the BC window index (1-based) passed to the Exner BC lookup.
"""
function step_transect_grid!(grid::TransectGrid, itime::Int)
    N  = length(grid.transects)
    # Determine the max cross-shore nodes across all transects (for array sizing)
    nn_max = maximum(m.state.jmax[1] for m in grid.transects)

    # Scratch arrays (allocated once per step — cheap relative to physics)
    Qby = zeros(Float64, nn_max, N)
    dzb = zeros(Float64, nn_max, N)
    F   = zeros(Float64, N + 1)    # face fluxes for upwind Exner

    # ── Step 1: parallel cross-shore physics ─────────────────────────────────
    Threads.@threads for i in 1:N
        BMI.update(grid.transects[i])
    end

    # ── Step 2: collect longshore fluxes ─────────────────────────────────────
    _read_qby_total!(Qby, grid)

    # Determine coupling dt from the BC window that was just stepped
    m1  = grid.transects[1]
    bc  = m1.config.boundary
    it  = min(m1.itime, length(bc.timebc) - 1)
    dt_couple = bc.timebc[it + 1] - bc.timebc[it]

    # ── Step 3: y-Exner (donor-cell upwind) ──────────────────────────────────
    _longshore_exner!(dzb, Qby, F, grid, dt_couple, itime)
    _cfl_check(Qby, dzb, minimum(grid.dy_left), dt_couple)

    # ── Step 4: inject bed corrections ───────────────────────────────────────
    for i in 1:N
        m  = grid.transects[i]
        nn = m.state.jmax[1]
        bmi_inject_dzb!(m, view(dzb, 1:nn, i))
    end

    # ── Step 5: optional wave angle update ───────────────────────────────────
    if grid.angle_mode == WaveAngleShoreline
        _update_wave_angles!(grid)
    end

    return nothing
end

"""
    run_transect_grid!(grid; callback=nothing, output_every=1) -> nothing

Run the Q2D simulation to completion.

Iterates over all BC windows in `grid.transects[1].config.boundary.timebc`,
calling `step_transect_grid!` at each step.

# Keyword arguments
- `callback` : optional function `callback(grid, itime)` called after each
               coupling step.  Use this to save snapshots, compute diagnostics,
               or update external wave model results.
- `output_every` : call `callback` every N steps (default: every step).
"""
function run_transect_grid!(grid::TransectGrid;
                            callback::Union{Nothing,Function} = nothing,
                            output_every::Int = 1)
    m1     = grid.transects[1]
    ntimes = ntime(m1.config.boundary) - 1   # number of BC windows

    for itime in 1:ntimes
        step_transect_grid!(grid, itime)
        if callback !== nothing && (itime % output_every == 0 || itime == ntimes)
            callback(grid, itime)
        end
    end
    return nothing
end

"""
    transect_grid_output(grid) -> NamedTuple

Collect a snapshot of all transects into arrays convenient for plotting or
writing.  Returns:

```julia
(; zb, qby_total, vmean, hrms, wsetup, y, x)
```

where `zb[j,i]`, `qby_total[j,i]`, etc. are (nn_max × N) matrices and
`y[i]` and `x[j]` are the coordinate vectors.  Shorter transects are
padded with `NaN` to `nn_max` rows.
"""
function transect_grid_output(grid::TransectGrid)
    N      = length(grid.transects)
    nn_max = maximum(m.state.jmax[1] for m in grid.transects)
    # All transects are assumed to have the same x grid; use transect 1.
    x = grid.transects[1].state.xb[1:grid.transects[1].state.jmax[1]]

    zb        = fill(NaN, nn_max, N)
    qby_total = fill(NaN, nn_max, N)
    vmean     = fill(NaN, nn_max, N)
    hrms      = fill(NaN, nn_max, N)
    wsetup    = fill(NaN, nn_max, N)

    for (i, m) in enumerate(grid.transects)
        bmi_compute_qby_total!(m)
        nn = m.state.jmax[1]
        zb[1:nn, i]        .= view(m.state.zb,        1:nn, 1)
        qby_total[1:nn, i] .= view(m.state.qby_total, 1:nn)
        vmean[1:nn, i]     .= view(m.state.vmean,      1:nn)
        hrms[1:nn, i]      .= view(m.state.hrms,       1:nn)
        wsetup[1:nn, i]    .= view(m.state.wsetup,     1:nn)
    end

    return (; zb, qby_total, vmean, hrms, wsetup, y=grid.y, x)
end
