# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
bmi.jl — Basic Model Interface (BMI v2.0) bindings for CSHORE.jl.

Implements the CSDMS BMI specification via `BasicModelInterface.jl` so that
CSHORE.jl can be driven by external couplers

Usage:

    import BasicModelInterface as BMI
    using CSHORE: CshoreBMI, build_config

    cfg = build_config(...)
    m = CshoreBMI(cfg)                            # programmatic config
    # or:
    # m = BMI.initialize(CshoreBMI, "config.toml")  # file-based (TODO)

    while BMI.get_current_time(m) < BMI.get_end_time(m)
        BMI.update(m)
        # Couple: read/write state via get_value_ptr / set_value
        zb = BMI.get_value_ptr(m, "zb")
        # ...
    end
    BMI.finalize(m)

Grid: a single 1D uniform-rectilinear grid with id `0`. All variables share
this grid. Per-fraction variables (`qbx`, `qsx`, `bed_mass`) have additional
rank beyond the base grid — BMI handles this via `get_var_nbytes` returning
the total byte count of the underlying array.

BMI step granularity: one "update" advances the simulation by one boundary-
condition window (`config.boundary.timebc[itime]` → `timebc[itime+1]`). Sub-
stepping within a BC window is handled internally by the adaptive CFL-based
Exner solver.
==============================================================================#

import BasicModelInterface as BMI

"""
    CshoreBMI

BMI wrapper around a `CshoreConfig` + `CshoreState` pair. Tracks the current
boundary-condition window index (`itime`) and owns an optional `NetcdfWriter`
for opt-in output.

Fields:
- `config` — immutable CshoreConfig (bathymetry, BCs, sediment, options)
- `state`  — mutable CshoreState (current solution)
- `itime`  — BC window index. `0` before the first `BMI.update`; after
             `update` it equals the window that was just advanced.
- `writer` — optional NetcdfWriter for per-step output
"""
mutable struct CshoreBMI
    config::CshoreConfig
    state::CshoreState
    itime::Int
    writer::Union{Nothing,NetcdfWriter}
end

"""
    CshoreBMI(config::CshoreConfig; outfile=nothing, output_interval_s=0.0)

Convenience programmatic constructor. Allocates state, applies initial
bathymetry, computes derived constants, and optionally opens a NetCDF
writer. Equivalent to `BMI.initialize(CshoreBMI, toml_file)` but takes
a `CshoreConfig` directly.
"""
function CshoreBMI(config::CshoreConfig;
    outfile::Union{Nothing,AbstractString}=nothing,
    output_interval_s::Real=0.0,
    outdir::AbstractString=".")
    state = initialize_state(config)
    apply_initial_bathymetry!(state, config)
    compute_derived_constants!(state, config)
    for l in 1:config.options.iline
        compute_bed_slope!(state, config, l)
    end
    writer = outfile === nothing ? nothing :
             open_netcdf(joinpath(outdir, outfile), config, state;
        output_interval_s=output_interval_s)
    m = CshoreBMI(config, state, 0, writer)
    # Record the initial state (t=0) if a writer is attached
    writer !== nothing && write_step!(writer, state, state.time)
    return m
end

# --------------------------------------------------------------------------
# Variable registry
# --------------------------------------------------------------------------
#
# Each entry:  name => (getter, units, standard_name, grid_id)
#
# The getter receives the `CshoreBMI` and returns a view into the underlying
# state array. `grid_id = 0` for all 1D fields. Per-fraction fields still
# use grid 0 for their first axis but carry extra data per node — BMI
# consumers use `get_var_nbytes` to size their buffers correctly.
#
# The BMI layer exposes state arrays directly. Couplers that mutate
# `bed_mass` directly bypass the grain-sorting machinery.

const _BMI_VAR_GETTERS = Dict{String,Function}()
const _BMI_VAR_UNITS = Dict{String,String}()
const _BMI_VAR_CFNAMES = Dict{String,String}()

function _register_bmi_var!(name::AbstractString, getter::Function,
    units::AbstractString, cfname::AbstractString)
    _BMI_VAR_GETTERS[name] = getter
    _BMI_VAR_UNITS[name] = units
    _BMI_VAR_CFNAMES[name] = cfname
    return nothing
end

# Register all output variables at module load
function _register_bmi_vars!()
    empty!(_BMI_VAR_GETTERS)
    empty!(_BMI_VAR_UNITS)
    empty!(_BMI_VAR_CFNAMES)

    # ── Cross-shore hydrodynamics ─────────────────────────────────────────────
    _register_bmi_var!("zb", m -> view(m.state.zb, 1:m.state.jmax[1], 1),
        "m", "sea_floor_depth_below_sea_surface")
    _register_bmi_var!("hrms", m -> view(m.state.hrms, 1:m.state.jmax[1]),
        "m", "sea_surface_wave_rms_height")
    _register_bmi_var!("h", m -> view(m.state.h, 1:m.state.jmax[1]),
        "m", "sea_water_depth")
    _register_bmi_var!("wsetup", m -> view(m.state.wsetup, 1:m.state.jmax[1]),
        "m", "sea_surface_wave_setup")
    _register_bmi_var!("sigma", m -> view(m.state.sigma, 1:m.state.jmax[1]),
        "m", "")
    _register_bmi_var!("umean", m -> view(m.state.umean, 1:m.state.jmax[1]),
        "m s-1", "sea_water_x_velocity")
    _register_bmi_var!("ustd", m -> view(m.state.ustd, 1:m.state.jmax[1]),
        "m s-1", "")
    _register_bmi_var!("qbreak", m -> view(m.state.qbreak, 1:m.state.jmax[1]),
        "1", "")

    # ── Longshore hydrodynamics (Q2D coupling) ────────────────────────────────
    # vmean and vstd are the time-averaged longshore current and its std dev.
    # These are the primary coupling quantities for Q2D momentum exchange.
    _register_bmi_var!("vmean", m -> view(m.state.vmean, 1:m.state.jmax[1]),
        "m s-1", "sea_water_y_velocity")
    _register_bmi_var!("vstd", m -> view(m.state.vstd, 1:m.state.jmax[1]),
        "m s-1", "")

    # ── Sediment transport ────────────────────────────────────────────────────
    _register_bmi_var!("qbx", m -> view(m.state.qbx, 1:m.state.jmax[1], :),
        "m2 s-1", "")
    _register_bmi_var!("qsx", m -> view(m.state.qsx, 1:m.state.jmax[1], :),
        "m2 s-1", "")
    _register_bmi_var!("q_total", m -> view(m.state.q_total, 1:m.state.jmax[1]),
        "m2 s-1", "")

    # Longshore transport per fraction — core Q2D coupling flux.
    # The orchestrator reads these after each BMI.update and computes ∂Qby/∂y.
    _register_bmi_var!("qby", m -> view(m.state.qby, 1:m.state.jmax[1], :),
        "m2 s-1", "")
    _register_bmi_var!("qsy", m -> view(m.state.qsy, 1:m.state.jmax[1], :),
        "m2 s-1", "")
    # Scalar sum over fractions: Qby_total = Σ_k (qby[:,k] + qsy[:,k]) / sporo
    # Pre-computed by the Q2D orchestrator or by set_value("qby_total", ...).
    # Exposed here as a read-back so the orchestrator can inspect it.
    _register_bmi_var!("qby_total", m -> view(m.state.qby_total, 1:m.state.jmax[1]),
        "m2 s-1", "")

    # ── Bed composition ───────────────────────────────────────────────────────
    _register_bmi_var!("bed_mass", m -> view(m.state.bed_mass, 1:m.state.jmax[1], 1, :),
        "kg m-2", "")

    # ── Cohesive (mud) sediment (CohesiveSedimentConfig) ─────────────────────
    # Zero-filled vectors when config.cohesive === nothing, so couplers can
    # always read them unconditionally.
    _register_bmi_var!("cohesive_bed_mass",
        m -> view(m.state.cohesive_bed_mass, 1:m.state.jmax[1]),
        "kg m-2", "")
    _register_bmi_var!("cohesive_concentration",
        m -> view(m.state.cohesive_concentration, 1:m.state.jmax[1]),
        "kg m-3", "")

    # ── Settable boundary forcing (Q2D angle + SWL injection) ─────────────────
    # These expose the *full* time series stored in BoundaryTimeSeries as
    # mutable Vector views. The Q2D orchestrator can overwrite individual
    # window values before calling BMI.update to inject neighbour-derived
    # wave angles or alongshore SWL gradients.
    # NOTE: BoundaryTimeSeries fields are Vector{Float64} (heap-allocated),
    # so mutating the view is safe even though CshoreConfig is declared immutable.
    _register_bmi_var!("wangbc", m -> m.config.boundary.wangbc,
        "rad", "sea_surface_wave_to_direction")
    _register_bmi_var!("swlbc", m -> m.config.boundary.swlbc,
        "m", "sea_surface_height_above_mean_sea_level")
    _register_bmi_var!("hrmsbc", m -> m.config.boundary.hrmsbc,
        "m", "sea_surface_wave_rms_height")

    return nothing
end

# --------------------------------------------------------------------------
# Q2D helpers: dzb injection and qby_total accumulation
# --------------------------------------------------------------------------

"""
    bmi_inject_dzb!(m, dzb)

Add a bed-elevation correction `dzb[j]` (m) to each node of transect `m`,
updating both `state.zb` and rescaling `state.bed_mass` to remain consistent.

This is the mechanism by which the Q2D orchestrator applies the y-direction
Exner equation result after each coupling step. Positive `dzb` means
deposition; negative means erosion.

The bed-mass rescaling preserves the *fractional composition* at each node —
the total mass changes proportionally but the grain-size distribution is
unchanged. This is the same update strategy used internally by `exner_step!`.
"""
function bmi_inject_dzb!(m::CshoreBMI, dzb::AbstractVector{Float64})
    l = 1
    nn = m.state.jmax[l]
    sed = m.config.sediment
    sporo = sed.sporo1           # sporo1 = 1 - porosity (CSHORE convention)
    nf = nfractions(m.config.multifraction)
    nlyr = m.config.multifraction.nlayers

    length(dzb) >= nn || error("dzb length $(length(dzb)) < jmax=$nn")

    @inbounds for j in 1:nn
        Δz = dzb[j]
        iszero(Δz) && continue

        # Update bed surface
        m.state.zb[j, l] += Δz

        # Rescale bed_mass proportionally in all layers to match new zb.
        # Strategy: the total column mass changes by Δz * ρs * sporo.
        # We spread the change uniformly across layers weighted by existing mass.
        total_mass = 0.0
        @inbounds for k in 1:nf, lyr in 1:nlyr
            total_mass += m.state.bed_mass[j, lyr, k]
        end

        if total_mass > 1e-10
            delta_mass = Δz * sed.sg * 1025.0 * sporo   # ρ_s ≈ sg × ρ_w
            scale = (total_mass + delta_mass) / total_mass
            @inbounds for k in 1:nf, lyr in 1:nlyr
                m.state.bed_mass[j, lyr, k] = max(0.0, m.state.bed_mass[j, lyr, k] * scale)
            end
        end
    end
    return nothing
end

"""
    bmi_compute_qby_total!(m)

Accumulate `state.qby_total[j]` as the depth-integrated longshore volume
flux at each node: `Qby[j] = Σ_k (qby[j,k] + qsy[j,k]) / sporo`.

Called by the Q2D orchestrator before reading "qby_total".  Also called
internally at the end of each `step_bc_window!` to keep the field current.
"""
function bmi_compute_qby_total!(m::CshoreBMI)
    sporo = m.config.sediment.sporo1   # sporo1 = 1 - porosity
    nn = m.state.jmax[1]
    nf = nfractions(m.config.multifraction)
    @inbounds for j in 1:nn
        s = 0.0
        for k in 1:nf
            s += m.state.qby[j, k] + m.state.qsy[j, k]
        end
        m.state.qby_total[j] = s / sporo
    end
    return nothing
end

# Populated on module load (called from CSHORE.jl at the end of include chain)
_register_bmi_vars!()

# --------------------------------------------------------------------------
# BMI v2.0 method implementations
# --------------------------------------------------------------------------

# ---- Component info ----
BMI.get_component_name(::CshoreBMI) = "CSHORE.jl"

# ---- Lifecycle ----

"""
    BMI.initialize(::Type{CshoreBMI}, config_file) -> CshoreBMI

Stub — programmatic `CshoreBMI(cfg)` is the recommended entry point.
File-based initialization (TOML / `.infile` parser) is not yet implemented.
"""
function BMI.initialize(::Type{CshoreBMI}, config_file::AbstractString)
    error("BMI.initialize(CshoreBMI, \"$config_file\"): file-based initialization " *
          "not yet implemented. Use `CshoreBMI(cfg)` with a programmatic config.")
end

function BMI.update(m::CshoreBMI)
    ntimes = ntime(m.config.boundary)
    next_itime = m.itime + 1
    if next_itime > ntimes - 1
        error("BMI.update: cannot advance past the last boundary window " *
              "(itime=$(m.itime), ntimes=$ntimes). Use get_end_time to check.")
    end
    step_bc_window!(m.state, m.config, next_itime)
    m.itime = next_itime
    m.writer !== nothing && maybe_write_step!(m.writer, m.state, m.state.time)
    return nothing
end

function BMI.update_until(m::CshoreBMI, t::Real)
    while m.state.time < t - 1e-9
        ntimes = ntime(m.config.boundary)
        if m.itime ≥ ntimes - 1
            break
        end
        BMI.update(m)
    end
    return nothing
end

function BMI.finalize(m::CshoreBMI)
    if m.writer !== nothing
        close_netcdf!(m.writer)
        m.writer = nothing
    end
    return nothing
end

# ---- Variable info ----

BMI.get_input_var_names(::CshoreBMI) = collect(keys(_BMI_VAR_GETTERS))
BMI.get_output_var_names(::CshoreBMI) = collect(keys(_BMI_VAR_GETTERS))
BMI.get_input_item_count(::CshoreBMI) = length(_BMI_VAR_GETTERS)
BMI.get_output_item_count(::CshoreBMI) = length(_BMI_VAR_GETTERS)

function BMI.get_var_grid(m::CshoreBMI, name::AbstractString)
    haskey(_BMI_VAR_GETTERS, name) || _unknown_var(name)
    return 0
end

function BMI.get_var_type(m::CshoreBMI, name::AbstractString)
    haskey(_BMI_VAR_GETTERS, name) || _unknown_var(name)
    return "float64"
end

function BMI.get_var_units(m::CshoreBMI, name::AbstractString)
    haskey(_BMI_VAR_UNITS, name) || _unknown_var(name)
    return _BMI_VAR_UNITS[name]
end

BMI.get_var_itemsize(m::CshoreBMI, name::AbstractString) = 8  # Float64

function BMI.get_var_nbytes(m::CshoreBMI, name::AbstractString)
    haskey(_BMI_VAR_GETTERS, name) || _unknown_var(name)
    arr = _BMI_VAR_GETTERS[name](m)
    return length(arr) * 8
end

function BMI.get_var_location(m::CshoreBMI, name::AbstractString)
    haskey(_BMI_VAR_GETTERS, name) || _unknown_var(name)
    return "node"
end

# ---- Time ----

BMI.get_current_time(m::CshoreBMI) = m.state.time
BMI.get_start_time(m::CshoreBMI) = first(m.config.boundary.timebc)
BMI.get_end_time(m::CshoreBMI) = last(m.config.boundary.timebc)
BMI.get_time_units(m::CshoreBMI) = "s"

"""
    BMI.get_time_step(m)

Returns the **last adaptive sub-step** taken by `exner_step!` (stored in
`state.delt`). The CSHORE.jl sub-step is chosen internally via a CFL
criterion, so this value reflects the most recent step rather than a
fixed interval. For a fixed coupling interval, BMI consumers should rely
on `get_start_time` / `get_end_time` and the BC window spacing in
`config.boundary.timebc`.
"""
BMI.get_time_step(m::CshoreBMI) = m.state.delt

# ---- Value access ----

function BMI.get_value(m::CshoreBMI, name::AbstractString, dest::AbstractArray)
    haskey(_BMI_VAR_GETTERS, name) || _unknown_var(name)
    arr = _BMI_VAR_GETTERS[name](m)
    copyto!(dest, vec(arr))
    return dest
end

"""
    BMI.get_value_ptr(m, name) -> AbstractArray

Returns a `view` into the underlying state array — **zero copy**. Mutating
the returned view mutates `state` directly; this is the BMI convention
for tight-coupling hot paths.
"""
function BMI.get_value_ptr(m::CshoreBMI, name::AbstractString)
    haskey(_BMI_VAR_GETTERS, name) || _unknown_var(name)
    return _BMI_VAR_GETTERS[name](m)
end

function BMI.get_value_at_indices(m::CshoreBMI, name::AbstractString,
    dest::AbstractArray, inds::AbstractArray)
    haskey(_BMI_VAR_GETTERS, name) || _unknown_var(name)
    arr = vec(_BMI_VAR_GETTERS[name](m))
    for (k, i) in enumerate(inds)
        dest[k] = arr[i]
    end
    return dest
end

function BMI.set_value(m::CshoreBMI, name::AbstractString, src::AbstractArray)
    haskey(_BMI_VAR_GETTERS, name) || _unknown_var(name)
    arr = _BMI_VAR_GETTERS[name](m)
    vec(arr) .= vec(src)
    return nothing
end

function BMI.set_value_at_indices(m::CshoreBMI, name::AbstractString,
    inds::AbstractArray, src::AbstractArray)
    haskey(_BMI_VAR_GETTERS, name) || _unknown_var(name)
    arr = vec(_BMI_VAR_GETTERS[name](m))
    for (k, i) in enumerate(inds)
        arr[i] = src[k]
    end
    return nothing
end

# ---- Grid metadata (1D uniform_rectilinear) ----

function _check_grid(grid::Int)
    grid == 0 || error("CshoreBMI only supports grid id 0 (got $grid)")
end

BMI.get_grid_rank(m::CshoreBMI, grid::Int) = (_check_grid(grid); 1)
BMI.get_grid_size(m::CshoreBMI, grid::Int) = (_check_grid(grid); m.state.jmax[1])
BMI.get_grid_type(m::CshoreBMI, grid::Int) = (_check_grid(grid); "uniform_rectilinear")

function BMI.get_grid_shape(m::CshoreBMI, grid::Int, shape::AbstractArray)
    _check_grid(grid)
    shape[1] = m.state.jmax[1]
    return shape
end

function BMI.get_grid_spacing(m::CshoreBMI, grid::Int, spacing::AbstractArray)
    _check_grid(grid)
    spacing[1] = m.config.grid.dx
    return spacing
end

function BMI.get_grid_origin(m::CshoreBMI, grid::Int, origin::AbstractArray)
    _check_grid(grid)
    origin[1] = m.state.xb[1]
    return origin
end

function BMI.get_grid_x(m::CshoreBMI, grid::Int, x::AbstractArray)
    _check_grid(grid)
    nn = m.state.jmax[1]
    copyto!(x, view(m.state.xb, 1:nn))
    return x
end

# 1D model → y and z are trivially zero-length
function BMI.get_grid_y(m::CshoreBMI, grid::Int, y::AbstractArray)
    _check_grid(grid)
    return y
end
function BMI.get_grid_z(m::CshoreBMI, grid::Int, z::AbstractArray)
    _check_grid(grid)
    return z
end

# Unstructured topology methods — not applicable to a 1D structured model
const _UNSTRUCT_ERR = "CshoreBMI uses a uniform_rectilinear grid; " *
                      "unstructured topology methods are not implemented."
BMI.get_grid_node_count(m::CshoreBMI, grid::Int) = error(_UNSTRUCT_ERR)
BMI.get_grid_edge_count(m::CshoreBMI, grid::Int) = error(_UNSTRUCT_ERR)
BMI.get_grid_face_count(m::CshoreBMI, grid::Int) = error(_UNSTRUCT_ERR)
BMI.get_grid_edge_nodes(m::CshoreBMI, grid, edges) = error(_UNSTRUCT_ERR)
BMI.get_grid_face_edges(m::CshoreBMI, grid, edges) = error(_UNSTRUCT_ERR)
BMI.get_grid_face_nodes(m::CshoreBMI, grid, nodes) = error(_UNSTRUCT_ERR)
BMI.get_grid_nodes_per_face(m::CshoreBMI, grid, nf) = error(_UNSTRUCT_ERR)

# ---- helpers ----

function _unknown_var(name::AbstractString)
    known = join(sort!(collect(keys(_BMI_VAR_GETTERS))), ", ")
    error("Unknown BMI variable \"$name\". Known variables: $known")
end
