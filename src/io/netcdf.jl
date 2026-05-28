# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
io/netcdf.jl — CF-compliant NetCDF output writer for CSHORE.jl.

Writes a single self-describing NetCDF-4 file per simulation with:
  - Unlimited `time` dimension (appended one slice per BC window or at a
    user-configurable interval)
  - `x` (cross-shore) coordinate dimension
  - `fraction` and `layer` dimensions for multifraction runs
  - 2D time-series fields: zb(time,x), hrms, h, wsetup, umean, ustd, q_total,
                           hrms_ig (infragravity RMS height; zeros when IgConfig absent),
                           cohesive_bed_mass / cohesive_concentration (mud;
                           zeros when CohesiveSedimentConfig absent)
  - 1D scalar time-series: runup_x, runup_z (= TWL), overtopping_flux
  - 3D per-fraction fields: qbx(time,x,fraction), qsx, bed_mass
  - CF-1.10 Conventions and provenance attrs

Usage via the driver:

    run_simulation!(config; outfile="out.nc", output_interval_s=3600.0)

Or manually:

    writer = open_netcdf("out.nc", config, state)
    write_step!(writer, state, 0.0)
    # ... time loop ...
    write_step!(writer, state, t)
    close_netcdf!(writer)
==============================================================================#

"""
    NetcdfWriter

Holds the open NCDataset handle plus bookkeeping for incremental writes
along the unlimited `time` dimension.

Fields:
- `ds`                  — NCDataset (open file handle)
- `it`                  — current time-slice index (0 before first write)
- `last_write_time`     — absolute simulation time of the last slice written
- `output_interval_s`   — minimum dt between slices; `0.0` writes every call
- `has_multifraction`   — whether the file has `bed_mass` / `d50_surface` vars
- `write_composition`   — write `bed_mass`, `d50_surface`, `d50_bulk` (default true)
- `write_T_profile`     — write full `T_profile(x, depth, time)` (default true;
                           set false to save ~5 GB/run; ALT and zb_hard still written)
- `write_transport`     — write per-fraction `qbx` / `qsx` arrays (default true)
"""
mutable struct NetcdfWriter
    ds::NCDataset
    it::Int
    last_write_time::Float64
    output_interval_s::Float64
    has_multifraction::Bool
    has_thermal::Bool
    write_composition::Bool
    write_T_profile::Bool
    write_transport::Bool
    has_provenance::Bool
end

function _attrs(units::String, longname::String, standard::String="")
    a = Pair{String,Any}[
        "units"     => units,
        "long_name" => longname,
    ]
    isempty(standard) || push!(a, "standard_name" => standard)
    return a
end

"""
    open_netcdf(path, config, state; output_interval_s=0.0) -> NetcdfWriter

Create a new NetCDF file at `path` with the CSHORE.jl variable layout and
return a writer handle. The file stays open until `close_netcdf!` is called.

`output_interval_s == 0.0` (default) writes every call to `maybe_write_step!`;
any positive value enforces a minimum dt between slices.
"""
function open_netcdf(path::AbstractString, config::CshoreConfig, state::CshoreState;
                     output_interval_s::Real=0.0,
                     write_composition::Bool=true,
                     write_T_profile::Bool=true,
                     write_transport::Bool=true,
                     provenance::Union{ProvenanceConfig,Nothing}=nothing,
                     prov_state::Union{ProvenanceState,Nothing}=nothing)
    ds = NCDataset(path, "c"; format=:netcdf4)
    nn = state.jmax[1]
    nf = nfractions(config.multifraction)
    nlayers = config.multifraction.nlayers

    # --- dimensions ---
    defDim(ds, "x", nn)
    defDim(ds, "fraction", nf)
    defDim(ds, "layer", nlayers)
    defDim(ds, "time", Inf)

    # --- global attributes ---
    ds.attrib["title"]       = "CSHORE.jl simulation output"
    ds.attrib["Conventions"] = "CF-1.10"
    ds.attrib["source"]      = "CSHORE.jl v0.1 (Julia port of CSHORE USACE)"
    ds.attrib["history"]     = string("Created ",
        Dates.format(now(UTC), "yyyy-mm-ddTHH:MM:SS"), " UTC by CSHORE.jl")
    ds.attrib["institution"] = "CSHORE.jl project"
    ds.attrib["model_dx"]            = config.grid.dx
    ds.attrib["model_gamma"]         = config.gamma
    ds.attrib["model_iprofl"]        = config.options.iprofl
    ds.attrib["model_iangle"]        = config.options.iangle
    ds.attrib["model_nfractions"]    = nf
    ds.attrib["model_nlayers"]       = nlayers
    ds.attrib["sediment_d50_m"]      = config.sediment.d50
    ds.attrib["sediment_porosity"]   = 1.0 - config.sediment.sporo1
    if nf > 1
        ds.attrib["grain_sizes_m"] = collect(config.multifraction.grain_sizes)
    end
    ds.attrib["max_dzb_per_step"] = config.max_dzb_per_step

    # --- coordinates ---
    xv = defVar(ds, "x", Float64, ("x",); attrib = [
        "units"         => "m",
        "long_name"     => "cross-shore distance",
        "standard_name" => "projection_x_coordinate",
        "axis"          => "X",
    ])
    xv[:] = view(state.xb, 1:nn)

    defVar(ds, "time", Float64, ("time",); attrib = [
        "units"         => "seconds since 2000-01-01 00:00:00",
        "long_name"     => "simulation time (model clock, seconds from t=0)",
        "standard_name" => "time",
        "axis"          => "T",
        "calendar"      => "standard",
    ])

    if nf > 1
        fv = defVar(ds, "fraction", Float64, ("fraction",); attrib = [
            "units"     => "m",
            "long_name" => "grain diameter per fraction",
        ])
        fv[:] = collect(config.multifraction.grain_sizes)
    end

    layer_v = defVar(ds, "layer", Int32, ("layer",); attrib = [
        "long_name" => "bed layer index (1 = surface)",
    ])
    layer_v[:] = collect(Int32, 1:nlayers)

    # --- time-series variables ---
    # NOTE: NCDatasets is column-major like Julia. Listing dims as
    # ("x","time") produces on-disk CF-conformant (time, x) layout.
    defVar(ds, "zb", Float64, ("x", "time"); attrib =
        _attrs("m", "bed elevation", "sea_floor_depth_below_sea_surface"))
    defVar(ds, "hrms", Float64, ("x", "time"); attrib =
        _attrs("m", "RMS wave height", "sea_surface_wave_rms_height"))
    defVar(ds, "h", Float64, ("x", "time"); attrib =
        _attrs("m", "total water depth", "sea_water_depth"))
    defVar(ds, "wsetup", Float64, ("x", "time"); attrib =
        _attrs("m", "wave setup"))
    defVar(ds, "sigma", Float64, ("x", "time"); attrib =
        _attrs("m", "wave standard deviation (Hrms / √8)"))
    defVar(ds, "umean", Float64, ("x", "time"); attrib =
        _attrs("m s-1", "mean cross-shore velocity (undertow)", "sea_water_x_velocity"))
    defVar(ds, "ustd", Float64, ("x", "time"); attrib =
        _attrs("m s-1", "cross-shore velocity standard deviation (orbital)"))
    defVar(ds, "vmean", Float64, ("x", "time"); attrib =
        _attrs("m s-1", "mean alongshore velocity (longshore current)", "sea_water_y_velocity"))
    defVar(ds, "qbreak", Float64, ("x", "time"); attrib =
        _attrs("1", "fraction of breaking waves"))
    defVar(ds, "q_total", Float64, ("x", "time"); attrib =
        _attrs("m2 s-1", "total cross-shore sediment transport (summed over fractions)"))
    defVar(ds, "hrms_ig", Float64, ("x", "time"); attrib = [
        "units"      => "m",
        "long_name"  => "infragravity RMS wave height (IgConfig Layer 1/2); " *
                        "zero everywhere when config.ig is nothing",
        "_FillValue" => NaN,
    ])
    defVar(ds, "gw_eta", Float64, ("x", "time"); attrib = [
        "units"      => "m",
        "long_name"  => "beach groundwater water-table elevation above datum " *
                        "(GroundwaterConfig Boussinesq model); " *
                        "zero everywhere when config.groundwater is nothing",
        "_FillValue" => NaN,
    ])
    defVar(ds, "theta", Float64, ("x", "time"); attrib = [
        "units"      => "1",
        "long_name"  => "volumetric surface moisture content (Van Genuchten retention curve); " *
                        "feeds aeolian transport threshold when GroundwaterConfig is active; " *
                        "zero everywhere when config.groundwater is nothing",
        "_FillValue" => NaN,
    ])
    defVar(ds, "cohesive_bed_mass", Float64, ("x", "time"); attrib = [
        "units"      => "kg m-2",
        "long_name"  => "cohesive (mud) bed mass per unit area " *
                        "(CohesiveSedimentConfig Partheniades-Krone); " *
                        "zero everywhere when config.cohesive is nothing",
        "_FillValue" => NaN,
    ])
    defVar(ds, "cohesive_concentration", Float64, ("x", "time"); attrib = [
        "units"      => "kg m-3",
        "long_name"  => "depth-averaged suspended cohesive (mud) concentration; " *
                        "zero everywhere when config.cohesive is nothing",
        "_FillValue" => NaN,
    ])

    # --- Runup / total water level (scalar per time step) ---
    # runup_x / runup_z locate the landward limit of the swash zone.
    # With iover=1 this is state.jdry (the node where swash depth → 0),
    # which is the true process-based runup position.  With iover=0 no
    # swash march is performed and the still-water shoreline (jswl) is
    # used as a proxy — the resulting runup_z then equals the bed
    # elevation at the SWL crossing, a lower bound on the true runup.
    # runup_z equals the total water level (TWL) at the swash limit
    # because water depth is zero there by construction.
    # overtopping_flux is the converged mean discharge per unit width
    # at the crest node (m²/s); zero when iover=0 or no overtopping.
    defVar(ds, "runup_x", Float64, ("time",); attrib = [
        "units"     => "m",
        "long_name" => "cross-shore position of wave runup / swash limit " *
                       "(jdry node when iover=1, jswl proxy when iover=0)",
        "_FillValue" => NaN,
    ])
    defVar(ds, "runup_z", Float64, ("time",); attrib = [
        "units"     => "m",
        "long_name" => "elevation of wave runup / swash limit = total water level (TWL) " *
                       "(bed elevation at runup position; NaN before first wave step)",
        "_FillValue" => NaN,
    ])
    defVar(ds, "overtopping_flux", Float64, ("time",); attrib = [
        "units"     => "m2 s-1",
        "long_name" => "mean wave overtopping discharge per unit width at crest (qotf); " *
                       "0.0 when iover=0 or when swash does not reach the crest",
    ])

    if write_transport
        defVar(ds, "qbx", Float64, ("x", "fraction", "time"); attrib =
            _attrs("m2 s-1", "per-fraction cross-shore bedload transport rate"))
        defVar(ds, "qsx", Float64, ("x", "fraction", "time"); attrib =
            _attrs("m2 s-1", "per-fraction cross-shore suspended transport rate"))
    end

    has_mf = nf > 1
    if has_mf && write_composition
        # bed_mass: stores the FULL (x, layer, fraction, time) — the
        # active layer (layer=1) is the top, deeper layers are buried
        # stratigraphy. Size: nn × nlayers × nf × nt.
        defVar(ds, "bed_mass", Float64, ("x", "layer", "fraction", "time");
               attrib = _attrs("kg m-2",
                               "bed mass per layer per fraction " *
                               "(layer 1 = active surface, higher indices = deeper)"))
        defVar(ds, "d50_surface", Float64, ("x", "time"); attrib =
            _attrs("m", "surface-layer (active) median grain diameter"))
        defVar(ds, "d50_bulk", Float64, ("x", "time"); attrib =
            _attrs("m", "vertically-integrated median grain diameter (all layers)"))
    end

    # --- thermal variables (ALT, zb_hard, surface temperature, full profile) ---
    has_therm = config.thermal !== nothing
    if has_therm
        nz_therm = config.thermal.nz
        dz_therm = config.thermal.dz

        defVar(ds, "ALT", Float64, ("x", "time"); attrib =
            _attrs("m", "active layer thickness (depth to permafrost)"))
        defVar(ds, "zb_hard", Float64, ("x", "time"); attrib =
            _attrs("m", "hardbottom / permafrost table elevation"))
        defVar(ds, "T_surface", Float64, ("x", "time"); attrib =
            _attrs("degC", "surface temperature (top thermal cell)"))

        if write_T_profile
            # Full thermal profile: T(x, thermal_depth, time) — large output,
            # disable with write_T_profile=false for faster calibration runs.
            defDim(ds, "thermal_depth", nz_therm)
            tdv = defVar(ds, "thermal_depth", Float64, ("thermal_depth",); attrib = [
                "units"     => "m",
                "long_name" => "depth below local bed surface (positive downward)",
                "positive"  => "down",
            ])
            tdv[:] = collect(Float64, range(dz_therm / 2, step=dz_therm, length=nz_therm))

            defVar(ds, "T_profile", Float64, ("x", "thermal_depth", "time"); attrib =
                _attrs("degC", "subsurface temperature profile (cell-centered, " *
                       "depth below local bed surface)"))
        end

        ds.attrib["thermal_nz"]     = nz_therm
        ds.attrib["thermal_dz_m"]   = dz_therm
        ds.attrib["thermal_depth_m"] = nz_therm * dz_therm
    end

    # --- Provenance variables ---
    # Written only when ProvenanceConfig is provided.  Dimensions: (x, fraction,
    # source, time).  The "source" dimension records one slot per named source.
    # Only the surface layer (layer 1) is tracked.
    has_prov = provenance !== nothing && prov_state !== nothing
    if has_prov
        ns = provenance.n_sources
        defDim(ds, "source", ns)
        srcv = defVar(ds, "source", String, ("source",); attrib = [
            "long_name" => "sediment source region label",
        ])
        for s in 1:ns
            srcv[s] = provenance.source_labels[s]
        end
        defVar(ds, "bed_mass_source", Float64,
               ("x", "fraction", "source", "time");
               attrib = _attrs("kg m-2",
                               "surface-layer bed mass by source region " *
                               "(layer 1 only; sum over source = bed_mass[:,1,:,it])"))
        ds.attrib["provenance_n_sources"]    = ns
        ds.attrib["provenance_source_labels"] = join(provenance.source_labels, ",")
    end

    return NetcdfWriter(ds, 0, -Inf, Float64(output_interval_s),
                        has_mf, has_therm,
                        write_composition, write_T_profile, write_transport,
                        has_prov)
end

"""
    maybe_write_step!(writer, state, t)

Write a slice only if the interval since the last write is ≥
`writer.output_interval_s`. When `output_interval_s == 0` this is
equivalent to `write_step!`.
"""
function maybe_write_step!(w::NetcdfWriter, state::CshoreState, t::Float64;
                           prov_state::Union{ProvenanceState,Nothing}=nothing)
    if t - w.last_write_time >= w.output_interval_s - 1e-9
        write_step!(w, state, t; prov_state=prov_state)
    end
    return nothing
end

"""
    write_step!(writer, state, t)

Unconditionally append one time slice at simulation time `t`. Increments
`writer.it` and writes every registered variable for that slice.

Called once at `t = 0` for the initial state, then at the end of each BC
window (if `output_interval_s == 0`) or at the configured interval.
"""
function write_step!(w::NetcdfWriter, state::CshoreState, t::Float64;
                     prov_state::Union{ProvenanceState,Nothing}=nothing)
    w.it += 1
    it = w.it
    nn = Int(w.ds.dim["x"])

    w.ds["time"][it]       = t
    w.ds["zb"][:, it]      = view(state.zb,      1:nn, 1)
    w.ds["hrms"][:, it]    = view(state.hrms,    1:nn)
    w.ds["h"][:, it]       = view(state.h,       1:nn)
    w.ds["wsetup"][:, it]  = view(state.wsetup,  1:nn)
    w.ds["sigma"][:, it]   = view(state.sigma,   1:nn)
    w.ds["umean"][:, it]   = view(state.umean,   1:nn)
    w.ds["ustd"][:, it]    = view(state.ustd,    1:nn)
    w.ds["vmean"][:, it]   = view(state.vmean,   1:nn)
    w.ds["qbreak"][:, it]  = view(state.qbreak,  1:nn)
    w.ds["q_total"][:, it]  = view(state.q_total, 1:nn)
    w.ds["hrms_ig"][:, it]                = view(state.hrms_ig,                1:nn)
    w.ds["gw_eta"][:, it]                 = view(state.gw_eta,                  1:nn, 1)
    w.ds["theta"][:, it]                  = view(state.theta,                    1:nn, 1)
    w.ds["cohesive_bed_mass"][:, it]      = view(state.cohesive_bed_mass,       1:nn)
    w.ds["cohesive_concentration"][:, it] = view(state.cohesive_concentration,  1:nn)

    # Runup / TWL: prefer jdry (swash limit, iover=1), fall back to jswl
    # (still-water shoreline proxy, iover=0).  Both are 0 at t=0 before
    # the first wave solve, in which case we write NaN.
    let j_run = state.jdry > 0 ? state.jdry :
                (!isempty(state.jswl) && state.jswl[1] > 0 ? state.jswl[1] : 0)
        if j_run > 0 && j_run <= nn
            w.ds["runup_x"][it] = state.xb[j_run]
            w.ds["runup_z"][it] = state.zb[j_run, 1]
        else
            w.ds["runup_x"][it] = NaN
            w.ds["runup_z"][it] = NaN
        end
    end
    w.ds["overtopping_flux"][it] = state.qotf

    if w.write_transport
        w.ds["qbx"][:, :, it] = view(state.qbx, 1:nn, :)
        w.ds["qsx"][:, :, it] = view(state.qsx, 1:nn, :)
    end

    if w.has_multifraction && w.write_composition
        # Write the FULL (x, layer, fraction) slab
        nf = size(state.bed_mass, 3)
        nlayers = size(state.bed_mass, 2)
        w.ds["bed_mass"][:, :, :, it] = view(state.bed_mass, 1:nn, :, :)

        # d50_surface: mass-weighted mean of fraction grain sizes over LAYER 1
        # d50_bulk:    mass-weighted mean over ALL layers
        grain_sizes = w.ds["fraction"][:]
        d50_surf = Vector{Float64}(undef, nn)
        d50_bulk = Vector{Float64}(undef, nn)
        @inbounds for j in 1:nn
            ts = 0.0; ws = 0.0
            tb = 0.0; wb = 0.0
            for k in 1:nf
                m_surf = max(0.0, state.bed_mass[j, 1, k])  # defensive clamp
                ts += m_surf
                ws += m_surf * grain_sizes[k]
                for ilay in 1:nlayers
                    m_l = max(0.0, state.bed_mass[j, ilay, k])
                    tb += m_l
                    wb += m_l * grain_sizes[k]
                end
            end
            d50_surf[j] = ts > 0 ? ws / ts : grain_sizes[1]
            d50_bulk[j] = tb > 0 ? wb / tb : grain_sizes[1]
        end
        w.ds["d50_surface"][:, it] = d50_surf
        w.ds["d50_bulk"][:, it]    = d50_bulk
    end

    # Thermal fields: ALT, zb_hard, surface temperature (always); full profile optional
    if w.has_thermal && state.thermal !== nothing
        thstate = state.thermal
        w.ds["ALT"][:, it]     = view(thstate.ALT, 1:nn)
        w.ds["zb_hard"][:, it] = view(state.zb_hard, 1:nn, 1)
        T_surf = Vector{Float64}(undef, nn)
        @inbounds for j in 1:nn
            T_surf[j] = thstate.T[1, j]  # T is (nz, nn)
        end
        w.ds["T_surface"][:, it] = T_surf
        if w.write_T_profile
            nz = size(thstate.T, 1)
            w.ds["T_profile"][:, :, it] = permutedims(view(thstate.T, 1:nz, 1:nn), (2, 1))
        end
    end

    # Provenance: bed_mass_source(x, fraction, source, time)
    if w.has_provenance && prov_state !== nothing
        nf_prov = size(prov_state.bed_mass_src, 2)
        ns_prov = size(prov_state.bed_mass_src, 3)
        # Write the (nn, nf, ns) slab for the current time slice
        w.ds["bed_mass_source"][:, :, :, it] =
            view(prov_state.bed_mass_src, 1:nn, :, :)
    end

    w.last_write_time = t
    NCDatasets.sync(w.ds)
    return nothing
end

"""
    close_netcdf!(writer)

Close the underlying NCDataset. Safe to call more than once.
"""
function close_netcdf!(w::NetcdfWriter)
    close(w.ds)
    return nothing
end
