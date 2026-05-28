# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
io/cshorejl.jl — Native CSHORE.jl TOML input format.

Reads a `.cshore` TOML file and builds a `CshoreConfig`. This format enables
all CSHORE.jl capabilities that neither the FORTRAN .infile nor XBeach
params.txt can express: spatial grain-size distributions, aeolian transport,
thermal/permafrost, current forcing, etc.

External file paths in the TOML are resolved relative to the `.cshore` file's
directory. CSV files use standard headers (column names on first row).

## Minimal example
```toml
[grid]
dx = 1.0

[bathymetry]
file = "bathy.csv"     # columns: x, z  (optionally: x, z, friction)

[boundary]
bc_file = "waves.csv"  # columns: time, hrms, tp, swl  (wangle optional)

[sediment]
d50 = 3e-4
```

## Full reference — see docs for field descriptions
==============================================================================#

using TOML

"""
    read_cshorejl(path::AbstractString) -> CshoreConfig

Parse a CSHORE.jl native `.cshore` TOML file and return a validated `CshoreConfig`.

All file references inside the TOML are resolved relative to the directory
containing `path`.

# Example
```julia
cfg = read_cshorejl("runs/storm/storm.cshore")
state = run_simulation!(cfg; outfile="storm_out.nc")
```
"""
function read_cshorejl(path::AbstractString)
    isfile(path) || throw(ArgumentError("read_cshorejl: file not found: $path"))
    d = TOML.parsefile(path)
    base = dirname(abspath(path))

    # ── [options] ──────────────────────────────────────────────────────────
    opt_d  = get(d, "options", Dict())
    options = _build_option_flags(opt_d)

    # ── [grid] ─────────────────────────────────────────────────────────────
    grid_d = get(d, "grid", Dict())
    dx = Float64(_req(grid_d, "dx", "grid.dx"))

    # ── [bathymetry] ───────────────────────────────────────────────────────
    bathy_d = get(d, "bathymetry", Dict())
    x_bathy, z_bathy, friction_vec = _load_bathymetry(bathy_d, base)
    friction_scalar = get(bathy_d, "friction", nothing)
    friction_arg = if friction_scalar !== nothing && friction_vec === nothing
        Float64(friction_scalar)
    elseif friction_vec !== nothing
        friction_vec
    else
        0.002   # CSHORE default fw
    end

    # Hardbottom / porous layer bottom
    hb_file = get(bathy_d, "hardbottom_file", nothing)
    hardbottom_z = hb_file !== nothing ? _load_1col(joinpath(base, hb_file), "hardbottom_z") : nothing

    # ── [sediment] ─────────────────────────────────────────────────────────
    sed_d  = get(d, "sediment", Dict())
    sediment = _build_sediment_config(sed_d)

    # ── [multifraction] ────────────────────────────────────────────────────
    mf_d   = get(d, "multifraction", Dict())
    multifraction, fractions_spatial = _build_multifraction(mf_d, base, length(x_bathy))

    # ── [boundary] ─────────────────────────────────────────────────────────
    bc_d   = get(d, "boundary", Dict())
    timebc, hrmsbc, tpbc, swlbc, wangbc, w10, wangle, windcd =
        _load_boundary(bc_d, base)

    # ── [vegetation] ───────────────────────────────────────────────────────
    veg_d  = get(d, "vegetation", nothing)
    vegetation = veg_d !== nothing ? _load_vegetation(veg_d, base, x_bathy) : nothing

    # ── [porous] ───────────────────────────────────────────────────────────
    por_d  = get(d, "porous", nothing)
    porous = por_d !== nothing ? _load_porous(por_d, base, x_bathy) : nothing

    # ── [thermal] ──────────────────────────────────────────────────────────
    therm_d = get(d, "thermal", nothing)
    thermal, thermal_bc = nothing, nothing
    if therm_d !== nothing && get(therm_d, "active", true)
        thermal    = _build_thermal_config(therm_d)
        thermal_bc = _load_thermal_bc(therm_d, base, timebc)
    end

    # ── [aeolian] ──────────────────────────────────────────────────────────
    aeol_d  = get(d, "aeolian", nothing)
    aeolian = nothing
    if aeol_d !== nothing && get(aeol_d, "active", true)
        aeolian = _build_aeolian_config(aeol_d)
        if options.iaeolian == 0
            options = OptionFlags(; (f => getfield(options, f) for f in fieldnames(OptionFlags) if f != :iaeolian)...,
                                    iaeolian=1)
        end
    end

    # ── [diffusion] ────────────────────────────────────────────────────────
    diff_d   = get(d, "diffusion", nothing)
    diffusion = diff_d !== nothing && get(diff_d, "active", true) ?
                _build_diffusion_config(diff_d) : nothing

    # ── [tidal] ────────────────────────────────────────────────────────────
    tide_d  = get(d, "tidal", nothing)
    tidal   = tide_d !== nothing ? _load_tidal_input(tide_d, base) : nothing

    # ── [current] ──────────────────────────────────────────────────────────
    curr_d  = get(d, "current", nothing)
    current = curr_d !== nothing ? _load_current_input(curr_d, base) : nothing

    # ── [overtopping] ──────────────────────────────────────────────────────
    ot_d  = get(d, "overtopping", Dict())
    rcrest = get(ot_d, "rcrest", NaN)

    # ── [gamma] shorthand ──────────────────────────────────────────────────
    gamma_val       = Float64(get(d, "gamma", get(grid_d, "gamma", 0.78)))
    gamma_method    = Symbol(get(d, "gamma_method", "constant"))
    gamma_a         = Float64(get(d, "gamma_a", 0.76))
    gamma_b         = Float64(get(d, "gamma_b", 0.29))
    gamma_min       = Float64(get(d, "gamma_min", 0.35))
    gamma_max       = Float64(get(d, "gamma_max", 0.90))
    facSK           = Float64(get(d, "facSK", 1.0))
    facAS           = Float64(get(d, "facAS", 0.0))
    max_dzb         = Float64(get(d, "max_dzb_per_step", 0.1))
    breaker_delay   = Float64(get(d, "breaker_delay",    0.0))

    # Apply spatial fractions if loaded
    mf_final = if !isempty(fractions_spatial)
        MultifractionConfig(;
            (f => getfield(multifraction, f) for f in fieldnames(MultifractionConfig)
             if f != :initial_fractions_spatial)...,
            initial_fractions_spatial=fractions_spatial,
        )
    else
        multifraction
    end

    build_config(;
        dx=dx,
        bathymetry_x=x_bathy,
        bathymetry_z=z_bathy,
        friction=friction_arg,
        hardbottom_z=hardbottom_z,
        timebc=timebc, hrmsbc=hrmsbc, tpbc=tpbc,
        swlbc=swlbc, wangbc=wangbc,
        w10=isempty(w10) ? nothing : w10,
        wangle=isempty(wangle) ? nothing : wangle,
        windcd=isempty(windcd) ? nothing : windcd,
        options=options,
        sediment=sediment,
        multifraction=mf_final,
        vegetation=vegetation,
        porous=porous,
        thermal=thermal,
        aeolian=aeolian,
        diffusion=diffusion,
        tidal=tidal,
        current=current,
        rcrest=Float64(rcrest),
        gamma=gamma_val, gamma_method=gamma_method,
        gamma_a=gamma_a, gamma_b=gamma_b,
        gamma_min=gamma_min, gamma_max=gamma_max,
        facSK=facSK, facAS=facAS,
        max_dzb_per_step=max_dzb,
        breaker_delay=breaker_delay,
        T_air=thermal_bc !== nothing ? thermal_bc.T_air : nothing,
        T_water=thermal_bc !== nothing ? thermal_bc.T_water : nothing,
        thermal_time=thermal_bc !== nothing ? thermal_bc.time : nothing,
        snow_depth=thermal_bc !== nothing && !isempty(thermal_bc.snow_depth) ?
                   thermal_bc.snow_depth : nothing,
    )
end

# ── Option flags ────────────────────────────────────────────────────────────

function _build_option_flags(d::Dict)
    _i(k, def) = Int(get(d, k, def))
    OptionFlags(;
        iprofl       = _i("iprofl",       1),
        iangle       = _i("iangle",       0),
        iroll        = _i("iroll",        0),
        iwind        = _i("iwind",        0),
        iperm        = _i("iperm",        0),
        iover        = _i("iover",        0),
        iwcint       = _i("iwcint",       0),
        isedav       = _i("isedav",       0),
        iwtran       = _i("iwtran",       0),
        ilab         = _i("ilab",         0),
        infilt       = _i("infilt",       0),
        ipond        = _i("ipond",        0),
        itide        = _i("itide",        0),
        iline        = _i("iline",        1),
        iqydy        = _i("iqydy",        0),
        iveg         = _i("iveg",         0),
        iclay        = _i("iclay",        0),
        ismooth      = _i("ismooth",      1),
        idiss        = _i("idiss",        0),
        ifv          = _i("ifv",          0),
        iweibull     = _i("iweibull",     0),
        iasym        = _i("iasym",        0),
        iwcint_along = _i("iwcint_along", 0),
        iv_transport = _i("iv_transport", 0),
        iaeolian     = _i("iaeolian",     0),
        icurrent     = _i("icurrent",     0),
    )
end

# ── Bathymetry ──────────────────────────────────────────────────────────────

function _load_bathymetry(d::Dict, base::String)
    if haskey(d, "file")
        path = joinpath(base, d["file"])
        return _read_bathy_csv(path)
    elseif haskey(d, "x") && haskey(d, "z")
        x = Float64.(d["x"])
        z = Float64.(d["z"])
        fw = haskey(d, "friction_vector") ? Float64.(d["friction_vector"]) : nothing
        return x, z, fw
    else
        throw(ArgumentError("[bathymetry] requires either 'file' or 'x'+'z' arrays"))
    end
end

function _read_bathy_csv(path::AbstractString)
    isfile(path) || throw(ArgumentError("bathymetry file not found: $path"))
    lines = filter(!isempty, strip.(readlines(path)))
    # Detect header
    header_line = lines[1]
    has_header = !all(c -> isdigit(c) || c in " \t.-+eE" , header_line)
    data_lines = has_header ? lines[2:end] : lines
    cols = _parse_csv_columns(data_lines)
    ncols = length(cols)
    ncols >= 2 || throw(ArgumentError("bathymetry file $path must have at least 2 columns (x, z)"))
    x  = Float64.(cols[1])
    z  = Float64.(cols[2])
    fw = ncols >= 3 ? Float64.(cols[3]) : nothing
    return x, z, fw
end

# ── Sediment ────────────────────────────────────────────────────────────────

function _build_sediment_config(d::Dict)
    _f(k, def) = Float64(get(d, k, def))
    make_sediment(;
        d50     = _f("d50",     0.3e-3),
        sg      = _f("sg",      2.65),
        rho_water = _f("rho_water", 1025.0),
        sporo   = _f("porosity", 0.4),
        shield  = _f("shield",  0.05),
        blp     = _f("blp",     2e-3),
        tanphi  = _f("tanphi",  0.63),
        effb    = _f("effb",    0.005),
        efff    = _f("efff",    0.005),
        slp     = _f("slp",     0.2),
        slpot   = _f("slpot",   0.1),
        wf      = haskey(d, "wf") ? Float64(d["wf"]) : nothing,
    )
end

# ── Multifraction ───────────────────────────────────────────────────────────

function _build_multifraction(d::Dict, base::String, n_bathy::Int)
    isempty(d) && return MultifractionConfig(), zeros(0, 0)

    _f(k, def) = Float64(get(d, k, def))
    _i(k, def) = Int(get(d, k, def))
    _b(k, def) = Bool(get(d, k, def))

    gs_raw = get(d, "grain_sizes", nothing)
    grain_sizes = gs_raw !== nothing ? Float64.(gs_raw) : [_f("d50", 0.3e-3)]
    nf = length(grain_sizes)

    frac_raw = get(d, "initial_fractions", nothing)
    init_fracs = frac_raw !== nothing ? Float64.(frac_raw) :
                 fill(1.0/nf, nf)

    tf_str = get(d, "transport_formula", "original")
    transport_formula = Symbol(tf_str)

    mf = MultifractionConfig(;
        grain_sizes = grain_sizes,
        nlayers     = _i("nlayers", 3),
        layer_thickness = _f("layer_thickness", 0.1),
        porosity    = _f("porosity", 0.4),
        initial_fractions = init_fracs,
        transport_formula = transport_formula,
        use_size_dependent_shields = _b("use_size_dependent_shields", true),
        use_hiding_exposure = _b("use_hiding_exposure", false),
        use_grainsize_tadapt = _b("use_grainsize_tadapt", true),
        n_face_flux_smooth = _i("n_face_flux_smooth", 3),
        n_pickup_smooth    = _i("n_pickup_smooth",    10),
        n_composition_smooth = _i("n_composition_smooth", 3),
    )

    # Spatial fractions: fractions_file = "fracs.csv"  (n_bathy × nf)
    fracs_spatial = zeros(0, 0)
    fracs_file = get(d, "fractions_file", nothing)
    if fracs_file !== nothing
        fp = joinpath(base, fracs_file)
        isfile(fp) || throw(ArgumentError("multifraction.fractions_file not found: $fp"))
        fracs_spatial = _read_fraction_csv(fp, n_bathy, nf)
    end

    return mf, fracs_spatial
end

function _read_fraction_csv(path::AbstractString, n_bathy::Int, nf::Int)
    lines = filter(!isempty, strip.(readlines(path)))
    has_header = !all(c -> isdigit(c) || c in " \t.-+eE", lines[1])
    data_lines = has_header ? lines[2:end] : lines
    cols = _parse_csv_columns(data_lines)
    length(cols) >= nf || throw(ArgumentError(
        "fractions_file $path has $(length(cols)) columns but nf=$nf fractions expected"))
    n = length(cols[1])
    mat = zeros(n, nf)
    for k in 1:nf
        mat[:, k] = Float64.(cols[k])
    end
    return mat
end

# ── Boundary conditions ─────────────────────────────────────────────────────

function _load_boundary(d::Dict, base::String)
    w10 = Float64[]; wangle = Float64[]; windcd = Float64[]

    if haskey(d, "bc_file")
        fp = joinpath(base, d["bc_file"])
        isfile(fp) || throw(ArgumentError("boundary.bc_file not found: $fp"))
        timebc, hrmsbc, tpbc, swlbc, wangbc = _read_bc_csv(fp)
    elseif haskey(d, "timebc")
        timebc  = Float64.(d["timebc"])
        hrmsbc  = Float64.(d["hrmsbc"])
        tpbc    = Float64.(d["tpbc"])
        swlbc   = haskey(d, "swlbc") ? Float64.(d["swlbc"]) : zeros(length(timebc))
        wangbc  = haskey(d, "wangbc") ? Float64.(d["wangbc"]) : zeros(length(timebc))
    else
        throw(ArgumentError("[boundary] requires either 'bc_file' or inline 'timebc'+'hrmsbc'+'tpbc' arrays"))
    end

    if haskey(d, "wind_file")
        fp = joinpath(base, d["wind_file"])
        isfile(fp) || throw(ArgumentError("boundary.wind_file not found: $fp"))
        w10, wangle, windcd = _read_wind_csv(fp)
    elseif haskey(d, "w10")
        w10     = Float64.(d["w10"])
        wangle  = haskey(d, "wangle") ? Float64.(d["wangle"]) : zeros(length(w10))
        windcd  = haskey(d, "windcd") ? Float64.(d["windcd"]) : fill(0.0015, length(w10))
    end

    return timebc, hrmsbc, tpbc, swlbc, wangbc, w10, wangle, windcd
end

function _read_bc_csv(path::AbstractString)
    lines = filter(!isempty, strip.(readlines(path)))
    has_header = !all(c -> isdigit(c) || c in " \t.-+eE", lines[1])
    header = has_header ? lowercase.(split(replace(lines[1], ","=>" "))) : nothing
    data_lines = has_header ? lines[2:end] : lines
    cols = _parse_csv_columns(data_lines)
    ncols = length(cols)
    ncols >= 3 || throw(ArgumentError("BC file $path must have at least 3 columns (time, hrms, tp)"))

    # Resolve columns by header name or position
    col_idx = _col_index(header, ncols)
    timebc = Float64.(cols[col_idx[:time]])
    hrmsbc = Float64.(cols[col_idx[:hrms]])
    tpbc   = Float64.(cols[col_idx[:tp]])
    swlbc  = ncols >= 4 && haskey(col_idx, :swl) ?
             Float64.(cols[col_idx[:swl]]) : zeros(length(timebc))
    wangbc = ncols >= 5 && haskey(col_idx, :wangle) ?
             Float64.(cols[col_idx[:wangle]]) : zeros(length(timebc))
    return timebc, hrmsbc, tpbc, swlbc, wangbc
end

function _col_index(header, ncols)
    d = Dict{Symbol,Int}()
    if header === nothing
        # positional: time hrms tp swl wangle
        d[:time] = 1; d[:hrms] = 2; d[:tp] = 3
        ncols >= 4 && (d[:swl] = 4)
        ncols >= 5 && (d[:wangle] = 5)
    else
        for (i, h) in enumerate(header)
            h = replace(h, ","=>"", " "=>"")
            if h in ("time", "t")           d[:time]   = i end
            if h in ("hrms", "hm0", "hs", "h") d[:hrms] = i end
            if h in ("tp", "period", "trep")   d[:tp]   = i end
            if h in ("swl", "eta", "zs")       d[:swl]  = i end
            if h in ("wangle", "angle", "dir") d[:wangle] = i end
        end
        haskey(d, :time) || (d[:time] = 1)
        haskey(d, :hrms) || (d[:hrms] = 2)
        haskey(d, :tp)   || (d[:tp]   = 3)
    end
    return d
end

function _read_wind_csv(path::AbstractString)
    lines = filter(!isempty, strip.(readlines(path)))
    has_header = !all(c -> isdigit(c) || c in " \t.-+eE", lines[1])
    data_lines = has_header ? lines[2:end] : lines
    cols = _parse_csv_columns(data_lines)
    ncols = length(cols)
    ncols >= 2 || throw(ArgumentError("wind file $path must have at least 2 columns (time, w10)"))
    # skip time column, read w10, wangle, windcd
    w10    = Float64.(cols[2])
    wangle = ncols >= 3 ? Float64.(cols[3]) : zeros(length(w10))
    windcd = ncols >= 4 ? Float64.(cols[4]) : fill(0.0015, length(w10))
    return w10, wangle, windcd
end

# ── Vegetation ──────────────────────────────────────────────────────────────

function _load_vegetation(d::Dict, base::String, x_bathy::Vector{Float64})
    vegcd  = Float64(get(d, "vegcd",  1.0))
    vegcdm = Float64(get(d, "vegcdm", 1.0))
    if haskey(d, "veg_file")
        fp = joinpath(base, d["veg_file"])
        isfile(fp) || throw(ArgumentError("vegetation.veg_file not found: $fp"))
        return _read_vegetation_file(fp, vegcd, vegcdm, x_bathy)
    else
        # Single-value vegetation over entire profile
        vegn_val = Float64(get(d, "vegn", 100.0))
        vegb_val = Float64(get(d, "vegb", 0.01))
        vegd_val = Float64(get(d, "vegd", 0.3))
        np = length(x_bathy)
        return VegetationInput(;
            vegcd=vegcd, vegcdm=vegcdm,
            vegn=fill(vegn_val, np, 1),
            vegb=fill(vegb_val, np, 1),
            vegd=fill(vegd_val, np, 1),
            vegh=fill(vegd_val, np, 1),
        )
    end
end

function _read_vegetation_file(path, vegcd, vegcdm, x_ref)
    lines = filter(!isempty, strip.(readlines(path)))
    has_header = !all(c -> isdigit(c) || c in " \t.-+eE", lines[1])
    data_lines = has_header ? lines[2:end] : lines
    cols = _parse_csv_columns(data_lines)
    ncols = length(cols)
    ncols >= 4 || throw(ArgumentError("veg_file $path must have at least 4 columns: x, vegn, vegb, vegd"))
    np = length(cols[1])
    vegn = reshape(Float64.(cols[2]), np, 1)
    vegb = reshape(Float64.(cols[3]), np, 1)
    vegd = reshape(Float64.(cols[4]), np, 1)
    vegh = ncols >= 5 ? reshape(Float64.(cols[5]), np, 1) : copy(vegd)
    return VegetationInput(; vegcd=vegcd, vegcdm=vegcdm,
                             vegn=vegn, vegb=vegb, vegd=vegd, vegh=vegh)
end

# ── Porous layer ─────────────────────────────────────────────────────────────

function _load_porous(d::Dict, base::String, x_bathy::Vector{Float64})
    porosity = Float64(get(d, "porosity", 0.4))
    dn50     = Float64(get(d, "stone_diameter", 0.02))
    if haskey(d, "porous_file")
        fp = joinpath(base, d["porous_file"])
        isfile(fp) || throw(ArgumentError("porous.porous_file not found: $fp"))
        lines = filter(!isempty, strip.(readlines(fp)))
        has_hdr = !all(c -> isdigit(c) || c in " \t.-+eE", lines[1])
        cols = _parse_csv_columns(has_hdr ? lines[2:end] : lines)
        length(cols) >= 2 || throw(ArgumentError("porous_file must have 2 columns: x, z_porous"))
        xp = Float64.(cols[1]); zp = Float64.(cols[2])
        return PorousInput(xp, zp; porosity=porosity, stone_diameter=dn50)
    else
        zp_uniform = Float64(get(d, "z_porous", -10.0))
        return PorousInput(x_bathy, fill(zp_uniform, length(x_bathy));
                           porosity=porosity, stone_diameter=dn50)
    end
end

# ── Thermal ──────────────────────────────────────────────────────────────────

function _build_thermal_config(d::Dict)
    _i(k, def) = Int(get(d, k, def))
    _f(k, def) = Float64(get(d, k, def))
    ThermalConfig(;
        nz           = _i("nz",        30),
        dz           = _f("dz",        0.1),
        k_frozen     = _f("k_frozen",  1.5),
        k_thawed     = _f("k_thawed",  0.8),
        C_frozen     = _f("C_frozen",  1.8e6),
        C_thawed     = _f("C_thawed",  3.0e6),
        L            = _f("L",         3.34e8),
        moisture     = _f("moisture",  0.35),
        T_init       = _f("T_init",   -5.0),
        T_lower      = _f("T_lower",  -5.0),
        z_bottom     = _f("z_bottom", -Inf),
        cfl_safety   = _f("cfl_safety", 0.4),
    )
end

function _load_thermal_bc(d::Dict, base::String, timebc_ref::Vector{Float64})
    if haskey(d, "thermal_bc_file")
        fp = joinpath(base, d["thermal_bc_file"])
        isfile(fp) || throw(ArgumentError("thermal.thermal_bc_file not found: $fp"))
        lines = filter(!isempty, strip.(readlines(fp)))
        has_hdr = !all(c -> isdigit(c) || c in " \t.-+eE", lines[1])
        cols = _parse_csv_columns(has_hdr ? lines[2:end] : lines)
        length(cols) >= 3 || throw(ArgumentError("thermal_bc_file must have 3 columns: time, T_air, T_water"))
        t = Float64.(cols[1]); Ta = Float64.(cols[2]); Tw = Float64.(cols[3])
        sd = length(cols) >= 4 ? Float64.(cols[4]) : Float64[]
        return ThermalBoundaryTimeSeries(; time=t, T_air=Ta, T_water=Tw, snow_depth=sd)
    elseif haskey(d, "T_air") && haskey(d, "T_water")
        Ta = Float64.(d["T_air"]); Tw = Float64.(d["T_water"])
        t  = haskey(d, "time") ? Float64.(d["time"]) : timebc_ref
        return ThermalBoundaryTimeSeries(; time=t, T_air=Ta, T_water=Tw)
    else
        @warn "read_cshorejl: [thermal] active but no thermal_bc_file or T_air/T_water provided — using T_init everywhere"
        return nothing
    end
end

# ── Aeolian ──────────────────────────────────────────────────────────────────

function _build_aeolian_config(d::Dict)
    _f(k, def) = Float64(get(d, k, def))
    z_contour    = _f("z_contour", NaN)
    decay_length = _f("decay_length", 0.0)
    vegetation = ContourVegetation(; z_contour=z_contour, decay_length=decay_length)
    AeolianConfig(;
        vegetation          = vegetation,
        Ck                  = _f("Ck",                 1.8),
        z_meas              = _f("z_meas",              10.0),
        karman              = _f("karman",              0.4),
        z0_factor           = _f("z0_factor",           2.0/30.0),
        rho_air             = _f("rho_air",             1.225),
        moisture            = _f("moisture",            0.0),
        runup_buffer        = _f("runup_buffer",        0.0),
        saturation_length   = _f("saturation_length",   5.0),
        dune_decay_length   = _f("dune_decay_length",   5.0),
        dt_aeolian_max      = _f("dt_aeolian_max",      Inf),
    )
end

# ── Diffusion ────────────────────────────────────────────────────────────────

function _build_diffusion_config(d::Dict)
    _f(k, def) = Float64(get(d, k, def))
    DiffusionConfig(;
        D_base           = _f("D_base",           0.01),
        critical_slope   = _f("critical_slope",   0.7),
        slope_exponent   = _f("slope_exponent",   2.0),
        critical_factor  = _f("critical_factor",  10.0),
        max_diffusivity  = _f("max_diffusivity",  1000.0),
        thermal_control  = Bool(get(d, "thermal_control", true)),
        thaw_threshold   = _f("thaw_threshold",   0.01),
    )
end

# ── Tidal input ──────────────────────────────────────────────────────────────

function _load_tidal_input(d::Dict, base::String)
    if haskey(d, "file")
        fp = joinpath(base, d["file"])
        lines = filter(!isempty, strip.(readlines(fp)))
        has_hdr = !all(c -> isdigit(c) || c in " \t.-+eE", lines[1])
        cols = _parse_csv_columns(has_hdr ? lines[2:end] : lines)
        t = Float64.(cols[1]); detady = Float64.(cols[2])
        dswldt = length(cols) >= 3 ? Float64.(cols[3]) : Float64[]
        return TidalInput(; time=t, detady=detady, dswldt=dswldt)
    else
        t      = Float64.(d["time"])
        detady = Float64.(d["detady"])
        dswldt = haskey(d, "dswldt") ? Float64.(d["dswldt"]) : Float64[]
        return TidalInput(; time=t, detady=detady, dswldt=dswldt)
    end
end

# ── Current input ────────────────────────────────────────────────────────────

function _load_current_input(d::Dict, base::String)
    if haskey(d, "file")
        fp = joinpath(base, d["file"])
        lines = filter(!isempty, strip.(readlines(fp)))
        has_hdr = !all(c -> isdigit(c) || c in " \t.-+eE", lines[1])
        cols = _parse_csv_columns(has_hdr ? lines[2:end] : lines)
        return CurrentInput(; time=Float64.(cols[1]), vbc=Float64.(cols[2]))
    else
        return CurrentInput(; time=Float64.(d["time"]), vbc=Float64.(d["vbc"]))
    end
end

# ── CSV / column parsing helpers ─────────────────────────────────────────────

"""Parse a vector of data line strings into column vectors (handles comma or whitespace separation)."""
function _parse_csv_columns(lines::Vector{<:AbstractString})
    isempty(lines) && return Vector{Float64}[]
    # Detect separator: comma present → CSV, else whitespace
    use_comma = any(l -> ',' in l, lines)
    splitter = use_comma ? r"[,\s]+" : r"\s+"
    parsed = [split(strip(l), splitter) for l in lines if !isempty(strip(l))]
    ncols = maximum(length(r) for r in parsed)
    cols = [Float64[] for _ in 1:ncols]
    for row in parsed
        for (i, v) in enumerate(row)
            isempty(v) && continue
            push!(cols[i], parse(Float64, v))
        end
    end
    return cols
end

function _load_1col(path::AbstractString, what::String)
    isfile(path) || throw(ArgumentError("$what file not found: $path"))
    vals = Float64[]
    for line in eachline(path)
        for tok in split(line)
            isempty(tok) && continue
            push!(vals, parse(Float64, tok))
        end
    end
    return vals
end

_req(d, key, label) = haskey(d, key) ? d[key] :
    throw(ArgumentError("read_cshorejl: missing required field '$label'"))
