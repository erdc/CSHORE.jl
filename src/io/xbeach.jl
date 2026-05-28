# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
io/xbeach.jl — XBeach params.txt → CshoreConfig adapter.

Reads an XBeach run directory (containing params.txt and associated grid /
boundary files) and returns a CshoreConfig that approximates the XBeach
setup as closely as CSHORE physics allow.

Supported:
  - 1D profiles (ny=0): depfile + xfile (uniform or variable dx)
  - Grain size: D50, D90 (D90 triggers two-fraction mode)
  - Wave BCs: stationary (instat=0/stat), JONSWAP file (instat=jons/jons_table),
    energy time series (wbctype=ts_1 gen.ezs format)
  - Morphology: sedtrans / morphology flags → iprofl
  - Friction: manning / chezy / quadratic → fw conversion
  - Breaking: gamma / alpha → CshoreConfig

Not mapped (silently ignored — XBeach-specific, no CSHORE equivalent):
  - 2D grids (ny>0 warns but takes first cross-shore transect)
  - morfac (morphological acceleration): warns if ≠ 1
  - XBeach flow/roller solver parameters (eps, ARC, nuhfac, order, …)
  - Non-JONSWAP spectral types (instat=2,3, …)
  - Output settings (outputformat, tintm, nglobalvar, …)
  - XBeach physics toggles with no CSHORE analogue (back, front, wbcEvarreduce, …)
  - Anything not explicitly mapped is simply discarded without warning
==============================================================================#

const _RHO_WATER = 1025.0   # kg/m³ for energy→Hrms conversion
const _GRAV_XB   = 9.81     # m/s²

"""
    read_xbeach_params(dir; dt_bc=nothing, strict=false) -> CshoreConfig

Read an XBeach run directory and return a `CshoreConfig`.

`dir` must contain a `params.txt` (or `params_*.txt` — if multiple are found
the user is asked to disambiguate via the file kwarg).

- `dt_bc`  — resample the wave BC time series to this interval (s).
  Defaults to ~100 evenly spaced points. Pass `nothing` to keep the original
  resolution (may be very dense for surfbeat ts_1 BCs).
- `strict` — if `true`, throw on any XBeach parameter that has no CSHORE
  equivalent; if `false` (default), warn and continue.
- `params_file` — override the params.txt filename within `dir`.

# Example
```julia
cfg = read_xbeach_params("examples/benchmarks/xbeach/data/Boers_1C";
                         params_file="params_BOI.txt")
state = run_simulation!(cfg; outfile="boers_out.nc")
```
"""
function read_xbeach_params(dir::AbstractString;
                             dt_bc::Union{Nothing,Float64}=nothing,
                             strict::Bool=false,
                             params_file::Union{Nothing,String}=nothing)
    isdir(dir) || throw(ArgumentError("read_xbeach_params: directory not found: $dir"))

    # --- locate params.txt ---
    pfile = if params_file !== nothing
        joinpath(dir, params_file)
    else
        candidates = filter(f -> startswith(f, "params") && endswith(f, ".txt"),
                            readdir(dir))
        if length(candidates) == 1
            joinpath(dir, candidates[1])
        elseif isempty(candidates)
            throw(ArgumentError("No params*.txt found in $dir"))
        else
            throw(ArgumentError("Multiple params*.txt in $dir: $candidates — pass params_file= to select one"))
        end
    end
    isfile(pfile) || throw(ArgumentError("params file not found: $pfile"))

    p = _parse_xbeach_params(pfile)

    # --- 2D guard ---
    ny = _xb_int(p, "ny", 0)
    if ny > 0
        _xb_warn_or_throw("XBeach ny=$ny — 2D grids not supported; only the first cross-shore transect will be used", strict)
    end

    # --- grid ---
    depfile = _xb_str(p, "depfile", "bed.dep")
    xfile   = _xb_str(p, "xfile",  nothing)
    posdwn  = _xb_int(p, "posdwn", -1)   # -1 = depth positive downward (common XBeach default)

    # Read grid files directly — use file token count rather than nx+1 so that
    # 2D XBeach setups (or params.txt files with stale nx) are handled gracefully.
    z_raw_all = _read_xb_array_all(joinpath(dir, depfile))
    nx     = _xb_int(p, "nx", length(z_raw_all) - 1)
    ny     = _xb_int(p, "ny", 0)
    n_nodes = nx + 1

    # For 2D setups take only the first cross-shore transect (first nx+1 values)
    if length(z_raw_all) > n_nodes
        @warn "read_xbeach_params: $(depfile) has $(length(z_raw_all)) values but nx=$(nx) ($(n_nodes) nodes) — using first $(n_nodes) values as cross-shore transect"
        z_raw = z_raw_all[1:n_nodes]
    else
        n_nodes = length(z_raw_all)
        z_raw = z_raw_all
    end

    # XBeach sign: posdwn=-1 → file stores positive depth (zb = -depth)
    #              posdwn= 0 → file stores bed elevation (positive up, same as CSHORE)
    #              posdwn= 1 → file stores depth positive downward
    z_bathy = if posdwn == -1 || posdwn == 1
        -z_raw    # flip to positive-up elevation
    else
        z_raw
    end

    vardx = _xb_int(p, "vardx", 0)
    x_bathy = if xfile !== nothing
        fp = joinpath(dir, xfile)
        isfile(fp) || throw(ArgumentError("XBeach xfile not found: $fp"))
        x_all = _read_xb_array_all(fp)
        if length(x_all) > n_nodes
            x_all[1:n_nodes]
        elseif length(x_all) < n_nodes
            @warn "xfile has $(length(x_all)) values but n_nodes=$n_nodes — padding with uniform spacing"
            vcat(x_all, x_all[end] .+ (1:(n_nodes - length(x_all))))
        else
            x_all
        end
    else
        xori = _xb_float(p, "xori", 0.0)
        dx_uniform = _xb_float(p, "dx", 1.0)
        xori .+ (0:(n_nodes-1)) .* dx_uniform
    end

    dx_val = if vardx == 1 || xfile !== nothing
        (x_bathy[end] - x_bathy[1]) / max(n_nodes - 1, 1)
    else
        dx_raw = _xb_float(p, "dx", nothing)
        dx_raw !== nothing ? dx_raw : (x_bathy[end] - x_bathy[1]) / max(n_nodes - 1, 1)
    end
    dx_val > 0 || (dx_val = 1.0; @warn "read_xbeach_params: computed dx ≤ 0, defaulting to 1.0 m")

    # Friction conversion
    bedfriction = lowercase(get(p, "bedfriction", "chezy"))
    bedfriccoef = _xb_float(p, "bedfriccoef", 0.01)
    fw = _xbeach_friction_to_fw(bedfriction, bedfriccoef, strict)

    # --- grain size ---
    d50 = _xb_float(p, "d50", 2.5e-4)
    d90 = _xb_float(p, "d90", nothing)

    multifrac = if d90 !== nothing && d90 > d50 * 1.05
        # Set up two-fraction: d50 as fine, d90 as coarse; 80/20 split
        f_fine = 0.8; f_coarse = 0.2
        MultifractionConfig(
            grain_sizes=[d50, d90],
            nlayers=3,
            initial_fractions=[f_fine, f_coarse],
            transport_formula=:size_adaptive,
        )
    else
        MultifractionConfig(grain_sizes=[d50], initial_fractions=[1.0])
    end

    sed = make_sediment(d50=d50)

    # --- SWL and timing ---
    zs0   = _xb_float(p, "zs0",   0.0)
    tstop = _xb_float(p, "tstop", nothing)
    tstop !== nothing || throw(ArgumentError("XBeach params: 'tstop' is required"))

    # --- wave breaking ---
    gamma_val = _xb_float(p, "gamma", 0.55)
    # XBeach alpha is a breaking dissipation coefficient (Roelvink), not directly
    # the CSHORE gamma. We use the XBeach gamma as the breaker index.

    # --- morphology ---
    sedtrans  = _xb_int(p, "sedtrans",  1)
    morphology = _xb_int(p, "morphology", 1)
    iprofl = (sedtrans == 1 && morphology == 1) ? 1 : 0

    morfac = _xb_float(p, "morfac", 1.0)
    if morfac != 1.0
        _xb_warn_or_throw("XBeach morfac=$morfac — morphological acceleration not supported in CSHORE (morfac will be ignored; run time is tstop=$tstop s)", strict)
    end

    # --- wave boundary conditions ---
    timebc, hrmsbc, tpbc, swlbc, wangbc = _xbeach_build_bc(p, dir, tstop, zs0, strict)

    # Optionally resample BC time series to dt_bc
    if dt_bc !== nothing && length(timebc) > 2
        timebc, hrmsbc, tpbc, swlbc, wangbc =
            _resample_bc(timebc, hrmsbc, tpbc, swlbc, wangbc, dt_bc)
    elseif dt_bc === nothing && length(timebc) > 200
        # Auto-thin very dense surfbeat BCs to ~100 points
        target = 100
        stride = max(1, div(length(timebc) - 1, target - 1))
        idx = unique(vcat(1:stride:length(timebc), length(timebc)))
        timebc  = timebc[idx]
        hrmsbc  = hrmsbc[idx]
        tpbc    = tpbc[idx]
        swlbc   = swlbc[idx]
        wangbc  = wangbc[idx]
    end

    bc = BoundaryTimeSeries(; timebc=timebc, tpbc=tpbc, hrmsbc=hrmsbc,
                              wsetbc=zeros(length(timebc)),
                              swlbc=swlbc, wangbc=wangbc)

    # --- build CshoreConfig ---
    opts = OptionFlags(; iprofl=iprofl)
    build_config(;
        dx=dx_val,
        bathymetry_x=collect(Float64, x_bathy),
        bathymetry_z=collect(Float64, z_bathy),
        friction=fw,
        timebc=timebc, hrmsbc=hrmsbc, tpbc=tpbc, swlbc=swlbc, wangbc=wangbc,
        options=opts,
        sediment=sed,
        multifraction=multifrac,
        gamma=gamma_val,
    )
end

# ─── Private helpers ────────────────────────────────────────────────────────

"""Parse an XBeach params.txt file into a Dict{String,String}."""
function _parse_xbeach_params(path::AbstractString)
    params = Dict{String,String}()
    for line in eachline(path)
        # Strip comments (% anywhere on line)
        idx = findfirst('%', line)
        content = idx === nothing ? line : line[1:idx-1]
        content = strip(content)
        isempty(content) && continue
        # Key = value  or  key value
        m = match(r"^(\w+)\s*=\s*(.+)$", content)
        if m !== nothing
            params[lowercase(strip(m[1]))] = strip(m[2])
        else
            # space-separated key value
            parts = split(content, limit=2)
            if length(parts) == 2
                params[lowercase(strip(parts[1]))] = strip(parts[2])
            end
        end
    end
    return params
end

_xb_str(p, key, default) = get(p, key, default === nothing ? nothing : string(default))
_xb_int(p, key, default) = haskey(p, key) ? parse(Int, p[key]) : default
_xb_float(p, key, default) = haskey(p, key) ? parse(Float64, p[key]) : default

function _xb_warn_or_throw(msg, strict)
    if strict
        throw(ArgumentError(msg))
    else
        @warn "read_xbeach_params: $msg"
    end
end

"""Read all whitespace/newline-separated numeric tokens from an XBeach data file."""
function _read_xb_array_all(path::AbstractString)
    isfile(path) || throw(ArgumentError("XBeach data file not found: $path"))
    vals = Float64[]
    for line in eachline(path)
        for tok in split(line)
            isempty(tok) && continue
            v = tryparse(Float64, tok)
            v !== nothing && push!(vals, v)
        end
    end
    return vals
end

# Kept for callers that need exactly expected_n values.
function _read_xb_array(path::AbstractString, expected_n::Int)
    vals = _read_xb_array_all(path)
    if length(vals) != expected_n
        @warn "read_xbeach_params: $path has $(length(vals)) values, expected $expected_n — truncating/padding to $expected_n"
        if length(vals) > expected_n
            resize!(vals, expected_n)
        else
            append!(vals, fill(last(vals), expected_n - length(vals)))
        end
    end
    return vals
end

"""Convert XBeach bedfriction type + coefficient to CSHORE fw (bottom friction factor)."""
function _xbeach_friction_to_fw(frtype::AbstractString, coef::Float64, strict::Bool)
    if frtype == "manning"
        # Manning n → fw: fw ≈ 2g n²/h^(1/3); representative depth h≈2m
        h_rep = 2.0
        return 2.0 * _GRAV_XB * coef^2 / h_rep^(1/3)
    elseif frtype == "chezy"
        # Chezy C → fw: fw = 2g/C²
        return 2.0 * _GRAV_XB / coef^2
    elseif frtype in ("quadratic", "constant", "white_colebrook")
        return coef
    else
        @warn "read_xbeach_params: unknown bedfriction='$frtype', using coef=$coef directly as fw"
        return coef
    end
end

"""
Build CSHORE boundary time series from XBeach params dict and BC files.
Returns (timebc, hrmsbc, tpbc, swlbc, wangbc).
"""
function _xbeach_build_bc(p::Dict, dir::AbstractString,
                           tstop::Float64, zs0::Float64, strict::Bool)
    instat  = get(p, "instat",  nothing)
    wbctype = get(p, "wbctype", nothing)

    # Derive the effective BC type token
    bc_type = if wbctype !== nothing
        lowercase(wbctype)
    elseif instat !== nothing
        lowercase(instat)
    else
        "stat"
    end

    trep = _xb_float(p, "trep", 8.0)
    hm0  = _xb_float(p, "hm0",  nothing)
    bcfile = _xb_str(p, "bcfile", nothing)
    rt     = _xb_float(p, "rt",    3600.0)   # record length per JONSWAP entry (s)
    mainang = _xb_float(p, "mainang", 270.0)  # default shore-normal in nautical degrees
    # Convert nautical wave direction to math angle (degrees from +x)
    # XBeach: 270° nautical = coming from west = toward east = shore-normal for west-facing coast
    wangle = 0.0   # CSHORE normal incidence default

    if bc_type in ("stat", "0", "stationary")
        # Stationary: single bulk parameters
        hrms_val = if hm0 !== nothing
            hm0 / sqrt(2.0)
        else
            _xb_float(p, "hrms", _xb_float(p, "hs", 1.0))
        end
        return [0.0, tstop], [hrms_val, hrms_val], [trep, trep],
               [zs0, zs0], [wangle, wangle]

    elseif bc_type in ("jons", "jons_table", "1")
        # JONSWAP file: each row = Hm0 fp mainang gammajsp s
        bcfile !== nothing || begin
            @warn "read_xbeach_params: instat=jons but no bcfile — using stationary fallback"
            hrms_val = hm0 !== nothing ? hm0/sqrt(2.0) : 1.0
            return [0.0, tstop], [hrms_val, hrms_val], [trep, trep],
                   [zs0, zs0], [wangle, wangle]
        end
        fp = joinpath(dir, bcfile)
        if !isfile(fp)
            # Try bcfile as directory of JONSWAP files (DUROS case)
            if isdir(fp)
                jfiles = sort(filter(f -> !startswith(f, "."), readdir(fp)))
                isempty(jfiles) && throw(ArgumentError("bcfile directory $fp is empty"))
                fp = joinpath(fp, jfiles[1])
                @info "read_xbeach_params: reading first JONSWAP file: $fp"
            else
                throw(ArgumentError("XBeach bcfile not found: $fp"))
            end
        end
        return _parse_jonswap_bc(fp, rt, tstop, zs0, strict)

    elseif bc_type in ("ts_1", "41", "jons_1d")
        # Surfbeat energy time series — gen.ezs format or ts format
        # Look for gen.ezs in bc/ subdirectory or current dir
        ezs_candidates = [
            joinpath(dir, "bc", "gen.ezs"),
            joinpath(dir, "gen.ezs"),
            bcfile !== nothing ? joinpath(dir, bcfile) : "",
        ]
        ezs_file = ""
        for c in ezs_candidates
            !isempty(c) && isfile(c) && (ezs_file = c; break)
        end
        if isempty(ezs_file)
            @warn "read_xbeach_params: ts_1 BC file not found (tried $(ezs_candidates)) — stationary fallback"
            hrms_val = hm0 !== nothing ? hm0/sqrt(2.0) : 1.0
            return [0.0, tstop], [hrms_val, hrms_val], [trep, trep],
                   [zs0, zs0], [wangle, wangle]
        end
        return _parse_gen_ezs_bc(ezs_file, trep, zs0, strict)

    else
        _xb_warn_or_throw("XBeach BC type '$bc_type' not supported — stationary fallback", strict)
        hrms_val = hm0 !== nothing ? hm0/sqrt(2.0) : 1.0
        return [0.0, tstop], [hrms_val, hrms_val], [trep, trep],
               [zs0, zs0], [wangle, wangle]
    end
end

"""
Parse a JONSWAP boundary condition file. Two formats are supported:

1. **Key=value** (single sea state): lines like `Hm0 = 1.5` or `fp = 0.1`
2. **Table** (time series): whitespace-separated rows of `Hm0 fp mainang gammajsp s`,
   one row per record interval `rt`.

Returns (timebc, hrmsbc, tpbc, swlbc, wangbc).
"""
function _parse_jonswap_bc(path::AbstractString, rt::Float64, tstop::Float64,
                            zs0::Float64, strict::Bool)
    # Detect format by checking if any line contains '='
    lines = filter(l -> !isempty(strip(l)) && !startswith(strip(l), '%') &&
                        !startswith(strip(l), '!') && !startswith(strip(l), '#'),
                   readlines(path))
    isempty(lines) && throw(ArgumentError("JONSWAP file $path is empty"))

    is_keyvalue = any(l -> '=' in l, lines)

    if is_keyvalue
        # Key=value single sea state
        kv = Dict{String,Float64}()
        for line in lines
            m = match(r"^\s*(\w+)\s*=\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)", line)
            m !== nothing && (kv[lowercase(m[1])] = parse(Float64, m[2]))
        end
        hm0 = get(kv, "hm0", get(kv, "hs", 1.0))
        fp  = get(kv, "fp",  get(kv, "fpeak", 0.1))
        hrms_val = hm0 / sqrt(2.0)
        tp_val   = 1.0 / fp
        return [0.0, tstop], [hrms_val, hrms_val], [tp_val, tp_val],
               [zs0, zs0], [0.0, 0.0]
    else
        # Table format: each row = Hm0 fp [mainang gammajsp s ...]
        rows = Float64[]
        ncols = 0
        for line in lines
            toks = split(strip(line))
            isempty(toks) && continue
            t1 = tryparse(Float64, toks[1])
            t1 === nothing && continue   # skip non-numeric rows (e.g. stray headers)
            ncols == 0 && (ncols = length(toks))
            for tok in toks[1:min(ncols, length(toks))]
                push!(rows, parse(Float64, tok))
            end
        end
        ncols == 0 && throw(ArgumentError("JONSWAP file $path: no numeric rows found"))
        nrows = div(length(rows), ncols)
        nrows == 0 && throw(ArgumentError("JONSWAP file $path: could not parse any rows"))
        mat = reshape(rows[1:nrows*ncols], ncols, nrows)'

        hm0_arr  = mat[:, 1]
        fp_arr   = ncols >= 2 ? mat[:, 2] : fill(0.1, nrows)
        hrms_arr = hm0_arr ./ sqrt(2.0)
        tp_arr   = 1.0 ./ fp_arr
        t_arr    = collect(Float64, (0:(nrows-1)) .* rt)
        if last(t_arr) < tstop
            push!(t_arr, tstop); push!(hrms_arr, hrms_arr[end]); push!(tp_arr, tp_arr[end])
        end
        return t_arr, collect(hrms_arr), collect(tp_arr), fill(zs0, length(t_arr)),
               zeros(length(t_arr))
    end
end

"""
Parse a gen.ezs surfbeat BC file.
Format: header lines (starting with * or row-count line), then data columns:
  t (s)  eta_LF (m)  E (J/m²)  eta_BI (m)  eta_F (m)

Derives Hrms from E: Hrms = sqrt(8E / (ρg)).
Uses `trep` as the representative wave period (Tp).
"""
function _parse_gen_ezs_bc(path::AbstractString, trep::Float64, zs0::Float64, strict::Bool)
    times  = Float64[]
    energies = Float64[]
    for line in eachline(path)
        s = strip(line)
        isempty(s) && continue
        # Skip header lines
        (startswith(s, '*') || startswith(s, '%') || startswith(s, '!') ||
         startswith(s, '#')) && continue
        toks = split(s)
        length(toks) < 3 && continue
        # Try parsing first token as time — if it fails it's a descriptor line
        t_val = tryparse(Float64, toks[1])
        t_val === nothing && continue
        e_val = tryparse(Float64, toks[3])
        e_val === nothing && continue
        push!(times, t_val)
        push!(energies, e_val)
    end
    isempty(times) && throw(ArgumentError("gen.ezs $path: no numeric rows found"))

    hrms_arr = sqrt.(max.(0.0, 8.0 .* energies ./ (_RHO_WATER * _GRAV_XB)))
    tp_arr   = fill(trep, length(times))
    swl_arr  = fill(zs0, length(times))
    wang_arr = zeros(length(times))

    return times, hrms_arr, tp_arr, swl_arr, wang_arr
end

"""Resample BC arrays to a uniform interval of dt_bc seconds."""
function _resample_bc(t, h, tp, swl, wang, dt_bc::Float64)
    t_new = collect(t[1]:dt_bc:t[end])
    isempty(t_new) && (t_new = [t[1], t[end]])
    last(t_new) < t[end] && push!(t_new, t[end])
    h_new    = [interp1(t, h,    ti) for ti in t_new]
    tp_new   = [interp1(t, tp,   ti) for ti in t_new]
    swl_new  = [interp1(t, swl,  ti) for ti in t_new]
    wang_new = [interp1(t, wang, ti) for ti in t_new]
    return t_new, h_new, tp_new, swl_new, wang_new
end
