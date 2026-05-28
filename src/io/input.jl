# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
io/input.jl — Input file parsing.

Parses a FORTRAN CSHORE `.infile` into a `CshoreConfig`. Uses a line-oriented
token stream that mirrors the FORTRAN READ statement order; FORMAT specifiers
are irregular so numeric tokens are split on whitespace and comment/header
lines are captured verbatim.
==============================================================================#

# ----------------------------------------------------------------------------
# InfileTokens — helper for reading FORTRAN free-format numeric tokens.
#
# The FORTRAN CSHORE `.infile` has two distinct regions:
# 1. A HEADER block (verbatim comment lines read with A70 format) — handled
#    separately by `_read_header`.
# 2. A NUMERIC BODY (read via `READ(11, *)`) where all tokens are space- or
#    newline-separated. Inline annotations like `  1  ->ILINE` are ignored by
#    FORTRAN's free-format reader.
#
# `_tokenize_infile` reads the whole file, splits the header out, and
# tokenizes the rest while stripping inline annotations (anything after `->`
# on a line is discarded). The parser then consumes the token stream via
# `_next_int` / `_next_float` / `_next_floats`.
# ----------------------------------------------------------------------------

mutable struct InfileTokens
    header::Vector{String}   # verbatim comment lines from the header
    tokens::Vector{String}   # remaining numeric tokens (annotation-stripped)
    pos::Int                 # next token to consume (1-based)
    erdc::Bool               # true if file uses `-->` ERDC-style annotations
                             # (implies all bathymetry rows have 3 columns)
end

"""
    _tokenize_infile(path) -> InfileTokens

Read `path`, separate the header block (NLINES + that many verbatim lines),
then tokenize the rest while stripping inline `->...` annotations.

The first numeric token on the first line is interpreted as `NLINES`; the
next `NLINES` lines are stored verbatim in `tokens.header`; everything after
that is broken into whitespace-separated tokens, with any token containing
`->` and everything following it on the same line dropped.
"""
function _tokenize_infile(path::AbstractString)
    raw_lines = readlines(path)
    isempty(raw_lines) && throw(ArgumentError("empty infile: $path"))

    # Line 1: number of header/comment lines (FORTRAN reads it with free format
    # so trailing tokens are OK).
    first_tok = _first_numeric_token(raw_lines[1])
    first_tok === nothing && throw(ArgumentError(
        "infile $path: expected NLINES on line 1, got: $(raw_lines[1])"))
    nlines = parse(Int, first_tok)
    nlines ≥ 0 || throw(ArgumentError("infile $path: NLINES must be non-negative, got $nlines"))

    # Lines 2..1+nlines: verbatim comment text
    1 + nlines ≤ length(raw_lines) || throw(ArgumentError(
        "infile $path: declared $nlines header lines but only $(length(raw_lines)-1) remain"))
    header = String[rstrip(raw_lines[1+i]) for i in 1:nlines]

    # Remaining lines → tokenized body
    # Also detect ERDC `-->` style (double-dash arrow): those files write
    # friction on the first bathymetry row too (3-column format throughout).
    body_lines = @view raw_lines[(2+nlines):end]
    tokens = String[]
    erdc = false
    for line in body_lines
        erdc = erdc || occursin("-->", line)
        cut = _annotation_start(line)
        trimmed = cut === nothing ? line : line[1:cut-1]
        for t in split(trimmed)
            isempty(t) && continue
            push!(tokens, String(t))
        end
        push!(tokens, "\n")   # line sentinel — preserves record boundaries
    end
    return InfileTokens(header, tokens, 1, erdc)
end

"""
Return the byte index of the start of an inline annotation (first `-` in `->`
or `-->` sequences), or `nothing` if the line has no annotation.
Handles both the CSHORE_USACE `->VAR` and ERDC CSHORE `-->VAR` styles.
"""
function _annotation_start(line::AbstractString)
    # Find `->` anywhere on the line; walk back over any leading dashes so
    # that `-->` is trimmed including the first `-`.
    m = findfirst("->", line)
    m === nothing && return nothing
    pos = first(m)           # index of the first `-` in `->`
    while pos > 1 && line[pos-1] == '-'
        pos -= 1
    end
    return pos
end

"""Return the first whitespace-separated token on `line`, or `nothing` if empty."""
function _first_numeric_token(line::AbstractString)
    cut = _annotation_start(line)
    trimmed = cut === nothing ? line : line[1:cut-1]
    for t in split(trimmed)
        isempty(t) || return String(t)
    end
    return nothing
end

function _next_token(t::InfileTokens)
    # Skip line sentinels ("\n") — normal reads are unaffected by record boundaries.
    while t.pos ≤ length(t.tokens) && t.tokens[t.pos] == "\n"
        t.pos += 1
    end
    t.pos ≤ length(t.tokens) || throw(ArgumentError(
        "infile parser: unexpected end of tokens at position $(t.pos)"))
    s = t.tokens[t.pos]
    t.pos += 1
    return s
end

"""
    _skip_rest_of_line(t)

Discard remaining tokens on the current line and advance past the line sentinel.
Mirrors FORTRAN sequential I/O record-advancing behavior: after a READ is
satisfied, leftover values on the current record are discarded.
"""
function _skip_rest_of_line(t::InfileTokens)
    while t.pos ≤ length(t.tokens)
        tok = t.tokens[t.pos]
        t.pos += 1
        tok == "\n" && return
    end
end

_next_int(t::InfileTokens)   = parse(Int,   _next_token(t))
_next_float(t::InfileTokens) = parse(Float64, _next_token(t))
_next_floats(t::InfileTokens, n::Int) = [_next_float(t) for _ in 1:n]

"""
Peek the next int-valued token without consuming it. Returns `nothing` if the
stream is exhausted or the next token does not parse as an integer (allowing
a `.0` fractional part). Used by optional fields whose presence depends on
the file's format flavor (e.g. ISWLSL after IWTRAN=1 in older infiles).
"""
function _peek_int(t::InfileTokens)
    p = t.pos
    while p ≤ length(t.tokens) && t.tokens[p] == "\n"
        p += 1
    end
    p > length(t.tokens) && return nothing
    s = t.tokens[p]
    val = tryparse(Int, s)
    if val === nothing
        f = tryparse(Float64, s)
        f === nothing && return nothing
        abs(f - round(f)) > 1e-9 && return nothing
        val = Int(round(f))
    end
    return val
end

"""
Count the number of numeric tokens remaining on the current line (i.e. up to
the next `"\\n"` sentinel) without consuming any tokens. Used to detect
variable-column rows like the bathymetry block where the first row may have
2, 3, or more columns depending on the writer.
"""
function _tokens_until_newline(t::InfileTokens)
    p = t.pos
    while p ≤ length(t.tokens) && t.tokens[p] == "\n"
        p += 1
    end
    n = 0
    while p ≤ length(t.tokens) && t.tokens[p] != "\n"
        n += 1
        p += 1
    end
    return n
end

"""
    _next_iprofl(t) -> (iprofl::Int, ismooth::Int)

Decode the ISMOOTH-in-IPROFL encoding convention. IPROFL is read as a float;
`1.1` encodes `iprofl=1, ismooth=0`; `1.0` encodes `iprofl=1, ismooth=1`.
"""
function _next_iprofl(t::InfileTokens)
    raw = _next_float(t)
    iprofl_int = floor(Int, raw + 1e-9)
    frac = raw - iprofl_int
    # tmp = nint(10 * frac); if tmp == 1, ismooth = 0; else ismooth = 1
    tmp = round(Int, 10 * frac)
    ismooth = tmp == 1 ? 0 : 1
    return iprofl_int, ismooth
end

"""
    read_infile(path::AbstractString; strict=true, extensions=nothing) -> CshoreConfig

Parse a FORTRAN CSHORE USACE `.infile` and return a validated `CshoreConfig`.

All reads use free-format (whitespace-separated) semantics matching the
FORTRAN `READ(11, *)` read order.

The parser is **tolerant** of common writer conventions:
- Inline annotations like `  1   ->ILINE` are stripped
- Extra trailing numeric tokens on a line flow into the next read
- Header lines are captured verbatim

## Keyword arguments

- `strict=true` — when `false`, flags that are out-of-scope in CSHORE.jl
  (e.g. `IPOND=1`, `IWTRAN=1`, `ISEDAV=2`, `IWEIBULL=1`) are zeroed with a
  `@warn` instead of throwing. Use this when parsing ERDC CSHORE infiles from
  the field that may enable features not supported here. Supported flags
  (IOVER, IPERM, IVEG, IROLL, IWCINT, ITIDE, IWIND) always work regardless.
- `extensions` — an optional `Dict{Symbol,Any}` that injects CSHORE.jl-native
  configuration not expressible in the FORTRAN .infile format:
  ```julia
  cfg = read_infile("run.infile"; strict=false,
        extensions=Dict(
            :multifraction => MultifractionConfig(
                grain_sizes=[0.15e-3, 0.30e-3, 0.60e-3],
                initial_fractions=[0.3, 0.5, 0.2],
                transport_formula=:size_adaptive),
            :aeolian  => AeolianConfig(vegetation=ContourVegetation(z_contour=3.0)),
            :thermal  => ThermalConfig(T_init=-5.0),
            :diffusion => DiffusionConfig(),
        ))
  ```
  Recognised keys: `:multifraction`, `:aeolian`, `:thermal`, `:diffusion`,
  `:tidal`, `:current`, `:vegetation`, `:porous`.

# Example

```julia
cfg = read_infile("test/fixtures/simple.infile")
state = run_simulation!(cfg; outfile="out.nc")

# Parse a field infile that uses unsupported flags — zero them, continue:
cfg = read_infile("field_run.infile"; strict=false)
```
"""
function read_infile(path::AbstractString; strict::Bool=true,
                     extensions::Union{Nothing,Dict}=nothing)
    isfile(path) || throw(ArgumentError("read_infile: file not found: $path"))
    t = _tokenize_infile(path)

    # -------- Flags --------
    iline = _next_int(t)
    iqydy = iline > 2 ? _next_int(t) : 0

    iprofl, ismooth = _next_iprofl(t)
    isedav = iprofl == 1 ? _next_int(t) : 0

    iperm = _next_int(t)
    iover = _next_int(t)
    iwtran = iover == 1 ? _next_int(t) : 0
    ipond  = (iover == 1 && iwtran == 0) ? _next_int(t) : 0
    infilt = (iover == 1 && iperm == 0 && iprofl == 1) ? _next_int(t) : 0

    iwcint = _next_int(t)
    iroll  = _next_int(t)
    iwind  = _next_int(t)
    itide  = _next_int(t)
    iveg   = _next_int(t)

    idiss = iveg == 3 ? _next_int(t) : 0
    ifv   = iveg == 3 ? _next_int(t) : 0

    iclay = (abs(isedav) == 1 && iperm == 0 && iveg == 0) ? _next_int(t) : 0

    # -------- Physical params (DX, GAMMA) --------
    dx    = _next_float(t)
    gamma = _next_float(t)

    # -------- Sediment block (if IPROFL==1) --------
    # Defaults for fixed-bed mode (not used by physics when iprofl=0 but
    # required for CshoreConfig construction).
    d50_mm = 0.3; wf = 0.0381; sg = 2.65
    effb = 0.005; efff = 0.01; slp = 0.2; slpot = 0.1
    tanphi = 0.63; blp = 0.002
    rwh = 0.0
    if iprofl == 1
        d50_mm = _next_float(t)
        wf     = _next_float(t)
        sg     = _next_float(t)
        effb = _next_float(t)
        efff = _next_float(t)
        slp  = _next_float(t)
        if iover == 1
            slpot = _next_float(t)
        end
        tanphi = _next_float(t)
        blp    = _next_float(t)
        if iperm == 1
            _next_float(t); _next_float(t); _next_float(t)  # SNP SDP CSTABN (skipped)
        end
        if iprofl == 2 || iclay == 1
            _next_float(t); _next_float(t)  # DEEB DEEF (skipped)
        end
    end
    # RWH is read whenever IOVER=1, independent of IPROFL
    if iover == 1
        rwh = _next_float(t)
    end

    # -------- Wave / SWL time series --------
    ilab  = _next_int(t)
    nwave = _next_int(t)
    nsurg = _next_int(t)

    local timebc::Vector{Float64}
    local tpbc::Vector{Float64}
    local hrmsbc::Vector{Float64}
    local wangbc::Vector{Float64}
    local wsetbc::Vector{Float64}
    local swlbc::Vector{Float64}

    if ilab == 1
        # Laboratory mode: NWAVE rows, each 6 cols (t, Tp, Hrms, Wsetup, SWL, angle)
        n = nwave
        timebc = Vector{Float64}(undef, n)
        tpbc   = Vector{Float64}(undef, n)
        hrmsbc = Vector{Float64}(undef, n)
        wsetbc = Vector{Float64}(undef, n)
        swlbc  = Vector{Float64}(undef, n)
        wangbc = Vector{Float64}(undef, n)
        for i in 1:n
            timebc[i] = _next_float(t)
            tpbc[i]   = _next_float(t)
            hrmsbc[i] = _next_float(t)
            wsetbc[i] = _next_float(t)
            swlbc[i]  = _next_float(t)
            wangbc[i] = _next_float(t)
        end
    else
        # Field mode: (NWAVE+1) rows of (t, Tp, Hrms, angle), then (NSURG+1) rows of (t, SWL)
        n = nwave + 1
        timebc = Vector{Float64}(undef, n)
        tpbc   = Vector{Float64}(undef, n)
        hrmsbc = Vector{Float64}(undef, n)
        wangbc = Vector{Float64}(undef, n)
        for i in 1:n
            timebc[i] = _next_float(t)
            tpbc[i]   = _next_float(t)
            hrmsbc[i] = _next_float(t)
            wangbc[i] = _next_float(t)
        end
        wsetbc = zeros(n)
        # Surge series on its own time base
        nsurg_total = nsurg + 1
        tsurg = Vector{Float64}(undef, nsurg_total)
        swlraw = Vector{Float64}(undef, nsurg_total)
        for i in 1:nsurg_total
            tsurg[i]  = _next_float(t)
            swlraw[i] = _next_float(t)
        end
        # Interpolate surge onto the wave time base so we have a single
        # timebc array (CshoreConfig uses one master time series).
        swlbc = Vector{Float64}(undef, n)
        for i in 1:n
            swlbc[i] = interp1(tsurg, swlraw, timebc[i])
        end
    end

    # -------- Wind (if IWIND==1) --------
    w10_arr     = Float64[]
    wangle_arr  = Float64[]
    windcd_arr  = Float64[]
    if iwind == 1
        nwind = _next_int(t)
        nwind_total = nwind + 1
        w10_arr    = Vector{Float64}(undef, nwind_total)
        wangle_arr = Vector{Float64}(undef, nwind_total)
        for i in 1:nwind_total
            _next_float(t)               # twind[i] (not stored in BoundaryTimeSeries)
            w10_arr[i]    = _next_float(t)
            wangle_arr[i] = _next_float(t)
        end
        windcd_arr = fill(0.0015, nwind_total)   # default drag coefficient
    end

    # -------- Landward SWL (IWTRAN==1) — consume to keep pointer aligned --------
    # ISWLSL is required by the canonical CSHORE_USACE format but absent from
    # some older / private infiles (e.g. the Vidal & Mansard 1995 structure
    # test). Peek and only consume if the next token looks like a valid
    # ISWLSL value (0, 1, or 2). Otherwise default to 0 (seaward SWL = landward
    # SWL) and leave the token for the bathymetry section.
    swl_landward_pairs = Tuple{Float64,Float64}[]
    if iwtran == 1
        iswlsl_peek = _peek_int(t)
        if iswlsl_peek !== nothing && iswlsl_peek in (0, 1, 2)
            iswlsl = _next_int(t)
            if iswlsl == 1
                nslan = _next_int(t)
                for _ in 1:(nslan + 1)
                    ts  = _next_float(t)
                    swl = _next_float(t)
                    push!(swl_landward_pairs, (ts, swl))
                end
            end
        end
    end

    # -------- Tidal gradient (ITIDE==1) --------
    if itide == 1
        ntide = _next_int(t)
        for _ in 1:(ntide + 1)
            _next_float(t); _next_float(t)
        end
    end

    # -------- Bathymetry loop (L = 1..ILINE) --------
    # First pass: read NBINP for each line so we can allocate (max_nbinp, iline)
    nbinp_vec = zeros(Int, iline)
    npinp_vec = zeros(Int, iline)
    yline_vec = zeros(Float64, iline)
    agline_vec = zeros(Float64, iline)

    # We do a single pass but need variable-length per line → read into
    # a vector of vectors, then pack into a matrix at the end.
    xbinp_cols = [Float64[] for _ in 1:iline]
    zbinp_cols = [Float64[] for _ in 1:iline]
    fbinp_cols = [Float64[] for _ in 1:iline]
    xpinp_cols = [Float64[] for _ in 1:iline]
    zpinp_cols = [Float64[] for _ in 1:iline]

    # Vegetation raw arrays (per line); populated when iveg >= 1
    vegcd_val  = 1.0
    vegcdm_val = 1.0
    vegn_cols  = [Float64[] for _ in 1:iline]  # stem density
    vegb_cols  = [Float64[] for _ in 1:iline]  # blade width
    vegd_cols  = [Float64[] for _ in 1:iline]  # stem height
    vegh_cols  = [Float64[] for _ in 1:iline]  # canopy height

    for l in 1:iline
        if iline > 1
            yline_vec[l]  = _next_float(t)
            agline_vec[l] = _next_float(t)
        end
        nb = _next_int(t)
        nbinp_vec[l] = nb
        if iperm == 1 || abs(isedav) ≥ 1
            npinp_vec[l] = _next_int(t)
        end

        # Bathymetry rows.
        # Column count on row 1 varies by writer:
        #   2 cols (x, z)             — private CSHORE_USACE format and the
        #                                Vidal & Mansard 1995 structure test
        #   3 cols (x, z, fw)         — most ERDC writers; fw on row 1 also
        #   ≥4 cols                   — wire-mesh / ISEDAV>1 with extra columns
        # FORTRAN READ(11,*) is record-based and discards leftover tokens on
        # the row, so we read x,z (and fw if present) and skip the rest.
        xbinp_cols[l] = Vector{Float64}(undef, nb)
        zbinp_cols[l] = Vector{Float64}(undef, nb)
        fbinp_cols[l] = Vector{Float64}(undef, nb)
        ncols_row1 = _tokens_until_newline(t)
        xbinp_cols[l][1] = _next_float(t)
        zbinp_cols[l][1] = _next_float(t)
        if ncols_row1 ≥ 3
            fbinp_cols[l][1] = _next_float(t)
        end
        _skip_rest_of_line(t)

        for j in 2:nb
            xbinp_cols[l][j]   = _next_float(t)
            zbinp_cols[l][j]   = _next_float(t)
            fbinp_cols[l][j-1] = _next_float(t)
            if abs(isedav) > 1
                _next_float(t)   # WMINP (wire mesh flag) — unused
            end
        end
        fbinp_cols[l][nb] = fbinp_cols[l][nb-1]   # pad last element for matrix packing

        # Porous / hardbottom layer (if applicable)
        if iperm == 1 || abs(isedav) ≥ 1
            np_l = npinp_vec[l]
            xpinp_cols[l] = Vector{Float64}(undef, np_l)
            zpinp_cols[l] = Vector{Float64}(undef, np_l)
            # First point is implicit: x=0, z=zbinp[1]
            xpinp_cols[l][1] = 0.0
            zpinp_cols[l][1] = zbinp_cols[l][1]
            for j in 2:np_l
                xpinp_cols[l][j] = _next_float(t)
                zpinp_cols[l][j] = _next_float(t)
                if iclay == 1
                    _next_float(t); _next_float(t)   # RCINP, FCINP — out of scope
                end
            end
        end

        # Vegetation (if IVEG>=1) — read into per-line raw arrays.
        # FORTRAN reads (nb-1) segment values; we assign each to its seaward
        # node and repeat the last value for the landward endpoint.
        if iveg >= 1
            vegcd_val  = _next_float(t)  # VEGCD
            if iveg == 3
                vegcdm_val = _next_float(t)  # VEGCDM
            end
            has4 = (iveg == 1 || iveg == 3)   # 4-column (SPVEG BVEG DVEG HVEG) vs 3-column
            vegn_l = Vector{Float64}(undef, nb)
            vegb_l = Vector{Float64}(undef, nb)
            vegd_l = Vector{Float64}(undef, nb)
            vegh_l = Vector{Float64}(undef, nb)
            for j in 1:(nb-1)
                vegn_l[j] = _next_float(t)
                vegb_l[j] = _next_float(t)
                vegd_l[j] = _next_float(t)
                vegh_l[j] = has4 ? _next_float(t) : vegd_l[j]
            end
            # Repeat last segment value for the landward endpoint
            vegn_l[nb] = vegn_l[nb-1]
            vegb_l[nb] = vegb_l[nb-1]
            vegd_l[nb] = vegd_l[nb-1]
            vegh_l[nb] = vegh_l[nb-1]
            vegn_cols[l] = vegn_l
            vegb_cols[l] = vegb_l
            vegd_cols[l] = vegd_l
            vegh_cols[l] = vegh_l
        end

        # Dike (IPROFL==2) — consume tokens
        if iprofl == 2
            for _ in 1:(nb - 1), _c in 1:3
                _next_float(t)
            end
        end
    end

    # -------- Pack bathymetry into matrices --------
    max_nb = maximum(nbinp_vec)
    xbinp = zeros(Float64, max_nb, iline)
    zbinp = zeros(Float64, max_nb, iline)
    fbinp = zeros(Float64, max_nb, iline)
    for l in 1:iline
        nb = nbinp_vec[l]
        xbinp[1:nb, l] = xbinp_cols[l]
        zbinp[1:nb, l] = zbinp_cols[l]
        fbinp[1:nb, l] = fbinp_cols[l]
    end

    # -------- Build VegetationInput (if IVEG >= 1) --------
    # vegfb = vegcd * vegn * vegb / fb2 must be precomputed here because
    # read_infile constructs CshoreConfig directly (bypassing build_config,
    # which is the other place that fills vegfb). Without this, the wave
    # solver hits a 0×0 BoundsError on the first vegetated step.
    vegetation_in = if iveg >= 1
        vegn_mat = zeros(Float64, max_nb, iline)
        vegb_mat = zeros(Float64, max_nb, iline)
        vegd_mat = zeros(Float64, max_nb, iline)
        vegh_mat = zeros(Float64, max_nb, iline)
        vegfb_mat = zeros(Float64, max_nb, iline)
        for l in 1:iline
            nb = nbinp_vec[l]
            vegn_mat[1:nb, l] = vegn_cols[l]
            vegb_mat[1:nb, l] = vegb_cols[l]
            vegd_mat[1:nb, l] = vegd_cols[l]
            vegh_mat[1:nb, l] = vegh_cols[l]
            for j in 1:nb
                fb2 = max(fbinp[j, l], 1e-12)
                vegfb_mat[j, l] = vegcd_val * vegn_mat[j, l] * vegb_mat[j, l] / fb2
            end
        end
        VegetationInput(; vegcd=vegcd_val, vegcdm=vegcdm_val,
                          vegn=vegn_mat, vegb=vegb_mat,
                          vegd=vegd_mat, vegh=vegh_mat,
                          vegfb=vegfb_mat)
    else
        nothing
    end

    bathy = BathyInput(
        xbinp=xbinp, zbinp=zbinp, fbinp=fbinp,
        nbinp=nbinp_vec,
        xs=[xbinp[1, l] for l in 1:iline],
        yline=iline > 1 ? yline_vec : Float64[],
        dyline=Float64[],
        agline=iline > 1 ? agline_vec : zeros(iline),
    )

    # -------- Tolerant-mode flag fixup --------
    # When strict=false, zero out-of-scope flags with warnings instead of
    # letting validate() throw. Flags that CSHORE.jl now supports are left alone.
    if !strict
        # IWTRAN=1 is now supported — leave it alone.
        if ipond == 1
            @warn "read_infile: IPOND=1 (ridge-runnel ponded water) is not supported — setting IPOND=0"
            ipond = 0
        end
        if isedav == 2
            @warn "read_infile: ISEDAV=2 (wire mesh) is not supported — setting ISEDAV=0 (unlimited sand)"
            isedav = 0
        end
    end

    # -------- Build CshoreConfig --------
    opt = OptionFlags(
        iprofl=iprofl, iangle=0, iroll=iroll, iwind=iwind, iperm=iperm,
        iover=iover, iwcint=iwcint, isedav=isedav, iwtran=iwtran,
        ivwall=Int[], ilab=ilab, infilt=infilt, ipond=ipond, itide=itide,
        iline=iline, iqydy=iqydy, iveg=iveg, iclay=iclay, ismooth=ismooth,
        idiss=idiss, ifv=ifv, iweibull=0,
    )
    sed = make_sediment(
        d50=d50_mm * 1e-3,   # infile stores mm; convert to m
        sg=sg,
        # Legacy FORTRAN infile format doesn't carry rho_water. Pin to 1000.0
        # (freshwater) so legacy infiles reproduce FORTRAN CSHORE numerics
        # exactly. Users wanting the seawater correction should migrate to the
        # native CSHORE.jl TOML format (read_cshorejl) or pass rho_water via
        # build_config directly.
        rho_water=1000.0,
        sporo=0.4,   # not read from infile; use default
        shield=0.05,
        blp=blp,
        tanphi=tanphi,
        effb=effb,
        efff=efff,
        slp=slp,
        slpot=slpot,
    )
    # Interpolate landward SWL series onto the wave timebc grid (IWTRAN=1)
    swl_landward_vec = if isempty(swl_landward_pairs)
        Float64[]
    else
        ts_land = [p[1] for p in swl_landward_pairs]
        swl_land = [p[2] for p in swl_landward_pairs]
        [interp1(ts_land, swl_land, t) for t in timebc]
    end
    bc = BoundaryTimeSeries(
        timebc=timebc, tpbc=tpbc, hrmsbc=hrmsbc,
        wsetbc=wsetbc, swlbc=swlbc, wangbc=wangbc,
        w10=w10_arr, wangle=wangle_arr, windcd=windcd_arr,
        swl_landward=swl_landward_vec,
    )
    grid = GridConfig(dx=dx, nn=max(20000, 2 * max_nb + 100), nl=iline)

    # Apply CSHORE.jl-native extensions from the `extensions` dict.
    # These cover capabilities the FORTRAN .infile cannot express.
    ext = extensions === nothing ? Dict() : extensions
    multifraction = get(ext, :multifraction, MultifractionConfig())
    aeolian       = get(ext, :aeolian,       nothing)
    thermal       = get(ext, :thermal,       nothing)
    diffusion     = get(ext, :diffusion,     nothing)
    tidal_ext     = get(ext, :tidal,         nothing)
    current_ext   = get(ext, :current,       nothing)
    vegetation_ext = get(ext, :vegetation,   nothing)
    porous_ext    = get(ext, :porous,        nothing)

    # Auto-set flags when extensions inject active objects
    if aeolian !== nothing && opt.iaeolian == 0
        opt = OptionFlags(; (f => getfield(opt, f) for f in fieldnames(OptionFlags) if f != :iaeolian)...,
                            iaeolian=1)
    end
    if tidal_ext !== nothing && opt.itide == 0
        opt = OptionFlags(; (f => getfield(opt, f) for f in fieldnames(OptionFlags) if f != :itide)...,
                            itide=1)
    end
    if current_ext !== nothing && opt.icurrent == 0
        opt = OptionFlags(; (f => getfield(opt, f) for f in fieldnames(OptionFlags) if f != :icurrent)...,
                            icurrent=1)
    end

    cfg = CshoreConfig(;
        options=opt,
        grid=grid,
        sediment=sed,
        multifraction=multifraction,
        boundary=bc,
        bathymetry=bathy,
        gamma=gamma,
        aeolian=aeolian,
        thermal=thermal,
        diffusion=diffusion,
        tidal=tidal_ext !== nothing ? tidal_ext : nothing,
        current=current_ext !== nothing ? current_ext : nothing,
        # Merge: extension-supplied vegetation takes precedence; fall back
        # to what was parsed from the infile (vegetation_in).
        vegetation=vegetation_ext !== nothing ? vegetation_ext : vegetation_in,
        porous=porous_ext,
    )
    validate(cfg)
    return cfg
end

"""
    build_config(; kwargs...) -> CshoreConfig

Convenience constructor for programmatic test inputs. Useful until
`read_infile` is ported. Wraps keyword arguments into the nested struct
hierarchy.
"""
function build_config(; dx::Float64, bathymetry_x::Vector{Float64},
                       bathymetry_z::Vector{Float64},
                       friction::Union{Float64, Vector{Float64}}=0.002,
                       manning::Union{Nothing, Float64, Vector{Float64}, VegetatedManningField}=nothing,
                       timebc::Vector{Float64},
                       hrmsbc::Vector{Float64},
                       tpbc::Vector{Float64},
                       swlbc::Vector{Float64},
                       wangbc::Vector{Float64}=fill(0.0, length(timebc)),
                       wsetbc::Vector{Float64}=fill(0.0, length(timebc)),
                       options::OptionFlags=OptionFlags(),
                       sediment::SedimentConfig=SedimentConfig(),
                       multifraction::MultifractionConfig=MultifractionConfig(),
                       hardbottom_z::Union{Nothing,Vector{Float64}}=nothing,
                       porous::Union{Nothing,PorousInput}=nothing,
                       porous_z::Union{Nothing,Vector{Float64}}=nothing,
                       porosity::Float64=0.4,
                       stone_diameter::Float64=0.02,
                       vegetation::Union{Nothing,VegetationInput}=nothing,
                       clay::Union{Nothing,ClayInput}=nothing,
                       dike::Union{Nothing,DikeErosionInput}=nothing,
                       tidal::Union{Nothing,TidalInput}=nothing,
                       current::Union{Nothing,CurrentInput}=nothing,
                       aeolian::Union{Nothing,AeolianConfig}=nothing,
                       windshear::Union{Nothing,WindShearConfig}=nothing,
                       ig::Union{Nothing,IgConfig}=nothing,
                       undertow::Union{Nothing,UndertowConfig}=nothing,
                       asymm::Union{Nothing,AsymmetryConfig}=nothing,
                       phase_lag::Union{Nothing,PhaseLagConfig}=nothing,
                       bailard::Union{Nothing,BailardConfig}=nothing,
                       wave_nonlinearity::Union{Nothing,WaveNonlinearityConfig}=nothing,
                       groundwater::Union{Nothing,GroundwaterConfig}=nothing,
                       rcrest::Float64=NaN,
                       diffusion::Union{Nothing,DiffusionConfig}=nothing,
                       cohesive::Union{Nothing,CohesiveSedimentConfig}=nothing,
                       thermal::Union{Nothing,ThermalConfig}=nothing,
                       snow::Union{Nothing,SnowConfig}=nothing,
                       snow_modifier::Union{Nothing,SnowSpatialModifier}=nothing,
                       thermal_time::Union{Nothing,Vector{Float64}}=nothing,
                       T_air::Union{Nothing,Vector{Float64}}=nothing,
                       T_water::Union{Nothing,Vector{Float64}}=nothing,
                       snow_depth::Union{Nothing,Vector{Float64}}=nothing,
                       w10::Union{Nothing,Vector{Float64}}=nothing,
                       wangle::Union{Nothing,Vector{Float64}}=nothing,
                       windcd::Union{Nothing,Vector{Float64}}=nothing,
                       gamma::Float64=0.78,
                       gamma_method::Symbol=:constant,
                       gamma_a::Float64=0.76,
                       gamma_b::Float64=0.29,
                       gamma_min::Float64=0.35,
                       gamma_max::Float64=0.90,
                       gamma_sr_slope::Float64=3.9,
                       morph_courant::Float64=0.3,
                       max_dzb_per_step::Float64=0.1,
                       morph_diffusion::Float64=0.0,
                       breaker_delay::Float64=0.0,
                       min_depth_wcint::Float64=0.10,
                       facSK::Float64=1.0,
                       facAS::Float64=0.0,
                       ur_sk_ref::Float64=0.20,
                       alpha_sk::Float64=0.0,
                       biphase_relax_L::Float64=0.0,
                       f_min::Float64=0.002,
                       f_sheet::Float64=0.015,
                       theta_sheet::Float64=1.0,
                       f_ripple_exp::Float64=0.5)
    np = length(bathymetry_x)
    np == length(bathymetry_z) || throw(DimensionMismatch("bathymetry_x and bathymetry_z must match"))

    # ── Wave-nonlinearity deprecation warnings (silent when wave_nonlinearity
    # is the source of truth or all legacy fields are at their defaults) ──
    if wave_nonlinearity === nothing
        legacy_used = (facSK != 1.0) || (facAS != 0.0) || (ur_sk_ref != 0.20) ||
                      (alpha_sk != 0.0) || (biphase_relax_L != 0.0) ||
                      (asymm !== nothing) || (phase_lag !== nothing) ||
                      (bailard !== nothing) ||
                      (options.iasym != 0) || (options.iskew_spatial != 0)
        if legacy_used
            Base.depwarn(
                "Legacy wave-nonlinearity kwargs (facSK / facAS / ur_sk_ref / " *
                "alpha_sk / biphase_relax_L / asymm / phase_lag / bailard) and " *
                "OptionFlags.iasym / iskew_spatial are deprecated. Pass " *
                "`wave_nonlinearity = WaveNonlinearityConfig(...)` to build_config " *
                "instead — see the WaveNonlinearityConfig docstring for fields.",
                :build_config)
        end
    end

    # Convenience: auto-build PorousInput from simple vector kwargs
    if porous === nothing && porous_z !== nothing
        length(porous_z) == np ||
            throw(DimensionMismatch("porous_z must match length of bathymetry_x"))
        porous = PorousInput(bathymetry_x, porous_z;
                             porosity=porosity, stone_diameter=stone_diameter)
    end

    nl = options.iline
    # Pack into the (max_nbinp, iline) matrices — for nl=1 just a single column.
    xbinp = reshape(copy(bathymetry_x), np, 1)
    zbinp = reshape(copy(bathymetry_z), np, 1)
    # Friction: scalar (uniform) or vector (per bathymetry node)
    fbinp = if friction isa Float64
        fill(friction, np, 1)
    else
        length(friction) == np ||
            throw(DimensionMismatch("friction vector must match length of bathymetry_x ($np)"))
        reshape(copy(friction), np, 1)
    end
    # Hardbottom: when iperm=1 and no explicit hardbottom_z is given,
    # auto-generate from the porous layer bottom (zpinp). The porous layer
    # thickness hp = zb - zb_hard, so setting zb_hard = zp gives the
    # correct porous depth. Also auto-set isedav=1 for hardbottom enforcement.
    effective_hb = hardbottom_z
    if effective_hb === nothing && porous !== nothing && options.iperm == 1
        # Interpolate zpinp onto the bathymetry grid
        np_por = size(porous.zpinp, 1)
        xp_raw = view(porous.xpinp, :, 1)
        zp_raw = view(porous.zpinp, :, 1)
        effective_hb = [interp1(xp_raw, zp_raw, bathymetry_x[i]) for i in 1:np]
    end
    zhinp = if effective_hb === nothing
        zeros(0, 0)
    else
        length(effective_hb) == np ||
            throw(DimensionMismatch("hardbottom_z must match length of bathymetry_x ($np)"))
        reshape(copy(effective_hb), np, 1)
    end
    # Manning's n: scalar, per-node, or vegetated field. When provided, sets
    # ifriction_spatial=2 so the driver recomputes fb2 from live water depth.
    #
    # Double-counting guard: a VegetatedManningField already encodes the
    # vegetation's frictional effect; combining it with a VegetationInput
    # (which adds Cd-based stem drag or fb2 multiplier) would count the
    # vegetation twice. Reject up front.
    if manning isa VegetatedManningField && vegetation !== nothing
        throw(ArgumentError(
            "Cannot combine a VegetatedManningField (categories=$(manning.categories)) " *
            "with a VegetationInput. The Manning field already encodes vegetation drag — " *
            "adding IVEG on top double-counts. Either (a) drop vegetation= and use the " *
            "Manning-only path, or (b) drop the VegetatedManningField and pass a plain " *
            "Vector{Float64} of bed-only Manning's n (e.g. bare-sand n ≈ 0.02) plus the " *
            "VegetationInput for Cd-based stem drag."))
    end
    if manning isa VegetatedManningField && options.iveg != 0
        throw(ArgumentError(
            "Cannot combine a VegetatedManningField with options.iveg=$(options.iveg). " *
            "Use iveg=0 with the Manning-only path, or drop the VegetatedManningField."))
    end
    manning_values = manning isa VegetatedManningField ? manning.values : manning
    manning_n_mat = if manning_values === nothing
        zeros(Float64, 0, 0)
    elseif manning_values isa Float64
        fill(manning_values, np, 1)
    else
        length(manning_values) == np ||
            throw(DimensionMismatch("manning vector must match length of bathymetry_x ($np)"))
        reshape(copy(manning_values), np, 1)
    end
    if snow_modifier !== nothing
        length(snow_modifier) == np || throw(DimensionMismatch(
            "snow_modifier length ($(length(snow_modifier))) must match bathymetry_x ($np)"))
        snow !== nothing || throw(ArgumentError(
            "snow_modifier requires a SnowConfig (pass snow=SnowConfig() with the modifier)"))
    end
    if manning_values !== nothing && options.ifriction_spatial == 0
        options = OptionFlags(; (k => getfield(options, k) for k in fieldnames(OptionFlags))...,
                                ifriction_spatial=2)
    end
    bathy = BathyInput(; xbinp=xbinp, zbinp=zbinp, fbinp=fbinp, zhinp=zhinp,
                         manning_n=manning_n_mat,
                         nbinp=[np], xs=[bathymetry_x[1]],
                         yline=Float64[], dyline=Float64[], agline=zeros(1))
    bc = BoundaryTimeSeries(; timebc=timebc, tpbc=tpbc, hrmsbc=hrmsbc,
                             wsetbc=wsetbc, swlbc=swlbc, wangbc=wangbc,
                             w10=w10 === nothing ? Float64[] : collect(w10),
                             wangle=wangle === nothing ? Float64[] : collect(wangle),
                             windcd=windcd === nothing ? Float64[] : collect(windcd))
    grid = GridConfig(; dx=dx, nn=max(20000, 2np + 100), nl=nl)

    # Auto-set flags from provided inputs
    opts = options
    # Auto-set iveg from vegetation presence if not already set
    if vegetation !== nothing && opts.iveg == 0
        opts = OptionFlags(; (f => getfield(opts, f) for f in fieldnames(OptionFlags) if f != :iveg)...,
                             iveg=1)
    end
    # Auto-set isedav=1 when iperm=1 and hardbottom is available (porous layer
    # bottom defines the non-erodible floor). Without isedav, hp stays at Inf
    # and the porous flow has no spatial extent information.
    if opts.iperm == 1 && opts.isedav == 0 && !isempty(zhinp)
        opts = OptionFlags(; (f => getfield(opts, f) for f in fieldnames(OptionFlags) if f != :isedav)...,
                             isedav=1)
    end

    # Optional thermal BCs: if any of thermal_time/T_air/T_water are
    # supplied, all three must be and must match. Falls back to the wave
    # BC time grid when `thermal_time` is omitted.
    thermal_bc = nothing
    if thermal !== nothing
        tt = thermal_time === nothing ? timebc : thermal_time
        Ta = T_air   === nothing ? fill(0.0, length(tt)) : T_air
        Tw = T_water === nothing ? fill(0.0, length(tt)) : T_water
        length(Ta) == length(tt) == length(Tw) ||
            throw(DimensionMismatch("thermal_time, T_air, T_water must all have the same length"))
        sd = snow_depth === nothing ? Float64[] : collect(snow_depth)
        if !isempty(sd)
            length(sd) == length(tt) ||
                throw(DimensionMismatch("snow_depth must match thermal_time length"))
        end
        thermal_bc = ThermalBoundaryTimeSeries(; time=collect(tt),
                                                 T_air=collect(Ta),
                                                 T_water=collect(Tw),
                                                 snow_depth=sd)
    end

    # Interpolate vegetation arrays onto the uniform grid and precompute vegfb.
    # The user supplies VegetationInput with raw-resolution arrays matching
    # bathymetry_x (np nodes). We interpolate each field to the uniform grid
    # (nn nodes) and compute vegfb = vegcd * vegn * vegb / fb2.
    veg_cfg = vegetation
    if vegetation !== nothing
        nn_grid = grid.nn
        x0 = bathymetry_x[1]
        nnodes = floor(Int, (bathymetry_x[end] - x0) / dx) + 1
        has_veg_arrays = size(vegetation.vegn, 1) > 0
        if has_veg_arrays
            # Interpolate each veg field from raw nodes to grid nodes
            vn = zeros(nn_grid, nl); vb = zeros(nn_grid, nl)
            vd = zeros(nn_grid, nl); vh = zeros(nn_grid, nl)
            vrd = zeros(nn_grid, nl); vrh = zeros(nn_grid, nl)
            vfb = zeros(nn_grid, nl)
            for l_i in 1:nl
                np_raw = size(vegetation.vegn, 1)
                # Raw arrays may be (np,1) or (np,nl) — use column l_i if available
                col = min(l_i, size(vegetation.vegn, 2))
                for j in 1:nnodes
                    xj = x0 + (j - 1) * dx
                    vn[j, l_i]  = interp1(bathymetry_x, view(vegetation.vegn, :, col), xj)
                    vb[j, l_i]  = interp1(bathymetry_x, view(vegetation.vegb, :, col), xj)
                    vd[j, l_i]  = interp1(bathymetry_x, view(vegetation.vegd, :, col), xj)
                    vh[j, l_i]  = interp1(bathymetry_x, view(vegetation.vegh, :, col), xj)
                    if size(vegetation.vegrd, 1) > 0
                        vrd[j, l_i] = interp1(bathymetry_x, view(vegetation.vegrd, :, col), xj)
                    end
                    if size(vegetation.vegrh, 1) > 0
                        vrh[j, l_i] = interp1(bathymetry_x, view(vegetation.vegrh, :, col), xj)
                    end
                    # Precompute vegfb = vegcd * vegn * vegb / fb2
                    fb2_j = if friction isa Float64
                        friction
                    else
                        interp1(bathymetry_x, friction, xj)
                    end
                    vfb[j, l_i] = vegetation.vegcd * vn[j, l_i] * vb[j, l_i] / max(fb2_j, 1e-12)
                end
            end
            veg_cfg = VegetationInput(; vegcd=vegetation.vegcd, vegcdm=vegetation.vegcdm,
                                        vegn=vn, vegb=vb, vegd=vd, vegh=vh,
                                        vegrd=vrd, vegrh=vrh, vegfb=vfb)
        end
    end

    # Swash config: set rcrest if provided
    swash_cfg = if !isnan(rcrest)
        SwashConfig(rcrest=[rcrest])
    else
        SwashConfig()
    end

    return CshoreConfig(; options=opts, grid=grid, sediment=sediment,
                         multifraction=multifraction, boundary=bc,
                         bathymetry=bathy, gamma=gamma,
                         gamma_method=gamma_method, gamma_a=gamma_a,
                         gamma_b=gamma_b, gamma_min=gamma_min, gamma_max=gamma_max,
                         gamma_sr_slope=gamma_sr_slope,
                         vegetation=veg_cfg,
                         porous=porous,
                         clay=clay,
                         dike=dike,
                         tidal=tidal,
                         current=current,
                         aeolian=aeolian,
                         windshear=windshear,
                         ig=ig,
                         undertow=undertow,
                         asymm=asymm,
                         phase_lag=phase_lag,
                         bailard=bailard,
                         wave_nonlinearity=wave_nonlinearity,
                         groundwater=groundwater,
                         diffusion=diffusion,
                         cohesive=cohesive,
                         swash=swash_cfg,
                         thermal=thermal, thermal_bc=thermal_bc,
                         snow=snow,
                         snow_modifier=snow_modifier,
                         morph_courant=morph_courant,
                         max_dzb_per_step=max_dzb_per_step,
                         morph_diffusion=morph_diffusion,
                         breaker_delay=breaker_delay,
                         min_depth_wcint=min_depth_wcint,
                         facSK=facSK, facAS=facAS,
                         ur_sk_ref=ur_sk_ref,
                         alpha_sk=alpha_sk,
                         biphase_relax_L=biphase_relax_L,
                         f_min=f_min, f_sheet=f_sheet,
                         theta_sheet=theta_sheet, f_ripple_exp=f_ripple_exp)
end
