#==============================================================================
CSHOREDataExt — HTTP / DataFrames-backed implementations of the CSHORE.jl
environmental-data-fetching stubs declared in `src/data.jl`.

Loaded automatically when the user has `HTTP`, `JSON3`, `CSV`, AND
`DataFrames` available in their environment (Julia 1.9+ package-extension
mechanism).

When any of those four packages is NOT loaded (e.g. inside the compiled
`cshore-julia` binary, where they're excluded from the build), the stubs
in `src/data.jl` raise a clear error pointing the user at the right
`using ...` line.

This split lets the bundled binary stay ~200 MB lighter (DataFrames + CSV
+ HTTP + JSON3 are heavy) and dramatically reduces sysimage compilation
peak memory for PackageCompiler builds on Windows / memory-constrained
machines.
==============================================================================#

module CSHOREDataExt

import CSHORE
import HTTP
import JSON3
import CSV
import DataFrames: DataFrame, rename!, select!, Not, nrow
import Dates
import Dates: DateTime, DateFormat, datetime2unix, unix2datetime, Second, Minute

# ---------------------------------------------------------------------------
# 1. fetch_ndbc_realtime
# ---------------------------------------------------------------------------

function CSHORE.fetch_ndbc_realtime(station_id::AbstractString; parameter::Symbol=:stdmet)
    url = "https://www.ndbc.noaa.gov/data/realtime2/$(station_id).$(parameter).txt"
    @info "Fetching NDBC realtime data" station_id url
    body = _http_get_text(url)
    return _parse_ndbc_text(body)
end

# ---------------------------------------------------------------------------
# 2. fetch_ndbc_historical
# ---------------------------------------------------------------------------

function CSHORE.fetch_ndbc_historical(station_id::AbstractString; year::Int)
    filename = "$(station_id)h$(year).txt.gz"
    url = "https://www.ndbc.noaa.gov/view_text_file.php?filename=$(filename)&dir=data/historical/stdmet/"
    @info "Fetching NDBC historical data" station_id year url
    body = _http_get_text(url)
    return _parse_ndbc_text(body)
end

# ---------------------------------------------------------------------------
# 3. fetch_tides
# ---------------------------------------------------------------------------

function CSHORE.fetch_tides(station_id::AbstractString;
                             start_date::AbstractString,
                             end_date::AbstractString,
                             datum::AbstractString="NAVD")
    @info "Fetching NOAA CO-OPS tides" station_id start_date end_date datum

    dfmt   = DateFormat("yyyymmdd")
    d_start = Dates.Date(start_date, dfmt)
    d_end   = Dates.Date(end_date,   dfmt)
    d_start > d_end && error("start_date ($start_date) is after end_date ($end_date)")

    frames = DataFrame[]
    chunk_start = d_start
    while chunk_start <= d_end
        chunk_end = min(chunk_start + Dates.Day(30), d_end)
        cs = Dates.format(chunk_start, dfmt)
        ce = Dates.format(chunk_end,   dfmt)
        df_chunk = _fetch_tides_chunk(station_id, cs, ce, datum)
        nrow(df_chunk) > 0 && push!(frames, df_chunk)
        chunk_start = chunk_end + Dates.Day(1)
    end

    length(frames) == 0 && return DataFrame(time = DateTime[], water_level = Float64[])
    return reduce(vcat, frames)
end

function _fetch_tides_chunk(station_id, start_date, end_date, datum)
    url = string(
        "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?",
        "begin_date=", start_date,
        "&end_date=",  end_date,
        "&station=",   station_id,
        "&product=water_level",
        "&datum=",     datum,
        "&units=metric",
        "&time_zone=gmt",
        "&format=json",
        "&application=CSHORE.jl",
    )
    body = _http_get_text(url)
    js   = JSON3.read(body)

    if haskey(js, :error)
        msg = get(get(js, :error, Dict()), :message, "unknown error")
        @warn "CO-OPS API error" station_id start_date end_date msg
        return DataFrame(time = DateTime[], water_level = Float64[])
    end

    data = get(js, :data, nothing)
    data === nothing && return DataFrame(time = DateTime[], water_level = Float64[])

    times  = DateTime[]
    levels = Float64[]
    tfmt   = DateFormat("yyyy-mm-dd HH:MM")
    for rec in data
        t_str = string(get(rec, :t, ""))
        v_str = string(get(rec, :v, ""))
        isempty(t_str) && continue
        isempty(v_str) && continue
        t = try DateTime(t_str, tfmt) catch; continue end
        v = try parse(Float64, v_str) catch; continue end
        push!(times,  t)
        push!(levels, v)
    end
    return DataFrame(time = times, water_level = levels)
end

# ---------------------------------------------------------------------------
# 4. ndbc_to_cshore_bc
# ---------------------------------------------------------------------------

function CSHORE.ndbc_to_cshore_bc(df::DataFrame; angle_convention::Symbol=:nautical)
    required = [:time, :WVHT, :DPD]
    for col in required
        hasproperty(df, col) || error("DataFrame missing required column: $col")
    end

    mask = .!ismissing.(df.WVHT) .& .!ismissing.(df.DPD)
    if hasproperty(df, :MWD)
        mask .&= .!ismissing.(df.MWD)
    end
    sub = df[mask, :]
    nrow(sub) == 0 && error("No valid rows after removing missing WVHT/DPD values")

    t0      = sub.time[1]
    timebc  = Float64[Dates.value(Dates.Second(t - t0)) for t in sub.time]
    hs      = Float64.(sub.WVHT)
    hrmsbc  = hs ./ sqrt(2.0)
    tpbc    = Float64.(sub.DPD)
    swlbc   = zeros(length(timebc))

    if hasproperty(sub, :MWD)
        mwd_raw = Float64.(sub.MWD)
        if angle_convention == :nautical
            wangbc = _nautical_to_math.(mwd_raw)
        else
            wangbc = mwd_raw
        end
    else
        wangbc = zeros(length(timebc))
    end

    return (timebc = timebc, hrmsbc = hrmsbc, tpbc = tpbc,
            swlbc = swlbc, wangbc = wangbc)
end

# ---------------------------------------------------------------------------
# 5. tides_to_swl
# ---------------------------------------------------------------------------

function CSHORE.tides_to_swl(df::DataFrame)
    required = [:time, :water_level]
    for col in required
        hasproperty(df, col) || error("DataFrame missing required column: $col")
    end
    nrow(df) == 0 && return (Float64[], Float64[])

    t0     = df.time[1]
    timebc = Float64[Dates.value(Dates.Second(t - t0)) for t in df.time]
    swlbc  = Float64.(df.water_level)
    return (timebc = timebc, swlbc = swlbc)
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

function _http_get_text(url::AbstractString)
    try
        resp = HTTP.get(url; readtimeout = 60, retry = true, retries = 2,
                        status_exception = true)
        return String(resp.body)
    catch e
        error("HTTP request failed for $url: $(sprint(showerror, e))")
    end
end

function _parse_ndbc_text(body::AbstractString)
    lines = split(body, '\n'; keepempty = false)
    length(lines) < 3 && error("NDBC data has fewer than 3 lines; check station ID")

    header_line = replace(strip(lines[1]), r"^#\s*" => "")
    colnames    = Symbol.(split(header_line))

    data_lines = String[]
    for i in 3:length(lines)
        ln = strip(lines[i])
        startswith(ln, '#') && continue
        isempty(ln) && continue
        push!(data_lines, ln)
    end
    length(data_lines) == 0 && error("No data rows found in NDBC response")

    ncols = length(colnames)
    columns = [Union{Float64,Missing}[] for _ in 1:ncols]

    for ln in data_lines
        parts = split(ln)
        length(parts) < ncols && continue
        for (ci, tok) in enumerate(parts[1:ncols])
            val = tryparse(Float64, tok)
            push!(columns[ci], val === nothing ? missing : val)
        end
    end

    df = DataFrame(columns, colnames)

    if all(c -> hasproperty(df, c), [:YY, :MM, :DD, :hh, :mm])
        df[!, :time] = [
            try
                DateTime(round(Int, r.YY), round(Int, r.MM), round(Int, r.DD),
                         round(Int, r.hh), round(Int, r.mm))
            catch
                missing
            end
            for r in eachrow(df)
        ]
        for c in [:YY, :MM, :DD, :hh, :mm]
            select!(df, Not(c))
        end
    end

    _ndbc_sentinels = Set([99.0, 999.0, 9999.0, 99.00, 999.00, 9999.00])
    for col in names(df)
        col == "time" && continue
        eltype(df[!, col]) <: Union{Missing, Number} || continue
        df[!, col] = [
            (v !== missing && v in _ndbc_sentinels) ? missing : v
            for v in df[!, col]
        ]
    end

    return df
end

function _nautical_to_math(mwd_deg::Real)
    math_deg = 270.0 - mwd_deg
    while math_deg > 180.0;  math_deg -= 360.0; end
    while math_deg < -180.0; math_deg += 360.0; end
    return math_deg
end

end # module CSHOREDataExt
