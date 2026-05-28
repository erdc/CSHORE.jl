# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
io/output.jl — FORTRAN-compatible ASCII output writers.

Reproduces the eight standard CSHORE_USACE output files so Julia runs can be
diffed against the FORTRAN reference. Each file has a header at the top of
each timestep record (time, jmax) and one data row per cross-shore node.

Files written (per BC window when enabled):
  OBPROF  — profile:      x, zb, swl_depth, bed_slope
  OSETUP  — wave setup:   x, h, sigma, wsetup, sigsta
  OCROSS  — cross-shore:  x, sxxsta, tbxsta, efsta, dfsta, dbsta
  OLONG   — longshore:    x, vmean, vstd, sxysta, tbysta   (iangle=1 only)
  OENERG  — energy:       x, hrms, sigma, efsta, dbsta, dfsta, dvegsta
  OTRANS  — transport:    x, qbx_tot, qsx_tot, q_total, vs_tot, ps_tot, pb_tot
  OWETDY  — swash zone:   x, pwet, hwd, uswd, umeawd, ustdwd   (jwd..jdry)
  ODOC    — scalar log:   time, hrms_bc, tp_bc, swl_bc, jr, xr, zr, qo

Format: fixed-width Fortran-style columns (F12.5 for floats), one record
header line per timestep ("# t=… jmax=… [aux]"), then per-node data rows.
Compatible with `numpy.loadtxt(comments='#')` and the original FORTRAN
post-processors. File names are uppercase with `.OUT` suffix; one set per
cross-shore line (suffix `_L{l}.OUT` when iline > 1).
==============================================================================#

"""
    CshoreAsciiWriter

Holds open IOStreams for the eight FORTRAN-compatible output files. Build
via [`open_ascii_outputs`](@ref); flush per timestep via
[`write_ascii_step!`](@ref); close via [`close_ascii_outputs!`](@ref).
"""
mutable struct CshoreAsciiWriter
    outdir::String
    prefix::String
    iline::Int
    # Per-line streams: streams[l] is a Dict{Symbol,IOStream}
    streams::Vector{Dict{Symbol,IOStream}}
    closed::Bool
end

# Standard output keys
const _ASCII_OUTPUTS = (:obprof, :osetup, :ocross, :olong, :oenerg,
                        :otrans, :owetdy, :odoc)

"""
    open_ascii_outputs(config; outdir=".", prefix="", outputs=_ASCII_OUTPUTS,
                       include_longshore=nothing) -> CshoreAsciiWriter

Open the ASCII output files for an upcoming simulation. Skips OLONG by
default unless `config.options.iangle == 1` (or `include_longshore=true`).

Files are opened in write-truncate mode and a `#`-prefixed header is
written identifying the columns. Each file is named
`{prefix}OBPROF.OUT`, `{prefix}OSETUP.OUT`, etc., suffixed with `_L{l}`
when `config.options.iline > 1`.
"""
function open_ascii_outputs(config::CshoreConfig;
                             outdir::AbstractString=".",
                             prefix::AbstractString="",
                             outputs=_ASCII_OUTPUTS,
                             include_longshore::Union{Nothing,Bool}=nothing)
    isdir(outdir) || mkpath(outdir)
    iline = config.options.iline
    want_longshore = include_longshore === nothing ?
                     (config.options.iangle == 1) : include_longshore

    streams = [Dict{Symbol,IOStream}() for _ in 1:iline]
    for l in 1:iline
        suffix = iline > 1 ? "_L$(l)" : ""
        for key in outputs
            key === :olong && !want_longshore && continue
            fname = uppercase(string(key)) * suffix * ".OUT"
            io = open(joinpath(outdir, prefix * fname), "w")
            _write_ascii_header(io, key)
            streams[l][key] = io
        end
    end
    return CshoreAsciiWriter(String(outdir), String(prefix), iline, streams, false)
end

function _write_ascii_header(io::IOStream, key::Symbol)
    headers = Dict(
        :obprof => "# CSHORE OBPROF — profile and bed slope\n" *
                   "# columns: x[m]  zb[m]  swldep[m]  bed_slope\n",
        :osetup => "# CSHORE OSETUP — wave setup\n" *
                   "# columns: x[m]  h[m]  sigma[m]  wsetup[m]  sigsta\n",
        :ocross => "# CSHORE OCROSS — cross-shore momentum\n" *
                   "# columns: x[m]  sxxsta[m]  tbxsta[m]  efsta[m³/s]  dfsta[m³/s²]  dbsta[m³/s²]\n",
        :olong  => "# CSHORE OLONG — longshore (iangle=1)\n" *
                   "# columns: x[m]  vmean[m/s]  vstd[m/s]  sxysta[m]  tbysta[m]\n",
        :oenerg => "# CSHORE OENERG — wave energy and dissipation\n" *
                   "# columns: x[m]  hrms[m]  sigma[m]  efsta[m³/s]  dbsta[m³/s²]  dfsta[m³/s²]  dvegsta[m³/s²]\n",
        :otrans => "# CSHORE OTRANS — sediment transport (summed over fractions)\n" *
                   "# columns: x[m]  qbx[m²/s]  qsx[m²/s]  q_total[m²/s]  vs[m]  ps  pb\n",
        :owetdy => "# CSHORE OWETDY — wet/dry swash zone (jwd..jdry)\n" *
                   "# columns: x[m]  pwet  hwd[m]  uswd[m/s]  umeawd[m/s]  ustdwd[m/s]\n",
        :odoc   => "# CSHORE ODOC — scalar log per BC window\n" *
                   "# columns: time[s]  hrms_bc[m]  tp_bc[s]  swl_bc[m]  jr  xr[m]  zr[m]  qo[m²/s]\n",
    )
    write(io, headers[key])
end

# Fortran-style float: 12-char field, 5 decimals. Handles NaN/Inf cleanly.
@inline _f12(x::Real) = isfinite(x) ? @sprintf("%12.5f", x) : @sprintf("%12s", string(x))
@inline _i8(x::Integer) = @sprintf("%8d", x)

"""
    write_ascii_step!(writer, state, config, time)

Write one record (per cross-shore line) to each open ASCII output file.
Should be called at each BC-window boundary; matches the FORTRAN cadence
of one record per wave-time-step boundary.
"""
function write_ascii_step!(writer::CshoreAsciiWriter, state::CshoreState,
                            config::CshoreConfig, time::Real)
    writer.closed && error("write_ascii_step! called on a closed CshoreAsciiWriter")
    iangle = config.options.iangle
    nf = length(config.multifraction.grain_sizes)
    bc = config.boundary
    itime = max(1, min(state.itime, length(bc.timebc)))

    for l in 1:writer.iline
        jmax_l = state.jmax[l]
        streams = writer.streams[l]

        # Per-step record header (skips ODOC, which is one-line-per-step)
        for key in keys(streams)
            key === :odoc && continue
            io = streams[key]
            @printf(io, "# t=%.5f jmax=%d\n", float(time), jmax_l)
        end

        # OBPROF
        if haskey(streams, :obprof)
            io = streams[:obprof]
            @inbounds for j in 1:jmax_l
                println(io, _f12(state.xb[j]), _f12(state.zb[j, l]),
                            _f12(state.swldep[j, l]),
                            _f12(state.bslope[j, l]))
            end
        end

        # OSETUP
        if haskey(streams, :osetup)
            io = streams[:osetup]
            @inbounds for j in 1:jmax_l
                println(io, _f12(state.xb[j]), _f12(state.h[j]),
                            _f12(state.sigma[j]),
                            _f12(state.wsetup[j]),
                            _f12(state.sigsta[j]))
            end
        end

        # OCROSS
        if haskey(streams, :ocross)
            io = streams[:ocross]
            @inbounds for j in 1:jmax_l
                println(io, _f12(state.xb[j]),
                            _f12(state.sxxsta[j]),
                            _f12(state.tbxsta[j]),
                            _f12(state.efsta[j]),
                            _f12(state.dfsta[j]),
                            _f12(state.dbsta[j]))
            end
        end

        # OLONG (iangle=1 only)
        if haskey(streams, :olong) && iangle == 1
            io = streams[:olong]
            @inbounds for j in 1:jmax_l
                println(io, _f12(state.xb[j]),
                            _f12(state.vmean[j]),
                            _f12(state.vstd[j]),
                            _f12(state.sxysta[j]),
                            _f12(state.tbysta[j]))
            end
        end

        # OENERG
        if haskey(streams, :oenerg)
            io = streams[:oenerg]
            @inbounds for j in 1:jmax_l
                println(io, _f12(state.xb[j]),
                            _f12(state.hrms[j]),
                            _f12(state.sigma[j]),
                            _f12(state.efsta[j]),
                            _f12(state.dbsta[j]),
                            _f12(state.dfsta[j]),
                            _f12(state.dvegsta[j]))
            end
        end

        # OTRANS — sum per-fraction → totals
        if haskey(streams, :otrans)
            io = streams[:otrans]
            @inbounds for j in 1:jmax_l
                qbx_t = 0.0; qsx_t = 0.0; vs_t = 0.0; ps_t = 0.0; pb_t = 0.0
                for k in 1:nf
                    qbx_t += state.qbx[j, k]
                    qsx_t += state.qsx[j, k]
                    vs_t  += state.vs[j, k]
                    ps_t  += state.ps[j, k]
                    pb_t  += state.pb[j, k]
                end
                println(io, _f12(state.xb[j]),
                            _f12(qbx_t), _f12(qsx_t),
                            _f12(state.q_total[j]),
                            _f12(vs_t), _f12(ps_t), _f12(pb_t))
            end
        end

        # OWETDY — only the swash band jwd..jdry
        if haskey(streams, :owetdy)
            io = streams[:owetdy]
            j0 = state.jwd
            j1 = state.jdry
            if j0 > 0 && j1 >= j0
                @inbounds for j in j0:j1
                    println(io, _f12(state.xb[j]),
                                _f12(state.pwet[j]),
                                _f12(state.hwd[j]),
                                _f12(state.uswd[j]),
                                _f12(state.umeawd[j]),
                                _f12(state.ustdwd[j]))
                end
            end
        end

        # ODOC — one line per step (no per-step header)
        if haskey(streams, :odoc)
            io = streams[:odoc]
            qo_l = isempty(state.qo) ? 0.0 : state.qo[l]
            println(io, _f12(float(time)),
                        _f12(bc.hrmsbc[itime]),
                        _f12(bc.tpbc[itime]),
                        _f12(bc.swlbc[itime]),
                        _i8(state.jr),
                        _f12(state.xr),
                        _f12(state.zr),
                        _f12(qo_l))
        end
    end
    return writer
end

"""
    close_ascii_outputs!(writer)

Close all open output streams. Idempotent.
"""
function close_ascii_outputs!(writer::CshoreAsciiWriter)
    writer.closed && return writer
    for streams in writer.streams
        for (_, io) in streams
            close(io)
        end
    end
    writer.closed = true
    return writer
end

"""
    write_outputs(state, config, time; outdir=".", prefix="",
                  writer=nothing, outputs=_ASCII_OUTPUTS)

One-shot ASCII output writer. If `writer` is supplied, appends a record;
otherwise opens fresh files (overwriting), writes one record, and closes.

For multi-step runs, prefer the persistent-writer pattern:

```julia
writer = open_ascii_outputs(config; outdir=".")
try
    for itime in 1:nsteps
        step_bc_window!(state, config, itime)
        write_ascii_step!(writer, state, config, state.time)
    end
finally
    close_ascii_outputs!(writer)
end
```
"""
function write_outputs(state::CshoreState, config::CshoreConfig, time::Real;
                       outdir::AbstractString=".", prefix::AbstractString="",
                       writer::Union{Nothing,CshoreAsciiWriter}=nothing,
                       outputs=_ASCII_OUTPUTS)
    if writer === nothing
        w = open_ascii_outputs(config; outdir=outdir, prefix=prefix, outputs=outputs)
        try
            write_ascii_step!(w, state, config, time)
        finally
            close_ascii_outputs!(w)
        end
    else
        write_ascii_step!(writer, state, config, time)
    end
    return nothing
end
