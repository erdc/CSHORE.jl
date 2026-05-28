# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
cli.jl — Command-line entry point for compiled CSHORE binary.

Invoked by PackageCompiler-built executable as:
    cshore-julia <config-file> [--outdir DIR] [--outfile NAME] [--interval SECS]

Dispatches to the appropriate reader based on file extension:
    *.cshore-julia        → read_cshore-juliajl    (native TOML)
    *.txt           → read_xbeach_params (XBeach params.txt; argument is the
                                          containing directory)
    anything else   → read_infile       (FORTRAN .infile)
==============================================================================#

const _USAGE = """
Usage: cshore-julia <config> [options]

  <config>            Path to a .cshore (TOML), .infile (FORTRAN), or
                      a directory containing an XBeach params.txt.

Options:
  --outdir DIR        Output directory (default: current directory)
  --outfile NAME      NetCDF output filename (default: derived from config)
  --interval SECS     Output write interval in seconds (default: 0 = every step)
  -h, --help          Show this message.
"""

function _parse_cli_args(args::Vector{String})
    config_path = ""
    outdir = "."
    outfile = nothing
    interval = 0.0
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "-h" || a == "--help"
            println(_USAGE)
            exit(0)
        elseif a == "--outdir"
            outdir = args[i+1]; i += 2
        elseif a == "--outfile"
            outfile = args[i+1]; i += 2
        elseif a == "--interval"
            interval = parse(Float64, args[i+1]); i += 2
        elseif startswith(a, "--")
            println(stderr, "Unknown option: $a")
            println(stderr, _USAGE)
            exit(2)
        else
            if isempty(config_path)
                config_path = a
            else
                println(stderr, "Unexpected positional argument: $a")
                exit(2)
            end
            i += 1
        end
    end
    if isempty(config_path)
        println(stderr, _USAGE)
        exit(2)
    end
    return config_path, outdir, outfile, interval
end

function _load_config(path::AbstractString)
    if isdir(path)
        return read_xbeach_params(path)
    end
    ext = lowercase(splitext(path)[2])
    if ext == ".cshore-julia" || ext == ".toml"
        return read_cshore-juliajl(path)
    elseif ext == ".txt"
        return read_xbeach_params(dirname(abspath(path)))
    else
        return read_infile(path)
    end
end

"""
    julia_main()::Cint

Entry point for the PackageCompiler-built executable. Reads `ARGS`, loads the
config, runs the simulation, and returns 0 on success / nonzero on failure.
"""
function julia_main()::Cint
    try
        config_path, outdir, outfile, interval = _parse_cli_args(copy(ARGS))
        @info "CSHORE.jl — loading config" path=config_path
        config = _load_config(config_path)

        if outfile === nothing
            base = splitext(basename(config_path))[1]
            outfile = base * ".nc"
        end
        isdir(outdir) || mkpath(outdir)

        @info "CSHORE.jl — running simulation" outdir=outdir outfile=outfile interval=interval
        run_simulation!(config;
                        outdir=outdir,
                        outfile=outfile,
                        output_interval_s=interval)
        @info "CSHORE.jl — done" output=joinpath(outdir, outfile)
        return Cint(0)
    catch err
        Base.invokelatest(Base.display_error, stderr, err, catch_backtrace())
        return Cint(1)
    end
end
