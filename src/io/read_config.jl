# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

"""
    read_config(path; kwargs...) -> CshoreConfig

Format-dispatched configuration reader. Routes `path` to the appropriate
underlying reader based on what the path points to:

| Path looks like…                      | Reader                |
|---------------------------------------|-----------------------|
| a directory                           | `read_xbeach_params`  |
| a file ending in `.toml` or `.cshore` | `read_cshorejl`       |
| anything else (e.g. `infile`, `*.in`) | `read_infile`         |

Keyword arguments are forwarded to the underlying reader. Each reader
accepts a different subset; unknown kwargs propagate the underlying
MethodError unchanged so the user sees which reader rejected them.

# Examples
```julia
cfg = read_config("test/fixtures/simple.infile")            # → read_infile
cfg = read_config("setups/run01.toml")                       # → read_cshorejl
cfg = read_config("examples/benchmarks/xbeach/Boers_1C")     # → read_xbeach_params
```
"""
function read_config(path::AbstractString; kwargs...)
    if isdir(path)
        return read_xbeach_params(path; kwargs...)
    end
    isfile(path) || throw(ArgumentError("read_config: path not found: $path"))
    ext = lowercase(splitext(path)[2])
    if ext in (".toml", ".cshore")
        return read_cshorejl(path; kwargs...)
    end
    return read_infile(path; kwargs...)
end
