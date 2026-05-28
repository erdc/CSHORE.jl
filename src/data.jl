# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
data.jl — Stubs for the optional environmental-data-fetching API.

Real implementations live in `ext/CSHOREDataExt.jl` and load automatically
when the user has `HTTP`, `JSON3`, `CSV`, AND `DataFrames` available in
their environment:

    using CSHORE
    using HTTP, JSON3, CSV, DataFrames    # triggers the extension
    fetch_ndbc_realtime("44013")

The compiled `cshore-julia` binary does not load these packages, so the
heavy DataFrames / HTTP / CSV / JSON3 stack (~200 MB precompile, hundreds
of MB bundle weight) is excluded from the bundle.
==============================================================================#

# Generic function declarations — concrete methods are added by the
# CSHOREDataExt extension when the four trigger packages are loaded.
function fetch_ndbc_realtime end
function fetch_ndbc_historical end
function fetch_tides end
function ndbc_to_cshore_bc end
function tides_to_swl end

const _DATA_NOT_LOADED_MSG = """
CSHORE environmental-data fetching requires the CSHOREDataExt extension.
Load it with:

    using HTTP, JSON3, CSV, DataFrames

(install with `using Pkg; Pkg.add(["HTTP","JSON3","CSV","DataFrames"])`
first if needed). The fetch_* / ndbc_to_cshore_bc / tides_to_swl functions
become available automatically once all four packages are loaded.
"""

# Catch-all fallbacks: only fire when no concrete extension method matches.
# When `using HTTP, JSON3, CSV, DataFrames` triggers the extension, the
# more specific methods defined there take precedence.
for _f in (:fetch_ndbc_realtime, :fetch_ndbc_historical, :fetch_tides,
           :ndbc_to_cshore_bc, :tides_to_swl)
    @eval $(_f)(args...; kwargs...) = error(_DATA_NOT_LOADED_MSG)
end
