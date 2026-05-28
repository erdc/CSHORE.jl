# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
plotting.jl — Stubs for the optional CairoMakie-backed plotting API.

Real implementations live in `ext/CSHOREMakieExt.jl` and load automatically
when the user has `CairoMakie` available in their environment:

    using CSHORE
    using CairoMakie       # triggers the extension to load
    plot_profile(x, z; ...)

The compiled `cshore-julia` binary (PackageCompiler create_app) does not load
CairoMakie, so the heavy graphics stack is excluded from the bundle.
==============================================================================#

# Generic function declarations — concrete methods are added by the
# CSHOREMakieExt extension when CairoMakie is loaded.
function plot_profile end
function plot_profile_evolution end
function plot_wave_field end
function plot_transport end
function plot_hovmoller end
function plot_mass_balance end
function plot_swash end
function plot_thermal end
function save_figure end

const _PLOTTING_NOT_LOADED_MSG = """
CSHORE plotting requires the CairoMakie extension. Load it with:

    using CairoMakie

(install with `using Pkg; Pkg.add("CairoMakie")` first if needed). The
plotting functions become available automatically once CairoMakie is loaded.
"""

# Catch-all fallback: only fires when no concrete extension method matches.
# When `using CairoMakie` triggers the extension, the more specific methods
# defined there take precedence and the fallback is never invoked.
for _f in (:plot_profile, :plot_profile_evolution, :plot_wave_field,
           :plot_transport, :plot_hovmoller, :plot_mass_balance,
           :plot_swash, :plot_thermal, :save_figure)
    @eval $(_f)(args...; kwargs...) = error(_PLOTTING_NOT_LOADED_MSG)
end
