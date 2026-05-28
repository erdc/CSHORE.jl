# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

# ThermalState struct definition only. Pulled out of thermal.jl so it can be
# included before state.jl, which lets CshoreState.thermal use the concrete
# Union{Nothing,ThermalState} type instead of Any. Methods that operate on
# ThermalState stay in thermal.jl.

"""
    ThermalState

Mutable runtime state for the thermal submodel. One `(nz, n_rep)`
temperature array, one `(nz, n_rep)` ice-fraction array (1=frozen, 0=thawed),
and the per-shore-node active-layer thickness vector `ALT`. `x_rep` holds
the representative stations' x-coordinates.
"""
mutable struct ThermalState
    T::Matrix{Float64}
    ice_frac::Matrix{Float64}
    ALT::Vector{Float64}
    nz_node::Vector{Int}
    z_top_last::Vector{Float64}
    R_insulation::Vector{Float64}
    Q_thermosyphon::Vector{Float64}
    snow_depth::Vector{Float64}
    R_sod_base::Vector{Float64}
    zb_hard_init::Vector{Float64}
    zb_hard_scour_floor::Vector{Float64}
end
