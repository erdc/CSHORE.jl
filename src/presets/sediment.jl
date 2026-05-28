# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
presets.jl — High-level helpers for common CSHORE setups.

Reduces the typical 30+ lines of build_config boilerplate to 3-5 lines.
Provides profile builders, sediment presets, vegetation/structure presets,
wave forcing builders, and a single-function `quick_run` entry point.
==============================================================================#
# ============================================================================
# Sediment presets — return MultifractionConfig
# ============================================================================

"""
    sediment_fine_sand(; d50=0.2e-3) -> MultifractionConfig
Single-fraction fine sand.
"""
sediment_fine_sand(; d50=0.2e-3) = MultifractionConfig(
    grain_sizes=[d50], initial_fractions=[1.0], nlayers=3, layer_thickness=0.1, porosity=0.4)

"""
    sediment_medium_sand(; d50=0.3e-3) -> MultifractionConfig
Single-fraction medium sand (CSHORE default).
"""
sediment_medium_sand(; d50=0.3e-3) = MultifractionConfig(
    grain_sizes=[d50], initial_fractions=[1.0], nlayers=3, layer_thickness=0.1, porosity=0.4)

"""
    sediment_coarse_gravel(; d50=5.0e-3) -> MultifractionConfig
Single-fraction coarse gravel/cobble.
"""
sediment_coarse_gravel(; d50=5.0e-3) = MultifractionConfig(
    grain_sizes=[d50], initial_fractions=[1.0], nlayers=3, layer_thickness=0.1, porosity=0.4)

"""
    sediment_arctic_mix(; nf=3) -> MultifractionConfig
Arctic permafrost bluff: fine silt + medium sand + coarse gravel.
"""
function sediment_arctic_mix(; nf::Int=3)
    if nf == 2
        return MultifractionConfig(
            grain_sizes=[0.063e-3, 0.3e-3],
            initial_fractions=[0.4, 0.6],
            nlayers=3, layer_thickness=0.1, porosity=0.4)
    elseif nf == 3
        return MultifractionConfig(
            grain_sizes=[0.063e-3, 0.3e-3, 2.0e-3],
            initial_fractions=[0.3, 0.5, 0.2],
            nlayers=3, layer_thickness=0.1, porosity=0.4)
    else
        # Linear interpolation between silt and gravel
        sizes = [0.063e-3 * (2.0e-3 / 0.063e-3)^((k-1)/(nf-1)) for k in 1:nf]
        fracs = fill(1.0 / nf, nf)
        return MultifractionConfig(
            grain_sizes=sizes, initial_fractions=fracs,
            nlayers=3, layer_thickness=0.1, porosity=0.4)
    end
end

"""
    sediment_custom(; grain_sizes, fractions, nlayers=3, layer_thickness=0.1, porosity=0.4) -> MultifractionConfig
Custom multifraction sediment.
"""
function sediment_custom(; grain_sizes::Vector{Float64}, fractions::Vector{Float64},
                           nlayers::Int=3, layer_thickness::Float64=0.1, porosity::Float64=0.4)
    return MultifractionConfig(; grain_sizes, initial_fractions=fractions,
                                 nlayers, layer_thickness, porosity)
end

# ============================================================================
# Spatial grain-size distribution
# ============================================================================

"""
    spatially_varying_fractions(x, z, grain_sizes;
        offshore_fracs, beach_fracs, bluff_fracs, transition_width=50.0) -> Matrix{Float64}

Create a spatially-varying initial fraction matrix `(np, nf)` for `init_bed_mass!`.

Uses smooth tanh transitions between offshore/beach/bluff zones based on
elevation. Returns a matrix suitable for the `initial_fractions_spatial`
argument of a future extended `build_config`.

For now, users can call this and then manually set `bed_mass` after initialization.
"""
function spatially_varying_fractions(x::Vector{Float64}, z::Vector{Float64},
                                      grain_sizes::Vector{Float64};
                                      offshore_fracs::Vector{Float64}=[0.2, 0.6, 0.2],
                                      beach_fracs::Vector{Float64}=[0.1, 0.3, 0.6],
                                      bluff_fracs::Vector{Float64}=[0.7, 0.2, 0.1],
                                      transition_width::Float64=50.0)
    nf = length(grain_sizes)
    np = length(x)
    fracs = zeros(np, nf)

    # Find approximate zone boundaries from elevation
    shore_idx = findfirst(zi -> zi ≥ 0.0, z)
    beach_idx = findfirst(zi -> zi ≥ 1.0, z)

    shore_x = shore_idx !== nothing ? x[shore_idx] : x[end] * 0.7
    beach_x = beach_idx !== nothing ? x[beach_idx] : x[end] * 0.85

    tw = transition_width
    for i in 1:np
        # Smooth blend: offshore → beach → bluff
        w_beach = 0.5 * (1.0 + tanh((x[i] - shore_x) / tw))
        w_bluff = 0.5 * (1.0 + tanh((x[i] - beach_x) / tw))
        w_offshore = 1.0 - w_beach

        for k in 1:nf
            fracs[i, k] = w_offshore * offshore_fracs[k] +
                          (w_beach - w_bluff) * beach_fracs[k] +
                          w_bluff * bluff_fracs[k]
        end
        # Normalize
        s = sum(view(fracs, i, :))
        if s > 0
            fracs[i, :] ./= s
        end
    end
    return fracs
end
