# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
presets.jl — High-level helpers for common CSHORE setups.

Reduces the typical 30+ lines of build_config boilerplate to 3-5 lines.
Provides profile builders, sediment presets, vegetation/structure presets,
wave forcing builders, and a single-function `quick_run` entry point.
==============================================================================#

# ============================================================================
# Profile builders — return (x::Vector, z::Vector)
# ============================================================================

"""
    planar_beach(; offshore_depth=6.0, slope=0.02, dx=1.0) -> (x, z)

Simple planar beach from `-offshore_depth` to `0` at the given slope.
"""
function planar_beach(; offshore_depth::Float64=6.0, slope::Float64=0.02, dx::Float64=1.0)
    x_end = offshore_depth / slope
    x = collect(0.0:dx:x_end)
    z = [-offshore_depth + slope * xi for xi in x]
    return x, z
end

"""
    beach_dune_profile(; kwargs...) -> (x, z)

Beach + dune profile with offshore shelf, beach face, dune crest, and upland.

Returns `(x, z)` vectors ready for `build_config(bathymetry_x=x, bathymetry_z=z)`.

# Keywords
- `offshore_depth=8.0`: depth at seaward boundary (m, positive down)
- `offshore_dist=600.0`: distance to seaward boundary (m)
- `beach_slope=0.04`: beach face slope
- `beach_width=30.0`: beach face width (m)
- `dune_height=3.0`: dune crest elevation above beach toe (m)
- `dune_width=20.0`: dune crest width (m)
- `upland_elev=nothing`: upland elevation (defaults to dune crest)
- `upland_width=50.0`: flat upland behind dune (m)
- `dx=1.0`: grid spacing (m)
"""
function beach_dune_profile(;
    offshore_depth::Float64=8.0,
    offshore_dist::Float64=600.0,
    beach_slope::Float64=0.04,
    beach_width::Float64=30.0,
    dune_height::Float64=3.0,
    dune_width::Float64=20.0,
    upland_elev::Union{Nothing,Float64}=nothing,
    upland_width::Float64=50.0,
    dx::Float64=1.0,
)
    beach_toe_z = beach_width * beach_slope
    crest_z = beach_toe_z + dune_height
    up_z = upland_elev === nothing ? crest_z : upland_elev

    x_pts = Float64[0.0, offshore_dist]
    z_pts = Float64[-offshore_depth, 0.0]

    # Beach face
    x_beach_end = offshore_dist + beach_width
    push!(x_pts, x_beach_end); push!(z_pts, beach_toe_z)

    # Dune crest
    x_crest = x_beach_end + dune_width / 2
    push!(x_pts, x_crest); push!(z_pts, crest_z)

    # Dune back slope → upland
    x_dune_end = x_beach_end + dune_width
    push!(x_pts, x_dune_end); push!(z_pts, up_z)

    # Upland
    x_end = x_dune_end + upland_width
    push!(x_pts, x_end); push!(z_pts, up_z)

    x = collect(0.0:dx:x_end)
    z = [_preset_interp(x_pts, z_pts, xi) for xi in x]
    return x, z
end

"""
    arctic_bluff_profile(; kwargs...) -> (x, z)

Arctic permafrost bluff with offshore shelf, beach, steep bluff face, and upland tundra.

# Keywords
- `offshore_depth=10.0`: depth at seaward boundary (m)
- `shelf_slope=0.01`: offshore shelf slope
- `beach_slope=0.04`: beach face slope
- `beach_width=25.0`: beach face width (m)
- `bluff_height=5.0`: bluff height above beach (m)
- `bluff_slope=0.25`: bluff face slope (~14°)
- `upland_width=100.0`: flat upland behind bluff (m)
- `dx=1.0`: grid spacing (m)
"""
function arctic_bluff_profile(;
    offshore_depth::Float64=10.0,
    shelf_slope::Float64=0.01,
    beach_slope::Float64=0.04,
    beach_width::Float64=25.0,
    bluff_height::Float64=5.0,
    bluff_slope::Float64=0.25,
    upland_width::Float64=100.0,
    dx::Float64=1.0,
)
    offshore_dist = offshore_depth / shelf_slope
    beach_toe_z = beach_width * beach_slope
    crest_z = beach_toe_z + bluff_height
    bluff_width = bluff_height / bluff_slope

    x_pts = Float64[0.0, offshore_dist,
                     offshore_dist + beach_width,
                     offshore_dist + beach_width + bluff_width,
                     offshore_dist + beach_width + bluff_width + upland_width]
    z_pts = Float64[-offshore_depth, 0.0, beach_toe_z, crest_z, crest_z]

    x = collect(0.0:dx:x_pts[end])
    z = [_preset_interp(x_pts, z_pts, xi) for xi in x]
    return x, z
end

"""
    rocky_shore_profile(; offshore_depth=6.0, platform_width=200.0,
        platform_depth=1.0, cliff_height=8.0, dx=1.0) -> (x, z)

Rocky shore with submerged platform and cliff.
"""
function rocky_shore_profile(;
    offshore_depth::Float64=6.0,
    offshore_dist::Float64=400.0,
    platform_width::Float64=200.0,
    platform_depth::Float64=1.0,
    cliff_height::Float64=8.0,
    upland_width::Float64=50.0,
    dx::Float64=1.0,
)
    x_pts = Float64[0.0, offshore_dist, offshore_dist + platform_width,
                     offshore_dist + platform_width + 10.0,
                     offshore_dist + platform_width + 10.0 + upland_width]
    z_pts = Float64[-offshore_depth, -platform_depth, -platform_depth,
                     cliff_height, cliff_height]
    x = collect(0.0:dx:x_pts[end])
    z = [_preset_interp(x_pts, z_pts, xi) for xi in x]
    return x, z
end

function _preset_interp(xp, zp, xi)
    xi ≤ xp[1] && return zp[1]
    xi ≥ xp[end] && return zp[end]
    for k in 1:(length(xp)-1)
        if xp[k] ≤ xi ≤ xp[k+1]
            t = (xi - xp[k]) / (xp[k+1] - xp[k])
            return zp[k] + t * (zp[k+1] - zp[k])
        end
    end
    return zp[end]
end

