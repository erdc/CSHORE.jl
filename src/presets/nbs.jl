# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
presets.jl — High-level helpers for common CSHORE setups.

Reduces the typical 30+ lines of build_config boilerplate to 3-5 lines.
Provides profile builders, sediment presets, vegetation/structure presets,
wave forcing builders, and a single-function `quick_run` entry point.
==============================================================================#
# ============================================================================
# Vegetation and structure presets — return VegetationInput or profile modifications
# ============================================================================

"""
    nbs_eelgrass(x, z; depth_range=(0.2, 2.0), cd=0.40, n=500.0, dia=0.003, ht=0.50) -> VegetationInput
Eelgrass (Zostera marina) meadow in the specified depth range.
"""
function nbs_eelgrass(x::Vector{Float64}, z::Vector{Float64};
                       depth_range::Tuple{Float64,Float64}=(0.2, 2.0),
                       cd::Float64=0.40, n::Float64=500.0, dia::Float64=0.003, ht::Float64=0.50)
    _make_depth_veg(x, z; depth_lo=depth_range[1], depth_hi=depth_range[2],
                    cd=cd, n=n, dia=dia, ht=ht)
end

"""
    nbs_kelp(x, z; depth_range=(2.0, 6.0), cd=0.30, n=10.0, dia=0.05, ht=2.0) -> VegetationInput
Kelp forest (Laminaria/Saccharina) in the specified depth range.
"""
function nbs_kelp(x::Vector{Float64}, z::Vector{Float64};
                   depth_range::Tuple{Float64,Float64}=(2.0, 6.0),
                   cd::Float64=0.30, n::Float64=10.0, dia::Float64=0.05, ht::Float64=2.0)
    _make_depth_veg(x, z; depth_lo=depth_range[1], depth_hi=depth_range[2],
                    cd=cd, n=n, dia=dia, ht=ht)
end

"""
    nbs_beach_grass(x, z; x_start, x_end, cd=1.0, n=100.0, dia=0.01, ht=0.25) -> VegetationInput
Beach/bluff grass vegetation over a specified cross-shore range.
"""
function nbs_beach_grass(x::Vector{Float64}, z::Vector{Float64};
                          x_start::Float64, x_end::Float64,
                          cd::Float64=1.0, n::Float64=100.0, dia::Float64=0.01, ht::Float64=0.25,
                          root_depth::Float64=0.15)
    _make_range_veg(x; x_start=x_start, x_end=x_end, cd=cd, n=n, dia=dia, ht=ht,
                    root_depth=root_depth)
end

"""
    nbs_dune_grass(x, z; x_start, x_end, cd=1.2, n=200.0, dia=0.004, ht=0.50) -> VegetationInput
Dense foredune grass (Leymus mollis) over a specified range.
"""
function nbs_dune_grass(x::Vector{Float64}, z::Vector{Float64};
                         x_start::Float64, x_end::Float64,
                         cd::Float64=1.2, n::Float64=200.0, dia::Float64=0.004, ht::Float64=0.50,
                         root_depth::Float64=0.30)
    _make_range_veg(x; x_start=x_start, x_end=x_end, cd=cd, n=n, dia=dia, ht=ht,
                    root_depth=root_depth)
end

"""
    nbs_log_jam(x, z; x_center, width=15.0, cd=0.5, n=3.0, dia=0.15, ht=0.30) -> VegetationInput
Large woody debris (LWD) jam or reef roughness elements.
Default Cd lowered to 0.5 for numerical stability; the KC-dependent
Cd correction in veg_dissipation! provides further dynamic reduction
at high wave heights (He et al. 2019, HESS 25:4825).
"""
function nbs_log_jam(x::Vector{Float64}, z::Vector{Float64};
                      x_center::Float64, width::Float64=15.0,
                      cd::Float64=0.5, n::Float64=3.0, dia::Float64=0.15, ht::Float64=0.30)
    _make_range_veg(x; x_start=x_center - width/2, x_end=x_center + width/2,
                    cd=cd, n=n, dia=dia, ht=ht)
end

"""
    nbs_breakwater(x, z; center_depth=3.0, width=30.0, crest_elev=0.5) -> (x_mod, z_mod, hardbottom_z)
Emergent or submerged breakwater. Returns modified profile + hardbottom.
"""
function nbs_breakwater(x::Vector{Float64}, z::Vector{Float64};
                         center_depth::Float64=3.0, width::Float64=30.0, crest_elev::Float64=0.5)
    # Find x position at the target depth
    center_x = 0.0
    for i in 1:length(x)-1
        if z[i] ≤ -center_depth && z[i+1] > -center_depth
            center_x = x[i] + (-center_depth - z[i]) / (z[i+1] - z[i]) * (x[i+1] - x[i])
            break
        end
    end

    x_mod = copy(x)
    z_mod = copy(z)
    zh = fill(-1e30, length(x))

    for i in eachindex(x)
        if abs(x[i] - center_x) ≤ width / 2
            z_mod[i] = max(z_mod[i], crest_elev)
            zh[i] = crest_elev  # non-erodible
        end
    end
    return x_mod, z_mod, zh
end

"""
    nbs_nourishment(x, z; width=100.0, crest_elev=1.0, ramp=20.0) -> (x_mod, z_mod)
Beach nourishment — widens the beach by adding a fill berm.
"""
function nbs_nourishment(x::Vector{Float64}, z::Vector{Float64};
                          width::Float64=100.0, crest_elev::Float64=1.0, ramp::Float64=20.0)
    # Find the shoreline (z=0)
    shore_x = x[end]
    for i in 1:length(x)-1
        if z[i] ≤ 0 && z[i+1] > 0
            shore_x = x[i]
            break
        end
    end

    z_mod = copy(z)
    fill_start = shore_x - width
    for i in eachindex(x)
        if x[i] ≥ fill_start && x[i] ≤ shore_x
            frac = (x[i] - fill_start) / width
            z_fill = crest_elev * min(frac * width / ramp, 1.0, (shore_x - x[i]) / ramp)
            z_fill = crest_elev * clamp(frac, 0.0, 1.0)
            z_mod[i] = max(z_mod[i], z_fill)
        end
    end
    return copy(x), z_mod
end

# Helpers
"""
    _taper(x, x_lo, x_hi, taper_width) -> Float64

Smooth tanh taper: 0 outside [x_lo, x_hi], 1 well inside, with smooth
transitions of width `taper_width` at each edge.  Prevents step-function
discontinuities in vegetation density that cause numerical instability
in the sediment transport.

    _taper(x, lo, hi, tw) = 0.5*(tanh((x-lo)/tw) - tanh((x-hi)/tw))

For tw → 0 this recovers a sharp box function; tw ≈ 5-10 m is typical.
"""
function _taper(xi::Float64, x_lo::Float64, x_hi::Float64, tw::Float64)
    tw = max(tw, 0.1)   # avoid division by zero
    return 0.5 * (tanh((xi - x_lo) / tw) - tanh((xi - x_hi) / tw))
end

# Default taper width: 5 m (5 grid cells at dx=1m). Can be overridden
# by passing taper_width to any preset function.
const DEFAULT_TAPER_WIDTH = 5.0

function _make_depth_veg(x, z; depth_lo, depth_hi, cd, n, dia, ht,
                          taper_width=DEFAULT_TAPER_WIDTH,
                          root_depth::Float64=0.0)
    np = length(x)
    vegn = zeros(np, 1); vegb = zeros(np, 1)
    vegd = zeros(np, 1); vegh = zeros(np, 1)
    vegrd = zeros(np, 1); vegrh = zeros(np, 1)
    for i in 1:np
        depth = -z[i]  # positive depth below MSL
        w = _taper(depth, depth_lo, depth_hi, taper_width)
        if w > 0.01
            vegn[i, 1] = n * w; vegb[i, 1] = dia
            vegd[i, 1] = ht * w; vegh[i, 1] = ht * w
            vegrd[i, 1] = root_depth * w
            vegrh[i, 1] = root_depth * w
        end
    end
    return VegetationInput(vegcd=cd, vegcdm=cd, vegn=vegn, vegb=vegb,
                           vegd=vegd, vegh=vegh, vegrd=vegrd, vegrh=vegrh)
end

function _make_range_veg(x; x_start, x_end, cd, n, dia, ht,
                          taper_width=DEFAULT_TAPER_WIDTH,
                          root_depth::Float64=0.0)
    np = length(x)
    vegn = zeros(np, 1); vegb = zeros(np, 1)
    vegd = zeros(np, 1); vegh = zeros(np, 1)
    vegrd = zeros(np, 1); vegrh = zeros(np, 1)
    for i in 1:np
        w = _taper(x[i], x_start, x_end, taper_width)
        if w > 0.01
            vegn[i, 1] = n * w; vegb[i, 1] = dia
            vegd[i, 1] = ht * w; vegh[i, 1] = ht * w
            vegrd[i, 1] = root_depth * w
            vegrh[i, 1] = root_depth * w
        end
    end
    return VegetationInput(vegcd=cd, vegcdm=cd, vegn=vegn, vegb=vegb,
                           vegd=vegd, vegh=vegh, vegrd=vegrd, vegrh=vegrh)
end

# ============================================================================
# Coastal vegetation roughness lookup
# ----------------------------------------------------------------------------
# Cd, stem density, blade width, and canopy height for common coastal NbS
# vegetation categories, indexed by (category, region). Cd values reflect
# published literature ranges (Anderson & Smith 2014; Augustin 2009; Dasgupta
# 2019; Jadhav 2012; Mudd 2010; Mendez & Losada 2004; Innocenti 2021; etc.);
# stem-scale parameters are mid-range picks from the same studies. Override
# via `nbs_vegetation` kwargs.
#
# CSHORE's vegetation drag (Mendez & Losada 2004) is Cd-based. Equivalent
# depth-averaged Manning's n ranges from the same literature are kept in
# `NBS_VEGETATION_MANNING_RANGE` for cross-comparison with depth-averaged
# 2D models — they can also be fed into CSHORE's bottom friction via the
# `manning=` kwarg on `build_config` (see ifriction_spatial=2).
#
# Region tags: :east_gulf, :west, :gulf.
# ============================================================================

"""
Representative `(cd, n, dia, ht, root_depth)` parameters for common coastal
vegetation categories. Mid-range literature values; override per-call via
`nbs_vegetation` kwargs.

`root_depth` (m) drives the root-biomass transport reduction in transport.jl
(reduces qbx/qsx/q_total by up to 30% via the Okin 2008 / Chen et al. 2012
soil-binding effect). Picks reflect typical effective rooting depths for the
dominant species in each category (Mudd 2010; Wigand 2014; Charbonneau 2017;
Feagin 2015 — and references therein for dune and mangrove rooting).
"""
const NBS_VEGETATION_PARAMS = Dict{Tuple{Symbol,Symbol}, NamedTuple{(:cd,:n,:dia,:ht,:root_depth),NTuple{5,Float64}}}(
    # Mangroves — Rhizophora mangle (E), Avicennia germinans + R. mangle (Gulf)
    # Cd 0.4–10 (Noarayanan 2012; Dasgupta 2019; Vanegas 2019; He 2019).
    # Prop-root + cable-root system extends 1–2 m below the substrate.
    (:mangrove,           :east_gulf) => (cd=1.5, n=2.0,    dia=0.10,  ht=2.0,  root_depth=1.50),
    (:mangrove,           :gulf)      => (cd=1.5, n=2.0,    dia=0.10,  ht=2.0,  root_depth=1.50),
    # Low marsh — Spartina alterniflora / Sporobolus alterniflorus.
    # Cd 0.1–4.3 (Augustin 2009; Anderson & Smith 2014; Jadhav 2012).
    # Dense fibrous root + rhizome mat; effective rooting depth ~0.3–0.5 m
    # for the tall form (Mudd 2010), ~0.2–0.4 m for the interior short form.
    (:low_marsh_tall,     :east_gulf) => (cd=1.5, n=200.0,  dia=0.006, ht=1.5,  root_depth=0.40),
    (:low_marsh_short,    :east_gulf) => (cd=1.5, n=600.0,  dia=0.004, ht=0.4,  root_depth=0.30),
    (:low_marsh,          :west)      => (cd=1.5, n=400.0,  dia=0.003, ht=1.0,  root_depth=0.30),
    # Shrubby low-marsh fringe — Iva fructescens, Baccharis halimifolia.
    # Cd 0.35–0.85 (Wunder 2011). Woody taproots reach 0.4–0.8 m.
    (:marsh_shrub,        :east_gulf) => (cd=0.6, n=50.0,   dia=0.020, ht=1.5,  root_depth=0.50),
    # Invasive emergent — Phragmites australis / Typha angustifolia.
    # Cd 0.61–26.24 (Kamali 2018; Zhao 2017). Deep rhizomes, 0.5–1.0 m.
    (:phragmites,         :east_gulf) => (cd=3.0, n=200.0,  dia=0.015, ht=2.5,  root_depth=0.70),
    # High marsh — Spartina patens, Distichlis spicata, Juncus spp. (E/Gulf).
    # Cd 0.1–4.3 (Augustin 2009; Jadhav 2012; Peruzzo 2018). Shallower mat
    # than low-marsh forms; ~0.2–0.4 m.
    (:high_marsh,         :east_gulf) => (cd=1.0, n=1000.0, dia=0.002, ht=0.3,  root_depth=0.25),
    # High marsh — Salicornia pacifica, Frankenia salina (West). Shallow
    # rooting (0.1–0.2 m) typical of succulent halophytes (Wigand 2014).
    (:high_marsh,         :west)      => (cd=0.5, n=800.0,  dia=0.003, ht=0.2,  root_depth=0.15),
    # Beach / foredune grass — Uniola paniculata, Ammophila breviligulata (E),
    # Ammophila arenaria, Leymus mollis (W), Uniola/Panicum (Gulf).
    # Cd 0.1–1.1 (Sepaskhah & Bondar 2002; Augustin 2009). Deep rhizome
    # networks (0.3–0.6 m) are central to dune stabilization (Charbonneau 2017).
    (:dune_grass,         :east_gulf) => (cd=0.5, n=300.0,  dia=0.003, ht=0.7,  root_depth=0.40),
    (:dune_grass,         :west)      => (cd=0.5, n=250.0,  dia=0.004, ht=0.8,  root_depth=0.40),
    (:dune_grass,         :gulf)      => (cd=0.5, n=300.0,  dia=0.003, ht=0.7,  root_depth=0.40),
    # Upper-dune shrub thicket — Lupinus arboreus, Baccharis pilularis.
    # Woody taproots 0.5–1.0 m; Cd from willow-analog (Chiaradia 2019).
    (:dune_shrub,         :west)      => (cd=1.0, n=30.0,   dia=0.025, ht=2.0,  root_depth=0.60),
    # Beach forbs / sand-binding herbs — Sesuvium, Cakile, Ipomoea, Solidago.
    # Cd 0.034–1.7 individual stem (Innocenti 2021); up to ~20 (Cantalice 2015).
    # Shallow taproots / runners, 0.05–0.15 m.
    (:beach_forbs,        :east_gulf) => (cd=0.5, n=100.0,  dia=0.005, ht=0.15, root_depth=0.10),
    (:beach_forbs,        :west)      => (cd=0.5, n=100.0,  dia=0.005, ht=0.15, root_depth=0.10),
    (:beach_forbs,        :gulf)      => (cd=0.5, n=100.0,  dia=0.005, ht=0.15, root_depth=0.10),
)

"""
Manning's n literature ranges for the same `(category, region)` keys as
`NBS_VEGETATION_PARAMS`. CSHORE's vegetation drag is Cd-based, but these
ranges can be (a) compared with depth-averaged 2D model parameterizations,
or (b) passed into `build_config(manning=...)` to drive bottom friction via
the dynamic Manning option (`ifriction_spatial=2`).
"""
const NBS_VEGETATION_MANNING_RANGE = Dict{Tuple{Symbol,Symbol}, Tuple{Float64,Float64}}(
    (:mangrove,        :east_gulf) => (0.124, 3.00),
    (:mangrove,        :gulf)      => (0.124, 3.00),
    (:low_marsh_tall,  :east_gulf) => (0.43,  4.45),
    (:low_marsh_short, :east_gulf) => (0.43,  4.45),
    (:low_marsh,       :west)      => (0.018, 0.024),
    (:marsh_shrub,     :east_gulf) => (0.043, 0.056),
    (:phragmites,      :east_gulf) => (0.043, 0.056),
    (:high_marsh,      :east_gulf) => (0.018, 0.024),
    (:high_marsh,      :west)      => (0.018, 0.024),
    (:dune_grass,      :east_gulf) => (0.043, 0.136),
    (:dune_grass,      :west)      => (0.043, 0.136),
    (:dune_grass,      :gulf)      => (0.043, 0.136),
    (:dune_shrub,      :west)      => (0.044, 0.067),
    (:beach_forbs,     :east_gulf) => (0.0196, 0.0238),
    (:beach_forbs,     :west)      => (0.0196, 0.0238),
    (:beach_forbs,     :gulf)      => (0.0196, 0.0238),
)

"""
    nbs_vegetation_params(category::Symbol; region::Symbol=:east_gulf) -> NamedTuple

Return the `(cd, n, dia, ht, root_depth)` defaults for a `(category, region)`
entry from `NBS_VEGETATION_PARAMS`. Throws if the pair is not in the lookup.
"""
function nbs_vegetation_params(category::Symbol; region::Symbol=:east_gulf)
    key = (category, region)
    haskey(NBS_VEGETATION_PARAMS, key) || throw(ArgumentError(
        "No vegetation entry for ($category, $region). Available: $(collect(keys(NBS_VEGETATION_PARAMS)))"))
    return NBS_VEGETATION_PARAMS[key]
end

"""
    nbs_vegetation(category, x, z; region=:east_gulf, placement=:auto,
                   depth_range=nothing, x_start=nothing, x_end=nothing,
                   cd=nothing, n=nothing, dia=nothing, ht=nothing,
                   taper_width=DEFAULT_TAPER_WIDTH,
                   root_depth=nothing) -> VegetationInput

Build a `VegetationInput` for a coastal vegetation category. Defaults come
from `NBS_VEGETATION_PARAMS` and can be overridden per-parameter. Pass
`root_depth=0.0` explicitly to opt out of the root-binding transport
reduction (transport.jl) — the lookup-derived default is non-zero for
every category.

`placement`:
- `:auto` (default) — depth-based for submerged categories (mangrove,
  low_marsh*, high_marsh, phragmites, marsh_shrub) and x-range-based for
  emergent dune/beach categories (dune_grass, dune_shrub, beach_forbs).
- `:depth` — uses `depth_range = (lo, hi)` (m below MSL, positive down).
- `:range` — uses `x_start`, `x_end` (m).

# Examples
```julia
veg = nbs_vegetation(:low_marsh_tall, x, z; region=:east_gulf)
veg = nbs_vegetation(:dune_grass, x, z; region=:west, x_start=180.0, x_end=210.0)
veg = nbs_vegetation(:mangrove, x, z; region=:gulf, depth_range=(0.0, 1.5), cd=2.5)
```
"""
function nbs_vegetation(category::Symbol, x::Vector{Float64}, z::Vector{Float64};
                  region::Symbol=:east_gulf,
                  placement::Symbol=:auto,
                  depth_range::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                  x_start::Union{Nothing,Float64}=nothing,
                  x_end::Union{Nothing,Float64}=nothing,
                  cd::Union{Nothing,Float64}=nothing,
                  n::Union{Nothing,Float64}=nothing,
                  dia::Union{Nothing,Float64}=nothing,
                  ht::Union{Nothing,Float64}=nothing,
                  taper_width::Float64=DEFAULT_TAPER_WIDTH,
                  root_depth::Union{Nothing,Float64}=nothing)
    p = nbs_vegetation_params(category; region=region)
    cd_         = something(cd,         p.cd)
    n_          = something(n,          p.n)
    dia_        = something(dia,        p.dia)
    ht_         = something(ht,         p.ht)
    root_depth_ = something(root_depth, p.root_depth)

    # Default placement zones — coastal zonation literature.
    # Depth = m below MSL (positive down).
    default_depth = Dict(
        :mangrove        => (0.0, 1.5),
        :low_marsh_tall  => (0.0, 1.0),
        :low_marsh_short => (-0.2, 0.6),
        :low_marsh       => (0.0, 1.0),
        :marsh_shrub     => (-0.5, 0.2),
        :phragmites      => (-0.5, 0.8),
        :high_marsh      => (-0.6, 0.1),
    )
    emergent_categories = (:dune_grass, :dune_shrub, :beach_forbs)

    use_depth = placement === :depth ||
                (placement === :auto && !(category in emergent_categories))

    if use_depth
        dr = something(depth_range, get(default_depth, category, (0.0, 1.5)))
        return _make_depth_veg(x, z; depth_lo=dr[1], depth_hi=dr[2],
                               cd=cd_, n=n_, dia=dia_, ht=ht_,
                               taper_width=taper_width, root_depth=root_depth_)
    else
        (x_start === nothing || x_end === nothing) && throw(ArgumentError(
            "nbs_vegetation(:$category) with :range placement requires x_start and x_end"))
        return _make_range_veg(x; x_start=x_start, x_end=x_end,
                               cd=cd_, n=n_, dia=dia_, ht=ht_,
                               taper_width=taper_width, root_depth=root_depth_)
    end
end

# ============================================================================
# Manning's-n path for vegetation (alternative to Cd-based vegetation drag)
# ----------------------------------------------------------------------------
# `nbs_vegetation_manning_field` returns a per-node Manning's n vector that
# already encodes the vegetation's frictional effect (per the literature
# ranges in NBS_VEGETATION_MANNING_RANGE). It is wrapped in a
# `VegetatedManningField` type so that `build_config` can refuse to combine
# it with a `VegetationInput`, which would double-count the vegetation drag.
#
# Use this when you want the depth-averaged 2D modeling approach (§2.2 of
# the lookup-table report): vegetation as bulk bed roughness, no IVEG.
# For the Cd-based stem-drag approach (§3.0), use `nbs_vegetation` with a
# bare-sand `friction` / `manning` baseline instead.
# ============================================================================

"""
    nbs_vegetation_manning_field(category, x, z;
        region=:east_gulf, placement=:auto,
        depth_range=nothing, x_start=nothing, x_end=nothing,
        n=nothing, n_bare=0.02, taper_width=DEFAULT_TAPER_WIDTH,
        base=nothing) -> VegetatedManningField

Build a per-node Manning's n field for a vegetation category. Inside the
placement zone the value blends to the category's representative n (mid of
`NBS_VEGETATION_MANNING_RANGE`, or user-supplied `n`); outside it falls back
to `n_bare`.

Pass `base=existing_field` to add another vegetation zone on top of a prior
field (the max of each node is taken). This is how to compose multiple
vegetation categories on the same profile.

Placement semantics match [`nbs_vegetation`](@ref): submerged categories
default to depth-based placement, emergent categories require `x_start` /
`x_end`.

The result must be passed to `build_config(manning=...)` AS-IS — do not
combine with `vegetation=` or `build_config` will throw.
"""
function nbs_vegetation_manning_field(category::Symbol,
                                       x::Vector{Float64}, z::Vector{Float64};
                                       region::Symbol=:east_gulf,
                                       placement::Symbol=:auto,
                                       depth_range::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                                       x_start::Union{Nothing,Float64}=nothing,
                                       x_end::Union{Nothing,Float64}=nothing,
                                       n::Union{Nothing,Float64}=nothing,
                                       n_bare::Float64=0.02,
                                       taper_width::Float64=DEFAULT_TAPER_WIDTH,
                                       base::Union{Nothing,VegetatedManningField}=nothing)
    key = (category, region)
    haskey(NBS_VEGETATION_MANNING_RANGE, key) || throw(ArgumentError(
        "No Manning's n range for ($category, $region). Available: $(collect(keys(NBS_VEGETATION_MANNING_RANGE)))"))
    n_lo, n_hi = NBS_VEGETATION_MANNING_RANGE[key]
    n_veg = something(n, 0.5 * (n_lo + n_hi))

    # Placement logic — same defaults as nbs_vegetation
    default_depth = Dict(
        :mangrove        => (0.0, 1.5),
        :low_marsh_tall  => (0.0, 1.0),
        :low_marsh_short => (-0.2, 0.6),
        :low_marsh       => (0.0, 1.0),
        :marsh_shrub     => (-0.5, 0.2),
        :phragmites      => (-0.5, 0.8),
        :high_marsh      => (-0.6, 0.1),
    )
    emergent_categories = (:dune_grass, :dune_shrub, :beach_forbs)
    use_depth = placement === :depth ||
                (placement === :auto && !(category in emergent_categories))

    np = length(x)
    out = base === nothing ? fill(n_bare, np) : copy(base.values)
    cats = base === nothing ? Tuple{Symbol,Symbol}[] : copy(base.categories)

    if use_depth
        dr = something(depth_range, get(default_depth, category, (0.0, 1.5)))
        for i in 1:np
            depth = -z[i]
            w = _taper(depth, dr[1], dr[2], taper_width)
            if w > 0.01
                n_blended = n_bare + w * (n_veg - n_bare)
                out[i] = max(out[i], n_blended)   # compose by max
            end
        end
    else
        (x_start === nothing || x_end === nothing) && throw(ArgumentError(
            "nbs_vegetation_manning_field(:$category) with :range placement requires x_start and x_end"))
        for i in 1:np
            w = _taper(x[i], x_start, x_end, taper_width)
            if w > 0.01
                n_blended = n_bare + w * (n_veg - n_bare)
                out[i] = max(out[i], n_blended)
            end
        end
    end
    push!(cats, key)
    return VegetatedManningField(out, cats)
end

# ============================================================================
# Porous layer presets — permeable revetments and gravel beaches
# ============================================================================

"""
    nbs_gravel_revetment(x, z; thickness=1.0, Dn50=0.05, porosity=0.4,
                          x_start=nothing, x_end=nothing) -> (porous_z, PorousInput)

Gravel or cobble dynamic revetment. Returns a porous layer bottom elevation
vector and a `PorousInput` struct ready to pass to `build_config`.

The porous layer extends from the bed surface downward by `thickness` metres
within the specified cross-shore range. Outside that range, `porous_z = zb`
(zero-thickness layer, no porous effect). If `x_start`/`x_end` are not
specified, the layer is applied wherever the bed is above the offshore flat.

# Arguments
- `thickness`: porous layer depth below the bed surface (m)
- `Dn50`: nominal stone diameter (m). Typical: 0.02 gravel, 0.05 cobble, 0.10 rubble
- `porosity`: pore fraction (0–1). Typical: 0.35–0.45 for gravel
- `x_start`, `x_end`: cross-shore extent of the revetment (m)

# Example
```julia
porous_z, por = nbs_gravel_revetment(x, z; thickness=0.8, Dn50=0.04)
cfg = build_config(..., iperm=1, porous=por, hardbottom_z=porous_z)
# Or more simply:
cfg = build_config(..., iperm=1, porous_z=porous_z, stone_diameter=0.04)
```
"""
function nbs_gravel_revetment(x::Vector{Float64}, z::Vector{Float64};
                               thickness::Float64=1.0,
                               Dn50::Float64=0.05,
                               porosity::Float64=0.4,
                               x_start::Union{Nothing,Float64}=nothing,
                               x_end::Union{Nothing,Float64}=nothing)
    np = length(x)
    porous_z = copy(z)   # default: zp = zb (hp=0, no porous layer)

    # Default range: wherever bed is above the offshore flat
    z_flat = minimum(z[1:min(50, np)])
    xs = x_start === nothing ? x[findfirst(zi -> zi > z_flat + 0.01, z)] : x_start
    xe = x_end === nothing ? x[end] : x_end

    for i in 1:np
        if x[i] >= xs && x[i] <= xe
            porous_z[i] = z[i] - thickness
        end
    end

    por = PorousInput(x, porous_z; porosity=porosity, stone_diameter=Dn50)
    return porous_z, por
end

# ============================================================================
# Thermal intervention presets — permafrost management via the heat equation
# ============================================================================

"""
    nbs_thermosyphon(x; x_start, x_end, Q_peak=100.0, spacing=3.0, radius=3.0) -> Vector{Float64}

Passive thermosyphon array for permafrost protection. Returns a per-node
heat extraction rate vector `Q_thermosyphon` (W/m²) to be assigned to
`state.thermal.Q_thermosyphon` after initialization.

Pipes are placed every `spacing` metres within `[x_start, x_end]`. Each
pipe's influence decays as a Gaussian with e-folding length `radius`:

    Q(x) = Σ_pipes  Q_peak · exp(-(x - x_pipe)² / radius²)

The thermosyphon is a one-way device: extraction is active only in winter
when the air is colder than the ground surface (`T_surface < T[1]`), and
inactive in summer. This physics is enforced inside `_step_heat_column!`.

# Typical values
- `Q_peak ≈ 50–150 W/m²` (peak extraction at the pipe)
- `spacing ≈ 2–5 m` (closer = more overlap = stronger aggregate cooling)
- `radius ≈ 2–5 m` (influence radius in cross-shore direction)

# Example
```julia
Q = nbs_thermosyphon(x; x_start=3850.0, x_end=3870.0)
state.thermal.Q_thermosyphon .= Q[1:state.jmax[1]]
```
"""
function nbs_thermosyphon(x::Vector{Float64};
                           x_start::Float64, x_end::Float64,
                           Q_peak::Float64=100.0,
                           spacing::Float64=3.0,
                           radius::Float64=3.0)
    np = length(x)
    Q = zeros(Float64, np)
    # Place pipes at regular intervals within the specified range
    pipe_positions = collect(x_start:spacing:x_end)
    for x_pipe in pipe_positions
        @inbounds for i in 1:np
            dx = x[i] - x_pipe
            Q[i] += Q_peak * exp(-dx * dx / (radius * radius))
        end
    end
    return Q
end

"""
    nbs_sod_insulation(x; x_start, x_end, thickness=0.15, k_sod=0.15) -> Vector{Float64}

Sod or moss insulation layer for permafrost protection. Returns a per-node
thermal resistance vector `R_insulation` (m²·K/W) to be assigned to
`state.thermal.R_insulation` after initialization.

The organic layer sits between the atmosphere and the mineral soil surface,
adding a series thermal resistance `R = thickness / k_sod` that slows heat
flux into the ground. This reduces summer thaw depth without any moving
parts or energy input — the insulation works year-round, slowing both
warming AND cooling. In practice, the net effect is permafrost-protective
because summer thaw is the limiting process.

# Typical values
- `thickness ≈ 0.10–0.30 m` (sod/moss mat depth)
- `k_sod ≈ 0.10–0.25 W/m/K` (organic soil thermal conductivity;
  lower = more insulating; mineral soil is ~0.8–1.5 for comparison)

# Example
```julia
R = nbs_sod_insulation(x; x_start=3850.0, x_end=3870.0)
state.thermal.R_insulation .= R[1:state.jmax[1]]
```
"""
function nbs_sod_insulation(x::Vector{Float64};
                             x_start::Float64, x_end::Float64,
                             thickness::Float64=0.15,
                             k_sod::Float64=0.15)
    np = length(x)
    R = zeros(Float64, np)
    R_val = thickness / k_sod   # m²·K/W
    @inbounds for i in 1:np
        if x[i] >= x_start && x[i] <= x_end
            R[i] = R_val
        end
    end
    return R
end

# ============================================================================
# Snow NbS presets — spatial snow management (fences, drifts, clearing)
# ----------------------------------------------------------------------------
# Each preset returns a `SnowSpatialModifier` that overrides snow depth in a
# spatial zone (per-node min/max applied inside `update_snow!`). Compose by
# passing `base=existing_modifier` to a subsequent call. Pass the final
# modifier to `build_config(snow_modifier=...)`.
# ============================================================================

# Snow modifier helpers
# Sentinel "no per-node cap" — large enough that snow_config.max_depth always
# wins the final clamp, but finite so taper blends don't overflow.
const SNOW_NO_CAP = 100.0

function _empty_snow_modifier(np::Int)
    SnowSpatialModifier(zeros(Float64, np), fill(SNOW_NO_CAP, np))
end

function _apply_min_in_zone!(depth_min::Vector{Float64}, x::Vector{Float64},
                              x_start, x_end, target_depth, tw)
    np = length(x)
    @inbounds for i in 1:np
        w = _taper(x[i], x_start, x_end, tw)
        if w > 0.01
            depth_min[i] = max(depth_min[i], w * target_depth)
        end
    end
end

function _apply_max_in_zone!(depth_max::Vector{Float64}, x::Vector{Float64},
                              x_start, x_end, cap_depth, tw, blend_range)
    # The cap fades from cap_depth (zone center, w=1) up to
    # `cap_depth + blend_range` at zone edges (w→0), then to the prior cap
    # outside the zone. `blend_range` should be comparable to snow_config.max_depth
    # so the cap is essentially inactive at the edges.
    np = length(x)
    @inbounds for i in 1:np
        w = _taper(x[i], x_start, x_end, tw)
        if w > 0.01
            zone_cap = cap_depth + (1.0 - w) * blend_range
            depth_max[i] = min(depth_max[i], zone_cap)
        end
    end
end

"""
    nbs_snow_fence(x; x_start, x_end, drift_depth=0.4,
                   taper_width=DEFAULT_TAPER_WIDTH,
                   base=nothing) -> SnowSpatialModifier

Snow fence — forces a minimum snow depth in the lee zone `[x_start, x_end]`
to represent engineered drift accumulation. `drift_depth` (m) is the
characteristic accumulated depth in the fence shadow (typical 0.3-1.0 m
for a 1-2 m fence). Tapers smoothly to zero at zone edges.

The forced minimum is applied each step inside `update_snow!`, so the
fence's drift persists even when the degree-day model would otherwise
melt off the snow (e.g. brief warm-ups). Final cap of
`snow_config.max_depth` still applies on top — bump `max_depth` if you
need drifts deeper than 0.3 m.
"""
function nbs_snow_fence(x::Vector{Float64};
                         x_start::Float64, x_end::Float64,
                         drift_depth::Float64=0.4,
                         taper_width::Float64=DEFAULT_TAPER_WIDTH,
                         base::Union{Nothing,SnowSpatialModifier}=nothing)
    np = length(x)
    mod = base === nothing ? _empty_snow_modifier(np) :
          SnowSpatialModifier(copy(base.depth_min), copy(base.depth_max))
    _apply_min_in_zone!(mod.depth_min, x, x_start, x_end, drift_depth, taper_width)
    return mod
end

"""
    nbs_snow_clearing(x; x_start, x_end, max_depth=0.0,
                      taper_width=DEFAULT_TAPER_WIDTH,
                      base=nothing) -> SnowSpatialModifier

Snow clearing — caps snow depth in zone `[x_start, x_end]` at `max_depth`
(default 0.0 = fully cleared). Use to model winter snow removal over
thermosyphon arrays or other engineered surfaces where you want to expose
the ground for maximum heat loss.
"""
function nbs_snow_clearing(x::Vector{Float64};
                            x_start::Float64, x_end::Float64,
                            max_depth::Float64=0.0,
                            taper_width::Float64=DEFAULT_TAPER_WIDTH,
                            blend_range::Float64=1.0,
                            base::Union{Nothing,SnowSpatialModifier}=nothing)
    np = length(x)
    mod = base === nothing ? _empty_snow_modifier(np) :
          SnowSpatialModifier(copy(base.depth_min), copy(base.depth_max))
    _apply_max_in_zone!(mod.depth_max, x, x_start, x_end, max_depth,
                         taper_width, blend_range)
    return mod
end

"""
    nbs_insulating_drift(x; x_start, x_end, target_depth=0.5,
                         taper_width=DEFAULT_TAPER_WIDTH,
                         base=nothing) -> SnowSpatialModifier

Engineered / maintained insulating snow drift — like a snow fence but
intended for sustained winter insulation of underlying permafrost or
sediment. Mechanically identical to `nbs_snow_fence` (forced minimum
depth) but semantically distinct. `target_depth` typically 0.4-1.0 m for
permafrost thermal protection applications.
"""
function nbs_insulating_drift(x::Vector{Float64};
                               x_start::Float64, x_end::Float64,
                               target_depth::Float64=0.5,
                               taper_width::Float64=DEFAULT_TAPER_WIDTH,
                               base::Union{Nothing,SnowSpatialModifier}=nothing)
    np = length(x)
    mod = base === nothing ? _empty_snow_modifier(np) :
          SnowSpatialModifier(copy(base.depth_min), copy(base.depth_max))
    _apply_min_in_zone!(mod.depth_min, x, x_start, x_end, target_depth, taper_width)
    return mod
end
