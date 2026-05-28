# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
presets.jl — High-level helpers for common CSHORE setups.

Reduces the typical 30+ lines of build_config boilerplate to 3-5 lines.
Provides profile builders, sediment presets, vegetation/structure presets,
wave forcing builders, and a single-function `quick_run` entry point.
==============================================================================#
# ============================================================================
# Quick-run — single-function simulation entry point
# ============================================================================

"""
    quick_run(; profile=:beach_dune, sediment=:medium_sand, waves=:constant,
               nbs=Symbol[], thermal=false, duration_days=1, kwargs...) -> CshoreState

Run a complete CSHORE simulation in one function call.

# Profiles
`:planar_beach`, `:beach_dune`, `:arctic_bluff`, `:rocky_shore`

# Sediment
`:fine_sand`, `:medium_sand`, `:coarse_gravel`, `:arctic_mix`

# Waves
`:constant`, `:storm`, `:seasonal`

# Vegetation / structures (vector of symbols, combined)
`:eelgrass`, `:kelp`, `:beach_grass`, `:dune_grass`, `:log_jam`, `:breakwater`, `:nourishment`,
`:gravel_revetment`, `:cobble_revetment`, `:rubble_mound`,
`:thermosyphon`, `:sod_insulation`, `:snow`

# Example
```julia
state = quick_run(profile=:arctic_bluff, sediment=:arctic_mix,
                   waves=:storm, nbs=[:eelgrass], thermal=true, duration_days=3)
```
"""
function quick_run(;
    profile::Symbol=:beach_dune,
    sediment::Symbol=:medium_sand,
    waves::Symbol=:constant,
    nbs::Vector{Symbol}=Symbol[],
    thermal::Bool=false,
    duration_days::Real=1,
    outfile::Union{Nothing,String}=nothing,
    kwargs...,
)
    # 1. Profile
    x, z = if profile == :planar_beach
        planar_beach()
    elseif profile == :beach_dune
        beach_dune_profile()
    elseif profile == :arctic_bluff
        arctic_bluff_profile()
    elseif profile == :rocky_shore
        rocky_shore_profile()
    else
        error("Unknown profile: $profile. Use :planar_beach, :beach_dune, :arctic_bluff, :rocky_shore")
    end

    # 2. Sediment
    mf = if sediment == :fine_sand;       sediment_fine_sand()
    elseif sediment == :medium_sand;      sediment_medium_sand()
    elseif sediment == :coarse_gravel;    sediment_coarse_gravel()
    elseif sediment == :arctic_mix;       sediment_arctic_mix()
    else error("Unknown sediment: $sediment")
    end

    # 3. Waves
    wf = if waves == :constant
        constant_waves(; duration_days=duration_days)
    elseif waves == :storm
        storm_sequence(; duration_hours=duration_days * 24.0)
    elseif waves == :seasonal
        seasonal_waves(; duration_years=duration_days / 365.25)
    else error("Unknown waves: $waves")
    end

    # 4. Vegetation / structures
    veg = nothing
    hardbottom_z = nothing
    iveg = 0
    idiss = 0
    iperm = 0
    porous_cfg = nothing
    porous_z_vec = nothing
    Q_thermo_vec = Float64[]
    R_insul_vec = Float64[]
    snow_cfg = nothing
    snow_depth_bc = nothing

    for n in nbs
        if n == :eelgrass
            veg = nbs_eelgrass(x, z); iveg = 3; idiss = 1
        elseif n == :kelp
            veg = nbs_kelp(x, z); iveg = 3; idiss = 1
        elseif n == :beach_grass
            shore_x = x[findfirst(zi -> zi ≥ 0.0, z)]
            bluff_x = x[findfirst(zi -> zi ≥ 2.0, z)]
            veg = nbs_beach_grass(x, z; x_start=shore_x, x_end=bluff_x); iveg = 1
        elseif n == :dune_grass
            crest_x = x[argmax(z)]
            veg = nbs_dune_grass(x, z; x_start=crest_x - 20.0, x_end=crest_x + 10.0); iveg = 1
        elseif n == :log_jam
            shore_x = x[findfirst(zi -> zi ≥ 0.0, z)]
            veg = nbs_log_jam(x, z; x_center=shore_x + 20.0); iveg = 1
        elseif n == :breakwater
            x, z, hardbottom_z = nbs_breakwater(x, z)
        elseif n == :nourishment
            x, z = nbs_nourishment(x, z)
        elseif n == :thermosyphon
            # Applied to the bluff face (above +2 m)
            bluff_start = x[findfirst(zi -> zi ≥ 2.0, z)]
            bluff_end = x[findlast(zi -> zi ≥ 2.0, z)]
            Q_thermo_vec = nbs_thermosyphon(x; x_start=bluff_start, x_end=bluff_end)
            thermal = true  # force thermal on
        elseif n == :sod_insulation
            # Applied to the bluff face (above +2 m)
            bluff_start = x[findfirst(zi -> zi ≥ 2.0, z)]
            bluff_end = x[findlast(zi -> zi ≥ 2.0, z)]
            R_insul_vec = nbs_sod_insulation(x; x_start=bluff_start, x_end=bluff_end)
            thermal = true  # force thermal on
        elseif n == :snow
            # Degree-day snow model with default Arctic settings
            snow_cfg = SnowConfig()
            thermal = true  # force thermal on
        elseif n == :gravel_revetment
            porous_z_vec, porous_cfg = nbs_gravel_revetment(x, z)
        elseif n == :cobble_revetment
            porous_z_vec, porous_cfg = nbs_gravel_revetment(x, z; Dn50=0.05, thickness=0.8)
        elseif n == :rubble_mound
            porous_z_vec, porous_cfg = nbs_gravel_revetment(x, z; Dn50=0.10, thickness=1.5)
        else
            @warn "Unknown NbS: $n — skipping"
        end
    end

    # 5. Thermal
    thermal_cfg = nothing
    thermal_time = nothing
    T_air = nothing
    T_water = nothing
    if thermal
        thermal_cfg = ThermalConfig(R_insulation=R_insul_vec, Q_thermosyphon=Q_thermo_vec)
        tf = seasonal_temperature(; duration_years=duration_days / 365.25)
        # Align thermal time with wave time
        thermal_time = wf.timebc
        T_air = [tf.T_air[clamp(round(Int, t / (tf.thermal_time[end] / length(tf.T_air))) + 1,
                                1, length(tf.T_air))] for t in wf.timebc]
        T_water = [tf.T_water[clamp(round(Int, t / (tf.thermal_time[end] / length(tf.T_water))) + 1,
                                    1, length(tf.T_water))] for t in wf.timebc]
    end

    # 6. Build config and run
    # Auto-set iperm if a porous revetment preset was selected
    if porous_cfg !== nothing; iperm = 1; end
    opts = OptionFlags(iprofl=1, iveg=iveg, idiss=idiss, iperm=iperm,
                        isedav=(hardbottom_z !== nothing || porous_z_vec !== nothing ? 1 : 0))

    cfg = build_config(
        dx=1.0, bathymetry_x=x, bathymetry_z=z, friction=0.002,
        timebc=wf.timebc, tpbc=wf.tpbc, hrmsbc=wf.hrmsbc,
        swlbc=wf.swlbc, wangbc=wf.wangbc,
        options=opts, sediment=make_sediment(),
        multifraction=mf,
        vegetation=veg,
        hardbottom_z=hardbottom_z,
        porous=porous_cfg,
        thermal=thermal_cfg,
        snow=snow_cfg,
        thermal_time=thermal_time,
        T_air=T_air,
        T_water=T_water,
        snow_depth=snow_depth_bc,
    )

    return run_simulation!(cfg; outfile=outfile)
end
