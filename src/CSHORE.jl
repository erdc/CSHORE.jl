# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

"""
    CSHORE

Julia implementation of CSHORE (USACE) with native multi-fraction bed composition.

Design:
- Immutable inputs are held in `CshoreConfig`; mutable per-timestep state in
  `CshoreState`. Multi-fraction state (`bed_mass[j,layer,k]`) is present
  throughout; single-grain runs use `nf=1` as a degenerate case with no
  branching in transport or Exner.

Feature scope (in): waves, sediment, Exner, hardbottom (ISEDAV=±1), wet/dry swash,
overtopping/runup (IOVER=1), vegetation (IVEG=1,2,3), porous beds (IPERM=1).

Feature scope (out): dike grass erosion (IPROFL=2), clay (ICLAY=1), wire mesh
(ISEDAV=2), ponded ridge-runnel (IPOND=1), tidal gradient (ITIDE=1), Weibull
spectra (IWEIBULL=1), measured spectrum (IDISS=3).
"""
module CSHORE

using LinearAlgebra
using Printf
using SpecialFunctions: erfc
using Statistics
using Dates
using NCDatasets
using FFTW

# Configuration (immutable inputs)
include("config.jl")

# ThermalState struct (needs to precede state.jl so CshoreState.thermal
# can use the concrete Union{Nothing,ThermalState} type instead of Any).
# Methods that operate on ThermalState live in thermal/thermal.jl below.
include("thermal/thermal_types.jl")

# Mutable simulation state
include("state.jl")

# Pure utilities (integration, smoothing, interpolation)
include("utils/utilities.jl")

# Hydrodynamics — waves, currents, swash, overtopping
include("hydrodynamics/waves.jl")           # wave transformation, breaking, dispersion
include("hydrodynamics/hydro.jl")           # wave-averaged hydrodynamics, roller
include("hydrodynamics/wetdry.jl")          # wet/dry front, swash probability
include("hydrodynamics/overtopping.jl")     # overtopping flux, runup
include("hydrodynamics/transmission.jl")    # landward wave transmission (IWTRAN=1)
include("hydrodynamics/infragravity.jl")    # infragravity wave energy (IgConfig)

# Vegetation (wave/flow drag)
include("vegetation/vegetation.jl")         # vegetation drag and wave dissipation

# Groundwater + porous flow
include("groundwater/porous.jl")            # porous bed flow
include("groundwater/groundwater.jl")       # beach groundwater + surface moisture (GroundwaterConfig)

# Aeolian (wind-driven) transport
include("aeolian/wind_shear.jl")            # Kroy shear perturbation over topography + lee mask
include("aeolian/aeolian.jl")               # wind-driven sediment transport (DRT-style)

# Sediment + morphology (bed evolution lives here)
include("sediment/fractions.jl")            # fall velocity, critical Shields, hiding/exposure
include("sediment/composition.jl")          # bed layer mass tracking, grain sorting
include("sediment/transport.jl")            # sediment transport per fraction
include("sediment/exner.jl")                # Exner equation + limiter
include("sediment/clay_dike_erosion.jl")     # grassed dike erosion (IPROFL=2) + sand-over-clay (ICLAY=1); FORTRAN EROSON port
include("sediment/diffusion.jl")            # hillslope diffusion (non-wave mass wasting)
include("sediment/avalanche.jl")            # underwater angle-of-repose enforcement
include("sediment/cohesive.jl")             # minimal-viable mud / Partheniades-Krone

# Provenance tracking (optional; included before driver.jl so driver can use types)
include("utils/provenance.jl")

# Thermal / permafrost (optional — drives zb_hard via active-layer thickness)
include("thermal/thermal.jl")

# I/O
include("io/input.jl")                # read_infile, build_config
include("io/output.jl")               # write_obprof, write_osetup (FORTRAN ASCII stubs)
include("io/netcdf.jl")               # NetcdfWriter — CF-1.10 output
include("io/xbeach.jl")               # read_xbeach_params — XBeach params.txt → CshoreConfig
include("io/cshorejl.jl")             # read_cshorejl — native TOML .cshore format
include("io/read_config.jl")          # read_config — format-dispatched single entry point

# Top-level driver (must come after io/netcdf.jl so it can use open_netcdf)
include("driver.jl")

# BMI bindings (must come after driver.jl — needs step_bc_window!)
include("bmi.jl")

# Quasi-2D transect-grid framework (must come after bmi.jl — uses CshoreBMI)
include("q2d/transect_grid.jl")

# High-level presets (must come after driver.jl + io/input.jl).
# Split by category — order matters: nbs.jl uses helpers from profiles.jl,
# quick_run.jl uses everything else.
include("presets/profiles.jl")    # planar_beach, beach_dune_profile, etc.
include("presets/sediment.jl")    # sediment_fine_sand, …, spatially_varying_fractions
include("presets/nbs.jl")         # nbs_eelgrass / kelp / breakwater / revetment / thermosyphon …
include("presets/forcing.jl")     # constant_waves, storm_sequence, seasonal_*
include("presets/quick_run.jl")   # quick_run — single-call simulation entry point

# Ensemble / scenario sweep runner (parallel run_simulation! over a population)
include("ensemble.jl")

# Visualization (CairoMakie-based plotting utilities)
include("plotting.jl")

# Environmental data fetching (NDBC, NOAA CO-OPS)
include("data.jl")

# Command-line entry point for compiled binary (PackageCompiler create_app)
include("cli.jl")

# Public API — core types
export CshoreConfig, CshoreState
export OptionFlags, GridConfig, SedimentConfig, MultifractionConfig
export BoundaryTimeSeries, BathyInput
export VegetationInput, SwashConfig, OvertoppingConfig, DiffusionConfig
export DikeErosionInput, ClayInput, TidalInput, CurrentInput
export AeolianConfig, AeolianVegetationModel, ContourVegetation, DensityVegetation
export IgConfig, compute_ig_field!
export UndertowConfig, AsymmetryConfig, PhaseLagConfig, BailardConfig
export WaveNonlinearityConfig
export GroundwaterConfig, step_groundwater!
export CohesiveSedimentConfig, cohesive_step!, initial_cohesive_bed_mass
export VegetationSpecies, MultiSpeciesVegetation
export species_dune_grass, species_shrub, species_forb
export step_aeolian!, update_vegetation_density!
export WindShearConfig, kroy_shear_perturbation, lee_separation_mask, compute_wind_shear!
export ThermalConfig, ThermalBoundaryTimeSeries, ThermalState, SnowConfig, SnowSpatialModifier
export read_infile, write_outputs, make_sediment, build_config, submerged_sgm1
export read_xbeach_params    # XBeach params.txt reader
export read_cshorejl          # Native CSHORE.jl TOML format reader
export read_config            # Format-dispatched single entry point
export run_simulation!, initialize_state, step_bc_window!
# Ensemble runner
export EnsembleResult, ensemble_run, successes, failures
# Provenance tracking
export ProvenanceConfig, ProvenanceState
export init_provenance, step_provenance!, provenance_fractions
export TransportMethod, WatanabeTransport, SvrTransport, OriginalCshoreTransport, SizeAdaptiveTransport
# NetCDF output
export NetcdfWriter, open_netcdf, write_step!, maybe_write_step!, close_netcdf!
# FORTRAN-compatible ASCII outputs
export CshoreAsciiWriter, open_ascii_outputs, write_ascii_step!, close_ascii_outputs!
# IWTRAN=1 — landward wave transmission
export transmission!, crest_transmission_coefficient
# BMI wrapper
export CshoreBMI
# Quasi-2D transect-grid framework
export TransectGrid, LongshoreBoundary, WaveAngleMode
export WaveAngleFixed, WaveAngleShoreline, WaveAngleSnell
export step_transect_grid!, run_transect_grid!, transect_grid_output
export bmi_inject_dzb!, bmi_compute_qby_total!
# Presets — profile builders
export planar_beach, beach_dune_profile, arctic_bluff_profile, rocky_shore_profile
# Presets — sediment
export sediment_fine_sand, sediment_medium_sand, sediment_coarse_gravel
export sediment_arctic_mix, sediment_custom
export spatially_varying_fractions
# Presets — vegetation and structures
export nbs_eelgrass, nbs_kelp, nbs_beach_grass, nbs_dune_grass
export nbs_log_jam, nbs_breakwater, nbs_nourishment
export nbs_vegetation, nbs_vegetation_params,
       nbs_vegetation_manning_field, VegetatedManningField,
       NBS_VEGETATION_PARAMS, NBS_VEGETATION_MANNING_RANGE
export nbs_gravel_revetment
export nbs_thermosyphon, nbs_sod_insulation
export nbs_snow_fence, nbs_snow_clearing, nbs_insulating_drift
# Presets — wave forcing
export constant_waves, storm_sequence, seasonal_waves
export seasonal_temperature, seasonal_snow
# Presets — quick run
export quick_run
# Plotting
export plot_profile, plot_profile_evolution, plot_wave_field, plot_transport
export plot_hovmoller, plot_mass_balance, plot_swash, plot_thermal, save_figure
# Data fetching
export fetch_ndbc_realtime, fetch_ndbc_historical, fetch_tides
export ndbc_to_cshore_bc, tides_to_swl
# CLI entry point (used by PackageCompiler create_app)
export julia_main

end # module CSHORE
