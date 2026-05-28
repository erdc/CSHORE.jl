#!/usr/bin/env julia
"""
coastal_vegetation_demo.jl — Demonstrates two parameterizations for coastal
nature-based features in CSHORE.jl:

  (1) Cd-based vegetation drag (Mendez & Losada 2004) via the
      `nbs_vegetation` preset and `NBS_VEGETATION_PARAMS` lookup.

  (2) Manning's n bottom friction (depth-averaged, Chezy-style quadratic
      drag) via `build_config(manning=...)` — dynamically recomputes the
      friction factor fb2 = g·n²/h^(1/3) each timestep from the live water
      depth (options.ifriction_spatial = 2).

Runs each vegetation category on a common beach-dune profile with
identical wave forcing, reports the wave-height attenuation, and compares
both parameterization paths for a representative species.
"""

using CSHORE
using CSHORE: build_config, run_simulation!, make_sediment, OptionFlags,
              VegetationInput, beach_dune_profile,
              nbs_vegetation, nbs_vegetation_params,
              nbs_vegetation_manning_field, VegetatedManningField,
              NBS_VEGETATION_PARAMS, NBS_VEGETATION_MANNING_RANGE,
              constant_waves

using Printf

# ----------------------------------------------------------------------------
# 1. Lookup table
# ----------------------------------------------------------------------------
println("Coastal vegetation parameter lookup")
println("  (category, region) → (Cd, N stems/m², b m, h m)  [Manning n range]")
for ((cat, region), p) in sort(collect(NBS_VEGETATION_PARAMS);
                                by = x -> (string(x[1][1]), string(x[1][2])))
    nrange = get(NBS_VEGETATION_MANNING_RANGE, (cat, region), (NaN, NaN))
    @printf("  %-22s %-10s Cd=%-4.2f  N=%-7.1f  b=%-6.4f  h=%-5.2f   n=[%.3f, %.3f]\n",
            cat, region, p.cd, p.n, p.dia, p.ht, nrange[1], nrange[2])
end
println()

# ----------------------------------------------------------------------------
# 2. Common profile + forcing
# ----------------------------------------------------------------------------
x, z = beach_dune_profile(; offshore_depth=4.0, offshore_dist=300.0,
                            beach_slope=0.04, beach_width=30.0,
                            dune_height=3.0, dune_width=20.0,
                            upland_width=40.0, dx=1.0)
dune_x = x[argmax(z)]

const HRMS = 1.2
const TP   = 8.0
const SWL  = 0.8
bc = constant_waves(; duration_days=1.0, dt_hours=1.0, hrms=HRMS, tp=TP, swl=SWL)

# ----------------------------------------------------------------------------
# 3. Vegetation scenarios — Cd-based drag (Mendez & Losada)
# ----------------------------------------------------------------------------
veg_scenarios = [
    (:bare,             ()  -> nothing),
    (:mangrove_gulf,    ()  -> nbs_vegetation(:mangrove,        x, z; region=:gulf)),
    (:low_marsh_tall,   ()  -> nbs_vegetation(:low_marsh_tall,  x, z; region=:east_gulf)),
    (:low_marsh_short,  ()  -> nbs_vegetation(:low_marsh_short, x, z; region=:east_gulf)),
    (:low_marsh_west,   ()  -> nbs_vegetation(:low_marsh,       x, z; region=:west)),
    (:high_marsh_east,  ()  -> nbs_vegetation(:high_marsh,      x, z; region=:east_gulf)),
    (:high_marsh_west,  ()  -> nbs_vegetation(:high_marsh,      x, z; region=:west)),
    (:phragmites,       ()  -> nbs_vegetation(:phragmites,      x, z; region=:east_gulf)),
    (:marsh_shrub,      ()  -> nbs_vegetation(:marsh_shrub,     x, z; region=:east_gulf)),
    (:dune_grass_east,  ()  -> nbs_vegetation(:dune_grass,      x, z; region=:east_gulf,
                                              x_start=dune_x - 15.0, x_end=dune_x + 10.0)),
    (:dune_grass_west,  ()  -> nbs_vegetation(:dune_grass,      x, z; region=:west,
                                              x_start=dune_x - 15.0, x_end=dune_x + 10.0)),
    (:dune_shrub_west,  ()  -> nbs_vegetation(:dune_shrub,      x, z; region=:west,
                                              x_start=dune_x - 5.0,  x_end=dune_x + 15.0)),
    (:beach_forbs_gulf, ()  -> nbs_vegetation(:beach_forbs,     x, z; region=:gulf,
                                              x_start=dune_x - 25.0, x_end=dune_x - 5.0)),
]

println("Cd-based vegetation drag (Mendez & Losada 2004)")
println("Scenario               Cd   Nstems   ht(m)   Hrms_off  Hrms_shore  Δ(%)")
println("--------------------------------------------------------------------------")
for (name, build_veg) in veg_scenarios
    veg = build_veg()
    iveg = veg === nothing ? 0 : 3
    cfg = build_config(
        dx = 1.0, bathymetry_x = x, bathymetry_z = z,
        sediment = make_sediment(d50=3e-4),
        timebc = bc.timebc, hrmsbc = bc.hrmsbc, tpbc = bc.tpbc, swlbc = bc.swlbc,
        vegetation = veg,
        options = OptionFlags(iprofl=1, iveg=iveg, idiss = iveg == 3 ? 1 : 0),
    )
    state = run_simulation!(cfg)
    wet = findall(>(0.0), state.hrms)
    shore_idx = isempty(wet) ? lastindex(state.hrms) : last(wet)
    hrms_off, hrms_shore = state.hrms[1], state.hrms[shore_idx]
    atten = 100.0 * (hrms_off - hrms_shore) / hrms_off

    if veg === nothing
        @printf("%-22s   -      -        -     %7.3f   %7.3f   %5.1f\n",
                name, hrms_off, hrms_shore, atten)
    else
        ks = findall(>(0.0), view(veg.vegn, :, 1))
        n_rep = isempty(ks) ? 0.0 : veg.vegn[ks[1], 1]
        h_rep = isempty(ks) ? 0.0 : veg.vegd[ks[1], 1]
        @printf("%-22s %4.2f  %6.1f   %5.2f  %7.3f   %7.3f   %5.1f\n",
                name, veg.vegcd, n_rep, h_rep, hrms_off, hrms_shore, atten)
    end
end
println()

# ----------------------------------------------------------------------------
# 4. Manning's n bottom friction (alternative path)
# ----------------------------------------------------------------------------
# Demonstrate a spatially-varying Manning's n: open water (low n) seaward,
# rough vegetated marsh (high n) shoreward of the surf zone.
np = length(x)
manning = fill(0.020, np)                  # bare sand baseline
for i in 1:np
    if z[i] > -0.5 && z[i] < 0.5           # marsh band straddling MSL
        manning[i] = 0.10                  # representative S. alterniflora
    end
end

println("Manning's n bottom friction (ifriction_spatial=2)")
println("Each scenario is the Manning-only path: vegetation effect is encoded")
println("in the spatially-varying Manning field; no VegetationInput passed.")
println()
println("Scenario               n_range          Hrms_off  Hrms_shore  Δ(%)")
println("---------------------------------------------------------------------")
n_scenarios = [
    (:bare_sand,         () -> fill(0.020, np)),
    (:low_marsh_tall,    () -> nbs_vegetation_manning_field(:low_marsh_tall, x, z;
                                  region=:east_gulf, n_bare=0.02, taper_width=0.3)),
    (:phragmites,        () -> nbs_vegetation_manning_field(:phragmites, x, z;
                                  region=:east_gulf, n_bare=0.02, taper_width=0.3)),
    (:dune_grass_east,   () -> nbs_vegetation_manning_field(:dune_grass, x, z;
                                  region=:east_gulf,
                                  x_start=dune_x - 15.0, x_end=dune_x + 10.0,
                                  n_bare=0.02)),
    (:marsh_plus_dune,   () -> begin
        f = nbs_vegetation_manning_field(:low_marsh_tall, x, z;
                                          region=:east_gulf, n_bare=0.02, taper_width=0.3)
        nbs_vegetation_manning_field(:dune_grass, x, z; region=:east_gulf,
                                      x_start=dune_x - 15.0, x_end=dune_x + 10.0,
                                      n_bare=0.02, base=f)
    end),
]
for (name, build_n) in n_scenarios
    n_field = build_n()
    n_values = n_field isa VegetatedManningField ? n_field.values : n_field
    cfg = build_config(
        dx = 1.0, bathymetry_x = x, bathymetry_z = z,
        sediment = make_sediment(d50=3e-4),
        timebc = bc.timebc, hrmsbc = bc.hrmsbc, tpbc = bc.tpbc, swlbc = bc.swlbc,
        manning = n_field,
        options = OptionFlags(iprofl=1),
    )
    state = run_simulation!(cfg)
    wet = findall(>(0.0), state.hrms)
    shore_idx = isempty(wet) ? lastindex(state.hrms) : last(wet)
    hrms_off, hrms_shore = state.hrms[1], state.hrms[shore_idx]
    atten = 100.0 * (hrms_off - hrms_shore) / hrms_off
    @printf("%-22s [%.3f, %.3f]  %7.3f   %7.3f   %5.1f\n",
            name, extrema(n_values)..., hrms_off, hrms_shore, atten)
end

# Demonstrate the double-counting guard
println()
println("Double-counting guard demonstration:")
mf = nbs_vegetation_manning_field(:low_marsh_tall, x, z; region=:east_gulf,
                                   n_bare=0.02, taper_width=0.3)
veg = nbs_vegetation(:low_marsh_tall, x, z; region=:east_gulf)
try
    build_config(dx=1.0, bathymetry_x=x, bathymetry_z=z,
        sediment=make_sediment(d50=3e-4),
        timebc=bc.timebc, hrmsbc=bc.hrmsbc, tpbc=bc.tpbc, swlbc=bc.swlbc,
        manning=mf, vegetation=veg,
        options=OptionFlags(iprofl=1, iveg=3, idiss=1))
    println("  ✗ Guard failed: should have rejected")
catch e
    @printf("  ✓ Rejected as expected:\n    %s\n", first(sprint(showerror, e), 220))
end
