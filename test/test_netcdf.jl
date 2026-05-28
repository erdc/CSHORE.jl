#==============================================================================
test_netcdf.jl — NetCDF writer + reader round-trip tests.
==============================================================================#

using Test
using CSHORE
using NCDatasets
using CSHORE: build_config, run_simulation!, OptionFlags, make_sediment,
              MultifractionConfig, open_netcdf, write_step!, close_netcdf!,
              initialize_state, apply_initial_bathymetry!, compute_derived_constants!

@testset "NetCDF output" begin
    @testset "single-grain round trip" begin
        tmp = joinpath(mktempdir(), "sg_roundtrip.nc")
        x = collect(0.0:1.0:200.0)
        z = range(-3.0, 1.0; length=length(x)) |> collect
        cfg = build_config(
            dx=1.0, bathymetry_x=x, bathymetry_z=copy(z), friction=0.002,
            timebc=[0.0, 3600.0], tpbc=[6.0, 6.0],
            hrmsbc=[0.5, 0.5], swlbc=[0.0, 0.0],
            options=OptionFlags(iprofl=0),
            sediment=make_sediment(d50=0.3e-3),
        )
        state = run_simulation!(cfg; outdir=dirname(tmp), outfile=basename(tmp))
        @test isfile(tmp)

        NCDataset(tmp, "r") do ds
            # Dimensions
            @test Int(ds.dim["x"]) == length(x)
            @test Int(ds.dim["time"]) == 2       # t=0 + one BC window
            @test Int(ds.dim["fraction"]) == 1
            @test Int(ds.dim["layer"]) == 3

            # CF conventions attribute
            @test ds.attrib["Conventions"] == "CF-1.10"
            @test occursin("CSHORE.jl", ds.attrib["source"])

            # Time coordinate — read with .var to skip CF decoding so we get
            # raw Float64 seconds instead of auto-parsed DateTime objects.
            @test ds["time"].var[1] == 0.0
            @test ds["time"].var[2] ≈ 3600.0 atol=1e-6
            time_units = ds["time"].attrib["units"]
            @test occursin("seconds since", time_units)

            # x coordinate matches bathymetry x
            @test ds["x"][:] ≈ x

            # Initial zb matches input
            @test ds["zb"][:, 1] ≈ z atol=1e-9
            # Fixed bed: final zb should equal initial
            @test ds["zb"][:, 2] ≈ z atol=1e-9

            # Wave transform results present
            @test all(ds["hrms"][:, 2] .≥ 0)
            @test any(ds["hrms"][:, 2] .> 0.1)   # some waves propagated
            @test ds["hrms"].attrib["units"] == "m"

            # No bed_mass variable in single-grain mode
            @test !haskey(ds, "bed_mass")
            @test !haskey(ds, "d50_surface")
        end
    end

    @testset "multifraction round trip" begin
        tmp = joinpath(mktempdir(), "mf_roundtrip.nc")
        x = collect(0.0:1.0:200.0)
        z = range(-3.0, 1.0; length=length(x)) |> collect
        mf = MultifractionConfig(
            grain_sizes=[0.15e-3, 0.25e-3, 0.50e-3],
            nlayers=3, layer_thickness=0.1, porosity=0.4,
            initial_fractions=[0.3, 0.5, 0.2],
            use_size_dependent_shields=true,
        )
        cfg = build_config(
            dx=1.0, bathymetry_x=x, bathymetry_z=copy(z), friction=0.002,
            timebc=collect(range(0.0, 86400.0; length=25)),  # hourly
            tpbc=fill(6.0, 25), hrmsbc=fill(0.5, 25), swlbc=fill(0.0, 25),
            options=OptionFlags(iprofl=1),
            sediment=make_sediment(d50=0.3e-3),
            multifraction=mf,
        )
        state = run_simulation!(cfg; outdir=dirname(tmp), outfile=basename(tmp))
        @test isfile(tmp)

        NCDataset(tmp, "r") do ds
            @test Int(ds.dim["fraction"]) == 3
            @test Int(ds.dim["time"]) == 25         # t=0 + 24 hourly slices
            @test haskey(ds, "bed_mass")
            @test haskey(ds, "d50_surface")

            # Grain-size coord
            @test ds["fraction"][:] ≈ [0.15e-3, 0.25e-3, 0.50e-3]

            # bed_mass dims: (x, layer, fraction, time) — nlayers=3 here
            @test size(ds["bed_mass"]) == (length(x), 3, 3, 25)

            # Initial active-layer bed_mass matches initial fractions
            # (uniform). ρs*(1-n)*thickness = 2650 * 0.6 * 0.1 = 159.0 kg/m² total per node
            total_init = 2650 * 0.6 * 0.1
            j = 100
            @test sum(ds["bed_mass"][j, 1, :, 1]) ≈ total_init atol=1e-6
            @test ds["bed_mass"][j, 1, 1, 1] / total_init ≈ 0.3 atol=1e-9
            @test ds["bed_mass"][j, 1, 2, 1] / total_init ≈ 0.5 atol=1e-9
            @test ds["bed_mass"][j, 1, 3, 1] / total_init ≈ 0.2 atol=1e-9

            # Subsurface layers have the same initial fractions
            @test sum(ds["bed_mass"][j, 2, :, 1]) ≈ total_init atol=1e-6
            @test sum(ds["bed_mass"][j, 3, :, 1]) ≈ total_init atol=1e-6

            # Initial d50_surface is the mass-weighted mean
            expected_d50 = 0.3 * 0.15e-3 + 0.5 * 0.25e-3 + 0.2 * 0.50e-3
            @test ds["d50_surface"][j, 1] ≈ expected_d50 atol=1e-9

            # d50_bulk should match initial d50_surface at t=0 (all layers identical)
            @test ds["d50_bulk"][j, 1] ≈ expected_d50 atol=1e-9

            # After 24 hours, zb should match in-memory state
            jmax = state.jmax[1]
            @test ds["zb"][1:jmax, end] ≈ state.zb[1:jmax, 1] atol=1e-9

            # qbx and qsx have per-fraction data
            @test size(ds["qbx"]) == (length(x), 3, 25)
            @test any(ds["qbx"][:, :, end] .!= 0)   # some transport occurred
        end
    end

    @testset "output_interval_s decimation" begin
        tmp = joinpath(mktempdir(), "decimated.nc")
        x = collect(0.0:1.0:100.0)
        z = range(-3.0, 1.0; length=length(x)) |> collect
        # 24 BC windows, but write only every 4 hours (every 4 windows) + initial.
        cfg = build_config(
            dx=1.0, bathymetry_x=x, bathymetry_z=copy(z), friction=0.002,
            timebc=collect(range(0.0, 86400.0; length=25)),
            tpbc=fill(6.0, 25), hrmsbc=fill(0.5, 25), swlbc=fill(0.0, 25),
            options=OptionFlags(iprofl=1),
            sediment=make_sediment(d50=0.3e-3),
        )
        run_simulation!(cfg; outdir=dirname(tmp), outfile=basename(tmp),
                         output_interval_s=4 * 3600.0)
        NCDataset(tmp, "r") do ds
            # Expected: t=0 + windows at 4h, 8h, 12h, 16h, 20h, 24h = 7 slices
            @test Int(ds.dim["time"]) == 7
            @test ds["time"].var[1] == 0.0
            @test ds["time"].var[end] ≈ 86400.0 atol=1.0
            # Differences between successive writes should be ≥ 4h
            ts = ds["time"].var[:]
            diffs = diff(ts)
            @test all(diffs .≥ 4 * 3600.0 - 1e-6)
        end
    end

    @testset "no output when outfile=nothing" begin
        # Default behavior: no NetCDF file written, state still returned
        x = collect(0.0:1.0:50.0)
        z = range(-2.0, 0.5; length=length(x)) |> collect
        cfg = build_config(
            dx=1.0, bathymetry_x=x, bathymetry_z=z, friction=0.002,
            timebc=[0.0, 3600.0], tpbc=[6.0, 6.0],
            hrmsbc=[0.5, 0.5], swlbc=[0.0, 0.0],
        )
        tmpdir = mktempdir()
        state = run_simulation!(cfg; outdir=tmpdir)  # no `outfile` kwarg
        @test state isa CshoreState
        # Directory should be empty (no NetCDF written)
        @test isempty(readdir(tmpdir))
    end
end
