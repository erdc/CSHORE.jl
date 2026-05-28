@testset "ensemble" begin
    # Tiny fixed-bed config (iprofl=0) so each ensemble member runs fast.
    n_nodes = 20
    x = collect(0.0:1.0:Float64(n_nodes - 1))
    z = collect(range(-3.0, 1.0; length=n_nodes))
    base_kwargs = (
        dx=1.0, bathymetry_x=x, bathymetry_z=copy(z), friction=0.002,
        timebc=[0.0, 60.0, 120.0], tpbc=[8.0, 8.0, 8.0],
        hrmsbc=[1.0, 1.0, 1.0], swlbc=[0.0, 0.0, 0.0],
        options=OptionFlags(iprofl=0),
        sediment=make_sediment(d50=0.25e-3),
    )
    base = build_config(; base_kwargs...)

    hrmss = [0.5, 1.0, 1.5]
    patches = [cfg -> begin
                    cfg.boundary.hrmsbc .= h
                    cfg
                end
               for h in hrmss]

    @testset "patch form" begin
        r = ensemble_run(base, patches; threaded=false)
        @test length(r) == 3
        @test all(r.succeeded)
        @test count(r.succeeded) == 3
        @test isempty(failures(r))
        @test length(successes(r)) == 3
        for s in successes(r)
            @test s.time > 0
        end
        @test r.runtimes_s isa Vector{Float64}
        @test all(t -> t >= 0, r.runtimes_s)
    end

    @testset "config-vector form" begin
        configs = [patches[i](deepcopy(base)) for i in 1:length(patches)]
        r = ensemble_run(configs; threaded=false)
        @test length(r) == 3
        @test all(r.succeeded)
    end

    @testset "patches don't share state" begin
        # Each member should see its own HRMS, not the last patch's value.
        r = ensemble_run(base, patches; threaded=false)
        for (i, h) in enumerate(hrmss)
            @test r.configs[i].boundary.hrmsbc[1] ≈ h
        end
        # Base config left untouched.
        @test base.boundary.hrmsbc[1] ≈ 1.0
    end

    @testset "failure captured, not raised" begin
        bad_patches = [
            cfg -> cfg,
            cfg -> error("synthetic ensemble member failure"),
            cfg -> cfg,
        ]
        r = ensemble_run(base, bad_patches; threaded=false)
        @test length(r) == 3
        @test r.succeeded == BitVector([true, false, true])
        fs = failures(r)
        @test length(fs) == 1
        @test fs[1][1] == 2
        @test fs[1][2] isa Exception
    end

    @testset "empty input" begin
        r = ensemble_run(CshoreConfig[])
        @test length(r) == 0
    end

    @testset "filename_pattern validation" begin
        @test_throws ArgumentError ensemble_run([base]; filename_pattern="no_index.nc")
    end

    @testset "threaded execution" begin
        # Smoke test: just confirms threaded=true doesn't crash and yields
        # the same success pattern. Cannot assert speedup with 3 tiny members.
        r = ensemble_run(base, patches; threaded=true)
        @test all(r.succeeded)
    end
end
