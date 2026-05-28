@testset "cohesive (mud) sediment" begin
    using CSHORE: build_config, OptionFlags, MultifractionConfig, make_sediment,
                  initialize_state, apply_initial_bathymetry!, run_simulation!,
                  CohesiveSedimentConfig, cohesive_step!, GRAV

    # Small grid; sloping bathymetry crossing SWL so some nodes are wet.
    n = 20
    x = collect(0.0:1.0:Float64(n - 1))
    z = collect(range(-3.0, 1.0; length=n))
    mf = MultifractionConfig(
        grain_sizes=[0.25e-3], initial_fractions=[1.0],
        nlayers=3, layer_thickness=0.1, porosity=0.4,
    )

    base_kwargs = (
        dx=1.0, bathymetry_x=x, bathymetry_z=copy(z), friction=0.002,
        timebc=[0.0, 60.0, 120.0], tpbc=[8.0, 8.0, 8.0],
        hrmsbc=[1.0, 1.0, 1.0], swlbc=[0.0, 0.0, 0.0],
        options=OptionFlags(iprofl=0),
        sediment=make_sediment(d50=0.25e-3),
        multifraction=mf,
    )

    @testset "config defaults are physical" begin
        c = CohesiveSedimentConfig()
        @test c.settling_velocity > 0
        @test 0 < c.tau_cd <= c.tau_ce
        @test c.M > 0
        @test c.initial_bed_mass >= 0
        @test c.rho_water > 0
    end

    @testset "initialize_state seeds bed mass" begin
        coh = CohesiveSedimentConfig(initial_bed_mass=50.0)
        cfg = build_config(; base_kwargs..., cohesive=coh)
        state = initialize_state(cfg)
        @test length(state.cohesive_bed_mass) > 0
        @test all(state.cohesive_bed_mass .≈ 50.0)
        @test all(state.cohesive_concentration .≈ 0.0)
    end

    @testset "no-op when config.cohesive is nothing" begin
        cfg = build_config(; base_kwargs...)
        state = initialize_state(cfg)
        @test all(state.cohesive_bed_mass .≈ 0.0)
        cohesive_step!(state, cfg, 1, 60.0)   # should not error
        @test all(state.cohesive_bed_mass .≈ 0.0)
        @test all(state.cohesive_concentration .≈ 0.0)
    end

    @testset "erosion: high shear erodes mud" begin
        # τ_ce very small so any shear erodes.
        coh = CohesiveSedimentConfig(tau_ce=0.001, tau_cd=0.0005,
                                     M=1e-3, initial_bed_mass=10.0,
                                     settling_velocity=1e-10)  # essentially no deposition
        cfg = build_config(; base_kwargs..., cohesive=coh)
        state = initialize_state(cfg)
        apply_initial_bathymetry!(state, cfg)
        # Force a non-zero bed shear stress at every wet node.
        # tbxsta is the normalized cross-shore bed shear (length units).
        # τ_b = ρ_w · g · |tbxsta|, so tbxsta = 0.1 gives τ_b ≈ 1005 Pa
        # — well above τ_ce.
        for j in 1:length(state.tbxsta)
            state.tbxsta[j] = 0.1
            state.h[j] = max(state.h[j], 1.0)  # ensure wet
        end
        bed_before = sum(state.cohesive_bed_mass)
        conc_before = sum(state.cohesive_concentration)
        cohesive_step!(state, cfg, 1, 10.0)
        bed_after = sum(state.cohesive_bed_mass)
        conc_after = sum(state.cohesive_concentration)
        @test bed_after < bed_before          # bed eroded
        @test conc_after > conc_before        # suspended pool grew
    end

    @testset "deposition: low shear deposits suspended" begin
        coh = CohesiveSedimentConfig(tau_ce=10.0, tau_cd=5.0,
                                     M=1e-3, initial_bed_mass=0.0,
                                     settling_velocity=1e-3)
        cfg = build_config(; base_kwargs..., cohesive=coh)
        state = initialize_state(cfg)
        apply_initial_bathymetry!(state, cfg)
        # Pre-seed suspended concentration; zero shear.
        for j in 1:length(state.tbxsta)
            state.tbxsta[j] = 0.0
            state.h[j] = max(state.h[j], 1.0)
            state.cohesive_concentration[j] = 1.0  # 1 kg/m³
        end
        conc_before = sum(state.cohesive_concentration)
        bed_before = sum(state.cohesive_bed_mass)
        cohesive_step!(state, cfg, 1, 60.0)
        @test sum(state.cohesive_concentration) < conc_before  # settled
        @test sum(state.cohesive_bed_mass)      > bed_before    # bed grew
    end

    @testset "mass conservation under erosion/deposition cycle" begin
        coh = CohesiveSedimentConfig(tau_ce=0.001, tau_cd=0.0005,
                                     M=1e-3, initial_bed_mass=5.0,
                                     settling_velocity=1e-3)
        cfg = build_config(; base_kwargs..., cohesive=coh)
        state = initialize_state(cfg)
        apply_initial_bathymetry!(state, cfg)
        for j in 1:length(state.tbxsta)
            state.h[j] = 1.0
        end

        # Total mass (bed + suspended·h) should be conserved across
        # an erosion-then-deposition cycle (each step explicit Euler,
        # so allow modest tolerance).
        function total_mass(s)
            tot = 0.0
            for j in 1:length(s.cohesive_bed_mass)
                tot += s.cohesive_bed_mass[j] + s.cohesive_concentration[j] * s.h[j]
            end
            tot
        end
        m0 = total_mass(state)

        # Erode
        for j in 1:length(state.tbxsta); state.tbxsta[j] = 0.1; end
        cohesive_step!(state, cfg, 1, 5.0)
        m1 = total_mass(state)
        @test isapprox(m1, m0; rtol=1e-9)

        # Deposit
        for j in 1:length(state.tbxsta); state.tbxsta[j] = 0.0; end
        cohesive_step!(state, cfg, 1, 5.0)
        m2 = total_mass(state)
        @test isapprox(m2, m0; rtol=1e-9)
    end

    @testset "bed cannot go negative" begin
        coh = CohesiveSedimentConfig(tau_ce=0.001, tau_cd=0.0005,
                                     M=1.0, initial_bed_mass=0.1)
        cfg = build_config(; base_kwargs..., cohesive=coh)
        state = initialize_state(cfg)
        apply_initial_bathymetry!(state, cfg)
        for j in 1:length(state.tbxsta)
            state.tbxsta[j] = 1.0
            state.h[j] = max(state.h[j], 1.0)
        end
        # Erode harder than bed can supply.
        for _ in 1:10
            cohesive_step!(state, cfg, 1, 100.0)
        end
        @test all(state.cohesive_bed_mass .>= 0.0)
        @test all(state.cohesive_concentration .>= 0.0)
    end

    @testset "integrates into run_simulation!" begin
        coh = CohesiveSedimentConfig(initial_bed_mass=10.0)
        cfg = build_config(; base_kwargs..., cohesive=coh)
        s = run_simulation!(cfg)
        @test s.time > 0
        @test length(s.cohesive_bed_mass) == length(s.cohesive_concentration)
    end
end
