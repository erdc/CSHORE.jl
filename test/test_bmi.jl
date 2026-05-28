#==============================================================================
test_bmi.jl — Basic Model Interface (BMI v2.0) compliance tests.
==============================================================================#

using Test
using CSHORE
using CSHORE: build_config, OptionFlags, make_sediment, MultifractionConfig,
              CshoreBMI
import BasicModelInterface as BMI

function _make_test_model(; nf::Int=1)
    x = collect(0.0:1.0:100.0)
    z = range(-3.0, 1.0; length=length(x)) |> collect
    mf = nf == 1 ? MultifractionConfig() :
        MultifractionConfig(
            grain_sizes=[0.15e-3, 0.25e-3, 0.50e-3],
            nlayers=3, layer_thickness=0.1, porosity=0.4,
            initial_fractions=[0.3, 0.5, 0.2],
            use_size_dependent_shields=true,
        )
    cfg = build_config(
        dx=1.0, bathymetry_x=x, bathymetry_z=copy(z), friction=0.002,
        timebc=collect(range(0.0, 3 * 3600.0; length=4)),  # 3 BC windows
        tpbc=fill(6.0, 4), hrmsbc=fill(0.5, 4), swlbc=fill(0.0, 4),
        options=OptionFlags(iprofl=1),
        sediment=make_sediment(d50=0.3e-3),
        multifraction=mf,
    )
    return CshoreBMI(cfg)
end

@testset "BMI v2.0 compliance" begin
    @testset "component info" begin
        m = _make_test_model()
        @test BMI.get_component_name(m) == "CSHORE.jl"
    end

    @testset "variable info" begin
        m = _make_test_model(nf=3)
        input = BMI.get_input_var_names(m)
        output = BMI.get_output_var_names(m)
        @test input isa Vector{<:AbstractString}
        @test output isa Vector{<:AbstractString}
        @test "zb" in output
        @test "hrms" in output
        @test "bed_mass" in output
        @test BMI.get_input_item_count(m) == length(input)
        @test BMI.get_output_item_count(m) == length(output)

        # units
        @test BMI.get_var_units(m, "zb") == "m"
        @test BMI.get_var_units(m, "hrms") == "m"
        @test BMI.get_var_units(m, "umean") == "m s-1"
        @test BMI.get_var_units(m, "qbx") == "m2 s-1"
        @test BMI.get_var_units(m, "bed_mass") == "kg m-2"

        # types
        @test BMI.get_var_type(m, "zb") == "float64"
        @test BMI.get_var_itemsize(m, "zb") == 8
        @test BMI.get_var_location(m, "zb") == "node"
        @test BMI.get_var_grid(m, "zb") == 0
        @test BMI.get_var_grid(m, "bed_mass") == 0

        # nbytes: zb is (nn,) so nn*8, bed_mass is (nn, nf)*8
        nn = m.state.jmax[1]
        @test BMI.get_var_nbytes(m, "zb") == nn * 8
        @test BMI.get_var_nbytes(m, "bed_mass") == nn * 3 * 8

        # Unknown variable errors cleanly
        @test_throws ErrorException BMI.get_var_units(m, "nonexistent")
        @test_throws ErrorException BMI.get_var_nbytes(m, "nonexistent")
    end

    @testset "time" begin
        m = _make_test_model()
        @test BMI.get_start_time(m) == 0.0
        @test BMI.get_end_time(m) == 3 * 3600.0
        @test BMI.get_current_time(m) == 0.0
        @test BMI.get_time_units(m) == "s"
        # get_time_step is 0 before any update (state.delt = 0 initially)
        @test BMI.get_time_step(m) isa Real
    end

    @testset "grid metadata" begin
        m = _make_test_model()
        @test BMI.get_grid_rank(m, 0) == 1
        @test BMI.get_grid_size(m, 0) == m.state.jmax[1]
        @test BMI.get_grid_type(m, 0) == "uniform_rectilinear"

        shape = zeros(Int, 1)
        BMI.get_grid_shape(m, 0, shape)
        @test shape[1] == m.state.jmax[1]

        spacing = zeros(1)
        BMI.get_grid_spacing(m, 0, spacing)
        @test spacing[1] ≈ 1.0

        origin = zeros(1)
        BMI.get_grid_origin(m, 0, origin)
        @test origin[1] ≈ 0.0

        xcoord = zeros(m.state.jmax[1])
        BMI.get_grid_x(m, 0, xcoord)
        @test xcoord[1] ≈ 0.0
        @test xcoord[end] ≈ 100.0

        # y, z are zero-length / no-op for 1D
        empty_y = Float64[]
        @test BMI.get_grid_y(m, 0, empty_y) === empty_y

        # Bad grid id errors
        @test_throws ErrorException BMI.get_grid_rank(m, 1)

        # Unstructured methods all error (CSHORE is structured)
        @test_throws ErrorException BMI.get_grid_node_count(m, 0)
    end

    @testset "value access (get / set)" begin
        m = _make_test_model()
        nn = m.state.jmax[1]

        # get_value into a destination buffer
        dest = zeros(nn)
        BMI.get_value(m, "zb", dest)
        @test dest[1] == m.state.zb[1, 1]
        @test dest[end] == m.state.zb[nn, 1]

        # get_value_ptr returns a view — zero-copy
        zb_view = BMI.get_value_ptr(m, "zb")
        @test zb_view isa AbstractArray
        @test length(zb_view) == nn
        # Mutate the view → state changes
        original = m.state.zb[50, 1]
        zb_view[50] = original + 0.001
        @test m.state.zb[50, 1] ≈ original + 0.001

        # set_value replaces the field
        new_zb = collect(Float64, range(-2.5, 0.5; length=nn))
        BMI.set_value(m, "zb", new_zb)
        @test view(m.state.zb, 1:nn, 1) ≈ new_zb

        # get_value_at_indices / set_value_at_indices
        inds = [10, 30, 60]
        buf = zeros(3)
        BMI.get_value_at_indices(m, "zb", buf, inds)
        @test buf ≈ new_zb[inds]

        new_vals = [-1.1, -0.7, -0.2]
        BMI.set_value_at_indices(m, "zb", inds, new_vals)
        buf2 = zeros(3)
        BMI.get_value_at_indices(m, "zb", buf2, inds)
        @test buf2 ≈ new_vals
    end

    @testset "lifecycle: update / update_until / finalize" begin
        m = _make_test_model()
        @test m.itime == 0
        @test BMI.get_current_time(m) == 0.0

        # One update advances one BC window (1 hour in our test config)
        BMI.update(m)
        @test m.itime == 1
        @test BMI.get_current_time(m) ≈ 3600.0 atol=1e-6

        # Another update
        BMI.update(m)
        @test m.itime == 2
        @test BMI.get_current_time(m) ≈ 7200.0 atol=1e-6

        # update_until
        BMI.update_until(m, 3 * 3600.0)
        @test m.itime == 3
        @test BMI.get_current_time(m) ≈ 3 * 3600.0 atol=1e-6

        # No-op finalize
        BMI.finalize(m)
        @test m.writer === nothing   # still nothing (never had one)
    end

    @testset "BMI.initialize error message is actionable" begin
        e = try
            BMI.initialize(CshoreBMI, "nonexistent.toml")
            nothing
        catch err
            err
        end
        @test e isa ErrorException
        @test occursin("CshoreBMI(cfg)", e.msg)
    end
end
