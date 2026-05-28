# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
presets.jl — High-level helpers for common CSHORE setups.

Reduces the typical 30+ lines of build_config boilerplate to 3-5 lines.
Provides profile builders, sediment presets, vegetation/structure presets,
wave forcing builders, and a single-function `quick_run` entry point.
==============================================================================#
# ============================================================================
# Wave forcing builders — return (timebc, hrmsbc, tpbc, swlbc [, wangbc])
# ============================================================================

"""
    constant_waves(; hrms=1.0, tp=8.0, swl=0.0, angle=0.0,
                     duration_days=1, dt_hours=1.0) -> NamedTuple

Constant wave forcing over the specified duration.
Returns a NamedTuple with `timebc`, `hrmsbc`, `tpbc`, `swlbc`, `wangbc`.
"""
function constant_waves(; hrms::Float64=1.0, tp::Float64=8.0, swl::Float64=0.0,
                          angle::Float64=0.0, duration_days::Real=1, dt_hours::Float64=1.0)
    ntimes = round(Int, duration_days * 24 / dt_hours) + 1
    timebc = collect(range(0.0, duration_days * 86400.0; length=ntimes))
    return (; timebc, hrmsbc=fill(hrms, ntimes), tpbc=fill(tp, ntimes),
              swlbc=fill(swl, ntimes), wangbc=fill(angle, ntimes))
end

"""
    storm_sequence(; hrms_peak=2.0, tp_peak=10.0, surge_peak=1.5, angle=0.0,
                     ramp_hours=6, duration_hours=24, dt_hours=0.5) -> NamedTuple

Synthetic storm: ramp-up → peak → ramp-down in wave height and surge.
"""
function storm_sequence(; hrms_peak::Float64=2.0, tp_peak::Float64=10.0,
                          surge_peak::Float64=1.5, angle::Float64=0.0,
                          ramp_hours::Real=6, duration_hours::Real=24,
                          dt_hours::Float64=0.5)
    total_s = duration_hours * 3600.0
    ntimes = round(Int, duration_hours / dt_hours) + 1
    timebc = collect(range(0.0, total_s; length=ntimes))
    ramp_s = ramp_hours * 3600.0
    mid_s = total_s / 2

    hrmsbc = [hrms_peak * _storm_envelope(t, ramp_s, total_s) for t in timebc]
    swlbc  = [surge_peak * _storm_envelope(t, ramp_s, total_s) for t in timebc]
    tpbc   = [tp_peak * (0.7 + 0.3 * _storm_envelope(t, ramp_s, total_s)) for t in timebc]
    return (; timebc, hrmsbc, tpbc, swlbc, wangbc=fill(angle, ntimes))
end

function _storm_envelope(t, ramp_s, total_s)
    mid = total_s / 2
    if t < ramp_s
        return t / ramp_s
    elseif t > total_s - ramp_s
        return (total_s - t) / ramp_s
    else
        return 1.0
    end
end

"""
    seasonal_waves(; hrms_summer=0.5, hrms_winter=2.0, tp=8.0,
                     duration_years=1, dt_hours=1.0) -> NamedTuple

Sinusoidal seasonal wave variation (min in summer, max in winter).
"""
function seasonal_waves(; hrms_summer::Float64=0.5, hrms_winter::Float64=2.0,
                          tp::Float64=8.0, duration_years::Real=1,
                          dt_hours::Float64=1.0)
    total_s = duration_years * 365.25 * 86400.0
    ntimes = round(Int, duration_years * 365.25 * 24 / dt_hours) + 1
    timebc = collect(range(0.0, total_s; length=ntimes))
    year_s = 365.25 * 86400.0

    hrms_mean = (hrms_summer + hrms_winter) / 2
    hrms_amp  = (hrms_winter - hrms_summer) / 2

    hrmsbc = [hrms_mean + hrms_amp * cos(2π * t / year_s) for t in timebc]
    return (; timebc, hrmsbc, tpbc=fill(tp, ntimes),
              swlbc=fill(0.0, ntimes), wangbc=fill(0.0, ntimes))
end

# ============================================================================
# Thermal forcing builder
# ============================================================================

"""
    seasonal_temperature(; T_summer=10.0, T_winter=-20.0,
        T_water_summer=8.0, T_water_winter=2.0,
        duration_years=1, dt_hours=1.0) -> NamedTuple

Sinusoidal air and water temperature over a seasonal cycle.
Returns `(thermal_time, T_air, T_water)` for `build_config`.
"""
function seasonal_temperature(;
    T_summer::Float64=10.0, T_winter::Float64=-20.0,
    T_water_summer::Float64=8.0, T_water_winter::Float64=2.0,
    duration_years::Real=1, dt_hours::Float64=1.0,
)
    total_s = duration_years * 365.25 * 86400.0
    ntimes = round(Int, duration_years * 365.25 * 24 / dt_hours) + 1
    thermal_time = collect(range(0.0, total_s; length=ntimes))
    year_s = 365.25 * 86400.0

    T_air_mean = (T_summer + T_winter) / 2
    T_air_amp  = (T_summer - T_winter) / 2
    T_water_mean = (T_water_summer + T_water_winter) / 2
    T_water_amp  = (T_water_summer - T_water_winter) / 2

    # Summer peak at t = year/2 (day 182)
    T_air   = [T_air_mean + T_air_amp * sin(2π * t / year_s - π/2) for t in thermal_time]
    T_water = [T_water_mean + T_water_amp * sin(2π * t / year_s - π/2) for t in thermal_time]

    return (; thermal_time, T_air, T_water)
end

"""
    seasonal_snow(; max_depth=0.3, duration_years=1, dt_hours=1.0) -> NamedTuple

Prescribed seasonal snow depth time series for Arctic conditions.
Snow builds linearly from 0 at the start of winter (day 0 / day 365) to
`max_depth` at mid-winter, then melts linearly back to 0 by summer.
Returns `(snow_time, snow_depth)` aligned to the same time grid as
`seasonal_temperature`.

# Example
```julia
sf = seasonal_snow(; max_depth=0.25, duration_years=1)
```
"""
function seasonal_snow(;
    max_depth::Float64=0.3,
    duration_years::Real=1,
    dt_hours::Float64=1.0,
)
    total_s = duration_years * 365.25 * 86400.0
    ntimes = round(Int, duration_years * 365.25 * 24 / dt_hours) + 1
    snow_time = collect(range(0.0, total_s; length=ntimes))
    year_s = 365.25 * 86400.0

    # Snow present in winter (T_air < 0), absent in summer.
    # Use a cosine that peaks at t=0 (mid-winter) and is zero at t=year/2.
    # snow_depth = max_depth * max(cos(2π t / year), 0)
    snow_depth = [max_depth * max(cos(2π * t / year_s), 0.0) for t in snow_time]

    return (; snow_time, snow_depth)
end
