#!/usr/bin/env julia
# =============================================================================
# generate_examples.jl -- writes the sample CSVs in this directory.
#
# Six files: three bathymetry profiles and three wave / SWL forcing scenarios
# in the formats the QML GUI expects:
#
#   bathy_*.csv    columns: x, z          (m)
#   waves_*.csv    columns: time, hrms, tp, swl  (s, m, s, m)
#
# Run from the repo root:
#     julia qml_gui/example_csvs/generate_examples.jl
# (no Project / package deps -- only the Julia stdlib is used.)
# =============================================================================

using Printf

const HERE = @__DIR__

function write_csv(path::AbstractString, header::AbstractString, rows; fmt::AbstractString)
    open(path, "w") do io
        println(io, header)
        for row in rows
            println(io, @eval @sprintf($fmt, $(row...)))
        end
    end
    println("wrote ", path, " (", length(rows), " rows)")
end

# ---------------------------------------------------------------------------
# Bathymetry CSVs (columns: x, z)
# ---------------------------------------------------------------------------

# 1. Planar beach: bed rises shoreward at slope 0.05 from offshore depth 8 m.
let
    path = joinpath(HERE, "bathy_planar.csv")
    dx = 1.0
    L = 300.0
    slope = 0.05
    depth = 8.0
    n = floor(Int, L / dx) + 1
    rows = [(i * dx - dx, -depth + slope * (i * dx - dx)) for i in 1:n]
    write_csv(path, "x,z", rows; fmt="%.3f,%.4f")
end

# 2. Beach + dune: planar slope offshore, Gaussian dune crest landward.
let
    path = joinpath(HERE, "bathy_beach_dune.csv")
    dx = 1.0
    L = 350.0
    slope = 0.04
    depth = 8.0
    dune_h = 5.0
    n = floor(Int, L / dx) + 1
    beach_end = depth / slope        # x at z = 0 (waterline at SWL=0)
    center = beach_end + 40.0
    sigma = 18.0
    rows = Tuple{Float64,Float64}[]
    for i in 1:n
        x = (i - 1) * dx
        z = -depth + slope * x
        z += dune_h * exp(-((x - center)^2) / (2 * sigma^2))
        push!(rows, (x, z))
    end
    write_csv(path, "x,z", rows; fmt="%.3f,%.4f")
end

# 3. Barred beach: planar slope + a single submerged sand bar.
let
    path = joinpath(HERE, "bathy_barred.csv")
    dx = 1.0
    L = 400.0
    slope = 0.025
    depth = 6.0
    bar_h = 1.2          # bar height above the underlying planar slope (m)
    bar_x = 120.0        # cross-shore position of bar crest (m)
    bar_sigma = 25.0
    n = floor(Int, L / dx) + 1
    rows = Tuple{Float64,Float64}[]
    for i in 1:n
        x = (i - 1) * dx
        z = -depth + slope * x
        z += bar_h * exp(-((x - bar_x)^2) / (2 * bar_sigma^2))
        push!(rows, (x, z))
    end
    write_csv(path, "x,z", rows; fmt="%.3f,%.4f")
end

# ---------------------------------------------------------------------------
# Wave / SWL CSVs (columns: time, hrms, tp, swl)
# Time in seconds; hrms in metres; tp in seconds; swl in metres.
# ---------------------------------------------------------------------------

# 1. Constant forcing: parity case for the form-driven 12 h run.
let
    path = joinpath(HERE, "waves_constant.csv")
    dt = 3600.0   # 1 hour
    nhr = 12
    rows = Tuple{Float64,Float64,Float64,Float64}[]
    for i in 0:nhr
        push!(rows, (i * dt, 1.0, 8.0, 0.5))
    end
    write_csv(path, "time,hrms,tp,swl", rows; fmt="%.0f,%.3f,%.3f,%.3f")
end

# 2. Storm: hrms ramps up to 3 m, peaks for ~6 h, then ramps down. Tp scales
#    with hrms (a rough Hrms-Tp empirical link). SWL elevated during peak
#    (storm surge).
let
    path = joinpath(HERE, "waves_storm.csv")
    dt = 3600.0
    nhr = 36              # 36-hour storm
    rows = Tuple{Float64,Float64,Float64,Float64}[]
    for i in 0:nhr
        t = i * dt
        # Triangular envelope: 1 m baseline, 3 m peak at hour 18.
        rise = clamp(i / 14, 0.0, 1.0)
        fall = clamp((nhr - i) / 14, 0.0, 1.0)
        hrms = 1.0 + 2.0 * min(rise, fall)
        tp   = 6.0 + 4.0 * (hrms - 1.0) / 2.0     # 6 → 10 s
        # Storm surge: smooth bump centred on peak Hrms hour.
        surge = 0.8 * exp(-((i - 18)^2) / (2 * 6 * 6))
        swl = 0.0 + surge
        push!(rows, (t, hrms, tp, swl))
    end
    write_csv(path, "time,hrms,tp,swl", rows; fmt="%.0f,%.3f,%.3f,%.3f")
end

# 3. Tide-dominated: gentle waves, semi-diurnal SWL with 1.5 m range.
let
    path = joinpath(HERE, "waves_tide.csv")
    dt = 1800.0   # 30 min — finer to resolve the tide
    n = 48        # 24 hours
    rows = Tuple{Float64,Float64,Float64,Float64}[]
    for i in 0:n
        t = i * dt
        hrms = 0.6 + 0.1 * sin(2π * (t / 3600) / 24)
        tp = 7.0
        # Two M2 tide cycles (~12.42 h period) over 24 h, ±0.75 m
        swl = 0.75 * sin(2π * t / (12.42 * 3600))
        push!(rows, (t, hrms, tp, swl))
    end
    write_csv(path, "time,hrms,tp,swl", rows; fmt="%.0f,%.3f,%.3f,%.3f")
end

# 4. Full year, hourly. Combines four physically-meaningful signals:
#    (a) Seasonal Hrms baseline (winter higher, summer lower).
#    (b) An idealised storm calendar — ~7 named storms placed across
#        fall/winter/spring with Gaussian envelopes in time. Peak Hrms
#        2.0–3.5 m, lasting 1.5–3 days each.
#    (c) Tp correlated with Hrms via a steepness-saturation rule plus
#        the storm peaks.
#    (d) M2 tide (12.42 h) modulated by a spring/neap envelope (14.77 d).
#    (e) Storm surge (~ 0.3–1.2 m) stacked on the tide, peaking with
#        each storm.
#
# Roughly 8760 rows; ~260 kB file. Suitable for testing year-long thermal
# runs with realistic wave + SWL forcing. Pairs nicely with
# `temps_arctic_year.csv` for a full annual coastal-permafrost demo.
let
    path = joinpath(HERE, "waves_year.csv")
    dt = 3600.0
    n_hr = 365 * 24

    # Storm calendar: (day_of_year, peak_Hrms_m, peak_Tp_s, peak_surge_m, full_width_days)
    storms = [
        ( 30, 3.0, 11.5, 1.0, 2.5),   # late January storm
        ( 55, 3.4, 12.0, 1.2, 3.0),   # mid-February nor'easter
        ( 85, 2.6, 10.5, 0.7, 2.0),   # late March
        (115, 1.8,  9.5, 0.4, 1.5),   # late April
        (260, 2.4, 11.0, 0.9, 2.0),   # mid-September hurricane
        (300, 2.9, 11.5, 1.0, 2.5),   # late October
        (340, 3.2, 11.8, 1.1, 2.5),   # early December storm
    ]

    rows = Tuple{Float64,Float64,Float64,Float64}[]
    for i in 0:n_hr
        t = i * dt
        doy = i / 24.0                          # day-of-year (fractional)

        # (a) Seasonal Hrms baseline — peaks near day 30 (winter), trough
        # near day 213 (mid-summer). Range ~ 0.4–1.0 m background.
        baseline_hrms = 0.7 + 0.3 * cos(2π * (doy - 30) / 365)

        # (b) + (c) + (e): Storm contributions
        storm_hrms = 0.0
        storm_tp_peak = 0.0
        storm_surge = 0.0
        for s in storms
            sd, sh, st, ssu, sdur = s
            sigma = sdur / 2.0                  # half-width = sdur / 2 days
            d = doy - sd
            env = exp(-(d * d) / (2.0 * sigma * sigma))
            # Add (peak − baseline) × env so the storm only ADDS height.
            storm_hrms     += max(0.0, sh - baseline_hrms) * env
            storm_tp_peak   = max(storm_tp_peak, st * env)
            storm_surge    += ssu * env
        end
        hrms = baseline_hrms + storm_hrms

        # Tp: a baseline steepness relationship (Tp ~ 6 + 1.5 √Hrms),
        # raised to the storm peak when a storm is active.
        tp_base = 6.0 + 1.5 * sqrt(max(hrms, 0.0))
        tp = max(tp_base, storm_tp_peak)

        # (d) Tide — M2 (12.42 h) modulated by spring/neap (14.77 d).
        tide_amp  = 0.6 + 0.2 * cos(2π * t / (14.77 * 86400.0))
        tide      = tide_amp * sin(2π * t / (12.42 * 3600.0))

        swl = tide + storm_surge

        push!(rows, (t, hrms, tp, swl))
    end
    write_csv(path, "time,hrms,tp,swl", rows; fmt="%.0f,%.3f,%.3f,%.4f")
end

# ---------------------------------------------------------------------------
# Thermal forcing CSVs (columns: time, T_air, T_water; optional snow_depth)
# Time in seconds; temperatures in °C; snow depth in metres.
# ---------------------------------------------------------------------------

# 1. Constant Arctic baseline: T_air = -5 °C, T_water = 0 °C for 1 day.
let
    path = joinpath(HERE, "temps_constant.csv")
    dt = 3600.0
    nhr = 24
    rows = Tuple{Float64,Float64,Float64}[]
    for i in 0:nhr
        push!(rows, (i * dt, -5.0, 0.0))
    end
    write_csv(path, "time,T_air,T_water", rows; fmt="%.0f,%.2f,%.2f")
end

# 2. Seasonal Arctic year: 365-day record with sinusoidal air temp swing
#    -25 °C (Feb) to +12 °C (Aug); water temp lags ~30 days and damps to
#    -1.5 → +6 °C. Daily samples; useful for ALT seasonal-thaw demos.
let
    path = joinpath(HERE, "temps_arctic_year.csv")
    nday = 365
    rows = Tuple{Float64,Float64,Float64,Float64}[]
    for i in 0:nday
        t = i * 86400.0
        # day-of-year offset so the minimum is around late Jan (day 30)
        phi = 2π * (i - 30) / 365
        T_air   = -6.5 + 18.5 * sin(phi)            # range  [-25, +12]
        # Water lags 30 days, smaller amplitude
        phi_w   = 2π * (i - 60) / 365
        T_water = 2.25 + 3.75 * sin(phi_w)          # range  [-1.5, +6]
        # Snow accumulates Oct-Apr; quick melt May. Crude triangular envelope
        # peaking late winter at ~30 cm.
        snow = if 270 <= i <= 365 || i <= 120
            d_from_peak = abs(i - (i <= 120 ? 60 : 330))
            max(0.0, 0.30 * (1.0 - d_from_peak / 60))
        else
            0.0
        end
        push!(rows, (t, T_air, T_water, snow))
    end
    write_csv(path, "time,T_air,T_water,snow_depth", rows;
              fmt="%.0f,%.3f,%.3f,%.4f")
end

println("\nAll example CSVs written to: $HERE")
