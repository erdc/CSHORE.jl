# Example CSVs for the CSHORE.jl QML GUI

Six sample input files in the formats the GUI's "Browse…" buttons expect.
All six are reproducible by running [`generate_examples.jl`](generate_examples.jl)
in the same folder.

## Bathymetry CSVs (columns: `x`, `z` in metres)

| File | Profile | Notes |
|---|---|---|
| `bathy_planar.csv`     | Straight slope, 1:20 (0.05) over 300 m, offshore depth 8 m | Sanity check / parity with the GUI's `planar_beach` preset. |
| `bathy_beach_dune.csv` | 1:25 slope + Gaussian dune crest at x ≈ 240 m, peak ≈ +6.6 m | Useful for testing morphology around a dune toe. |
| `bathy_barred.csv`     | 1:40 slope + submerged sand bar at x ≈ 120 m, ~1.2 m above the trend | Triggers the surf-zone breaking and undertow that bars are known for. |

## Wave / SWL CSVs (columns: `time` (s), `hrms` (m), `tp` (s), `swl` (m))

| File | Forcing | Notes |
|---|---|---|
| `waves_constant.csv` | Hrms=1.0, Tp=8.0, SWL=0.5 for 12 h hourly | Parity check vs running with the GUI form. |
| `waves_storm.csv`    | 36 h triangular envelope: Hrms 1 → 3 → 1 m, Tp scaling 6 → 10 s, surge bump up to 0.8 m at peak | Stress-test for morphological response. |
| `waves_tide.csv`     | 24 h gentle waves (~0.6 m) with semi-diurnal M2 tide ±0.75 m at 30-min cadence | Checks that SWL variation drives swash / runup. |
| `waves_year.csv`     | **Full year, hourly (8760 rows)**. Combines a seasonal Hrms baseline (winter ~1 m, summer ~0.4 m), 7 idealised storms (peak Hrms 1.8–3.4 m, lasting 1.5–3 days each, placed across fall/winter/spring), Tp correlated with Hrms (6–12 s), an M2 tide (~12.42 h) modulated by a spring/neap envelope (~14.77 d), and storm surge (0.3–1.2 m peaks) stacked on the tide. Realistic SWL range: ±0.8 to +1.9 m. | Year-long thermal/morphological demos. Pairs with `temps_arctic_year.csv` for a full annual coastal-permafrost simulation. |

## Thermal forcing CSVs (columns: `time` (s), `T_air` (°C), `T_water` (°C); optional `snow_depth` (m))

Loaded via the **Thermal CSV** browse button when the *Enable thermal model* box is checked. When no CSV is loaded but thermal is enabled, the GUI uses the constant *Air temperature* / *Water temperature* fields broadcast across the wave time grid.

| File | Forcing | Notes |
|---|---|---|
| `temps_constant.csv`     | T_air = -5 °C, T_water = 0 °C, hourly for 24 h | Parity check vs constant thermal mode (no CSV). |
| `temps_arctic_year.csv`  | 365-day record: T_air sin(±18.5 °C ± offset), T_water lagged & damped, plus a snow-depth column | Drives a full seasonal active-layer / permafrost cycle. ALT plot shows the seasonal thaw envelope. |

## How to use

In the GUI:

1. Click **Browse…** next to *Bathy CSV* and pick one of the `bathy_*.csv` files.
2. Click **Browse…** next to *Waves CSV* and pick one of the `waves_*.csv` files.
3. The Geometry / Wave-forcing GroupBoxes gray out (they're now overridden by the CSVs).
4. Click **Run simulation**. The progress bar fills as BC windows complete; the plot panel populates at the end.

You can mix and match — e.g., the storm waves on the barred profile is a sensible morphology test.

## Format details

The CSV reader is permissive about column ordering and case, but the names
must match (case-insensitive):

- **Bathymetry**: requires `x` and `z` (or aliases like `x_m`, `elevation`, `depth`).
- **Waves**: requires `time`, `hrms`, `tp`, `swl`. Optional `wangle` (radians).

Times must start at 0 and increase monotonically. Hrms / Tp / SWL can vary
freely between rows; the simulator linearly interpolates between them.
