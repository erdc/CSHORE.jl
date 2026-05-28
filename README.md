# CSHORE.jl

Cross-shore wave, hydrodynamic, sediment-transport, and morphology model in Julia.

Julia port of the original FORTRAN CSHORE
([erdc/cshore](https://github.com/erdc/cshore)), with extensions for Arctic
thermal / active-layer processes and multi-fraction sediment transport.

The code is under early development and therefore all features may not work or be fully tested.

---

## Installation

Install via [juliaup](https://github.com/JuliaLang/juliaup):

```bash
curl -fsSL https://install.julialang.org | sh
juliaup add 1.10
juliaup default 1.10
```

Clone the repo and instantiate the project:

```bash
git clone https://github.com/ncohn/CSHORE.jl.git
cd CSHORE.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

---

## Running the GUI

The QML.jl GUI under [`qml_gui/`](qml_gui/) is the user-facing front-end. It drives
`run_simulation!` directly from a form, supports CSV imports for bathymetry / waves
/ thermal forcing, multi-fraction sediment, optional thermal-permafrost coupling,
and produces a multi-panel result figure.

### Launch

| Platform | Command |
| -------- | ------- |
| macOS    | Double-click [`qml_gui/run_qml.command`](qml_gui/run_qml.command) (right-click → Open the first time, Gatekeeper) |
| Windows  | Double-click [`qml_gui/run_qml.bat`](qml_gui/run_qml.bat)               |
| Linux    | `./qml_gui/run_qml.sh` from a terminal                                  |

---
## Running CSHORE Benchmark Cases

```bash
julia --project=. examples/benchmarks/erdc_usace/run_benchmarks.jl
julia --project=. examples/benchmarks/erdc_usace/plot_benchmarks.jl
```

## Analytical Benchmarks

Seven simple cases with closed-form expected solutions (linear shoaling, Snell's
refraction, Dean equilibrium profile, undertow balance, longshore current, etc.):

```bash
for f in examples/benchmarks/analytical/*.jl; do
    julia --project=. "$f"
done
```

---

## Setting up a new case

Three input paths are supported. Pick whichever fits your workflow.

### 1. FORTRAN-compatible `.infile` for legacy USACE workflows

```julia
using CSHORE
cfg = read_infile("path/to/run.infile")
run_simulation!(cfg; outfile="run.nc", outdir="./out")
```

### 2. CSHORE.jl native `.cshore` TOML

A TOML file plus separate CSVs for bathymetry and wave time series. Easier to edit
than the FORTRAN format:

```toml
# example.cshore
[grid]
dx = 1.0

[bathymetry]
file = "bathy.csv"          # columns: x, z

[boundary]
bc_file = "waves.csv"       # columns: time, hrms, tp, swl  (wangle optional)

[sediment]
d50 = 3e-4
```

```julia
cfg = read_cshorejl("example.cshore")
run_simulation!(cfg; outfile="example.nc")
```

### 3. Programmatic config

```julia
cfg = build_config(;
    dx = 1.0,
    bathymetry_x = collect(0.0:1.0:300.0),
    bathymetry_z = [-8.0 + 0.05*x for x in 0.0:1.0:300.0],
    timebc       = [0.0, 43200.0],
    hrmsbc       = [1.0, 1.0],
    tpbc         = [8.0, 8.0],
    swlbc        = [0.5, 0.5],
    sediment     = make_sediment(d50 = 3e-4, effb = 0.005),
    multifraction = MultifractionConfig(
        grain_sizes      = [0.15e-3, 0.30e-3, 0.60e-3],
        initial_fractions = [0.3, 0.5, 0.2],
    ),
    options      = OptionFlags(iprofl = 1),
)
run_simulation!(cfg; outfile = "out.nc")
```
---

## Output format

`run_simulation!(...; outfile="run.nc")` writes a CF-1.10 compliant NetCDF:

| Dimension      | Meaning                          |
| -------------- | -------------------------------- |
| `x`            | Cross-shore distance (m)         |
| `time`         | Time (unlimited)                 |
| `fraction`     | Grain-size class (m)             |
| `layer`        | Bed layer index (1 = active top) |

**2D fields** `(x, time)`: `zb`, `zb_hard`, `hrms`, `wsetup`, `sigma`, `umean`,
`ustd`, `qbreak`, `q_total`, `pwet`, `hwd`. With thermal: `ALT`, `T_surface`.
**3D per-fraction**: `qbx`, `qsx`, `bed_mass`. 

---

## Running the test suite

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
---

## Repository layout

```
src/
├── CSHORE.jl
├── cli.jl
├── config.jl
├── state.jl
├── driver.jl
├── bmi.jl
├── presets.jl
├── data.jl
├── plotting.jl
│
├── hydrodynamics/
│   ├── waves.jl
│   ├── hydro.jl
│   ├── wetdry.jl
│   ├── infragravity.jl
│   └── overtopping.jl
│
├── sediment/
│   ├── transport.jl
│   ├── exner.jl
│   ├── composition.jl
│   ├── fractions.jl
│   ├── erosion.jl
│   ├── diffusion.jl
│   └── avalanche.jl
│   └── cohesive.jl
│   └── clay_dike_erosion.jl
|
├── aeolian/
│   ├── aeolian.jl
│   └── wind_shear.jl
│
├── groundwater/
│   ├── groundwater.jl
│   └── porous.jl
│
├── vegetation/
│   └── vegetation.jl
│
├── utils/
│   ├── utilities.jl
│   └── provenance.jl
│
├── thermal/
│   └── thermal.jl
│
├── q2d/
│   └── transect_grid.jl
│
└── io/
    ├── input.jl
    ├── cshorejl.jl
    ├── xbeach.jl
    ├── netcdf.jl
    └── output.jl

ext/
├── CSHOREMakieExt.jl
└── CSHOREDataExt.jl

qml_gui/
examples/
test/
```
