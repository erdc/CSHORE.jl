# CSHORE.jl Benchmark Suite

A collection of test cases for validating CSHORE.jl against known solutions,
laboratory experiments, and other coastal morphodynamic models.

## Structure

```
benchmarks/
├── analytical/    # Closed-form solutions (no data required)
├── field/         # Real-world field data (Duck FRF, Egmond, etc.)
├── xbeach/        # Published XBeach test cases (LIP 11D, etc.)
└── output/        # Generated plots (gitignored)
```

## How to run

Each benchmark is a standalone Julia script that uses CSHORE.jl and
CairoMakie to produce a comparison figure in `output/`. To run:

```bash
cd /path/to/CSHORE.jl
julia --project=examples examples/benchmarks/analytical/01_linear_shoaling.jl
```

Or run them all:

```bash
julia --project=examples -e '
  for f in readdir("examples/benchmarks/analytical"; join=true)
      endswith(f, ".jl") || continue
      println("=== $f ===")
      include(f)
  end
'
```

## Analytical benchmarks

Each test sets up a controlled scenario where the answer is known in
closed form. CSHORE.jl should match these to within discretization error.

| Test | Validates | Expected result |
|------|-----------|-----------------|
| 01_linear_shoaling | LWAVE, Green's Law | Hrms follows `H·(c0/c)^(1/2)` in non-breaking waves |
| 02_snell_refraction | Wave angle tracking | `sin(θ)/c = const` through shoaling |
| 03_longshore_current | Longuet-Higgins 1970 | Steady LSC peaks at the break point |
| 04_dean_equilibrium | Morphodynamic stability | h=Ax^(2/3) profile remains stationary under steady waves |
| 05_undertow_balance | Mass conservation | Cross-shore mean flow balances Stokes drift |

## Field benchmarks (to be added)

- **Duck FRF / SandyDuck 1997** — storm response, N.C., USA
- **Egmond aan Zee** — Roelvink et al. 2009, Netherlands

See `field/README_data.md` for data download instructions.

## XBeach comparison (to be added)

- **LIP 11D test 1C** — Delft large wave flume, the de facto bar migration benchmark
- **Barrier overwash** — from XBeach toolbox examples

See `xbeach/README.md` for test-case setup.

## Design principles

1. **Each script is self-contained** — no external data files for analytical tests
2. **Always produce a figure** — visual inspection catches issues numbers hide
3. **Show the exact solution overlaid** on CSHORE output when possible
4. **Print diagnostic metrics** (RMSE, bias, skill) at the end
5. **Keep runtime < 1 minute** for analytical tests; field tests may take longer
