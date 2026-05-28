# ERDC USACE Testbed Benchmarks

Cross-shore profile comparisons between **Julia CSHORE.jl** and pre-computed
**FORTRAN CSHORE_USACE** reference outputs for 8 standard testbed cases from
the ERDC USACE distribute bundle.

## Cases

| Slug | Description | ILAB | IPROFL |
|------|-------------|------|--------|
| `agate_sep` | Agate Beach, OR — September 2013 storm | 0 | 1 |
| `agate_oct` | Agate Beach, OR — October 2013 storm | 0 | 1 |
| `frf_runup_2` | FRF BathyDuck runup event #2 | 1 | 0 |
| `frf_runup_3` | FRF BathyDuck runup event #3 | 1 | 0 |
| `frf_morpho_070` | FRF general morphology run 070 | 0 | 1 |
| `frf_morpho_071` | FRF general morphology run 071 | 0 | 1 |
| `gee` | GEE laboratory flume | 1 | 1 |
| `lstf` | LSTF laboratory basin | 1 | 1 |

## Directory layout

```
erdc_usace/
├── compare_utils.jl        # shared parser, BenchmarkCase, run_and_compare
├── run_benchmarks.jl       # standalone runner — prints RMS table
├── plot_benchmarks.jl      # generates erdc_comparison.png
├── README.md
└── cases/
    └── <slug>/
        ├── infile           # original FORTRAN infile
        ├── ref_OSETUP.txt   # final time-step block from FORTRAN OSETUP output
        ├── ref_OBPROF.txt   # final time-step block from FORTRAN OBPROF output
        └── ref_OXVELO.txt   # final time-step block from FORTRAN OXVELO output
```

Reference files are compact single-block extracts (the last time step only)
from the full FORTRAN output files, formatted as:

```
L  NJ  TIME
x1  val1 ...
x2  val2 ...
...
```

## Running

**RMS comparison table:**
```bash
julia --project=../../.. examples/benchmarks/erdc_usace/run_benchmarks.jl
```

**Multi-panel comparison figure:**
```bash
julia --project=../../.. examples/benchmarks/erdc_usace/plot_benchmarks.jl
# → examples/benchmarks/erdc_usace/erdc_comparison.png
```

**As part of the test suite:**
```bash
julia --project=. test/runtests.jl
```

## RMS thresholds (test suite)

| Quantity | Threshold |
|----------|-----------|
| Hrms | < 0.05 m |
| Δzb (bed-level change) | < 0.10 m |
| umean (cross-shore) | < 0.15 m/s |

Wave setup is not threshold-tested because the reference SWL at the final
time step must be estimated from the infile time series, introducing
systematic offsets for field cases with time-varying tides.

## Reference column layout

| File | Col 1 | Col 2 | Col 3 | Col 4 |
|------|-------|-------|-------|-------|
| `ref_OSETUP.txt` | x (m) | MWS = SWL + setup (m) | h (m) | σ (Hrms/√8) |
| `ref_OBPROF.txt` | x (m) | zb (m) | | |
| `ref_OXVELO.txt` | x (m) | umean (m/s) | ustd (m/s) | |
