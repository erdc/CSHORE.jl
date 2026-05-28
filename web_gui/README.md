# CSHORE.jl Web GUI

Browser-based front-end for CSHORE.jl quick-runs. Companion to
`qml_gui/` (desktop QML app) — same physics, same `run_simulation!`
call, different front-end.

The stack is intentionally minimal:

- **Server:** `HTTP.jl` in `app.jl` (~370 lines). Two routes — `GET /`
  serves the static HTML page, `POST /run` builds a `CshoreConfig` from
  form JSON, runs the simulation, reads the resulting NetCDF, and
  returns the plot arrays as JSON.
- **Frontend:** one static HTML file (`public/index.html`) with vanilla
  JS. Plotly.js (loaded from CDN) renders the charts client-side.
- **Deps:** `HTTP`, `JSON3`, `NCDatasets`, `CSHORE`. That's it.

No Genie / Stipple / Alpine / CairoMakie. Cold image fits inside Fly's
2 GB builder.

## Run locally

```bash
julia --project=web_gui --threads=auto web_gui/app.jl
```

Then open <http://localhost:8000>. First launch downloads + compiles
deps (~5 min); subsequent launches are fast.

To skip the startup warmup (useful when iterating on `app.jl` and you
don't need the JIT pre-pay):

```bash
WEB_GUI_SKIP_WARMUP=1 julia --project=web_gui --threads=auto web_gui/app.jl
```

## Deploy

See [`DEPLOY.md`](DEPLOY.md) for the Fly.io walkthrough. Single
`Dockerfile`, single `fly deploy` command.

## What the form covers

| Tab          | Inputs                                                                                 |
| ------------ | -------------------------------------------------------------------------------------- |
| Case         | name, profile preset (planar / beach-dune), depth, slope, backshore, dune, grid dx     |
| Env          | Hrms / Tp / SWL / duration; optional sinusoidal tide; optional sin² storm recurrence   |
| Sediment     | grain sizes + mass fractions (multifraction supported; fractions auto-normalized)      |
| Thermal      | T_air / T_water / T_init (Arctic permafrost / active-layer model)                      |
| Structures   | marsh / dune grass / oyster reef / breakwater; per-type fields only                    |
| Params       | free CSHORE coefficients (`effb`, `efff`, `blp`, `slp`, `tanphi`)                      |

Output panel renders three plots via Plotly:

- Bed elevation snapshots over time
- Δz heatmap (bed change vs. time and cross-shore distance)
- Multifraction-only: surface composition + d50 heatmap

Each plot's hover toolbar includes a camera icon that exports a 2× PNG.

## Architecture notes

### Static asset

`public/index.html` is the entire UI. Vanilla JS, no framework, no
build step. The form fields use plain `name=` attributes;
`readForm()` scrapes them into a plain object and POSTs to `/run`.

If a CDN dep (Plotly) fails to load or any JS error fires, a red
banner appears at the top of the page surfacing the error — no need
to open devtools to diagnose "the button doesn't do anything".

### Server warmup

`app.jl::_warmup()` runs one tiny single-fraction and one tiny
multifraction simulation before `HTTP.serve` binds. This pre-pays JIT
specialization for both code paths so the user's first real `/run` is
just the simulation work (~10–60 s) instead of multi-minute compile.

Set `WEB_GUI_SKIP_WARMUP=1` to bypass during local iteration.

### Cost cap

Each request is rejected if `duration_h / dx_m² > 50_000` (`MAX_COST_PROXY`
in `app.jl`). Tune to your VM size. Default is calibrated for
shared-cpu-2x — roughly 10 minutes worst case.

### Multifraction normalization

The form-side validation auto-normalizes mass fractions to sum to 1.0
within ±0.05, then `app.jl` renormalizes again to floating-point
exact before calling `MultifractionConfig` (whose `validate` requires
sum == 1.0 at atol=1e-9).

## Resource limits to add before public deploy

The defaults are fine for a small trusted audience. Before posting the
URL widely:

1. **Rate limiter** — Fly's edge doesn't include one. Either token-
   bucket inside `handler()` keyed by IP or put Cloudflare in front.
2. **`runs/` janitor** — runs accumulate inside the machine until the
   next restart. A small periodic cleanup keeps disk usage bounded.
3. **One-job-at-a-time queue** — if you expect concurrent users on a
   small VM, wrap `run_simulation_json` in a `lock(::ReentrantLock)`
   so Julia's JIT doesn't thrash on contended specialization.
