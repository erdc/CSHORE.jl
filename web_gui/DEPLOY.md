# Deploying the CSHORE.jl Web GUI on Fly.io

About 10 minutes of setup + ~$25/mo for a shared-cpu-4x / 8 GB machine.
You get a public HTTPS URL and a 24/7 instance.

The stack is intentionally lean: HTTP.jl serves a single static HTML
page plus a `/run` endpoint that drives `run_simulation!` and replies
with JSON. Plots render client-side in Plotly.js. No Genie / Stipple /
CairoMakie. Cold-boot fits inside Fly's default 2 GB remote builder.

## One-time setup

1. **Install flyctl**
   ```bash
   brew install flyctl                            # macOS
   curl -L https://fly.io/install.sh | sh         # Linux / WSL
   ```

2. **Create a Fly account** and add a card on file.
   ```bash
   fly auth signup    # or `fly auth login`
   ```

3. **Edit `fly.toml`** at the repo root. Change `app = 'cshore-jl'` to
   a globally-unique name (and update the matching `[build] image =
   'registry.fly.io/<name>:latest'` line), and set `primary_region`
   to whatever's closest to your users (`iad` = Virginia, `sjc` = San
   Jose, `ams` = Amsterdam, `lhr` = London).

## Deploy

Deploys are automated: pushing to `main` triggers
`.github/workflows/docker-publish.yml`, which builds the Docker image
on GitHub's runners, pushes it to Fly's own private registry
(`registry.fly.io/cshore-jl`), then runs `flyctl deploy --image …` so
Fly just pulls the prebuilt image. Wall-clock time is typically 5–15 min
on the first run (cold build cache) and 3–8 min thereafter.

For an ad-hoc deploy outside of CI:

```bash
fly auth docker                            # one-time, sets up Docker creds
fly deploy                                 # pulls :latest from registry.fly.io
fly deploy --dockerfile web_gui/Dockerfile # forces a local rebuild
                                            # (slow path, only if needed)
```

`fly logs` shows progress — look for the `web_gui: warmup complete`
line. After that the page is ready.

### Required GitHub Actions secret

The workflow expects one repository secret, set under
**Settings → Secrets and variables → Actions → New repository secret**:

| Name | What it is | How to mint |
| --- | --- | --- |
| `FLY_API_TOKEN` | Fly deploy token, scoped to the `cshore-jl` app. Used both to push the image to `registry.fly.io/cshore-jl` and to run `flyctl deploy`. | `fly tokens create deploy -a cshore-jl` |

We use Fly's own registry (not GHCR or Docker Hub) because Fly's
machine infra auto-authenticates to it — no second pull credential
needed at runtime, no PAT to rotate, and recent `flyctl` removed the
`--image-username` / `--image-password` flags that would have been
needed to deploy from a private external registry.

## Optional: GitHub Pages mirror of the UI (infrastructure ready, workflow deferred)

The plumbing to serve `web_gui/public/index.html` from GitHub Pages while
the model continues to run on Fly is already in place:

- `web_gui/public/index.html` defines `const API_BASE = window.CSHORE_API_BASE || ''`.
  All `fetch(API_BASE + '/run', ...)` calls go to relative URLs by
  default (same-origin Fly deploy), or to a remote backend when a
  static-host deployment injects `window.CSHORE_API_BASE` before the
  main script tag.
- `app.jl` returns `Access-Control-Allow-Origin: *` on every response
  (including OPTIONS preflights for cross-origin POSTs from a static
  Pages-hosted UI).

The Pages publish workflow (`.github/workflows/pages.yml`) was added and
then removed, pending repo-owner action: GitHub Pages must first be
enabled under **Settings → Pages → Source = "GitHub Actions"**. To
re-enable Pages hosting later, restore that workflow file from git
history (it copied `web_gui/public/` to a staging dir, injected the
backend URL via an awk one-liner, and used `actions/deploy-pages@v4`).

## Day-to-day operations

```bash
fly logs                       # tail server logs
fly status                     # machine state, IPs, region
fly machine restart            # bounce after a config-only change
fly scale memory 4096          # bump RAM if runs OOM
fly machine stop / start       # pause / resume billing
fly releases                   # list previous deploys
fly releases revert <version>  # roll back
```

## Cost expectations

| Tier                 | Hourly    | Monthly (24/7) |
| --------------------- | --------- | -------------- |
| `shared-cpu-2x` 4 GB *(minimum — runtime precompile thrashes here, expect 30-60+ min boots)* | ~$0.007 | ~$5 |
| `shared-cpu-2x` 8 GB *(workable; 2 parallel precompile workers fit)* | ~$0.013 | ~$10 |
| `shared-cpu-4x` 8 GB *(matches `fly.toml` — 4 workers fit, ~5-10 min cold boot. 8 GB is the cap for this VM size.)* | ~$0.034 | ~$25 |
| `performance-2x` 8 GB *(dedicated CPU; ~3-5 min cold boot)* | ~$0.05 | ~$35 |
| `shared-cpu-8x` 16 GB *(if 8 GB ever feels tight)* | ~$0.07 | ~$50 |

Fly bills per second. Stop the machine when idle if you want to pay
only for time it's running.

## Runtime caveats

1. **`runs/` is ephemeral.** Every restart wipes outputs. Each `/run`
   writes its NetCDF + plots into `/app/runs/<case>_<stamp>_lean/` for
   the lifetime of the machine.

2. **Auto-stop is disabled on purpose.** Cold-boot is ~5–10 min
   (precompile + warmup). Scaling to zero would make every fresh
   visitor wait that long.

3. **One process, multi-thread.** A second concurrent `/run` runs in
   parallel via `--threads=auto`, but Julia's JIT is global so the
   second request blocks compilation. After warmup this is rarely
   noticeable; you can lower `http_service.concurrency.hard_limit` if
   you ever see contention.

4. **No persistent state.** `runs/` and `.julia/` are inside the
   container. To survive deploys, mount a Fly Volume — see Fly's
   volumes docs.

## Hardening before exposing the URL publicly

The default config is fine for a small trusted audience. Before
posting the URL widely:

- **Cost cap is already enforced.** `app.jl` rejects any run with
  `duration_h / dx_m² > 50_000` (~10 min worst case on shared-cpu-2x).
  Tune `MAX_COST_PROXY` in `app.jl` if you want a tighter bound.
- **Add a rate limiter.** Fly Proxy doesn't include one out of the
  box; either use a token bucket inside `handler()` or put a Cloudflare
  Worker in front.
- **Janitor** — `runs/` grows forever otherwise. A cron job over `fly
  ssh console` cleaning entries older than N days is enough.

## If the deploy fails

| Symptom                                         | Likely cause                | Fix                                                                  |
| ----------------------------------------------- | --------------------------- | -------------------------------------------------------------------- |
| `fatal: lookup app: Could not find App`         | App name not registered     | `fly apps create <name>` then `fly deploy`                           |
| Build hangs on `Pkg.instantiate` for > 25 min   | Network slow                | Wait. First-time dep download is genuinely slow.                     |
| App boots but never serves                      | Stuck in warmup             | `fly logs` — look for the `warmup complete` line.                    |
| HTTP 500 on every request                       | App crashed at boot         | `fly logs` — look for ERROR / LoadError near the top.                |
| OOM at runtime                                  | 4 GB too small for your run | `fly scale memory 8192` and bump the VM size.                        |
| Page renders but Run button is dead             | JS / CDN failure            | Check the red banner at the top of the page — it surfaces JS errors. |
