# CSHORE.jl QML GUI (prototype)

A native-desktop GUI for CSHORE.jl quick-runs, built with [QML.jl](https://github.com/JuliaGraphics/QML.jl)
on top of the Qt 6 framework. Calls `run_simulation!` directly from the
running Julia process — same code path as the Pluto notebook in
`../gui/cshore_pluto.jl`, but in a real desktop window with native widgets,
file dialogs, and menus.

## Status: prototype

This is a working first cut. What's in:

- Form fields for the same parameters as the Python / Pluto GUIs
  (case name, profile preset, geometry, wave forcing, sediment, output dir)
- Two-way binding between QML form fields and a Julia `JuliaPropertyMap`
- A "Run simulation" button that calls `run_quick_sim` in Julia
- A status pane that updates from a Julia `Observable`
- A `BusyIndicator` while the simulation runs
- Per-platform double-click launchers (`run_qml.bat / .command / .sh`)

What's **not** in the prototype yet:

- Inline plot of the result (would require Makie + image embedding via
  QML's `Image` element — separate follow-up)
- Async simulation execution — currently the UI blocks while `run_simulation!`
  runs. For 12 h × dx=1 m runs this is sub-second so it's fine, but multi-day
  fine-grid runs will freeze the window. Fix is `Threads.@spawn` + Observable
  notifications; ~10 lines, see `cshore_qml.jl :: run_quick_sim`.
- Custom-CSV bathymetry import (currently only the two presets)
- PackageCompiler bundle target
- Code-signing for distribution

## Requirements

- **Julia 1.10 LTS** on PATH (juliaup recommended).
- A working desktop environment. On Linux you'll also need the system
  XCB libraries Qt links against — see comments in `run_qml.sh`.

The first launch will resolve and install QML.jl + Qt6 binaries
(~300–400 MB download), so initial startup takes a few minutes. Subsequent
launches start in seconds.

## How to launch

| Platform | Command                          |
| -------- | -------------------------------- |
| Windows  | Double-click `run_qml.bat`       |
| macOS    | Double-click `run_qml.command` (right-click → Open the first time, Gatekeeper) |
| Linux    | `./run_qml.sh` from a terminal   |

You can also invoke Julia directly from the repo root:

```bash
julia --project=qml_gui qml_gui/cshore_qml.jl
```

## How it's wired

The Julia file (`cshore_qml.jl`):

1. Activates and instantiates the local `qml_gui` environment, dev-installing
   the parent CSHORE.jl package on first launch.
2. Defines a `run_quick_sim(params::JuliaPropertyMap)` function that reads
   form values out of `params`, calls `build_config` + `run_simulation!`,
   and writes status updates to an `Observable`.
3. Registers that function via `@qmlfunction` so QML can call it.
4. Constructs a `JuliaPropertyMap` of default parameter values and an
   Observable for the status text.
5. Calls `loadqml(...)` to bind both into the QML scene, then `exec()` to
   start the Qt event loop.

The QML file (`cshore_main.qml`):

- Declares an `ApplicationWindow` containing four `GroupBox` sections
  (Case / Geometry / Wave forcing / Sediment & output) plus the Run row
  and status pane.
- Each `TextField` reads from `params.<key>` (the JuliaPropertyMap) and
  writes back via `onEditingFinished`. The two-way binding is automatic.
- The Run button calls `Julia.run_quick_sim(params)`. The `params` object
  is the same object the form fields write to, so Julia sees the latest
  values.
- The `statusMsg` Observable is bound directly to the read-only
  `TextArea.text` — Julia setting it to a new string updates the UI live.

## Trade-offs vs the other GUIs in this repo

|                              | Python tkinter (`gui/`) | Pluto (`gui/`)        | QML (`qml_gui/`)       |
| ---------------------------- | ----------------------- | --------------------- | ---------------------- |
| User needs Julia installed   | No (uses compiled exe)  | Yes                   | Yes                    |
| User needs Python installed  | Yes                     | No                    | No                     |
| Calls run_simulation! directly | No (launches subprocess) | Yes                | Yes                    |
| Native-feeling UI            | No (Tk-style)           | No (browser tab)      | Yes (Qt widgets)       |
| First launch download size   | ~0 MB                   | ~150 MB (Pluto deps)  | ~300 MB (Qt6_jll)      |
| Live reactive parameters     | No (manual Run)         | Yes (slider → re-run) | No (manual Run)        |
| Inline plot of result        | No                      | Yes                   | Not yet (todo)         |

## Known issues

1. **macOS Gatekeeper warning on first launch.** `run_qml.command` is an
   unsigned shell script; you'll get an "unidentified developer" prompt
   on first run. Right-click → Open to bypass.
2. **Long simulations freeze the UI.** Synchronous run path; use the
   async refactor mentioned above for multi-day runs.
3. **Qt6_jll download is large.** First launch can take 5–10 min on slow
   connections.
4. **Display scaling on some Linux DEs** can render too small or too
   large; Qt env var `QT_SCALE_FACTOR=1.5` (or whichever) fixes it.

## Next steps (roughly in priority order)

1. **Async simulation:** wrap the body of `run_quick_sim` in
   `Threads.@spawn`, push status updates from the worker via the existing
   Observable, re-enable the Run button via `RUNNING[]=false` in the
   `finally`.
2. **Inline plot:** generate a PNG with CairoMakie after each run, expose
   the path through a separate Observable, bind a QML `Image` element
   to it. This is the largest feature delta vs the Python / Pluto GUIs.
3. **Custom-CSV bathymetry:** add a "Custom CSV…" entry to the profile
   ComboBox and a `FileDialog` invocation.
4. **PackageCompiler bundle:** mirror `compile/build.jl` with a
   `compile/build_qml.jl` that calls `create_app` with this entry point.
   Result is a single-folder bundle the user double-clicks. Big lift on
   macOS due to notarization; modest on Windows / Linux.
