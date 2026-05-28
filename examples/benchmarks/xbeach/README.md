# XBeach Comparison Benchmarks

Tests based on published cases from the XBeach validation suite and
the broader coastal engineering literature. These use well-documented
wave conditions and compare CSHORE.jl's qualitative response to what
XBeach (and Delft3D, CSHORE-FORTRAN) produce for the same inputs.

## Tests

### `lip_11d_test_1c.jl` — LIP 11D Delta Flume, Test 1C

The **standard cross-shore morphodynamic benchmark**. Delta Flume
experiments by Roelvink & Reniers (1995) covering stable (1A),
erosive (1B), and accretive (1C) beach states under irregular waves.

**Test 1C is the accretive case**: moderate waves migrate a
pre-existing bar (from the end of 1B) onshore.

**Conditions:**
- Flume: 200 m × 4.1 m working section
- Initial profile: final state of LIP-1B (barred, erosive)
- Waves: Hm0 = 0.6 m, Tp = 8.0 s, JONSWAP γ=3.3
- Duration: 13 hours
- Sediment: d50 = 0.22 mm

**Expected response:**
- Bar migrates **onshore** by ~5-10 m over the experiment
- Post-breaking wave heights decrease across bar
- Undertow peaks over the bar
- Classic "dipole" bed-change pattern: erosion on seaward face,
  deposition on shoreward face

**What this tests in CSHORE.jl:**
- Wave transform over a barred bathymetry (LWAVE + DBREAK)
- Cross-shore transport balance (undertow offshore vs. asymmetry onshore)
- The IASYM=1 (Ruessink 2012 skewness) option — onshore-biased bedload
  is the mechanism that reproduces bar migration in accretive cases
- Mass conservation (total volume should not drift)

**Data provenance:** The initial bathymetry (`h_1C.dep`) and JONSWAP
boundary file (`jonswap.inp`) in `data/DeltaflumeLIP11D_1C/` are pulled
directly from the XBeach public skillbed repository:

    https://svn.oss.deltares.nl/repos/xbeach/skillbed/input/DeltaflumeLIP11D/1C/

The `params_original.txt` file is the original XBeach parameter file
shipped with the skillbed case (for reference; CSHORE.jl does not read
XBeach param files — the Julia script translates the relevant inputs
directly).

**Expected observation (from skillbed outputs):** Bar migrates onshore
~5–10 m over 13 hours. The Julia benchmark reports a centroid-based
bar-migration metric that should fall within this range.

## Related benchmarks

See `../field/egmond_1998_storm.jl` for a field-scale storm case
using XBeach skillbed data from the Egmond aan Zee Oct 1998 campaign.

## Planned additions

- **LIP 11D Test 1B** (erosive case) — offshore bar migration under storm waves
- **XBeach barrier_1d test** — overwash with IOVER=1 + IPERM=1
- **Deltares dune erosion test** — storm-driven scarping

## References

- Roelvink, J.A. & Reniers, A.J.H.M. (1995). *LIP 11D Delta Flume
  experiments: A dataset for profile model validation*. Delft
  Hydraulics report H2130.
- Ruessink, B.G., Ramaekers, G., & van Rijn, L.C. (2012). On the
  parameterization of the free-stream non-linear wave orbital motion
  in nearshore morphodynamic models. *Coastal Engineering*, 65, 56–63.
- [XBeach documentation](https://oss.deltares.nl/web/xbeach/)
- [CROCO test cases — sandbar (LIP-1C)](https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.test_cases.sandbar.html)

## Running

```bash
cd /path/to/CSHORE.jl
julia --project=examples examples/benchmarks/xbeach/lip_11d_test_1c.jl
```

Output figures are saved to `examples/benchmarks/output/`.
