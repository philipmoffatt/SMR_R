# SMR CPU Reference Release and JAX Rewrite Plan

## Status and purpose

The CPU reference path is the behavior-preservation target for a later Python
and JAX implementation.  The target is the 30-day Mica Creek run for
2015-01-01 through 2015-01-30, with the original GRASS/Perl engine retained as
the scientific and numerical reference.  The JAX work is a translation first;
changes to model science, units, or routing semantics are separate decisions.

The CPU-optimized reference completed all 30 days.  Its final `saturation` and
`storage_amt` rasters were identical to the initial reference at all 68,121
common valid cells.  The 30 daily records were present in both outputs.  The
text CSVs are not byte-identical because the optimized reporting path retains
more precision than the legacy GRASS report, not because the final states
differ.

## Frozen Mica Creek reference capture

`reference/mica-creek-jan2015/manifest.json` is the canonical inventory for
the frozen January 2015 JAX target. It hashes the compact replay bundle:

- `model_inputs.zip`: 86 model-ready GRASS/ASCII inputs;
- `gridmet_weather_2015-01.csv`: the forcing series;
- `daily_mass_balance_2015-01.csv`: the 30 legacy daily records;
- `final_state_rasters.zip`: final `saturation` and `storage_amt` rasters;
- reference configuration and the successful legacy simulation manifest.

The capture excludes raw SSURGO/source downloads and the GRASS workspace: both
are larger, reconstructible from the captured inputs, and unnecessary to run
the JAX parity case. The final-state rasters matched the CPU-optimized result
exactly at 68,121 valid cells.

The completed legacy run did not retain daily internal radiation, snow, soil,
or lateral-flow rasters. The next verification-preparation task must create a
short instrumented replay to capture those named daily/substep fixtures; it
must not substitute newly generated values for the frozen final target.

## Target package

The JAX engine will be introduced beside, rather than in place of, the
existing R/Perl/GRASS workflow.

```text
src/smr_jax/
  state.py              # immutable physical state groups
  parameters.py         # static rasters, land-cover properties, constants
  forcing.py            # daily and four six-hour forcing fields
  initialize.py         # reference-compatible initial conditions
  engine.py             # daily and subdaily composition via lax.scan
  diagnostics.py        # mass balance, annual counters, legacy reports
  fluxes/
    precipitation.py    # lapse correction and rain/snow partition
    canopy.py           # interception, evaporation, throughfall
    radiation.py        # shortwave provider boundary and longwave flux
    turbulent_heat.py   # vapor, latent, and sensible heat fluxes
    snow.py              # refreeze, melt, ice/liquid and thermal state
    soil_water.py       # input and A/B layer allocation
    lateral_flow.py     # conductivity, outflow and neighborhood inflow
    evapotranspiration.py
    percolation.py
    surface_runoff.py   # road and saturation-excess runoff
    baseflow.py         # watershed reservoir recurrence
  io/
    legacy_bundle.py    # weather, ASCII rasters, watershed properties
    raster.py
    checkpoints.py
    reports.py
```

Each flux module will be a pure function with an explicit contract:

```python
next_state, named_fluxes = flux(state, static_parameters, forcing)
```

Units and the sign convention belong in each module's public docstring and
test.  Reading ASCII/GRASS data, date handling, file writing, and report
formatting stay outside compiled functions.  This is deliberately a compact
module layout, not a plugin framework.

## State and execution contract

Static arrays are dense `[row, col]` fields with an explicit basin mask,
georeferencing metadata, and eight routing weights.  Dynamic state is grouped
as follows:

- Soil: `storage_amt`, `storage_amt_A`, and `storage_amt_B`.
- Canopy: `canopy_storage_amt`.
- Snow: ice, liquid water, SWE and previous SWE, snow age, total and surface
  energy, and surface temperature.
- Basin: reservoir volume and baseflow for each watershed.
- Diagnostics: annual runoff, saturation counts, and daily accumulators.

Initially, redundant legacy variables remain explicit.  In particular, total
soil storage must not be reconstructed from layers until parity tests establish
that this does not alter intermediate behavior.

`engine.py` must preserve the source order in `R/smr_buildup.pl`:

1. Precipitation and temperature lapse correction, then rain/snow partition.
2. Canopy interception and evaporation, which consumes PET.
3. Snow age/albedo, radiation, turbulent heat, refreeze, and melt.
4. Road runoff and the corresponding adjustment to throughfall and melt.
5. Four six-hour soil steps: input, layer allocation, lateral routing,
   reallocation, ET, percolation, reallocation, and saturation runoff.
6. Daily runoff aggregation, annual diagnostics, basin reservoir recurrence,
   and reports.

The four substeps should be a fixed-length `jax.lax.scan` nested inside a daily
`jax.lax.scan`.  CPU eager, CPU compiled, and GPU runs must call the same flux
functions and use float64 for the first parity release.

## Translation invariants and hazards

The port must reproduce the engine before attempting to improve it.

- **Radiation:** GRASS `r.sun` uses daily evolving albedo.  Recorded
  beam/diffuse/reflected fields are valid only as a replay provider for parity
  testing.  A reusable terrain/date decomposition or native JAX radiation
  implementation must be demonstrated against GRASS before it becomes the
  production provider.
- **Lateral routing:** this is an eight-neighbor, slope-weighted
  multiple-flow-direction style redistribution, not canonical D-infinity.
  The source expression gathers shifted neighbor lateral-flow fields and
  multiplies them by directional weights at the receiving expression.  A 3x3
  impulse fixture must establish its literal row/column and direction behavior
  before a conventional donor-weighted formulation is considered.
- **Sequential updates:** snow energy, refreeze, ice, liquid water, soil
  layers, ET, and percolation are sequentially overwritten in Perl.  JAX must
  snapshot old values and advance them in the same order.
- **GRASS conditional behavior:** `jnp.where` evaluates both operands.  ET
  cells with `root_zone = 0`, GRASS NULLs, and denominator thresholds require
  explicit safe-denominator and mask rules, proven with fixtures.
- **Legacy behavior is evidence, not endorsement:** negative moisture,
  dewpoint lapse algebra, slope units, differing runoff/percolation conversions,
  and February's 29-day convention must be preserved and reported initially,
  not silently corrected.
- **Reporting is currently coupled to basin state:** the Perl reservoir/baseflow
  recurrence occurs in the reporting block.  Moving it outside reporting,
  including for `final-only` mode, is a behavior change that needs its own
  compatibility decision and tests.

## Delivery decision: one build, then verification

The JAX engine will be built in one isolated implementation packet rather than
one flux module per implementation packet.  The build will create the complete
package layout above, all flux functions, state and forcing contracts, the
ordered engine, legacy I/O, a runnable CLI, and the replay-radiation provider.
It must leave the R/Perl/GRASS engine unchanged.

Reference capture is the only prerequisite: freeze hashes of Mica Creek inputs,
masks, weather, selected GRASS intermediate rasters, and final outputs, plus
small fixtures for NULLs, land cover, thresholds, directional routing, and
borders.  This capture establishes the comparison target; it is not an
incremental JAX implementation.

The one-shot build may be code-complete without being scientifically accepted.
It must import, compile, load the frozen bundle, and execute a two-day replay
run.  Corrections found during verification stay within the established module
boundaries unless an explicit architecture decision is recorded.

## Verification packet sequence after the build

1. **JAX-VERIFY-01 — Initialization and surface/snow fluxes.** Verify static
   arrays, masks, initial state, precipitation, canopy, longwave/turbulent heat,
   and snow against the captured cell and substep fixtures.
2. **JAX-VERIFY-02 — Soil vertical fluxes.** Verify allocation, ET,
   percolation, saturation runoff, and road runoff with synthetic threshold
   cases and representative reference cells.
3. **JAX-VERIFY-03 — Lateral routing and basin recurrence.** Verify the literal
   eight-neighbor gather, orientation, zero-flow borders, reservoir volume, and
   baseflow recurrence.
4. **JAX-VERIFY-04 — Integrated CPU parity.** Compare two days, then all 30
   January days, against named fluxes and daily states as well as terminal
   rasters.
5. **JAX-VERIFY-05 — General radiation provider.** Verify native or decomposed
   radiation under changing snow albedo, then repeat end-to-end parity without
   replay radiation fields.
6. **JAX-VERIFY-06 — GPU parity and benchmark.** Run the identical float64
   functions on GPU.  Report compilation, host/device transfer, steady-state
   computation, and full workflow time independently.  Only after parity should
   float32, ensembles, calibration, or autodiff be evaluated.

## Acceptance gates

Each translation packet must keep the R/Perl/GRASS reference runnable and add
small focused tests before refactoring a new flux group.  The integrated gate
requires matching date sequence, raster shape and mask; finite-value agreement
for every tracked state and flux using documented absolute and relative
tolerances; water-balance diagnostics; and explicit reports of any legacy
nonphysical values.  Terminal agreement alone is insufficient.

Cross-device bit identity is not required, but any CPU/GPU difference must be
within predeclared variable-specific tolerances.  January is an important
reference window, not complete coverage; fixtures must also cover dry soil,
rain, accumulation, melt/refreeze, zero-flow cells, all land-cover classes,
and routing edges.

## Non-goals of the first JAX release

- No changes to calibrated scientific equations or units.
- No replacement of the current R spatial or weather preparation pipeline.
- No removal of the Perl/GRASS engine.
- No float32 performance claim before float64 parity.
- No claim that replay radiation is a general native radiation solution.
