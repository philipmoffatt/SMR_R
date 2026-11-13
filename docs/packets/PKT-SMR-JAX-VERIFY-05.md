# VERIFY-05 - General radiation with evolving snow albedo

## Result

The decomposed radiation provider passes its source-precision checks. Actual
two-day and thirty-day runs use evolving JAX albedo without completed shortwave
replay, preserve required support and retain the previously deferred tiny water
differences. Strict integrated legacy parity is explicitly **not passed**.

## Boundary and operator decision

The operator explicitly authorized moving on from the tiny VERIFY-04 water
differences on 2026-11-13. Strict default-JAX month parity remains deferred and
its checker still fails. This packet implements the planned **GRASS radiation
decomposition**, not a native JAX replacement for r.sun. No native FFI or Python
callback is used by the new radiation path.

Preparation computes beam, diffuse and reflected irradiation at unit ground
albedo for the requested terrain, dates, cloud and Linke turbidity. During each
actual JAX daily step, snow age and albedo are updated first; shortwave is then
computed from those components and the evolving albedo. Completed `q.srad` maps
and reference snow states are not model inputs. The existing replay provider
remains available for reference diagnostics.

## Equations and source precision

GRASS's [r.sun source](https://github.com/OSGeo/grass/blob/main/raster/r.sun/rsunlib.c)
applies local ground albedo to reflected irradiation. Beam and diffuse
irradiation are independent of its value on a fixed support mask. The
[r.sun input/output code](https://github.com/OSGeo/grass/blob/main/raster/r.sun/main.c)
reads albedo and writes these radiation components as FCELL. The
[manual](https://grass.osgeo.org/grass-stable/manuals/r.sun.html) specifies daily
component units as Wh/m2/day.

The JAX kernel preserves those boundaries:

```text
reflected = FCELL(FCELL(model_albedo) * reflected_at_unit_albedo)
incident  = FCELL(FCELL(beam + diffuse) + reflected)
q.srad    = (1 - canopy_cover) * (1 - model_albedo) * incident * 3600 / 1000
```

The outside absorption factor retains DCELL model albedo. Resulting shortwave
is in kJ/m2/day. The unit reflected basis has already been rounded by GRASS,
so this decomposition has a declared float32 quantization error rather than
exact equality to a newly evaluated r.sun at every possible albedo.

Before comparing values, the radiation budget was defined as eight float32
epsilons times the absolute component sum and absorption/unit-conversion
factor, plus the existing 512-float64-epsilon arithmetic allowance. This is a
conservative forward-error bound for the rounded basis, scaling and additions;
it was not fitted to the failed reference comparisons. It is a radiation
budget, not a replacement soil/snow acceptance tolerance.

## Terrain support correction

The first preparation used constant unit albedo across the full terrain.
An independent full-grid comparison rejected it: maximum January 1 shortwave
error was about 286 kJ/m2/day, far beyond the declared bound. GRASS removes
NULL-albedo cells from its shadow-casting elevation array, so an extra valid
albedo cell can change a neighboring cell's beam irradiation.

The production preparation now constructs unit albedo only where the source
elevation and land use are finite. r.sun applies the remaining terrain/input
NULL rules. This mask is derived from static inputs, not captured state. Cache
schema 2 names this convention and rejects the former maskless cache. The
corrected January 1 maximum is `0.00011500488290039357 kJ/m2/day`, within the
original radiation bound. All required model cells are preserved.

## Evidence

The small independent fixture contains **60 cases x 128 representative cells**:
days 1, 80, 172 and 355; cloud 0, 0.5 and 1; albedos 0, 0.2, 0.7383, 0.95 and 1.
Each albedo is independently evaluated by r.sun; beam/diffuse independence is
asserted during capture. Maximum absorbed-shortwave error is
`0.00009200390627484012 kJ/m2/day`, using at most 12.0% of the declared bound.
The fixture preserves its original full-region albedo support as an algebraic
test. Its exact original capture/preparation scripts are retained in packet
commit `96138d3`; full-model static-support verification is a separate check.

The source-mask-corrected January cache is independently anchored at the
reference's albedos over **2,089,830 full-domain cell-days**. Maximum shortwave
error is `0.0004811334983969573 kJ/m2/day`, using at most 12.45% of the declared
bound. Reference states are used only in that comparison, never injected into
the actual simulation.

The actual 30-day default legacy-snow run uses only prepared components.
Its maximum storage and percolation differences remain respectively
`2.2823668732740998e-8 cm` and `4.427795019523728e-8 cm`, the same maxima already
deferred in VERIFY-04. Daily snow-presence and annual saturation-count
decisions match. Maximum SWE difference is `5.573519423762718e-9 cm` and maximum
surface-temperature difference is `1.9092599123915477e-7 C`.

Strict radiation/thermal and other legacy comparisons remain visible in the
JSON audit. Source-precision component acceptance is reported separately from
strict integrated parity. Structural failures, incomplete support and corrupt
inputs still fail verification. The optional next-day residue mode remains a
distinct operator-selected behavior; comparisons with that mode use a separate
same-policy replay control to distinguish policy effects from radiation error.

Over thirty days, basin discharge/baseflow differ by at most
`4.667761403093218e-16 cms`, and terminal reservoir volume by
`2.278559285695536e-15` in inherited units. These pass the existing strict basin
budgets. Radiation/thermal reporting sums and several terminal fields exceed
their old strict budgets and are reported as failures, not hidden by the
component pass.

The selected next-day residue policy deliberately differs from legacy snow
physics: its maximum SWE difference from the legacy reference is
`0.13759970718978032 cm`. With that policy held constant in a separate replay
control, the radiation-provider change contributes at most
`5.573519423762718e-9 cm` SWE, `2.1394686466891244e-8 cm` soil storage and
`4.183503365950969e-8 cm` percolation on the required 68,121-cell domain. Policy
effects are therefore not attributed to the new provider. The default legacy
snow mode is unchanged.

The full suite passed **119 tests in 170.524 seconds**. After adding a negative
native-evidence check and tightening malformed/missing cache cases, all **11
radiation tests passed in 8.749 seconds**. These runs cover 120 distinct tests.
The runnable CLI completed a prepared-radiation two-day smoke test. Report-only
scope and `git diff --check` passed.
After final import cleanup, the 11 radiation tests passed again in 12.485
seconds and all 16 surface/replay tests passed in 36.316 seconds.

## Preparation and use

Create a new disposable GRASS project in the source terrain CRS. For the tested
Mica Creek case that is EPSG:32611:

```powershell
C:\OSGeo4W\bin\grass85.bat --tmp-project EPSG:32611 --exec C:\OSGeo4W\bin\python3.exe -B tools/prepare_jax_radiation.py --bundle reference/mica-creek-jan2015/model_inputs.zip --weather reference/mica-creek-jan2015/gridmet_weather_2015-01.csv --output <NEW-output-directory-outside-Git>
```

The tool refuses an existing output directory or populated GRASS mapset.
It preserves r.sun's terrain shadowing, default half-hour integration and
literal legacy cloud equations, including diffuse coefficients above one.
The cache records projection, grid, source bundle/weather hashes, preparation
script hash, component hash and each date's cloud/Linke/day-of-year inputs.
The provider verifies grid, bundle, component inventory/hash, FCELL type,
nonnegative finite values, required support and forcing identity. Grid-origin
serialization uses the existing TIFF reader's `1e-8`-cell alignment floor.
Changed cloud, turbidity, dates, terrain or land use require a new cache.
Weather fields that do not affect radiation may change without rebuilding it.

Run the model using the prepared cache:

```powershell
python -B -m smr_jax.cli --days 30 --radiation-dir <prepared-directory> --summary <summary.json>
# Optional existing selected snow convention:
python -B -m smr_jax.cli --days 30 --radiation-dir <prepared-directory> --melt-small-snow-next-day --summary <summary.json>
```

`--weather-file` accepts another legacy-format forcing file. Library users can
construct `PreparedRadiationProvider(directory, static.grid, bundle_path,
required_mask=static.active_mask)` and pass it to `load_daily_forcings`. The
result contains components and explicitly absent replay shortwave. Supplying
both providers, missing required days, or missing both providers fails.

GRASS is required for preparation. An already prepared model run uses JAX
arrays only; GPU verification remains VERIFY-06. Other terrain/CRS cases and
held-out observations have not been used for tuning or claimed as validated.

## Reproduction and external artifacts

```powershell
python -B -m smr_jax.verify_radiation
python -B -m smr_jax.verify_radiation --radiation-dir C:\projects\git\worktrees\smr-verify05-radiation-masked --capture-dir C:\projects\git\worktrees\smr-verify04-reference --days 30 --report-path <NEW-report.json>
```

The independently prepared cache and large reports stay outside Git at
`C:\projects\git\worktrees\smr-verify05-radiation-masked`. The original maskless
preparation and its rejected comparison are retained as diagnostics at
`smr-verify05-radiation`; neither is accepted by the schema-2 provider.
The small versioned radiation fixture records source, bundle, capture-tool
and binary artifact hashes. Original Perl and frozen reference files are
unchanged. No board, plant repository, default branch or shared GRASS project
was modified.

Accepted prepared-cache witnesses:

- Component archive: `21,981,988 bytes`, SHA-256
  `e4bf0394008ff425f6c30e082490cf619d8b20605f28cafd415e1ca3e845dce9`.
- Preparation source SHA-256:
  `c3b806cf2642af7a49e3ab02863b32a1134820b0fd5db3adf7228d625f0c6698`.
- Final reports: `two-day-final.json`, `thirty-day-final.json`,
  `selected-snow-control-final.json` in the prepared-cache directory.
