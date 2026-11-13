# PKT-SMR-JAX-VERIFY-04 — Integrated CPU verification

## Result: strict default pure-JAX month parity deferred by operator

**Operator decision, 2026-11-13:** "THRESHOLDS ARE TINY MOVE ON" authorizes
proceeding to VERIFY-05 while deferring strict VERIFY-04 closure. The measured
default-JAX maximum differences are `2.2823668732740998e-8 cm` in soil storage
and `4.427795019523728e-8 cm` in percolation, below `4.5e-7 mm` of water.
The 167 failing strict checks remain recorded. Existing budgets and nonzero
strict-gate failure exits are unchanged; neither default JAX nor native FFI is
relabelled as a default-JAX strict pass. This explicit sequencing decision
replaces the requirement to hold VERIFY-05 pending strict closure below.

**Status correction, 2026-11-13, PKT-SMR-JAX-VERIFY-04A:** the completion claim
landed in `0c0152e` was incorrect. The requested default pure-JAX engine has not
passed strict thirty-day parity. Passing an optional native CPU exponential
profile does not close that requirement. VERIFY-04 is **active and incomplete**;
VERIFY-05 must not advance on the withdrawn completion claim.

The code changes landed on the JAX feature branch; this is not evidence that
the default JAX path achieved one-to-one parity. Keep the native profile as
diagnostic evidence identifying arithmetic causes. Closure requires the default
pure-JAX continuous run to pass the existing two-day and thirty-day criteria,
without replacing its exponential with a native FFI or relaxing the budgets.

## Native CPU diagnostic result (not VERIFY-04 acceptance)

The optional native CPU profile passes all **960 daily field comparisons**, full
daily reporting, actual terminal carry and both original golden rasters, with
the original tolerances and all **68,121 required cells**. Every daily soil,
snow, thermal and flux field matches exactly; 31 of the 32 field families have
zero maximum error. The remaining annual runoff counter differs by at most
`1.3322676295501878e-15 cm`. No tolerance, source physics, six-hour step, input,
reference artifact, support mask or threshold was changed.

This diagnostic uses an **explicit Windows x64 MinGW scalar-exp FFI profile**.
The ordinary pure JAX exponential remains the model default and continues to
fail the strict month gate. There is no default JAX month-parity acceptance.
The new native library is optional, built outside Git, and needs no GRASS
runtime. General radiation and GPU verification remain VERIFY-05 and VERIFY-06.
The original two-day gate also passes with the default JAX math profile.

See the native diagnostic and reproduction evidence below. The earlier failure
evidence is retained as a historical record rather than recast as a pass.

## Historical audit before the source arithmetic corrections

At packet revision `b17dcc6`, the actual JAX daily scan passed all **40 existing two-day daily-map comparisons**
and the CPU report checks. The full **30-day, 960-map-comparison audit** is now
available and failed the unchanged strict numerical budget. The packet remained
active at that revision, with accepted feature base `47b517e`.

The thirty-day audit has 173 failed checks: 167 daily map comparisons plus a
reporting check and five terminal checks. Its first failing daily field is
**January 5 total soil storage** (17 of 68,121 cells beyond budget). Maximum
differences over the entire month are:

| Quantity | Maximum absolute difference |
| --- | ---: |
| Total storage | 2.82600165490976e-8 cm |
| A storage | 2.3085569456071653e-8 cm |
| B storage | 2.82600165490976e-8 cm |
| Percolation | 4.427795019523728e-8 cm |
| Actual ET | 2.0274703027567398e-9 cm |
| Last-substep lateral flow | 2.21389745442283e-8 cm |
| Runoff | 3.0531133177191805e-15 cm |
| Saturation | 1.5716224766038067e-7 percentage points |
| SWE | 7.105427357601002e-15 cm |
| Total heat flux | 4.729372449219227e-11 inherited heat units |
| Terminal storage versus original golden raster | 2.0096048380935372e-8 cm |

The largest water-flux difference is `4.427795019523728e-7 mm`. Its physical
smallness does not make it a pass against the declared `1e-12 cm` absolute
floor and float64 relative budget. No numerical tolerance was relaxed, no cells
were dropped and no reference states or branch masks were inserted into the
actual model evolution.

## Independent thirty-day reference

`tools/capture_jax_integrated_reference.py` ran the frozen Perl engine through
all January 1–30 dates in a new temporary GRASS 8.5.0 EPSG:32611 project. It
imports the original hashed ASCII inputs. A temporary Perl copy adds output
exports only: the scientific expressions, weather sequence, original mask,
four six-hour steps and full-report recurrence are unchanged. It also captures
initial total/A storage and first-day soil checkpoints for diagnosis.

The resulting **335,131,312-byte** lossless archive contains all original 20
daily maps plus 12 additional state/flux maps for every day. The large generated
files are retained outside Git at:

```text
C:\projects\git\worktrees\smr-verify04-reference\
  manifest.json
  daily_trace.zip
  daily.csv
  capture-engine.pl
  integrated-audit.json
  integrated-audit.log
```

The versioned witness is
`tests/fixtures/jax_verify04_capture.json`. It records the original-source,
instrumented-source, capture-tool, seven frozen-input and generated-artifact
hashes, sizes, GRASS build and map inventory. The capture tool refuses an
existing output directory. Original frozen files and `R/smr_buildup.pl` are
unchanged. Existing GRASS runtimes were not used or modified.

Before any new shortwave field is prescribed, the extension is anchored against
every original two-day map over its entire finite footprint, both original
terminal rasters and all 30 original legacy reports. The original zonal
reporting path rounded sums through FCELL and six decimal places, whereas the
optimized CPU trace retains double sums. The validator explicitly reproduces
that presentation and the corresponding legacy reservoir recurrence; it does
not equate the differently formatted rows or widen the state tolerances.
All reference-extension anchors passed before the JAX month comparison.

Only shortwave is prescribed in model evolution. All remaining surface,
snow/thermal, canopy, soil, runoff, routing, annual and basin quantities come
from the actual JAX scan. Native/general radiation and GPU acceptance remain
VERIFY-05 and VERIFY-06.

## Demonstrated corrections and instrumentation

1. The original two-day integrated run differed in A/B storage by about
   `2.24e-8 cm`. GRASS rounds the initialization moisture difference as FCELL
   before the `0.2` multiplication. Correcting both initialization entry points
   and the checker removed that discrepancy; two-day maximum soil-water error
   is now below `2e-15 cm`.
2. The source longwave coefficient is `9.2*(10^-6.0)`, which differs by one
   float64 rounding unit from the transcribed `9.2e-6`. The literal source
   coefficient is preserved. First-day longwave is now checked for exact
   agreement over the full captured footprint.
3. Six-hour rate arithmetic retains the source `*6/24.0` ordering, and the four
   temperatures are accumulated sequentially as in Perl. These choices alone
   did not resolve the wet-period mismatch.
4. The source annual runoff counter accumulates soil runoff before road runoff
   is added to the daily report. The engine previously included the road term
   in this counter. It now preserves the source distinction, with a warm/rainy
   synthetic road case checking that behavior.
5. Additive diagnostics expose the last lateral substep separately from the
   existing daily lateral sum, and expose road-adjusted soil input, canopy and
   complete snow/thermal state, annual counters and basin carry. Original
   diagnostics remain available. These are measured model values, not values
   copied from the reference.

## Required support and strict failures

Every day requires the same independently captured 68,121-cell soil footprint.
Each map additionally compares its full original finite footprint: surface
reporting maps can cover more cells than the final soil state. Every required
soil cell must remain present in every compared field. Reporting sums use
independent per-field footprints rather than silently using the soil mask for
all quantities. Condensation support is derived from the source conditional's
captured snow, previous SWE and latent-flux operands, not from actual finite
values. The original two-day fixture lacks canopy-storage/rain-heat rasters;
those extra reporting fields are checked in the explicitly extended month run.

Comparison budgets are unchanged: water `atol=1e-12 cm`, heat `1e-9`,
saturation/temperature `1e-10`, and `rtol=512*eps(float64)`. Shortwave, snow age
and saturation-day counts are exact. SWE must also preserve exact snow-presence
decisions. Aggregate floors sum the existing per-cell floors; basin flows use
`atol=1e-15` in inherited cms units. The full actual terminal carry is checked,
not just its last history entry.

The default CLI fails at the first mismatch and names its date. `--audit-all`
collects numerical disagreements while retaining exit code **1**. Missing
support, missing fields, nonfinite values and malformed/checksum-invalid
evidence still fail. `--report-path` requires a new file and preserves existing
reports. A diagnostic report cannot be mistaken for a successful acceptance.

## Investigation at the initial failing revision

All surface/snow/thermal maps pass the month budget, including exact snow
presence. The failures are small soil-state/transport differences accumulating
from wet conditions. A separate original-GRASS soil trace, seeded with actual
January 4 states and January 5 surface input, agrees with the same JAX stages
to about `2e-15 cm`. The continuously evolved JAX state instead diverges at the
strict threshold budget. This supports threshold-sensitive accumulated
rounding rather than a large error in the isolated routing equation.

An explicit January 5 fourth-substep branch witness occurs at cell `(74,79)`:
JAX B is `-0x1.ac00220000001p-1`, while GRASS B and the frozen FC value are
`-0x1.ac00220000000p-1`. They select different conductivity arms, producing
six-hour effective values `5.103468288619541e-6` and
`5.113047980945727e-6`. The frozen case includes such signed legacy properties;
they were preserved. This demonstrates sensitivity but does not claim to have
localized every one of the first 17 failing cells to that exact branch.

Compiler-only vector-width (1/64) and platform-dependent-math experiments did
not remove the mismatch. A floating-power experiment also did not resolve it
and is not retained. No NumPy callbacks, expected-value substitutions, FC snaps,
new infiltration physics or threshold smoothing were added.

At that revision, closure required a source-faithful arithmetic treatment or an
operator decision changing acceptance or conductivity behavior. The operator
selected continued strict verification. The measured maxima above were retained
as evidence, not fitted replacement tolerances.

## Native diagnostic on 2026-11-13; closure claim withdrawn

Two source arithmetic details jointly produced the native diagnostic pass:

1. The snow surface temperature expression is `u.surface / 24.0 / capacity`.
   The translation grouped its denominator as `u.surface / (24.0 * capacity)`.
   Restoring the sequential divisions preserves the original rounding and
   first-day internal temperature exactly. A regression checks this distinction
   on the actual integrated state, including cells where regrouping differs.
2. XLA's CPU float64 exponential and the frozen GRASS executable's scalar
   exponential can differ by one rounding unit. The source conductivity
   threshold can amplify that difference during continuous wet-period
   evolution. The optional `legacy_cpu_math` profile binds the same model
   functions to a typed JAX FFI exponential, leaving the default `jnp.exp`
   binding available. All other state evolution remains in the actual JAX
   daily/four-substep scans.

The C++ handler applies native MinGW scalar `exp` directly to XLA float64 CPU
buffers. Its x87 precision guard restores the calling thread's control word
after each call. Native C starts at 64-bit x87 precision, while CPython starts
at 53-bit precision; an early Python/GRASS-DLL diagnostic without this guard
was therefore not an equivalent test of the executable. A separate 100,001-point
comparison over `[-700,700]` found the installed compiler's native exponential
and GRASS's bundled exponential identical under the native precision setting.
The model does not depend on the GRASS DLLs.

The FFI path contains no Python/NumPy callback, copied reference state, FC snap,
threshold smoothing, alternative infiltration equation or global XLA flag.
JIT, batching and an elementwise exponential JVP are supported. The builder
requires an explicit compiler/output, refuses existing library/build-record
paths and records source/library hashes, JAX version, compiler version and the
exact command. Loading verifies these hashes and JAX version, retains the DLL
for executable lifetimes and refuses changing the library within a process.
Other compiler/platform profiles require their own successful strict gate.

The implementation follows the public [JAX typed FFI API](https://docs.jax.dev/en/latest/ffi.html).
The distinction from XLA's intrinsic exponential is visible in the
[CPU elemental emitter](https://github.com/openxla/xla/blob/main/xla/service/cpu/elemental_ir_emitter.cc)
and [intrinsic implementation](https://github.com/openxla/xla/blob/main/xla/codegen/intrinsic/exp.cc).

The native diagnostic full audit is retained outside Git as
`integrated-native-ffi-v2-audit.json`. It exits **0** with zero failed checks.
The report includes the selected math profile and build provenance:

- Profile: `windows-mingw64-scalar-exp`, JAX `0.12.0.1`.
- Compiler: MinGW-W64 x86_64-ucrt-posix-seh GCC `13.2.0`, r8.
- C++ source SHA-256: `b456ae671ac5ec98ab00fa52090df492548029b9f9bb1ed09e66d84a3b69b32d`.
- Library SHA-256: `63a180c1d413d143747a0ad7c40f9f6596fa5944277af808db9b5cea78abad4b`.
- Largest water-map difference: `1.3322676295501878e-15 cm`, annual runoff only.
- Largest daily reporting sum difference: `4.291534423828125e-6` in inherited
  aggregate units, within the unchanged field-specific aggregate budget.

The full suite passed **108 tests in 118.191 seconds** with the native option
enabled. After adding a binding-immutability regression and defensively copying
returned build metadata, all **8 native boundary/execution checks passed in
12.726 seconds**, including the added test and the original two-day gate. These
runs cover 109 distinct tests. Report-only scope and `git diff --check` passed.

Build a NEW library outside Git to reproduce the native diagnostic:

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PYTHONPATH="$PWD\src;$env:TEMP\smr-jax-plant-testdeps"
# Choose a NEW output name; this is the already-built diagnostic artifact:
$mathLibrary='C:\projects\git\worktrees\smr-verify04-reference\smr_legacy_exp_v2.dll'
# To reproduce the build, use a different new output path:
# python -B -m smr_jax.legacy_cpu_math --compiler C:\Strawberry\c\bin\g++.exe --output <new-library.dll>
& C:\projects\git\plant-game-of-life\.venv\Scripts\python.exe -B -m smr_jax.verify_integrated --days 30 --capture-dir C:\projects\git\worktrees\smr-verify04-reference --native-math-library $mathLibrary
$env:SMR_JAX_NATIVE_MATH_LIBRARY=$mathLibrary
& C:\projects\git\plant-game-of-life\.venv\Scripts\python.exe -B -m unittest discover -s tests -p 'test_jax_*.py' -q
```

For application use, load this library explicitly with
`smr_jax.legacy_cpu_math.load_library(path)` and bind
`compile_native_cpu(run_days, static_argnames=('melt_small_snow_next_day',))`.
The ordinary `run_days` default and GPU path still use JAX math. This packet
does not claim GPU numerical parity or qualify a different math library.

## Initial validation before the native diagnostic

- Preflight and report-only scope checks passed from feature base `47b517e`.
- All **100 JAX tests passed in 141.913 seconds**, including the earlier gates,
  selected snow policy, plant interface, integrated two-day checks, deliberate
  corruption, terminal-state/history disagreement and audit-failure exits.
- After applying the same initialization correction to the foundation entry
  point, its three checks passed in 0.625 seconds; the final two-day CLI also
  passed and its JSON report is retained as `integrated-two-day.json` beside
  the month audit. The integrated engine's initializer and month numerics did
  not change in that final foundation adjustment.
- The strict 30-day CLI exited **1**, first on January 5 total storage.
- The complete audit evaluated **960 daily fields** and exited **1** with the
  173 failures recorded above. Its first aggregate failure is January 17 storage:
  `1.34285074e-7 cell-cm`, against the unchanged aggregate budget.

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PYTHONPATH="$PWD\src;$env:TEMP\smr-jax-plant-testdeps"
& C:\projects\git\plant-game-of-life\.venv\Scripts\python.exe -B -m smr_jax.verify_integrated --days 2
& C:\projects\git\plant-game-of-life\.venv\Scripts\python.exe -B -m smr_jax.verify_integrated --days 30 --capture-dir C:\projects\git\worktrees\smr-verify04-reference
# Optional complete numerical audit; still exits 1 on failures; output must be NEW:
& C:\projects\git\plant-game-of-life\.venv\Scripts\python.exe -B -m smr_jax.verify_integrated --days 30 --capture-dir C:\projects\git\worktrees\smr-verify04-reference --audit-all --report-path <new-report-path>
```

To rebuild the reference, run the versioned capture tool under a new temporary
GRASS project and supply a NEW output directory outside Git:

```text
grass --tmp-project EPSG:32611 --exec python -B tools/capture_jax_integrated_reference.py --output <new-output-directory>
```
