---
title: Initialization and two-day surface snow legacy parity
status: active
owner: SMR_R
last_reviewed: 2026-11-12
---

# PKT-SMR-JAX-VERIFY-01A completion evidence

**The VERIFY-01 numerical gate passes for the explicit float64 CPU legacy
compilation profile**, covering initialization and all 24 named surface/snow
comparisons for January 1–2, 2015, on all 68,121 required cells. Fifteen
comparisons are exact; the others differ only at double-precision rounding
scale. This is not a claim of bitwise identity for every energy expression,
optimized execution, arbitrary dates, general radiation or GPU parity.

Implementation is isolated on the packet branch from JAX feature base
`722151c3c38b14d0e961e20695e2da34102f1dcd`; the operator accepted feature
landing on 2026-11-12. No
equation, source fixture, soil/routing module, selected residue policy, original
VERIFY worktree or default branch was changed. The original packet's earlier
failure records remain historical evidence.

## Cause and correction

The ordinary JAX lowering changes operation rounding relative to GRASS's
successive raster expressions. A direct primitive check showed 24,303 active
cells differed after division by 130, with maximum error 2.22e-16. Those small
differences propagate into day-one residual SWE and its strict `SWE > 0`
decision, causing material day-two heat flux differences.

Independent eager NumPy evaluation of the literal surface/snow expressions
reproduced day-one SWE and its decisions exactly. JAX optimization barriers
alone did not. Per-compilation controls isolated the required combination:

| Disabled passes (fast math off) | Day-one snow-decision mismatches | Day-two maximum SWE error (cm) |
|---|---:|---:|
| `algsimp` | 2,005 | 0.09923476034 |
| `fusion` | 1,646 | 0.09834262584 |
| `algsimp,fusion` | 0 | 7.7716e-16 |

`fluxes/legacy_math.py::compile_legacy_cpu` makes the accepted CPU profile
explicit and reusable. It turns off fast math and disables those two passes
for that compilation only, checks float64/CPU inputs, and uses the actual
`engine.run_days` implementation. There is no expected-state injection,
reference-derived branch mask, host NumPy callback, global XLA environment
change or tolerance fitted to these failures.

The ordinary optimized `run_days`/replay CLI remains unchanged and is **not**
the accepted legacy numerical profile. For a legacy calculation requiring
this parity contract, use:

```python
from smr_jax.engine import run_days
from smr_jax.fluxes.legacy_math import compile_legacy_cpu

run_reference_cpu = compile_legacy_cpu(
    run_days, static_argnames=("melt_small_snow_next_day",))
final, history = run_reference_cpu(initial, forcings, static)
```

This is a host compilation entry point, not a wrapper to nest inside another
JIT. The installed JAX must support per-compilation `compiler_options`.
Performance and GPU profiles require their own verification; no benchmark or
cross-device claim is made here.

## Fail-closed verification

`python -m smr_jax.verify_surface` is now an assertion-based gate. It checks:

- All seven frozen artifact sizes/SHA256 values before loading the target.
- Initial soil total/A/B amounts and their sum; canopy storage; snow stores,
  age, temperature and thermal initial constants; canopy/partition constants
  for the loaded land-cover classes.
- Exact agreement of model support with each day's independently captured
  finite saturation footprint and the manifest's 68,121-cell requirement.
  That footprint selects comparison cells; it is never used in model physics.
- Every required field on every required cell: missing diagnostics, nonfinite
  values, empty/changed support, wrong shapes or loss of float64 fail.
- Exact day-one `SWE > 0` decisions, even if a changed value lies below the
  ordinary numerical tolerance. This protects the transition into day two.
- Exact replayed shortwave and explicit float64 numerical budgets for other
  fields: water absolute tolerance `1e-12 cm`, energy absolute tolerance
  `1e-9` in inherited diagnostic units, plus relative tolerance
  `512 * eps(float64)` (approximately 1.14e-13). These fixed budgets represent
  numerical operation accumulation, not scientific acceptance thresholds.

The operator's `1e-3 cm` next-day melt convention is a separate model choice,
not this numerical comparison tolerance. The unmodified frozen trace used
strict zero. Therefore the changed convention requires
`--report-only --melt-small-snow-next-day`; that output explicitly declares
no numerical acceptance. The chosen convention and all its regression tests
remain available unchanged.

The raster reader now rejects fractional-cell shifts, resolution mismatches
and missing maps, accounts for nonzero tiepoint pixels, and translates encoded
GDAL nodata into NaN. It no longer silently rounds a half-cell shift or treats
a numeric nodata sentinel as valid evidence.

## Results on 2026-11-12

Windows CPU, JAX `0.12.0.1`, float64. The existing plant virtual environment
was used read-only with the isolated temporary tifffile dependency; no plant
source or environment was installed into or modified.

- Baseline: 41 inherited tests passed.
- The new frozen parity regression initially failed on SWE, latent heat,
  sensible heat and total energy on day two, reproducing the defect.
- Final full suite: **57 tests passed in 36.159 seconds**.
- Strict verification CLI: exit 0, initialization PASS, all 24 comparisons
  PASS on 68,121 cells apiece. All seven frozen inventory hashes match.
- Negative tests exercise corrupt manifests, missing diagnostics/cells,
  nonfinite values, sub-tolerance changes to the snow decision, reader geometry
  and nodata, a genuinely wrong flux, and nonzero CLI failure status.
- Synthetic checks cover the rain/snow thresholds, canopy water/PET ledgers,
  and independently prescribed melt/refreeze amounts. All ten previously
  added next-day residue-policy tests still pass.

| Quantity | January 1 max absolute error | January 2 max absolute error |
|---|---:|---:|
| Precipitation/rain/snow (cm) | 0 | 0 |
| Throughfall (cm) | 0 | 1.3878e-17 |
| Canopy evaporation (cm) | 0 | 0 |
| SWE (cm) | 0 | 7.7716e-16 |
| Snowmelt (cm) | 0 | 1.1102e-16 |
| Replayed shortwave | 0 | 0 |
| Longwave energy | 7.2760e-12 | 3.2742e-11 |
| Latent energy | 0 | 2.2510e-11 |
| Sensible energy | 0 | 1.4552e-11 |
| Total energy | 7.2760e-12 | 4.9113e-11 |

The strict CLI regenerates the complete per-map RMSE/max table. No fields or
cells were dropped, and no frozen expected values were rewritten. No held-out
observations, productivity parameters or evaluation splits were used.

## EOS scope and closure handoff

Live toolkit create/preflight passed before edits on the isolated packet
branch. This packet owns only the numerical profile, reader/verifier, tests
and local evidence; report-only scope passes with seven in-scope paths.
No strict check was added to a Git hook or commit/push path.

The adjacent completed packet definition explicitly names `PKT-SMR-JAX-VERIFY-01`
as the gate it closes. The local parent record now records the passing result
on this branch. The operator accepted feature landing on 2026-11-12; the EOS
run record records the actual landing commit after integration. Default-branch
publication remains a separate action.

The accepted landing reconciles the parent gate's operational completion record
with `PKT-SMR-JAX-VERIFY-01A`, the feature landing and this evidence as provenance.
The archived `engineering-os` source definition's historical active status is
preserved in its read-only repository; it is not the current completion evidence.
No global board or registry identity change is needed for this numerical result.
Soil (VERIFY-02) is the next verification gate; routing/basin, 30-day integration,
general radiation and GPU remain open. No worktree removal or push is included.

## Integrated initialization follow-up

VERIFY-04 subsequently compared the integrated engine with actual GRASS daily
states. It found that the original source-only initialization check had treated
`fieldcap_mc_A - wiltpt_mc_A` as a double operation. Both operands are imported
FCELL rasters, so GRASS rounds this difference before multiplying by the DCELL
constant `0.2`. Both initialization entry points and the initialization checker
now preserve that source arithmetic. The surface/snow gate remains passing.
See [VERIFY-04](PKT-SMR-JAX-VERIFY-04.md) for the two-day pass and the separate,
still-open thirty-day strict-parity result. Earlier evidence above is historical.
