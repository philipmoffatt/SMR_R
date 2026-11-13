---
title: Frozen surface trace tooling preservation
status: active
owner: SMR_R
last_reviewed: 2026-11-12
---

# PKT-SMR-JAX-VERIFY-01 preservation evidence

The operator requested a focused commit of the existing uncommitted trace
reader, comparison reporter and alignment test. Source behavior is preserved
verbatim; only this evidence record is new. Base:
`4b5733ec1a9da2a86572975f3b1e1564cdbd8616`.

The patch extracts `FrozenTraceProvider` from the existing radiation replay
reader and keeps `ReplayRadiationProvider.shortwave` as a wrapper. It adds
`verify_surface.py` to report named two-day errors and
`test_jax_trace_contract.py` to check trace raster alignment. It changes no
flux equations, source fixtures, science choices or acceptance tolerances.

Validation on Windows CPU, 2026-11-12:

- Live EOS preflight passed, with one warning for the three existing dirty
  paths; no writes preceded that preflight in this pass.
- All 9 JAX tests passed in 6.587 seconds using the existing plant JAX
  interpreter and isolated temporary tifffile dependency, with bytecode writes
  disabled. The plant interpreter/environment was not modified.
- `python -m smr_jax.verify_surface` ran successfully and reported all 24
  date/map comparisons on 68,121 common finite cells apiece.

The numerical verification gate is **still open**. Day one water fields agree
to approximately 2.22e-16 cm maximum error; day one total energy max error is
1.4552e-11 in the inherited diagnostic units. Day two has material differences:

| Map, 2015-01-02 | RMSE | Maximum absolute error |
|---|---:|---:|
| SWE (cm) | 0.0136112959859 | 0.0964109100625 |
| q.latent (inherited energy units) | 385.880241200 | 2733.24930027 |
| q.sensible (inherited energy units) | 2207.56637903 | 16723.2041799 |
| q.total (inherited energy units) | 2572.57087660 | 17575.4683617 |

The reporter uses the finite intersection of actual/expected arrays and has no
acceptance thresholds or assertions. Its zero exit code and the alignment test
are not a surface/snow parity pass. Missing/mismatched comparison support must
be investigated alongside the day-two differences before accepting this gate.
The complete reporter can reproduce the per-map table; no tolerance was fitted
to these discrepancies.

This preservation commit stays on the VERIFY packet branch. It is excluded
from the finished JAX feature landing until its gate is resolved or an explicit
tooling-only integration is separately accepted. No claim is made about soil,
routing, reservoir, 30-day, general radiation or GPU parity. No plant files,
frozen inputs or other worktrees were edited by this preservation pass.

## Integration update — 2026-11-12

The operator accepted the `PKT-SMR-JAX-SNOW-01` feature landing, which carries
these preserved trace tools alongside the explicit next-day residual-melt
policy. The earlier exclusion statement describes the preservation handoff;
tooling integration is now explicitly accepted through that packet. The
surface/snow numerical verification gate remains open, including comparison
under the agreed convention. This note does not declare legacy parity.

## Numerical gate closure — 2026-11-12

`PKT-SMR-JAX-VERIFY-01A` establishes passing initialization and all 24 named
daily surface/snow comparisons on the full 68,121-cell support, using an
explicit strict float64 CPU compilation profile. It resolves the day-one
rounding/branch mismatch without modifying equations or the frozen reference.
The source of the discrepancy and per-field errors are recorded in
`PKT-SMR-JAX-VERIFY-01A.md`; all 57 tests pass.

This closes the numerical VERIFY-01 gate for the frozen legacy two-day CPU
target on that packet branch, pending feature integration. It does not accept
the ordinary optimized profile or claim bitwise equality for every energy
operation. The selected next-day residue-melt convention remains separately
tested and explicitly report-only against the unchanged legacy trace.
Subsequent soil, routing/reservoir, 30-day, radiation and GPU gates remain open.
Earlier open-gate statements above describe the preserved history.
