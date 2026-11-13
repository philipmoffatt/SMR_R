# SMR-JAX feature closure and PR boundary

The agreed SMR-side prototype feature is complete under the operator-selected
practical acceptance profile. This record reconciles current release status;
older packet records remain historical evidence of their execution boundary.
Completion does not claim plant-side integration, a new Richards solver,
arbitrary soil layers, general native r.sun replacement or field validation.

## Delivered and accepted

- Float64 two-layer daily engine with surface/snow, soil, legacy receiver-weighted
  routing and basin recurrence, plus frozen source verification tools.
- Bounded external interception, per-layer transpiration and soil evaporation;
  actual/unmet extraction, supply-aware refusal and independent water ledgers.
- Adaptive transport substeps for the external plant-water path; legacy ordering
  remains separately selectable. Layer-boundary sensitivity is documented.
- Explicit optional next-day release of positive SWE <=0.001 cm; the literal
  legacy policy remains available and is not mislabeled as the selected policy.
- Prepared terrain/weather radiation components with evolving model snow albedo.
  GRASS is required for preparation; completed shortwave maps are not inputs to
  the prepared-provider simulation.
- Actual RTX 3070 Ti CUDA float64 verification and synchronized benchmark
  evidence, with strict diagnostics retained alongside practical acceptance.

VERIFY-01/02/03 source-stage evidence is recorded in their packet notes.
The operator explicitly deferred tiny strict default-JAX month differences
after VERIFY-04 and rejected strict parity as the GPU release gate. VERIFY-05
passes its declared source-precision radiation checks. VERIFY-06 passes
practical acceptance for two and thirty days on all 68,121 active cells,
complete histories, terminal state and basin scalars. It preserves exact
mask/counter/finite-footprint checks and independent conservation limits.

Accepted maximum month water difference is 0.00017706047951770998 cm against
the explicit 0.001 cm water limit. Strict arithmetic identity remains a
diagnostic failure. GPU verification uses explicit compiler settings preserving
source FCELL rounding; arbitrary default JIT configurations are not qualified.
See [VERIFY-06](PKT-SMR-JAX-VERIFY-06.md) for all unit limits and source/data hashes.

## Final validation and provenance

- Implementation feature revision: `f30fc7a56cc87b9db0c3b766d6fbe6f41402c350`.
- Exact executed GPU source: `089c1eea1b0183aabf2ba59fbe00c5fbee859e62`.
- Final feature regression: 137 tests, successful, with five optional native
  CPU-library tests skipped. No native FFI is used by the accepted GPU path.
- Final harness tests: 17 successful on Windows CPU and CUDA-enabled WSL;
  CUDA-only foundation tests: 3 successful. Independent plant-water ledgers
  pass on both devices for ten cases.
- Accepted month warmed execution: CPU 2.02236 s / GPU 0.48445 s (4.17x).
  Measured prepared-input workflow with full download: 3.29849 s / 2.59985 s
  (1.27x). These are prototype observations on the recorded hardware.
- GitHub's actual default/PR base is `master`, not a branch named `main`.
  Remote base checked: `27774fc095da63d33adae694d669e0a7b894a5e4`, CPU PR #10.
  A three-way `git merge-tree --write-tree` check is clean; the default branch
  is not advanced by this release packet.

This packet changes only current README guidance and its own Markdown/YAML.
It leaves implementation, reference bytes, dependencies, historical packet
records, global boards, other repositories and existing uncommitted work intact.
Preflight passed; report-only scope covers the three declared paths.

## Separate follow-up work

- Plant repository: adopt actual supplied uptake before committing growth,
  replace estimated flow accounting and review the documented root/growth issues.
- New scientific backend: arbitrary soil layers, Green-Ampt/Richards decisions,
  restriction/bottom-boundary laws and any broader convergence study.
- Complete automatic surface/snow/basin integration for the external plant
  interface; the delivered boundary accepts explicit surface/snow bookkeeping.
- Additional terrain, forcing, GPU and field evaluation. No held-out observation
  was used to tune equations or acceptance, and no evaluation split changed.

The feature can be reviewed and merged as this bounded prototype update while
those follow-ups remain explicit. Historical strict failures are not rewritten
as passes, and no plant physiology or soil backend replacement is implied.
