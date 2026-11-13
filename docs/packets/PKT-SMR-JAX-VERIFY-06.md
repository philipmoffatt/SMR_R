# VERIFY-06 - CUDA float64 verification and benchmark audit

## Accepted result

VERIFY-06 is complete under the operator-authorized practical profile.
Both two-day and thirty-day CUDA runs pass, alongside both strict synthetic
full-model cases and all ten strict plant-water/independent ledger cases.
The original strict diagnostic still fails; strict arithmetic parity is not
claimed. The existing selected `<=0.001 cm` next-day snow-release policy and
explicit precision-preserving compiler settings define this accepted run.

Maximum active month water difference is `0.00017706047951770998 cm`, or
17.7% of the declared `0.001 cm` practical limit. No active cell, day, terminal
field or basin scalar is omitted. Masks, counters and finite footprints match.

| Accepted practical benchmark | Two-day CPU | Two-day GPU | Month CPU | Month GPU |
|---|---:|---:|---:|---:|
| Median of five synchronized resident executions, seconds | 0.13194 | 0.04248 | 2.02236 | 0.48445 |
| Measured prepared-input workflow including full download, seconds | 1.02299 | 0.94154 | 3.29849 | 2.59985 |

Observed resident compute ratios are 3.11x for two days and 4.17x for the
month; the measured month workflow ratio is 1.27x. These are qualified by
practical active-watershed acceptance, not strict parity. Full histories are
still downloaded. GPU allocator peaks remain 270,532,608 and 1,608,254,720
bytes. Timing variation and environment limitations from the original audit
still apply; these measurements do not predict other hardware or workloads.

The final 17 harness tests pass on Windows CPU and in the existing CUDA-enabled
WSL environment; CUDA-only foundation tests also pass (3). Acceptance tests
explicitly reject water errors above the limit, counter changes (including
outside the numeric domain), required NaNs and every unit's excessive error.
The original strict comparer and independent conservation checks are unchanged.

Accepted executed source: `089c1eea1b0183aabf2ba59fbe00c5fbee859e62`.
Source archive SHA-256:
`f519dd738e86701bec563048b8ca64aadfb18114d8b2a4852f45b6bf9250e47b`.
The report is `C:\projects\git\worktrees\smr-verify06-snapshot-089c1ee\verify06.json`,
SHA-256 `a5e75a072627b8a2f163a6535f994558c97a121a490bf90aac3f20cd76ee9366`.
The Linux runtime evidence is under
`/home/c829439661/src/smr-jax-verify-06-089c1ee-run1/runtime/`.
It records `passed: true`, `strict_passed: false` and
`qualification: practical-acceptance-passed` for both benchmarks. The CLI exits 0.

Reproduce from that bounded exact-commit snapshot using the same radiation
cache and environment described below:

```bash
python -B -m smr_jax.verify_gpu --radiation-dir runtime/radiation \
  --snapshot-manifest snapshot.json --days 2 30 --repeats 5 \
  --acceptance-profile practical --melt-small-snow-next-day \
  --report-path runtime/verify06.json
```

## Acceptance amendment - 2026-11-13

The operator explicitly rejected strict parity as the release gate: "the strict
parity is too strict". Practical acceptance now checks every active watershed
cell, every simulated day, the actual terminal state and every basin scalar.
The original full-array strict result is retained separately as
`strict_comparison`/`strict_passed`, including extra-domain snow differences.
It does not block practical acceptance. This changes verification policy;
model physics, compiler precision settings and snow-release behavior are unchanged.

The practical profile uses absolute engineering limits, with no relative
allowance, averaging, percentiles or exclusion of difficult active cells:

| Quantity | Absolute CPU/GPU difference limit |
|---|---:|
| Water stores and flux depths | 0.001 cm (0.01 mm) |
| Energy quantities | 0.001 in the model's energy units |
| Temperature / saturation | 0.001 C / percentage points |
| Discharge / baseflow | 1e-10 cms |
| Legacy reservoir carry | 1e-9 in its retained legacy units |
| Albedo | 1e-6 |

The water limit adopts the existing operator-selected 0.001 cm resolution for
insignificant water residues. These are explicit prototype engineering
acceptance choices, not a claimed literature-derived uncertainty bound or
calibration to a newly fitted parameter. They are declared in the verifier
and checked against the complete output maxima. Larger errors still fail.

Masks, types, shapes, counters and complete finite footprints retain their
original exact checks, including outside the active numeric domain. Missing
or nonfinite required values fail. External plant-water comparisons and
independent conservation checks retain their original strict limits. Raw
snow-presence disagreements remain visible. Unknown units fail explicitly.

`--acceptance-profile practical` is the default for watershed runs;
`--acceptance-profile strict` retains the original strict gate. Synthetic
full-model and plant-water cases still have to pass their strict checks.
Timing qualification explicitly names the profile, so practical acceptance
never establishes strict arithmetic identity or default-JIT compatibility.

## Original strict audit and boundary

The original audit below records GPU execution and synthetic cases passing
on the actual RTX device. Strict complete-array watershed parity failed. Its checker
returns a nonzero exit for failed comparisons and does not publish a validated
speedup for them. VERIFY-04's separate Windows/default-JAX legacy deferral is
unchanged. At that audit, the packet remained on its packet branch pending
acceptance. The acceptance amendment above supersedes strict-only blocking;
the historical numbers and rejected comparisons below remain unchanged.

The same `run_days` and external plant-water functions run on CPU and CUDA in
the same existing WSL Python environment. No physics, native FFI, callbacks,
driver, dependency installations or existing Linux checkout were changed.
Windows remains the canonical source; Linux receives a new exact-commit
execution replica with a per-file SHA-256 manifest.

## Compiler precision finding

The initial CUDA compiler profile removed explicit float64-to-float32-to-float64
round trips used to reproduce GRASS FCELL raster boundaries. Float64 output
types alone did not establish that these intermediate roundings survived.
For example, the radiation fixture's rounded albedo differed at 4,608 of
7,680 entries, and shortwave differed by up to `0.0009937133781932062`
kJ/m2/day. Disabling selected HLO passes or using optimization level zero did
not fix this diagnostic.

Setting `xla_allow_excess_precision=False` on both compiled executables preserves
these source boundaries. All seven diagnostic radiation stages then match
bit for bit, including final shortwave. OpenXLA describes the option as
allowing increased instruction output precision in its
[debug-option implementation](https://github.com/openxla/xla/blob/main/xla/debug_options_flags.cc).
The checked-in regression compares compiled casts with explicit NumPy rounding
on every device returned by the test process.

The verification executable also disables CPU fast math and fast min/max on
both backends. These are explicit compiler settings in the report, not a
replacement for any scientific function. The ordinary application's default
GPU JIT configuration has not been qualified by these results.

## Comparison contract

The primary result compares every finite CPU output value, the complete daily
history and the actual terminal state. Shapes and types must match; floating
outputs must be float64; masks, boolean/integer values and counters are exact.
Required inputs define the active watershed independently of GPU output.
Missing or nonfinite required values fail. The full CPU/GPU finite footprints
must match, including outside the active domain; their intersection is never
used to hide missing support.

Absolute budgets remain `1e-12 cm` for water, `1e-9` for energy,
`1e-10` for temperature/percent, `1e-15 cms`, `1e-12` for reservoir carry,
`1e-10 m3` for domain ledger residuals and `1e-12` for albedo, plus
`512 * float64 epsilon` relative error. Unknown fields fail rather than inherit
an implicit tolerance. These budgets were not fitted to device differences.

An additional `required_domain` result applies the same numeric budgets to
every active cell and all scalar outputs. It does **not** replace primary
acceptance or qualify timing. Snow-presence disagreements are also counted
inside the active domain, even when their absolute water difference is tiny.

Both seven-day synthetic full-model runs pass: literal legacy snow behavior
and the existing selected `<=1e-3 cm` next-day release. They include dry/wet,
snow/melt/refreeze, all configured land covers, a mask hole, zero flow and
domain edges. Maximum synthetic water difference is `2.1316282072803006e-14 cm`;
maximum energy difference is `4.3655745685100555e-11`.

All ten external plant-water cases pass full CPU/GPU comparison and independent
layer, cell and domain water-ledger checks on both devices. Five physical
cases cover dry/excess demand, a wet pulse, ordinary demand, no input/demand
and consistent snow transfer; each runs with closed and legacy-drainage bottom
boundaries. These are synthetic coupling verification cases, not validation
of every legacy watershed parameter as a plant-interface input. A separate
read-only real-grid coupling probe rejected negative/nonfinite legacy
`fieldcap_b_cm` values through the existing input validator; no parameter
repair or successful real-grid coupling benchmark is claimed.

## Actual hardware and environment

- Intel Core i7-12800H, 20 logical CPUs visible in WSL.
- NVIDIA GeForce RTX 3070 Ti Laptop GPU, 8 GB physical VRAM, driver 566.24.
- Ubuntu 24.04 under WSL2; Linux `6.6.87.2-microsoft-standard-WSL2`,
  glibc 2.39, Python 3.12.3.
- JAX, jaxlib, CUDA12 plugin and PJRT `0.11.1`; NumPy `2.5.3`.
- Existing CUDA runtime `12.9.79`, nvcc `12.9.86`, cuDNN `9.26.0.51`,
  cuBLAS `12.9.2.10`; the JSON report inventories the remaining packages.
- `JAX_ENABLE_X64=true`, `JAX_PLATFORMS=cuda,cpu`,
  `XLA_PYTHON_CLIENT_PREALLOCATE=false`; GPU allocator limit 6 GiB.
- Existing interpreter: `/home/c829439661/venvs/smr-jax/bin/python`.

CUDA-only foundation tests pass (3 tests). The final harness tests pass on
the CUDA-enabled environment (13 tests), including explicit GPU cast rounding.
Windows CPU regression ran 131 tests successfully, with five optional native
library tests skipped; the subsequent 13 targeted harness tests also pass.

## Watershed comparison results

Both reference runs use 68,121 independently defined active cells and the
existing selected next-day snow-release policy. Complete histories, scalar
basin outputs and actual terminal state are retained and compared.

| Result | Two days | Thirty days |
|---|---:|---:|
| Strict full-array comparison | Fail: 16 fields | Fail: 34 fields |
| Additional active-domain comparison | Pass | Fail: 18 fields |
| Largest active water difference, cm | 7.1054e-15 | 1.7706e-4 |
| Active terminal total-soil-water difference, cm | 7.1054e-15 | 1.5144e-5 |
| Active terminal layer-water difference, cm | 7.1054e-15 | 1.4114e-4 |
| Largest active snow-water difference, cm | 8.8818e-16 | 1.0658e-14 |
| Largest active energy difference | 5.8208e-11 | 1.1642e-10 |
| Active discharge difference, cms | 0 | 7.5184e-14 |
| Raw active snow-presence disagreements, cell-days | 2,370 | 2,370 |

The month water maximum is approximately 1.77 micrometres of water depth;
it remains above the existing strict tolerance. It occurs in a daily layer
storage map, rather than the terminal total store. Month active snow/heat,
albedo and exact counters pass. The soil discrepancy's first causal operation
has not been isolated; no physics correction or new tolerance is inferred.

The large selected-policy snow/temperature failures (`0.08192291949773434 cm`,
`1.9575634040258594 C`) are confined to 58 cells **outside** the active domain.
The selected release deliberately applies only to active cells. Raw daily SWE
still differs at the sign boundary inside the watershed: a first-day trace
shows zero versus `5.551115123125783e-17 cm`. The next-day policy removes this
residue before the next active-cell snow-presence decision. The raw presence
count is reported even though the active snow/heat histories pass numerically.

For comparison, the strict literal-legacy audit at revision `c04bcb6` fails
both two and thirty days, with up to `0.08458120608596698 cm` snow difference
and `2.028215874817433 C` temperature difference. Those results are retained;
selected-policy results do not establish literal-legacy GPU parity.

Closing strict arithmetic parity would require resolving the month soil difference and
the extra-domain full-array snow differences, or an explicit operator decision
to accept a narrower/toleranced boundary. The operator subsequently chose
practical acceptance above; strict diagnostics are never silently marked passed.

## Timing method and limitations

Every device transfer, first execution and repeated execution is synchronized
with `block_until_ready`, following
[JAX benchmarking guidance](https://docs.jax.dev/en/latest/benchmarking.html).
The harness records input loading/initialization, host-to-device transfer,
lowering, compilation, first execution, complete output download and allocator
memory separately. A cold phase sum is labeled as a sum, not an independently
measured wall-clock workflow. Successful cases additionally receive five
warmed executions and a separately measured warm workflow with fresh local
input loading, initialization, transfer and full host output download.

By default, failed full-array comparisons receive no warmed timing. Explicit
`--diagnostic-timing` permits synchronized warmed measurements while marking
them `unqualified-full-array-parity-failed`. It never adds a validated speedup
or changes the failed result/nonzero exit. The final audit uses this option.
The complete thirty-day output contract contains about 1.236 GB of arrays;
downloading that history can dominate GPU runtime. GRASS/r.sun preparation is
outside the prepared-input measurement boundary. No optimized terminal-only
contract, other GPU, power setting or concurrent-load isolation is claimed.

Final selected-policy timing observations, in seconds, are **unqualified**:

| Measurement | 2-day CPU | 2-day GPU | 30-day CPU | 30-day GPU |
|---|---:|---:|---:|---:|
| Host-to-device/input placement | 0.00210 | 0.13571 | 0.00307 | 0.29582 |
| Lowering | 0.19755 | 0.24280 | 0.17620 | 0.17769 |
| Compilation | 1.75334 | 1.72913 | 1.28850 | 1.39232 |
| First synchronized execution | 0.19287 | 0.14974 | 2.43052 | 0.78305 |
| Median of five synchronized warmed executions | 0.12432 | 0.09089 | 2.12029 | 0.48458 |
| First complete output download | 0.00235 | 0.13644 | 0.00240 | 2.34241 |
| Measured warm prepared-input workflow | 1.28632 | 1.01263 | 3.57545 | 3.19540 |
| Cold executable phase sum, including loading | 4.31842 | 4.56403 | 5.39921 | 6.48982 |

Loading/initialization used in those phase sums was `2.17021 s` for two days
and `1.49851 s` for thirty days. The workflow measures include fresh input
loading/initialization and full host output transfer; they exclude compilation,
JAX import and initial device-context setup. The cold executable phase sum also
excludes process startup and device-context setup, which already occurred in
the synthetic checks. Timing samples vary substantially; the JSON retains
every sample rather than implying a controlled steady-state benchmark.

The GPU allocator's cumulative process peak was `270,532,608 bytes` after
the two-day benchmark, then `1,608,254,720 bytes` (about 1.50 GiB) after thirty
days, below its 6 GiB limit. Output sizes were `94,911,355` and `1,236,217,339`
bytes. The faster resident GPU execution is substantially offset by moving
the full history back to the host. Neither timing column establishes an
accepted speedup while full-array parity fails.

## Original strict-audit reproduction and evidence

```powershell
python -B tools/package_jax_gpu_snapshot.py --revision df6586086b2f20313ad2cf36e85807a28bf19245 --output <NEW-directory>
```

Extract the bounded archive into a new Linux directory, copy the manifest and
prepared radiation cache, then run using the unchanged existing environment:

```bash
export PYTHONPATH="$PWD/src" JAX_ENABLE_X64=true JAX_PLATFORMS=cuda,cpu
export XLA_PYTHON_CLIENT_PREALLOCATE=false
python -B -m smr_jax.verify_gpu --radiation-dir runtime/radiation \
  --snapshot-manifest snapshot.json --days 2 30 --repeats 5 \
  --audit-all --melt-small-snow-next-day --diagnostic-timing \
  --report-path runtime/verify06.json
```

Omit `--melt-small-snow-next-day` to audit the literal legacy policy. Without
`--audit-all`, a failed synthetic or two-day comparison stops the larger audit.
A real CUDA device is mandatory; CPU fallback fails. Report paths and snapshot
directories must be new. Neither command installs anything.

Strict-audit executed source: `df6586086b2f20313ad2cf36e85807a28bf19245`.
Source archive SHA-256:
`338b6916aecd3c3ddb6c67d4bec85ae7ee7567db5a4651dcc70453b277534cc7`.
Prepared radiation archive SHA-256:
`e4bf0394008ff425f6c30e082490cf619d8b20605f28cafd415e1ca3e845dce9`.
The provider checks the source bundle, grid, weather identity and component
hash before either device runs.

External evidence is preserved outside Git:

- Windows snapshot/launch scripts:
  `C:\projects\git\worktrees\smr-verify06-snapshot-df65860`.
- Linux snapshot/results:
  `/home/c829439661/src/smr-jax-verify-06-df65860-run1/runtime/`.
- `verify06.json`, `verify06.log`, `foundation.log`, `harness-tests.log`,
  `nvidia-smi.txt` and `radiation/` identify the final execution.
- A review copy of `verify06.json` is in the Windows snapshot directory;
  SHA-256 `3cf2ceb9c12e44a96ed3f2b6af82b0473ecdf3b7bcd39d27396f309d5e0f5f92`.
- Earlier snapshots `08bcb36`, `c04bcb6`, `b630d9c` and `87ce2a7` retain rejected profiles,
  the literal-legacy audit, the selected two-day audit and the snow diagnostic.

The source manifest distinguishes the executed revision from later evidence-only
documentation commits. Existing worktrees, uncommitted work, reference archives,
the plant repository, EOS protected chats and the global board are preserved.
The report-only scope check found all five changed paths within this packet's
declared boundary. No strict scope check or new Git hook gates commits/pushes.
