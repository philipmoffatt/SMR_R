---
title: CPU feature integration evidence
status: active
owner: SMR_R
last_reviewed: 2026-11-12
---

# CPU reference feature landing

The operator authorized consolidating the completed CPU work onto a CPU
feature and opening a PR to the repository's default branch. GitHub's repository
metadata confirms the default is `master`; there is no `main` branch. This
does not authorize merging the feature into that default branch.

Feature: `feature/smr-cpu-reference`. Source boundary: `6af1a5e`, the frozen
two-day verification trace capture. This includes the existing reference-run
and CPU optimization history through reporting modes, the single-watershed
reporting fast path, temporary-map lifetime, CPU configuration/timing, and
frozen Mica Creek inputs, daily records, final states and trace maps.

Original commits remain unchanged and reachable. The subsequent JAX foundation
(`23bee09`), GPU setup note (`90c546f`), modular engine (`4b5733e`) and external
plant interface (`cd42e36`) belong to the JAX feature and are excluded from
this CPU PR. The packet name CPU-06 spans reference capture and JAX preparation;
the implementation boundary here follows the actual files, not its name alone.

No implementation files were changed by this integration packet. Its only
new files are this evidence and the adjacent packet YAML. The original
CPU-06 design document describes both the completed CPU evidence and future
JAX work; it remains historical source context.

## Validation on 2026-11-12

All five existing CPU checks passed in this isolated integration worktree:

- `tests/test_benchmark_smr.ps1`: date completeness and comparison fixtures.
- `tests/test_output_modes.ps1`: output-mode contracts and Perl syntax.
- `tests/test_single_watershed_reporting.ps1`: reporting dispatch and Perl syntax.
- `tests/test_map_lifetime.ps1`: map lifetime contracts and Perl syntax.
- `tests/test_cpu_configuration.ps1`: CPU/timing contracts and Perl syntax.

All seven artifact SHA256 values in
`reference/mica-creek-jan2015/manifest.json` were independently recalculated
and match. No fixture was regenerated. R/Perl host-locale fallback warnings do
not affect these checks. This pass did not rerun live downloads or a 30-day
GRASS simulation. The previously recorded 30-day comparison at 68,121 common
valid cells remains historical evidence, not a new run or field validation.

Live EOS create/preflight passed with no preflight warnings. Report-only scope
covers exactly the two evidence/definition files. No new Git hooks or commit
checks were installed. No default branch was advanced, no user files were
discarded, and no worktrees/branches were removed.

## Related preservation and JAX landing

FIRST-RUN changes are preserved in focused commits on its existing branch.
All eight source/config files exactly match successor commit `f626fbb`, already
in this CPU history, so the historical snapshot branch is not merged over
the optimized CPU code. Its acceleration plan remains a historical proposal.

VERIFY's trace reader/comparison tools are preserved on its own branch.
Although its nine tests pass, the reporter finds day-two SWE and energy
discrepancies, so that packet's numerical gate is still open and it is not
part of the finished JAX landing. The completed JAX engine and bounded plant
interface land separately into `feature/smr-jax-plant-water`.

## Publication and reconciliation boundary

The CPU PR should describe the R/Perl/GRASS reference and CPU tooling, target
`master`, and retain original commit identities. Do not claim a JAX engine,
new scientific equations, productivity tuning or field validation. Feature-to-
default merge remains operator-owned. Public repository metadata is readable,
but initial publication probes found no GitHub CLI login and no usable SSH
key; publication requires restoring repository authentication.

The live EOS run record records the actual feature landing commit. These local
completion documents do not update central scheduling or claim default-branch
integration. After the PR is merged, propose any necessary global board/status
reconciliation as a separately bounded follow-up; no shared EOS board is owned
by this product packet.
