---
title: Bounded SMR plant water compatibility completion evidence
status: active
owner: SMR_R
last_reviewed: 2026-11-12
---

# PKT-SMR-JAX-PLANT-01

Implementation is complete on the isolated packet branch; integration has not
been performed. This record does not claim plant-side integration or field
validation. The execution definition is the adjacent YAML, identical to the
controller bootstrap at
`C:/projects/git/worktrees/packet-controls/smr-plant-20261112/PKT-SMR-JAX-PLANT-01.yaml`.

## Boundary and provenance

- Owning repository: `C:/projects/git/SMR_R`.
- Feature target: `feature/smr-jax-plant-water`, created at
  `4b5733ec1a9da2a86572975f3b1e1564cdbd8616`; not advanced or merged.
- Implementation branch:
  `packet/PKT-SMR-JAX-PLANT-01/add-bounded-external-plant-water-couplin`.
- Worktree: `C:/projects/git/worktrees/PKT-SMR-JAX-PLANT-01`.
- Base: `4b5733ec1a9da2a86572975f3b1e1564cdbd8616`, matching both committed
  JAX reference tips. Existing dirty VERIFY work was inspected and preserved.
- Scope: exactly `src/smr_jax/plant.py`, `tests/test_jax_plant.py`,
  `examples/plant_water.py`, `docs/smr-plant-interface.md`, and this packet's
  Markdown/YAML. No changes to existing engine, fluxes, tests or fixtures.
- No merge, push, branch deletion, worktree removal, default-branch advance,
  board/roadmap mutation, repository restructure, or plant-code edit.

The obsolete engineering-os scripts now refuse execution and direct callers
to `engineering-os-toolkit`; active data instructions resolve EOS data to
`eos-data`. Live toolkit scripts were used for create, preflight and scope.
The checked-out EOS status is newer than the historical startup snapshot;
this product packet neither changes nor resumes any EOS structural packet.

## Delivered

| Output | Result |
|---|---|
| Compatibility findings and decisions | `docs/smr-plant-interface.md`: actual code vs proposals, units, time order, boundaries and follow-ups. |
| External vegetation-water API | Exclusive external interception/ET; requested, actual and unmet water; proportional competing layer extraction; checked entry refuses unmet demands by default. |
| Water diagnostics | Explicit A/B downward/upward movements, per-layer routing out/in/export, bottom drainage, return flow, evaporation/extraction, snow handoff and substep bounds. |
| Independent verification | NumPy layer/cell/domain budgets and explicit donor/receiver-edge reconstruction; absent or nonfinite required diagnostics fail. |
| Two-layer preservation | Original code and fixtures unchanged; synthetic agreement with literal old soil sequence when bounds are inactive. |
| Arbitrary-layer assessment | Explicitly unsupported by current backend. Documented selectable-backend boundary and unresolved science; no stand-in port. |
| Interface example | `examples/plant_water.py`, synthetic properties and tomorrow's contract kwargs, no dependency on plant files. |
| Synthetic tests | 22 new tests plus all 8 original JAX tests pass. |

## Validation and numerical evidence

Executed on Windows CPU with Python from the existing plant virtual
environment, **without installing into or changing that environment**:

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PYTHONPATH="$env:TEMP\smr-jax-plant-testdeps"
& 'C:\projects\git\plant-game-of-life\.venv\Scripts\python.exe' -B -m unittest discover -s tests -p 'test_jax_*.py' -v
& 'C:\projects\git\plant-game-of-life\.venv\Scripts\python.exe' -B examples\plant_water.py
```

JAX version: `0.12.0.1`, device `CpuDevice(id=0)`, float64 enabled. The shell's
default interpreter lacked JAX. The available JAX interpreter lacked tifffile,
so the initial baseline had 7 passes and one environment error. Installed only
the existing project's tifffile dependency (`2026.11.8`, no dependencies) into
the temporary `smr-jax-plant-testdeps` directory using the default interpreter's
pip. Sandbox access to that temporary dependency required elevated execution
of the final suite. No reference virtual environment or source was modified.
With the dependency accessible, the final run had **30 tests pass in 29.620 seconds**.
Python syntax checks passed for the three new Python files.

Coverage includes zero/dry soil, zero/excessive requests, proportional
realization, same-day wetting, interception bounds, supersaturation return
flow, donor flow bounds, upward redistribution, all eight routing orientations,
open exports, inactive holes, duplicate-donor refusal, closed/draining bottoms,
snow accumulation/melt/vapor, invalid properties, missing/bad diagnostics,
packed-cell mapping, JIT and multi-day scan. Layer-boundary sensitivity uses
10/90, 30/70 and 50/50 cm at conserved total depth and soil stocks; it does not
claim arbitrary-layer convergence.

Additional synthetic runs used `verify_water_balance` to recalculate budgets
with NumPy rather than reuse the kernel's residual values. Uniform synthetic
profile: 30/70 cm, WP=.10, FC=.28, saturation=.45, cell area 900 m2; standard
requests .08/.12 cm transpiration, .04/0 cm soil evaporation and .02 cm
interception; standard rain .4 cm. Residuals below are independent ledger
values, rounded for presentation:

| Scenario | Max absolute layer residual (cm) | Max absolute cell residual (cm) | Absolute domain residual (m3) |
|---|---:|---:|---:|
| Standard example, initial theta .20 | 2.554e-15 | 5.829e-16 | 4.934e-15 |
| Dry, theta 0, no rain | 0 | 0 | 0 |
| Excess demand, theta .11, no rain, closed bottom, 10/5 cm requests per layer | 0 | 0 | 4.996e-16 |
| Wet pulse, theta .40, 40 cm rain, zero ET, closed bottom | 0 | 0 | 0 |
| 3x1 routed grid, theta .20, 10-degree legacy slope, north weights 1 for lower rows | 5.523e-15 | 2.665e-15 | 5.367e-14 |
| Snow handoff: SWE 1 -> .72, melt .3, vapor +.02, rain .2 | 6.550e-15 | 7.959e-15 | 7.138e-14 |

The wet-pulse case exports 35 cm as surface runoff/return flow. The routed
case exports approximately .000078 summed cell-cm through unreceived lateral
fractions. The standard example realizes all .20 cm transpiration, .04 cm soil
evaporation and .02 cm interception; its kernel-reported maximum layer
residual is 2.554e-15 cm and domain residual 4.746e-15 m3. Slight differences
between kernel and NumPy totals reflect summation order and remain far below
the fixed `1e-10 cm` numerical tolerance. No field observations, evaluation
splits or tuned parameters contributed to acceptance.

## EOS evidence

The toolkit created the branch/worktree and required run record. The first
sandbox preflight could not resolve the feature branch's shared metadata;
the elevated read-access rerun passed **with zero warnings before any
implementation writes**. It verified the clean tree and current feature base
and echoed the six allowed paths and forbidden references.

The toolkit warned that its operational run file
`SMR_R/.engineering/runs/PKT-SMR-JAX-PLANT-01.json` is not ignored. That file
is not staged or committed. A proposed repository-local follow-up is adding
`.engineering/runs/` to SMR's ignore policy; this packet does not expand its
scope to change that policy in an existing reference checkout. The EOS
controller's run pointer is operational state, not a central board change.

Final report-only `scope.ps1 PKT-SMR-JAX-PLANT-01`: all six changed paths in
scope, zero violations. No strict check or hook was introduced. The packet
YAML remains `active` pending operator integration; this Markdown records the
delivered implementation without claiming that integration happened.

## Remaining work and scientific uncertainties

- **Plant integration open, separate repository:** replace copied substeps and
  estimated flows; separate candidate growth from commitment and use actual
  uptake. Default unmet-demand refusal prevents a silent supply assumption,
  but this packet cannot finalize growth on the plant's behalf.
- **Arbitrary-layer backend deferred:** vertical constitutive laws, infiltration
  timing/capacity, restricted bottom, cross-profile lateral delivery and
  upward flow require explicit scientific decisions. Campbell drainage and
  the stand-in's .01 bottom multiplier were not adopted.
- **Surface/basin integration open:** explicit snow handoff is supported and
  conserved; automatic legacy energy-snow/radiation bridging and full-reservoir
  budgets were not implemented or claimed. Runoff/drainage are exported at the
  soil/snow domain boundary.
- **Plant process review remains plant-owned:** partial-layer root access,
  profile-scaled water index, zero-water growth floor and reserve mobilization
  vs new assimilation remain visible issues, with no SMR physiology workaround.
- **Legacy scientific uncertainty retained:** slope degrees/100, conductivity
  factors, A/B redistribution and receiver-weight routing remain inherited
  choices. Water conservation is not validation of those equations or nitrate
  transport, and two-day/30-day/GPU hydrologic parity is not newly claimed.
- **Reconciliation after integration:** propose a local manifest/interface-index
  update if those surfaces are adopted, and any authorized central scheduling
  update separately. No global surface is worker-owned here.
