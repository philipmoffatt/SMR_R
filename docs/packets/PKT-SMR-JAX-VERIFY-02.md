# PKT-SMR-JAX-VERIFY-02 — Legacy soil vertical fluxes

## Result and boundary

The fixed six-hour allocation, redistribution, ET, percolation,
saturation-excess runoff and road-runoff calculations have an independent,
strict float64 CPU verification gate: `python -m smr_jax.verify_soil`.

The packet starts at `62b36a8b6258c2e96652730099a5668d74a7ceae` on
`feature/smr-jax-plant-water`. Implementation is isolated in the
`PKT-SMR-JAX-VERIFY-02` worktree. Its YAML records the allowed write paths.

This gate retains four six-hour steps. Adaptive stepping, arbitrary soil layers,
matric-potential physics, lateral routing, reservoir recurrence, integrated
two-/30-day acceptance, general radiation and GPU acceptance remain separate
work. In particular, these results do not establish infiltration-front accuracy.

## Independent evidence

`verify_soil.SourceOracle` extracts the original assignments directly from
`R/smr_buildup.pl`, including the land-cover setup and five soil/surface blocks.
An allowlisted arithmetic/conditional interpreter evaluates those expressions on
the CPU without calling the JAX flux implementations. Scalar conditionals are
lazy, as required for zero-root branches. Hand-calculated ET, road partitioning
and drainage anchors check the interpreter itself. Missing assignments,
unsupported expressions and changed six-hour/30-metre constants fail explicitly.
This is a source-expression oracle, not a fresh GRASS execution or a new frozen
substep capture.

Before reading frozen evidence, the gate checks every artifact size and SHA-256
in `reference/mica-creek-jan2015/manifest.json`. Both daily saturation footprints
must exactly match the model's entire required support of **68,121 cells**.
Root depth, wilting storage and ET-reduction moisture are independently checked
on every required cell.

The numerical comparisons use the strict CPU compilation profile established by
VERIFY-01A: float64, CPU fast math disabled, and `algsimp`/`fusion` disabled at the
explicit host compilation boundary. Water comparisons use `atol=1e-12 cm`;
saturation uses `atol=1e-10` percentage points; relative tolerance is
`512 * eps(float64)`. These budgets were declared before evaluating the cases.
Missing, wrong-shaped, non-float64 and nonfinite required values fail; no finite
intersection or report-only success path exists.

The case matrix contains:

- **918 synthetic cases** spanning all six defined land uses plus default
  codes `0`, `41` and `43`; both sides and equality at wilting/stress/FC/layer
  saturation/profile saturation boundaries; dry, wet and large inputs; zero
  and strong PET; zero and large drainage capacity; road thresholds and both
  water sources; and cold/asymmetric melt-allocation weights.
- **98 actual frozen cells**, selected deterministically across every observed
  class (`1, 2, 3, 4, 5, 41, 43`). January 1 and 2 daily states and weather,
  together with controlled dry/wet/ET probes, produce **588 representative
  cases**. The daily snapshots are starting conditions for isolated column
  checks; they are not presented as captured vertical-substep truth.

Every case compares 23 values: three land-cover properties, road source
partitioning, soil input, the total/A/B storage at each stage, ET, percolation,
runoff and saturation. Independent road, allocation, ET, drainage and complete
column ledgers check total-water conservation.

## Defects found and corrected

The original 57-test baseline passed before the new gate existed. The new gate
then failed on the unmodified implementation: **61,171 of 68,121** wilting-storage
values differed beyond the fixed tolerance, with a maximum difference of
**41.75 cm**. The corrections follow the source expressions rather than changing
the scientific model:

1. Perl overwrites imported `wiltpt_amt` during land-cover setup. The loader now
   reproduces that calculation from layer moisture and depth: grass uses A;
   forest/shrub/crop use A+B; water, urban and unknown classes use zero. The old
   imported raster included both land-cover differences and precomputed rounding.
2. A zero root depth does not force ET to zero in Perl. Positive root storage
   selects its high-storage arm. JAX now preserves that branch while protecting
   unselected divisions. The earlier safety test incorrectly required zero ET;
   it now checks the literal high-storage flux and finite resulting storage.
3. The frozen bundle contains 171 cells with landuse `41` and 111 with `43`.
   Perl's default ET coefficient is `cc_water`, and its default A redistribution
   is zero. JAX previously clipped the coefficient index to crop and applied the
   ordinary nonforest A redistribution. Both defaults now match the source.

No frozen artifact, Perl source, plant coupling implementation or lateral/basin
algorithm was modified. Allocation, percolation, saturation runoff and road
runoff needed no formula correction.

## Water representation retained

The input `sat_amt` raster and independently computed A+B saturation capacities
carry slightly different rounding. Perl clips total storage and layer storage
against those separate quantities. Saturated representative probes therefore
retain a maximum **7.450580596923828e-9 cm** difference between layer sum and
total storage. The gate reports that representation gap and checks it against
the source explicitly; it does not relax the zero-residual total-water ledgers.
Synthetic properties have a gap of at most `3.552713678800501e-15 cm`.

Source parity also preserves signed inputs and inherited algebraic behavior. It
does not certify nonnegative physical states for every frozen snapshot or add a
new bounded-transfer model. The separately bounded external plant path remains
subject to its existing checks.

## Validation

- Packet preflight passed from the clean, current feature base.
- The initial 57 JAX tests passed in 40.839 seconds.
- All 14 focused soil tests passed in 18.530 seconds, including deliberately
  wrong fluxes, incorrect static properties, lost support, nonfinite values,
  missing source assignments, changed timestep and nonzero CLI failure.
- All **70 JAX tests passed in 92.808 seconds**, including the strict two-day
  VERIFY-01 surface/snow gate, the selected residue policy, and the bounded
  external plant interface.
- Final soil CLI passed: 68,121 required cells, 918 synthetic cases and 588
  representative cases. Maximum source-comparison/total-ledger error was
  **1.4210854715202004e-14** in both case groups.
- `git diff --check` and the report-only EOS scope check passed; all seven changed
  paths are inside the declared boundary. Perl and frozen-reference diffs are empty.

Runtime: the existing read-only
`plant-game-of-life/.venv/Scripts/python.exe` (JAX 0.12.0.1 CPU), with the existing
temporary `tifffile` dependency directory on `PYTHONPATH`; no environment changes.

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PYTHONPATH="$PWD\src;$env:TEMP\smr-jax-plant-testdeps"
& C:\projects\git\plant-game-of-life\.venv\Scripts\python.exe -B -m smr_jax.verify_soil
& C:\projects\git\plant-game-of-life\.venv\Scripts\python.exe -B -m unittest discover -s tests -p 'test_jax_*.py' -q
```

## Follow-up

VERIFY-03 owns lateral gather/orientation, the provisional slope convention and
basin recurrence. VERIFY-04 owns whole-engine parity. A later adaptive-step
packet must first replace embedded six-hour rate factors with explicit elapsed
time in a separate selectable path, then verify conservative exchange and
fixed-step refinement under identical physical profiles and forcing. Shared
steps for interacting layers and synchronization at forcing/coupling events are
the proposed starting design, not implemented behavior in this packet.

## VERIFY-03 follow-up: source raster precision

Fresh GRASS 8.5 evidence captured by `PKT-SMR-JAX-VERIFY-03` exposed a limitation
of this packet's original source-only oracle: it evaluated imported raster
subexpressions uniformly in double precision. GRASS retains FCELL arithmetic
inside an expression until a DCELL operand promotes the result. The original
results above remain historical evidence, not proof of GRASS raster precision.

VERIFY-03 corrects the wilting-storage products/sum, layer saturation products,
and the percolation `Ksubsurface * 6 / 24.0` intermediate accordingly. The source
oracle now models FCELL/integer/DCELL promotion, and synthetic imported
properties are represented as FCELL values. Grass stress and bottom-layer-full
threshold probes are placed at those represented thresholds. The updated gate
retains the same tolerances and passes all 918 synthetic and 588 representative
cases. Fresh captured GRASS root depth, wilting storage and ET-reduction moisture
match exactly on all 68,121 required cells. See
[VERIFY-03](PKT-SMR-JAX-VERIFY-03.md) for capture provenance and validation.
