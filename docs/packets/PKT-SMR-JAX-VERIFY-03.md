# PKT-SMR-JAX-VERIFY-03 — Lateral routing and basin recurrence

## Scope and result

This packet establishes `python -m smr_jax.verify_routing`, a strict CPU gate
for the literal receiver-weighted lateral rule, GRASS slope defaults and the
reporting-coupled reservoir recurrence. It starts at
`5dbdd71d5a762809c1f86b3288c47fe2cd6c81e2` on
`feature/smr-jax-plant-water`, in its own worktree and declared YAML boundary.

No new routing physics, adaptive stepping or basin unit correction is introduced.
VERIFY-04 still owns whole-engine two-/30-day parity; general radiation and GPU
acceptance remain later gates. The selected next-day snow residue convention is
still separate from unmodified legacy parity.

## Independent GRASS capture

`tools/capture_jax_routing_reference.py` ran against installed **GRASS 8.5.0** in
a **new temporary EPSG:32611 location**, importing the original frozen ASCII
rasters with `r.in.gdal` and applying the source watershed mask. Existing GRASS
runtimes, Perl source and all original frozen artifacts were untouched.

The script checks all seven frozen artifact hashes/sizes before capture and
refuses to overwrite an existing output. It runs the actual Perl map-algebra
assignments, not a NumPy reimplementation, to capture:

- Default slope and the root-depth, wilting-storage and ET-reduction rasters.
- Eight lateral fields (two conductivities, two layer outflows, total outflow,
  and resulting total/A/B storage) from each of the two frozen daily starting
  states, sampled at **386 deterministic cells across all seven observed land
  uses**. These are isolated extra lateral steps from existing daily states,
  not a claim to have recovered missing original substep snapshots.

The compressed fixture is **422,240 bytes**. Its adjacent JSON records the
GRASS version/build, source/capture-script/input/binary SHA-256 hashes, raster
types, reservoir initialization text, support size and sample count. The gate
verifies that provenance, disallows pickle loading, and requires exact equality
between captured support, both frozen daily footprints and the model's entire
**68,121-cell** support. No finite-intersection fallback is used.

Normal tests consume this fixture without requiring GRASS. To reproduce it, run
the capture script inside a new temporary GRASS project and supply a NEW output
stem; the existing fixture is intentionally not overwritten automatically.

```text
grass --tmp-project EPSG:32611 --exec python -B tools/capture_jax_routing_reference.py --output <new-output-stem>
```

## Defects demonstrated and corrected

1. **Slope precision and NULL centers.** The default slope is in degrees and
   stored as FCELL. The old emulation differed on every required cell, by up to
   `1.9010091776294757e-6` degrees. Rounding through float32 before retaining
   float64 engine arrays now matches captured GRASS values exactly. A synthetic
   NULL-center test also exposed that the old validity predicate checked only
   the eight neighbors; it now requires a finite center too. Frozen active
   support remains unchanged.
2. **Masked donors.** Previously, `legacy_receive_gather` masked receivers but
   could read finite backing-array values from inactive donors. A new fixture
   failed with a spurious `0.4666666666666667 cm` reception. Donors are now masked
   before reading neighbors, matching GRASS MASK behavior.
3. **Inherited raster arithmetic.** Fresh GRASS evidence showed up to
   `6.854534149169922e-7 cm` error in the prior wilting-storage calculation.
   The layer products and sum must each round as FCELL before the conditional
   promotes the final raster to DCELL. Corrected root, wilting and ET-reduction
   properties now match GRASS exactly on every required cell.
4. **Conductivity and soil intermediates.** The old double arithmetic differed
   from captured effective conductivity by up to `1.024577942481536e-5`, outgoing
   flow by `2.8061147051649084e-8 cm`, and gathered storage by
   `2.178474023395438e-8 cm`. Raster-only conductivity differences/ratios and
   saturation products now retain the source FCELL rounding. The corresponding
   layer caps and `Ksubsurface * 6` percolation intermediate were corrected too.
   The reusable source oracle now models integer/FCELL/DCELL promotion, and
   VERIFY-02 was rerun without changing its tolerances. Its original evidence
   has an additive clarification rather than a rewritten history.

The default units and type are documented in the
[GRASS slope manual](https://grass.osgeo.org/grass85/manuals/r.slope.aspect.html).
Raster expression promotion is documented in
[GRASS map algebra](https://grass.osgeo.org/grass85/manuals/r.mapcalc.html#data-types-and-their-precision).
The capture is the numerical evidence for the corrections above.

## Lateral acceptance

The gate reads the ten actual Perl lateral assignments and their explicit
neighbor offsets into the bounded source interpreter. It compares all eight
fields on every required cell for both daily starting states, and independently
compares the 386 captured GRASS cells on each day. Synthetic checks include
below/at/above FC and saturation, conductivity interpolation, storage-limited
outflow, zero flowunits, unequal receiver weights and an inactive finite donor.

Another **64 gather cases** cover all eight orientations on `1x1`, `1x4`, `4x1`
and `4x5` grids, with and without mask holes. They use independent explicit
receiver indexing and pass exactly. Existing NaN-neighbor fixtures remain.

A separate domain ledger transposes receiver weights to count each donor's
accepted fraction, then compares the predicted gain/loss to the actual storage
change. It does not assume donor fractions sum to one. This matters: the
synthetic receiver rule can duplicate water, and zero weights can discard it.
Those are retained legacy semantics, not a conservation guarantee for arbitrary
weights. The bounded external plant path retains its existing stricter checks.

The strict CPU profile is unchanged: float64 arrays, fast math off, and
`algsimp`/`fusion` disabled at the explicit compilation boundary. Water/flow
comparisons use `atol=1e-12 cm`, `rtol=512*eps(float64)`. The aggregate domain
ledger uses `atol=1e-10 cell-cm` for summation across the raster; this is separate
from per-cell tolerances. No tolerance was increased to admit a mismatch.

## Basin acceptance and limits

The loader's IDs, area counts, initial reservoir and coefficients are checked
against independently captured initialization text: `1 75218 0 0.05`.
Initial baseflow must be zero. The actual JAX recurrence then consumes all
**30 original full-report aggregate rows**, and its discharge and baseflow are
compared with the printed frozen outputs. Reservoir volume is independently
inferred from reported baseflow and the known nonzero coefficient.

Hand-calculated tests verify subtraction of PREVIOUS baseflow before calculating
today's baseflow, zero forcing, zero coefficients, masks and separate basin
arithmetic. They verify that runoff includes cell area while percolation does
not, and that configured area counts do not replace source cell sums. Reported
aggregates are supplied as aggregate inputs, preserving the legacy reporting
precision; this does not establish the model's daily flux aggregation parity.

The original recurrence's units and update order are preserved, not repaired.
Captured acceptance covers the single watershed with reporting enabled. The
source final-only path skips this reporting-coupled update; JAX does not expose
that reporting mode. Synthetic multiple-ID arithmetic is tested but does not
claim parity with the source's multi-watershed reporting path, which also
contains the unconditional watershed-to-1 conversion and zonal-reporting quirks.

Discharge/baseflow comparisons use `atol=1e-15` in inherited cms units;
reservoir comparisons use `atol=1e-12` in inherited reservoir units, both with
the existing float64 relative budget.

## Validation

Packet preflight passed from a clean, current feature base. The preceding feature
had 70 passing JAX tests. The new checks first demonstrated the slope, masked
donor, raster arithmetic and NULL-center failures before their corrections.

- **All 84 JAX tests passed in 120.156 seconds**, including strict VERIFY-01,
  VERIFY-02, selected residue policy and the bounded external plant interface.
- After adding the final exact-property/configuration assertions, all **16
  routing tests passed in 14.997 seconds**. Failure tests deliberately corrupt
  capture/source provenance, slope, reservoir coefficients, masks, gather and
  discharge. The CLI returns nonzero for verification failures and missing files.
- Final VERIFY-03 CLI passed. The 64 gather cases agree exactly. The two full
  source comparisons have maximum errors `3.552713678800501e-15` and
  `5.684341886080802e-14`; captured GRASS comparisons have maximum errors
  `8.881784197001252e-16` and `5.684341886080802e-14`.
- The independent donor-accounting domain ledger has a maximum residual of
  `2.2026824808563106e-13 cell-cm` over the frozen grids. Synthetic receiver gain
  is retained and accounted for, rather than silently normalized.
- Across 30 frozen reports, maximum discharge/baseflow/reservoir differences
  are `1.214306433183765e-17`, `5.692061405548898e-19` and
  `1.1275702593849246e-17`, respectively.
- Final VERIFY-02 CLI passed its 918 synthetic and 588 representative cases
  with maximum source/total-ledger error `1.4210854715202004e-14 cm`.
- `git diff --check` and the report-only scope check passed. All 12 changed
  paths are inside the packet boundary. Perl and original frozen-artifact diffs
  remain empty. No shared/default branch or existing GRASS runtime was changed.

Runtime: existing read-only `plant-game-of-life/.venv/Scripts/python.exe`,
JAX 0.12.0.1 CPU with the existing temporary TIFF dependency directory. GRASS
was required only for the independently recorded capture, not normal tests.

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PYTHONPATH="$PWD\src;$env:TEMP\smr-jax-plant-testdeps"
& C:\projects\git\plant-game-of-life\.venv\Scripts\python.exe -B -m smr_jax.verify_routing
& C:\projects\git\plant-game-of-life\.venv\Scripts\python.exe -B -m smr_jax.verify_soil
& C:\projects\git\plant-game-of-life\.venv\Scripts\python.exe -B -m unittest discover -s tests -p 'test_jax_*.py' -q
```

Completion closes VERIFY-03 for the frozen single-watershed CPU profile. The
next gate is VERIFY-04 integrated parity, starting with the two captured days;
the remaining January days still require radiation fixtures or a verified
general radiation provider.
