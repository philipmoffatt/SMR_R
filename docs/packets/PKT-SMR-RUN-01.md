# PKT-SMR-RUN-01 completion record

## Purpose

Capture the source workflow used for the successful Mica Creek, Idaho,
January 2015 reference simulation so later CPU work has a fixed, reviewable
target. This packet preserves source and configuration only; it does not add
downloaded inputs, GRASS locations, or model output to version control.

## Captured behavior

- `scripts/run_smr.ps1` runs spatial, weather, import, simulation, and
  validation stages independently, with explicit runtime root, configuration,
  UTM CRS, and location parameters.
- The R stages consume `SMR_CONFIG` and `SMR_ROOT_DIR`, allowing the same source
  checkout to use isolated runtime directories.
- The Mica Creek configuration uses HUC12 `170103030403`, UTM 11N, and the
  QA-selected basin mask 37.
- Current external-source compatibility is explicit: HTTPS SSURGO cache
  retrieval, an NLCD COG fallback, and the current GRIDMET variable interface.
- Where SSURGO lacks saturated water content, the historical `1.1 * field
  capacity` fallback is used only for missing values. This is an explicit
  legacy assumption, not a calibration result.

## Reference evidence

The completed reference runtime remains at
`C:\projects\git\worktrees\SMR-FIRST-RUN-01\scratch\mica-creek-jan2015`.
Its simulation manifest records success for 2015-01-01 through 2015-01-30.
It is read-only evidence for CPU regression checks.

## Out of scope

No scientific parameters, equations, or generated runtime data were changed.
CPU reporting and GRASS process reductions begin in `PKT-SMR-CPU-01` and later
packets only after this captured source baseline is reviewed.
