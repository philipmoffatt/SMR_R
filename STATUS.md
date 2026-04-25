# STATUS — SMR_R Development State
<!-- AI-optimized: load alongside CONCEPT.md to understand current state without re-reading code.
     Last verified against repo: 2026-04-24. Update this file when work progresses. -->

## Summary
All four pipeline stages exist and are code-complete. The project is **functional but unvalidated** — no confirmed end-to-end run is documented, one known gap in post-processing, and several integration prerequisites are undocumented.

---

## Component Status

| Component | File | State | Notes |
|-----------|------|-------|-------|
| Spatial preprocessing | `R/scratch/SMR_data_setup.R` | ✅ Complete | Last touched Feb 2023; unit conversions and file renaming added at end of script |
| Weather builder | `R/scratch/wxData.R` | ✅ Complete | Added Mar 2023 |
| Simulation engine | `R/scratch/smr_buildup.pl` | ✅ Complete | 884 lines, added Mar 2023 |
| Validation / plots | `R/scratch/functional_PP.R` | ⚠️ Mostly complete | 230 lines; `annual_mass_balance` function explicitly flagged as incomplete (see below) |

---

## Known Issues

### 1. `annual_mass_balance` is incomplete — `functional_PP.R`
**Source**: DuncanJurayj commit message, 2023-03-30:
> "annual_mass_balance (which is very likely missing components but is simple to improve)"

**Impact**: Annual water balance plots may undercount flux components. Other plot functions (NSE, KGE, time series, SAM diagnostics, radiation budget) are not affected.

### 2. No observed streamflow data in repo
`functional_PP.R` requires a USGS gauge observation file to run comparisons. No such file is committed. The path/format expected by the script must be inferred from the code before validation can run.

### 3. No example run or output files
No `processed_data/`, `imitate_smr_setup/`, or simulation CSV output is committed. A new user has no reference run to diff against.

### 4. GRASS GIS workspace setup not documented
`smr_buildup.pl` calls `r.mapcalc`, `r.sun`, and `r.slope.aspect` inside a GRASS session. The required GRASS location/mapset initialization steps are not written down anywhere in the repo.

### 5. Stage linkage is implicit
No script or config file ties the four stages together. Output paths from Stage 1 must match the input paths hardcoded in Stage 3 (Perl). This linkage lives only in the developer's head.

---

## Contributor History

| Contributor | Role | Last active |
|-------------|------|-------------|
| DuncanJurayj (`duncan_jurayj@brown.edu`) | Primary developer — all R scripts and Perl engine | 2023-03-30 |
| ames-fowler | Merge / sync (Dropbox ↔ GitHub) | 2023-03-30 |
| philipmoffatt | Repo owner (GitHub org) | — |

---

## Git Timeline

| Date | Commit | What changed |
|------|--------|-------------|
| 2023-02-10 | `74f5e83` | `SMR_data_setup.R`: watershed `.ini` folder creation + unit conversion stubs |
| 2023-02-14 | `f5a31ff` | `SMR_data_setup.R`: unit conversions, file renaming, map calculations (+172 lines) |
| 2023-02-19 | `ac5e5d3` | `SMR_data_setup.R`: flow direction, flipped to match Perl output convention |
| 2023-02-22 | `8d316a8` | Merge: ames-fowler syncs DuncanJurayj's setup work |
| 2023-03-30 | `1f4633c` | `functional_PP.R` added by DuncanJurayj (230 lines) |
| 2023-03-30 | `cf2e3e7` | `functional_PP.R` deleted (same day — likely replaced by correct version) |
| 2023-03-30 | `cb9c6cf` | `wxData.R` + `smr_buildup.pl` added |
| 2023-03-30 | `576b366` | Merge: ames-fowler brings in all three new files |

**Last commit: 2023-03-30. No activity since.**

---

## What a New Session Should Do First

Before starting any new work, verify:

1. **Can Stage 1 run?** — Check that `FedData`, `whitebox`, `terra`, `sf` are installed and WhiteBox Tools binary is accessible via `whitebox::wbt_init()`.
2. **Can Stage 2 run?** — Check `climateR`, `rnoaa`, `chillR` install cleanly; NOAA API key may be required.
3. **Is the observed streamflow file available?** — `functional_PP.R` needs it; locate or obtain a USGS gauge file for Mica Fork Creek.
4. **Is GRASS GIS installed and on PATH?** — Required by `smr_buildup.pl`.
5. **Do output paths match?** — Confirm that paths hardcoded in `SMR_data_setup.R` output section match what `smr_buildup.pl` expects as input.

---

## Suggested Next Work (priority order)

1. Document the GRASS GIS workspace setup (location, projection, mapset) needed to run `smr_buildup.pl`.
2. Locate / commit the observed streamflow file so `functional_PP.R` can run.
3. Fix `annual_mass_balance` in `functional_PP.R` — confirm which flux terms are missing.
4. Add a `run_pipeline.R` or shell script that calls all four stages in order with consistent paths.
5. Commit a small example output CSV so output format can be verified without a full run.
