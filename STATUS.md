# STATUS — SMR_R Development State
<!-- AI-optimized: load alongside CONCEPT.md to understand current state without re-reading code.
     Last verified against repo: 2026-06-20. Update this file when work progresses. -->

## Summary
All four pipeline stages have been rewritten and bug-fixed into `R/` (canonical versions). Original scratch scripts preserved in `R/scratch/` for reference. Environment is set up. Stage 1 is in progress — blocked on a FedData v4 / SSURGO temp-path bug (see Known Issues). Stages 2–4 are ready to run once Stage 1 clears.

---

## Component Status

| Component | Canonical file | State | Notes |
|-----------|---------------|-------|-------|
| Spatial preprocessing | `R/SMR_data_setup.R` | 🔄 Running | Stage 1 blocked on SSURGO temp-path bug (see Known Issues #1) |
| Weather builder | `R/wxData.R` | ✅ Ready | Rewritten: GRIDMET-only, FedData v4 / terra compatible, row_crop DOY bug fixed |
| Simulation engine | `R/smr_buildup.pl` | ✅ Ready | Hardcoded Mac paths removed; takes CLI args; annual map init now dynamic |
| Validation / plots | `R/functional_PP.R` | ✅ Ready | `annual_mass_balance` fixed (ET + perc + ΔS added); `radiation_ts` bug fixed |
| Pipeline runner | `R/run_pipeline.R` | ✅ New | Ties all four stages together with consistent paths |
| Config | `scratch/config.yaml` | ✅ New | Bull Run OR (USGS 14138870), GRIDMET WY2015–2019 |
| GRASS workspace | `scratch/grass_workspace/BullRun/PERMANENT` | ✅ Created | UTM Zone 10N (EPSG:32610) |

---

## Known Issues

### 1. ⚠️ SSURGO download fails — FedData v4 temp-path bug (ACTIVE BLOCKER)
**Symptom** (2026-06-20):
```
Error in purrr::map2():
Cannot open "C:\Users\C82943~1\AppData\Local\Temp\...\OR051/spatial";
The file doesn't seem to exist.
```
**Root cause**: FedData v4 generates a temp path using the Windows 8.3 short name (`C82943~1`) but then tries to open it with the full username path — the two don't resolve to the same location under some Windows configurations.
**Workaround to try**: Set `TMPDIR`/`TMP`/`TEMP` to a short path with no spaces or special chars (e.g., `C:\tmp`) before running Stage 1:
```r
Sys.setenv(TMPDIR = "C:/tmp", TMP = "C:/tmp", TEMP = "C:/tmp")
```
Or call `get_ssurgo` with `extraction.dir` if a newer FedData patch restores that param.
**Status**: Unresolved — Stage 1 cannot complete until this clears.

### 2. No observed streamflow data in repo
`functional_PP.R` expects a USGS gauge CSV at `scratch/raw_data/observed/BullRun_discharge.csv`.
Format: columns `date` (MM/DD/YYYY) and `mean_discharge_cfs`. Download from USGS NWIS for gauge 14138870.

### 3. No example run or output files
No simulation CSV is committed. Once Stage 3 runs, commit `scratch/smr_output/MFC_mass_balance_1.csv` as the reference output.

### ~~4. GRASS GIS workspace setup not documented~~ ✅ Fixed
Documented in `SETUP.md`. GRASS 8.5 location created at `scratch/grass_workspace/BullRun/PERMANENT`.

### ~~5. Stage linkage is implicit~~ ✅ Fixed
`config.yaml` centralises all paths. `run_pipeline.R` connects all four stages.

### ~~6. `annual_mass_balance` incomplete~~ ✅ Fixed
ET + percolation + ΔS terms added in `R/functional_PP.R`.

### ~~7. Hardcoded paths in Perl script~~ ✅ Fixed
`smr_buildup.pl` now takes `<weather_file>` and `<output_file>` as CLI arguments.

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

1. **Fix SSURGO temp-path bug** (Known Issue #1) — try `Sys.setenv(TMPDIR="C:/tmp")` workaround and re-run Stage 1.
2. **Verify Stage 1 outputs** — check `scratch/processed_data/imitate_smr_setup/` has all `.asc` files including `n/ne/e/.../nw.asc` (flow direction maps).
3. **Confirm `basin_mask_value`** — re-check the `wbt_basins` plot to ensure ID `1` covers Bull Run watershed.
4. **Run Stage 2** — `source("../R/wxData.R")` from `scratch/`.
5. **Import rasters into GRASS** — commands in `SETUP.md`.
6. **Download observed discharge** — USGS NWIS gauge 14138870, save as `scratch/raw_data/observed/BullRun_discharge.csv`.
7. **Run Stage 3** — Perl via GRASS `--exec` (exact command printed by `run_pipeline.R`).
8. **Run Stage 4** — `source("../R/functional_PP.R")`.

---

## Environment (as of 2026-06-20)

| Tool | Location | Status |
|------|----------|--------|
| R 4.5.1 | `C:\Program Files\R\R-4.5.1\` | ✅ |
| WhiteBox Tools | `C:\Users\c829439661\AppData\Roaming\R\data\R\whitebox\WBT\` | ✅ |
| Strawberry Perl | `C:\Strawberry\perl\bin\perl.exe` | ✅ (added to User PATH) |
| GRASS GIS 8.5 | `C:\OSGeo4W\apps\grass\grass85\` | ✅ (needs OSGeo4W env vars — see SETUP.md) |