# SMR Modernization Plan
<!-- Built from planning Q&A session 2026-04-24. Update as decisions change. -->

---

## Core Architecture Decisions

**Q: What language stack should the rewrite use?**
A: Python for data ingestion and geospatial APIs; Julia for the simulation engine. Python's `dataretrieval` (USGS NWIS), `rasterio`, and `osmnx` (road networks) are too mature to replace. Julia solves the two-language problem — the same syntax runs fast on CPU and GPU without calling C extensions. The boundary between the two is clean: Python writes GeoTIFF / NetCDF; Julia reads them and runs physics.

**Q: Is GPU a nice-to-have or a hard requirement?**
A: Hard requirement. Use `KernelAbstractions.jl` so the same kernel code runs on CUDA (NVIDIA), ROCm (AMD), and Metal (Apple) without rewriting. The spatial raster loop — thousands of cells × 365 days × years — is the primary target. Start CPU-Julia, add GPU module by module.

**Q: Should the system be locked to one watershed or general-purpose?**
A: General-purpose. Any HUC watershed (2/4/8/12) must work. CAMELS is used only to select the initial 10 benchmark catchments — it is not a hard dependency. The data pipeline takes a HUC code as input and builds everything from scratch.

---

## Phase 0 — Working Repository
*Goal: get the existing R/Perl code running end-to-end with documented setup. This is the reference truth before any rewrite.*

**Q: What is broken or undocumented right now?**
A: Four things (from STATUS.md):
1. GRASS GIS workspace initialization is not documented anywhere.
2. The observed streamflow file needed by `functional_PP.R` is not in the repo.
3. Output paths from Stage 1 must manually match Stage 3 input paths — no config ties them.
4. `annual_mass_balance()` in `functional_PP.R` is explicitly flagged as missing flux components.

**Q: What does "done" look like for Phase 0?**
A: A single `run_pipeline.sh` (or `run_pipeline.R`) that, given a HUC code and date range, executes all four stages in order with consistent paths. A `SETUP.md` documents the GRASS GIS location/mapset init steps. `annual_mass_balance` is fixed. An example output CSV is committed so future runs can diff against it.

**Tasks:**
- [ ] Write `SETUP.md`: GRASS GIS install, `grass --text` location creation, projection, mapset
- [ ] Locate or pull USGS observed flow for MFC; commit to `data/observed/`
- [ ] Audit path constants across all 4 scripts; centralize into one config block or `.env`
- [ ] Fix `annual_mass_balance` — read function, identify missing terms, add them
- [ ] Write `run_pipeline.R` or shell script connecting all 4 stages
- [ ] Commit one example simulation output CSV to `data/example_output/`
- [ ] Confirm end-to-end run produces NSE > 0 on MFC test period (1965–1970)

---

## Phase 1 — Architecture Translation
*Goal: rewrite as a clean, modular, testable Python + Julia codebase.*

**Q: How should the repository be structured?**
A:
```
smr/
  ingest/          # Python — all data fetching and preprocessing
    dem.py         # DEM download, WhiteBox conditioning, culvert burning
    soils.py       # SSURGO download and processing
    landcover.py   # NLCD download and reclassification
    weather.py     # GRIDMET + USGS GHCND + elevation corrections
    usgs.py        # streamflow from dataretrieval
    huc.py         # HUC boundary fetch, any level
  engine/          # Julia — simulation physics
    snow.jl        # SAM: energy balance, albedo decay, SWE
    soil.jl        # two-layer moisture, Darcy lateral flow, percolation
    et.jl          # PET, Kc, actual ET, soil moisture reduction
    radiation.jl   # Bristow-Campbell, Unsworth-Monteith
    baseflow.jl    # linear recession reservoir
    routing.jl     # D8 lateral flow accumulation
    freeze.jl      # frost depth, Ksat impedance, freeze-thaw cycling
    run.jl         # daily timestep loop, calls all modules
  gpu/             # Julia — GPU kernels (added per module)
    kernels.jl     # KernelAbstractions.jl kernel definitions
  benchmark/       # Python — CAMELS selection, validation metrics
    select.py      # choose N watersheds from CAMELS or NWIS
    metrics.py     # NSE, KGE, PBIAS, FHV (peak flow bias), FLV (low flow bias)
    report.py      # generate comparison plots
  tests/
    test_snow.jl
    test_soil.jl
    test_et.jl
    test_freeze.jl
    test_ingest.py
```

**Q: What is the file interchange format between Python and Julia?**
A: GeoTIFF for spatial rasters (readable by both `rasterio` and `Rasters.jl`). NetCDF for time-series weather forcing (standard, CF-conventions). Julia reads these at simulation start; no runtime Python dependency.

**Q: How should modules be tested?**
A: Each physics module gets unit tests with synthetic inputs — e.g., feed `snow.jl` a known energy budget and assert melt volume. Use Julia's built-in `Test` stdlib. Python ingest modules use `pytest`. A CI workflow runs all tests on push. The Phase 0 example output CSV becomes the integration test target: full run must reproduce it within tolerance.

**Key Julia packages:**
- `Rasters.jl` — raster I/O and spatial operations
- `NCDatasets.jl` — NetCDF weather file I/O
- `CUDA.jl` / `KernelAbstractions.jl` — GPU
- `DifferentialEquations.jl` — available for future ODE-based process modules
- `Test.jl` — unit testing

**Key Python packages:**
- `dataretrieval` — USGS NWIS streamflow and site metadata
- `rasterio`, `geopandas`, `xarray` — raster and vector I/O
- `whitebox` — DEM conditioning
- `osmnx` — OpenStreetMap road network
- `climateR` (or direct GRIDMET API) — weather data
- `pytest` — testing

---

## Phase 2 — GPU Acceleration
*Goal: run the per-cell daily loop on GPU, targeting 10–100x speedup on large watersheds.*

**Q: Which computations benefit most from GPU parallelism?**
A: The per-cell daily physics loop is embarrassingly parallel in space — each cell's state update (snow, soil moisture, ET) is independent of neighbors except for lateral flow routing. Priority order:
1. Snow energy balance (`snow.jl`) — 6 sub-steps × all cells, compute-heavy
2. Soil moisture update (`soil.jl`) — matrix-wide θ update
3. ET calculation (`et.jl`) — simple but high volume
4. Lateral flow routing (`routing.jl`) — sequential along flow path; GPU gains are smaller here

**Q: How do we avoid rewriting everything twice (CPU and GPU)?**
A: `KernelAbstractions.jl` — write one kernel, pass `CUDABackend()` or `CPU()` as a parameter. Start all modules on CPU. Add `@kernel` annotations module by module once the CPU version is validated. GPU becomes a runtime flag, not a separate code path.

**Q: What memory layout should rasters use?**
A: Store each parameter as a 2D Julia `Array{Float32}` (rows = northings, cols = eastings). Float32 halves GPU memory vs Float64 with negligible precision loss for hydrology. Stack daily state as a 3D array `[rows, cols, vars]` to allow coalesced GPU memory access. Avoid Julia's `DataFrame` for cell-level state — stay in arrays.

---

## Phase 3 — Benchmark Network (CAMELS / Any HUC)
*Goal: automated validation against 10 reference watersheds, expandable to any HUC.*

**Q: How do we select the initial 10 watersheds?**
A: Pull from CAMELS (531 minimally-disturbed CONUS catchments with pre-computed attributes). Select 10 by:
- Spread across Level I ecoregions (at least 5 different)
- Drainage area 10–500 km² (compatible with 30 m grid feasibility)
- Snow-dominated, rain-dominated, and mixed regimes represented
- USGS record length ≥ 20 years continuous

Script this selection in `benchmark/select.py` using the CAMELS attributes CSV — reproducible, not hand-picked.

**Q: How does the system generalize to any HUC?**
A: `ingest/huc.py` takes a HUC code string, queries the USGS WBD (Watershed Boundary Dataset) via `dataretrieval` or direct REST, returns the boundary polygon. All downstream ingest steps (DEM, soils, NLCD, weather, observed flow) clip to that polygon. A single entry point:
```python
run_watershed(huc="17060208", start="1990-10-01", end="2000-09-30")
```
No hardcoded MFC paths anywhere in the new codebase.

**Q: What validation metrics beyond NSE and KGE?**
A: Add PBIAS (percent bias — catches systematic over/underestimation), FHV (top 2% flow volume bias — peak flow performance), FLV (bottom 30% flow volume bias — baseflow performance). Report all five per watershed. Target: NSE > 0.6, KGE > 0.6, |PBIAS| < 15% as passing thresholds.

---

## Phase 4 — DEM Depressions and Culverts
*Goal: distinguish real storage depressions from routing artifacts; correctly route flow through road culverts and under bridges.*

**Q: How do we identify real depressions vs. DEM artifacts?**
A: After WhiteBox breach (which handles routing artifacts), depressions that remain in the original unbreached DEM and are spatially coincident with flat or concave terrain signatures are candidates for real storage. Algorithm:
1. Compute fill-minus-breach difference raster — cells with large difference are likely real depressions.
2. Apply minimum area and depth thresholds (tunable, e.g., > 0.5 m depth, > 4 cells area).
3. Assign each qualifying depression a maximum surface storage depth (cm) written to the `depression_storage` parameter raster.
4. In the simulation: precipitation fills depression storage before generating runoff; evaporation draws it down.

**Q: How do we detect and burn culverts?**
A: Road-stream intersection approach using OpenStreetMap via `osmnx`:
1. Download road network for the watershed.
2. Extract stream network from DEM (same WhiteBox threshold used elsewhere).
3. Find road-stream crossing points (intersection of road linestrings and stream polylines, buffered 15 m for DEM imprecision).
4. At each crossing, sample DEM elevation profile across the road — a culvert shows as an elevation ridge interrupting the stream valley.
5. **Burn**: lower the DEM by the ridge height at crossing points, leaving a channel at the estimated culvert invert elevation.
6. **Bridge**: if road segment crosses a stream on a span > ~10 m (detectable from OSM `bridge=yes` tag), treat as transparent — do not modify DEM, do not assign storage.
7. Re-run WhiteBox flow conditioning on the burned DEM.

This is implemented in `ingest/dem.py`. The culvert burn is optional (flag `burn_culverts=True`) so the original WhiteBox-only workflow still works.

---

## Phase 5 — Plant Growth Dynamics
*Goal: move from static Kc lookup tables to dynamic vegetation that responds to climate.*

**Stage 5a — Improved Kc curves (implement with Phase 1)**
Replace the hardcoded seasonal Kc values in `wxData.R` with parameterized curves per land cover class and ecoregion. Parameters (green-up DOY, peak Kc, senescence DOY) are stored in a config file, not hardcoded. This is a small lift and improves all benchmarks immediately.

**Stage 5b — LAI-driven canopy model (Phase 5, future)**
Replace Kc curves with a leaf area index (LAI) model:
- LAI driven by GDD (growing degree days) accumulation with land-cover-specific base temperatures.
- Transpiration = PET × f(LAI) instead of PET × Kc.
- Canopy interception storage scales with LAI (replaces fixed 0.31 cm max).
- Stomatal resistance from Ball-Berry or simplified Jarvis model.

**Stage 5c — C/N cycling (long-term roadmap)**
Add litter decomposition, soil organic matter pools, and N mineralization. This feeds back to plant growth rate (N availability → LAI) and soil hydraulic properties (OM → water retention). Scope is significant — plan as a separate module with its own development cycle when 5b is validated.

---

## Phase 6 — Cold Season Hydrology and Freeze-Thaw
*Goal: correct infiltration and runoff generation in frozen or thawing soils; accurate rain-on-snow.*

**Q: How do we model soil frost depth?**
A: Stefan equation — frost depth scales with the square root of cumulative freezing degree days (CFDD). Frost front advances downward daily when T_soil < 0°C. Implement in `freeze.jl`:
- Track `frost_depth` (m) as a state variable per cell.
- When frost depth > 0, effective Ksat of Layer A is reduced by a freezing impedance factor (typically 0.01–0.1× unfrozen Ksat, tunable).
- Full infiltration block when soil water content at freezing is above a threshold (ice-saturated condition).

**Q: How do we handle diurnal freeze-thaw cycling?**
A: Use the existing 6-hour temperature sub-steps already in the model. At each sub-step:
- If T_6hr crosses 0°C, toggle near-surface Ksat between frozen and thawed value.
- Track latent heat of fusion consumed/released — this is already partially handled in SAM; extend it to the soil surface layer.
- Net daily effect: reduced infiltration capacity on mornings when surface refreezes overnight, even if afternoon T > 0°C.

**Q: How is rain-on-snow handled?**
A: Rain falling on an existing snowpack adds energy via:
1. Sensible heat: Q = ρ_water × Cp_water × (T_rain − 0°C) × rain_depth
2. Latent heat if rain refreezes (net heat sink — important for cold snowpacks)
This energy term is added to the SAM energy balance in `snow.jl`. The current SAM model already has the energy accounting structure (`u.total`); rain-on-snow is an additional input term.

**Q: What about porting to a full energy-balance snow model?**
A: The current SAM is a hybrid (energy balance for melt, but simplified albedo and no explicit snow layer stratigraphy). A full model (e.g., SNOWPACK or a simplified version of DHSVM's snow routine) would add snow layer physics, density evolution, and liquid water retention. Plan this as a swap-in replacement for `snow.jl` — same inputs/outputs, better physics. Sequence: implement Phase 6 freeze-thaw first, validate cold-season performance, then evaluate whether the snow model is the remaining bottleneck.

---

## Sequencing Summary

| Phase | Depends on | Deliverable |
|-------|-----------|-------------|
| 0 — Working repo | nothing | R/Perl pipeline runs end-to-end, documented |
| 1 — Translation | Phase 0 as reference | Python ingest + Julia engine, modular, tested |
| 2 — GPU | Phase 1 CPU validated | Same engine, GPU flag, KernelAbstractions kernels |
| 3 — Benchmark | Phase 1 + Python ingest | 10-watershed CAMELS run, any-HUC capable |
| 4 — Culverts/depressions | Phase 1 ingest | `dem.py` with OSM burn + depression storage raster |
| 5a — Kc curves | Phase 1 | Config-driven Kc, replaces hardcoded values |
| 5b — LAI model | Phase 5a validated | LAI state variable, GDD accumulation |
| 6 — Freeze-thaw | Phase 1 | `freeze.jl`, frost depth state, rain-on-snow energy term |
| 5c — C/N cycling | Phase 5b + Phase 6 | Separate module, long-term |

Phases 1, 3, 4, 5a, and 6 can proceed in parallel once Phase 0 is done. Phase 2 (GPU) waits for Phase 1 CPU code to be validated — no point optimizing unvalidated physics.
