# CONCEPT — SMR_R Architecture Reference
<!-- AI-optimized: load this file to understand the system without reading source code. -->

## Identity
- **Project**: Soil Moisture Routing (SMR) + Snow Accumulation and Melt (SAM) watershed model
- **Language stack**: R (preprocessing) + Perl (simulation engine) + GRASS GIS (raster algebra)
- **Study area**: Mica Fork Creek (MFC), Moscow Mountain / Pullman, WA–ID
- **Purpose**: Simulate daily watershed hydrology (snow, soil moisture, ET, streamflow) from spatial inputs

---

## File Map

| File | Role | Reads | Writes |
|------|------|-------|--------|
| `R/scratch/SMR_data_setup.R` | Stage 1: spatial preprocessing | NED, SSURGO, NLCD (remote APIs) | `processed_data/ASC/MFC/*.asc`, `imitate_smr_setup/` |
| `R/scratch/wxData.R` | Stage 2: weather file builder | GRIDMET or NOAA GHCND CSV | space-delimited weather text file (30 cols) |
| `R/scratch/smr_buildup.pl` | Stage 3: simulation engine | `.asc` grids + weather file | daily CSV mass balance per cell |
| `R/scratch/functional_PP.R` | Stage 4: validation & plots | simulation CSV + USGS observed flow | PNG diagnostic plots, NSE/KGE metrics |

---

## Pipeline Sequence

```
[APIs] NED + SSURGO + NLCD
        ↓  SMR_data_setup.R
ASCII parameter rasters (30 m)
        ↓  wxData.R
Weather forcing file (30 cols, daily)
        ↓  smr_buildup.pl (GRASS GIS)
Daily mass balance CSV
        ↓  functional_PP.R
NSE / KGE metrics + diagnostic PNGs
```

---

## Spatial Data Processing (SMR_data_setup.R)

### Inputs (downloaded via FedData)
- DEM: USGS NED 1 m → reproject + resample to UTM 30 m
- Soils: USDA SSURGO (tabular + spatial)
- Land cover: NLCD 2019

### WhiteBox Tools DEM workflow (in order)
1. `wbt_feature_preserving_smoothing` — smooth DEM
2. `wbt_breach_depressions` — remove sinks
3. `wbt_d8_pointer` — flow direction
4. `wbt_d8_flow_accumulation` — accumulation grid
5. `wbt_extract_streams` (threshold: 1000 cells ≈ 0.09 km²) — stream network
6. `wbt_watershed` — basin delineation
7. `wbt_isobasins` (threshold: 1000 cells) — sub-watershed isolation

### Soil processing
- Select majority component per MUKEY (component.comppct_r)
- Two depth layers:
  - **Layer A**: 0–30 cm (surface)
  - **Layer B**: 30 cm – restrictive horizon depth
- Per-layer depth-weighted averages of: organic matter (OM), bulk density (BD), Ksat (cm/hr), θ_sat, θ_fc, θ_wp
- Rasterized via spatial join to SSURGO polygon map

### Land cover reclassification (NLCD → 6 classes)
| Class ID | Label | NLCD codes |
|----------|-------|------------|
| 1 | Water | 11 |
| 2 | Urban / Barren | 21–24, 31 |
| 3 | Forest | 41–43 |
| 4 | Shrub | 51–52 |
| 5 | Grass / Pasture | 71–74, 81–82 |
| 6 | Row crop (wheat) | 82 (subset) |

### Output format
- ASCII grid (`.asc`), GRASS GIS compatible
- CRS: UTM, 30 m cells
- One file per parameter: `elev`, `aspect`, `ksat_a`, `ksat_b`, `bd_a`, `bd_b`, `om_a`, `om_b`, `theta_sat_a`, `theta_sat_b`, `theta_fc_a`, `theta_fc_b`, `theta_wp_a`, `theta_wp_b`, `lc` (land cover), `depth_b` (restrictive horizon), `flow_dir`, `streams`, `basins`

---

## Weather File Format (wxData.R)

### 30-column schema (space-separated, one row per day)
| Col # | Variable | Units | Notes |
|-------|----------|-------|-------|
| 1 | year | — | |
| 2 | doy | — | day of year |
| 3–6 | T6hr_1..4 | °C | 6-hourly temps (0, 6, 12, 18 h) |
| 7 | precip | cm | total daily |
| 8 | rain | cm | precip × rain fraction |
| 9 | snow | cm | precip × snow fraction |
| 10 | PET | mm | potential ET |
| 11 | dew_pt | °C | |
| 12 | cloud_frac | 0–1 | obs/potential solar radiation |
| 13 | Kc | — | crop coefficient by land cover |
| 14 | resist | s/m | heat transfer resistance |
| 15–30 | (additional flux/radiation vars) | mixed | site-specific |

### Data sources
| Period | Source | Station / Dataset |
|--------|--------|-------------------|
| 2000–2022 | GRIDMET (U. Idaho) | gridded 4 km |
| 1960–1979 | NOAA GHCND | USC00456789 (Pullman WA) |

### Derived variables (historical period)
- **6-hourly T**: `chillR::make_hourly_temps()` from daily Tmin / Tmax
- **PET**: Hamon method — f(saturation vapor pressure, day length)
- **Dew point**: Gentilli 1955 threshold formula
- **Cloud fraction**: T_obs / T_potential (ratio of measured to clear-sky solar)

### Elevation corrections (applied to Pullman station data)
- Temperature: −3.20 °C / km
- Precipitation: +0.09195 cm / m
- PET: linear model — `PET_site = 1.8281 × PET_pullman − 0.2486`

### Crop coefficient (Kc) by land cover
| Class | Kc behavior |
|-------|-------------|
| Water | Constant 0.65 |
| Forest | Constant 1.0 |
| Shrub | Seasonal curve |
| Grass | Seasonal curve, active May 1 – Sep 19 |
| Row crop (wheat) | Spring planting through harvest |

---

## Simulation Engine (smr_buildup.pl)

### Time stepping
- Main step: 24 hours
- Snow sub-step: 6 hours (4 per day, requires T6hr cols)

### Process modules (executed each cell each day)

#### 1. Precipitation partitioning
- Threshold temperature splits precip into rain vs. snow
- Canopy interception: land-cover max storage (0–0.31 cm); overflow becomes throughfall

#### 2. Snow Accumulation and Melt (SAM) — energy balance
| Component | Model |
|-----------|-------|
| Shortwave radiation | Bristow–Campbell transmissivity (Campbell 1977; Bristow & Campbell 1985) |
| Longwave radiation | Unsworth–Monteith cloud emissivity (Unsworth & Monteith 1975) |
| Latent heat | Vapor pressure gradient → condensation / evaporation |
| Sensible heat | Convective: f(T_air − T_snow, wind) |
| Albedo | Decays with snow age: 0.7383 → 0.2 |
| Phase change tracking | Internal energy `u.total`; melt when u.total > 0 |

#### 3. Soil moisture — two-layer
- Layer A (0–30 cm) fills first; excess percolates to Layer B
- Lateral flow: Darcy's law with Bresler unsaturated hydraulic conductivity
- At saturation: macropore flow replaces matrix flow
- Saturation-excess runoff generated when θ > θ_sat

#### 4. Evapotranspiration
- Actual ET = PET × Kc × f(θ)
- f(θ): linear reduction from θ_ETred (= wilting + fraction) to θ_wp → ET = 0
- ET reduction coefficients by land cover: 0.7–0.8 (Thornthwaite–Mather approach)

#### 5. Baseflow
- Deep percolation feeds a linear reservoir
- Recession coefficient set per watershed

#### 6. Road runoff
- Impervious road cells: 1.5 cm storage, excess → direct runoff

### Mass balance output (daily CSV)
Columns include: `wshed_id`, `date`, `year`, `runoff`, `precip`, `rain`, `ET`, `streamflow`, `SWE`, radiation flux components, soil moisture states

---

## Validation (functional_PP.R)

### Metrics
- **NSE**: Nash–Sutcliffe Efficiency (`hydroGOF::NSE`) — target > 0.5
- **KGE**: Kling–Gupta Efficiency (`hydroGOF::KGE`) — target > 0.5

### Standard plot set
1. Observed vs. modeled Q (log scale, daily)
2. Precip decomposed into rain + snow
3. SWE time series + SAM energy components
4. Basin-averaged daily water flux
5. Annual accumulated fluxes (bar chart)
6. Radiation budget (SW, LW, latent, sensible)

---

## Key Numerical Parameters

| Parameter | Value |
|-----------|-------|
| Cell size | 30 m |
| Stream threshold | 1000 cells (0.09 km²) |
| Sub-watershed threshold | 1000 cells |
| Soil Layer A depth | 0–30 cm |
| Max canopy storage (forest) | 0.31 cm |
| Road surface storage | 1.50 cm |
| Snow albedo (fresh) | 0.7383 |
| Snow albedo (aged) | 0.2 |
| T lapse rate | −3.20 °C / km |
| P lapse rate | +0.09195 cm / m |
| PET lapse (intercept) | −0.2486 mm |
| PET lapse (slope) | 1.8281 |

---

## R Package Dependencies

| Group | Packages |
|-------|----------|
| Geospatial | `terra`, `raster`, `sf`, `mapview`, `FedData` |
| Climate | `climateR`, `rnoaa`, `chillR`, `humidity`, `bigleaf` |
| Wrangling | `tidyverse`, `dplyr`, `lubridate` |
| Hydrology | `hydroGOF` |
| Visualization | `ggplot2`, `RColorBrewer`, `gridExtra`, `cowplot` |
| Terrain | `whitebox` |
| Spatial tools | `AOI` |

## External Tool Dependencies
- **WhiteBox Tools** — DEM conditioning (called via `whitebox` R package)
- **GRASS GIS 7+** — raster map algebra during Perl simulation (`r.mapcalc`, `r.sun`, `r.slope.aspect`)
- **Perl 5** — simulation engine runtime
