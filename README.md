# SMR_R — Soil Moisture Routing with Snow Accumulation and Melt

A spatial watershed hydrology system that preprocesses geospatial data and drives a physics-based simulation of snowmelt, soil moisture, evapotranspiration, and streamflow. The model is tested on the Mica Fork Creek (MFC) watershed in the Moscow Mountain / Pullman region of Washington and Idaho.

---

## What This Does

The project chains three stages:

1. **Spatial data preparation** — downloads and processes elevation, soils, and land cover into 30 m ASCII raster grids that the simulation engine can read.
2. **Weather file generation** — builds daily/6-hourly climate forcing files from GRIDMET (2000–2022) or NOAA station data (1960–1979), including derived variables like PET, cloud fraction, and dew point.
3. **Simulation and validation** — runs a Perl-based hydrological model inside GRASS GIS and evaluates output against observed streamflow using NSE and KGE efficiency metrics.

---

## Repository Layout

```
R/scratch/
  SMR_data_setup.R   — spatial data download, DEM processing, soil & land cover rasterization
  wxData.R           — weather data download, disaggregation, PET calculation, file formatting
  functional_PP.R    — post-processing: observed vs. modeled streamflow, diagnostic plots
  smr_buildup.pl     — Perl simulation engine (daily timestep, 6-hour sub-steps for snow)
```

Output files land in `processed_data/ASC/MFC/` (ASCII rasters) and `imitate_smr_setup/` (consolidated model inputs).

---

## The Hydrological Model

The Perl script (`smr_buildup.pl`) runs a physically-based simulation on a 30 m grid. Each cell tracks:

- **Snow accumulation and melt (SAM)** — energy-balance approach: shortwave radiation (Bristow–Campbell transmissivity), longwave radiation (Unsworth–Monteith cloud emissivity), latent heat from vapor pressure gradients, and sensible heat from convective transfer. Albedo decays with snow age.
- **Canopy interception** — land-cover-specific maximum storage (0–0.31 cm); partitions precipitation into rain and snow using temperature thresholds.
- **Two-layer soil moisture** — Layer A: 0–30 cm; Layer B: 30 cm to restrictive horizon. Lateral flow uses Darcy's law with Bresler's formula for unsaturated conductivity. Saturation-excess runoff when soil moisture exceeds capacity.
- **Evapotranspiration** — actual ET is controlled by soil moisture, a crop coefficient (Kc), and PET. Kc varies by land cover class and season following FAO guidelines.
- **Baseflow** — linear recession reservoir fed by deep percolation.

The main time step is 24 hours; snow energy calculations sub-step at 6 hours.

---

## Spatial Data Pipeline (`SMR_data_setup.R`)

| Step | Tool | Output |
|------|------|--------|
| Download DEM (1 m NED), soils (SSURGO), land cover (NLCD 2019) | `FedData` | Raw rasters & shapefiles |
| Smooth DEM, breach depressions, compute D8 flow | WhiteBox Tools | Hydrologically conditioned DEM |
| Delineate watershed, extract streams (1000-cell threshold) | WhiteBox Tools | Basin & stream rasters |
| Select majority soil component per map unit, compute depth-weighted properties | R / `dplyr` | Ksat, bulk density, OM, θ_sat / θ_fc / θ_wp |
| Reclassify NLCD into 6 classes (water, urban, forest, shrub, grass, row crop) | `terra` | Land cover raster |
| Reproject everything to UTM 30 m, export ASCII grids | `terra` | Model-ready `.asc` files |

---

## Weather Data Pipeline (`wxData.R`)

**Modern period (2000–2022):** Downloaded from GRIDMET (University of Idaho). PET taken directly from the dataset.

**Historical period (1960–1979):** Downloaded from NOAA GHCND station USC00456789 (Pullman, WA). Missing variables are derived:
- **6-hourly temperature** — disaggregated from daily min/max using `chillR::make_hourly_temps()`
- **PET** — Hamon method (saturation vapor pressure × day length)
- **Dew point** — threshold-based formula (Gentilli 1955)
- **Cloud fraction** — observed / potential solar radiation

**Elevation corrections** applied to the station data before use:
- Temperature lapse rate: −3.20 °C / km
- Precipitation lapse rate: +0.09195 cm / m
- PET: linear regression model fitted to site

Each weather output file has 30 space-separated columns covering temperature (6-hourly), precipitation, PET, radiation, cloud fraction, crop coefficient, and heat-transfer resistance.

---

## Model Validation (`functional_PP.R`)

Simulated daily streamflow is compared to USGS gauge observations using:
- **NSE** — Nash–Sutcliffe Efficiency (`hydroGOF::NSE`)
- **KGE** — Kling–Gupta Efficiency (`hydroGOF::KGE`)

Diagnostic plots produced:
- Observed vs. modeled streamflow (log scale)
- Precipitation broken into rain vs. snow
- Snow water equivalent and SAM energy components
- Daily basin-averaged water fluxes
- Annual accumulated fluxes (bar charts)
- Radiation budget (shortwave, longwave, latent, sensible)

---

## Dependencies

**R packages**

| Category | Packages |
|----------|----------|
| Geospatial | `terra`, `raster`, `sf`, `mapview`, `FedData` |
| Climate / weather | `climateR`, `rnoaa`, `chillR`, `humidity`, `bigleaf` |
| Data wrangling | `tidyverse`, `dplyr`, `lubridate` |
| Hydrology metrics | `hydroGOF` |
| Visualization | `ggplot2`, `RColorBrewer`, `gridExtra`, `cowplot` |
| Terrain processing | `whitebox` |

**External tools**
- [WhiteBox Tools](https://www.whiteboxgeo.com/) — DEM hydrological conditioning
- [GRASS GIS 7+](https://grass.osgeo.org/) — raster map algebra during simulation (`r.mapcalc`, `r.sun`)
- Perl 5 — simulation engine runtime

**Remote data APIs**
- USGS National Elevation Dataset (NED)
- USDA SSURGO Web Soil Survey
- USGS National Land Cover Database (NLCD 2019)
- NOAA GHCND (station USC00456789, Pullman WA)
- University of Idaho GRIDMET archive

---

## JAX Prototype and NVIDIA GPU Setup

`src/smr_jax/` provides a float64 two-layer hydrology engine, an explicit
plant-water interface, prepared terrain/weather radiation and CPU/CUDA
verification. The agreed SMR-side prototype feature is complete under
[practical GPU acceptance](docs/packets/PKT-SMR-JAX-VERIFY-06.md).
The GRASS/Perl production path remains available; field validation, plant-side
adoption, arbitrary soil layers and Richards physics are separate work.

See the [plant-water interface](docs/smr-plant-interface.md) for requested,
realized and unmet extraction, independently checked water ledgers, adaptive
soil transport substeps and the [synthetic example](examples/plant_water.py).
The [feature closure](docs/packets/PKT-SMR-JAX-LAND-01.md) records acceptance
and deliberate deferrals. No plant repository changes are included.

For an NVIDIA GPU on Windows, use a Linux distribution under WSL2.  The native
Windows JAX build is CPU-only.  From the WSL checkout, create and activate an
isolated environment, then install the CUDA 12 JAX extra and this project:

```bash
python3 -m venv ~/venvs/smr-jax
source ~/venvs/smr-jax/bin/activate
python -m pip install --upgrade pip
python -m pip install --upgrade "jax[cuda12]"
python -m pip install -e .
```

Use float64 to match the legacy reference calculations, avoid JAX reserving
all GPU memory, and confirm that the GPU is visible:

```bash
export JAX_PLATFORMS=cuda
export JAX_ENABLE_X64=true
export XLA_PYTHON_CLIENT_PREALLOCATE=false
python -c "import jax; print(jax.devices())"
python -m unittest discover -s tests -p test_jax_foundation.py -v
```

The expected device listing includes `CudaDevice(id=0)`. The foundation test
loads the frozen reference grid and checks the explicit two-layer initial
moisture state. Flux and integrated checks are in the remaining JAX tests.

The modular engine's first runnable case is intentionally a two-day replay,
not a radiation replacement:

```bash
python -m smr_jax.cli --reference-dir reference/mica-creek-jan2015 --days 2
```

This diagnostic reads frozen `q.srad` trace fields while retaining the legacy
daily and four six-hour-step order. The Horn slope and receiver-weighted
routing have passed the independent Mica Creek source checks; that evidence
does not establish general terrain equivalence.

For simulation with evolving snow albedo, use the
[prepared radiation provider](docs/packets/PKT-SMR-JAX-VERIFY-05.md). GRASS
computes terrain/weather components during preparation; the JAX daily scan
applies the simulated albedo. An already prepared model run needs no GRASS
subprocess or completed shortwave replay. After preparing the requested dates:

```bash
python -m smr_jax.cli --days 30 --radiation-dir <prepared-directory> \
  --melt-small-snow-next-day --summary <summary.json>

JAX_PLATFORMS=cuda,cpu python -m smr_jax.verify_gpu \
  --radiation-dir <prepared-directory> --days 2 30 \
  --melt-small-snow-next-day --acceptance-profile practical \
  --report-path <new-report.json>
```

The optional snow policy releases positive start-of-day SWE <=0.001 cm as
melt on the next day. Practical GPU acceptance checks all active cells and
basin scalars at declared tolerances, retaining strict parity as a diagnostic.
The verifier explicitly preserves FCELL roundings during compilation; these
GPU results do not qualify an arbitrary default JIT configuration. Independent
conservation, exact masks/counters and required finite values remain checked.

---

## Running the Model

1. **Prepare spatial inputs** — open `R/scratch/SMR_data_setup.R`, set your watershed HUC code and output directory, then source the script. Expect several minutes of download and DEM processing time.
2. **Build the weather file** — open `R/scratch/wxData.R`, confirm the date range and station ID, then source it. The script writes a space-delimited text file ready for the Perl engine.
3. **Run the simulation** — execute `smr_buildup.pl` with GRASS GIS initialized and the parameter/weather files in place. Output is a daily CSV with one row per watershed cell.
4. **Validate results** — open `R/scratch/functional_PP.R`, point it at the simulation output and the observed streamflow file, then source it to generate efficiency metrics and plots.

---

## Key Parameters at a Glance

| Parameter | Value |
|-----------|-------|
| Grid cell size | 30 m |
| Main time step | 24 hours |
| Snow sub-step | 6 hours |
| Soil layer A depth | 0–30 cm |
| Soil layer B depth | 30 cm – restrictive horizon |
| Stream extraction threshold | 1000 cells (~0.09 km²) |
| Temperature lapse rate | −3.20 °C / km |
| Precipitation lapse rate | +0.09195 cm / m |
