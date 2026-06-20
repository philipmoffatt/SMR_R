# SMR_R — Environment Setup

## What You Need Before Running

| Tool | Stage | Status |
|------|-------|--------|
| R 4.5+ | 1, 2, 4 | installed |
| WhiteBox Tools binary | 1 | installed (auto-detected by `whitebox` R pkg) |
| Strawberry Perl | 3 | **needs install** |
| GRASS GIS 8 | 3 | **needs install** |

---

## 1. Install Strawberry Perl (Windows)

1. Download from https://strawberryperl.com/ — grab the 64-bit MSI installer.
2. Run the installer with all defaults. It adds `perl` to your PATH.
3. Verify: open a new PowerShell and run `perl --version`. You should see `v5.3x.x`.

---

## 2. Install GRASS GIS 8 (Windows)

1. Download the standalone installer from https://grass.osgeo.org/download/windows/
   - Use the **OSGeo4W** installer for the easiest setup, or the standalone GRASS 8 MSI.
   - OSGeo4W path: https://trac.osgeo.org/osgeo4w/ → click "OSGeo4W Network Installer"
   - In the installer choose "Express Install" → select **GRASS GIS**.
2. After install, GRASS is at something like `C:\OSGeo4W\bin\grass.bat`.
3. Add `C:\OSGeo4W\bin` to your system PATH (or use the OSGeo4W Shell).
4. Verify: `grass --version` should print `GRASS GIS 8.x.x`.

---

## 3. GRASS Location (already created)

The GRASS location `scratch/grass_workspace/BullRun/PERMANENT` was created automatically
in UTM Zone 10N (EPSG:32610). No manual setup needed.

### How to run GRASS commands on Windows (OSGeo4W)

GRASS on Windows requires OSGeo4W environment variables. Use this pattern from PowerShell:

```powershell
$o4w = "C:\OSGeo4W"
$env:OSGEO4W_ROOT    = $o4w
$env:GISBASE         = "$o4w\apps\grass\grass85"
$env:GRASS_PYTHON    = "$o4w\bin\python3.exe"
$env:PYTHONHOME      = "$o4w\apps\Python312"
$env:PYTHONPATH      = "$o4w\apps\grass\grass85\etc\python;$o4w\apps\Python312\Lib;$o4w\apps\Python312\DLLs"
$env:GRASS_PROJSHARE = "$o4w\share\proj"
$env:PATH            = "$o4w\bin;$o4w\apps\grass\grass85\bin;$o4w\apps\Python312;$env:PATH"

# Then call GRASS:
& "$o4w\bin\python3.exe" "$o4w\apps\grass\grass85\etc\grass85.py" `
    --text scratch/grass_workspace/BullRun/PERMANENT --exec <command>
```

Or use the OSGeo4W Shell (Start Menu → OSGeo4W Shell) and run `grass85` directly.

### Import the parameter rasters before running the Perl script

Once `SMR_data_setup.R` has produced the `.asc` files in `processed_data/imitate_smr_setup/`,
import them into GRASS from inside the GRASS shell:

```bash
# Start GRASS in the BullRun location
grass --text C:\projects\git\SMR_R\scratch\grass_workspace\BullRun\PERMANENT

# Set region from the watershed mask
r.in.gdal input=processed_data/imitate_smr_setup/watershed.asc output=watershed --overwrite
g.region rast=watershed

# Import all parameter rasters (run from the scratch directory)
r.in.gdal input=processed_data/imitate_smr_setup/el.asc              output=el              --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/strms_30m.asc       output=strms_30m       --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/landuse.asc         output=landuse         --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/Ksat_matrix_A.asc   output=Ksat_matrix_A   --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/Ksat_matrix_B.asc   output=Ksat_matrix_B   --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/Ksat_mpore_A.asc    output=Ksat_mpore_A    --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/Ksat_mpore_B.asc    output=Ksat_mpore_B    --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/Kfc_A.asc           output=Kfc_A           --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/Kfc_B.asc           output=Kfc_B           --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/sat_mc_A.asc        output=sat_mc_A        --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/sat_mc_B.asc        output=sat_mc_B        --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/sat_amt.asc         output=sat_amt         --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/wiltpt_mc_A.asc     output=wiltpt_mc_A     --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/wiltpt_mc_B.asc     output=wiltpt_mc_B     --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/wiltpt_amt.asc      output=wiltpt_amt      --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/fieldcap_mc_A.asc   output=fieldcap_mc_A   --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/fieldcap_mc_B.asc   output=fieldcap_mc_B   --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/fieldcap_amt_A.asc  output=fieldcap_amt_A  --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/fieldcap_amt_B.asc  output=fieldcap_amt_B  --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/fieldcap_amt.asc    output=fieldcap_amt    --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/soil_depth_A.asc    output=soil_depth_A    --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/soil_depth_B.asc    output=soil_depth_B    --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/soil_depth.asc      output=soil_depth      --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/ETreduction_mc.asc  output=ETreduction_mc  --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/Ksubsurface.asc     output=Ksubsurface     --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/n.asc               output=north           --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/ne.asc              output=northeast       --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/e.asc               output=east            --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/se.asc              output=southeast       --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/s.asc               output=south           --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/sw.asc              output=southwest       --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/w.asc               output=west            --overwrite
r.in.gdal input=processed_data/imitate_smr_setup/nw.asc              output=northwest       --overwrite

exit
```

---

## 4. Run the Pipeline

Once everything above is installed and the GRASS location is populated:

```r
# From R, with working directory set to scratch/
source("../R/scratch/SMR_data_setup.R")   # Stage 1 — ~10-20 min download + DEM
source("../R/scratch/wxData.R")            # Stage 2 — ~2 min
# Stage 3 — run from the GRASS shell:
#   cd C:\projects\git\SMR_R\scratch
#   grass --text grass_workspace/BullRun/PERMANENT --exec perl ../R/scratch/smr_buildup.pl
source("../R/scratch/functional_PP.R")     # Stage 4 — plots + NSE/KGE
```

---

## Key Parameter Choices for Bull Run

| Parameter | Value | Notes |
|-----------|-------|-------|
| Watershed | Bull Run River near Bull Run, OR | USGS 14138870, CAMELS |
| HUC12 | 170900120402 | verify with nhdplusTools output |
| CRS | EPSG:32610 (UTM Zone 10N) | Oregon Cascades |
| Period | WY2015–WY2019 | GRIDMET only, no historical station files |
| `basin_mask_value` | 1 | **confirm after `wbt_basins` runs** — check plot output |
