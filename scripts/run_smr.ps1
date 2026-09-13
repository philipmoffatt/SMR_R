[CmdletBinding()]
param(
  [ValidateSet("spatial", "weather", "import", "simulate", "validate")]
  [string]$Stage = "spatial",
  [string]$GridmetStart,
  [string]$GridmetEnd,
  [string]$TempDir = "C:\tmp\smr-first-run",
  [string]$ConfigPath,
  [string]$RuntimeRoot,
  [string]$Epsg = "32610",
  [string]$LocationName = "BullRun_utm10",
  [ValidateSet("full", "daily", "final-only")]
  [string]$OutputMode = "full",
  [switch]$QuietGrass
)

$ErrorActionPreference = "Stop"
$repoRoot = Split-Path -Parent $PSScriptRoot
$scratch = if ($RuntimeRoot) { $RuntimeRoot } else { Join-Path $repoRoot "scratch" }
$config = if ($ConfigPath) { $ConfigPath } else { Join-Path $repoRoot "config.example.yaml" }
$rscript = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe"
$perl = "C:\Strawberry\perl\bin\perl.exe"
$grassPython = "C:\OSGeo4W\bin\python3.exe"
$grassScript = "C:\OSGeo4W\apps\grass\grass85\etc\grass85.py"
$o4wEnv = "C:\OSGeo4W\bin\o4w_env.bat"
$grassEnv = "C:\OSGeo4W\apps\grass\grass85\etc\env.bat"

foreach ($path in @($config, $rscript, $perl, $grassPython, $grassScript, $o4wEnv, $grassEnv)) {
  if (!(Test-Path -LiteralPath $path)) { throw "Required path is missing: $path" }
}

New-Item -ItemType Directory -Force -Path $scratch, $TempDir | Out-Null
$env:SMR_CONFIG = $config
$env:SMR_ROOT_DIR = $scratch
$env:SMR_TEMP_DIR = $TempDir
# The machine also has PostgreSQL's older PROJ database on PATH.  GDAL/terra
# must use the matching OSGeo4W database for raster reprojection.
$env:PROJ_LIB = "C:\OSGeo4W\share\proj"
$env:SMR_OUTPUT_MODE = $OutputMode
if ($GridmetStart) { $env:SMR_GRIDMET_START = $GridmetStart }
if ($GridmetEnd) { $env:SMR_GRIDMET_END = $GridmetEnd }

function Write-RunManifest([string]$Status, [int]$ExitCode) {
  $manifestDir = Join-Path $scratch "run-manifests"
  New-Item -ItemType Directory -Force -Path $manifestDir | Out-Null
  [ordered]@{
    timestamp_utc = (Get-Date).ToUniversalTime().ToString("o")
    stage = $Stage
    status = $Status
    exit_code = $ExitCode
    repo_root = $repoRoot
    runtime_root = $scratch
    config = $config
    gridmet_start = $env:SMR_GRIDMET_START
    gridmet_end = $env:SMR_GRIDMET_END
    output_mode = $OutputMode
  } | ConvertTo-Json | Set-Content -Encoding utf8 (Join-Path $manifestDir ((Get-Date -Format "yyyyMMddTHHmmss") + "_$Stage.json"))
}

function Invoke-Grass([string[]]$GrassArgs) {
  $quotedArgs = ($GrassArgs | ForEach-Object { '"' + $_.Replace('"', '\"') + '"' }) -join ' '
  $command = 'call "' + $o4wEnv + '" && call "' + $grassEnv + '" && "' +
    $grassPython + '" "' + $grassScript + '" ' + $quotedArgs
  # GRASS writes ordinary progress notices to stderr on Windows. Merge it with
  # stdout so ErrorActionPreference does not misclassify a successful launch.
  $previousErrorAction = $ErrorActionPreference
  $ErrorActionPreference = "Continue"
  try {
    if ($QuietGrass) {
      & cmd.exe /d /s /c $command 2>&1 | Out-Null
      if ($LASTEXITCODE -ne 0) { throw "GRASS command failed: $($GrassArgs -join ' ')" }
    } else {
      & cmd.exe /d /s /c $command 2>&1
    }
  } finally {
    $ErrorActionPreference = $previousErrorAction
  }
}

try {
  switch ($Stage) {
    "spatial" {
      & $rscript (Join-Path $repoRoot "R\run_pipeline.R") spatial
    }
    "weather" {
      & $rscript (Join-Path $repoRoot "R\run_pipeline.R") weather
    }
    "import" {
      $inputDir = Join-Path $scratch "processed_data\imitate_smr_setup"
      $locationRoot = Join-Path $scratch (Join-Path "grass_workspace" $LocationName)
      $location = Join-Path $locationRoot "PERMANENT"
      if (!(Test-Path -LiteralPath $inputDir)) { throw "Run spatial before import: $inputDir" }
      if (!(Test-Path -LiteralPath (Join-Path $location "DEFAULT_WIND"))) {
        New-Item -ItemType Directory -Force -Path (Split-Path -Parent $locationRoot) | Out-Null
        Invoke-Grass @("-c", "EPSG:$Epsg", "-e", $locationRoot)
        if ($LASTEXITCODE -ne 0) { throw "GRASS location creation failed." }
      }
      $maps = @(
        "watershed", "el", "strms_30m", "landuse",
        "Ksat_matrix_A", "Ksat_matrix_B", "Ksat_mpore_A", "Ksat_mpore_B",
        "Kfc_A", "Kfc_B", "sat_mc_A", "sat_mc_B", "sat_amt",
        "wiltpt_mc_A", "wiltpt_mc_B", "wiltpt_amt",
        "fieldcap_mc_A", "fieldcap_mc_B", "fieldcap_amt_A", "fieldcap_amt_B", "fieldcap_amt",
        "soil_depth_A", "soil_depth_B", "soil_depth", "ETreduction_mc", "Ksubsurface",
        "flowunits", "north", "northeast", "east", "southeast", "south", "southwest", "west", "northwest"
      )
      $sourceAliases = @{
        north = "n"; northeast = "ne"; east = "e"; southeast = "se"
        south = "s"; southwest = "sw"; west = "w"; northwest = "nw"
      }
      foreach ($map in $maps) {
        $sourceMap = if ($sourceAliases.ContainsKey($map)) { $sourceAliases[$map] } else { $map }
        $input = Join-Path $inputDir "$sourceMap.asc"
        if (!(Test-Path -LiteralPath $input)) { throw "Missing spatial input: $input" }
        Invoke-Grass @("--text", $location, "--exec", "r.in.gdal", "-o", "input=$input", "output=$map", "--overwrite")
        if ($LASTEXITCODE -ne 0) { break }
      }
    }
    "simulate" {
      $inputDir = Join-Path $scratch "processed_data\imitate_smr_setup"
      $location = Join-Path $scratch (Join-Path (Join-Path "grass_workspace" $LocationName) "PERMANENT")
      $weather = Join-Path $scratch "raw_data\weather\gridmet_weather.csv"
      $outputDir = Join-Path $scratch "smr_output"
      $output = Join-Path $outputDir "MFC_mass_balance_1.csv"
      foreach ($path in @($inputDir, $location, $weather)) {
        if (!(Test-Path -LiteralPath $path)) { throw "Required simulation input is missing: $path" }
      }
      New-Item -ItemType Directory -Force -Path $outputDir | Out-Null
      # smr_buildup.pl appends output. A stage run must begin with exactly the
      # requested weather interval, never rows retained from a prior attempt.
      if (Test-Path -LiteralPath $output) { Remove-Item -LiteralPath $output -Force }
      Push-Location $inputDir
      try {
        Invoke-Grass @("--text", $location, "--exec", $perl, (Join-Path $repoRoot "R\smr_buildup.pl"), $weather, $output)
      } finally {
        Pop-Location
      }
    }
    "validate" {
      & $rscript (Join-Path $repoRoot "R\run_pipeline.R") validate
    }
  }
  $code = $LASTEXITCODE
  if ($null -eq $code) { $code = 0 }
  Write-RunManifest $(if ($code -eq 0) { "success" } else { "failed" }) $code
  exit $code
} catch {
  Write-RunManifest "failed" 1
  throw
}
