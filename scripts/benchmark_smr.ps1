[CmdletBinding()]
param(
  [Parameter(Mandatory)] [string]$SeedRuntime,
  [Parameter(Mandatory)] [string]$RuntimeRoot,
  [string]$GoldenRuntime,
  [ValidateSet('smoke', 'full')] [string]$Scenario = 'smoke',
  [string]$StartDate = '2015-01-01',
  [string]$EndDate,
  [switch]$SkipCompare
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot
$perlScript = Join-Path $repoRoot 'R\smr_buildup.pl'
$config = Join-Path $repoRoot 'config.mica-creek.yaml'
$o4wEnv = 'C:\OSGeo4W\bin\o4w_env.bat'
$grassEnv = 'C:\OSGeo4W\apps\grass\grass85\etc\env.bat'
$grassPython = 'C:\OSGeo4W\bin\python3.exe'
$grassScript = 'C:\OSGeo4W\apps\grass\grass85\etc\grass85.py'

if (!$EndDate) { $EndDate = if ($Scenario -eq 'full') { '2015-01-30' } else { '2015-01-02' } }
foreach ($path in @($SeedRuntime, $perlScript, $config, $o4wEnv, $grassEnv, $grassPython, $grassScript)) {
  if (!(Test-Path -LiteralPath $path)) { throw "Required benchmark path is missing: $path" }
}
if (Test-Path -LiteralPath $RuntimeRoot) {
  $existing = @(Get-ChildItem -LiteralPath $RuntimeRoot -Force -ErrorAction Stop)
  if ($existing.Count -gt 0) { throw "Benchmark runtime must be new or empty: $RuntimeRoot" }
}

function Invoke-TimedStage {
  param([Parameter(Mandatory)] [string]$Name, [Parameter(Mandatory)] [scriptblock]$Action)
  $watch = [Diagnostics.Stopwatch]::StartNew()
  & $Action
  $watch.Stop()
  [pscustomobject]@{ name = $Name; seconds = [math]::Round($watch.Elapsed.TotalSeconds, 3); exit_code = 0 }
}

function Invoke-Grass {
  param([Parameter(Mandatory)] [string[]]$GrassArguments)
  $quotedArgs = ($GrassArguments | ForEach-Object { '"' + $_.Replace('"', '\"') + '"' }) -join ' '
  $command = 'call "' + $o4wEnv + '" && call "' + $grassEnv + '" && "' + $grassPython + '" "' + $grassScript + '" ' + $quotedArgs
  $priorErrorAction = $ErrorActionPreference
  $ErrorActionPreference = 'Continue'
  try {
    $ignoredOutput = & cmd.exe /d /s /c $command 2>&1
    $exitCode = $LASTEXITCODE
    if ($exitCode -ne 0) { throw "GRASS command failed: $($GrassArguments -join ' ')" }
  } finally {
    $ErrorActionPreference = $priorErrorAction
  }
}

New-Item -ItemType Directory -Force -Path $RuntimeRoot | Out-Null
$tracePath = Join-Path $RuntimeRoot 'benchmark-trace.log'
Set-Content -LiteralPath $tracePath -Encoding utf8 -Value 'benchmark start'
$stages = @()
$stages += Invoke-TimedStage -Name 'seed_inputs' -Action {
  Copy-Item -LiteralPath (Join-Path $SeedRuntime 'processed_data\imitate_smr_setup') -Destination (Join-Path $RuntimeRoot 'processed_data\imitate_smr_setup') -Recurse
  $weatherDirectory = Join-Path $RuntimeRoot 'raw_data\weather'
  New-Item -ItemType Directory -Force -Path $weatherDirectory | Out-Null
  $sourceWeather = Join-Path $SeedRuntime 'raw_data\weather\gridmet_weather.csv'
  if (!(Test-Path -LiteralPath $sourceWeather)) { throw "Seed runtime has no weather input: $sourceWeather" }
  $start = [datetime]::ParseExact($StartDate, 'yyyy-MM-dd', $null)
  $end = [datetime]::ParseExact($EndDate, 'yyyy-MM-dd', $null)
  $selected = Get-Content -LiteralPath $sourceWeather | Where-Object {
    $parts = $_ -split '\s+'
    $parts.Count -gt 0 -and [datetime]::ParseExact($parts[0], 'yyyy-MM-dd', $null) -ge $start -and [datetime]::ParseExact($parts[0], 'yyyy-MM-dd', $null) -le $end
  }
  if (@($selected).Count -eq 0) { throw "No seed weather rows fall within $StartDate to $EndDate" }
  # Weather is numeric ASCII; ASCII avoids PowerShell 5.1's UTF-8 BOM, which
  # the Perl reader would otherwise retain on the first date token.
  Set-Content -LiteralPath (Join-Path $weatherDirectory 'gridmet_weather.csv') -Value $selected -Encoding ascii
}
Add-Content -LiteralPath $tracePath -Value 'seed complete'
Add-Content -LiteralPath $tracePath -Value 'import start'
$stages += Invoke-TimedStage -Name 'import' -Action {
  $inputDirectory = Join-Path $RuntimeRoot 'processed_data\imitate_smr_setup'
  $locationRoot = Join-Path $RuntimeRoot 'grass_workspace\MicaCreek_utm11'
  $location = Join-Path $locationRoot 'PERMANENT'
  if (!(Test-Path -LiteralPath (Join-Path $location 'DEFAULT_WIND'))) {
    New-Item -ItemType Directory -Force -Path (Split-Path -Parent $locationRoot) | Out-Null
    Invoke-Grass -GrassArguments @('-c', 'EPSG:32611', '-e', $locationRoot)
  }
  $maps = @('watershed', 'el', 'strms_30m', 'landuse', 'Ksat_matrix_A', 'Ksat_matrix_B', 'Ksat_mpore_A', 'Ksat_mpore_B', 'Kfc_A', 'Kfc_B', 'sat_mc_A', 'sat_mc_B', 'sat_amt', 'wiltpt_mc_A', 'wiltpt_mc_B', 'wiltpt_amt', 'fieldcap_mc_A', 'fieldcap_mc_B', 'fieldcap_amt_A', 'fieldcap_amt_B', 'fieldcap_amt', 'soil_depth_A', 'soil_depth_B', 'soil_depth', 'ETreduction_mc', 'Ksubsurface', 'flowunits', 'north', 'northeast', 'east', 'southeast', 'south', 'southwest', 'west', 'northwest')
  $aliases = @{ north = 'n'; northeast = 'ne'; east = 'e'; southeast = 'se'; south = 's'; southwest = 'sw'; west = 'w'; northwest = 'nw' }
  foreach ($map in $maps) {
    $sourceMap = if ($aliases.ContainsKey($map)) { $aliases[$map] } else { $map }
    $input = Join-Path $inputDirectory "$sourceMap.asc"
    if (!(Test-Path -LiteralPath $input)) { throw "Missing prepared spatial input: $input" }
    Invoke-Grass -GrassArguments @('--text', $location, '--exec', 'r.in.gdal', '-o', "input=$input", "output=$map", '--overwrite')
  }
}
Add-Content -LiteralPath $tracePath -Value 'import complete'
Add-Content -LiteralPath $tracePath -Value 'simulate start'
$stages += Invoke-TimedStage -Name 'simulate' -Action {
  $inputDirectory = Join-Path $RuntimeRoot 'processed_data\imitate_smr_setup'
  $location = Join-Path $RuntimeRoot 'grass_workspace\MicaCreek_utm11\PERMANENT'
  $weather = Join-Path $RuntimeRoot 'raw_data\weather\gridmet_weather.csv'
  $outputDirectory = Join-Path $RuntimeRoot 'smr_output'
  $output = Join-Path $outputDirectory 'MFC_mass_balance_1.csv'
  New-Item -ItemType Directory -Force -Path $outputDirectory | Out-Null
  Push-Location $inputDirectory
  try {
    Invoke-Grass -GrassArguments @('--text', $location, '--exec', 'C:\Strawberry\perl\bin\perl.exe', $perlScript, $weather, $output)
  } finally {
    Pop-Location
  }
}
Add-Content -LiteralPath $tracePath -Value 'simulate complete'

$manifest = [ordered]@{
  generated_utc = (Get-Date).ToUniversalTime().ToString('o')
  status = 'success'
  scenario = $Scenario
  start_date = $StartDate
  end_date = $EndDate
  seed_runtime = (Resolve-Path $SeedRuntime).Path
  runtime_root = (Resolve-Path $RuntimeRoot).Path
  source_revision = (git -C $repoRoot rev-parse HEAD).Trim()
  config_sha256 = (Get-FileHash -Algorithm SHA256 $config).Hash
  engine_sha256 = (Get-FileHash -Algorithm SHA256 $perlScript).Hash
  stages = @($stages)
  comparison = $null
}
if ($GoldenRuntime -and !$SkipCompare) {
  $comparisonPath = Join-Path $RuntimeRoot 'qc\comparison.json'
  & powershell.exe -NoProfile -ExecutionPolicy Bypass -File (Join-Path $PSScriptRoot 'compare_smr_outputs.ps1') -GoldenRuntime $GoldenRuntime -CandidateRuntime $RuntimeRoot -OutputPath $comparisonPath -ExpectedStart $StartDate -ExpectedEnd $EndDate
  if ($LASTEXITCODE -ne 0) { throw "Comparison stage failed with exit code $LASTEXITCODE" }
  $manifest.comparison = $comparisonPath
}
$manifestPath = Join-Path $RuntimeRoot 'benchmark-manifest.json'
$manifest | ConvertTo-Json -Depth 6 | Set-Content -Encoding utf8 -LiteralPath $manifestPath
Add-Content -LiteralPath $tracePath -Value 'manifest complete'
Write-Host "Benchmark complete: $manifestPath"
