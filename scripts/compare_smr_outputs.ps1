[CmdletBinding()]
param(
  [string]$GoldenRuntime,
  [string]$CandidateRuntime,
  [string]$OutputPath,
  [string]$ExpectedStart,
  [string]$ExpectedEnd,
  [string[]]$Maps = @('saturation', 'storage_amt')
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Get-SmrCsvSummary {
  param(
    [Parameter(Mandatory)] [string]$Path,
    [string]$ExpectedStart,
    [string]$ExpectedEnd
  )

  if (!(Test-Path -LiteralPath $Path)) { throw "Missing SMR output CSV: $Path" }
  $dates = New-Object System.Collections.Generic.List[string]
  $badRows = New-Object System.Collections.Generic.List[int]
  $lineNumber = 0
  foreach ($line in Get-Content -LiteralPath $Path) {
    $lineNumber++
    if ([string]::IsNullOrWhiteSpace($line)) { continue }
    $fields = $line -split '\s+'
    if ($fields.Count -lt 2 -or $fields[1] -notmatch '^\d{4}-\d{2}-\d{2}$') {
      $badRows.Add($lineNumber)
      continue
    }
    $dates.Add($fields[1])
  }

  $distinct = @($dates | Sort-Object -Unique)
  $missing = New-Object System.Collections.Generic.List[string]
  if ($ExpectedStart -and $ExpectedEnd) {
    $cursor = [datetime]::ParseExact($ExpectedStart, 'yyyy-MM-dd', $null)
    $end = [datetime]::ParseExact($ExpectedEnd, 'yyyy-MM-dd', $null)
    while ($cursor -le $end) {
      $date = $cursor.ToString('yyyy-MM-dd')
      if ($date -notin $dates) { $missing.Add($date) }
      $cursor = $cursor.AddDays(1)
    }
  }

  [pscustomobject]@{
    path = $Path
    rows = $dates.Count
    first_date = if ($dates.Count) { $dates[0] } else { $null }
    last_date = if ($dates.Count) { $dates[$dates.Count - 1] } else { $null }
    duplicate_dates = @($dates | Group-Object | Where-Object Count -gt 1 | ForEach-Object Name)
    malformed_rows = @($badRows)
    missing_dates = @($missing)
    complete = ($badRows.Count -eq 0 -and $missing.Count -eq 0 -and $dates.Count -eq $distinct.Count)
    date_sequence = @($dates)
  }
}

function Compare-SmrCsvSummaries {
  param([Parameter(Mandatory)]$Golden, [Parameter(Mandatory)]$Candidate)
  [pscustomobject]@{
    golden = $Golden
    candidate = $Candidate
    same_date_sequence = (@($Golden.date_sequence) -join ',') -eq (@($Candidate.date_sequence) -join ',')
  }
}

function Get-GrassLocation {
  param([Parameter(Mandatory)] [string]$Runtime)
  $workspace = Join-Path $Runtime 'grass_workspace'
  $location = Get-ChildItem -LiteralPath $workspace -Directory -ErrorAction Stop |
    Where-Object { Test-Path -LiteralPath (Join-Path $_.FullName 'PERMANENT') } |
    Select-Object -First 1
  if ($null -eq $location) { throw "No GRASS location with a PERMANENT mapset under $workspace" }
  Join-Path $location.FullName 'PERMANENT'
}

function Export-GrassRaster {
  param([Parameter(Mandatory)] [string]$Location, [Parameter(Mandatory)] [string]$Map, [Parameter(Mandatory)] [string]$Output)
  $o4w = 'C:\OSGeo4W\bin\o4w_env.bat'
  $grassEnv = 'C:\OSGeo4W\apps\grass\grass85\etc\env.bat'
  $python = 'C:\OSGeo4W\bin\python3.exe'
  $grass = 'C:\OSGeo4W\apps\grass\grass85\etc\grass85.py'
  foreach ($path in @($o4w, $grassEnv, $python, $grass)) {
    if (!(Test-Path -LiteralPath $path)) { throw "Required GRASS path missing: $path" }
  }
  $command = 'call "' + $o4w + '" && call "' + $grassEnv + '" && "' + $python + '" "' + $grass +
    '" --text "' + $Location + '" --exec r.out.gdal input=' + $Map + ' output="' + $Output +
    '" format=GTiff createopt=COMPRESS=DEFLATE --overwrite'
  & cmd.exe /d /s /c $command 2>&1 | Out-Null
  if ($LASTEXITCODE -ne 0) { throw "GRASS export failed for $Map from $Location" }
}

function Get-RasterComparison {
  param([Parameter(Mandatory)] [string]$GoldenTif, [Parameter(Mandatory)] [string]$CandidateTif, [Parameter(Mandatory)] [string]$OutputCsv)
  $rscript = 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe'
  if (!(Test-Path -LiteralPath $rscript)) { throw "Rscript is missing: $rscript" }
  $code = @'
library(terra)
args <- commandArgs(trailingOnly = TRUE)
a <- rast(args[[1]])
b <- rast(args[[2]])
if (!compareGeom(a, b, stopOnError = FALSE)) stop('Raster geometry differs')
av <- values(a, mat = FALSE)
bv <- values(b, mat = FALSE)
keep <- is.finite(av) & is.finite(bv)
d <- bv[keep] - av[keep]
stats <- function(x, prefix) c(setNames(c(min(x), max(x), mean(x), sum(is.finite(x))), paste0(prefix, c('_min','_max','_mean','_valid_cells'))))
row <- c(stats(av, 'golden'), stats(bv, 'candidate'),
         difference_min=if(length(d)) min(d) else NA_real_,
         difference_max=if(length(d)) max(d) else NA_real_,
         difference_mean=if(length(d)) mean(d) else NA_real_,
         difference_rmse=if(length(d)) sqrt(mean(d^2)) else NA_real_,
         common_valid_cells=sum(keep))
write.csv(data.frame(metric=names(row), value=as.numeric(row)), args[[3]], row.names=FALSE)
'@
  & $rscript -e $code $GoldenTif $CandidateTif $OutputCsv
  if ($LASTEXITCODE -ne 0) { throw "Raster comparison failed for $GoldenTif and $CandidateTif" }
  Import-Csv -LiteralPath $OutputCsv
}

if ($MyInvocation.InvocationName -ne '.') {
  foreach ($argument in @(
    [pscustomobject]@{ name = 'GoldenRuntime'; value = $GoldenRuntime },
    [pscustomobject]@{ name = 'CandidateRuntime'; value = $CandidateRuntime },
    [pscustomobject]@{ name = 'OutputPath'; value = $OutputPath }
  )) {
    if ([string]::IsNullOrWhiteSpace($argument.value)) { throw "$($argument.name) is required." }
  }
  $outDirectory = Split-Path -Parent $OutputPath
  New-Item -ItemType Directory -Force -Path $outDirectory | Out-Null
  $goldenCsv = Get-SmrCsvSummary -Path (Join-Path $GoldenRuntime 'smr_output\MFC_mass_balance_1.csv') -ExpectedStart $ExpectedStart -ExpectedEnd $ExpectedEnd
  $candidateCsv = Get-SmrCsvSummary -Path (Join-Path $CandidateRuntime 'smr_output\MFC_mass_balance_1.csv') -ExpectedStart $ExpectedStart -ExpectedEnd $ExpectedEnd
  $result = [ordered]@{
    generated_utc = (Get-Date).ToUniversalTime().ToString('o')
    csv = Compare-SmrCsvSummaries -Golden $goldenCsv -Candidate $candidateCsv
    rasters = @()
  }
  $goldenLocation = Get-GrassLocation -Runtime $GoldenRuntime
  $candidateLocation = Get-GrassLocation -Runtime $CandidateRuntime
  foreach ($map in $Maps) {
    $goldenTif = Join-Path $outDirectory "golden_$map.tif"
    $candidateTif = Join-Path $outDirectory "candidate_$map.tif"
    $statsCsv = Join-Path $outDirectory "$map`_comparison.csv"
    Export-GrassRaster -Location $goldenLocation -Map $map -Output $goldenTif
    Export-GrassRaster -Location $candidateLocation -Map $map -Output $candidateTif
    $result.rasters += [ordered]@{ map = $map; metrics = @(Get-RasterComparison -GoldenTif $goldenTif -CandidateTif $candidateTif -OutputCsv $statsCsv) }
  }
  $result | ConvertTo-Json -Depth 8 | Set-Content -Encoding utf8 -LiteralPath $OutputPath
  Write-Host "Comparison written to $OutputPath"
}
