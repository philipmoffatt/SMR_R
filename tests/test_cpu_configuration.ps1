$ErrorActionPreference = 'Stop'
$engine = Join-Path (Split-Path -Parent $PSScriptRoot) 'R\smr_buildup.pl'
& 'C:\Strawberry\perl\bin\perl.exe' -c $engine
if ($LASTEXITCODE -ne 0) { throw 'Perl syntax check failed.' }
$text = Get-Content $engine -Raw
foreach ($required in @('SMR_RSUN_NPROCS', 'nprocs=$rsun_nprocs', 'timing.tsv', 'radiation_seconds', 'soil_seconds')) {
  if ($text -notmatch [regex]::Escape($required)) { throw "Missing CPU timing contract: $required" }
}
$runner = Get-Content (Join-Path (Split-Path -Parent $PSScriptRoot) 'scripts\run_smr.ps1') -Raw
if ($runner -notmatch 'RSunNProcs' -or $runner -notmatch 'rsun_nprocs') { throw 'Runner CPU configuration is missing.' }
'CPU configuration checks OK'
