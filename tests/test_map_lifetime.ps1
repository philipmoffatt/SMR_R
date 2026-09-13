$ErrorActionPreference = 'Stop'
$engine = Join-Path (Split-Path -Parent $PSScriptRoot) 'R\smr_buildup.pl'
& 'C:\Strawberry\perl\bin\perl.exe' -c $engine
if ($LASTEXITCODE -ne 0) { throw 'Perl syntax check failed.' }
$text = Get-Content $engine -Raw
foreach ($required in @('daily_ephemeral_maps', 'g.remove -f type=raster', 'storage_amt_A_tmp', 'q.total')) {
  if ($text -notmatch [regex]::Escape($required)) { throw "Missing cleanup contract: $required" }
}
$cleanupBlock = [regex]::Match($text, 'my @daily_ephemeral_maps = qw\([\s\S]*?\);').Value
foreach ($state in @('storage_amt_A', 'storage_amt_B', 'storage_amt', 'swe.yesterday', 'saturation')) {
  if ($cleanupBlock -match "\b$([regex]::Escape($state))\b") { throw "Carry-over state listed for cleanup: $state" }
}
'map-lifetime checks OK'
