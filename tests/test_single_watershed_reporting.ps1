$ErrorActionPreference = 'Stop'
$engine = Join-Path (Split-Path -Parent $PSScriptRoot) 'R\smr_buildup.pl'
& 'C:\Strawberry\perl\bin\perl.exe' -c $engine
if ($LASTEXITCODE -ne 0) { throw 'Perl syntax check failed.' }
$text = Get-Content $engine -Raw
foreach ($required in @('my $single_watershed', 'r.univar -g map=$var', 'r.stats.zonal base=watershed')) {
  if ($text -notmatch [regex]::Escape($required)) { throw "Missing reporting route: $required" }
}
if ($text.IndexOf('if ($single_watershed)') -lt 0 -or $text.IndexOf('} else {') -lt 0) {
  throw 'Single and multi-watershed reporting routes are not both present.'
}
'single-watershed routing checks OK'
