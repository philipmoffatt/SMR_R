$ErrorActionPreference = 'Stop'
. (Join-Path (Split-Path -Parent $PSScriptRoot) 'scripts\compare_smr_outputs.ps1')
$testRoot = Join-Path ([System.IO.Path]::GetTempPath()) ('smr-benchmark-test-' + [guid]::NewGuid())
try {
  New-Item -ItemType Directory -Force -Path $testRoot | Out-Null
  $complete = Join-Path $testRoot 'complete.csv'
  $gapped = Join-Path $testRoot 'gapped.csv'
  Set-Content -LiteralPath $complete -Value @('1 2015-01-01 2015 0', '1 2015-01-02 2015 0') -Encoding utf8
  Set-Content -LiteralPath $gapped -Value @('1 2015-01-01 2015 0', '1 2015-01-03 2015 0') -Encoding utf8
  $good = Get-SmrCsvSummary -Path $complete -ExpectedStart '2015-01-01' -ExpectedEnd '2015-01-02'
  if (!$good.complete -or $good.rows -ne 2) { throw 'Expected complete two-day CSV summary.' }
  $bad = Get-SmrCsvSummary -Path $gapped -ExpectedStart '2015-01-01' -ExpectedEnd '2015-01-03'
  if ($bad.complete -or '2015-01-02' -notin $bad.missing_dates) { throw 'Expected missing-date detection.' }
  $comparison = Compare-SmrCsvSummaries -Golden $good -Candidate $good
  if (!$comparison.same_date_sequence) { throw 'Expected equal date sequence.' }
  'benchmark CSV checks OK'
} finally {
  if (Test-Path -LiteralPath $testRoot) { Remove-Item -LiteralPath $testRoot -Recurse -Force }
}
