$ErrorActionPreference = 'Stop'
$perl = 'C:\Strawberry\perl\bin\perl.exe'
$engine = Join-Path (Split-Path -Parent $PSScriptRoot) 'R\smr_buildup.pl'

& $perl -c $engine
if ($LASTEXITCODE -ne 0) { throw 'Perl syntax check failed.' }

$previousErrorAction = $ErrorActionPreference
$ErrorActionPreference = 'Continue'
try {
  $invalid = & cmd.exe /d /s /c "set SMR_OUTPUT_MODE=invalid&& `"$perl`" `"$engine`" weather output" 2>&1
  $invalidExitCode = $LASTEXITCODE
} finally {
  $ErrorActionPreference = $previousErrorAction
}
if ($invalidExitCode -eq 0 -or ($invalid -join "`n") -notmatch 'SMR_OUTPUT_MODE must be') {
  throw 'Invalid SMR_OUTPUT_MODE did not fail before GRASS execution.'
}

$text = Get-Content -LiteralPath $engine -Raw
foreach ($required in @('SMR_OUTPUT_MODE', 'final-only', 'output_mode eq ''full''', 'final_state_maps=storage_amt,saturation')) {
  if ($text -notmatch [regex]::Escape($required)) { throw "Missing reporting-mode contract: $required" }
}

$runner = Get-Content -LiteralPath (Join-Path (Split-Path -Parent $PSScriptRoot) 'scripts\run_smr.ps1') -Raw
if ($runner -notmatch 'OutputMode' -or $runner -notmatch 'output_mode') { throw 'Runner does not expose and manifest OutputMode.' }
'reporting-mode checks OK'
