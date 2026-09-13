# PKT-SMR-CPU-01A benchmark harness

`scripts/benchmark_smr.ps1` is the CPU reference runner. It creates a new
runtime, seeds it from a preserved runtime's already-prepared spatial inputs,
filters its weather file to the requested interval, and times the GRASS import
and simulation stages. It does not use the completed reference runtime as a
write target.

Example two-day smoke benchmark:

```powershell
scripts\benchmark_smr.ps1 `
  -SeedRuntime C:\projects\git\worktrees\SMR-FIRST-RUN-01\scratch\mica-creek-jan2015 `
  -RuntimeRoot C:\projects\git\worktrees\PKT-SMR-CPU-01A\scratch\bench-smoke-20150101 `
  -Scenario smoke
```

For a regression comparison, use `-Scenario full` and supply the same completed
reference runtime as `-GoldenRuntime`. The comparison reports raw state values;
it does not mask known legacy negative saturation values or claim scientific
validation. The JSON benchmark manifest records the committed source revision,
configuration hash, runner hash, runtime paths, dates, and per-stage seconds.

## Completion evidence

The isolated two-day smoke run (2015-01-01 through 2015-01-02) completed with
two output rows. Its runtime-only manifest recorded 0.235 s to seed prepared
inputs, 34.866 s for GRASS import, and 97.659 s for simulation. The runtime is
in the system temporary directory and is intentionally not versioned.
