# PKT-SMR-CPU-05 CPU configuration and timing

`SMR_RSUN_NPROCS` controls `r.sun` workers and defaults to four. The runner
validates `-RSunNProcs` and records it in the stage manifest. Each simulated
day writes a TSV record for radiation, soil, reporting, and total wall time.
These records support a later repeated benchmark matrix; they do not claim an
optimal thread count from one run.
