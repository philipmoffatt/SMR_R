# PKT-SMR-CPU-02 reporting modes

`SMR_OUTPUT_MODE` is explicit and defaults to `full`.

- `full` preserves the legacy daily 24-field mass-balance record and all 19
  daily raster aggregation pairs.
- `daily` writes a compact daily state/flux row: watershed, date, year,
  runoff, precipitation, actual ET, storage, percolation, and Q. It avoids the
  diagnostic snow/canopy/energy aggregates.
- `final-only` creates no daily CSV rows and makes no daily zonal-statistic
  calls. It writes `<output>.summary.txt` with run dates and final-map names;
  final `storage_amt` and `saturation` maps remain in GRASS.

The PowerShell runner accepts `-OutputMode` and includes the selected mode in
its stage manifest. These modes change reporting cost only; state evolution and
scientific parameters are unchanged.
