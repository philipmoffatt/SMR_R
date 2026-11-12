"""Report two-day frozen-trace errors for the first JAX verification gate."""

import argparse
from pathlib import Path

import numpy as np

from .engine import run_days
from .initialization import initialize_legacy_state
from .io import FrozenTraceProvider, ReplayRadiationProvider, load_daily_forcings
from .legacy import load_legacy_bundle


_MAPS = {
    "precip": "precip_cm", "rain": "rain_cm", "snow": "snow_cm",
    "throughfall": "throughfall_cm", "canopy_evap": "canopy_evap_cm",
    "swe": "swe_cm", "snowmelt": "snowmelt_cm", "q.srad": "q_srad",
    "q.lw": "q_lw", "q.latent": "q_latent", "q.sensible": "q_sensible",
    "q.total": "q_total",
}


def comparison_rows(reference_dir: Path, days: int = 2):
    static = load_legacy_bundle(reference_dir / "model_inputs.zip")
    trace = FrozenTraceProvider(reference_dir / "two_day_internal_trace.zip", static.grid)
    forcings = load_daily_forcings(reference_dir / "gridmet_weather_2015-01.csv",
                                   ReplayRadiationProvider(reference_dir / "two_day_internal_trace.zip", static.grid), days)
    _, history = run_days(initialize_legacy_state(static), forcings, static)
    rows = []
    for day in range(days):
        date = f"2015-01-{day + 1:02d}"
        for trace_name, diagnostic_name in _MAPS.items():
            actual = np.asarray(history[diagnostic_name][day])
            expected = trace.raster(date, trace_name)
            mask = np.asarray(static.active_mask) & np.isfinite(expected) & np.isfinite(actual)
            error = actual[mask] - expected[mask]
            rows.append((date, trace_name, int(mask.sum()), float(np.sqrt(np.mean(error ** 2))), float(np.max(np.abs(error)))))
    return rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Compare first-gate JAX surface/snow fields to frozen trace maps.")
    parser.add_argument("--reference-dir", type=Path, default=Path("reference/mica-creek-jan2015"))
    args = parser.parse_args(argv)
    print("date\tmap\tcells\trmse\tmax_abs")
    for row in comparison_rows(args.reference_dir):
        print("%s\t%s\t%d\t%.17g\t%.17g" % row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
