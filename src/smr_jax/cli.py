"""Runnable two-day frozen-trace replay entry point."""

import argparse
import json
from pathlib import Path

import numpy as np

from .engine import run_days
from .initialization import initialize_legacy_state
from .io import ReplayRadiationProvider, load_daily_forcings
from .legacy import load_legacy_bundle


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Run the SMR JAX frozen replay-radiation case.")
    parser.add_argument("--reference-dir", type=Path, default=Path("reference/mica-creek-jan2015"))
    parser.add_argument("--days", type=int, default=2, help="Rows to run; only days with replay q.srad are supported.")
    parser.add_argument("--summary", type=Path, default=Path("smr-jax-replay-summary.json"))
    args = parser.parse_args(argv)
    if args.days < 1 or args.days > 2:
        parser.error("the frozen replay provider contains exactly two days; use --days 1 or 2")
    static = load_legacy_bundle(args.reference_dir / "model_inputs.zip")
    forcings = load_daily_forcings(args.reference_dir / "gridmet_weather_2015-01.csv",
                                   ReplayRadiationProvider(args.reference_dir / "two_day_internal_trace.zip", static.grid), args.days)
    final, history = run_days(initialize_legacy_state(static), forcings, static)
    summary = {"days": args.days, "active_cells": int(np.asarray(static.active_mask).sum()),
               "storage_sum_cm": float(np.nansum(np.asarray(final.soil.total_water_cm))),
               "runoff_sum_cm": float(np.nansum(np.asarray(final.annual_runoff_cm))),
               "baseflow_cms": np.asarray(final.basin.baseflow_cms).tolist()}
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
