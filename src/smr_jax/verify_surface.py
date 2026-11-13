"""Fail-closed two-day initialization/surface/snow gate on an explicit CPU profile."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import jax
import numpy as np

from .engine import run_days
from .fluxes.legacy_math import compile_legacy_cpu
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

# Fixed float64 numerical budgets, not tolerances fitted to failing fields.
# The relative budget allows up to 512 double-precision rounding units through
# the surface expressions; the unit-specific absolute floors cover near zero.
RELATIVE_TOLERANCE = 512 * np.finfo(np.float64).eps
WATER_ATOL_CM = 1e-12
ENERGY_ATOL = 1e-9
_run_legacy_cpu = compile_legacy_cpu(run_days, static_argnames=("melt_small_snow_next_day",))


def verify_reference_manifest(reference_dir):
    """Verify the frozen inventory before using any of its artifacts as evidence."""
    root = Path(reference_dir).resolve()
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8-sig"))
    if manifest.get("reference_id") != "mica-creek-jan2015-v1":
        raise ValueError("VERIFY-01 requires the identified frozen Mica Creek reference")
    artifacts = manifest.get("artifacts", [])
    names = [item["path"] for item in artifacts]
    required = {"model_inputs.zip", "gridmet_weather_2015-01.csv", "two_day_internal_trace.zip"}
    if len(names) != len(set(names)) or not required.issubset(names):
        raise ValueError("frozen manifest has duplicate or missing required artifacts")
    for item in artifacts:
        path = (root / item["path"]).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise ValueError(f"missing or invalid frozen artifact: {item['path']}")
        if path.stat().st_size != item["bytes"] or hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"frozen artifact checksum/size failed: {item['path']}")
    count = manifest.get("parity_evidence", {}).get("common_valid_cells")
    if not isinstance(count, int) or count <= 0:
        raise ValueError("frozen manifest is missing its required comparison support")
    return count


def compare_field(actual, expected, mask, *, name, assert_close=True):
    """Compare every required cell; never use a finite-intersection fallback."""
    actual, expected, mask = np.asarray(actual), np.asarray(expected), np.asarray(mask)
    if mask.dtype != np.bool_ or mask.ndim != 2 or not mask.any():
        raise ValueError(f"{name}: required support must be a nonempty boolean grid")
    if actual.shape != mask.shape or expected.shape != mask.shape:
        raise ValueError(f"{name}: actual/reference shapes must equal required grid {mask.shape}")
    if actual.dtype != np.float64 or expected.dtype != np.float64:
        raise ValueError(f"{name}: comparison requires float64 arrays")
    if not np.isfinite(actual[mask]).all() or not np.isfinite(expected[mask]).all():
        raise ValueError(f"{name}: missing/nonfinite value on required support")
    if name not in _MAPS:
        raise ValueError(f"unknown required surface diagnostic: {name}")
    error = actual[mask] - expected[mask]
    if assert_close:
        atol = 0. if name == "q.srad" else ENERGY_ATOL if name.startswith("q.") else WATER_ATOL_CM
        rtol = 0. if name == "q.srad" else RELATIVE_TOLERANCE
        np.testing.assert_allclose(actual[mask], expected[mask], rtol=rtol, atol=atol, err_msg=name)
    return int(mask.sum()), float(np.sqrt(np.mean(error ** 2))), float(np.max(np.abs(error)))


def verify_initialization(static, state):
    """Check the literal initial expressions and surface land-cover constants."""
    mask = np.asarray(static.active_mask)
    if not np.array_equal(mask, np.asarray(state.soil.active_mask)):
        raise ValueError("initial soil support differs from static support")

    def check(name, actual, expected):
        actual = np.asarray(actual)
        expected = np.broadcast_to(expected, mask.shape)
        if actual.shape != mask.shape or actual.dtype != np.float64 or not np.isfinite(actual[mask]).all():
            raise ValueError(f"invalid initialization field: {name}")
        np.testing.assert_allclose(actual[mask], expected[mask], rtol=RELATIVE_TOLERANCE,
                                   atol=WATER_ATOL_CM, err_msg=f"initial {name}")

    total = np.asarray(static.sat_amt_cm) * .4
    a = (np.asarray(static.wiltpt_mc_a) + (np.asarray(static.fieldcap_mc_a) -
         np.asarray(static.wiltpt_mc_a)) * .2) * np.asarray(static.soil_depth_a_cm)
    check("total water", state.soil.total_water_cm, total)
    check("A water", state.soil.layer_water_cm[0], a)
    check("B water", state.soil.layer_water_cm[1], total - a)
    check("layer sum", np.asarray(state.soil.layer_water_cm).sum(axis=0), total)
    check("canopy storage", state.canopy_storage_cm, 0.)
    for name, expected in (("swe_cm", 0.), ("previous_swe_cm", 0.), ("ice_cm", 0.),
                           ("liquid_cm", 0.), ("age_days", 1.), ("surface_temp_c", -2.),
                           ("total_energy", -10080.), ("surface_energy", -2016.)):
        check(name, getattr(state.snow, name), expected)
    lu = np.asarray(static.landuse)
    for name, values in (("max_canopy_cm", (0., .1, .31, .15, .15, .147)),
                          ("canopy_cover", (0., 0., .8, .6, .6, .5)),
                          ("tmax_rain_c", (-.5, -.5, 3.1, 1.8, -.5, 1.8)),
                          ("tmin_snow_c", (-3., -3., .9, -1., -3., -1.))):
        expected = np.select([lu == k for k in range(1, 7)], values, default=0.)
        check(name, getattr(static, name), expected)


def comparison_rows(reference_dir: Path, days: int = 2, *, melt_small_snow_next_day=False,
                    assert_close=True):
    if days not in (1, 2):
        raise ValueError("the frozen daily trace contains exactly January 1 and 2")
    if melt_small_snow_next_day and assert_close:
        raise ValueError("the selected snow convention is not the frozen legacy target; request report-only")
    reference_dir = Path(reference_dir)
    expected_count = verify_reference_manifest(reference_dir)
    with jax.default_device(jax.devices("cpu")[0]):
        static = load_legacy_bundle(reference_dir / "model_inputs.zip")
        trace = FrozenTraceProvider(reference_dir / "two_day_internal_trace.zip", static.grid)
        forcings = load_daily_forcings(reference_dir / "gridmet_weather_2015-01.csv",
                                       ReplayRadiationProvider(reference_dir / "two_day_internal_trace.zip", static.grid), days)
        initial = initialize_legacy_state(static)
        verify_initialization(static, initial)
        _, history = _run_legacy_cpu(initial, forcings, static,
                                     melt_small_snow_next_day=melt_small_snow_next_day)
    rows = []
    for day in range(days):
        date = f"2015-01-{day + 1:02d}"
        # The frozen engine's finite saturation footprint independently owns
        # comparison support. It is never supplied to the actual model physics.
        mask = np.isfinite(trace.raster(date, "saturation"))
        if int(mask.sum()) != expected_count or not np.array_equal(mask, np.asarray(static.active_mask)):
            raise ValueError(f"{date}: model/frozen support differs from the required {expected_count} cells")
        for trace_name, diagnostic_name in _MAPS.items():
            if diagnostic_name not in history:
                raise ValueError(f"missing required diagnostic: {diagnostic_name}")
            actual = np.asarray(history[diagnostic_name][day])
            expected = trace.raster(date, trace_name)
            rows.append((date, trace_name, *compare_field(actual, expected, mask, name=trace_name,
                                                         assert_close=assert_close)))
            if assert_close and day == 0 and trace_name == "swe":
                np.testing.assert_array_equal(actual[mask] > 0., expected[mask] > 0.,
                                               err_msg="day-one SWE decisions carried into day two")
    return rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Verify initialization and two-day surface/snow parity using the strict float64 CPU profile.")
    parser.add_argument("--reference-dir", type=Path, default=Path("reference/mica-creek-jan2015"))
    parser.add_argument("--melt-small-snow-next-day", action="store_true",
                        help="Report the explicit 1e-3 cm next-day-melt departure from frozen legacy behavior.")
    parser.add_argument("--report-only", action="store_true",
                        help="Print differences without asserting numerical parity; required for the changed snow convention.")
    args = parser.parse_args(argv)
    if args.melt_small_snow_next_day and not args.report_only:
        parser.error("--melt-small-snow-next-day requires --report-only against this unmodified legacy trace")
    try:
        rows = comparison_rows(args.reference_dir, melt_small_snow_next_day=args.melt_small_snow_next_day,
                                assert_close=not args.report_only)
    except (AssertionError, ValueError, OSError, KeyError, RuntimeError) as exc:
        print(f"FAIL: VERIFY-01: {exc}", file=sys.stderr)
        return 1
    print("date\tmap\tcells\trmse\tmax_abs")
    for row in rows:
        print("%s\t%s\t%d\t%.17g\t%.17g" % row)
    print("REPORT ONLY: no numerical parity acceptance" if args.report_only else
          "PASS: initialization and 24 daily surface/snow comparisons; strict float64 CPU profile")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
