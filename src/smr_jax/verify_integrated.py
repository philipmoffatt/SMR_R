"""Fail-closed integrated CPU parity: existing two-day or anchored 30-day trace."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from zipfile import ZipFile

import jax
import numpy as np

from .engine import run_days
from .fluxes.legacy_math import compile_legacy_cpu
from .initialization import initialize_legacy_state
from .io import FrozenTraceProvider, ReplayRadiationProvider, load_daily_forcings
from .legacy import load_legacy_bundle
from .legacy_cpu_math import compile_native_cpu, load_library
from .verify_surface import (RELATIVE_TOLERANCE, WATER_ATOL_CM, ENERGY_ATOL,
                             _MAPS as SURFACE_MAPS, verify_reference_manifest, verify_initialization)


DAILY_MAPS = {**SURFACE_MAPS, "throughfall": "soil_throughfall_cm", "snowmelt": "soil_snowmelt_cm",
    "storage_amt": "storage_amt_cm", "storage_amt_A": "storage_amt_a_cm", "storage_amt_B": "storage_amt_b_cm",
    "actualET_daily_flow": "actual_et_cm", "perc_daily_flow": "percolation_cm", "runoff_daily_flow": "runoff_cm",
    "saturation": "saturation_percent", "lateral_flow": "last_lateral_flow_cm"}
EXTRA_MAPS = {"canopy_storage_amt": "canopy_storage_cm", "ice.content": "ice_cm", "liquid.water": "liquid_cm",
    "snow.age": "snow_age_days", "u.total": "snow_total_energy", "u.surface": "snow_surface_energy",
    "tsnow_surf": "snow_surface_temp_c", "swe.yesterday": "previous_swe_cm", "road_runoff": "road_runoff_cm",
    "q.rain.ground": "q_rain_ground", "runoff_2015": "annual_runoff_cm", "Psat_2015": "annual_saturation_count"}
REPORT_COLUMNS = {3: "runoff_daily_flow", 4: "precip", 5: "rain", 6: "actualET_daily_flow", 7: "canopy_evap",
    8: "snowmelt", 9: "storage_amt", 10: "throughfall", 11: "canopy_storage_amt", 12: "perc_daily_flow",
    14: "swe", 15: "condens", 16: "snow", 18: "q.srad", 19: "q.latent", 20: "q.sensible", 21: "q.lw",
    22: "q.rain.ground", 23: "q.total"}
REPORT_MAPS = {**DAILY_MAPS, **EXTRA_MAPS, "condens": "condens_cm"}
_run_cpu = compile_legacy_cpu(run_days, static_argnames=("melt_small_snow_next_day",))
_run_native_cpu = compile_native_cpu(run_days, static_argnames=("melt_small_snow_next_day",))


def _budget(name):
    if name == "q.srad" or name in ("snow.age", "Psat_2015"):
        return 0., 0.
    if name.startswith("q.") or name in ("u.total", "u.surface"):
        return ENERGY_ATOL, RELATIVE_TOLERANCE
    if name in ("saturation", "tsnow_surf"):
        return 1e-10, RELATIVE_TOLERANCE
    return WATER_ATOL_CM, RELATIVE_TOLERANCE


def compare_field(name, actual, expected, mask):
    actual, expected, mask = np.asarray(actual), np.asarray(expected), np.asarray(mask)
    if mask.dtype != np.bool_ or mask.ndim != 2 or not mask.any():
        raise ValueError(f"{name}: required nonempty boolean support")
    if actual.shape != mask.shape or expected.shape != mask.shape or actual.dtype != np.float64 or expected.dtype != np.float64:
        raise ValueError(f"{name}: missing/incorrect shape or float64 type")
    if not np.isfinite(actual[mask]).all() or not np.isfinite(expected[mask]).all():
        raise ValueError(f"{name}: nonfinite required value")
    atol, rtol = _budget(name)
    np.testing.assert_allclose(actual[mask], expected[mask], rtol=rtol, atol=atol, err_msg=name)
    if name == "swe":
        np.testing.assert_array_equal(actual[mask] > 0., expected[mask] > 0., err_msg="snow-presence decisions")
    return float(np.max(np.abs(actual[mask]-expected[mask])))


def read_reports(text, days):
    rows = [line.split() for line in text.splitlines() if line.strip()]
    if len(rows) != days or any(len(row) != 24 or row[0] != "1" or row[1] != f"2015-01-{i+1:02d}" for i, row in enumerate(rows)):
        raise ValueError(f"expected {days} ordered full-report rows for watershed 1")
    values = np.asarray([[float(value) for j, value in enumerate(row) if j != 1] for row in rows], dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError("nonfinite full-report value")
    return rows


def validate_extension(capture_dir, reference_dir, static):
    """Anchor an independent capture BEFORE prescribing its shortwave fields."""
    capture_dir, reference_dir = Path(capture_dir), Path(reference_dir)
    metadata = json.loads((capture_dir/"manifest.json").read_text(encoding="utf-8"))
    root = Path(__file__).resolve().parents[2]
    if metadata.get("schema_version") != 1 or metadata.get("days") != 30:
        raise ValueError("capture must identify a complete thirty-day extension")
    for key, path in (("source_sha256", root/"R/smr_buildup.pl"),
                      ("capture_script_sha256", root/"tools/capture_jax_integrated_reference.py")):
        if sha256(path.read_bytes()).hexdigest() != metadata[key]:
            raise ValueError(f"capture provenance differs: {key}")
    frozen = json.loads((reference_dir/"manifest.json").read_text(encoding="utf-8-sig"))
    if metadata["frozen_artifacts"] != {item["path"]: item["sha256"] for item in frozen["artifacts"]}:
        raise ValueError("capture/frozen artifact provenance differs")
    if set(metadata["artifacts"]) != {"daily_trace.zip", "daily.csv", "capture-engine.pl"}:
        raise ValueError("incomplete capture artifact inventory")
    for name, item in metadata["artifacts"].items():
        path = capture_dir/name
        if path.stat().st_size != item["bytes"] or sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"capture artifact checksum/size failed: {name}")
    if set(metadata["daily_maps"]) != set(DAILY_MAPS)|set(EXTRA_MAPS):
        raise ValueError("incomplete or changed daily capture map inventory")
    provider = FrozenTraceProvider(capture_dir/"daily_trace.zip", static.grid)
    frozen_trace = FrozenTraceProvider(reference_dir/"two_day_internal_trace.zip", static.grid)
    mask = np.asarray(static.active_mask)
    for day in range(2):
        date = f"2015-01-{day+1:02d}"
        for name in DAILY_MAPS:
            actual, expected = provider.raster(date, name), frozen_trace.raster(date, name)
            support = np.isfinite(expected)
            if not np.array_equal(support, np.isfinite(actual)):
                raise ValueError(f"capture anchor support differs: {date} {name}")
            compare_field(name, actual, expected, support)
    terminal = FrozenTraceProvider(reference_dir/"final_state_rasters.zip", static.grid)
    for name in ("storage_amt", "saturation"):
        compare_field(name, provider.raster("2015-01-30", name), terminal.raster_member(f"golden_{name}.tif"), mask)
    cpu_rows = read_reports((capture_dir/"daily.csv").read_text(encoding="utf-8"), 30)
    legacy_rows = read_reports((reference_dir/"daily_mass_balance_2015-01.csv").read_text(encoding="utf-8"), 30)
    # The original zonal raster/report rounded aggregate values through FCELL
    # and six decimal places. The optimized CPU trace retains double sums.
    # Anchor that reporting distinction explicitly, without widening state or
    # daily-field tolerances or changing the actual model's basin state.
    volume, previous_baseflow = 0., 0.
    for cpu, legacy in zip(cpu_rows, legacy_rows):
        rounded = {column: float(format(float(np.float32(float(cpu[column]))), ".6f")) for column in REPORT_COLUMNS}
        for column, value in rounded.items():
            if value != float(legacy[column]):
                raise ValueError(f"legacy report anchor differs: {legacy[1]} column {column}")
        volume += rounded[12]/100./86400.-previous_baseflow
        previous_baseflow = .05*volume
        np.testing.assert_allclose([previous_baseflow+rounded[3]/100.*900./86400., previous_baseflow],
                                   [float(legacy[13]), float(legacy[17])],
                                   rtol=RELATIVE_TOLERANCE, atol=1e-15, err_msg="legacy report recurrence anchor")
    return provider, metadata


def compare_reports(history, rows, provider, soil_mask, initial_mask, *, extended):
    errors = []
    previous_swe = np.where(initial_mask, 0., np.nan)
    for day, row in enumerate(rows):
        available = {**DAILY_MAPS, **EXTRA_MAPS} if extended else DAILY_MAPS
        snow = provider.raster(row[1], "snow")
        latent = provider.raster(row[1], "q.latent")
        # Literal condens conditional: NULL condition stays NULL; the latent
        # operand is needed only on the snow-present arm. This support is
        # derived from independent trace operands, never actual finite values.
        condens_support = (np.isfinite(previous_swe) & np.isfinite(snow) &
                           ((previous_swe+snow <= 0.) | np.isfinite(latent)))
        for column, name in REPORT_COLUMNS.items():
            if name not in available and name != "condens":
                continue  # These maps exist only in the explicitly extended capture.
            mask = condens_support if name == "condens" else np.isfinite(provider.raster(row[1], name))
            if not mask.any() or np.any(soil_mask & ~mask):
                raise ValueError(f"missing required reporting support: {name}")
            # The two-day frozen trace does not include all state maps, but the
            # actual engine now exposes them for full-report comparison.
            if REPORT_MAPS[name] not in history:
                raise ValueError(f"missing required report diagnostic {REPORT_MAPS[name]}")
            value = np.asarray(history[REPORT_MAPS[name]][day])
            if value.shape != mask.shape or not np.isfinite(value[mask]).all():
                raise ValueError(f"nonfinite or missing report support: {name}")
            actual, expected = float(np.sum(value[mask], dtype=np.float64)), float(row[column])
            atol, _ = _budget(name)
            # Aggregate budget is the sum of declared cell floors; relative
            # allowance covers reduction and printed-report rounding.
            np.testing.assert_allclose(actual, expected, rtol=RELATIVE_TOLERANCE,
                                       atol=mask.sum()*atol, err_msg=f"{row[1]} report {name}")
            errors.append(abs(actual-expected))
        for key, column in (("discharge_cms", 13), ("baseflow_cms", 17)):
            actual = np.asarray(history[key][day])
            if actual.shape != (1,) or actual.dtype != np.float64 or not np.isfinite(actual).all():
                raise ValueError(f"invalid required basin diagnostic: {key}")
            np.testing.assert_allclose(actual, [float(row[column])], rtol=RELATIVE_TOLERANCE, atol=1e-15,
                                       err_msg=f"{row[1]} {key}")
        previous_swe = provider.raster(row[1], "swe")
    return float(max(errors))


def verification_report(reference_dir, *, days=2, capture_dir=None, collect_failures=False, native_math_library=None):
    if days not in (2, 30):
        raise ValueError("integrated gate requires exactly two or thirty days")
    reference_dir = Path(reference_dir)
    required_count = verify_reference_manifest(reference_dir)
    static = load_legacy_bundle(reference_dir/"model_inputs.zip")
    mask = np.asarray(static.active_mask)
    if int(mask.sum()) != required_count:
        raise ValueError("model required support differs from frozen inventory")
    if days == 30:
        if capture_dir is None:
            raise ValueError("thirty-day verification requires --capture-dir; missing days cannot be skipped")
        provider, metadata = validate_extension(capture_dir, reference_dir, static)
        trace_path = Path(capture_dir)/"daily_trace.zip"
        rows = read_reports((Path(capture_dir)/"daily.csv").read_text(encoding="utf-8"), days)
        maps = {**DAILY_MAPS, **EXTRA_MAPS}
    else:
        trace_path = reference_dir/"two_day_internal_trace.zip"
        provider = FrozenTraceProvider(trace_path, static.grid)
        with ZipFile(trace_path) as archive:
            rows = read_reports(archive.read("MFC_mass_balance_1.csv").decode("utf-8"), days)
        maps = DAILY_MAPS
    # Validate complete required support every day before using the trace.
    for day in range(days):
        date = f"2015-01-{day+1:02d}"
        support = np.isfinite(provider.raster(date, "saturation"))
        if not np.array_equal(mask, support):
            raise ValueError(f"{date}: required model/reference support differs")
    forcings = load_daily_forcings(reference_dir/"gridmet_weather_2015-01.csv", ReplayRadiationProvider(trace_path, static.grid), days)
    initial = initialize_legacy_state(static)
    verify_initialization(static, initial)
    math_profile = {"profile": "jax-exp"}
    run = _run_cpu
    if native_math_library is not None:
        math_profile = load_library(native_math_library)
        run = _run_native_cpu
    final, history = run(initial, forcings, static, melt_small_snow_next_day=False)
    if days == 30:
        for name, value in (("storage_amt", initial.soil.total_water_cm), ("storage_amt_A", initial.soil.layer_water_cm[0])):
            expected = provider.raster("initial", name)
            compare_field(name, np.asarray(value), expected, np.isfinite(expected))
    errors, failures = {}, []

    def observe(label, name, actual, expected, support):
        try:
            return compare_field(name, actual, expected, support)
        except AssertionError as exc:
            if not collect_failures:
                raise AssertionError(f"{label}: {exc}") from exc
            actual, expected = np.asarray(actual), np.asarray(expected)
            maximum = float(np.max(np.abs(actual[support]-expected[support])))
            failures.append({"field": label, "max_abs": maximum, "detail": str(exc)[:1000]})
            return maximum
    for day in range(days):
        date = f"2015-01-{day+1:02d}"
        for name, diagnostic in maps.items():
            if diagnostic not in history or history[diagnostic].shape != (days,)+static.grid.shape:
                raise ValueError(f"missing/incomplete integrated diagnostic: {diagnostic}")
            expected = provider.raster(date, name)
            field_mask = np.isfinite(expected)
            if np.any(mask & ~field_mask):
                raise ValueError(f"{date} {name}: missing required soil support")
            errors[f"{date} {name}"] = observe(f"{date} {name}", name, np.asarray(history[diagnostic][day]), expected, field_mask)
    try:
        report_error = compare_reports(history, rows, provider, mask, np.isfinite(np.asarray(static.watershed)), extended=days==30)
    except AssertionError as exc:
        if not collect_failures:
            raise
        report_error = None
        failures.append({"field": "daily reporting", "detail": str(exc)[:1000]})
    # Terminal state must be the actual final scan state, not just history.
    terminal_state = {"storage_amt": final.soil.total_water_cm, "storage_amt_A": final.soil.layer_water_cm[0],
        "storage_amt_B": final.soil.layer_water_cm[1], "swe": final.snow.swe_cm}
    if days == 30:
        terminal_state.update({"canopy_storage_amt": final.canopy_storage_cm, "ice.content": final.snow.ice_cm,
            "liquid.water": final.snow.liquid_cm, "snow.age": final.snow.age_days, "u.total": final.snow.total_energy,
            "u.surface": final.snow.surface_energy, "tsnow_surf": final.snow.surface_temp_c,
            "swe.yesterday": final.snow.previous_swe_cm, "runoff_2015": final.annual_runoff_cm,
            "Psat_2015": final.annual_saturation_count})
    if not np.array_equal(mask, np.asarray(final.soil.active_mask)):
        raise ValueError("terminal soil support changed")
    for name, value in terminal_state.items():
        observe(f"terminal {name}", name, np.asarray(value), provider.raster(f"2015-01-{days:02d}", name), mask)
    for name, actual, expected, atol in (("baseflow", final.basin.baseflow_cms, [float(rows[-1][17])], 1e-15),
        ("reservoir", final.basin.reservoir_volume, np.asarray([float(rows[-1][17])])/np.asarray(static.basin_reservoir_coeff), 1e-12)):
        try:
            np.testing.assert_allclose(np.asarray(actual), expected, atol=atol, rtol=RELATIVE_TOLERANCE, err_msg=f"terminal {name}")
        except AssertionError as exc:
            if not collect_failures:
                raise
            failures.append({"field": f"terminal {name}", "detail": str(exc)[:1000]})
    if days == 30:
        terminal = FrozenTraceProvider(reference_dir/"final_state_rasters.zip", static.grid)
        observe("original golden storage", "storage_amt", np.asarray(final.soil.total_water_cm), terminal.raster_member("golden_storage_amt.tif"), mask)
        observe("original golden saturation", "saturation", np.asarray(final.soil.total_water_cm/static.sat_amt_cm*100.), terminal.raster_member("golden_saturation.tif"), mask)
    return {"passed": not failures, "failures": failures, "math_profile": math_profile,
            "days": days, "required_cells": required_count, "daily_field_comparisons": len(errors),
            "max_water_error_cm": max(error for name, error in errors.items() if not name.split(" ", 1)[1].startswith("q.")
                                      and name.split(" ", 1)[1] not in ("u.total", "u.surface", "tsnow_surf", "saturation")),
            "max_field_errors": {name: max(error for key, error in errors.items() if key.endswith(" "+name)) for name in maps},
            "max_report_sum_error": report_error}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, default=Path("reference/mica-creek-jan2015"))
    parser.add_argument("--days", type=int, choices=(2, 30), default=2)
    parser.add_argument("--capture-dir", type=Path)
    parser.add_argument("--native-math-library", type=Path,
                        help="Explicit optional Windows/MinGW scalar exp FFI library and build record.")
    parser.add_argument("--audit-all", action="store_true", help="Collect numerical failures; exit remains nonzero if any check fails.")
    parser.add_argument("--report-path", type=Path, help="Write JSON to a NEW file; existing files are preserved.")
    args = parser.parse_args(argv)
    try:
        report = verification_report(args.reference_dir, days=args.days, capture_dir=args.capture_dir,
                                     collect_failures=args.audit_all, native_math_library=args.native_math_library)
        if args.report_path is not None:
            with args.report_path.open("x", encoding="utf-8", newline="\n") as stream:
                json.dump(report, stream, indent=2)
                stream.write("\n")
    except (AssertionError, ValueError, OSError, KeyError, RuntimeError, IndexError) as exc:
        print(f"FAIL: VERIFY-04: {exc}", file=sys.stderr)
        return 1
    print(("PASS" if report["passed"] else "FAIL")+": VERIFY-04 " + json.dumps(report, sort_keys=True))
    print("Strict legacy CPU profile with prescribed shortwave replay; general radiation and GPU acceptance remain separate.")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
