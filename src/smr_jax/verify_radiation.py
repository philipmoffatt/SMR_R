"""VERIFY-05 component evidence and end-to-end dynamic-albedo diagnostics.

Radiation acceptance uses a predeclared FCELL forward-error bound. Integrated
strict legacy differences remain visible under the operator's VERIFY-04
deferral and are never recast as a strict pass.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

from .contracts import RadiationComponents
from .engine import run_days
from .fluxes.legacy_math import compile_legacy_cpu
from .initialization import initialize_legacy_state
from .io import ReplayRadiationProvider, load_daily_forcings
from .legacy import load_legacy_bundle
from .radiation import PreparedRadiationProvider, shortwave_from_components, radiation_error_bound
from .verify_integrated import (DAILY_MAPS, EXTRA_MAPS, compare_field, compare_reports,
                                read_reports, validate_extension)
from .verify_surface import verify_reference_manifest, RELATIVE_TOLERANCE


_shortwave = compile_legacy_cpu(shortwave_from_components)
_run = compile_legacy_cpu(run_days, static_argnames=("melt_small_snow_next_day",))


def component_report(fixture_dir):
    fixture_dir = Path(fixture_dir)
    metadata = json.loads((fixture_dir/"jax_verify05_radiation.json").read_text(encoding="utf-8"))
    path = fixture_dir/"jax_verify05_radiation.npz"
    if metadata.get("schema_version") != 1 or path.stat().st_size != metadata["artifact_bytes"] or sha256(path.read_bytes()).hexdigest() != metadata["artifact_sha256"]:
        raise ValueError("radiation evidence checksum/inventory differs")
    with np.load(path, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    shape = (len(metadata["cases"]), metadata["representatives"])
    combinations = {(c["doy"], c["cloud"], c["albedo"]) for c in metadata["cases"]}
    required_cases = {(day, cloud, albedo) for day in (1, 80, 172, 355)
                      for cloud in (0., .5, 1.) for albedo in (0., .2, .7383, .95, 1.)}
    if combinations != required_cases or len(metadata["cases"]) != len(required_cases):
        raise ValueError("radiation evidence is missing required season/cloud/albedo cases")
    required = {"beam", "diffuse", "reflected_unit", "albedo", "canopy", "expected_q", "expected_reflected", "indices"}
    if set(arrays) != required or metadata["representatives"] < 1:
        raise ValueError("radiation evidence fields/representatives differ")
    for name, value in arrays.items():
        if name != "indices" and (value.shape != shape or value.dtype != np.float64 or not np.isfinite(value).all()):
            raise ValueError(f"invalid required radiation evidence: {name}")
    components = RadiationComponents(*(jnp.asarray(arrays[k]) for k in ("beam", "diffuse", "reflected_unit")))
    albedo, canopy = arrays["albedo"], arrays["canopy"]
    actual = np.asarray(_shortwave(components, jnp.asarray(albedo), jnp.asarray(canopy)))
    errors = np.abs(actual-arrays["expected_q"])
    bound = radiation_error_bound(components, albedo, canopy)
    if np.any(errors > bound):
        raise AssertionError("radiation exceeds declared FCELL forward-error bound")
    return {"passed": True, "cases": shape[0], "representatives_per_case": shape[1],
            "max_shortwave_error_kj_m2_day": float(errors.max()),
            "max_fraction_of_bound": float(np.max(errors/bound))}


def compare_reference_albedo(reference, forcings, static, days):
    """Separate component validation at independently captured albedos.

    These reference states are used in this comparison only. They are never
    supplied to run_days or substituted into its evolving state.
    """
    previous_swe = np.zeros(static.grid.shape, dtype=np.float64)
    maximum, fraction, cells = 0., 0., 0
    for day in range(days):
        when = f"2015-01-{day+1:02d}"
        age, snow = reference.raster(when, "snow.age"), reference.raster(when, "snow")
        albedo = np.where(previous_swe+snow > 0., np.minimum(.95, .7383*age**(-.1908)), .2)
        components = jax.tree_util.tree_map(lambda value: value[day], forcings.radiation_components)
        expected = reference.raster(when, "q.srad")
        support = np.isfinite(expected)
        actual = np.asarray(_shortwave(components, jnp.asarray(albedo), static.canopy_cover))
        bound = radiation_error_bound(components, albedo, static.canopy_cover)
        if not np.isfinite(actual[support]).all() or not np.isfinite(bound[support]).all():
            raise ValueError(f"{when}: missing component comparison support")
        errors = np.abs(actual[support]-expected[support])
        if np.any(errors > bound[support]):
            raise AssertionError(f"{when}: radiation components exceed FCELL forward-error bound")
        maximum = max(maximum, float(errors.max()))
        fraction = max(fraction, float(np.max(errors/bound[support])))
        cells += int(support.sum())
        previous_swe = reference.raster(when, "swe")
    return {"passed": True, "cell_days": cells, "max_shortwave_error_kj_m2_day": maximum,
            "max_fraction_of_bound": fraction}


def integrated_report(reference_dir, radiation_dir, capture_dir, *, days=30, melt_small_snow_next_day=False,
                      compare_replay_control=False):
    if days not in (2, 30):
        raise ValueError("radiation integration requires two or thirty days")
    reference_dir = Path(reference_dir)
    required_cells = verify_reference_manifest(reference_dir)
    static = load_legacy_bundle(reference_dir/"model_inputs.zip")
    reference, _ = validate_extension(capture_dir, reference_dir, static)
    prepared = PreparedRadiationProvider(radiation_dir, static.grid, reference_dir/"model_inputs.zip", required_mask=static.active_mask)
    forcings = load_daily_forcings(reference_dir/"gridmet_weather_2015-01.csv", prepared, days)
    if forcings.replay_shortwave is not None or forcings.radiation_components is None:
        raise ValueError("integrated radiation must not prescribe completed shortwave")
    component_comparison = compare_reference_albedo(reference, forcings, static, days)
    final, history = _run(initialize_legacy_state(static), forcings, static,
                          melt_small_snow_next_day=melt_small_snow_next_day)
    failures, maxima, snow_presence = [], {}, 0
    mask = np.asarray(static.active_mask)
    maps = {**DAILY_MAPS, **EXTRA_MAPS}
    for day in range(days):
        when = f"2015-01-{day+1:02d}"
        for name, key in maps.items():
            expected = reference.raster(when, name)
            support = np.isfinite(expected)
            actual = np.asarray(history[key][day])
            if np.any(mask & ~support) or actual.shape != support.shape or not np.isfinite(actual[support]).all():
                raise ValueError(f"{when} {name}: missing/nonfinite required support")
            error = float(np.max(np.abs(actual[support]-expected[support])))
            maxima[name] = max(maxima.get(name, 0.), error)
            try:
                compare_field(name, actual, expected, support)
            except AssertionError:
                failures.append({"date": when, "field": name, "max_abs": error})
            if name == "swe":
                snow_presence += int(np.count_nonzero((actual[support] > 0.) != (expected[support] > 0.)))
    # The stored terminal carry must still be the final actual history state.
    for state, key in ((final.soil.total_water_cm, "storage_amt_cm"), (final.snow.swe_cm, "swe_cm")):
        np.testing.assert_array_equal(np.asarray(state)[mask], np.asarray(history[key][-1])[mask])
    rows = read_reports((Path(capture_dir)/"daily.csv").read_text(encoding="utf-8"), 30)[:days]
    reporting_failure = None
    try:
        compare_reports(history, rows, reference, mask, np.isfinite(np.asarray(static.watershed)), extended=True)
    except AssertionError as exc:
        reporting_failure = str(exc)[:1000]
    basin_errors, basin_failures = {}, []
    for key, column in (("discharge_cms", 13), ("baseflow_cms", 17)):
        actual = np.asarray(history[key])
        expected = np.asarray([[float(row[column])] for row in rows])
        if actual.shape != expected.shape or not np.isfinite(actual).all():
            raise ValueError(f"incomplete/nonfinite basin history: {key}")
        basin_errors[key] = float(np.max(np.abs(actual-expected)))
        try:
            np.testing.assert_allclose(actual, expected, atol=1e-15, rtol=RELATIVE_TOLERANCE)
        except AssertionError:
            basin_failures.append(key)
    terminal = {"storage_amt": final.soil.total_water_cm, "storage_amt_A": final.soil.layer_water_cm[0],
                "storage_amt_B": final.soil.layer_water_cm[1], "swe": final.snow.swe_cm,
                "canopy_storage_amt": final.canopy_storage_cm, "ice.content": final.snow.ice_cm,
                "liquid.water": final.snow.liquid_cm, "snow.age": final.snow.age_days,
                "u.total": final.snow.total_energy, "u.surface": final.snow.surface_energy,
                "tsnow_surf": final.snow.surface_temp_c, "swe.yesterday": final.snow.previous_swe_cm,
                "runoff_2015": final.annual_runoff_cm, "Psat_2015": final.annual_saturation_count}
    terminal_errors, terminal_failures = {}, []
    for name, value in terminal.items():
        actual = np.asarray(value)
        expected = reference.raster(f"2015-01-{days:02d}", name)
        if not np.isfinite(actual[mask]).all() or not np.isfinite(expected[mask]).all():
            raise ValueError(f"nonfinite terminal state: {name}")
        terminal_errors[name] = float(np.max(np.abs(actual[mask]-expected[mask])))
        try:
            compare_field(name, actual, expected, mask)
        except AssertionError:
            terminal_failures.append(name)
    for state, key in ((final.basin.baseflow_cms, "baseflow_cms"), (final.basin.reservoir_volume, "basin_reservoir")):
        np.testing.assert_array_equal(np.asarray(state), np.asarray(history[key][-1]))
    expected_volume = np.asarray([float(rows[-1][17])])/np.asarray(static.basin_reservoir_coeff)
    basin_errors["terminal_reservoir"] = float(np.max(np.abs(np.asarray(final.basin.reservoir_volume)-expected_volume)))
    try:
        np.testing.assert_allclose(np.asarray(final.basin.reservoir_volume), expected_volume, atol=1e-12, rtol=RELATIVE_TOLERANCE)
    except AssertionError:
        basin_failures.append("terminal_reservoir")
    control_report = None
    if compare_replay_control:
        # Independent diagnostic run with the SAME snow policy. Its state is
        # never fed to the target model; this separates provider changes from
        # the already-authorized policy's deliberate differences from Perl.
        replay = load_daily_forcings(reference_dir/"gridmet_weather_2015-01.csv",
            ReplayRadiationProvider(Path(capture_dir)/"daily_trace.zip", static.grid), days)
        _, control = _run(initialize_legacy_state(static), replay, static,
                         melt_small_snow_next_day=melt_small_snow_next_day)
        control_errors = {}
        for name, key in maps.items():
            # Some thermal fields have wider first-day support than later
            # days. This optional control compares the fixed, independently
            # established model domain; the primary comparison above retains
            # every reference field's complete day-specific footprint.
            support = mask
            target, baseline = np.asarray(history[key]), np.asarray(control[key])
            if not np.isfinite(baseline[:, support]).all():
                raise ValueError(f"nonfinite same-policy replay control: {name}")
            control_errors[name] = float(np.max(np.abs(target[:, support]-baseline[:, support])))
        control_report = {"same_snow_policy": True, "required_cells": int(mask.sum()),
                          "max_field_errors": control_errors,
                          "note": "Separate diagnostic control only; target state evolution uses prepared components."}
    return {"executed_without_shortwave_replay": True, "days": days, "required_cells": required_cells,
            "radiation_at_reference_albedo": component_comparison,
            "snow_policy": "next_day_melt_1e-3_cm" if melt_small_snow_next_day else "legacy",
            "strict_legacy_parity_passed": not (failures or reporting_failure or basin_failures or terminal_failures),
            "strict_failure_count": len(failures),
            "strict_failures": failures, "max_field_errors": maxima,
            "strict_reporting_passed": reporting_failure is None, "strict_reporting_first_failure": reporting_failure,
            "basin_max_errors": basin_errors, "terminal_max_errors": terminal_errors,
            "strict_basin_failures": basin_failures, "strict_terminal_failures": terminal_failures,
            "same_policy_replay_control": control_report,
            "snow_presence_disagreements_cell_days": snow_presence,
            "note": "Numerical legacy differences are reported under the operator-authorized VERIFY-04 deferral."}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture-dir", type=Path, default=Path("tests/fixtures"))
    parser.add_argument("--reference-dir", type=Path, default=Path("reference/mica-creek-jan2015"))
    parser.add_argument("--radiation-dir", type=Path)
    parser.add_argument("--capture-dir", type=Path)
    parser.add_argument("--days", type=int, choices=(2, 30), default=30)
    parser.add_argument("--melt-small-snow-next-day", action="store_true")
    parser.add_argument("--compare-replay-control", action="store_true",
                        help="Separate diagnostic run to isolate radiation differences under the same snow policy.")
    parser.add_argument("--report-path", type=Path)
    args = parser.parse_args(argv)
    try:
        report = {"radiation_components": component_report(args.fixture_dir)}
        if args.radiation_dir is not None:
            if args.capture_dir is None:
                raise ValueError("integration requires the independent --capture-dir")
            report["integration"] = integrated_report(args.reference_dir, args.radiation_dir, args.capture_dir,
                days=args.days, melt_small_snow_next_day=args.melt_small_snow_next_day,
                compare_replay_control=args.compare_replay_control)
        if args.report_path is not None:
            with args.report_path.open("x", encoding="utf-8") as stream:
                json.dump(report, stream, indent=2)
                stream.write("\n")
    except (AssertionError, ValueError, KeyError, OSError, RuntimeError) as exc:
        print(f"FAIL: VERIFY-05: {exc}")
        return 1
    print(json.dumps(report, sort_keys=True))
    print("Radiation component acceptance is distinct from the reported, deferred strict integrated legacy comparison.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
