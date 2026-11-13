"""Strict VERIFY-03: lateral routing, GRASS slope defaults and basin recurrence."""

import argparse
import ast
from dataclasses import fields, replace
from hashlib import sha256
import json
from io import StringIO
from pathlib import Path
import re
import sys
from zipfile import ZipFile

import jax.numpy as jnp
import numpy as np

from .fluxes.basin import basin_recurrence
from .fluxes.legacy_math import compile_legacy_cpu
from .fluxes.routing import legacy_receive_gather
from .fluxes.soil import _effective_conductivity, lateral_route
from .initialization import initialize_legacy_state
from .io import FrozenTraceProvider
from .legacy import load_legacy_bundle
from .state import BasinState, GridSpec, TwoLayerState
from .verify_soil import _evaluate, check_values
from .verify_surface import verify_reference_manifest


_DIRECTIONS = ((-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1))
_SOURCE_FIELDS = {"slope": "slope_percent", "flowunits": "flowunits", "soil_depth_A": "soil_depth_a_cm",
    "soil_depth_B": "soil_depth_b_cm", "fieldcap_amt_A": "fieldcap_a_cm", "fieldcap_amt_B": "fieldcap_b_cm",
    "sat_mc_A": "sat_mc_a", "sat_mc_B": "sat_mc_b", "Ksat_matrix_A": "ksat_matrix_a", "Ksat_matrix_B": "ksat_matrix_b",
    "Ksat_mpore_A": "ksat_mpore_a", "Ksat_mpore_B": "ksat_mpore_b", "Kfc_A": "kfc_a", "Kfc_B": "kfc_b"}
_run_lateral = compile_legacy_cpu(lateral_route)
_run_basin = compile_legacy_cpu(basin_recurrence)
_run_gather = compile_legacy_cpu(legacy_receive_gather)
_run_conductivity = compile_legacy_cpu(_effective_conductivity)


def _neighbor(values, dr, dc):
    # Reference indexing is deliberately independent of JAX pad/slice routing.
    h, w = values.shape
    result = np.full((h, w), np.nan, dtype=np.float64)
    rows, cols = np.indices((h, w))
    inside = (rows+dr >= 0) & (rows+dr < h) & (cols+dc >= 0) & (cols+dc < w)
    result[inside] = values[rows[inside]+dr, cols[inside]+dc]
    return result


class RoutingOracle:
    """Read the ten actual Perl assignments, including literal neighbor offsets."""

    def __init__(self, source):
        text = Path(source).read_text(encoding="utf-8")
        start, end = "# Lateral flow (Darcy + Bresler)", "# Re-distribute between A and B"
        if text.count(start) != 1 or text.count(end) != 1:
            raise ValueError("missing/ambiguous source lateral block")
        block = text.split(start, 1)[1].split(end, 1)[0]
        expressions = re.findall(r'expression=([^"\r\n]+)', block)
        if len(expressions) != 10:
            raise ValueError("source lateral block must contain ten assignments")
        self.neighbors, self.assignments = {}, []

        def neighbor(match):
            key = (match[1], int(match[2]), int(match[3]))
            alias = f"neighbor_{len(self.neighbors)}"
            self.neighbors[alias] = key
            return alias

        for expression in expressions:
            target, value = expression.split(" = ", 1)
            value = value.replace("$temp_time_step", "6.0").replace("$gridsize", "30.0")
            value = re.sub(r"(lateral_flow(?:_[AB])?)\[(-?\d+),(-?\d+)\]", neighbor, value)
            value = re.sub(r"\bif\(", "choose(", value)
            self.assignments.append((target, ast.parse(value, mode="eval").body))

    def run(self, soil, static):
        mask = np.asarray(static.active_mask)
        env = {name: np.where(mask, np.asarray(getattr(static, key), dtype=np.float32), np.float32(np.nan)) for name, key in _SOURCE_FIELDS.items()}
        env.update(storage_amt=np.where(mask, np.asarray(soil.total_water_cm), np.nan),
                   storage_amt_A=np.where(mask, np.asarray(soil.layer_water_cm[0]), np.nan),
                   storage_amt_B=np.where(mask, np.asarray(soil.layer_water_cm[1]), np.nan))
        for name, value in zip(("north", "northeast", "east", "southeast", "south", "southwest", "west", "northwest"),
                                np.asarray(static.routing_weights)):
            env[name] = value
        with np.errstate(invalid="ignore", divide="ignore"):
            for i, (target, expression) in enumerate(self.assignments):
                if i == 7:
                    env.update({alias: _neighbor(env[name], dr, dc) for alias, (name, dr, dc) in self.neighbors.items()})
                env[target] = _evaluate(expression, env)
        return env


def _lateral_values(soil, static):
    after, flow = _run_lateral(soil, static)
    values = {"storage_amt": after.total_water_cm, "storage_amt_A": after.layer_water_cm[0],
              "storage_amt_B": after.layer_water_cm[1], "lateral_flow": flow}
    for i, letter in enumerate(("A", "B")):
        lower = letter.lower()
        eff = _run_conductivity(soil.layer_water_cm[i], getattr(static, f"fieldcap_{lower}_cm"),
              getattr(static, f"sat_mc_{lower}"), getattr(static, f"soil_depth_{lower}_cm"),
              getattr(static, f"ksat_matrix_{lower}"), getattr(static, f"ksat_mpore_{lower}"), getattr(static, f"kfc_{lower}"))
        out = jnp.minimum(eff * static.slope_percent / 100. * getattr(static, f"soil_depth_{lower}_cm") / 30. / 100., soil.layer_water_cm[i])
        values[f"effK_{letter}"] = eff
        values[f"lateral_flow_{letter}"] = jnp.where(static.flowunits == 0., 0., out)
    return {name: np.asarray(value) for name, value in values.items()}


def check_lateral(oracle, soil, static):
    mask = np.asarray(static.active_mask)
    actual, expected = _lateral_values(soil, static), oracle.run(soil, static)
    errors = {name: check_values(f"source lateral {name}", value[mask], np.asarray(expected[name])[mask]) for name, value in actual.items()}
    # Account for the receiver rule's actual domain gain/loss. Receiver weights
    # are not assumed to be a conservative donor partition.
    out = np.where(mask, expected["lateral_flow"], 0.)
    accepted_fraction = np.zeros(mask.shape, dtype=np.float64)
    for weight, (dr, dc) in zip(np.asarray(static.routing_weights), _DIRECTIONS):
        # Transpose the receiver rule to count each donor's accepted fractions.
        receiver = np.where(mask, weight, 0.)
        accepted_fraction += np.nan_to_num(_neighbor(receiver, -dr, -dc), nan=0.)
    delta = np.sum(actual["storage_amt"][mask] - np.asarray(soil.total_water_cm)[mask])
    gain_loss = np.sum((out * (accepted_fraction-1.))[mask])
    np.testing.assert_allclose(delta, gain_loss, rtol=512*np.finfo(float).eps,
                               atol=1e-10, err_msg="receiver-rule domain ledger")
    return {"max_error_cm": max(errors.values()), "domain_change_cell_cm": float(delta),
            "domain_ledger_error_cell_cm": float(abs(delta-gain_loss))}


def check_gather():
    cases, maximum = 0, 0.
    # Differing receiver weights distinguish receiver and donor conventions;
    # finite inactive donors verify the GRASS MASK, beyond NaN-only fixtures.
    for shape in ((1, 1), (1, 4), (4, 1), (4, 5)):
        flow = np.arange(1, np.prod(shape)+1, dtype=np.float64).reshape(shape)
        for direction in range(8):
            for holes in (False, True):
                mask = np.ones(shape, dtype=bool)
                if holes:
                    mask.flat[0] = False
                weights = np.zeros((8,)+shape, dtype=np.float64)
                weights[direction] = np.linspace(.1, 1.2, np.prod(shape)).reshape(shape)
                expected = np.zeros(shape, dtype=np.float64)
                dr, dc = _DIRECTIONS[direction]
                for row, col in np.argwhere(mask):
                    r, c = row+dr, col+dc
                    if 0 <= r < shape[0] and 0 <= c < shape[1] and mask[r, c] and np.isfinite(flow[r, c]):
                        expected[row, col] = flow[r, c] * weights[direction, row, col]
                actual = _run_gather(jnp.asarray(flow), jnp.asarray(weights), jnp.asarray(mask))
                maximum = max(maximum, check_values("orientation/MASK/border gather", np.asarray(actual), expected))
                cases += 1
    return {"cases": cases, "max_error_cm": maximum}


def synthetic_lateral(template, oracle):
    shape = (5, 6)
    f = lambda value: jnp.full(shape, value, dtype=jnp.float64)
    updates = {field.name: f(0.) for field in fields(template) if np.ndim(getattr(template, field.name)) == 2}
    mask = np.ones(shape, dtype=bool)
    mask[0, 0] = False
    updates.update(active_mask=jnp.asarray(mask), slope_percent=jnp.asarray(np.arange(30.).reshape(shape)),
        flowunits=jnp.asarray((np.arange(30).reshape(shape) % 4 != 0).astype(np.float64)),
        soil_depth_a_cm=f(30.), soil_depth_b_cm=f(70.), fieldcap_a_cm=f(8.4), fieldcap_b_cm=f(19.6),
        sat_mc_a=f(.45), sat_mc_b=f(.45), ksat_matrix_a=f(50.), ksat_matrix_b=f(20.),
        ksat_mpore_a=f(1e7), ksat_mpore_b=f(1e7), kfc_a=f(5.), kfc_b=f(2.))
    weights = np.arange(8*30, dtype=np.float64).reshape((8,)+shape) / 400.
    static = replace(template, **updates, routing_weights=jnp.asarray(weights), grid=GridSpec(*shape, 30., 0., 0.))
    static = replace(static, **{key: getattr(static, key).astype(jnp.float32).astype(jnp.float64) for key in _SOURCE_FIELDS.values()})
    fc_a, fc_b = float(static.fieldcap_a_cm[0, 0]), float(static.fieldcap_b_cm[0, 0])
    a = np.resize(np.array([0., fc_a-1e-8, fc_a, fc_a+1e-8, 13.5-1e-8, 13.5, 13.5+1e-8]), shape)
    b = np.resize(np.array([0., fc_b-1e-8, fc_b, fc_b+1e-8, 31.5-1e-8, 31.5, 31.5+1e-8]), shape)
    soil = TwoLayerState(jnp.asarray(a+b), jnp.asarray(np.stack((a, b))), static.active_mask)
    return check_lateral(oracle, soil, static)


def load_capture(path, source_path, reference_dir):
    path = Path(path)
    metadata = json.loads(path.with_suffix(".json").read_text(encoding="utf-8"))
    binary = path.with_suffix(".npz")
    if metadata.get("schema_version") != 1 or binary.stat().st_size != metadata["binary_bytes"] or sha256(binary.read_bytes()).hexdigest() != metadata["binary_sha256"]:
        raise ValueError("corrupt GRASS capture inventory")
    if sha256(Path(source_path).read_bytes()).hexdigest() != metadata["source_sha256"]:
        raise ValueError("source differs from GRASS capture provenance")
    capture_script = Path(source_path).resolve().parents[1]/"tools/capture_jax_routing_reference.py"
    if sha256(capture_script.read_bytes()).hexdigest() != metadata["capture_script_sha256"]:
        raise ValueError("capture script differs from GRASS provenance")
    manifest = json.loads((Path(reference_dir)/"manifest.json").read_text(encoding="utf-8-sig"))
    if {item["path"]: item["sha256"] for item in manifest["artifacts"]} != metadata["frozen_artifacts"]:
        raise ValueError("capture/frozen artifact provenance differs")
    with np.load(binary, allow_pickle=False) as arrays:
        return metadata, {name: arrays[name] for name in arrays.files}


def check_frozen_basin(reference_dir, static):
    # Feed the actual Perl-reported aggregates into the isolated recurrence.
    # Their printed precision is part of this reporting-coupled legacy path.
    rows = [line.split() for line in (Path(reference_dir)/"daily_mass_balance_2015-01.csv").read_text().splitlines() if line.strip()]
    if len(rows) != 30 or any(len(row) != 24 or row[0] != "1" or row[1] != f"2015-01-{i+1:02d}" for i, row in enumerate(rows)):
        raise ValueError("expected all 30 frozen full-report rows for watershed 1")
    initial = initialize_legacy_state(static).basin
    check_values("initial reservoir", np.asarray(initial.reservoir_volume), np.asarray(static.basin_reservoir_initial))
    check_values("initial baseflow", np.asarray(initial.baseflow_cms), np.zeros_like(np.asarray(static.basin_ids)))
    point = replace(static, grid=GridSpec(1, 1, 30., 0., 0.), active_mask=jnp.ones((1, 1), dtype=bool),
                    watershed=jnp.ones((1, 1), dtype=jnp.float64))
    state, max_q, max_base, max_volume = initial, 0., 0., 0.
    for row in rows:
        state, discharge = _run_basin(state, point, jnp.asarray([[float(row[3])]]), jnp.asarray([[float(row[12])]]))
        expected_q, expected_base = np.asarray([float(row[13])]), np.asarray([float(row[17])])
        max_q = max(max_q, check_values("frozen discharge", np.asarray(discharge), expected_q, atol=1e-15))
        max_base = max(max_base, check_values("frozen baseflow", np.asarray(state.baseflow_cms), expected_base, atol=1e-15))
        max_volume = max(max_volume, check_values("frozen reservoir", np.asarray(state.reservoir_volume),
                              expected_base/np.asarray(point.basin_reservoir_coeff)))
    return {"days": len(rows), "max_discharge_error_cms": max_q, "max_baseflow_error_cms": max_base,
            "max_reservoir_error": max_volume}


def verification_report(reference_dir, source_path, capture_path):
    reference_dir = Path(reference_dir)
    count = verify_reference_manifest(reference_dir)
    metadata, capture = load_capture(capture_path, source_path, reference_dir)
    static = load_legacy_bundle(reference_dir/"model_inputs.zip")
    reservoir = np.loadtxt(StringIO(metadata["reservoir_ini"]), dtype=np.float64, ndmin=2)
    if reservoir.shape != (1, 4) or reservoir[0, 0] != 1.:
        raise ValueError("this frozen reporting gate requires the captured single watershed")
    for i, name in enumerate(("basin_ids", "basin_area_cells", "basin_reservoir_initial", "basin_reservoir_coeff")):
        check_values(f"captured reservoir configuration {name}", np.asarray(getattr(static, name)), reservoir[:, i], atol=0.)
    mask = np.asarray(static.active_mask)
    if capture["support"].dtype != np.bool_ or not np.array_equal(mask, capture["support"]) or int(mask.sum()) != count or metadata["required_cells"] != count:
        raise ValueError("GRASS/model required support differs")
    if len(capture["sample_flat"]) != metadata["representative_cells"] or not mask.ravel()[capture["sample_flat"]].all():
        raise ValueError("invalid required GRASS representative cells")
    for actual, expected in (("slope_percent", "slope"), ("root_zone_cm", "root_zone"),
                             ("wiltpt_amt_cm", "wiltpt_amt"), ("et_reduction_mc", "ETreduction_mc")):
        observed, target = np.asarray(getattr(static, actual))[mask], capture[expected][mask]
        check_values(f"captured GRASS {expected}", observed, target, atol=0.)
        np.testing.assert_array_equal(observed, target, err_msg=f"captured GRASS {expected}")
    gather = check_gather()
    oracle = RoutingOracle(source_path)
    synthetic = synthetic_lateral(static, oracle)
    trace = FrozenTraceProvider(reference_dir/"two_day_internal_trace.zip", static.grid)
    days = []
    for date in ("2015-01-01", "2015-01-02"):
        if not np.array_equal(mask, np.isfinite(trace.raster(date, "saturation"))):
            raise ValueError("frozen lateral support differs")
        soil = TwoLayerState(jnp.asarray(trace.raster(date, "storage_amt")),
            jnp.asarray(np.stack((trace.raster(date, "storage_amt_A"), trace.raster(date, "storage_amt_B")))), static.active_mask)
        report = check_lateral(oracle, soil, static)
        actual = _lateral_values(soil, static)
        capture_error = 0.
        for name, value in actual.items():
            capture_error = max(capture_error, check_values(f"captured GRASS {date} {name}",
                value.ravel()[capture["sample_flat"]], capture[f"{date}_{name}"]))
        report["max_grass_error_cm"] = capture_error
        days.append(report)
    basin = check_frozen_basin(reference_dir, static)
    return {"required_cells": count, "grass_representatives": metadata["representative_cells"],
            "gather": gather, "synthetic": synthetic, "lateral_days": days, "basin": basin}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, default=Path("reference/mica-creek-jan2015"))
    parser.add_argument("--source", type=Path, default=Path("R/smr_buildup.pl"))
    parser.add_argument("--capture", type=Path, default=Path("tests/fixtures/jax_verify03_grass"))
    args = parser.parse_args(argv)
    try:
        report = verification_report(args.reference_dir, args.source, args.capture)
    except (AssertionError, ValueError, OSError, KeyError, RuntimeError, IndexError, SyntaxError) as exc:
        print(f"FAIL: VERIFY-03: {exc}", file=sys.stderr)
        return 1
    print("PASS: VERIFY-03 " + json.dumps(report, sort_keys=True))
    print("Legacy receiver rule and asymmetric reservoir units retained; whole-engine parity remains VERIFY-04.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
