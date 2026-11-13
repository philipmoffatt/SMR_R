"""VERIFY-02: fixed six-hour soil stages against the original Perl expressions.

This deliberately excludes lateral transfers and basin routing. The oracle reads
the source, evaluates scalar map algebra, and never calls the JAX soil routines.
Frozen daily states supply representative starting conditions, not substep truth.
"""

import argparse
import ast
from dataclasses import fields, replace
import operator
from pathlib import Path
import re
import sys
from zipfile import ZipFile

import jax.numpy as jnp
import numpy as np

from .fluxes.legacy_math import compile_legacy_cpu
from .fluxes.soil import _reallocate, allocate_input, evapotranspiration, percolation_and_runoff
from .fluxes.surface import road_runoff
from .io import FrozenTraceProvider
from .legacy import _landcover_properties, _read_esri_ascii, load_legacy_bundle
from .state import GridSpec, TwoLayerState
from .verify_surface import RELATIVE_TOLERANCE, WATER_ATOL_CM, verify_reference_manifest


_BINARY = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul, ast.Div: operator.truediv}
_COMPARE = {ast.Lt: operator.lt, ast.Gt: operator.gt, ast.LtE: operator.le, ast.GtE: operator.ge, ast.Eq: operator.eq}
_STATIC = {
    "landuse": "landuse", "soil_depth_a_cm": "soil_depth_A", "soil_depth_b_cm": "soil_depth_B",
    "sat_amt_cm": "sat_amt", "fieldcap_amt_cm": "fieldcap_amt", "fieldcap_a_cm": "fieldcap_amt_A",
    "fieldcap_b_cm": "fieldcap_amt_B", "sat_mc_a": "sat_mc_A", "sat_mc_b": "sat_mc_B",
    "wiltpt_mc_a": "wiltpt_mc_A", "ksubsurface": "Ksubsurface",
}
_RAW = set(_STATIC.values()) | {"soil_depth", "wiltpt_mc_B", "wiltpt_amt"}


def _evaluate(node, env):
    """Small allowlisted evaluator; scalar conditionals evaluate one arm only."""
    if isinstance(node, ast.Constant) and type(node.value) in (int, float):
        return node.value
    if isinstance(node, ast.Name):
        return env[node.id]
    if isinstance(node, ast.BinOp) and type(node.op) in _BINARY:
        return _BINARY[type(node.op)](_evaluate(node.left, env), _evaluate(node.right, env))
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        value = _evaluate(node.operand, env)
        return -value if isinstance(node.op, ast.USub) else value
    if isinstance(node, ast.Compare) and len(node.ops) == 1 and type(node.ops[0]) in _COMPARE:
        return _COMPARE[type(node.ops[0])](_evaluate(node.left, env), _evaluate(node.comparators[0], env))
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and not node.keywords:
        if node.func.id == "choose" and len(node.args) == 3:
            condition = _evaluate(node.args[0], env)
            if np.ndim(condition) == 0:
                return _evaluate(node.args[1] if condition else node.args[2], env)
            return np.where(condition, _evaluate(node.args[1], env), _evaluate(node.args[2], env))
        if node.func.id in ("min", "max") and len(node.args) == 2:
            function = np.minimum if node.func.id == "min" else np.maximum
            return function(_evaluate(node.args[0], env), _evaluate(node.args[1], env))
    raise ValueError(f"unsupported source expression: {ast.dump(node)}")


class SourceOracle:
    """Extract only the five bounded source blocks owned by this gate."""

    def __init__(self, path):
        source = Path(path).read_text(encoding="utf-8")
        for name, expected in (("gridsize", 30.), ("temp_time_step", 6.)):
            values = re.findall(r"my \$" + name + r"\s*=\s*([0-9.]+)\s*;", source)
            if len(values) != 1 or float(values[0]) != expected:
                raise ValueError(f"unsupported source {name}; this gate requires {expected}")
        self.blocks = {}
        bounds = {
            "properties": ("# Land cover maps", "# Initial conditions", 7),
            "road": ("# Road runoff", "# Daily accumulators", 4),
            "allocation": ("# Water distribution to soil layers", "# Lateral flow", 5),
            "redistribution": ("# Re-distribute between A and B", "# ET (", 3),
            "et": ("# ET (", "# Percolation", 7),
            "drainage": ("# Percolation", 'expression=runoff_$year = runoff_flow', 11),
        }
        for name, (start, end, count) in bounds.items():
            if source.count(start) != 1 or source.count(end) != 1:
                raise ValueError(f"ambiguous/missing source block: {name}")
            block = source.split(start, 1)[1].split(end, 1)[0]
            expressions = re.findall(r'expression=([^"\r\n]+)', block)
            if len(expressions) != count:
                raise ValueError(f"source block {name}: expected {count} assignments, got {len(expressions)}")
            parsed = []
            for expression in expressions:
                target, value = expression.split(" = ", 1)
                value = re.sub(r"\$(\w+)", r"\1", value)
                value = re.sub(r"\bif\(", "choose(", value)
                parsed.append((target.strip(), ast.parse(value, mode="eval").body))
            self.blocks[name] = parsed

    def run(self, block, env):
        for target, expression in self.blocks[block]:
            env[target] = _evaluate(expression, env)
        return env


def check_values(name, actual, expected, *, atol=WATER_ATOL_CM):
    actual, expected = np.asarray(actual), np.asarray(expected)
    if actual.shape != expected.shape or actual.size == 0:
        raise ValueError(f"{name}: missing values or different shapes")
    if actual.dtype != np.float64 or expected.dtype != np.float64:
        raise ValueError(f"{name}: float64 values required")
    if not np.isfinite(actual).all() or not np.isfinite(expected).all():
        raise ValueError(f"{name}: nonfinite required value")
    np.testing.assert_allclose(actual, expected, rtol=RELATIVE_TOLERANCE, atol=atol, err_msg=name)
    return float(np.max(np.abs(actual - expected)))


def synthetic_cases():
    """A branch matrix including all six land uses and zero-root/PET cases."""
    cases = []
    for landuse in (0, 1, 2, 3, 4, 5, 6, 41, 43):
        base = dict(landuse=float(landuse), soil_depth=100., soil_depth_A=30., soil_depth_B=70.,
                    sat_amt=45., sat_mc_A=.45, sat_mc_B=.45, fieldcap_amt=28.,
                    fieldcap_amt_A=8.4, fieldcap_amt_B=19.6, wiltpt_mc_A=.1, wiltpt_mc_B=.1,
                    wiltpt_amt=10., temp_time_step=6., gridsize=30.,
                    temp=1., temp_sum=4., snowmelt=0., throughfall=0.,
                    cc_water=.65, cc_urban=.3, cc_forest=1., cc_shrub=.45, cc_grass=.3, cc_row_crop=1.15)
        # Include each side of the allocation/FC/saturation boundaries, and the
        # forest A-only root storage versus other classes' total storage.
        for total in (0., 3.-1e-8, 3., 3.+1e-8, 10.-1e-8, 10., 10.+1e-8,
                      8.4-1e-8, 8.4, 8.4+1e-8, 41.5-1e-8, 41.5, 41.5+1e-8,
                      19.6-1e-8, 19.6, 19.6+1e-8, 22.4-1e-8, 22.4, 22.4+1e-8,
                      28.-1e-8, 28., 28.+1e-8, 39.9-1e-8, 39.9, 39.9+1e-8,
                      45.-1e-8, 45., 45.+1e-8):
            for water, pet, ksub in ((0., 0., 0.), (.04, .8, .5), (60., 8., 1000.)):
                env = dict(base, storage_amt=total, storage_amt_A=min(total, 8.4),
                           storage_amt_B=total-min(total, 8.4), throughfall=water*4,
                           pet=pet, Ksubsurface=ksub)
                cases.append(env)
        # Road threshold, both source types, and the dry-source division guard.
        for amount in (0., 1.5-1e-8, 1.5, 1.5+1e-8, 12.):
            for melt_fraction in (0., .4, 1.):
                cases.append(dict(base, storage_amt=20., storage_amt_A=6., storage_amt_B=14.,
                                  snowmelt=amount*melt_fraction, throughfall=amount*(1-melt_fraction),
                                  pet=.8, Ksubsurface=.5))
        for temperature, temperature_sum in ((-1., 0.), (0., 0.), (2., 3.)):
            cases.append(dict(base, storage_amt=20., storage_amt_A=6., storage_amt_B=14.,
                              snowmelt=1., throughfall=.4, pet=.8, Ksubsurface=.5,
                              temp=temperature, temp_sum=temperature_sum))
    return cases


def _packed_static(template, cases):
    shape = (1, len(cases))
    column = lambda key: np.asarray([[c[key] for c in cases]], dtype=np.float64)
    props = _landcover_properties(*(column(k) for k in ("landuse", "soil_depth", "soil_depth_A", "soil_depth_B",
                                  "fieldcap_amt", "fieldcap_amt_A", "wiltpt_mc_A", "wiltpt_mc_B")))
    # Unused raster fields have the same shape, avoiding accidental broadcasting.
    updates = {f.name: jnp.zeros(shape, dtype=jnp.float64) for f in fields(template)
               if f.name != "active_mask" and np.ndim(getattr(template, f.name)) == 2}
    updates.update({name: jnp.asarray(column(key)) for name, key in _STATIC.items()})
    updates.update(zip(("max_canopy_cm", "canopy_cover", "root_zone_cm", "et_reduction_mc",
                        "tmax_rain_c", "tmin_snow_c", "wiltpt_amt_cm"), map(jnp.asarray, props)))
    return replace(template, **updates, active_mask=jnp.ones(shape, dtype=bool),
                   grid=GridSpec(1, len(cases), 30., 0., 0.))


def _jax_stages(soil, static, forcing, snowmelt, throughfall, pet, temp, temp_sum):
    melt, rain, road = road_runoff(static, snowmelt, throughfall)
    water = melt * jnp.maximum(temp, 0.) / jnp.where(temp_sum == 0., 1., temp_sum) + rain * .25
    allocated = allocate_input(soil, static, water)
    a, b = _reallocate(allocated.total_water_cm, allocated.layer_water_cm[0], static)
    redistributed = TwoLayerState(allocated.total_water_cm, jnp.stack((a, b)), soil.active_mask)
    evaporated, et = evapotranspiration(redistributed, static, forcing, pet)
    drained, perc, runoff, saturation = percolation_and_runoff(evaporated, static)
    return melt, rain, road, water, allocated, redistributed, evaporated, et, drained, perc, runoff, saturation


_run_stages = compile_legacy_cpu(_jax_stages)


def verify_cases(oracle, template, forcing, cases):
    """Compare every stage and independently reconstruct each water ledger."""
    if not cases:
        raise ValueError("no required soil cases")
    column = lambda key: np.asarray([[c[key] for c in cases]], dtype=np.float64)
    static = _packed_static(template, cases)
    # The forcing coefficient vector is deliberately nonuniform in the branch matrix.
    coeff = [cases[0][k] for k in ("cc_water", "cc_urban", "cc_forest", "cc_shrub", "cc_grass", "cc_row_crop")]
    if any(any(case[key] != coeff[i] for i, key in enumerate(
               ("cc_water", "cc_urban", "cc_forest", "cc_shrub", "cc_grass", "cc_row_crop"))) for case in cases):
        raise ValueError("soil cases must share one land-use coefficient vector")
    forcing = replace(forcing, et_coeff_by_landuse=jnp.asarray(coeff, dtype=jnp.float64))
    soil = TwoLayerState(jnp.asarray(column("storage_amt")),
                         jnp.asarray(np.stack((column("storage_amt_A"), column("storage_amt_B")))), static.active_mask)
    actual = _run_stages(soil, static, forcing, *(jnp.asarray(column(k)) for k in
                          ("snowmelt", "throughfall", "pet", "temp", "temp_sum")))
    expected = []
    for case in cases:
        env = dict(case, actualET_daily_flow=0., perc_daily_flow=0., runoff_daily_flow=0.)
        oracle.run("properties", env)
        properties = [env[k] for k in ("root_zone", "wiltpt_amt", "ETreduction_mc")]
        oracle.run("road", env)
        road = [env[k] for k in ("snowmelt", "throughfall", "road_runoff")]
        oracle.run("allocation", env)
        state = lambda: [env[k] for k in ("storage_amt", "storage_amt_A", "storage_amt_B")]
        allocation, water = state(), env["input"]
        oracle.run("redistribution", env)
        redistribution = state()
        oracle.run("et", env)
        evaporation, et = state(), env["actualET_flow"]
        oracle.run("drainage", env)
        expected.append(properties + road + [water] + allocation + redistribution + evaporation +
                        [et] + state() + [env[k] for k in ("perc", "runoff_flow", "saturation")])
    expected = np.asarray(expected, dtype=np.float64).T[:, None, :]
    melt, rain, road, water, allocated, redistributed, evaporated, et, drained, perc, runoff, saturation = actual
    flatten_state = lambda s: [s.total_water_cm, *s.layer_water_cm]
    observed = [static.root_zone_cm, static.wiltpt_amt_cm, static.et_reduction_mc, melt, rain, road, water,
                *flatten_state(allocated), *flatten_state(redistributed), *flatten_state(evaporated), et,
                *flatten_state(drained), perc, runoff, saturation]
    names = ["root zone", "wilting storage", "ET reduction", "retained melt", "retained rain", "road runoff", "soil input"]
    names += [f"{stage} {part}" for stage in ("allocation", "redistribution", "ET") for part in ("total", "A", "B")]
    names += ["actual ET", "drainage total", "drainage A", "drainage B", "percolation", "saturation runoff", "saturation percent"]
    errors = {name: check_values(name, np.asarray(value), expected[i],
                                atol=1e-10 if name == "saturation percent" else WATER_ATOL_CM)
              for i, (name, value) in enumerate(zip(names, observed))}
    initial = np.asarray(soil.total_water_cm)
    zeros = np.zeros_like(initial)
    ledgers = {
        "road": column("snowmelt") + column("throughfall") - np.asarray(melt + rain + road),
        "allocation": initial + np.asarray(water) - np.asarray(allocated.total_water_cm),
        "ET": np.asarray(redistributed.total_water_cm) - np.asarray(evaporated.total_water_cm + et),
        "drainage": np.asarray(evaporated.total_water_cm) - np.asarray(drained.total_water_cm + perc + runoff),
        "column": initial + np.asarray(water) - np.asarray(drained.total_water_cm + et + perc + runoff),
    }
    for name, residual in ledgers.items():
        errors[f"ledger {name}"] = check_values(f"ledger {name}", residual, zeros)
    # The frozen sat_amt raster is not exactly the sum of the independently
    # rounded layer capacities. Above saturation, Perl clips layers and total
    # against those separate quantities. Check that inherited representation
    # gap against the source rather than hiding it in a relaxed water budget.
    layer_gap = np.asarray(drained.layer_water_cm).sum(axis=0) - np.asarray(drained.total_water_cm)
    expected_gap = expected[18] + expected[19] - expected[17]
    errors["source layer representation"] = check_values("source layer representation", layer_gap, expected_gap)
    return {"cases": len(cases), "max_error": max(errors.values()), "errors": errors,
            "max_layer_representation_gap_cm": float(np.max(np.abs(layer_gap)))}


def verification_report(reference_dir, source_path):
    reference_dir = Path(reference_dir)
    required = verify_reference_manifest(reference_dir)
    oracle = SourceOracle(source_path)
    static = load_legacy_bundle(reference_dir / "model_inputs.zip")
    mask = np.asarray(static.active_mask)
    trace = FrozenTraceProvider(reference_dir / "two_day_internal_trace.zip", static.grid)
    for date in ("2015-01-01", "2015-01-02"):
        if not np.array_equal(mask, np.isfinite(trace.raster(date, "saturation"))) or int(mask.sum()) != required:
            raise ValueError("frozen soil support differs from the full required support")
    with ZipFile(reference_dir / "model_inputs.zip") as bundle:
        raw = {key: _read_esri_ascii(bundle, f"{key}.asc")[1] for key in _RAW}
    props = oracle.run("properties", dict(raw))
    for key, name in (("root_zone", "root_zone_cm"), ("wiltpt_amt", "wiltpt_amt_cm"), ("ETreduction_mc", "et_reduction_mc")):
        check_values(f"all frozen {key}", np.asarray(getattr(static, name))[mask], props[key][mask])
    # Read the real weather coefficient order without using the model's loader.
    weather = [line.split() for line in (reference_dir / "gridmet_weather_2015-01.csv").read_text().splitlines() if line.strip()]
    from .io import ReplayRadiationProvider, load_daily_forcings
    import jax
    forcing = jax.tree_util.tree_map(lambda x: x[0], load_daily_forcings(
        reference_dir / "gridmet_weather_2015-01.csv",
        ReplayRadiationProvider(reference_dir / "two_day_internal_trace.zip", static.grid), 1))
    synthetic = verify_cases(oracle, static, forcing, synthetic_cases())
    representatives = []
    classes = sorted(int(x) for x in np.unique(raw["landuse"][mask]))
    for lu in classes:
        cells = np.argwhere(mask & (raw["landuse"] == lu))
        representatives.extend(tuple(cells[i]) for i in np.unique(np.linspace(0, len(cells)-1, min(16, len(cells)), dtype=int)))
    cases = []
    for day, date in enumerate(("2015-01-01", "2015-01-02")):
        names = ("storage_amt", "storage_amt_A", "storage_amt_B", "snowmelt", "throughfall", "canopy_evap")
        snapshots = {key: trace.raster(date, key) for key in names}
        for key, value in snapshots.items():
            check_values(f"frozen required {date} {key}", value[mask], value[mask])
        for row, col in representatives:
            base = {key: float(value[row, col]) for key, value in raw.items()}
            base.update({key: float(value[row, col]) for key, value in snapshots.items()})
            base.update(temp_time_step=6., gridsize=30., temp=1., temp_sum=4.,
                        pet=max(0., float(weather[day][10])-base.pop("canopy_evap")))
            base.update(zip(("cc_water", "cc_urban", "cc_forest", "cc_shrub", "cc_grass", "cc_row_crop"), map(float, weather[day][17:23])))
            # Actual state/forcing plus controlled dry, wet and strong-ET probes.
            for throughfall, pet in ((base["throughfall"], base["pet"]), (0., 0.), (4., .8)):
                cases.append(dict(base, throughfall=throughfall, pet=pet))
    representative = verify_cases(oracle, static, forcing, cases)
    return {"required_cells": required, "landuse_classes": classes, "representative_cells": len(representatives),
            "synthetic": synthetic, "representative": representative}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, default=Path("reference/mica-creek-jan2015"))
    parser.add_argument("--source", type=Path, default=Path("R/smr_buildup.pl"))
    args = parser.parse_args(argv)
    try:
        report = verification_report(args.reference_dir, args.source)
    except (AssertionError, ValueError, OSError, KeyError, RuntimeError, ZeroDivisionError, SyntaxError) as exc:
        print(f"FAIL: VERIFY-02: {exc}", file=sys.stderr)
        return 1
    print(f"PASS: VERIFY-02; properties on {report['required_cells']} frozen cells; classes {report['landuse_classes']}")
    for group in ("synthetic", "representative"):
        result = report[group]
        print(f"{group}: {result['cases']} cases; max comparison/ledger error {result['max_error']:.17g}")
        print(f"  inherited layer/total representation gap: {result['max_layer_representation_gap_cm']:.17g} cm")
    print("Fixed six-hour column stages only; lateral routing, integrated parity and adaptive stepping remain separate.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
