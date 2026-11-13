"""Float64 CPU/CUDA verification and synchronized prepared-workflow benchmarks.

No physics functions are replaced. Strict comparisons use the existing unit
budgets; audit mode retains failures and a nonzero exit. Timing ratios are
qualified only when the corresponding complete-output comparison passes.
"""

import argparse
from dataclasses import replace
from hashlib import sha256
import importlib.metadata
import json
from pathlib import Path
import platform
import statistics
import time

import jax
import jax.numpy as jnp
import numpy as np

from .contracts import StaticInputs, DailyForcing, RadiationComponents
from .engine import run_days
from .initialization import initialize_legacy_state
from .io import load_daily_forcings
from .legacy import load_legacy_bundle, _landcover_properties
from .plant import (VegetationDemand, SurfaceWater, external_day_kernel,
                    validate_external_inputs, verify_water_balance)
from .radiation import PreparedRadiationProvider
from .state import GridSpec, TwoLayerState


RTOL = 512*np.finfo(np.float64).eps
# FCELL round trips are part of the source model, not redundant conversions.
# XLA's excess-precision optimization can remove them on CUDA, changing snow
# branches and subsequent fluxes despite float64 inputs and outputs.
CPU_OPTIONS = {"xla_cpu_enable_fast_math": False, "xla_cpu_enable_fast_min_max": False,
               "xla_allow_excess_precision": False}
GPU_OPTIONS = {"xla_gpu_enable_fast_min_max": False, "xla_allow_excess_precision": False}
EXACT_FIELDS = {"age_days", "snow_age_days", "annual_saturation_count"}
TEMP_FIELDS = {"tavg_c", "tdew_c", "surface_temp_c", "snow_surface_temp_c", "saturation_percent"}


def field_budget(name):
    if name in EXACT_FIELDS:
        return 0., 0., "exact"
    if name.startswith("q_") or "energy" in name:
        return 1e-9, RTOL, "energy"
    if name in TEMP_FIELDS:
        return 1e-10, RTOL, "temperature_or_percent"
    if name.endswith("_cms"):
        return 1e-15, RTOL, "cms"
    if name in ("reservoir_volume", "basin_reservoir"):
        return 1e-12, RTOL, "reservoir"
    if name == "domain_residual_m3":
        return 1e-10, RTOL, "m3"
    if name == "albedo":
        return 1e-12, RTOL, "dimensionless"
    if name.endswith("_cm"):
        return 1e-12, RTOL, "cm"
    raise ValueError(f"no declared comparison budget for {name}")


def compare_outputs(expected, actual, required_mask):
    """Compare complete pytrees against independent CPU support, never intersection."""
    expected_leaves, expected_tree = jax.tree_util.tree_flatten_with_path(expected)
    actual_leaves, actual_tree = jax.tree_util.tree_flatten_with_path(actual)
    if expected_tree != actual_tree:
        raise ValueError("CPU/GPU output structures differ")
    required_mask = np.asarray(required_mask)
    if required_mask.dtype != np.bool_ or required_mask.ndim != 2 or not required_mask.any():
        raise ValueError("need nonempty independent boolean model support")
    failures, errors, units = [], {}, {}
    for (path, cpu), (_, gpu) in zip(expected_leaves, actual_leaves):
        label = jax.tree_util.keystr(path)
        name = getattr(path[-1], "key", getattr(path[-1], "name", None))
        cpu, gpu = np.asarray(cpu), np.asarray(gpu)
        if cpu.shape != gpu.shape or cpu.dtype != gpu.dtype:
            raise ValueError(f"{label}: CPU/GPU shape or type differs")
        if cpu.dtype.kind in "biu":
            if not np.array_equal(cpu, gpu):
                failures.append({"field": label, "reason": "exact boolean/integer mismatch"})
            continue
        if cpu.dtype != np.float64:
            raise ValueError(f"{label}: floating outputs must be float64")
        required = (np.broadcast_to(required_mask, cpu.shape) if cpu.ndim >= 2 and cpu.shape[-2:] == required_mask.shape
                    else np.ones(cpu.shape, dtype=bool))
        if not np.isfinite(cpu[required]).all() or not np.isfinite(gpu[required]).all():
            raise ValueError(f"{label}: nonfinite required CPU/GPU output")
        support = np.isfinite(cpu)
        if not np.array_equal(support, np.isfinite(gpu)):
            raise ValueError(f"{label}: CPU/GPU finite footprints differ")
        atol, rtol, unit = field_budget(name)
        delta = np.abs(gpu[support]-cpu[support])
        maximum = float(delta.max()) if delta.size else 0.
        errors[label] = maximum
        units[unit] = max(units.get(unit, 0.), maximum)
        bad = delta > atol+rtol*np.abs(cpu[support])
        if bad.any():
            failures.append({"field": label, "max_abs": maximum, "outside_budget": int(bad.sum()),
                             "atol": atol, "rtol": rtol})
    presence = 0
    if isinstance(expected, tuple) and len(expected) == 2 and isinstance(expected[1], dict) and "swe_cm" in expected[1]:
        cpu, gpu = np.asarray(expected[1]["swe_cm"]), np.asarray(actual[1]["swe_cm"])
        presence = int(np.count_nonzero((cpu[..., required_mask] > 0.) != (gpu[..., required_mask] > 0.)))
    return {"passed": not failures, "failed_fields": len(failures), "failures": failures,
            "max_errors_by_unit": units, "max_field_errors": errors,
            "snow_presence_disagreements_cell_days": presence}


def require_devices():
    if not jax.config.jax_enable_x64:
        raise ValueError("GPU verification requires float64 enabled before initialization")
    gpu, cpu = jax.devices("gpu"), jax.devices("cpu")
    if not gpu or gpu[0].platform != "gpu":
        raise RuntimeError("CUDA GPU required; CPU fallback is not GPU verification")
    return cpu[0], gpu[0]


def verify_snapshot(root, metadata_path):
    root = Path(root).resolve()
    metadata = json.loads(Path(metadata_path).read_text(encoding="utf-8"))
    if metadata.get("schema_version") != 1 or not metadata.get("files"):
        raise ValueError("incomplete source snapshot witness")
    for relative, expected in metadata["files"].items():
        path = (root/relative).resolve()
        if not path.is_relative_to(root) or sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"source snapshot checksum differs: {relative}")
    return metadata["source_revision"]


def _host_tree(value):
    return jax.tree_util.tree_map(lambda v: np.asarray(v), jax.device_get(value))


def _put(value, device):
    placed = jax.device_put(value, device)
    jax.block_until_ready(placed)
    for leaf in jax.tree_util.tree_leaves(placed):
        if isinstance(leaf, jax.Array) and leaf.devices() != {device}:
            raise ValueError("unexpected device placement")
    return placed


def _memory(device):
    values = device.memory_stats()
    return None if values is None else {k: int(v) for k, v in values.items()}


def execute(function, host_args, device):
    """First execution plus a resident executable for subsequent qualified timing."""
    started = time.perf_counter()
    args = _put(host_args, device)
    transfer = time.perf_counter()-started
    flat_args, args_tree = jax.tree_util.tree_flatten(args)
    # JAX lower() reconstructs its public argument metadata using ArgInfo
    # objects, which are not arrays accepted by TwoLayerState's validator.
    # Stage array leaves and reconstruct the unchanged contract inside tracing.
    def flat_function(*leaves):
        return function(*jax.tree_util.tree_unflatten(args_tree, leaves))
    options = CPU_OPTIONS if device.platform == "cpu" else GPU_OPTIONS
    compiled_function = jax.jit(flat_function, compiler_options=dict(options))
    started = time.perf_counter()
    lowered = compiled_function.lower(*flat_args)
    lowering = time.perf_counter()-started
    started = time.perf_counter()
    executable = lowered.compile()
    compilation = time.perf_counter()-started
    started = time.perf_counter()
    result = executable(*flat_args)
    jax.block_until_ready(result)
    first = time.perf_counter()-started
    for leaf in jax.tree_util.tree_leaves(result):
        if isinstance(leaf, jax.Array) and leaf.devices() != {device}:
            raise ValueError("compiled output is on the wrong device")
    started = time.perf_counter()
    host_result = _host_tree(result)
    download = time.perf_counter()-started
    metrics = {"device": str(device), "device_kind": device.device_kind, "compiler_options": options,
        "input_transfer_s": transfer, "lowering_s": lowering, "compilation_s": compilation,
        "first_execution_s": first, "warmed_execution_s": [], "warmed_median_s": None,
        "output_download_s": download, "memory": _memory(device),
        "input_bytes": sum(v.nbytes for v in jax.tree_util.tree_leaves(host_args)),
        "output_bytes": sum(v.nbytes for v in jax.tree_util.tree_leaves(host_result))}
    return host_result, metrics, (executable, flat_args, args_tree)


def warmed_timing(resident, device, metrics, repeats, loader=None):
    executable, args, args_tree = resident
    times = []
    for _ in range(repeats):
        started = time.perf_counter()
        value = executable(*args)
        jax.block_until_ready(value)
        times.append(time.perf_counter()-started)
        del value
    metrics.update(warmed_execution_s=times, warmed_median_s=statistics.median(times))
    if loader is not None:
        # A separately measured complete warm workflow: fresh local input load,
        # initialization, transfer, cached executable, sync, and full download.
        started = time.perf_counter()
        fresh = loader()
        fresh_args = _put(fresh, device)
        fresh_leaves, fresh_tree = jax.tree_util.tree_flatten(fresh_args)
        if fresh_tree != args_tree:
            raise ValueError("warm workflow changed input contract")
        output = executable(*fresh_leaves)
        jax.block_until_ready(output)
        downloaded = _host_tree(output)
        metrics["measured_warm_prepared_workflow_s"] = time.perf_counter()-started
        del downloaded, output, fresh_args, fresh
    metrics["memory"] = _memory(device)


def synthetic_model():
    """Explicit small physical properties; all land covers, a mask hole and edges."""
    shape = (5, 7)
    field = lambda v: jnp.full(shape, v, dtype=jnp.float64)
    lu = (np.arange(35).reshape(shape) % 8)+1.
    lu[lu == 7.] = 41.
    lu[lu == 8.] = 43.
    props = _landcover_properties(lu, np.full(shape, 100.), np.full(shape, 30.), np.full(shape, 70.),
        np.full(shape, 28.), np.full(shape, 8.4), np.full(shape, .1), np.full(shape, .1))
    # Choose field-capacity conductivity from the stated matrix curve, rather
    # than injecting an arbitrary jump between its two source branches.
    kfc_a = 50.*np.exp((-13./.45)*(.45-.28))
    kfc_b = 20.*np.exp((-13./.45)*(.45-.28))
    mask = jnp.ones(shape, dtype=bool).at[2, 3].set(False)
    static = StaticInputs(GridSpec(*shape, 30., 0., 0.), mask, field(784.), field(15.),
        jnp.asarray(lu), field(1.), field(1.).at[0, 0].set(0.), jnp.full((8,)+shape, .125),
        field(45.), field(28.), field(8.4), field(19.6), jnp.asarray(props[-1]), field(30.), field(70.),
        field(.45), field(.45), field(.28), field(.28), field(.1), field(50.), field(20.), field(200.),
        field(80.), field(kfc_a), field(kfc_b), field(.5), *(jnp.asarray(v) for v in props[:-1]),
        jnp.asarray([1.]), jnp.asarray([34.]), jnp.asarray([0.]), jnp.asarray([.05]))
    conditions = [(-5., 0., 0.), (-5., 1., .05), (8., 1., .1), (-8., .2, 0.),
                  (3., 50., .5), (20., 0., 2.), (1., 0., 0.)]
    days = []
    for index, (temperature, precip, pet) in enumerate(conditions):
        scalar = lambda v: jnp.asarray(v, dtype=jnp.float64)
        components = RadiationComponents(field(1000. if temperature > 0 else 100.), field(100.), field(20.))
        days.append(DailyForcing(scalar(2015.), scalar(1.), scalar(index+1.), scalar(index+1.),
            scalar(temperature), scalar(temperature-2.), scalar(precip), scalar(pet),
            jnp.full((4,), temperature), scalar(.5), scalar(100.), scalar(10.), jnp.ones((6,)), None, components))
    forcing = jax.tree_util.tree_map(lambda *v: jnp.stack(v), *days)
    return initialize_legacy_state(static), forcing, static


def coupling_cases(static):
    """Host-validated conservation cases with explicitly consistent snow sources."""
    shape = static.grid.shape
    field = lambda v: jnp.full(shape, v, dtype=jnp.float64)
    static = replace(static, wiltpt_amt_cm=field(10.))
    for moisture, rain, demand, snow_melt in ((.1, 0., 100., False), (.44, 50., .2, False),
            (.2, .4, .02, False), (.2, 0., 0., False), (.2, .3, .02, True)):
        water = jnp.stack((field(moisture*30.), field(moisture*70.)))
        soil = TwoLayerState(water.sum(axis=0), water, static.active_mask)
        want = VegetationDemand(jnp.stack((field(demand), field(demand))),
            jnp.stack((field(demand*.2), field(0.))), field(.01 if demand else 0.))
        surface = SurfaceWater(field(rain+(.1 if snow_melt else 0.)), field(rain), field(.1 if snow_melt else 0.),
            field(.2 if snow_melt else 0.), field(1. if snow_melt else 0.), field(.91 if snow_melt else 0.),
            field(.01 if snow_melt else 0.))
        for boundary in ("legacy_drainage", "closed"):
            validate_external_inputs(soil, static, want, surface, boundary)
            yield (soil, static, want, surface), boundary


def verify_case(function, args, cpu, gpu, mask, *, repeats=0, loader=None):
    host_args = _host_tree(args)
    expected, cpu_timing, cpu_resident = execute(function, host_args, cpu)
    actual, gpu_timing, gpu_resident = execute(function, host_args, gpu)
    comparison = compare_outputs(expected, actual, mask)
    if comparison["passed"] and repeats:
        warmed_timing(cpu_resident, cpu, cpu_timing, repeats, loader)
        warmed_timing(gpu_resident, gpu, gpu_timing, repeats, loader)
        comparison["validated_compute_speedup"] = cpu_timing["warmed_median_s"]/gpu_timing["warmed_median_s"]
    return {"comparison": comparison, "cpu": cpu_timing, "gpu": gpu_timing}, expected, actual


def environment_report(cpu, gpu):
    versions = {d.metadata["Name"]: d.version for d in importlib.metadata.distributions()
                if d.metadata["Name"].lower().startswith(("jax", "numpy", "nvidia-", "tifffile"))}
    return {"system": platform.platform(), "python": platform.python_version(), "packages": versions,
            "cpu": str(cpu), "gpu": str(gpu), "gpu_kind": gpu.device_kind,
            "float64": jax.config.jax_enable_x64}


def verification_report(reference_dir, radiation_dir, *, days=(2, 30), repeats=5, audit_all=False):
    cpu, gpu = require_devices()
    report = {"environment": environment_report(cpu, gpu), "budgets": {"rtol": RTOL, "water_atol_cm": 1e-12,
        "energy_atol": 1e-9, "temperature_percent_atol": 1e-10, "cms_atol": 1e-15}, "synthetic": [], "cases": []}
    passed = True
    with jax.default_device(cpu):
        initial, forcing, static = synthetic_model()
    for policy in (False, True):
        function = lambda state, weather, constants: run_days(state, weather, constants, melt_small_snow_next_day=policy)
        result, _, _ = verify_case(function, (initial, forcing, static), cpu, gpu, np.asarray(static.active_mask))
        result.update(kind="synthetic_model", snow_policy="next_day_melt_1e-3_cm" if policy else "legacy")
        report["synthetic"].append(result)
        passed &= result["comparison"]["passed"]
    with jax.default_device(cpu):
        coupling = list(coupling_cases(static))
    for args, boundary in coupling:
        function = lambda soil, constants, want, surface: external_day_kernel(soil, constants, want, surface, bottom_boundary=boundary)
        result, expected, actual = verify_case(function, args, cpu, gpu, np.asarray(static.active_mask))
        result.update(kind="external_water", bottom_boundary=boundary)
        for label, output in (("cpu_ledger", expected), ("gpu_ledger", actual)):
            residuals = verify_water_balance(args[0], output[0], output[1], args[1])
            result[label] = {name: float(np.max(np.abs(value))) for name, value in residuals.items()}
        report["synthetic"].append(result)
        passed &= result["comparison"]["passed"]
    if not passed and not audit_all:
        report.update(passed=False, stopped_before_reference_cases=True)
        return report
    reference_dir = Path(reference_dir)
    for count in days:
        started = time.perf_counter()
        def load():
            with jax.default_device(cpu):
                constants = load_legacy_bundle(reference_dir/"model_inputs.zip")
                radiation = PreparedRadiationProvider(radiation_dir, constants.grid, reference_dir/"model_inputs.zip", required_mask=constants.active_mask)
                weather = load_daily_forcings(reference_dir/"gridmet_weather_2015-01.csv", radiation, count)
                return _host_tree((initialize_legacy_state(constants), weather, constants))
        args = load()
        static = args[2]
        loading = time.perf_counter()-started
        function = lambda state, weather, constants: run_days(state, weather, constants)
        result, _, _ = verify_case(function, args, cpu, gpu, np.asarray(static.active_mask), repeats=repeats, loader=load)
        qualifies = result["comparison"]["passed"]
        result.update(days=count, loading_initialization_s=loading,
            radiation_components_sha256=json.loads((Path(radiation_dir)/"manifest.json").read_text())["components_sha256"],
            workflow_boundary="prepared local inputs to complete host state/history outputs; excludes GRASS preparation and verification")
        for device in ("cpu", "gpu"):
            timing = result[device]
            timing["cold_phase_sum_s"] = loading+sum(timing[k] for k in ("input_transfer_s", "lowering_s", "compilation_s", "first_execution_s", "output_download_s"))
        report["cases"].append(result)
        passed &= qualifies
        if not qualifies and not audit_all:
            break
    report["passed"] = bool(passed)
    report["note"] = "Same WSL JAX version on CPU/CUDA; Windows strict legacy parity remains separately deferred."
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, default=Path("reference/mica-creek-jan2015"))
    parser.add_argument("--radiation-dir", required=True, type=Path)
    parser.add_argument("--days", type=int, nargs="+", choices=(2, 30), default=[2, 30])
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--audit-all", action="store_true")
    parser.add_argument("--snapshot-manifest", type=Path)
    parser.add_argument("--report-path", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.repeats < 1 or args.days[0] != 2:
        parser.error("positive repeats and the two-day case first are required")
    if args.report_path.exists():
        parser.error("report output must be NEW")
    try:
        revision = verify_snapshot(Path(__file__).resolve().parents[2], args.snapshot_manifest) if args.snapshot_manifest else None
        report = verification_report(args.reference_dir, args.radiation_dir, days=args.days, repeats=args.repeats, audit_all=args.audit_all)
        report["source_revision"] = revision
        with args.report_path.open("x", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2)
            stream.write("\n")
    except (AssertionError, ValueError, RuntimeError, OSError) as exc:
        print(f"FAIL: VERIFY-06: {exc}")
        return 1
    print(json.dumps({"passed": report["passed"], "reference_cases": len(report["cases"]), "report": str(args.report_path)}))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
