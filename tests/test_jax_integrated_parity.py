"""Whole-engine maps/reporting, terminal carry, and strict gate failures."""

import contextlib
from dataclasses import replace
from io import StringIO
import json
import math
import os
import jax
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import jax.numpy as jnp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

from smr_jax import verify_integrated as gate
from smr_jax.initialization import initialize_legacy_state
from smr_jax.legacy import load_legacy_bundle


class IntegratedParityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = ROOT/"reference/mica-creek-jan2015"
        cls.static = load_legacy_bundle(cls.reference/"model_inputs.zip")
        original = gate._run_cpu

        def retain(*args, **kwargs):
            cls.final, cls.history = original(*args, **kwargs)
            return cls.final, cls.history

        with patch.object(gate, "_run_cpu", side_effect=retain):
            cls.report = gate.verification_report(cls.reference)

    def test_two_days_all_twenty_maps_and_native_cpu_reports(self):
        self.assertEqual(self.report["days"], 2)
        self.assertEqual(self.report["required_cells"], 68121)
        self.assertEqual(self.report["daily_field_comparisons"], 40)
        self.assertLessEqual(self.report["max_water_error_cm"], 1e-12)
        self.assertEqual(self.report["max_field_errors"]["q.srad"], 0.)

    def test_first_day_longwave_matches_literal_reference_exactly(self):
        provider = gate.FrozenTraceProvider(self.reference/"two_day_internal_trace.zip", self.static.grid)
        expected = provider.raster("2015-01-01", "q.lw")
        mask = np.isfinite(expected)
        np.testing.assert_array_equal(np.asarray(self.history["q_lw"][0])[mask], expected[mask])

    def test_first_day_surface_temperature_preserves_source_division_order(self):
        # This internal state is in the extended capture, rather than the
        # original twenty-map fixture. Check the independent source expression
        # against the actual integrated state without regrouping its divisions.
        swe = np.asarray(self.history["swe_cm"][0])
        energy = np.asarray(self.history["snow_surface_energy"][0])
        mask = np.asarray(self.static.active_mask)
        capacity = 2100.*np.minimum(swe/100., .02)+2100.*np.maximum(0., .02-swe/100.)
        expected = energy/24./capacity
        actual = np.asarray(self.history["snow_surface_temp_c"][0])
        np.testing.assert_array_equal(actual[mask], expected[mask])
        regrouped = energy/(24.*capacity)
        self.assertTrue(np.any(expected[mask] != regrouped[mask]))

    def test_initialization_preserves_inner_fcell_subtraction(self):
        mask = np.asarray(self.static.active_mask)
        initial = initialize_legacy_state(self.static)
        wilt = np.asarray(self.static.wiltpt_mc_a)
        fc = np.asarray(self.static.fieldcap_mc_a)
        depth = np.asarray(self.static.soil_depth_a_cm)
        expected = (wilt+(fc.astype(np.float32)-wilt.astype(np.float32)).astype(np.float64)*.2)*depth
        np.testing.assert_array_equal(np.asarray(initial.soil.layer_water_cm[0])[mask], expected[mask])
        old = (wilt+(fc-wilt)*.2)*depth
        self.assertTrue(np.any(old[mask] != expected[mask]))

    def test_last_substep_flow_is_distinct_from_daily_accumulation(self):
        mask = np.asarray(self.static.active_mask)
        final_flow = np.asarray(self.history["last_lateral_flow_cm"])
        daily = np.asarray(self.history["lateral_flow_cm"])
        self.assertTrue(np.any(daily[:, mask] != final_flow[:, mask]))
        self.assertLessEqual(self.report["max_field_errors"]["lateral_flow"], 1e-12)

    def test_missing_diagnostic_cannot_pass(self):
        bad = dict(self.history)
        bad.pop("storage_amt_b_cm")
        with patch.object(gate, "_run_cpu", return_value=(self.final, bad)):
            with self.assertRaisesRegex(ValueError, "missing/incomplete integrated diagnostic"):
                gate.verification_report(self.reference)

    def test_changed_terminal_state_cannot_hide_behind_correct_history(self):
        soil = replace(self.final.soil, total_water_cm=self.final.soil.total_water_cm+.01)
        with patch.object(gate, "_run_cpu", return_value=(replace(self.final, soil=soil), self.history)):
            with self.assertRaisesRegex(AssertionError, "storage_amt"):
                gate.verification_report(self.reference)

    def test_changed_basin_carry_is_detected(self):
        basin = replace(self.final.basin, baseflow_cms=self.final.basin.baseflow_cms+.001)
        with patch.object(gate, "_run_cpu", return_value=(replace(self.final, basin=basin), self.history)):
            with self.assertRaisesRegex(AssertionError, "terminal baseflow"):
                gate.verification_report(self.reference)

    def test_lost_reference_support_cannot_be_silently_intersected(self):
        original = gate.FrozenTraceProvider.raster

        def lost(provider, date, name):
            values = original(provider, date, name)
            if name == "saturation":
                values.ravel()[np.flatnonzero(np.asarray(self.static.active_mask).ravel())[0]] = np.nan
            return values

        with patch.object(gate.FrozenTraceProvider, "raster", new=lost):
            with self.assertRaisesRegex(ValueError, "required model/reference support differs"):
                gate.verification_report(self.reference)

    def test_nonfinite_types_shapes_and_small_branch_errors_fail(self):
        mask = np.ones((1, 2), dtype=bool)
        expected = np.zeros((1, 2), dtype=np.float64)
        for actual in (expected+np.nan, expected.astype(np.float32), expected[:, :1]):
            with self.assertRaises(ValueError):
                gate.compare_field("storage_amt", actual, expected, mask)
        actual = expected.copy()
        actual[0, 0] = np.nextafter(0., 1.)
        with self.assertRaisesRegex(AssertionError, "snow-presence decisions"):
            gate.compare_field("swe", actual, expected, mask)

    def test_report_rows_require_all_dates_columns_and_finite_values(self):
        row = "1 2015-01-01 "+" ".join(["0"]*22)
        self.assertEqual(len(gate.read_reports(row, 1)), 1)
        for text in (row.replace("2015-01-01", "2015-01-02"), row+"\n"+row,
                     row.rsplit(" ", 1)[0], row.rsplit(" ", 1)[0]+" nan"):
            with self.assertRaises(ValueError):
                gate.read_reports(text, 1)

    def test_thirty_days_requires_complete_independent_capture(self):
        with self.assertRaisesRegex(ValueError, "requires --capture-dir"):
            gate.verification_report(self.reference, days=30)
        with TemporaryDirectory() as directory:
            (Path(directory)/"manifest.json").write_text(json.dumps({"schema_version": 1, "days": 2}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "complete thirty-day extension"):
                gate.validate_extension(directory, self.reference, self.static)

    def test_missing_extended_provenance_and_incomplete_days_fail(self):
        with TemporaryDirectory() as directory:
            (Path(directory)/"manifest.json").write_text(json.dumps({"schema_version": 1, "days": 30,
                "source_sha256": "changed"}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "capture provenance differs"):
                gate.validate_extension(directory, self.reference, self.static)
        with self.assertRaises(ValueError):
            gate.verification_report(self.reference, days=3)

    def test_failed_verification_cli_returns_nonzero(self):
        with patch.object(gate, "verification_report", side_effect=AssertionError("bad daily state")):
            with contextlib.redirect_stderr(StringIO()):
                self.assertEqual(gate.main([]), 1)
        with contextlib.redirect_stderr(StringIO()):
            self.assertEqual(gate.main(["--reference-dir", "missing"]), 1)

    def test_audit_collects_failures_and_still_returns_nonzero(self):
        bad = dict(self.history)
        bad["storage_amt_cm"] = bad["storage_amt_cm"].at[0].add(.01)
        with patch.object(gate, "_run_cpu", return_value=(self.final, bad)):
            report = gate.verification_report(self.reference, collect_failures=True)
        self.assertFalse(report["passed"])
        self.assertTrue(any(f["field"] == "2015-01-01 storage_amt" for f in report["failures"]))
        with patch.object(gate, "verification_report", return_value=report):
            with contextlib.redirect_stdout(StringIO()):
                self.assertEqual(gate.main(["--audit-all"]), 1)

    def test_existing_report_file_is_preserved(self):
        with TemporaryDirectory() as directory:
            report = Path(directory)/"report.json"
            report.write_text("preserve", encoding="utf-8")
            with patch.object(gate, "verification_report", return_value=self.report):
                with contextlib.redirect_stderr(StringIO()):
                    self.assertEqual(gate.main(["--report-path", str(report)]), 1)
            self.assertEqual(report.read_text(), "preserve")

    def test_annual_runoff_excludes_road_loss(self):
        sys.path.insert(0, str(ROOT/"examples"))
        from plant_water import synthetic_case
        _, static, _, _ = synthetic_case()
        static = replace(static, landuse=jnp.full((1, 1), 2., dtype=jnp.float64))
        forcing = jax.tree_util.tree_map(lambda x: x[0], gate.load_daily_forcings(
            self.reference/"gridmet_weather_2015-01.csv",
            gate.ReplayRadiationProvider(self.reference/"two_day_internal_trace.zip", self.static.grid), 1))
        forcing = replace(forcing, tavg_c=jnp.asarray(10.), tdew_c=jnp.asarray(5.), precip_cm=jnp.asarray(5.),
                          pet_cm=jnp.asarray(.5), temperatures_c=jnp.full((4,), 10.), replay_shortwave=jnp.zeros((1, 1)))
        initial = replace(initialize_legacy_state(static), annual_runoff_cm=jnp.full((1, 1), 7.))
        final, history = gate._run_cpu(initial, jax.tree_util.tree_map(lambda x: x[None], forcing), static)
        self.assertGreater(float(history["road_runoff_cm"][0, 0, 0]), 0.)
        np.testing.assert_allclose(history["runoff_cm"], history["road_runoff_cm"], rtol=0., atol=1e-12)
        self.assertEqual(float(final.annual_runoff_cm[0, 0]), 7.)


class NativeMathBoundaryTest(unittest.TestCase):
    def test_unloaded_library_and_wrong_dtype_fail_without_fallback(self):
        from smr_jax import legacy_cpu_math as native
        with self.assertRaisesRegex(ValueError, "float64"):
            native.native_exp(jnp.asarray([0.], dtype=jnp.float32))
        with patch.object(native, "_binding", None):
            with self.assertRaisesRegex(ValueError, "load a native"):
                native.native_exp(jnp.asarray([0.], dtype=jnp.float64))

    def test_build_refuses_existing_artifacts_before_invoking_compiler(self):
        from smr_jax import legacy_cpu_math as native
        with TemporaryDirectory() as directory:
            output = Path(directory)/"preserve.dll"
            for path in (output, output.with_suffix(".dll.json")):
                path.write_text("preserve", encoding="utf-8")
                with patch.object(native.subprocess, "run") as run:
                    with self.assertRaisesRegex(ValueError, "new library and record"):
                        native.build_library("missing-compiler", output)
                    run.assert_not_called()
                self.assertEqual(path.read_text(encoding="utf-8"), "preserve")

    def test_invalid_provenance_is_rejected_before_loading_native_code(self):
        from smr_jax import legacy_cpu_math as native
        with TemporaryDirectory() as directory:
            output = Path(directory)/"invalid.dll"
            output.write_bytes(b"invalid")
            output.with_suffix(".dll.json").write_text(json.dumps({"schema_version": 1}), encoding="utf-8")
            with patch.object(native.ctypes, "CDLL") as load:
                with self.assertRaisesRegex(ValueError, "provenance differs"):
                    native.load_library(output)
                load.assert_not_called()


@unittest.skipUnless(os.environ.get("SMR_JAX_NATIVE_MATH_LIBRARY"),
                     "optional native CPU profile: set SMR_JAX_NATIVE_MATH_LIBRARY")
class NativeMathExecutionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from smr_jax import legacy_cpu_math as native
        cls.native = native
        cls.path = Path(os.environ["SMR_JAX_NATIVE_MATH_LIBRARY"])
        cls.metadata = native.load_library(cls.path)

    def test_jitted_scalar_exp_special_values_empty_arrays_and_accuracy(self):
        from smr_jax.fluxes.legacy_math import compile_legacy_cpu
        run = compile_legacy_cpu(self.native.native_exp)
        samples = np.linspace(-700., 700., 4097, dtype=np.float64)
        actual = np.asarray(run(jnp.asarray(samples)))
        expected = np.asarray([math.exp(float(v)) for v in samples])
        # This checks scalar math accuracy, not the unchanged integrated water
        # budgets. Platform libm results can differ by one rounding unit.
        np.testing.assert_allclose(actual, expected, rtol=3*np.finfo(np.float64).eps, atol=0.)
        special = np.asarray(run(jnp.asarray([-np.inf, np.inf, np.nan, 0.], dtype=jnp.float64)))
        self.assertEqual(special[0], 0.)
        self.assertEqual(special[1], np.inf)
        self.assertTrue(np.isnan(special[2]))
        self.assertEqual(special[3], 1.)
        self.assertEqual(run(jnp.empty((0,), dtype=jnp.float64)).shape, (0,))
        self.assertEqual(float(run(jnp.asarray(0., dtype=jnp.float64))), 1.)

    def test_native_execution_contains_no_python_callback(self):
        values = jnp.linspace(-13., 2., 1024, dtype=jnp.float64)
        lowered = jax.jit(self.native.native_exp).lower(values).as_text()
        self.assertIn("smr_legacy_exp_", lowered)
        self.assertNotIn("callback", lowered.lower())
        # Assert that selecting the profile actually changes the math backend.
        actual = np.asarray(jax.jit(self.native.native_exp)(values))
        self.assertTrue(np.any(actual != np.asarray(jax.jit(jnp.exp)(values))))

    def test_batched_values_and_gradient_keep_elementwise_exp_semantics(self):
        values = jnp.asarray([[-2., 0.], [1., 2.]], dtype=jnp.float64)
        expected = np.asarray([[math.exp(float(v)) for v in row] for row in values])
        actual = jax.jit(jax.vmap(self.native.native_exp))(values)
        derivative = jax.jit(jax.grad(lambda v: jnp.sum(self.native.native_exp(v))))(values)
        np.testing.assert_allclose(actual, expected, rtol=3*np.finfo(np.float64).eps, atol=0.)
        np.testing.assert_allclose(derivative, expected, rtol=3*np.finfo(np.float64).eps, atol=0.)

    def test_native_profile_passes_original_two_day_integrated_gate(self):
        report = gate.verification_report(ROOT/"reference/mica-creek-jan2015", native_math_library=self.path)
        self.assertTrue(report["passed"])
        self.assertEqual(report["math_profile"]["library_sha256"], self.metadata["library_sha256"])
        self.assertEqual(report["daily_field_comparisons"], 40)

    def test_provenance_and_loaded_library_binding_cannot_be_changed(self):
        from hashlib import sha256
        returned = self.native.load_library(self.path)
        returned["command"].append("must not mutate the registered build record")
        self.assertEqual(self.native.load_library(self.path), self.metadata)
        with TemporaryDirectory() as directory:
            other = Path(directory)/"other.dll"
            other.write_bytes(self.path.read_bytes())
            metadata = dict(self.metadata)
            record = other.with_suffix(".dll.json")
            record.write_text(json.dumps(metadata), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "already bound"):
                self.native.load_library(other)
            other.write_bytes(b"corrupt")
            self.assertNotEqual(sha256(other.read_bytes()).hexdigest(), metadata["library_sha256"])
            with self.assertRaisesRegex(ValueError, "provenance differs"):
                self.native.load_library(other)


if __name__ == "__main__":
    unittest.main()
