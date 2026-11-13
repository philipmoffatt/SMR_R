"""Whole-engine maps/reporting, terminal carry, and strict gate failures."""

import contextlib
from dataclasses import replace
from io import StringIO
import json
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


if __name__ == "__main__":
    unittest.main()
