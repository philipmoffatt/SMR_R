"""Soil source parity, threshold coverage, budgets, and gate failure detection."""

import ast
import contextlib
from dataclasses import replace
from io import StringIO
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import jax
import jax.numpy as jnp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from smr_jax import verify_soil as gate
from smr_jax.io import ReplayRadiationProvider, load_daily_forcings
from smr_jax.legacy import load_legacy_bundle


class SoilParityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = ROOT / "reference/mica-creek-jan2015"
        cls.source = ROOT / "R/smr_buildup.pl"
        cls.report = gate.verification_report(cls.reference, cls.source)
        cls.oracle = gate.SourceOracle(cls.source)
        cls.static = load_legacy_bundle(cls.reference / "model_inputs.zip")
        cls.forcing = jax.tree_util.tree_map(lambda x: x[0], load_daily_forcings(
            cls.reference / "gridmet_weather_2015-01.csv",
            ReplayRadiationProvider(cls.reference / "two_day_internal_trace.zip", cls.static.grid), 1))

    def test_all_required_frozen_properties_and_classes(self):
        self.assertEqual(self.report["required_cells"], 68121)
        self.assertEqual(self.report["landuse_classes"], [1, 2, 3, 4, 5, 41, 43])
        self.assertEqual(self.report["representative_cells"], 98)
        self.assertEqual(self.report["representative"]["cases"], 588)

    def test_synthetic_thresholds_and_all_stage_ledgers(self):
        self.assertEqual(self.report["synthetic"]["cases"], 918)
        for group in ("synthetic", "representative"):
            errors = self.report[group]["errors"]
            for stage in ("allocation", "redistribution", "ET"):
                for part in ("total", "A", "B"):
                    self.assertLessEqual(errors[f"{stage} {part}"], 1e-12)
            for name in ("road", "allocation", "ET", "drainage", "column"):
                self.assertLessEqual(errors[f"ledger {name}"], 1e-12)

    def test_capacity_rounding_gap_is_reported_and_matches_source(self):
        result = self.report["representative"]
        self.assertGreater(result["max_layer_representation_gap_cm"], 1e-12)
        self.assertLessEqual(result["errors"]["source layer representation"], 1e-12)
        self.assertLessEqual(result["errors"]["ledger column"], 1e-12)

    def test_source_oracle_has_independent_hand_calculated_anchors(self):
        env = dict(gate.synthetic_cases()[0], landuse=2., snowmelt=1.2, throughfall=1.8)
        self.oracle.run("properties", env)
        self.assertEqual(env["root_zone"], 0.)
        self.assertEqual(env["wiltpt_amt"], 0.)
        self.oracle.run("road", env)
        self.assertAlmostEqual(env["road_runoff"], .05)
        self.assertAlmostEqual(env["snowmelt"], 1.18)
        self.assertAlmostEqual(env["throughfall"], 1.77)
        env.update(storage_amt=30., storage_amt_A=8.4, storage_amt_B=21.6,
                   Ksubsurface=.5, perc_daily_flow=0., runoff_daily_flow=0.)
        self.oracle.run("drainage", env)
        self.assertAlmostEqual(env["perc"], .125)
        self.assertAlmostEqual(env["storage_amt"], 29.875)
        self.assertEqual(env["runoff_flow"], 0.)

    def test_zero_root_positive_storage_uses_source_high_et_arm(self):
        env = dict(gate.synthetic_cases()[0], landuse=2., storage_amt=1., storage_amt_A=.3,
                   storage_amt_B=.7, pet=1., actualET_daily_flow=0.)
        self.oracle.run("properties", env)
        self.oracle.run("et", env)
        self.assertAlmostEqual(env["actualET_flow"], .075)
        self.assertAlmostEqual(env["storage_amt"], .925)

    def test_unknown_landuse_uses_water_coefficient_and_zero_a(self):
        env = dict(gate.synthetic_cases()[0], landuse=41., storage_amt=20., storage_amt_A=6.,
                   storage_amt_B=14., pet=1., actualET_daily_flow=0.)
        self.oracle.run("properties", env)
        self.oracle.run("et", env)
        self.assertAlmostEqual(env["actualET_flow"], .1625)
        self.assertEqual(env["storage_amt_A"], 0.)

    def test_source_conditional_is_lazy_and_rejects_other_code(self):
        self.assertEqual(gate._evaluate(ast.parse("choose(1, 2., 1./0.)", mode="eval").body, {}), 2.)
        with self.assertRaises(ValueError):
            gate._evaluate(ast.parse("object.method()", mode="eval").body, {})

    def test_missing_source_assignment_and_changed_timestep_fail(self):
        text = self.source.read_text(encoding="utf-8")
        variants = [text.replace('expression=perc = ', 'removed=perc = ', 1),
                    text.replace('my $temp_time_step = 6.0;', 'my $temp_time_step = 12.0;', 1)]
        with TemporaryDirectory() as directory:
            path = Path(directory) / "source.pl"
            for value in variants:
                path.write_text(value, encoding="utf-8")
                with self.assertRaises(ValueError):
                    gate.SourceOracle(path)

    def test_bad_flux_cannot_pass_source_parity(self):
        original = gate._run_stages

        def changed(*args):
            result = list(original(*args))
            result[7] = result[7] + .01  # Realized ET; leave the state untouched.
            return tuple(result)

        with patch.object(gate, "_run_stages", side_effect=changed):
            with self.assertRaisesRegex(AssertionError, "actual ET"):
                gate.verify_cases(self.oracle, self.static, self.forcing, gate.synthetic_cases())

    def test_wrong_landcover_input_is_detected_on_full_support(self):
        corrupt = replace(self.static, wiltpt_amt_cm=self.static.wiltpt_amt_cm + .01)
        with patch.object(gate, "load_legacy_bundle", return_value=corrupt):
            with self.assertRaisesRegex(AssertionError, "all frozen wiltpt_amt"):
                gate.verification_report(self.reference, self.source)

    def test_nonfinite_wrong_shape_float32_and_empty_values_fail(self):
        expected = np.ones((1, 2), dtype=np.float64)
        for actual in (expected * np.nan, expected * np.inf, expected[:, :1],
                       expected.astype(np.float32), np.empty((1, 0), dtype=np.float64)):
            with self.assertRaises(ValueError):
                gate.check_values("required soil", actual, expected)
        with self.assertRaises(ValueError):
            gate.check_values("required soil", expected, expected * np.nan)

    def test_missing_or_wrong_frozen_support_fails(self):
        original = gate.FrozenTraceProvider.raster

        def missing_cell(provider, date, name):
            value = original(provider, date, name)
            if name == "saturation":
                value.ravel()[np.flatnonzero(np.isfinite(value).ravel())[0]] = np.nan
            return value

        with patch.object(gate.FrozenTraceProvider, "raster", new=missing_cell):
            with self.assertRaisesRegex(ValueError, "full required support"):
                gate.verification_report(self.reference, self.source)

    def test_cli_has_nonzero_exit_for_failed_gate(self):
        with patch.object(gate, "verification_report", side_effect=AssertionError("bad soil flux")):
            with contextlib.redirect_stderr(StringIO()):
                self.assertEqual(gate.main([]), 1)
        with contextlib.redirect_stderr(StringIO()):
            self.assertEqual(gate.main(["--reference-dir", str(self.reference / "missing")]), 1)


if __name__ == "__main__":
    unittest.main()
