"""Independent GRASS routing evidence, basin recurrence and failure detection."""

import ast
import contextlib
from dataclasses import replace
from io import StringIO
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import jax.numpy as jnp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
sys.path.insert(0, str(ROOT/"examples"))

from plant_water import synthetic_case
from smr_jax import verify_routing as gate
from smr_jax.fluxes.routing import legacy_receive_gather
from smr_jax.legacy import _grass_default_slope_degrees
from smr_jax.state import BasinState
from smr_jax.verify_soil import _evaluate


class RoutingParityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = ROOT/"reference/mica-creek-jan2015"
        cls.source = ROOT/"R/smr_buildup.pl"
        cls.capture = ROOT/"tests/fixtures/jax_verify03_grass"
        cls.report = gate.verification_report(cls.reference, cls.source, cls.capture)

    def test_required_grass_static_fields_and_representative_flows(self):
        self.assertEqual(self.report["required_cells"], 68121)
        self.assertEqual(self.report["grass_representatives"], 386)
        self.assertEqual(len(self.report["lateral_days"]), 2)
        for day in self.report["lateral_days"]:
            self.assertLessEqual(day["max_grass_error_cm"], 1e-12)
            self.assertLessEqual(day["max_error_cm"], 1e-12)
            self.assertLessEqual(day["domain_ledger_error_cell_cm"], 1e-10)

    def test_all_gather_orientations_masks_and_skinny_borders(self):
        self.assertEqual(self.report["gather"]["cases"], 64)
        self.assertEqual(self.report["gather"]["max_error_cm"], 0.)

    def test_finite_inactive_donor_cannot_enter_active_receiver(self):
        flow = jnp.asarray([[10., 0.]], dtype=jnp.float64)
        weights = jnp.zeros((8, 1, 2), dtype=jnp.float64).at[6, 0, 1].set(1.)
        mask = jnp.asarray([[False, True]])
        np.testing.assert_array_equal(legacy_receive_gather(flow, weights, mask), np.zeros((1, 2)))

    def test_receiver_weights_can_duplicate_or_export_water(self):
        flow = jnp.asarray([[0., 4., 0.]], dtype=jnp.float64)
        weights = jnp.zeros((8, 1, 3), dtype=jnp.float64).at[2, 0, 0].set(.75).at[6, 0, 2].set(.75)
        mask = jnp.ones((1, 3), dtype=bool)
        result = gate._run_gather(flow, weights, mask)
        np.testing.assert_array_equal(result, [[3., 0., 3.]])
        self.assertEqual(float(result.sum()-flow.sum()), 2.)
        exported = gate._run_gather(flow, weights*0., mask)
        self.assertEqual(float(exported.sum()-flow.sum()), -4.)

    def test_synthetic_lateral_budget_accepts_literal_receiver_gain(self):
        synthetic = self.report["synthetic"]
        self.assertGreater(synthetic["domain_change_cell_cm"], 0.)
        self.assertLessEqual(synthetic["domain_ledger_error_cell_cm"], 1e-10)

    def test_flat_planar_slope_units_rounding_and_null_edges(self):
        y, x = np.indices((5, 5), dtype=np.float64)
        slope = _grass_default_slope_degrees(100.+30.*x, 30.)
        np.testing.assert_array_equal(slope[1:-1, 1:-1], np.full((3, 3), 45.))
        self.assertTrue(np.isnan(slope[0]).all())
        flat = _grass_default_slope_degrees(np.ones((5, 5)), 30.)
        np.testing.assert_array_equal(flat[1:-1, 1:-1], np.zeros((3, 3)))
        curved = _grass_default_slope_degrees(100.+x*x+2*y, 30.)
        np.testing.assert_array_equal(curved, curved.astype(np.float32).astype(np.float64))
        dem = 100.+30.*x
        dem[2, 2] = np.nan
        self.assertTrue(np.isnan(_grass_default_slope_degrees(dem, 30.))[1:4, 1:4].all())

    def test_raster_only_products_round_before_double_promotion(self):
        expression = ast.parse("choose(1., a*b+a*c, 0.)", mode="eval").body
        env = {"a": np.float32(.1234567), "b": np.float32(30.), "c": np.float32(70.)}
        expected = np.float64(np.float32(env["a"]*env["b"]+env["a"]*env["c"]))
        self.assertEqual(_evaluate(expression, env), expected)
        mixed = _evaluate(ast.parse("a*6/24.0", mode="eval").body, env)
        self.assertEqual(mixed, np.float64(np.float32(env["a"]*np.float32(6)))/24.)

    def test_all_thirty_frozen_discharge_baseflow_and_volume_values(self):
        result = self.report["basin"]
        self.assertEqual(result["days"], 30)
        self.assertLessEqual(result["max_discharge_error_cms"], 1e-15)
        self.assertLessEqual(result["max_baseflow_error_cms"], 1e-15)
        self.assertLessEqual(result["max_reservoir_error"], 1e-12)

    def test_basin_update_uses_previous_baseflow_and_asymmetric_units(self):
        _, static, _, _ = synthetic_case(shape=(1, 2))
        static = replace(static, basin_reservoir_coeff=jnp.asarray([.2]),
                         basin_area_cells=jnp.asarray([999999.]))
        before = BasinState(jnp.asarray([10.]), jnp.asarray([.4]))
        state, discharge = gate._run_basin(before, static, jnp.asarray([[1., 3.]]), jnp.asarray([[2., 6.]]))
        volume = 10.+8./100./86400.-.4
        np.testing.assert_allclose(state.reservoir_volume, [volume], rtol=0., atol=1e-14)
        np.testing.assert_allclose(state.baseflow_cms, [.2*volume], rtol=0., atol=1e-14)
        np.testing.assert_allclose(discharge, [.2*volume+4./100.*900./86400.], rtol=0., atol=1e-14)
        # Empty inflow still subtracts yesterday's baseflow before the update.
        next_state, _ = gate._run_basin(state, static, jnp.zeros((1, 2)), jnp.zeros((1, 2)))
        np.testing.assert_allclose(next_state.reservoir_volume, [.8*volume], rtol=0., atol=1e-14)

    def test_separate_basin_math_masks_and_zero_coefficient(self):
        _, static, _, _ = synthetic_case(shape=(2, 2))
        static = replace(static, basin_ids=jnp.asarray([3., 8., 99.]),
            basin_reservoir_coeff=jnp.asarray([.1, 0., .2]), watershed=jnp.asarray([[3., 8.], [3., 8.]]),
            active_mask=jnp.asarray([[True, True], [False, True]]))
        before = BasinState(jnp.asarray([2., 4., 6.]), jnp.asarray([.1, .2, .3]))
        runoff = jnp.asarray([[1., 2.], [jnp.nan, 4.]])
        perc = jnp.asarray([[2., 3.], [jnp.nan, 5.]])
        after, q = gate._run_basin(before, static, runoff, perc)
        volume = np.asarray([2.+2./100./86400.-.1, 4.+8./100./86400.-.2, 6.-.3])
        base = volume*np.asarray([.1, 0., .2])
        np.testing.assert_allclose(after.reservoir_volume, volume, rtol=0., atol=1e-14)
        np.testing.assert_allclose(after.baseflow_cms, base, rtol=0., atol=1e-14)
        np.testing.assert_allclose(q, base+np.asarray([1., 6., 0.])/100.*900./86400., rtol=0., atol=1e-14)

    def test_corrupt_capture_binary_or_changed_source_fails(self):
        metadata = self.capture.with_suffix(".json").read_bytes()
        with TemporaryDirectory() as directory:
            path = Path(directory)/"capture"
            path.with_suffix(".json").write_bytes(metadata)
            path.with_suffix(".npz").write_bytes(b"corrupt")
            with self.assertRaisesRegex(ValueError, "capture inventory"):
                gate.load_capture(path, self.source, self.reference)
        with TemporaryDirectory() as directory:
            source = Path(directory)/"changed.pl"
            source.write_text("different source", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "source differs"):
                gate.load_capture(self.capture, source, self.reference)

    def test_changed_slope_and_missing_support_cannot_pass(self):
        original = gate.load_legacy_bundle

        def changed(*args):
            static = original(*args)
            return replace(static, slope_percent=static.slope_percent+.01)

        with patch.object(gate, "load_legacy_bundle", side_effect=changed):
            with self.assertRaisesRegex(AssertionError, "captured GRASS slope"):
                gate.verification_report(self.reference, self.source, self.capture)
        static = original(self.reference/"model_inputs.zip")
        corrupt = replace(static, basin_reservoir_coeff=static.basin_reservoir_coeff+.01)
        with patch.object(gate, "load_legacy_bundle", return_value=corrupt):
            with self.assertRaisesRegex(AssertionError, "captured reservoir configuration"):
                gate.verification_report(self.reference, self.source, self.capture)
        meta, arrays = gate.load_capture(self.capture, self.source, self.reference)
        arrays["support"] = arrays["support"].copy()
        arrays["support"].ravel()[np.flatnonzero(arrays["support"].ravel())[0]] = False
        with patch.object(gate, "load_capture", return_value=(meta, arrays)):
            with self.assertRaisesRegex(ValueError, "required support"):
                gate.verification_report(self.reference, self.source, self.capture)

    def test_changed_gather_or_recurrence_fails_gate(self):
        original = gate._run_gather
        with patch.object(gate, "_run_gather", side_effect=lambda *args: original(*args)+.01):
            with self.assertRaises(AssertionError):
                gate.check_gather()
        static = gate.load_legacy_bundle(self.reference/"model_inputs.zip")
        basin = gate._run_basin

        def wrong(*args):
            state, discharge = basin(*args)
            return state, discharge+.001

        with patch.object(gate, "_run_basin", side_effect=wrong):
            with self.assertRaisesRegex(AssertionError, "frozen discharge"):
                gate.check_frozen_basin(self.reference, static)

    def test_missing_capture_or_failed_verification_cli_is_nonzero(self):
        with patch.object(gate, "verification_report", side_effect=AssertionError("bad routing")):
            with contextlib.redirect_stderr(StringIO()):
                self.assertEqual(gate.main([]), 1)
        with contextlib.redirect_stderr(StringIO()):
            self.assertEqual(gate.main(["--capture", str(self.capture.parent/"missing")]), 1)


if __name__ == "__main__":
    unittest.main()
