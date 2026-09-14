"""Frozen-bundle static contracts, including explicit GRASS NULL propagation."""

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from smr_jax.io import ReplayRadiationProvider
from smr_jax.legacy import load_legacy_bundle


class LegacyBundleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = ROOT / "reference" / "mica-creek-jan2015"
        cls.static = load_legacy_bundle(cls.reference / "model_inputs.zip")

    def test_full_engine_mask_is_legacy_operand_intersection(self):
        self.assertEqual(self.static.grid.shape, (311, 381))
        self.assertEqual(int(np.asarray(self.static.active_mask).sum()), 68121)
        self.assertTrue(np.isfinite(np.asarray(self.static.slope_percent)[np.asarray(self.static.active_mask)]).all())

    def test_replay_field_aligns_to_full_frozen_grid(self):
        field = ReplayRadiationProvider(self.reference / "two_day_internal_trace.zip", self.static.grid).shortwave("2015-01-01")
        self.assertEqual(field.shape, self.static.grid.shape)
        self.assertTrue(np.isfinite(field).any())


if __name__ == "__main__":
    unittest.main()
