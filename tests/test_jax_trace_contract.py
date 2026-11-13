"""Frozen trace maps are aligned to the model grid before parity assertions."""

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from smr_jax.io import FrozenTraceProvider
from smr_jax.legacy import load_legacy_bundle


class TraceContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        reference = ROOT / "reference" / "mica-creek-jan2015"
        cls.static = load_legacy_bundle(reference / "model_inputs.zip")
        cls.trace = FrozenTraceProvider(reference / "two_day_internal_trace.zip", cls.static.grid)

    def test_each_surface_trace_map_is_grid_aligned(self):
        for name in ("precip", "rain", "snow", "throughfall", "canopy_evap", "swe", "snowmelt", "q.srad", "q.lw", "q.latent", "q.sensible", "q.total"):
            values = self.trace.raster("2015-01-01", name)
            self.assertEqual(values.shape, self.static.grid.shape, name)
            self.assertTrue(np.isfinite(values[np.asarray(self.static.active_mask)]).any(), name)


if __name__ == "__main__":
    unittest.main()
