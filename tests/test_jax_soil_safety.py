"""Focused fixtures for legacy conditional hazards that jnp.where would evaluate eagerly."""

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import jax.numpy as jnp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from smr_jax.fluxes.soil import evapotranspiration
from smr_jax.state import TwoLayerState


class SoilSafetyTest(unittest.TestCase):
    def test_zero_root_zone_preserves_high_storage_et_without_nan(self):
        field = jnp.ones((1, 1), dtype=jnp.float64)
        static = SimpleNamespace(landuse=field * 2, root_zone_cm=jnp.zeros_like(field),
                                 et_reduction_mc=jnp.zeros_like(field), wiltpt_amt_cm=jnp.zeros_like(field),
                                 fieldcap_amt_cm=field, fieldcap_a_cm=field)
        forcing = SimpleNamespace(et_coeff_by_landuse=jnp.ones((6,), dtype=jnp.float64))
        soil = TwoLayerState(field, jnp.stack((field * .5, field * .5)), jnp.ones_like(field, dtype=bool))
        next_soil, et = evapotranspiration(soil, static, forcing, field)
        self.assertEqual(float(et[0, 0]), .25)
        self.assertEqual(float(next_soil.total_water_cm[0, 0]), .75)
        self.assertTrue(np.isfinite(np.asarray(next_soil.total_water_cm)).all())


if __name__ == "__main__":
    unittest.main()
