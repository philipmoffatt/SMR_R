"""Two-layer JAX foundation checks against the frozen Mica Creek input bundle."""

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from smr_jax import initialize_two_layer_reference_state, load_two_layer_legacy_inputs


class TwoLayerFoundationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        bundle = ROOT / "reference" / "mica-creek-jan2015" / "model_inputs.zip"
        cls.inputs = load_two_layer_legacy_inputs(bundle)
        cls.state = initialize_two_layer_reference_state(cls.inputs)

    def test_grid_contract(self) -> None:
        self.assertEqual(self.inputs.grid.shape, (311, 381))
        self.assertEqual(self.inputs.grid.cell_area_m2, 900.0)

    def test_two_layer_state_contract(self) -> None:
        self.assertEqual(self.state.n_layers, 2)
        self.assertEqual(self.state.layer_water_cm.shape, (2, 311, 381))
        self.assertEqual(self.state.total_water_cm.shape, (311, 381))
        # `storage_amt = sat_amt * 0.4` is the legacy initialization
        # expression.  Its valid input domain is wider than the 68,121-cell
        # final-saturation comparison domain, which is shaped by later fluxes.
        self.assertEqual(int(np.asarray(self.state.active_mask).sum()), 72730)

    def test_initial_total_is_preserved_separately_but_balanced_at_start(self) -> None:
        np.testing.assert_allclose(
            np.asarray(self.state.total_water_cm),
            np.asarray(self.state.layer_water_cm).sum(axis=0),
            equal_nan=True,
        )


if __name__ == "__main__":
    unittest.main()
