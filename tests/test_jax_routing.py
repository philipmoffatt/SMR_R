"""Fixtures that prove the legacy receiver-weighted routing orientation."""

import sys
import unittest
from pathlib import Path

import jax.numpy as jnp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from smr_jax.fluxes.routing import legacy_receive_gather


class LegacyReceiveGatherTest(unittest.TestCase):
    def test_each_named_receiver_weight_reads_its_literal_neighbor(self):
        # [N, NE, E, SE, S, SW, W, NW], in row/column coordinates.
        offsets = ((-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1))
        active = jnp.ones((3, 3), dtype=bool)
        for index, (row, col) in enumerate(offsets):
            flow = jnp.zeros((3, 3), dtype=jnp.float64).at[1 + row, 1 + col].set(7.0)
            weights = jnp.zeros((8, 3, 3), dtype=jnp.float64).at[index, 1, 1].set(0.25)
            received = legacy_receive_gather(flow, weights, active)
            self.assertEqual(float(received[1, 1]), 1.75, msg=f"direction {index}")

    def test_null_and_border_donors_are_zero_before_weighting(self):
        flow = jnp.array([[jnp.nan, 2.0], [3.0, 4.0]], dtype=jnp.float64)
        weights = jnp.ones((8, 2, 2), dtype=jnp.float64)
        received = legacy_receive_gather(flow, weights, jnp.ones((2, 2), dtype=bool))
        self.assertTrue(np.isfinite(np.asarray(received)).all())
        self.assertEqual(float(received[0, 0]), 9.0)  # east + south + southeast


if __name__ == "__main__":
    unittest.main()
