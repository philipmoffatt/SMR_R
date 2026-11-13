"""Operator-selected 1e-3 cm snow cleanup: delayed release, no lost water."""

from dataclasses import replace
from functools import partial
from pathlib import Path
from types import SimpleNamespace
import sys
import unittest
from unittest.mock import patch

import jax
import jax.numpy as jnp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "examples")]

from plant_water import synthetic_case
from smr_jax.contracts import DailyForcing
from smr_jax.engine import daily_step, run_days
from smr_jax.fluxes.snow import snow_energy
from smr_jax.initialization import initialize_legacy_state


class SnowResidueTest(unittest.TestCase):
    threshold = 1e-3

    def case(self, swe=0., shape=(1, 1)):
        _, static, _, _ = synthetic_case(shape=shape)
        static = replace(static, ksubsurface=jnp.zeros(shape))
        state = initialize_legacy_state(static)
        field = lambda value: jnp.full(shape, value, dtype=jnp.float64)
        snow = replace(state.snow, swe_cm=field(swe), previous_swe_cm=field(swe),
                       ice_cm=field(swe), liquid_cm=field(0.), total_energy=field(0.),
                       surface_energy=field(0.), surface_temp_c=field(0.))
        # At 0 C all vapor terms cancel. Shortwave balances the explicit
        # longwave loss, retaining a small fresh snowfall for the timing test.
        cover = np.asarray(static.canopy_cover)
        emissivity = np.maximum(.72, (9.2e-6 * 273.15**2) * (1 - cover) + .92 * cover)
        shortwave = 5.67e-8 * 86400 / 1000 * (.98 - emissivity) * 273.15**4
        forcing = DailyForcing(
            *[jnp.asarray(v) for v in (2015., 1., 1., 1., 0., 0., 0., .2486 / 1.8281)],
            jnp.zeros(4), jnp.asarray(0.), jnp.asarray(100.), jnp.asarray(10.),
            jnp.zeros(6), jnp.asarray(shortwave))
        return replace(state, snow=snow), static, forcing, field

    def snow_step(self, state, static, forcing, field, snowfall=0., enabled=True):
        return snow_energy(state, static, forcing, field(0.), field(0.), field(snowfall), field(0.),
                           melt_small_snow_next_day=enabled)

    def test_threshold_below_at_and_above(self):
        for amount in (0., np.nextafter(1e-3, 0.), 1e-3, np.nextafter(1e-3, np.inf), .1):
            with self.subTest(amount=amount):
                state, static, forcing, field = self.case(amount)
                after, d = self.snow_step(state, static, forcing, field)
                expected = amount if 0 < amount <= self.threshold else 0.
                self.assertEqual(float(d["residual_snowmelt_cm"][0, 0]), expected)
                # Independent snow mass ledger, including the released residue.
                np.testing.assert_allclose(after.swe_cm + d["snowmelt_cm"],
                                           state.snow.swe_cm + d["condens_cm"], rtol=0, atol=1e-14)

    def test_new_residue_waits_until_following_day(self):
        state, static, forcing, field = self.case()
        first, d1 = self.snow_step(state, static, forcing, field, snowfall=.5e-3)
        self.assertEqual(float(d1["residual_snowmelt_cm"][0, 0]), 0.)
        self.assertGreater(float(first.swe_cm[0, 0]), 0.)
        self.assertLessEqual(float(first.swe_cm[0, 0]), self.threshold)
        second, d2 = self.snow_step(replace(state, snow=first), static, forcing, field)
        self.assertEqual(float(d2["residual_snowmelt_cm"][0, 0]), float(first.swe_cm[0, 0]))
        self.assertEqual(float(second.swe_cm[0, 0]), 0.)
        third, d3 = self.snow_step(replace(state, snow=second), static, forcing, field)
        self.assertEqual(float(d3["residual_snowmelt_cm"][0, 0]), 0.)
        self.assertAlmostEqual(float(second.swe_cm.sum() + d1["snowmelt_cm"].sum() + d2["snowmelt_cm"].sum()), .5e-3)

    def test_released_components_do_not_reappear_in_snow(self):
        state, static, forcing, field = self.case(.9e-3)
        state = replace(state, snow=replace(state.snow, ice_cm=field(.7e-3), liquid_cm=field(.2e-3)))
        after, d = self.snow_step(state, static, forcing, field)
        for value in (after.swe_cm, after.previous_swe_cm, after.ice_cm, after.liquid_cm):
            np.testing.assert_allclose(value, 0., atol=1e-14)
        self.assertAlmostEqual(float(d["snowmelt_cm"].sum()), .9e-3)

    def test_fresh_snow_does_not_absorb_yesterdays_residue(self):
        state, static, forcing, field = self.case(.8e-3)
        after, d = self.snow_step(state, static, forcing, field, snowfall=.5)
        self.assertEqual(float(d["residual_snowmelt_cm"].sum()), .8e-3)
        np.testing.assert_allclose(after.swe_cm + d["snowmelt_cm"],
                                   state.snow.swe_cm + .5 + d["condens_cm"], rtol=0, atol=1e-14)
        self.assertGreater(float(after.swe_cm.sum()), self.threshold)

    def test_tiny_positive_residue_uses_bare_ground_resistance(self):
        tiny, static, forcing, field = self.case(1e-16)
        empty, _, _, _ = self.case(0.)
        forcing = replace(forcing, replay_shortwave=field(0.))
        _, cleared = snow_energy(tiny, static, forcing, field(2.), field(0.), field(0.), field(0.),
                                  melt_small_snow_next_day=True)
        _, bare = snow_energy(empty, static, forcing, field(2.), field(0.), field(0.), field(0.))
        np.testing.assert_array_equal(cleared["q_sensible"], bare["q_sensible"])
        np.testing.assert_array_equal(cleared["q_latent"], bare["q_latent"])

    def test_negative_and_inactive_values_are_not_cleanup_sources(self):
        state, static, forcing, field = self.case(.8e-3, shape=(1, 2))
        static = replace(static, active_mask=jnp.asarray([[True, False]]))
        state = replace(state, snow=replace(state.snow, swe_cm=jnp.asarray([[-1e-16, .8e-3]])))
        _, d = self.snow_step(state, static, forcing, field)
        np.testing.assert_array_equal(d["residual_snowmelt_cm"], 0.)

    def test_legacy_mode_has_no_deferred_release(self):
        state, static, forcing, field = self.case(.8e-3)
        _, d = self.snow_step(state, static, forcing, field, enabled=False)
        np.testing.assert_array_equal(d["residual_snowmelt_cm"], 0.)

    def test_cold_day_still_delivers_all_residue_to_soil(self):
        state, static, forcing, field = self.case(.8e-3)
        forcing = replace(forcing, tavg_c=jnp.asarray(-5.), temperatures_c=jnp.full(4, -5.))
        after, d = daily_step(state, forcing, static, melt_small_snow_next_day=True)
        self.assertEqual(float(d["residual_snowmelt_to_soil_cm"].sum()), .8e-3)
        self.assertAlmostEqual(float((after.soil.total_water_cm - state.soil.total_water_cm).sum()), .8e-3)
        np.testing.assert_allclose(after.snow.swe_cm, 0., atol=1e-14)

    def test_road_partition_preserves_the_residue_tracer(self):
        state, static, forcing, field = self.case(.8e-3)
        static = replace(static, landuse=field(2.))
        snow_flux = {"snowmelt_cm": field(2.), "residual_snowmelt_cm": field(.8e-3)}
        with patch("smr_jax.engine.snow_energy", return_value=(state.snow, snow_flux)):
            _, d = daily_step(state, forcing, static, melt_small_snow_next_day=True)
        np.testing.assert_allclose(d["residual_snowmelt_to_soil_cm"] + d["residual_snowmelt_road_runoff_cm"],
                                   .8e-3, rtol=0, atol=1e-16)
        self.assertGreater(float(d["residual_snowmelt_road_runoff_cm"].sum()), 0.)

    def test_multiday_jit_releases_once(self):
        state, static, forcing, field = self.case(.8e-3)
        forcing = replace(forcing, tavg_c=jnp.asarray(-5.), temperatures_c=jnp.full(4, -5.))
        stacked = jax.tree.map(lambda a: jnp.stack((a, a, a)), forcing)
        go = jax.jit(partial(run_days, static=static, melt_small_snow_next_day=True))
        final, history = go(state, stacked)
        np.testing.assert_allclose(history["residual_snowmelt_cm"][:, 0, 0], [.8e-3, 0., 0.], rtol=0, atol=1e-16)
        self.assertAlmostEqual(float((final.soil.total_water_cm - state.soil.total_water_cm).sum()), .8e-3)


if __name__ == "__main__":
    unittest.main()
