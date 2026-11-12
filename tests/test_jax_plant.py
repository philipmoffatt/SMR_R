"""Synthetic coupling, independent edge ledgers, and legacy preservation checks."""

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
from smr_jax.fluxes import soil as legacy
from smr_jax.plant import (SurfaceWater, VegetationDemand, demand_from_plant, external_day_kernel,
                           external_water_step, require_realized_uptake, soil_to_plant,
                           validate_external_inputs, verify_water_balance)
from smr_jax.state import TwoLayerState


class PlantWaterTest(unittest.TestCase):
    def run_case(self, soil, static, demand, surface, **kwargs):
        after, d = external_water_step(soil, static, demand, surface, allow_unmet=True, **kwargs)
        # Test-owned arithmetic uses only initial/final storage and measured
        # inputs/outputs. No model residual or inferred lateral input is used.
        delta = np.asarray(after.layer_water_cm - soil.layer_water_cm)
        vertical = np.asarray(d["vertical_down_cm"] - d["vertical_up_cm"])[0]
        sources = np.stack((np.asarray(d["infiltration_cm"]) - vertical, vertical))
        sinks = np.stack((np.asarray(d["runoff_cm"]), np.asarray(d["bottom_drainage_cm"])))
        expected = sources - sinks + np.asarray(d["lateral_in_cm"] - d["lateral_out_cm"] - d["transp_cm"] - d["soil_evap_cm"])
        mask = np.asarray(static.active_mask)
        np.testing.assert_allclose(delta[:, mask], expected[:, mask], atol=1e-10, rtol=0)
        return after, d

    def zero_demand(self, demand):
        return jax.tree.map(jnp.zeros_like, demand)

    def no_rain(self, surface):
        return jax.tree.map(jnp.zeros_like, surface)

    def test_external_replaces_canopy_and_et_without_duplication(self):
        soil, static, demand, surface = synthetic_case()
        with patch("smr_jax.fluxes.surface.canopy", side_effect=AssertionError("internal canopy")), patch(
                "smr_jax.fluxes.soil.evapotranspiration", side_effect=AssertionError("internal ET")):
            after, d = self.run_case(soil, static, demand, surface, bottom_boundary="closed")
        np.testing.assert_allclose(d["transp_cm"], demand.transp_cm, atol=1e-12)
        np.testing.assert_allclose(d["soil_evap_cm"], demand.soil_evap_cm, atol=1e-12)
        self.assertAlmostEqual(float(after.total_water_cm[0, 0]), 20 + .4 - .02 - .2 - .04)
        require_realized_uptake(d)

    def test_zero_and_dry_soil(self):
        for moisture in (0., .05, .10):
            with self.subTest(moisture=moisture):
                soil, static, demand, surface = synthetic_case(moisture=moisture)
                after, d = self.run_case(soil, static, demand, self.no_rain(surface), bottom_boundary="closed")
                np.testing.assert_array_equal(d["transp_cm"], 0.)
                np.testing.assert_array_equal(d["soil_evap_cm"], 0.)
                self.assertAlmostEqual(float(after.total_water_cm[0, 0]), float(soil.total_water_cm[0, 0]))
                with self.assertRaisesRegex(ValueError, "reconcile plant"):
                    require_realized_uptake(d)

    def test_zero_demand_preserves_dry_storage(self):
        soil, static, demand, surface = synthetic_case(moisture=.05)
        after, d = external_water_step(soil, static, self.zero_demand(demand), self.no_rain(surface))
        np.testing.assert_allclose(after.layer_water_cm, soil.layer_water_cm, atol=1e-12)

    def test_excess_demand_is_bounded_proportionally_and_refused_by_default(self):
        soil, static, demand, surface = synthetic_case(moisture=.11)
        demand = replace(demand, transp_cm=jnp.full_like(demand.transp_cm, 10.),
                         soil_evap_cm=jnp.full_like(demand.soil_evap_cm, 5.), interception_cm=surface.rain_cm * 0)
        surface = self.no_rain(surface)
        after, d = self.run_case(soil, static, demand, surface, bottom_boundary="closed")
        np.testing.assert_allclose(d["transp_cm"], 2 * d["soil_evap_cm"], atol=1e-12)
        np.testing.assert_allclose(after.layer_water_cm[:, 0, 0], [3., 7.], atol=1e-12)
        np.testing.assert_allclose(d["transp_cm"] + d["soil_evap_cm"], soil.layer_water_cm - after.layer_water_cm)
        with self.assertRaisesRegex(ValueError, "unmet transp"):
            external_water_step(soil, static, demand, surface, bottom_boundary="closed")

    def test_wetting_pulse_supplies_same_day_demand(self):
        soil, static, demand, surface = synthetic_case(moisture=.10)
        after, d = self.run_case(soil, static, demand, surface, bottom_boundary="closed")
        self.assertGreater(float(d["transp_cm"].sum()), 0.)
        self.assertGreater(float(d["unmet_transp_cm"][1, 0, 0]), 0.)
        self.assertGreaterEqual(float(after.layer_water_cm.min()), 0.)

    def test_interception_is_limited_to_liquid_not_snowfall(self):
        soil, static, demand, surface = synthetic_case()
        demand = replace(demand, interception_cm=surface.rain_cm * 10)
        _, d = self.run_case(soil, static, demand, surface)
        np.testing.assert_array_equal(d["interception_evap_cm"], surface.rain_cm)
        np.testing.assert_array_equal(d["infiltration_cm"], 0.)
        self.assertGreater(float(d["unmet_interception_cm"].sum()), 0.)

    def test_wet_pulse_has_explicit_returnflow_and_transfers(self):
        soil, static, demand, surface = synthetic_case(moisture=.40)
        surface = replace(surface, rain_cm=surface.rain_cm * 100, precip_cm=surface.precip_cm * 100)
        after, d = self.run_case(soil, static, self.zero_demand(demand), surface, bottom_boundary="closed")
        np.testing.assert_allclose(after.layer_water_cm[:, 0, 0], [13.5, 31.5], atol=1e-12)
        self.assertAlmostEqual(float(d["runoff_cm"][0, 0]), 35.)
        self.assertGreater(float(d["vertical_down_cm"].sum()), 0.)

    def test_closed_and_legacy_draining_bottom(self):
        soil, static, demand, surface = synthetic_case(moisture=.40)
        demand, surface = self.zero_demand(demand), self.no_rain(surface)
        closed, dc = self.run_case(soil, static, demand, surface, bottom_boundary="closed")
        drained, dd = self.run_case(soil, static, demand, surface)
        self.assertEqual(float(dc["bottom_drainage_cm"].sum()), 0.)
        self.assertAlmostEqual(float(dd["bottom_drainage_cm"].sum()), .5)
        self.assertAlmostEqual(float(closed.total_water_cm.sum() - drained.total_water_cm.sum()), .5)

    def test_flow_cap_and_upward_returnflow_under_strong_routing(self):
        soil, static, demand, surface = synthetic_case(shape=(2, 1), moisture=.45)
        # Only the first row exports water to the second. Large conductivities
        # stress donor limits without changing the source profile properties.
        static = replace(static, slope_percent=jnp.full((2, 1), 45.),
                         routing_weights=static.routing_weights.at[0, 1, 0].set(1.),
                         flowunits=static.flowunits.at[1, 0].set(0.),
                         ksat_mpore_a=jnp.full((2, 1), 1e8), ksat_mpore_b=jnp.full((2, 1), 1e8))
        after, d = self.run_case(soil, static, self.zero_demand(demand), self.no_rain(surface), bottom_boundary="closed")
        self.assertGreater(float(d["runoff_cm"].sum()), 0.)
        self.assertGreater(float(d["vertical_up_cm"].sum()), 0.)
        self.assertGreaterEqual(float(after.layer_water_cm.min()), 0.)

    def test_equal_fieldcapacity_and_saturation_is_supported(self):
        soil, static, demand, surface = synthetic_case()
        static = replace(static, sat_mc_a=static.fieldcap_mc_a, sat_mc_b=static.fieldcap_mc_b,
                         sat_amt_cm=static.fieldcap_amt_cm)
        self.run_case(soil, static, self.zero_demand(demand), surface)

    def test_literal_legacy_soil_parity_when_bounds_are_inactive(self):
        soil, static, demand, surface = synthetic_case(moisture=.35)
        demand = self.zero_demand(demand)
        after, _ = self.run_case(soil, static, demand, surface)
        expected = soil
        for _ in range(4):
            expected = legacy.allocate_input(expected, static, surface.rain_cm * .25)
            expected, _ = legacy.lateral_route(expected, static)
            a, b = legacy._reallocate(expected.total_water_cm, expected.layer_water_cm[0], static)
            expected = TwoLayerState(expected.total_water_cm, jnp.stack((a, b)), expected.active_mask)
            expected, _, _, _ = legacy.percolation_and_runoff(expected, static)
        np.testing.assert_allclose(after.layer_water_cm, expected.layer_water_cm, atol=1e-12, rtol=0)

    def test_routing_edges_and_outlet_are_independent(self):
        soil, static, demand, surface = synthetic_case(shape=(3, 1), moisture=.35)
        weights = static.routing_weights.at[0, 1:, 0].set(1.)
        static = replace(static, slope_percent=jnp.full((3, 1), 10.), routing_weights=weights)
        after, d = self.run_case(soil, static, self.zero_demand(demand), self.no_rain(surface), bottom_boundary="closed")
        out, incoming, export = (np.asarray(d[n]) for n in ("lateral_out_cm", "lateral_in_cm", "lateral_export_cm"))
        np.testing.assert_array_equal(incoming[:, 0], 0.)
        np.testing.assert_allclose(incoming[:, 1:], out[:, :-1], atol=1e-12)
        np.testing.assert_allclose(export[:, -1], out[:, -1], atol=1e-12)
        np.testing.assert_array_equal(export[:, :-1], 0.)
        self.assertAlmostEqual(float(soil.total_water_cm.sum() - after.total_water_cm.sum()), float(export.sum()))
        corrupted = {**d, "lateral_in_cm": d["lateral_in_cm"] + .001}
        with self.assertRaises(FloatingPointError):
            verify_water_balance(soil, after, corrupted, static)

    def test_each_routing_orientation_conserves(self):
        offsets = ((-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1))
        for k, (dr, dc) in enumerate(offsets):
            soil, static, demand, surface = synthetic_case(shape=(3, 3), moisture=.35)
            weights = static.routing_weights.at[k, 1, 1].set(.7)
            static = replace(static, slope_percent=jnp.full((3, 3), 10.), routing_weights=weights)
            _, d = self.run_case(soil, static, self.zero_demand(demand), self.no_rain(surface))
            np.testing.assert_allclose(d["lateral_in_cm"][:, 1, 1], .7 * d["lateral_out_cm"][:, 1 + dr, 1 + dc])

    def test_routing_water_duplication_is_rejected(self):
        soil, static, demand, surface = synthetic_case(shape=(3, 3))
        static = replace(static, routing_weights=jnp.ones((8, 3, 3)))
        with self.assertRaisesRegex(ValueError, "duplicate donor water"):
            external_water_step(soil, static, demand, surface)

    def test_masked_nan_holes_are_not_sources_or_receivers(self):
        soil, static, demand, surface = synthetic_case(shape=(2, 2), moisture=.35)
        mask = static.active_mask.at[0, 1].set(False)
        static = replace(static, active_mask=mask, routing_weights=jnp.ones((8, 2, 2)) * .1,
                         slope_percent=jnp.where(mask, 10., jnp.nan))
        soil = replace(soil, active_mask=mask, layer_water_cm=jnp.where(mask, soil.layer_water_cm, jnp.nan),
                       total_water_cm=jnp.where(mask, soil.total_water_cm, jnp.nan))
        _, d = self.run_case(soil, static, demand, surface)
        for name in ("lateral_in_cm", "lateral_out_cm", "lateral_export_cm"):
            np.testing.assert_array_equal(d[name][:, 0, 1], 0.)

    def test_snow_accumulation_melt_and_signed_vapor_accounting(self):
        soil, static, demand, surface = synthetic_case()
        demand = self.zero_demand(demand)
        f = lambda v: jnp.full(static.grid.shape, v)
        for snow, melt, vapor, start, end in ((1., 0., 0., 0., 1.), (0., .3, .02, 1., .72), (0., 0., -.1, 1., .9)):
            surface = SurfaceWater(f(.2 + snow), f(.2), f(snow), f(melt), f(start), f(end), f(vapor))
            self.run_case(soil, static, demand, surface)
        bad = replace(surface, swe_end_cm=f(.8))
        with self.assertRaisesRegex(FloatingPointError, "water balance failed"):
            external_water_step(soil, static, demand, bad)

    def test_missing_and_nonfinite_diagnostics_cannot_pass(self):
        soil, static, demand, surface = synthetic_case()
        after, d = self.run_case(soil, static, demand, surface)
        for key in ("vertical_down_cm", "lateral_in_cm", "bottom_drainage_cm", "interception_evap_cm", "snow_vapor_cm", "domain_residual_m3"):
            bad = dict(d)
            del bad[key]
            with self.assertRaisesRegex(ValueError, "missing required"):
                verify_water_balance(soil, after, bad, static)
        with self.assertRaises(ValueError):
            verify_water_balance(soil, after, {**d, "domain_residual_m3": jnp.asarray(jnp.nan)}, static)
        with self.assertRaisesRegex(FloatingPointError, "plant drainage view"):
            verify_water_balance(soil, after, {**d, "drain_cm": d["drain_cm"] + .1}, static)

    def test_property_state_and_demand_validation(self):
        soil, static, demand, surface = synthetic_case()
        for bad in (replace(static, soil_depth_a_cm=static.soil_depth_a_cm * 0),
                    replace(static, wiltpt_amt_cm=static.wiltpt_amt_cm * 20),
                    replace(static, sat_amt_cm=static.sat_amt_cm + 1),
                    replace(static, ksubsurface=static.ksubsurface * -1)):
            with self.assertRaises(ValueError):
                external_water_step(soil, bad, demand, surface)
        for bad in (replace(demand, transp_cm=-demand.transp_cm), replace(demand, soil_evap_cm=demand.soil_evap_cm * jnp.nan)):
            with self.assertRaises(ValueError):
                external_water_step(soil, static, bad, surface)
        with self.assertRaisesRegex(ValueError, "total/layer mismatch"):
            external_water_step(replace(soil, total_water_cm=soil.total_water_cm + 1), static, demand, surface)
        with self.assertRaises(ValueError):
            external_water_step(soil, static, demand, surface, bottom_boundary="campbell")

    def test_layer_resolution_sensitivity_with_conserved_profile(self):
        totals = []
        for depths in ((10., 90.), (30., 70.), (50., 50.)):
            soil, static, demand, surface = synthetic_case(depths=depths)
            depth = jnp.stack((static.soil_depth_a_cm, static.soil_depth_b_cm))
            demand = replace(demand, transp_cm=depth * .002, soil_evap_cm=depth * 0, interception_cm=surface.rain_cm * 0)
            after, _ = self.run_case(soil, static, demand, self.no_rain(surface), bottom_boundary="closed")
            totals.append(float(after.total_water_cm.sum()))
        np.testing.assert_allclose(totals, 19.8, atol=1e-12)
        # Hydrologic redistribution itself IS resolution sensitive: conserved
        # profile stocks do not imply matching A/B mixing or plant access.
        upper_fraction = []
        for depths in ((10., 90.), (30., 70.), (50., 50.)):
            soil, static, demand, surface = synthetic_case(depths=depths, moisture=.30)
            after, _ = self.run_case(soil, static, self.zero_demand(demand), self.no_rain(surface), bottom_boundary="closed")
            upper_fraction.append(float(after.layer_water_cm[0, 0, 0] / after.total_water_cm[0, 0]))
        self.assertGreater(max(upper_fraction) - min(upper_fraction), .1)

    def test_plant_contract_mapping_and_previous_actual_flows(self):
        soil, static, demand, surface = synthetic_case()
        after, d = self.run_case(soil, static, demand, surface)
        view = soil_to_plant(after, static, jnp.full((1, 1), 15.), d)
        self.assertEqual(view["depth_cm"].shape, (2, 1))
        np.testing.assert_array_equal(view["drain_cm"], d["drain_cm"].reshape(2, 1))
        plant = SimpleNamespace(transp_cm=demand.transp_cm.reshape(2, 1),
                                soil_evap_cm=demand.soil_evap_cm.reshape(2, 1),
                                interception_cm=demand.interception_cm.reshape(1), ground_light_frac=jnp.asarray([.5]))
        converted = demand_from_plant(plant, static)
        np.testing.assert_array_equal(converted.transp_cm, demand.transp_cm)
        with self.assertRaisesRegex(ValueError, "two layers"):
            demand_from_plant(SimpleNamespace(**{**vars(plant), "transp_cm": jnp.ones((10, 1))}), static)
        with self.assertRaisesRegex(ValueError, "previous-day"):
            soil_to_plant(after, static, jnp.full((1, 1), 15.), {})

    def test_packed_active_cell_mapping_has_explicit_row_major_order(self):
        soil, static, demand, surface = synthetic_case(shape=(2, 2))
        mask = static.active_mask.at[0, 1].set(False)
        static = replace(static, active_mask=mask)
        plant = SimpleNamespace(transp_cm=jnp.asarray([[.1, .2, .3], [.4, .5, .6]]),
                                soil_evap_cm=jnp.zeros((2, 3)), interception_cm=jnp.asarray([.01, .02, .03]),
                                ground_light_frac=jnp.ones(3))
        converted = demand_from_plant(plant, static)
        np.testing.assert_array_equal(converted.transp_cm[:, 0, 1], 0.)
        np.testing.assert_array_equal(converted.transp_cm[:, mask], plant.transp_cm)

    def test_jit_and_multiday_scan(self):
        soil, static, demand, surface = synthetic_case()
        validate_external_inputs(soil, static, demand, surface, "closed")
        step = jax.jit(partial(external_day_kernel, static=static, demand=demand, surface=surface, bottom_boundary="closed"))
        after, d = step(soil)
        verify_water_balance(soil, after, d, static)
        final, days = jax.jit(lambda initial: jax.lax.scan(lambda s, _: step(s), initial, None, length=3))(soil)
        self.assertEqual(days["transp_cm"].shape, (3, 2, 1, 1))
        self.assertTrue(np.isfinite(final.layer_water_cm).all())


if __name__ == "__main__":
    unittest.main()
