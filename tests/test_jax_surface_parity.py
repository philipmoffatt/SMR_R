"""Strict VERIFY-01 coverage: frozen fields, initial state, and failure paths."""

from dataclasses import replace
from io import BytesIO, StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch
from zipfile import ZipFile
import contextlib
import sys
import unittest

import jax.numpy as jnp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from smr_jax.verify_surface import comparison_rows


class SurfaceParityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = ROOT / "reference/mica-creek-jan2015"
        cls.rows = comparison_rows(cls.reference)

    def test_both_days_all_surface_fields_match_on_full_support(self):
        self.assertEqual(len(self.rows), 24)
        self.assertEqual({r[0] for r in self.rows}, {"2015-01-01", "2015-01-02"})
        for date, name, cells, rmse, maximum in self.rows:
            with self.subTest(date=date, field=name):
                self.assertEqual(cells, 68121)
                self.assertLessEqual(maximum, 1e-9 if name.startswith("q.") else 1e-12)

    def test_bad_support_and_nonfinite_values_cannot_pass(self):
        from smr_jax.verify_surface import compare_field
        mask = np.ones((2, 2), dtype=bool)
        expected = np.ones((2, 2), dtype=np.float64)
        for corrupt in (np.asarray([[np.nan, 1.], [1., 1.]]), np.full((2, 2), np.inf), np.ones((1, 2))):
            with self.subTest(shape=corrupt.shape):
                with self.assertRaises(ValueError):
                    compare_field(corrupt, expected, mask, name="swe")
        with self.assertRaises(ValueError):
            compare_field(expected, expected * np.nan, mask, name="swe")
        with self.assertRaises(ValueError):
            compare_field(expected, expected, mask * False, name="swe")
        with self.assertRaises(ValueError):
            compare_field(expected.astype(np.float32), expected, mask, name="swe")

    def test_wrong_flux_fails_and_report_only_remains_honest(self):
        from smr_jax.verify_surface import compare_field
        expected = np.ones((1, 2), dtype=np.float64)
        mask = np.ones((1, 2), dtype=bool)
        with self.assertRaises(AssertionError):
            compare_field(expected + .01, expected, mask, name="swe")
        cells, rmse, maximum = compare_field(expected + .01, expected, mask, name="swe", assert_close=False)
        self.assertEqual(cells, 2)
        self.assertAlmostEqual(maximum, .01)

    def test_replayed_shortwave_is_compared_exactly(self):
        from smr_jax.verify_surface import compare_field
        expected = np.ones((1, 1), dtype=np.float64)
        with self.assertRaises(AssertionError):
            compare_field(np.nextafter(expected, np.inf), expected, np.ones((1, 1), dtype=bool), name="q.srad")

    def test_strict_cpu_division_preserves_nearest_double_results(self):
        from smr_jax.fluxes.legacy_math import compile_legacy_cpu
        values = np.linspace(.1, 300., 1001, dtype=np.float64)
        actual = compile_legacy_cpu(lambda a: a / 130.)(jnp.asarray(values))
        np.testing.assert_array_equal(actual, values / 130.)

    def test_rain_snow_partition_at_both_thresholds(self):
        from smr_jax.fluxes.legacy_math import compile_legacy_cpu
        from smr_jax.fluxes.surface import precipitation_and_temperature
        field = lambda v: jnp.full((1, 1), v, dtype=jnp.float64)
        static = SimpleNamespace(elevation_m=field(784.), tmin_snow_c=field(-1.), tmax_rain_c=field(1.))

        def evaluate(tavg):
            forcing = SimpleNamespace(tavg_c=tavg, tdew_c=jnp.asarray(0.),
                                      precip_cm=jnp.asarray(1.), pet_cm=jnp.asarray(.3))
            return precipitation_and_temperature(static, forcing)

        run = compile_legacy_cpu(evaluate)
        for temperature, rain_share in ((-2., 0.), (-1., 0.), (0., .5), (1., 1.), (2., 1.)):
            values = run(jnp.asarray(temperature))
            np.testing.assert_allclose(values["rain_cm"], rain_share * np.asarray(values["precip_cm"]), rtol=0, atol=1e-14)
            np.testing.assert_allclose(values["rain_cm"] + values["snow_cm"], values["precip_cm"], rtol=0, atol=1e-14)

    def test_canopy_water_and_pet_ledgers(self):
        from smr_jax.fluxes.legacy_math import compile_legacy_cpu
        from smr_jax.fluxes.surface import canopy
        water = jnp.asarray([[.02, .10, 0., .05]])
        rain = jnp.asarray([[0., .10, .02, 0.]])
        pet = jnp.asarray([[.01, .02, .04, 0.]])
        static = SimpleNamespace(max_canopy_cm=jnp.full((1, 4), .1))
        run = compile_legacy_cpu(lambda stored, r, p: canopy(SimpleNamespace(canopy_storage_cm=stored), static, r, p))
        after, soil_pet, d = run(water, rain, pet)
        np.testing.assert_allclose(after, [[.01, .08, 0., .05]], rtol=0, atol=1e-14)
        np.testing.assert_allclose(after - water - rain + d["throughfall_cm"] + d["canopy_evap_cm"], 0., rtol=0, atol=1e-14)
        np.testing.assert_allclose(soil_pet + d["canopy_evap_cm"], pet, rtol=0, atol=1e-14)

    def test_prescribed_energy_melt_and_refreeze(self):
        from smr_jax.fluxes.legacy_math import compile_legacy_cpu
        from smr_jax.fluxes.snow import snow_energy
        from smr_jax.state import SnowState
        field = lambda v: jnp.full((1, 1), v, dtype=jnp.float64)
        static = SimpleNamespace(active_mask=jnp.ones((1, 1), dtype=bool), canopy_cover=field(1.))
        # At 0 C vapor and sensible heat cancel. Supply radiation that cancels
        # the known longwave loss plus a prescribed latent-energy amount.
        longwave_loss = 5.67e-8 * 86400. / 1000. * (.98 - .92) * 273.15**4
        for ice, liquid, latent_cm in ((1., 0., .25), (.5, .1, -.05)):
            with self.subTest(latent_cm=latent_cm):
                before = SnowState(field(ice + liquid), field(ice + liquid), field(ice), field(liquid),
                                   field(1.), field(0.), field(0.), field(0.))
                forcing = SimpleNamespace(roughness_urban=jnp.asarray(100.), roughness_crop=jnp.asarray(10.),
                                          cloud=jnp.asarray(0.), replay_shortwave=field(longwave_loss + 3350. * latent_cm))
                run = compile_legacy_cpu(lambda snow: snow_energy(SimpleNamespace(snow=snow), static, forcing,
                                          field(0.), field(0.), field(0.), field(0.)))
                after, d = run(before)
                expected_ice = ice - latent_cm
                expected_liquid = min(.03 / .97 * expected_ice, liquid + latent_cm)
                np.testing.assert_allclose(after.ice_cm, expected_ice, rtol=0, atol=1e-12)
                np.testing.assert_allclose(after.liquid_cm, expected_liquid, rtol=0, atol=1e-12)
                np.testing.assert_allclose(after.swe_cm - before.swe_cm - d["condens_cm"] + d["snowmelt_cm"],
                                           0., rtol=0, atol=1e-12)

    def test_initialization_state_and_properties_are_checked(self):
        from smr_jax.initialization import initialize_legacy_state
        from smr_jax.legacy import load_legacy_bundle
        from smr_jax.verify_surface import verify_initialization
        static = load_legacy_bundle(self.reference / "model_inputs.zip")
        state = initialize_legacy_state(static)
        verify_initialization(static, state)
        with self.assertRaises(AssertionError):
            verify_initialization(static, replace(state, snow=replace(state.snow, swe_cm=state.snow.swe_cm + .01)))
        with self.assertRaises(AssertionError):
            verify_initialization(replace(static, canopy_cover=static.canopy_cover + .01), state)

    def test_corrupt_frozen_manifest_cannot_pass(self):
        from smr_jax.verify_surface import verify_reference_manifest
        with TemporaryDirectory() as folder:
            path = Path(folder)
            (path / "manifest.json").write_text('{"reference_id":"wrong","artifacts":[]}', encoding="utf-8")
            with self.assertRaises(ValueError):
                verify_reference_manifest(path)

    def test_missing_diagnostic_and_changed_reference_support_fail(self):
        from smr_jax.io import FrozenTraceProvider
        with patch("smr_jax.verify_surface._run_legacy_cpu", return_value=(None, {})):
            with self.assertRaisesRegex(ValueError, "missing required diagnostic"):
                comparison_rows(self.reference)
        original = FrozenTraceProvider.raster

        def missing_cell(provider, date, name):
            value = original(provider, date, name)
            if name == "saturation":
                value.ravel()[np.flatnonzero(np.isfinite(value).ravel())[0]] = np.nan
            return value

        with patch("smr_jax.verify_surface._run_legacy_cpu", return_value=(None, {})), patch.object(
                FrozenTraceProvider, "raster", new=missing_cell):
            with self.assertRaisesRegex(ValueError, "required.*cells"):
                comparison_rows(self.reference)

    def test_subtolerance_swe_cannot_hide_a_changed_day_two_branch(self):
        from smr_jax.io import FrozenTraceProvider
        original = FrozenTraceProvider.raster

        def altered_zero(provider, date, name):
            value = original(provider, date, name)
            if date == "2015-01-01" and name == "swe":
                value.ravel()[np.flatnonzero((value == 0.).ravel())[0]] = np.nextafter(0., 1.)
            return value

        with patch.object(FrozenTraceProvider, "raster", new=altered_zero):
            with self.assertRaisesRegex(AssertionError, "SWE decisions"):
                comparison_rows(self.reference)

    def test_cli_nonzero_on_failure_and_policy_requires_explicit_report(self):
        from smr_jax.verify_surface import main
        with patch("smr_jax.verify_surface.comparison_rows", side_effect=AssertionError("bad flux")):
            with contextlib.redirect_stdout(StringIO()), contextlib.redirect_stderr(StringIO()):
                self.assertEqual(main([]), 1)
        with contextlib.redirect_stdout(StringIO()), contextlib.redirect_stderr(StringIO()):
            with self.assertRaises(SystemExit) as error:
                main(["--melt-small-snow-next-day"])
            self.assertEqual(error.exception.code, 2)
        with patch("smr_jax.verify_surface.comparison_rows", return_value=self.rows) as report:
            output = StringIO()
            with contextlib.redirect_stdout(output):
                self.assertEqual(main(["--report-only", "--melt-small-snow-next-day"]), 0)
            self.assertIn("REPORT ONLY", output.getvalue())
            self.assertNotIn("PASS:", output.getvalue())
            self.assertFalse(report.call_args.kwargs["assert_close"])


class TraceReaderFailureTest(unittest.TestCase):
    def archive(self, folder, *, origin=100., pixel=30., nodata=False):
        import tifffile
        from smr_jax.state import GridSpec
        from smr_jax.io import FrozenTraceProvider
        buffer = BytesIO()
        values = np.ones((2, 2), dtype=np.float64)
        tags = [(33550, "d", 3, (pixel, pixel, 0.), False),
                (33922, "d", 6, (0., 0., 0., origin, 260., 0.), False)]
        if nodata:
            values[0, 0] = -9999.
            tags.append((42113, "s", 0, "-9999", False))
        tifffile.imwrite(buffer, values, extratags=tags)
        path = Path(folder) / "trace.zip"
        with ZipFile(path, "w") as archive:
            archive.writestr("trace\\2015-01-01--precip.tif", buffer.getvalue())
        return FrozenTraceProvider(path, GridSpec(2, 2, 30., 100., 200.))

    def test_missing_field_has_an_explicit_error(self):
        with TemporaryDirectory() as folder:
            with self.assertRaises(KeyError):
                self.archive(folder).raster("2015-01-02", "swe")

    def test_fractional_cell_shift_and_resolution_mismatch_fail(self):
        for options in ({"origin": 115.}, {"pixel": 31.}):
            with self.subTest(options=options), TemporaryDirectory() as folder:
                with self.assertRaises(ValueError):
                    self.archive(folder, **options).raster("2015-01-01", "precip")

    def test_encoded_nodata_is_not_counted_as_finite_support(self):
        with TemporaryDirectory() as folder:
            values = self.archive(folder, nodata=True).raster("2015-01-01", "precip")
            self.assertTrue(np.isnan(values[0, 0]))
            self.assertEqual(np.isfinite(values).sum(), 3)


if __name__ == "__main__":
    unittest.main()
